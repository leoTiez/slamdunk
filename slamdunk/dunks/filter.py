#!/usr/bin/env python

# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
#
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Date located in: -
from __future__ import print_function
import pysam, random, os

from slamdunk.version import __version__, __bam_version__  # @UnresolvedImport

from slamdunk.utils.BedReader import bedToIntervallTree  # @UnresolvedImport
from slamdunk.utils.misc import checkStep, run, removeFile, getBinary, pysamIndex, SlamSeqInfo, md5  # @UnresolvedImport


def bamSort(outputBAM, log, newHeader, paired, verbose):
    tmp = outputBAM + "_tmp"
    if newHeader is not None:
        pyOutputBAM = pysam.AlignmentFile(outputBAM, "rb")
        pyTmp = pysam.AlignmentFile(tmp, "wb", header=newHeader)
        for read in pyOutputBAM:
            pyTmp.write(read)
        pyOutputBAM.close()
        pyTmp.close()
    else:
        os.rename(outputBAM, tmp)

    if not paired:
        run("samtools sort %s -o %s" % (tmp, outputBAM), log, verbose=verbose, dry=False)
    else:
        run("samtools sort -n %s -o %s" % (tmp, outputBAM), log, verbose=verbose)
    removeFile(tmp)


def dumpBufferToBam (buffer, multimapList, outbam, infile):
    # Randomly write hit from read
    read = list(buffer.values()).pop().pop()
    read.set_tag("RD", multimapList.rstrip(" "), "Z")
    read.is_secondary = False
    read.is_supplementary = False
    outbam.write(read)


def multimapUTRRetainment (infile, outfile, bed, minIdentity, NM, MQ, log):
    mappedReads, unmappedReads, filteredReads, mqFiltered, idFiltered, nmFiltered = 0, 0, 0, 0, 0, 0
    utrIntervallTreeDict = bedToIntervallTree(bed)  # Is interpreted now to represent any entity in the bed file

    # Buffers for multimappers
    multimapBuffer = {}
    prevRead = ""
    # If read maps to another than previously recorded UTR -> do not dump reads to file
    dumpBuffer = True
    # This string tracks all multiple alignments
    multimapList = ""

    for read in infile:  # infile is AlignedFile according to pysam definition. Read is AlignedSegment
        if not read.is_secondary and not read.is_supplementary:
            if read.is_unmapped:
                unmappedReads += 1
            else:
                mappedReads += 1

        # First pass general filters
        if read.is_unmapped:
            continue
        if float(read.get_tag("XI")) < minIdentity:
            idFiltered += 1
            continue
        if -1 < NM < int(read.get_tag("NM")):
            nmFiltered += 1
            continue
        if read.mapping_quality < MQ:
            # Previous read was also multimapper
            if read.query_name != prevRead and prevRead != "":
                if dumpBuffer and len(multimapBuffer) > 0:
                    dumpBufferToBam(multimapBuffer, multimapList, outfile, infile)
                    filteredReads += 1

                dumpBuffer = True
                multimapList = ""
                multimapBuffer = {}

            # Query Intervall tree for given chromosome for UTRs
            chr = infile.get_reference_name(read.reference_id)
            start = read.reference_start
            end = read.reference_end

            if chr in utrIntervallTreeDict:
                query = utrIntervallTreeDict[chr][start:end]  # This makes sure that mapping is in a bed region
            else:
                query = set()

            if len(query) > 0:
                # First UTR hit is recorded without checks
                if len(multimapBuffer) == 0:
                    for result in query:
                        if result.data not in multimapBuffer:
                            multimapBuffer[result.data] = []
                        multimapBuffer[result.data].append(read)
                # Second UTR hit looks at previous UTR hits -> no dump if hit on different UTR
                else:
                    for result in query:
                        if result.data not in multimapBuffer:
                            multimapBuffer[result.data] = []
                            multimapBuffer[result.data].append(read)
                            dumpBuffer = False
                        else:
                            multimapBuffer[result.data].append(read)

            multimapList = multimapList + chr + ":" + str(start) + "-" + str(end) + " "
            prevRead = read.query_name
        else:  # If read.mapping_quality > mq
            # Dump any multimappers before a unique mapper
            if len(multimapBuffer) > 0:
                if dumpBuffer:
                    dumpBufferToBam(multimapBuffer, multimapList, outfile, infile)
                    filteredReads += 1
                multimapBuffer = {}
                dumpBuffer = True
                multimapList = ""

            # Record all unique mappers
            prevRead = read.query_name
            outfile.write(read)
            filteredReads += 1

    # Dump last portion if it was multimapper
    if dumpBuffer and len(multimapBuffer) > 0:
        dumpBufferToBam(multimapBuffer, multimapList, outfile, infile)
        filteredReads += 1

    multimapper = mappedReads - filteredReads - idFiltered - nmFiltered

    print("Criterion\tFiltered reads", file=log)
    print("MQ < 0\t0", file=log)
    print("ID < %s\t%s" % (minIdentity, idFiltered), file=log)
    print("NM > %s\t%s" % (NM, nmFiltered), file=log)
    print("MM\t%s" % multimapper, file=log)

    return mappedReads, unmappedReads, filteredReads, mqFiltered, idFiltered, nmFiltered, multimapper


def Filter(
        inputBAM,
        outputBAM,
        log,
        bed,
        MQ=2,
        minIdentity=0.8,
        NM=-1,
        printOnly=False,
        verbose=True,
        force=False,
        paired=False
):
    inputBAM = os.path.expanduser(inputBAM)
    outputBAM = os.path.expanduser(outputBAM)
    if printOnly or checkStep([inputBAM], [outputBAM], force):
        (mappedReads, unmappedReads, filteredReads,
         mqFiltered, idFiltered, nmFiltered, multimapper) = 0, 0, 0, 0, 0, 0, 0

        infile = pysam.AlignmentFile(inputBAM, "rb")
        outfile = pysam.AlignmentFile(outputBAM, "wb", template=infile)
        # Default filtering without bed
        if bed is None:
            print("#No bed-file supplied. Running default filtering on " + inputBAM + ".", file=log)
            if paired:
                read1 = None
                read2 = None
            for read in infile:
                if paired:
                    if not read.is_paired or read.mate_is_unmapped or read.is_duplicate:
                        unmappedReads += 1
                        continue
                    if read.is_read2:
                        read2 = read
                    else:
                        read1 = read
                        read2 = None
                        continue

                if not read.is_secondary and not read.is_supplementary:
                    if read.is_unmapped:
                        unmappedReads += 1
                        continue
                    else:
                        mappedReads += 1

                if not paired:
                    if read.mapping_quality < MQ:
                        mqFiltered += 1
                        continue
                    if float(read.get_tag("XI")) < minIdentity:
                        idFiltered += 1
                        continue
                    if -1 < NM < int(read.get_tag("NM")):
                        nmFiltered += 1
                        continue

                    filteredReads += 1
                    outfile.write(read)
                else:
                    if read1 is None or read2 is None:
                        continue
                    if read1.query_name != read2.query_name:
                        continue

                    if read1.mapping_quality < MQ and read2.mapping_quality < MQ:
                        mqFiltered += 1
                        continue
                    if float(read1.get_tag("XI")) < minIdentity and float(read2.get_tag("XI")) < minIdentity:
                        idFiltered += 1
                        continue
                    if -1 < NM < int(read1.get_tag("NM")) and -1 < NM < int(read2.get_tag("NM")):
                        nmFiltered += 1
                        continue
                    filteredReads += 1
                    outfile.write(read1)
                    outfile.write(read2)

            print("Criterion\tFiltered reads", file=log)
            print("MQ < 0\t0", file=log)
            print("ID < %s\t%s" % (minIdentity, idFiltered), file=log)
            print("NM > %s\t%s" % (NM, nmFiltered), file=log)
            print("MM\t0", file=log)
        else:
            # Multimap retention strategy filtering when bed is supplied
            print("#Bed-file supplied. Running multimap retention filtering strategy on " + inputBAM + ".", file=log)
            (
                mappedReads,
                unmappedReads,
                filteredReads,
                mqFiltered,
                idFiltered,
                nmFiltered,
                multimapper
            ) = multimapUTRRetainment(infile, outfile, bed, minIdentity, NM, MQ, log)

        # Add number of sequenced and number of mapped reads to the read group description
        # Used for creating summary file
        inFileBamHeader = outfile.header
        if "RG" in inFileBamHeader and len(inFileBamHeader["RG"]) > 0:
            slamseqInfo = SlamSeqInfo()
            slamseqInfo.SequencedReads = mappedReads + unmappedReads
            slamseqInfo.MappedReads = mappedReads
            slamseqInfo.FilteredReads = filteredReads
            slamseqInfo.MQFilteredReads = mqFiltered
            slamseqInfo.IdFilteredReads = idFiltered
            slamseqInfo.NmFilteredReads = nmFiltered
            slamseqInfo.MultimapperReads = multimapper

            if bed:
                slamseqInfo.AnnotationName = os.path.basename(bed)
                slamseqInfo.AnnotationMD5 = md5(bed)
            else:
                slamseqInfo.AnnotationName = ""
                slamseqInfo.AnnotationMD5 = ""

            if not isinstance(inFileBamHeader, dict):
                inFileBamHeader = inFileBamHeader.to_dict()
            inFileBamHeader["RG"][0]["DS"] = str(slamseqInfo)

        slamDunkPG = {"ID": "slamdunk", "PN": "slamdunk filter v" + __version__, "VN": __bam_version__}
        if "PG" in inFileBamHeader:
            inFileBamHeader["PG"].append(slamDunkPG)
        else:
            inFileBamHeader["PG"] = [slamDunkPG]

        infile.close()
        outfile.close()

        # Sort afterwards
        bamSort(outputBAM, log, inFileBamHeader, paired=False, verbose=verbose)
        if not paired:
            pysamIndex(outputBAM)
    else:
        print("Skipped filtering for " + inputBAM, file=log)
