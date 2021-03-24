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

#########################################################################
# Main routine for the SLAMdunk analyzer
#########################################################################
# Imports
#########################################################################
from __future__ import print_function
import sys, os, random

from time import sleep
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS

from os.path import basename

from joblib import Parallel, delayed
from slamdunk.dunks import tcounter, mapper, filter, deduplicator, snps
from slamdunk.utils.misc import replaceExtension, estimateMaxReadLength
from slamdunk.version import __version__

########################################################################
# Global variables
########################################################################

printOnly = False
verbose = False

mainOutput = sys.stderr

logToMainOutput = False

########################################################################
# Routine definitions
########################################################################


def getLogFile(path):
    if(logToMainOutput):
        return mainOutput
    else:
        log = open(path, "a")
        return log


def closeLogFile(log):
    if(not logToMainOutput):
        log.close()


def message(msg):
    print(msg, file=mainOutput)


def error(msg, code=-1):
    print(msg, file=mainOutput)
    sys.exit(code)


def stepFinished():
    print(".", end="", file=mainOutput)


def dunkFinished():
    print("", file=mainOutput)


def createDir(directory):
    if not os.path.exists(directory):
        message("Creating output directory: " + directory)
        os.makedirs(directory)


def readSampleFile(fileName):
    samples = []
    infos = []

    with open(fileName, "r") as ins:
        for line in ins:
            line = line.strip()
            if(len(line) > 1):
                if(fileName.endswith(".tsv")):
                    cols = line.split("\t")
                elif(fileName.endswith(".csv")):
                    cols = line.split(",")
                else:
                    raise RuntimeError("Unknown file extension found: " + fileName)
                if(len(cols) < 4):
                    raise RuntimeError("Invalid sample file found: " + fileName)
                samples.append(cols[0])
                infos.append(cols[1] + ":" + cols[2] + ":" + cols[3])

    return samples, infos


def getSamples(bams, runOnly=-1):
    samples = []
    samplesInfos = []
    if len(bams) == 1 and (bams[0].endswith(".tsv") or bams[0].endswith(".csv")):
        # Sample file specified
        samples, samplesInfos = readSampleFile(bams[0])
    else:
        # List of BAM files specified
        samples = bams
        samplesInfos = [""] * len(samples)

    if runOnly > 0:
        if runOnly > len(samples):
            raise RuntimeError("Sample index out of range. %s > %s. Check -i/--sample-index" % (runOnly, len(samples)))
        message("Running only job " + str(runOnly))
        samples = [samples[runOnly - 1]]
        samplesInfos = [samplesInfos[runOnly - 1]]
    elif runOnly == 0:
        raise RuntimeError("Sample index (" + str(runOnly) + ") out of range. Starts with 1. Check -i/--sample-index")

    return samples, samplesInfos


def processPrepare(args):
    if args.sampleIndex > -1:
        sec = random.randrange(200, 2000) / 1000.0
        message("Waiting " + str(sec) + " seconds")
        sleep(sec)

    # Setup slamdunk run folder
    outputDirectory = args.outputDir
    createDir(outputDirectory)

    n = args.threads
    referenceFile = args.referenceFile
    samples, samplesInfos = getSamples(args.files, runOnly=args.sampleIndex)
    return outputDirectory, n, referenceFile, samples, samplesInfos


def processMap(
        samples,
        paired,
        sampleName,
        sampleType,
        sampleTime,
        sampleIndex,
        trim5,
        maxPolyA,
        quantseq,
        endtoend,
        topn,
        skipSAM,
        outputDirectory,
        samplesInfos,
        referenceFile,
        n
):
    mapper.checkNextGenMapVersion()
    dunkPath = os.path.join(outputDirectory, "map")
    createDir(dunkPath)

    message("Running slamDunk map for " + str(len(samples)) + " files (" + str(n) + " threads)")
    if paired:
        message("Assuming files to be paired, i.e. representing forward and reverse reads.")
        if len(samples) != 2:
            raise ValueError("Files are required to be paired, but not 2 samples passed."
                             "Several paired processing is not supported. Please pass only 2 files.")
        if not sampleName:
            sampleName = replaceExtension(basename(samples[0]), "", "")
        else:
            sampleName = sampleName

        sampleInfo = samplesInfos[0]
        if sampleInfo == "":
            sampleInfo = "%s:%s:%s" % (sampleName, sampleType, str(sampleTime))

        runMap(
            None,
            samples,
            referenceFile,
            n,
            trim5,
            maxPolyA,
            quantseq,
            endtoend,
            topn,
            sampleInfo,
            dunkPath,
            skipSAM,
            isPaired=True
        )

    else:
        message("Assuming files to be separate, i.e. they are mapped sequentially.")
        for i in range(0, len(samples)):
            bam = samples[i]

            if not sampleName or len(samples) > 1:
                sampleName = replaceExtension(basename(bam), "", "")
            else:
                sampleName = sampleName

            sampleInfo = samplesInfos[i]
            if sampleInfo == "":
                sampleInfo = "%s:%s:%s" % (sampleName, sampleType, str(sampleTime))
            tid = i
            if sampleIndex > -1:
                tid = sampleIndex
            runMap(
                tid,
                bam,
                referenceFile,
                n,
                trim5,
                maxPolyA,
                quantseq,
                endtoend,
                topn,
                sampleInfo,
                dunkPath,
                skipSAM
            )

    dunkFinished()
    dunkbufferIn = []
    if not paired:
        if not skipSAM:
            message("Running slamDunk sam2bam for %s files (%s threads)" % (len(samples), n))
            _ = Parallel(n_jobs=1, verbose=verbose)(delayed(runSam2Bam)(
                samples[tid],
                n,
                dunkPath,
            ) for tid in range(0, len(samples)))
            dunkFinished()

        for file in samples:
            dunkbufferIn.append(os.path.join(dunkPath, replaceExtension(basename(file), ".bam", "_slamdunk_mapped")))
    else:
        if not skipSAM:
            message("Running slamDunk sam2bam for the paired sam file (%s threads)" % n)
            runSam2Bam(samples[0], n, dunkPath, isPaired=True)
            dunkFinished()
        dunkbufferIn.append(os.path.join(
            dunkPath,
            replaceExtension(basename(samples[0]), ".bam", "_slamdunk_mapped_paired")
        ))

    return dunkbufferIn


def processFilter(bam, mq, identity, nm, bed, paired, outputDirectory, n):
    dunkPath = os.path.join(outputDirectory, "filter")
    createDir(dunkPath)
    message("Running slamDunk filter for %s files (%s threads)" % (len(bam), n))
    _ = Parallel(n_jobs=n, verbose=verbose)(
        delayed(runFilter)(bam[tid], bed, mq, identity, nm, paired, dunkPath) for tid in range(0, len(bam)))
    dunkFinished()
    return dunkPath


def processSNP(bam, fasta, minCov, minVarFreq, minQual, outputDirectory, n):
    if n > 1:
        n = n // 2
    dunkPath = os.path.join(outputDirectory, "snp")
    createDir(dunkPath)
    message("Running slamDunk SNP for %s files (%s threads)" % (len(bam), n))
    _ = Parallel(n_jobs=n, verbose=verbose)(
        delayed(runSnp)(fasta, minCov, minVarFreq, minQual, bam[tid], dunkPath) for tid in range(0, len(bam)))
    dunkFinished()
    return dunkPath


def processCount(
        bams,
        ref,
        bed,
        maxLength,
        minQual,
        conversionThreshold,
        snpDirectory,
        vcfFile,
        outputDirectory,
        n
):
    dunkPath = os.path.join(outputDirectory, "count")
    createDir(dunkPath)

    message("Running slamDunk count for %s files (%s threads)" % (len(bams), n))
    _ = Parallel(n_jobs=n, verbose=verbose)(
        delayed(runCount)(
            bams[tid],
            ref,
            bed,
            maxLength,
            minQual,
            conversionThreshold,
            dunkPath,
            snpDirectory,
            vcfFile
        ) for tid in range(0, len(bams))
    )
    dunkFinished()


def runMap(
        tid,
        inputBAM,
        referenceFile,
        threads,
        trim5p,
        maxPolyA,
        quantseqMapping,
        endtoendMapping,
        topn,
        sampleDescription,
        outputDirectory,
        skipSAM,
        isPaired=False
):
    extension = "_slamdunk_mapped%s" % "" if not isPaired else "_paired"
    base_name = basename(inputBAM) if not isPaired else basename(inputBAM[0])

    if skipSAM:
        outputSAM = os.path.join(outputDirectory, replaceExtension(base_name, ".bam", extension))
    else:
        outputSAM = os.path.join(outputDirectory, replaceExtension(base_name, ".sam", extension))
    outputLOG = os.path.join(outputDirectory, replaceExtension(base_name, ".log", extension))

    sampleName = replaceExtension(base_name, ".bam", "")
    sampleType = "NA"
    sampleTime = "-1"
    if sampleDescription != "":
        sampleDescriptions = sampleDescription.split(":")
        if len(sampleDescriptions) >= 1:
            sampleName = sampleDescriptions[0]
        if len(sampleDescriptions) >= 2:
            typeDict = {"p": "pulse", "c": "chase", "pulse": "pulse", "chase": "chase", "": "NA"}
            if sampleDescriptions[1] in typeDict:
                sampleType = typeDict[sampleDescriptions[1]]
            else:
                sampleType = sampleDescriptions[1]
        if len(sampleDescriptions) >= 3:
            sampleTime = sampleDescriptions[2]

    mapper.Map(
        inputBAM,
        referenceFile,
        outputSAM,
        getLogFile(outputLOG),
        quantseqMapping,
        endtoendMapping,
        threads=threads,
        trim5p=trim5p,
        maxPolyA=maxPolyA,
        topn=topn,
        sampleId=tid,
        sampleName=sampleName,
        sampleType=sampleType,
        sampleTime=sampleTime,
        printOnly=printOnly,
        verbose=verbose,
        isPaired=isPaired
    )
    stepFinished()


def runSam2Bam(bam, threads, outputDirectory, isPaired=False):
    inputSAM = os.path.join(
        outputDirectory,
        replaceExtension(basename(bam), ".sam", "_slamdunk_mapped%s" % "" if not isPaired else "_paired")
    )
    outputBAM = os.path.join(
        outputDirectory,
        replaceExtension(basename(bam), ".bam", "_slamdunk_mapped%s" % "" if not isPaired else "_paired")
    )
    outputLOG = os.path.join(
        outputDirectory,
        replaceExtension(basename(bam), ".log", "_slamdunk_mapped%s" % "" if not isPaired else "_paired")
    )
    mapper.sort(inputSAM, outputBAM, getLogFile(outputLOG), threads, False, printOnly, verbose, isPaired)
    stepFinished()


def runDedup(bam, outputDirectory) :
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_dedup"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_dedup"))
    log = getLogFile(outputLOG)
    deduplicator.Dedup(bam, outputBAM, log)
    closeLogFile(log)
    stepFinished()


def runFilter(bam, bed, mq, minIdentity, maxNM, paired, outputDirectory):
    outputBAM = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bam", "_filtered"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_filtered"))
    filter.Filter(bam, outputBAM, getLogFile(outputLOG), bed, mq, minIdentity, maxNM, printOnly, verbose, paired=paired)
    stepFinished()


def runSnp(referenceFile, minCov, minVarFreq, minQual, inputBAM, outputDirectory) :
    outputSNP = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".vcf", "_snp"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(inputBAM), ".log", "_snp"))
    snps.SNPs(
        inputBAM,
        outputSNP,
        referenceFile,
        minVarFreq,
        minCov,
        minQual,
        getLogFile(outputLOG),
        printOnly,
        verbose,
        False
    )
    stepFinished()


def runCount(
        bam,
        ref,
        bed,
        maxLength,
        minQual,
        conversionThreshold,
        outputDirectory,
        snpDirectory,
        vcfFile,
):
    outputCSV = os.path.join(outputDirectory, replaceExtension(basename(bam), ".tsv", "_tcount"))
    outputBedgraphPlus = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bedgraph", "_tcount_plus"))
    outputBedgraphMinus = os.path.join(outputDirectory, replaceExtension(basename(bam), ".bedgraph", "_tcount_mins"))
    outputLOG = os.path.join(outputDirectory, replaceExtension(basename(bam), ".log", "_tcount"))

    if vcfFile is not None:
        inputSNP = vcfFile
    elif snpDirectory is not None:
        inputSNP = os.path.join(snpDirectory, replaceExtension(basename(bam), ".vcf", "_snp"))
    else:
        inputSNP = None

    if maxLength is None:
        maxLength = estimateMaxReadLength(bam)
    if maxLength < 0:
        print("Difference between minimum and maximum read length is > 10. "
              "Please specify --max-read-length parameter.")
        sys.exit(0)

    log = getLogFile(outputLOG)
    print("Using " + str(maxLength) + " as maximum read length.", file=log)
    if bed is not None:
        message("Bed file detected.")
        tcounter.computeTconversions(
            ref,
            bed,
            inputSNP,
            bam,
            maxLength,
            minQual,
            outputCSV,
            outputBedgraphPlus,
            outputBedgraphMinus,
            conversionThreshold,
            log
        )
    else:
        message("No bed file passed. Count w.r.t. the full genome.")
        outputBedgraphPlusNew = os.path.join(
            outputDirectory,
            replaceExtension(basename(bam), ".bedgraph", "_tcount_plus_new")
        )
        outputBedgraphMinusNew = os.path.join(
            outputDirectory,
            replaceExtension(basename(bam), ".bedgraph", "_tcount_mins_new")
        )
        tcounter.computeTconversionsAll(
            ref,
            inputSNP,
            bam,
            outputBedgraphPlus,
            outputBedgraphPlusNew,
            outputBedgraphMinus,
            outputBedgraphMinusNew,
            conversionThreshold,
            minQual,
            log
        )
    stepFinished()
    return outputCSV


def runAll(args):
    message("slamdunk all")
    outputDirectory, n, referenceFile, samples, samplesInfos = processPrepare(args)
    # Run mapper dunk
    dunkbufferIn = processMap(
        samples,
        args.paired,
        args.sampleName,
        args.sampleType,
        args.sampleTime,
        args.sampleIndex,
        args.trim5,
        args.maxPolyA,
        args.quantseq,
        args.endtoend,
        args.topn,
        args.skipSAM,
        outputDirectory,
        samplesInfos,
        referenceFile,
        n
    )

    # Run filter dunk
    bed = args.bed
    if args.filterbed:
        bed = args.filterbed
        args.multimap = True
    if not args.multimap:
        bed = None
    bam = args.bam
    mq = args.mq
    identity = args.identity
    nm = args.nm
    dunkPath = processFilter(bam, mq, identity, nm, bed, args.paired, outputDirectory, n)
    dunkbufferOut = []

    for file in dunkbufferIn:
        dunkbufferOut.append(os.path.join(dunkPath, replaceExtension(basename(file), ".bam", "_filtered")))

    dunkbufferIn = dunkbufferOut
    dunkbufferOut = []
    dunkFinished()

    # Run snps dunk only if vcf not specified
    snpDirectory = None
    vcfFile = None

    if not "vcfFile" in args:
        snpDirectory = processSNP(bam, referenceFile, args.cov, args.var, args.minQual, outputDirectory, n)
        dunkFinished()
    else:
        vcfFile = args.vcfFile

    # Run count dunk
    processCount(
        dunkbufferIn,
        referenceFile,
        bed,
        args.maxLength,
        args.minQual,
        args.conversionThreshold,
        snpDirectory,
        vcfFile,
        outputDirectory,
        n
    )


def run():
    ########################################################################
    # Argument parsing
    ########################################################################

    # Info
    usage = "SLAMdunk software for analyzing SLAM-seq data"

    # Main Parsers
    parser = ArgumentParser(description=usage, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    # Initialize Subparsers
    subparsers = parser.add_subparsers(help="", dest="command")

    # map command
    mapparser = subparsers.add_parser("map", formatter_class=ArgumentDefaultsHelpFormatter,
                                      help="Map SLAM-seq read data",)
    mapparser.add_argument("files", action="store", nargs="+",
                           help="Single csv/tsv file (recommended) containing all sample files and sample info or a "
                                "list of all sample BAM/FASTA(gz)/FASTQ(gz) files")
    mapparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile", default=SUPPRESS,
                           help="Reference fasta file")
    mapparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", default=SUPPRESS,
                           help="Output directory for mapped BAM files.")
    mapparser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", default=12,
                           help="Number of bp removed from 5' end of all reads.")
    mapparser.add_argument("-n", "--topn", type=int, required=False, dest="topn", default=1,
                           help="Max number of alignments to report per read")
    mapparser.add_argument("-a", "--max-polya", type=int, required=False, dest="maxPolyA", default=4,
                           help="Max number of As at the 3' end of a read.")
    mapparser.add_argument("-t", "--threads", type=int, required=False, dest="threads", default=1,
                           help="Thread number")
    mapparser.add_argument("-q", "--quantseq", dest="quantseq", action="store_true", required=False,
                           help="Run plain Quantseq alignment without SLAM-seq scoring")
    mapparser.add_argument("-e", "--endtoend", action="store_true", dest="endtoend",
                           help="Use a end to end alignment algorithm for mapping.")
    mapparser.add_argument("-sn", "--sampleName", type=str, dest="sampleName", required=False,
                           help="Use this sample name for all supplied samples")
    mapparser.add_argument("-sy", "--sampleType", type=str, dest="sampleType", required=False, default="pulse",
                           help="Use this sample type for all supplied samples")
    mapparser.add_argument("-st", "--sampleTime", type=int, dest="sampleTime", required=False, default=0,
                           help="Use this sample time for all supplied samples")
    mapparser.add_argument("-i", "--sample-index", type=int, required=False, default=-1, dest="sampleIndex",
                           help="Run analysis only for sample <i>. Use for distributing slamdunk analysis on a cluster"
                                " (index is 1-based).")
    mapparser.add_argument("-ss", "--skip-sam", action="store_true", dest="skipSAM",
                           help="Output BAM while mapping. Slower but, uses less hard disk.")
    
    mapparser.add_argument("--paired", action="store_true", dest="paired",
                           help="Set this flag if the data files represent paired data, i.e. forward and reverse reads")
    
    # filter command
    filterparser = subparsers.add_parser("filter", help="Filter SLAM-seq aligned data")
    filterparser.add_argument("bam", action="store", nargs="+", help="Bam file(s)")
    filterparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir",
                              help="Output directory for mapped BAM files.")
    filterparser.add_argument("-b", "--bed", type=str, required=False, dest="bed",
                              help="BED file, overrides MQ filter to 0")
    filterparser.add_argument("-mq", "--min-mq", type=int, required=False, default=2, dest="mq",
                              help="Minimum mapping quality (default: %(default)d)")
    filterparser.add_argument("-mi", "--min-identity", type=float, required=False, default=0.95, dest="identity",
                              help="Minimum alignment identity (default: %(default)s)")
    filterparser.add_argument("-nm", "--max-nm", type=int, required=False, default=-1, dest="nm",
                              help="Maximum NM for alignments (default: %(default)d)")
    filterparser.add_argument("-t", "--threads", type=int, required=False, dest="threads", default=1,
                              help="Thread number (default: %(default)d)")

    filterparser.add_argument("--paired", action="store_true", dest="paired",
                              help="Set this flag if the data files represent paired data, "
                                   "i.e. forward and reverse reads")

    # snp command

    snpparser = subparsers.add_parser("snp", formatter_class=ArgumentDefaultsHelpFormatter,
                                      help="Call SNPs on SLAM-seq aligned data")
    snpparser.add_argument("bam", action="store", nargs="+", help="Bam file(s)")
    snpparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", default=SUPPRESS,
                           help="Output directory for mapped BAM files.")
    snpparser.add_argument("-r", "--reference", required=True, dest="fasta", type=str, default=SUPPRESS,
                           help="Reference fasta file")
    snpparser.add_argument("-c", "--min-coverage", required=False, dest="cov", type=int, default=10,
                           help="Minimimum coverage to call variant")
    # snpparser.add_argument("-q", "--min-base-qual", type=int, default=13, required=False, dest="minQual",
    # help="Min base quality for T -> C conversions (default: %(default)d)")
    snpparser.add_argument("-f", "--var-fraction", required=False, dest="var", type=float, default=0.8,
                           help="Minimimum variant fraction to call variant")
    snpparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads", help="Thread number")

    # count command

    countparser = subparsers.add_parser("count", help="Count T/C conversions in SLAM-seq aligned data")
    countparser.add_argument("bam", action="store", nargs="+", help="Bam file(s)")
    countparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir", default=SUPPRESS,
                             help="Output directory for mapped BAM files.")
    countparser.add_argument("-s", "--snp-directory", type=str, required=False, dest="snpDir", default=SUPPRESS,
                             help="Directory containing SNP files.")
    countparser.add_argument("-v", "--vcf", type=str, required=False, dest="vcfFile", default=SUPPRESS,
                             help="Externally provided custom variant file.")
    countparser.add_argument("-r", "--reference", type=str, required=True, dest="ref", default=SUPPRESS,
                             help="Reference fasta file")
    countparser.add_argument("-b", "--bed", type=str, required=False, dest="bed", default=SUPPRESS,
                             help="BED file")
    countparser.add_argument("-c", "--conversion-threshold", type=int, dest="conversionThreshold", required=False,
                             default=1,
                             help="Number of T>C conversions required to count read as T>C read (default: %(default)d)")
    countparser.add_argument("-l", "--max-read-length", type=int, required=False, dest="maxLength",
                             help="Max read length in BAM file")
    countparser.add_argument("-q", "--min-base-qual", type=int, default=27, required=False, dest="minQual",
                             help="Min base quality for T -> C conversions (default: %(default)d)")
    countparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads",
                             help="Thread number (default: %(default)d)")

    # all command
    allparser = subparsers.add_parser("all", help="Run entire SLAMdunk analysis")
    allparser.add_argument("files", action="store", nargs="+",
                           help="Single csv/tsv file (recommended) containing all sample files and sample info or a"
                                " list of all sample BAM/FASTA(gz)/FASTQ(gz) files")
    allparser.add_argument("-r", "--reference", type=str, required=True, dest="referenceFile",
                           help="Reference fasta file")
    allparser.add_argument("-b", "--bed", type=str, required=False, dest="bed",
                           help="BED file with 3'UTR coordinates")
    allparser.add_argument("-fb", "--filterbed", type=str, required=False, dest="filterbed",
                           help="BED file with 3'UTR coordinates to filter multimappers (activates -m)")
    allparser.add_argument("-v", "--vcf", type=str, required=False, dest="vcfFile", default=SUPPRESS,
                           help="Skip SNP step and provide custom variant file.")
    allparser.add_argument("-o", "--outputDir", type=str, required=True, dest="outputDir",
                           help="Output directory for slamdunk run.")
    allparser.add_argument("-5", "--trim-5p", type=int, required=False, dest="trim5", default=12,
                           help="Number of bp removed from 5' end of all reads (default: %(default)s)")
    allparser.add_argument("-a", "--max-polya", type=int, required=False, dest="maxPolyA", default=4,
                           help="Max number of As at the 3' end of a read (default: %(default)s)")
    allparser.add_argument("-n", "--topn", type=int, required=False, dest="topn", default=1,
                           help="Max number of alignments to report per read (default: %(default)s)")
    allparser.add_argument("-t", "--threads", type=int, required=False, default=1, dest="threads",
                           help="Thread number (default: %(default)s)")
    allparser.add_argument("-q", "--quantseq", dest="quantseq", action="store_true", required=False,
                           help="Run plain Quantseq alignment without SLAM-seq scoring")
    allparser.add_argument("-e", "--endtoend", action="store_true", dest="endtoend",
                           help="Use a end to end alignment algorithm for mapping.")
    allparser.add_argument("-m", "--multimap", action="store_true", dest="multimap",
                           help="Use reference to resolve multimappers (requires -n > 1).")
    allparser.add_argument("-mq", "--min-mq", type=int, required=False, default=2, dest="mq",
                           help="Minimum mapping quality (default: %(default)s)")
    allparser.add_argument("-mi", "--min-identity", type=float, required=False, default=0.95, dest="identity",
                           help="Minimum alignment identity (default: %(default)s)")
    allparser.add_argument("-nm", "--max-nm", type=int, required=False, default=-1, dest="nm",
                           help="Maximum NM for alignments (default: %(default)s)")
    allparser.add_argument("-mc", "--min-coverage", required=False, dest="cov", type=int, default=10,
                           help="Minimimum coverage to call variant (default: %(default)s)")
    allparser.add_argument("-mv", "--var-fraction", required=False, dest="var", type=float, default=0.8,
                           help="Minimimum variant fraction to call variant (default: %(default)s)")
    allparser.add_argument("-c", "--conversion-threshold",
                           type=int, dest="conversionThreshold", required=False, default=1,
                           help="Number of T>C conversions required to count read as T>C read (default: %(default)d)")
    allparser.add_argument("-rl", "--max-read-length", type=int, required=False, dest="maxLength",
                           help="Max read length in BAM file")
    allparser.add_argument("-mbq", "--min-base-qual", type=int, default=27, required=False, dest="minQual",
                           help="Min base quality for T -> C conversions (default: %(default)d)")
    allparser.add_argument("-sn", "--sampleName", type=str, dest="sampleName", required=False,
                           help="Use this sample name for all supplied samples")
    allparser.add_argument("-sy", "--sampleType", type=str, dest="sampleType", required=False, default="pulse",
                           help="Use this sample type for all supplied samples")
    allparser.add_argument("-st", "--sampleTime", type=int, dest="sampleTime", required=False, default=0,
                           help="Use this sample time for all supplied samples")
    allparser.add_argument("-i", "--sample-index", type=int, required=False, default=-1, dest="sampleIndex",
                           help="Run analysis only for sample <i>. Use for distributing slamdunk analysis on a cluster"
                                " (index is 1-based).")
    allparser.add_argument("-ss", "--skip-sam", action="store_true", dest="skipSAM",
                           help="Output BAM while mapping. Slower but, uses less hard disk.")

    allparser.add_argument("--paired", action="store_true", dest="paired",
                           help="Set this flag if the data files represent paired data, i.e. forward and reverse reads")
    args = parser.parse_args()

    ########################################################################
    # Routine selection
    ########################################################################

    command = args.command

    if command == "map":
        outputDirectory, n, referenceFile, samples, samplesInfos = processPrepare(args)
        _ = processMap(
            samples,
            args.paired,
            args.sampleName,
            args.sampleType,
            args.sampleTime,
            args.sampleIndex,
            args.trim5,
            args.maxPolyA,
            args.quantseq,
            args.endtoend,
            args.topn,
            args.skipSAM,
            outputDirectory,
            samplesInfos,
            referenceFile,
            n
        )

    elif command == "filter":
        outputDirectory = args.outputDir
        n = args.threads
        bed = args.bed
        bam = args.bam
        mq = args.mq
        identity = args.identity
        nm = args.nm
        paired = args.paired

        createDir(outputDirectory)
        _ = processFilter(bam, mq, identity, nm, bed, paired, outputDirectory, n)

    elif command == "snp":
        bam = args.bam
        outputDirectory = args.outputDir
        fasta = args.fasta
        minCov = args.cov
        minVarFreq = args.var
        minQual = 15
        n = args.threads

        createDir(outputDirectory)
        processSNP(bam, fasta, minCov, minVarFreq, minQual, outputDirectory, n)

    elif command == "count":
        outputDirectory = args.outputDir
        createDir(outputDirectory)
        n = args.threads
        processCount(
            args.bam,
            args.ref,
            args.bed if "bed" in args else None,
            args.maxLength,
            args.minQual,
            args.conversionThreshold,
            args.snpDir if "snpDir" in args else None,
            args.vcfFile if "vcfFile" in args else None,
            outputDirectory,
            n
        )

    elif command == "all":
        runAll(args)
        dunkFinished()

    else:
        parser.error("Too few arguments.")


if __name__ == "__main__":
    run()
