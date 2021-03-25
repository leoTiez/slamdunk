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

from __future__ import print_function
import os, re

from slamdunk.utils.misc import files_exist, checkStep, run, pysamIndex, removeFile, replaceExtension, shellerr  # @UnresolvedImport
from slamdunk.version import __ngm_version__  # @UnresolvedImport

projectPath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sort(inputSAM, outputBAM, log, threads=1, keepSam=True, dry=False, verbose=True, isPaired=False):
    if files_exist(inputSAM) and checkStep([inputSAM], [outputBAM + ".bai"]):
        runSam2bam(
            inputSAM,
            outputBAM,
            log,
            index=False,
            sort="name" if isPaired else None,
            delinFile=not keepSam,
            onlyProperPaired=True,
            threads=threads,
            dry=dry,
            verbose=verbose
        )
    else:
        print("Skipped sorting for " + inputSAM, file=log)


def checkNextGenMapVersion():
    ngmHelp = shellerr("ngm", raiseError=False)
    matchObj = re.match(r'.*([0-9]+\.[0-9]+\.[0-9]+).*', str(ngmHelp), re.M|re.I)
    if matchObj:
        version = matchObj.group(1)
        if version != __ngm_version__:
            raise RuntimeError("NextGenMap version expected: %s  but found %s. Please reinstall slamdunk package." % (
                __ngm_version__, version))
    else:
        raise RuntimeError("Could not get NextGenMap version. Please reinstall slamdunk package.")


def runSam2bam(
        inFile,
        outFile,
        log,
        index=True,
        sort=None,
        delinFile=False,
        onlyUnique=False,
        onlyProperPaired=False,
        filterMQ=0,
        L=None,
        threads=1,
        verbose=False,
        dry=False
):
    if delinFile and files_exist(outFile) and not files_exist(inFile):
        print("Skipping sam2bam for %s" % outFile, file=log)
    else:
        if onlyUnique and filterMQ == 0:
            filterMQ = 1

        success = True
        cmd = ["samtools view", "-@", str(threads), "-Sb", "-o", outFile, inFile]
        if filterMQ > 0:
            cmd += ["-q", str(filterMQ)]
        if onlyProperPaired:
            cmd += ["-f", "2"]
        if L is not None:
            cmd += ["-L", L]
        run(" ".join(cmd), log, verbose=verbose, dry=dry)

        if sort is not None:
            tmp = outFile + "_tmp"
            if not dry:
                os.rename(outFile, tmp)
            if sort.lower() == "index":
                run(" ".join(["samtools sort", "-@", str(threads), "-o",  outFile, tmp]), log, verbose=verbose, dry=dry)
            elif sort.lower() == "name":
                run(
                    " ".join(["samtools sort -n", "-@", str(threads), "-o", outFile, tmp]),
                    log,
                    verbose=verbose,
                    dry=dry
                )
            if success:
                removeFile(tmp)
        if success and delinFile:
            if not dry:
                removeFile(inFile)

    if index:
        pysamIndex(outFile)


def Map(
        inputBAM,
        inputReference,
        outputSAM,
        log,
        quantseqMapping,
        endtoendMapping,
        threads=1,
        parameter="--no-progress --slam-seq 2",
        outputSuffix="_ngm_slamdunk",
        trim5p=0,
        maxPolyA=-1,
        topn=1,
        sampleId=None,
        sampleName="NA",
        sampleType="NA",
        sampleTime=0,
        printOnly=False,
        verbose=True,
        force=False,
        isPaired=False
):
    if quantseqMapping:
        parameter = "--no-progress"

    if trim5p > 0:
        parameter = parameter + " -5 " + str(trim5p)

    if maxPolyA > -1:
        parameter = parameter + " --max-polya " + str(maxPolyA)

    if endtoendMapping:
        parameter = parameter + " -e "
    else:
        parameter = parameter + " -l "

    if sampleId is not None:
        parameter = parameter + " --rg-id " + str(sampleId)
        if sampleName != "":
            parameter = parameter + " --rg-sm " + sampleName + ":" + sampleType + ":" + str(sampleTime)

    if topn > 1:
        parameter = parameter + " -n " + str(topn) + " --strata "

    files = [inputReference]
    files.append(inputBAM) if not isPaired else files.extend(inputBAM)
    files = [os.path.expanduser(p) for p in files]
    if checkStep(files, [replaceExtension(outputSAM, ".bam")], force):
        cmd = "ngm %s -r %s %s -t %s %s -o %s" % (
            "" if outputSAM.endswith(".sam") else "-b",
            files[0],
            "-q %s" % files[1] if not isPaired else "-1 %s -2 %s" % (files[1], files[2]),
            threads,
            parameter,
            outputSAM
        )
        run(cmd, log, verbose=verbose, dry=printOnly)
    else:
        print("Skipped mapping for " + inputBAM if not isPaired else inputBAM[0], file=log)
