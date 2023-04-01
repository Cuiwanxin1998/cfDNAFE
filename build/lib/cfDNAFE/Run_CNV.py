import pandas as pd
import shutil
from .Base import Base
from .utils import commonError, maxCore
import os
import math

__metaclass__ = type


class runCNV(Base):
    def __init__(
            self,
            pathTorunIchorCNA=None,
            wig=None,
            ploidy="'c(2)'",
            normal="'c(0.95, 0.99, 0.995, 0.999)'",
            maxCN=3,
            gcWig=None,
            ID=None,
            mapWig=None,
            centromere=None,
            normalPanel=None,
            includeHOMD=False,
            chrs="'c(1:22)'",
            chrTrain="'c(1:22)'",
            estimateNormal=True,
            estimatePloidy=True,
            estimateScPrevalence=False,
            scStates="'c()'",
            txnE=0.9999999,
            txnStrength=1000000,
            outputdir=None,
            threads=None,
            **kwargs
    ):
        """
                This function is used for analysising Copy number.
                runCNV(pathToreadCounter=None, bamInput=None, outputdir=None, window_size=1000000, quality=20, chromosome=None, threads=None,)
                {P}arameters:
                    pathTorunIchorCNA: str, /path/to/ichorCNA/scripts/runIchorCNA.R.
                    wig: list, Path to tumor WIG file
                    ploidy: str, Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired.
                    normal: str, Initial normal contamination; can be more than one value if additional normal initializations are desired.
                    maxCN: int, (default 3). Total clonal CN states.
                    gcWig: str, Path to GC-content WIG file.
                    ID: list, Patient ID.
                    mapWig: str, Path to mappability score WIG file.
                    centromere:str, File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package.
                    normalPanel: str, Median corrected depth from panel of normals.
                    includeHOMD: bool, If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb).
                    chrs: str, Specify chromosomes to analyze.
                    chrTrain: str, Specify chromosomes to estimate params.
                    estimateNormal: bool, Estimate normal.
                    estimatePloidy: bool, Estimate tumour ploidy.
                    estimateScPrevalence, bool, Estimate subclonal prevalence.
                    scStates, str, Subclonal states to consider.
                    txnE, int, Self-transition probability. Increase to decrease number of segments
                    txnStrength, int, Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE.
                    outputdir: str, output result folder, None means the same folder as input files.
                    threads: int (default None) Whether to use multithreading.
                {O}utput:
                    outputfile1: sample.cna.seg


                """
        super(runCNV, self).__init__()
        if pathTorunIchorCNA is None:
            commonError("Path to R file runIchorCNA; Required!")

        if len(wig) != len(ID):
            commonError(
                "the length of wig file does not equal the length of patient ID! Please Check!")
        if wig is None:
            commonError(
                "Wig Path to sample WIG file. Required! Generating read count coverage information using readCounter()")
        else:
            self.setInput("wigInput", wig)

        if ID is None:
            commonError(
                "patient ID. Required! ")
        else:
            self.setParam("ID", ID)
        if gcWig is None:
            commonError("gcWig Path to GC-content WIG file; Required!")

        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True


        if outputdir is None:
            self.setOutput(
                "outputdir",
                [os.path.dirname(os.path.abspath(self.getInput("wigInput")[0]))],
            )
        else:
            self.setOutput("outputdir", outputdir)
        self.setParam('ploidy', ploidy)
        self.setParam('normal', normal)
        self.setParam('maxCN', maxCN)
        self.setParam('gcWig', gcWig)

        self.setParam('mapWig', mapWig)
        self.setParam('centromere', centromere)
        self.setParam('normalPanel', normalPanel)
        self.setParam('includeHOMD', includeHOMD)
        self.setParam('chrs', chrs)
        self.setParam('chrTrain', chrTrain)
        self.setParam('estimateNormal', estimateNormal)
        self.setParam('estimatePloidy', estimatePloidy)
        self.setParam('estimateScPrevalence', estimateScPrevalence)
        self.setParam('scStates', scStates)
        self.setParam('txnE', txnE)
        self.setParam('txnStrength', txnStrength)
        self.setParam('threads', threads)

        dirs = []
        for x in self.getParam("ID"):
            newdir = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x))
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            dirs.append(newdir)
        self.setOutput(
            "sampleOutputdir", dirs,
        )
        multi_run_len = len(self.getInput("wigInput"))
        cmd = []

        for i in range(multi_run_len):
            cmd.append(
                self.cmdCreate(
                    [
                        'Rscript',
                        pathTorunIchorCNA,
                        "--id",
                        self.getParam("ID")[i],
                        "--WIG",
                        self.getInput('wigInput')[i],
                        "--ploidy",
                        self.getParam('ploidy'),
                        "--normal",
                        self.getParam('normal'),
                        "--maxCN",
                        self.getParam('maxCN'),
                        "--gcWig",
                        self.getParam('gcWig'),
                        "--mapWig",
                        'NULL' if self.getParam('mapWig') is None else self.getParam('mapWig'),
                        "--centromere",
                        'NULL' if self.getParam('centromere') is None else self.getParam('centromere'),
                        "--normalPanel ",
                        'NULL' if self.getParam('normalPanel') is None else self.getParam('normalPanel'),
                        "--includeHOMD",
                        self.getParam('includeHOMD'),
                        "--chrs",
                        self.getParam('chrs'),
                        "--chrTrain",
                        self.getParam('chrTrain'),
                        "--estimateNormal",
                        self.getParam('estimateNormal'),
                        "--estimatePloidy",
                        self.getParam('estimatePloidy'),
                        "--estimateScPrevalence",
                        self.getParam('estimateScPrevalence'),
                        "--scStates",
                        self.getParam('scStates'),
                        "--txnE",
                        self.getParam('txnE'),
                        "--txnStrength",
                        self.getParam('txnStrength'),
                        "--outDir",
                        self.getOutput('sampleOutputdir')[i]
                    ]
                )
            )

        if verbose:
            self.run(cmd)
        else:
            self.multiRun(args=cmd, func=None, nCore=maxCore(self.getParam("threads")))

        for i in range(multi_run_len):
            outputdir = self.getOutput('outputdir')
            path = self.getOutput('sampleOutputdir')[i]
            ID = self.getParam('ID')[i]
            cnvfile = os.path.join(path, ID + '.cna.seg')
            data = pd.read_table(cnvfile, usecols=[0, 1, 2, 5])
            outputfile = os.path.join(outputdir, ID+'.CNV')
            data.to_csv(outputfile, sep='\t', na_rep='nan')
            shutil.rmtree(path)
