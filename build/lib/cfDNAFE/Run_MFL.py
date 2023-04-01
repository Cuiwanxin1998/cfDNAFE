from .Base import Base
from .utils import commonError, calc_MFL, maxCore
import os
import math

__metaclass__ = type


class runMFL(Base):
    def __init__(
        self,
        bamInput=None,
        markInput=None,
        outputdir=None,
        mapQuality=30,
        minCpG=3,
        methyThreshold=0.75,
        unmethyThreshold=0.25,
        threads=None,
        **kwargs
    ):
        """
        This function is used to process the sorted bam file.
        runBamProcess(bamInput=None, blacklistInput=None, outputdir=None, genome_reference=None, CHR=None, mapQuality=30, fragFilter=False, minLen=None, maxLen=None, k_mer=6, threads=None,)
        {P}arameters:
            bamInput: list, input .bam files.
            markInput: str, regions of mark file. including chrom start end.
            outputdir: str, output result folder, None means the same folder as input files.
            mapQuality: int, Min fragment quality.
            minCpG: int, Min fragment CpG number.
            methyThreshold: float, thresholds of more than or equal to 75% methylated CpGs for most methylated reads.
            unmethyThreshold: float, thresholds of less than or equal to 25% methylated CpGs for most unmethylated reads.
            threads: int (default None) Whether to use multithreading
        {O}utput:
            outputfile : sample.bed.gz         chr                Mfragment       Xfragment       Ufragment
                                       chr17:77127220-77127487     0.3333          0.3333          0.33333
                                                      ......

        """


        super(runMFL, self).__init__()
        # set bam input
        if bamInput is None:
            raise commonError("Parameter bamInput must require.")
        else:
            self.setInput("bamInput", bamInput)



        # set outputdir
        if outputdir is None:
            self.setOutput(
                "outputdir", os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
            )
        else:
            self.setOutput("outputdir", outputdir)

        self.setParam("markInput",markInput)
        self.setParam("mapQuality", mapQuality)
        self.setParam("minCpG", minCpG)
        self.setParam("methyThreshold", methyThreshold)
        self.setParam("unmethyThreshold", unmethyThreshold)

        # set threads
        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True
        outputFile = []
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        for x in self.getInput("bamInput"):

            file = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".tsv"

            outputFile.append(file)
        self.setOutput(
            "outputfile", outputFile,
        )
        multi_run_len = len(self.getInput("bamInput"))
        if verbose:
            for i in range(multi_run_len):
                calc_MFL(

                    bamInput=self.getInput("bamInput")[i],
                    markInput=self.getParam("markInput"),
                    outputFile=self.getOutput("outputfile")[i],
                    mapQuality=self.getParam("mapQuality"),
                    minCpG=self.getParam('minCpG'),
                    methyThreshold=self.getParam('methyThreshold'),
                    unmethyThreshold=self.getParam('unmethyThreshold'),

                )
        else:
            args = [
                [
                    self.getInput("bamInput")[i],
                    self.getParam("markInput"),
                    self.getOutput("outputfile")[i],
                    self.getParam("mapQuality"),
                    self.getParam('minCpG'),
                    self.getParam('methyThreshold'),
                    self.getParam('unmethyThreshold'),
                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=calc_MFL, nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
            )
