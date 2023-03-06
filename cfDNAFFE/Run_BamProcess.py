from .Base import Base
from .utils import commonError, bam_process, maxCore
import os
import math

__metaclass__ = type


class runBamProcess(Base):
    def __init__(
        self,
        bamInput=None,
        blacklistInput=None,
        outputdir=None,
        genome_reference=None,
        CHR=None,
        mapQuality=30,
        fragFilter=False,
        minLen=None,
        maxLen=None,
        k_mer=6,
        threads=None,
        **kwargs
    ):
        """
        This function is used to process the sorted bam file.
        runBamProcess(bamInput=None, blacklistInput=None, outputdir=None, genome_reference=None, CHR=None, mapQuality=30, fragFilter=False, minLen=None, maxLen=None, k_mer=6, threads=None,)
        {P}arameters:
            bamInput: list, input .bam files.
            blacklistInput: str, regions of blacklist file. download from this site https://github.com/Boyle-Lab/Blacklist/.
            outputdir: str, output result folder, None means the same folder as input files.
            genome_reference: str, input .fa or fa.gz file.
            CHR: list, Chromosomes to be processed(default ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr18', 'chr19', 'chr20', 'chr21',
               'chr22', 'chrX']).
            empty: bool, keep files of empty blocks (default False).
            fragFilter: bool, (default False). Whether filter fragment by length, only for paired data.
            minLen: int, Min fragment length.
            maxLen: int, Max fragment length.
            k_mer: int, (default 6) The number of motifs to analyze.
            threads: int (default None) Whether to use multithreading
        {O}utput:
            outputfile1: sample.bed.gz chr      start      end      gc_content
                                       chr1     13160     13334     0.597701149
                                       chr1     13294     13461     0.550898203
                                       chr1     14555     14730     0.611428571
                                       chr1     15072     15255     0.644808746
                                                      ......
            outputflie2: Endmotif.txt  Endmotif     frequency
                                       AAAAAA       0.00141453
                                       AAAAAC       0.00132175
                                       AAAAAT       0.00130899
                                       AAAAAG       0.00098876
                                             ......
            outputflie3: Breakpointmotif.txt  Breakpointmotif     frequency
                                                   AAAAAA       0.00108359
                                                   AAAAAC       0.00041196
                                                   AAAAAT       0.00054755
                                                   AAAAAG       0.00032743
                                                          ......
            outputfile4: MDS.txt    outputfile2_path    MDS
                                    path               0.964121

        """


        super(runBamProcess, self).__init__()
        # set bam input
        if bamInput is None:
            raise commonError("Parameter bamInput must require.")
        else:
            self.setInput("bamInput", bamInput)

        # set blacklist input
        if blacklistInput is None:
            raise commonError("Parameter blacklistInput must be not None!")
        else:
            self.setParam("blacklistInput", blacklistInput)


        # set outputdir
        if outputdir is None:
            self.setOutput(
                "outputdir", os.path.dirname(os.path.abspath(self.getInput("bedgzInput")[0])),
            )
        else:
            self.setOutput("outputdir", outputdir)

        # set genome_reference input
        if genome_reference is None:
            raise commonError("Parameter genome_reference must be not None!")
        else:
            self.setParam("genome_reference", genome_reference)

        #set chromosome input
        if CHR is None:
            self.setParam("CHR", ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr18', 'chr19', 'chr20', 'chr21',
               'chr22', 'chrX'])
        else:
            self.setParam("CHR", CHR)
        self.setParam("mapQuality", mapQuality)
        self.setParam("fragFilter", fragFilter)
        self.setParam("minLen", minLen)
        self.setParam("maxLen", maxLen)
        self.setParam("k_mer", k_mer)

        # set threads
        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True
        sample_Bed = []
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        for x in self.getInput("bamInput"):
            outputfile = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".bed"

            sample_Bed.append(outputfile)
        self.setOutput(
            "sampleOutput", sample_Bed,
        )
        multi_run_len = len(self.getInput("bamInput"))
        if verbose:
            for i in range(multi_run_len):
                bam_process(

                    bamInput=self.getInput("bamInput")[i],
                    blacklistInput=self.getParam("blacklistInput"),
                    bedOutput=self.getOutput("sampleOutput")[i],
                    genome_reference=self.getParam("genome_reference"),
                    CHR=self.getParam("CHR"),
                    mapQuality=self.getParam("mapQuality"),
                    fragFilter=self.getParam("fragFilter"),
                    minLen=self.getParam("minLen"),
                    maxLen=self.getParam("maxLen"),
                    k_mer=self.getParam("k_mer"),

                )
        else:
            args = [
                [
                    self.getInput("bamInput")[i],
                    self.getParam("blacklistInput"),
                    self.getOutput("sampleOutput")[i],
                    self.getParam("genome_reference"),
                    self.getParam("CHR"),
                    self.getParam("mapQuality"),
                    self.getParam("fragFilter"),
                    self.getParam("minLen"),
                    self.getParam("maxLen"),
                    self.getParam("k_mer"),

                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=bam_process, nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
            )
