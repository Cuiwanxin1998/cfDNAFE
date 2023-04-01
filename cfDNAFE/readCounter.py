from .Base import Base
from .utils import commonError, maxCore
import os
import math

__metaclass__ = type


class readCounter(Base):
    def __init__(
            self,
            bamInput=None,
            outputdir=None,
            window_size=1000000,
            quality=20,
            chromosome=None,
            threads=None,
            **kwargs
    ):
        """
                This function is used for Generating read count coverage information using readCounter from the HMMcopy suite.
                readCounter(pathToreadCounter=None, bamInput=None, outputdir=None, window_size=1000000, quality=20, chromosome=None, threads=None,)
                {P}arameters:
                    bamInput: list, input .bam files
                    outputdir: str, output result folder, None means the same folder as input files.
                    window_size: int, (default 1000000). input window size.
                    quality: int, (default 20). remove reads with lower quality
                    chromosome: str, Chromosomes to be processed(default 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22').
                    threads: int (default None) Whether to use multithreading
                {O}utput:
                    outputfile1: sample.wig  fixedStep chrom=chr1 start=1 step=1000000 span=1000000
                                             5051
                                             21833
                                             17901
                                             28555
                                             25501
                                             24853
                                             24833
                                             24551
                                             21182


                """
        super(readCounter, self).__init__()

        if bamInput is None:
            raise commonError("Parameter bamInput must require. the input bam must be sorted")
        else:
            self.setInput("bamInput", bamInput)

        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True

        if outputdir is None:
            self.setOutput(
                "outputdir",
                os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
            )
        else:
            self.setOutput("outputdir", outputdir)

        self.setOutput("wigOutput",
                       [
                           os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x))
                           + '.wig'
                           for x in self.getInput("bamInput")
                       ])
        if chromosome is None:
            self.setParam('chromosome',
                          'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22')
        else:
            self.setParam('chromosome', chromosome)
        self.setParam("window", window_size)
        self.setParam('quality', quality)
        multi_run_len = len(self.getInput("bamInput"))
        cmd = []
        for i in range(multi_run_len):
            cmd.append(
                self.cmdCreate(
                    [
                        'readCounter',
                        "--window",
                        self.getParam("window"),
                        "--quality",
                        self.getParam("quality"),
                        "--chromosome",
                        self.getParam('chromosome'),
                        self.getInput("bamInput")[i],
                        '>' + self.getOutput("wigOutput")[i]
                    ]
                )
            )

        if verbose:
            self.run(cmd)
        else:
            self.multiRun(args=cmd, func=None, nCore=maxCore(math.ceil(self.getParam("threads") / 4)))
