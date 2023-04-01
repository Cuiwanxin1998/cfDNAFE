from .Base import Base
from .utils import commonError, maxCore
import os
import math

__metaclass__ = type


class BamToVcf(Base):
    def __init__(
            self,
            bamInput=None,
            genome_reference=None,
            outputdir=None,
            baseQ=30,
            mapQ=60,
            threads=None,
            **kwargs
    ):
        """
                This function is used for Converting the bam file into a vcf file for mutation detection.
                readCounter(pathToreadCounter=None, bamInput=None, outputdir=None, window_size=1000000, quality=20, chromosome=None, threads=None,)
                {P}arameters:
                    bamInput: list, input .bam files
                    outputdir: str, output result folder, None means the same folder as input files.
                    genome_reference: str, Human reference genome storage path, input .fa or fa.gz file.
                    baseQ: int, (default 30) skip bases with baseQ/BAQ smaller than INT
                    mapQ: int, (default 60) skip alignments with mapQ smaller than INT
                    threads: int (default None) Whether to use multithreading
                {O}utput:
                    outputfile1: sample.vcf.gz


                """
        super(BamToVcf, self).__init__()
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

        # set genome_reference input
        if genome_reference is None:
            raise commonError("Parameter genome_reference must be not None!")
        else:
            self.setParam("genome_reference", genome_reference)

        if outputdir is None:
            self.setOutput(
                "outputdir",
                os.path.dirname(os.path.abspath(self.getInput("bamInput")[0])),
            )
        else:
            self.setOutput("outputdir", outputdir)

        self.setOutput("vcfOutput",
                       [
                           os.path.join(self.getOutput('outputdir'), self.getMaxFileNamePrefixV2(x))
                           + '.vcf'
                           for x in self.getInput("bamInput")
                       ])

        self.setParam("baseQ", baseQ)
        self.setParam('mapQ', mapQ)
        multi_run_len = len(self.getInput("bamInput"))
        cmd = []
        for i in range(multi_run_len):
            cmd.append(
                self.cmdCreate(
                    [
                        'bcftools mpileup',
                        self.getInput("bamInput")[i],
                        "-f",
                        self.getParam("genome_references"),
                        "-q",
                        self.getParam("mapQ"),
                        "-Q",
                        self.getParam('baseQ'),
                        "|",
                        "bcftools call",
                        "-mv"
                        '-o',
                        self.getOutput("vcfOutput")[i],


                    ]
                )
            )

        if verbose:
            self.run(cmd)
        else:
            self.multiRun(args=cmd, func=None, nCore=maxCore(math.ceil(self.getParam("threads") / 4)))
