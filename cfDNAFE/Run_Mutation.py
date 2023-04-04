from .Base import Base
from .utils import commonError, maxCore
import os
import math

__metaclass__ = type


class runMutation(Base):
    def __init__(
            self,
            pathToRunMutation=None,
            bamInput=None,
            genome_reference=None,
            outputdir=None,
            refGenome='hg38',
            ID='test',
            baseQ=30,
            mapQ=60,
            threads=None,
            **kwargs
    ):
        """
                This function is used for Converting the bam file into a vcf file for mutation detection.
                runMutation(pathToRunMutation=None, bamInput=None,  genome_reference=None, outputdir=None, refGenome='hg38',ID='test',baseQ=30, mapQ=60, threads=None,)
                {P}arameters:
                    pathToRunMutation: str, /path/to/RunMutation.R.
                    bamInput: list, input .bam files
                    ID:str, group id.
                    refGenome: str, human genome reference 'hg19' or 'hg38'.
                    outputdir: str, output result folder, None means the same folder as input files.
                    genome_reference: str, Human reference genome storage path, input .fa or fa.gz file.
                    baseQ: int, (default 30) skip bases with baseQ/BAQ smaller than INT
                    mapQ: int, (default 60) skip alignments with mapQ smaller than INT
                    threads: int (default None) Whether to use multithreading

                {O}utput:
                    outputfile1: sample.vcf.gz


                """
        super(runMutation, self).__init__()
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
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        vcfPath = os.path.join(self.getOutput('outputdir'), 'vcfpath')
        if not os.path.exists(vcfPath):
            os.mkdir(vcfPath)
        self.setOutput("vcfOutput",
                       [
                           os.path.join(vcfPath, self.getMaxFileNamePrefixV2(x))
                           + '.vcf'
                           for x in self.getInput("bamInput")
                       ])

        self.setParam("baseQ", baseQ)
        self.setParam('mapQ', mapQ)
        self.setParam("ID", ID)
        self.setParam("refGenome", refGenome)
        sample_signature = []
        mutation_profile = []

        for x in self.getParam("ID"):
            outputfile1 = os.path.join(self.getOutput("outputdir"), x) + ".96.mutation.profle"
            outputfile2 = os.path.join(self.getOutput("outputdir"), x) + ".signatures"
            mutation_profile.append(outputfile1)
            sample_signature.append(outputfile2)
        self.setOutput(
            "sampleOutput1", mutation_profile,
        )
        self.setOutput(
            "sampleOutput2", sample_signature,
        )
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
            self.multiRun(args=cmd, func=None, nCore=maxCore(self.getParam("threads")))



        multi_run_len = len(vcfPath)
        cmd = []
        for i in range(multi_run_len):
            cmd.append(
                self.cmdCreate(
                    [
                        'Rscript',
                        pathToRunMutation,
                        "--inputFile",
                        self.getInput("vcfPath")[i],
                        "--refGenome",
                        self.getParam("refGenome"),
                        "--outputfile1",
                        self.getOutput("sampleOutput1")[i],
                        "--outputfile2",
                        self.getOutput("sampleOutput2")[i]
                    ]
                )
            )

        if verbose:
            self.run(cmd)
        else:
            self.multiRun(args=cmd, func=None, nCore=maxCore(self.getParam("threads")))
