from .Base import Base
from .utils import commonError, maxCore
import os
import math

__metaclass__ = type


class runMutation(Base):
    def __init__(
            self,
            pathToRunMutation=None,
            vcfInput=None,
            outputdir=None,
            id=['test'],
            cutoff=0.85,
            refGenome='hg38',
            threads=None,
            **kwargs
    ):
        """
                       This function is used for Generating mutation signature.
                       runMutation(pathToRunMutation=None, vcfInput=None, outputdir=None, id='test', cutoff=None, refGenome='hg38', threads=None,)
                       {P}arameters:
                           pathToRunMutation: str, /path/to/RunMutation.R.
                           vcfInput: str, Path to .vcf files.
                           outputdir: str, output result folder, None means the same folder as input files.
                           id: list, Files ID
                           cutoff: int, a cosine similarity of more than cutoff with an existing COSMIC signature
                           refGenome: str, human genome reference 'hg19' or 'hg38'
                           threads: int (default None) Whether to use multithreading
                       {O}utput:
                           outputfile1: ID.mutation_signature.txt


                       """
        super(runMutation, self).__init__()
        if pathToRunMutation is None:
            raise commonError("Parameter pathToRunMutation must require. /path/to/RunMutation.R.")
        if vcfInput is None:
            raise commonError("Parameter vcfInput must require. /path/to/.vcf")
        else:
            self.setInput("vcfInput", vcfInput)

        if id is None:
            raise commonError("Parameter id must require.")
        else:
            self.setInput("id", id)


        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True

        if outputdir is None:
            self.setOutput(
                "outputdir",
                os.path.dirname(os.path.abspath(self.getParam("vcfInput"))),
            )
        else:
            self.setOutput("outputdir", outputdir)


        self.getParam("refGenome", refGenome)
        self.getParam("cutoff", cutoff)

        multi_run_len = len(self.getInput("vcfInput"))
        cmd = []
        for i in range(multi_run_len):
            cmd.append(
                self.cmdCreate(
                    [
                        'Rscript',
                        pathToRunMutation,
                        "--inputFile",
                        self.getInput("vcfInput")[i],
                        "--refGenome",
                        self.getParam("refGenome"),
                        "--cutoff",
                        self.getParam('cutoff'),
                        "--id",
                        self.getInput("id")[i],
                        "--outputdir",
                        self.getOutput('outputdir')
                    ]
                )
            )

        if verbose:
            self.run(cmd)
        else:
            self.multiRun(args=cmd, func=None, nCore=maxCore(math.ceil(self.getParam("threads") / 4)))
