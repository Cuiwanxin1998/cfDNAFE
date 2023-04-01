from .Base import Base
from .utils import commonError, maxCore
import os
import math

__metaclass__ = type


class runMutation(Base):
    def __init__(
            self,
            pathToRunMutation=None,
            vcfPath=None,
            outputdir=None,
            refGenome='hg38',
            ID=None,
            threads=None,
            **kwargs
    ):
        """
                       This function is used for Generating mutation signature.
                       runMutation(pathToRunMutation=None, vcfInput=None, outputdir=None, id='test', cutoff=None, refGenome='hg38', threads=None,)
                       {P}arameters:
                           pathToRunMutation: str, /path/to/RunMutation.R.
                           vcfPath: str, Path to .vcf or .vcf.gz files.
                           outputdir: str, output result folder, None means the same folder as input files.
                           ID:list, group id.
                           refGenome: str, human genome reference 'hg19' or 'hg38'.
                           threads: int (default None) Whether to use multithreading.
                       {O}utput:
                           outputfile1: ID.signatures


                       """
        super(runMutation, self).__init__()
        if pathToRunMutation is None:
            raise commonError("Parameter pathToRunMutation must require. /path/to/RunMutation.R.")
        if isinstance(vcfPath, list):
            for path in vcfPath:
                if os.path.exists(path):
                    vcfInput = []
                    files = os.listdir(path)
                    for file in files:
                        if file.endswith('.vcf'):
                            vcfInput.append(os.path.join(path, file))
                        if file.endswith('.vcf.gz'):
                            vcfInput.append(os.path.join(path, file))
                    if len(vcfInput) == 0:
                        raise commonError(
                            "Check that the vcf path " + path + "you entered contains files ending in.vcf or.vcf.gz")
                else:
                    raise commonError("this path" + path + 'not exist')
            self.setInput("vcfPath", vcfPath)
        else:
            if os.path.exists(vcfPath):
                vcfInput = []
                files = os.listdir(vcfPath)
                for file in files:
                    if file.endswith('.vcf'):
                        vcfInput.append(os.path.join(vcfPath, file))
                    if file.endswith('.vcf.gz'):
                        vcfInput.append(os.path.join(vcfPath, file))
                if len(vcfInput) == 0:
                    raise commonError(
                        "Check that the vcf path " + vcfPath + "you entered contains files ending in.vcf or.vcf.gz")
            else:
                raise commonError("this path" + vcfPath + 'not exist')
            self.setInput("vcfPath", [vcfPath])
        if ID is None:
            ID = []
            for i in range(len(self.getInput("vcfPath"))):
                ID.append("run" + str(i))
            self.setInput("ID", ID)
        else:
            if isinstance(ID, list):
                if len(ID) != len(self.getInput("vcfPath")):
                    raise commonError("Number of input ID and vcfPath not equal.")
            else:
                if len(self.getInput("vcfPath")) != 1:
                    raise commonError("Number of vcfPath not equal 1, you must input a list ID")
            self.setInput("ID", ID)

        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True

        if outputdir is None:
            self.setOutput(
                "outputdir",
                os.path.dirname(os.path.abspath(self.getInput("vcfPath")[0])),
            )
        else:
            self.setOutput("outputdir", outputdir)

        self.setParam("refGenome", refGenome)
        sample_signature = []
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        for x in self.getInput("ID"):
            outputfile = os.path.join(self.getOutput("outputdir"), x) + ".signatures"
            sample_signature.append(outputfile)
        self.setOutput(
            "sampleOutput", sample_signature,
        )
        multi_run_len = len(self.getInput("vcfPath"))
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
                        "--outputfile",
                        self.getOutput("sampleOutput")[i]
                    ]
                )
            )

        if verbose:
            self.run(cmd)
        else:
            self.multiRun(args=cmd, func=None, nCore=maxCore(math.ceil(self.getParam("threads") / 4)))
