from .Base import Base
from .utils import commonError, calc_FSR, maxCore
import os
import math

__metaclass__ = type


class runFSR(Base):
    def __init__(
        self,
        bedgzInput=None,
        binInput=None,
        windows=100000,
        outputdir=None,
        threads=None,
        **kwargs
    ):
        """
        This function is used for computing FSC.
        It was generated using the coverages of short (65-150bp), intermediate (151-260bp), long (261-400bp), and total (65-400bp) cfDNA fragments.
        runFSR(bedgzInput=None, binInput=None, outputdir=None, outputdir=None, threads=None,)
        {P}arameters:
            bedgzInput: list, input bed.gz files.
            binInput: str, regions of chromosome N kb bin file.
            windows:int, the window size.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, Whether to use multithreading
        {O}utput:
        outputfile1: sample.FSC.txt chr      start           end                     ratio(short,intermediate,long)
                                   chr1        0	        4999999 	    0.17272139343023732,0.7420107578411599,0.07155221677046829
                                   chr1	    5000000	        9999999	        0.17621836169747568,0.7299949941017261,0.07902522243317145
                                   chr1	    10000000	    14999999	    0.1732550594006567,0.7378802361663213,0.07380528089811797
                                   chr1	    15000000	    19999999	    0.1767141959043928,0.7277467257505905,0.08046175989872899
                                                                ......

        """
        super(runFSR, self).__init__()
        # set bedgz input
        if bedgzInput is None:
            raise commonError("Parameter bedgzInput must require.")
        else:
            self.setInput("bedgzInput", bedgzInput)

        # set bin input
        if binInput is None:
            raise commonError("binInput must be not None!")
        else:
            self.setParam("binInput", binInput)


        # set outputdir
        if outputdir is None:
            self.setOutput(
                "outputdir", os.path.dirname(os.path.abspath(self.getInput("bedgzInput")[0])),
            )
        else:
            self.setOutput("outputdir", outputdir)

        # set threads
        if threads is not None:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", 1)
            verbose = True
        self.setParam("windows", windows)
        sample_FSC = []
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        for x in self.getInput("bedgzInput"):
            outputfile = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".FSR.txt"

            sample_FSC.append(outputfile)
        self.setOutput(
            "sampleOutput", sample_FSC,
        )
        multi_run_len = len(self.getInput("bedgzInput"))
        if verbose:
            for i in range(multi_run_len):
                calc_FSR(
                    bedgzInput=self.getInput("bedgzInput")[i],
                    binInput=self.getParam("binInput"),
                    windows=self.getParam("windows"),
                    outputfile=self.getOutput("sampleOutput")[i],

                )
        else:
            args = [
                [
                    self.getInput("bedgzInput")[i],
                    self.getParam("binInput"),
                    self.getParam("windows"),
                    self.getOutput("sampleOutput")[i],

                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=calc_FSR, nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
            )
