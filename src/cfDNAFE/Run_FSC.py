from .Base import Base
from .utils import commonError, calc_FSC, maxCore
import os
import math

__metaclass__ = type


class runFSC(Base):
    def __init__(
        self,
        bedgzInput=None,
        binInput=None,
        windows=100000,
        outputdir=None,
        continue_N=50,
        threads=1,
        **kwargs
    ):
        """
        This function is used for computing FSC.
        It was generated using the coverages of short (65-150bp), intermediate (151-260bp), long (261-400bp), and total (65-400bp) cfDNA fragments.
        runFSC(bedgzInput=None, binInput=None, outputdir=None, outputdir=None, threads=None,)
        {P}arameters:
            bedgzInput: list, input bed.gz files.
            binInput: str, regions of chromosome N kb bin file.
            windows:int, the window size.
            continue_N:int, continues window number.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, Whether to use multithreading
        {O}utput:
        outputfile1: sample.FSC.txt chr      start           end                     z-score(short,intermediate,long,total)
                                   chr1        0	        4999999 	    -1.0    -1.2    -0.5    -0.7    -0.8    .....
                                                                ......

        """
        super(runFSC, self).__init__()
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
        if threads != 1:
            self.setParam("threads", threads)
            verbose = False
        else:
            self.setParam("threads", threads)
            verbose = True
        self.setParam("windows", windows)
        self.setParam("continue_N", continue_N)
        sample_FSC = []
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        for x in self.getInput("bedgzInput"):
            outputfile = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".FSC"

            sample_FSC.append(outputfile)
        self.setOutput(
            "sampleOutput", sample_FSC,
        )
        multi_run_len = len(self.getInput("bedgzInput"))
        if verbose:
            for i in range(multi_run_len):
                calc_FSC(
                    bedgzInput=self.getInput("bedgzInput")[i],
                    binInput=self.getParam("binInput"),
                    windows=self.getParam("windows"),
                    continue_N=self.getParam("continue_N"),
                    outputfile=self.getOutput("sampleOutput")[i],

                )
        else:
            args = [
                [
                    self.getInput("bedgzInput")[i],
                    self.getParam("binInput"),
                    self.getParam("windows"),
                    self.getParam("continue_N"),
                    self.getOutput("sampleOutput")[i],

                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=calc_FSC, nCore=maxCore(self.getParam("threads")),
            )
