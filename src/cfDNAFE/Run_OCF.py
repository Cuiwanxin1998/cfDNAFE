from .Base import Base
from .utils import commonError, calc_OCF, maxCore
import os
import math

__metaclass__ = type


class runOCF(Base):
    def __init__(
        self,
        bedgzInput=None,
        ocrInput=None,
        outputdir=None,
        threads=1,
        **kwargs
    ):
        """
        This function is used for computing windowed protection score.
        runWPS(bedgzInput=None, ocrInput=None, outputdir=None, threads=None)
        {P}arameters:
            bedgzInput: list, input bed.gz files.
            ocrInput: str, regions of open chromosome file.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int (default None) Whether to use multithreading.
        """
        super(runOCF, self).__init__()
        # set bedgz input
        if bedgzInput is None:
            raise commonError("Parameter bedgzInput must require.")
        else:
            self.setInput("bedgzInput", bedgzInput)

        # set tsv input
        if ocrInput is None:
            raise commonError("ocrInput must be not None!")
        else:
            self.setParam("ocrInput", ocrInput)


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

        dirs = []
        for x in self.getInput("bedgzInput"):
            newdir = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x))
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            dirs.append(newdir)
        self.setOutput(
            "sampleOutputdir", dirs,
        )
        multi_run_len = len(self.getInput("bedgzInput"))
        if verbose:
            for i in range(multi_run_len):
                calc_OCF(
                    bedgzInput=self.getInput("bedgzInput")[i],
                    ocrInput=self.getParam("ocrInput"),
                    outputdir=self.getOutput("sampleOutputdir")[i],

                )
        else:
            args = [
                [
                    self.getInput("bedgzInput")[i],
                    self.getParam("ocrInput"),
                    self.getOutput("sampleOutputdir")[i],

                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=calc_OCF, nCore=maxCore(self.getParam("threads")),
            )
