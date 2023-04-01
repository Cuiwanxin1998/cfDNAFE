from .Base import Base
from .utils import commonError, calc_WPS, maxCore
import os
import math

__metaclass__ = type


class runWPS(Base):
    def __init__(
        self,
        bedgzInput=None,
        tsvInput=None,
        outputdir=None,
        protectInput=120,
        empty=False,
        minSize=120,
        maxSize=180,
        threads=None,
        **kwargs
    ):
        """
        This function is used for computing windowed protection score.
        runWPS(bedgzInput=None, tsvInput=None, outputdir=None, protectInput=120, empty=None, minSize=120, maxSize=180, threads=None, )
        {P}arameters:
            bedgzInput: list, input bed.gz files.
            tsvInput: str, regions of transcript file.
            outputdir: str, output result folder, None means the same folder as input files.
            protectInput: int, base pair protection assumed for elements (default 120).
            empty: bool, keep files of empty blocks (default False).
            minSize: int, Min fragment length (default 120).
            maxSize: int, Max fragment length (default 180).
            threads: int (default None) Whether to use multithreading


        {O}utput:
            sample.ENSG00000000938.12.tsv     chrom: chromatin, pos: position in the genome, covCount:how many reads span this site, startCount: how many reads end point located, WPS
                                                1	1722097	1	0	-2
                                                1	1722096	1	0	-2
                                                1	1722095	1	0	-2
                                                1	1722094	1	0	-2
                                                1	1722093	1	0	-2
                                                1	1722092	1	0	-2
                                                1	1722091	1	0	-2
                                                1	1722090	1	0	-2
                                                1	1722089	1	0	-2
        """
        super(runWPS, self).__init__()
        print(outputdir)
        # set bedgz input
        if bedgzInput is None:
            raise commonError("Parameter bedgzInput must require.")
        else:
            self.setInput("bedgzInput", bedgzInput)

        # set tsv input
        if tsvInput is None:
            raise commonError("tsvInput must be not None!")
        else:
            self.setParam("tsvInput", tsvInput)


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

        # set protect
        self.setParam("protectInput", protectInput)

        # set empty
        self.setParam("empty", empty)

        # set insertsize
        self.setParam("minSize", minSize)
        self.setParam("maxSize", maxSize)
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
                calc_WPS(
                    bedgzInput=self.getInput("bedgzInput")[i],
                    tsvInput=self.getParam("tsvInput"),
                    protectInput=self.getParam("protectInput"),
                    outputfile=os.path.join(
                        self.getOutput("sampleOutputdir")[i],
                        self.getMaxFileNamePrefixV2(self.getInput("bedgzInput")[i]),
                    )
                    + "_%s.tsv.gz",
                    empty=self.getParam("empty"),
                    minSize=self.getParam("minSize"),
                    maxSize=self.getParam("maxSize"),
                )
        else:
            args = [
                [
                    self.getInput("bedgzInput")[i],
                    self.getParam("tsvInput"),
                    os.path.join(
                        self.getOutput("sampleOutputdir")[i],
                        self.getMaxFileNamePrefixV2(self.getInput("bedgzInput")[i]),
                    )
                    + "_%s.tsv.gz",
                    self.getParam("empty"),
                    self.getParam("protectInput"),
                    self.getParam("minSize"),
                    self.getParam("maxSize"),
                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=calc_WPS, nCore=maxCore(math.ceil(self.getParam("threads") / 4)),
            )
