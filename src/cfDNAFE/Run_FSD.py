from .Base import Base
from .utils import commonError, calc_FSD, maxCore
import os
import math

__metaclass__ = type


class runFSD(Base):
    def __init__(
        self,
        bedgzInput=None,
        armsInput=None,
        outputdir=None,
        threads=1,
        **kwargs
    ):
        """
        This function is used for FSD.
        The FSD examined fragment length patterns at a high resolution by grouping cfDNA fragments into length bins of 5bp ranging from 65bp and 400bp and calculating the ratio of fragments in each bin at arm level for each chromosome.
        runFSD(bedgzInput=None, armsInput=None, outputdir=None, threads=None,)
        {P}arameters:
            bedgzInput: list, input bed.gz files.
            armsInput: str, regions of chromosome arms file.
            outputdir: str, output result folder, None means the same folder as input files.
            threads: int, Whether to use multithreading

        {O}utput:
        outputfile1: sample.FSD                     65bin_ratio(65~69bp,70~74bp,...,396-400bp)
                                   0.0,0.0,0.006622516556291391,0.0016556291390728477,0.0033112582781456954,0.0033112582781456954...
                                   0.0,0.0,0.0011111111111111111,0.001851851851851852,0.0014814814814814814,0.0033333333333333335,0.0014814814814814814...
                                   0.0,0.0,0.0013169446883230904,0.001755926251097454,0.0013169446883230904,0.0021949078138718174,0.0013169446883230904...
                                                                ......

        """
        super(runFSD, self).__init__()
        # set bedgz input
        if bedgzInput is None:
            raise commonError("Parameter bedgzInput must require.")
        else:
            self.setInput("bedgzInput", bedgzInput)

        # set tsv input
        if armsInput is None:
            raise commonError("armsInput must be not None!")
        else:
            self.setParam("armsInput", armsInput)


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

        sample_FSD = []
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
        for x in self.getInput("bedgzInput"):
            outputfile = os.path.join(self.getOutput("outputdir"), self.getMaxFileNamePrefixV2(x)) + ".FSD"

            sample_FSD.append(outputfile)
        self.setOutput(
            "sampleOutput", sample_FSD,
        )
        multi_run_len = len(self.getInput("bedgzInput"))
        if verbose:
            for i in range(multi_run_len):
                calc_FSD(
                    bedgzInput=self.getInput("bedgzInput")[i],
                    armsfile=self.getParam("armsInput"),
                    outputfile=self.getOutput("sampleOutput")[i],

                )
        else:
            args = [
                [
                    self.getInput("bedgzInput")[i],
                    self.getParam("armsInput"),
                    self.getOutput("sampleOutput")[i],

                ]
                for i in range(multi_run_len)
            ]
            self.multiRun(
                args=args, func=calc_FSD, nCore=maxCore(self.getParam("threads")),
            )
