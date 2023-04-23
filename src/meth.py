#!/usr/bin/python3 -u

import os.path as op
import os
import argparse
from pathlib import Path
from cfDNAFE import *


def main():
    args = parse_args()
    dpath = str(Path(op.realpath(__file__)).parent.parent)
    bamPath = args.bamPath
    bamInput = []
    files = os.listdir(bamPath)
    for file in files:
        if file.endswith('.bam'):
            bamInput.append(os.path.join(bamPath, file))

    if args.markInput is None:
        markInput = op.join(dpath, 'data/MethMark/Atlas.U25.l4.hg19.bed')
    else:
        markInput = args.markInput
    outputdir = args.outputdir
    mapQuality = args.mapQuality
    minCpG = args.minCpG
    methyThreshold = args.methyThreshold
    unmethyThreshold = args.unmethyThreshold

    threads = args.threads

    runMeth(
        bamInput=bamInput,
        markInput=markInput,
        outputdir=outputdir,
        mapQuality=mapQuality,
        minCpG=minCpG,
        methyThreshold=methyThreshold,
        unmethyThreshold=unmethyThreshold,
        threads=threads
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bamPath', required=True, type=str,
                        help='path to bam file')
    parser.add_argument('-m', '--markInput', type=str,
                        help="regions of mark file, including chrom start end.")
    parser.add_argument('-o', '--outputdir', type=str,
                        help='output result folder')
    parser.add_argument('-mq', '--mapQuality', type=int, default=30,
                        help='min fragment quality')
    parser.add_argument('-mCpG', '--minCpG', type=int, default=4,
                        help='min fragment CpG number')
    parser.add_argument('-mT', '--methyThreshold', type=float, default=0.75,
                        help='thresholds of more than or equal to T percentage methylated CpGs for most methylated '
                             'reads.')
    parser.add_argument('-umT', '--unmethyThreshold', type=float, default=0.25,
                        help='thresholds of less than or equal to T percentage methylated CpGs for most unmethylated '
                             'reads.')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Whether to use multithreading')
    return parser.parse_args()


if __name__ == '__main__':
    main()
