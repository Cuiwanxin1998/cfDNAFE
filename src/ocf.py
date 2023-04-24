#!/usr/bin/python3 -u

import os
import argparse
import os.path as op
from pathlib import Path
from cfDNAFE import *


def main():
    args = parse_args()
    dpath = str(Path(op.realpath(__file__)).parent.parent)
    bedgzPath = args.bedgzPath
    bedgzInput = []
    files = os.listdir(bedgzPath)
    for file in files:
        if file.endswith('.bed.gz'):
            bedgzInput.append(os.path.join(bedgzPath, file))

    if args.ocrInput is None:
        ocrInput = op.join(dpath, 'data/OpenChromatinRegion/7specificTissue.all.OC.bed')
    else:
        ocrInput = args.ocrInput
    outputdir = args.output
    threads = args.threads
    runOCF(
        bedgzInput=bedgzInput,
        ocrInput=ocrInput,
        outputdir=outputdir,
        threads=threads
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bedgzPath',  required=True, type=str,
            help='this is a file folder, path to bedgz file')
    parser.add_argument('-ocr', '--ocrInput', type=str,
            help="regions of open chromosome file")
    parser.add_argument('-o', '--output', type=str,
            help='output result folder')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Whether to use multithreading')
    return parser.parse_args()



if __name__ == '__main__':
    main()
