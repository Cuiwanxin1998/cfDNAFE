#!/usr/bin/python3 -u
from cfDNAFE import *
import os
import sys
import argparse
import os.path as op
from pathlib import Path

def main():
    args = parse_args()
    dpath = str(Path(op.realpath(__file__)).parent.parent)
    bedgzPath = args.bedgzPath
    bedgzInput = []
    files = os.listdir(bedgzPath)
    for file in files:
        if file.endswith('.bed.gz'):
            bedgzInput.append(os.path.join(bedgzPath, file))

            # You can get ocrInput in /cfDNAFE/data/OpenChromatinRegion/
    if args.tsvInput == None:
        tsvInput = op.join(dpath, 'data/TranscriptAnno/transcriptAnno-hg38-1kb.tsv')
    else:
        tsvInput = args.tsvInput

    if args.wpstype == 'L':
        protectInput = 120
        minSize = 120
        maxSize = 180
    else:
        protectInput = 16
        minSize = 35
        maxSize = 80
    outputdir = args.output
    empty = args.empty
    threads = args.threads
    runWPS(
        bedgzInput=bedgzInput,
        tsvInput=tsvInput,
        outputdir=outputdir,
        protectInput=protectInput,
        empty=empty,
        minSize=minSize,
        maxSize=maxSize,
        threads=threads
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bedgzPath', required=True, type=str,
            help='this is a file folder, path to bedgz file')
    parser.add_argument('-tsv', '--tsvInput',  type=str,
            help="regions of open chromosome file")
    parser.add_argument('-o', '--output', type=str,
            help='output result folder')
    parser.add_argument('-w', '--wpstype', type=str, default='L', choices=['L', 'S'],
                        help='output result folder')
    parser.add_argument('-empty', '--empty',  type=bool, default=False, choices=[True, False],
                        help='keep files of empty blocks')
    parser.add_argument( '-t', '--threads', type=int, default=1,
                        help='Whether to use multithreading')
    return parser.parse_args()



if __name__ == '__main__':
    main()
