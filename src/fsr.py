#!/usr/bin/python3 -u
from cfDNAFE import *
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


    if args.binInput == None:
        binInput = op.join(dpath, 'data/ChormosomeBins/hg38_window_100kb.bed')
    else:
        binInput = args.binInput


    outputdir = args.output
    windows = args.windows
    threads = args.threads
    continue_N = args.continue_N
    runFSR(
        bedgzInput=bedgzInput,
        binInput=binInput,
        windows=windows,
        continue_N =continue_N,
        outputdir=outputdir,
        threads=threads
    )

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bedgzPath',required=True, type=str,
            help='this is a file folder, path to bedgz file')
    parser.add_argument('-b', '--binInput', type=str,
            help="regions of chromosome N kb bin file")
    parser.add_argument('-w', '--windows', default=100000, type=int,
            help="the window size.")
    parser.add_argument('--continue_N', '-c', default=50, type=int,
                        help="continue window number.")
    parser.add_argument( '-o', '--output', type=str,
            help='output result folder')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Whether to use multithreading, (default None) ')
    return parser.parse_args()



if __name__ == '__main__':
    main()