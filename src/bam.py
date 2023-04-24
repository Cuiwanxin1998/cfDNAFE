#!/usr/bin/env python3
import os
import argparse
import os.path as op
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
    if args.blacklist == None:
        blacklistInput = op.join(dpath, 'data/BlackList/hg38-blacklist.v2.bed')
    else:
        blacklistInput = args.blacklist



    genome_reference = args.genome_reference

    outputdir = args.output

    chr = args.chr
    fragFilter = args.fragFilter
    minlen = args.minLen
    maxlen = args.maxLen
    k_mer = args.k_mer
    threads = args.threads
    runBamProcess(
        bamInput=bamInput,
        blacklistInput=blacklistInput,
        outputdir=outputdir,
        genome_reference=genome_reference,
        chr=chr,
        fragFilter=fragFilter,
        minlen=minlen,
        maxlen=maxlen,
        k_mer=k_mer,
        threads=threads
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bamPath', required=True, type=str,
            help='this is a file folder, path to bam file')
    parser.add_argument('-b', '--blacklist', type=str,
            help="regions of blacklist file")
    parser.add_argument('-g', '--genome_reference', required=True, type=str,
            help="genome reference .fa file")
    parser.add_argument('-o', '--output',  type=str,
            help='output result folder')
    parser.add_argument('-c', '--chr',  type=list, default=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr18', 'chr19', 'chr20', 'chr21','chr22', 'chrX']
                        , help="Chromosomes to be processed")
    parser.add_argument('-f', '--fragFilter',  action='store_true',
                        help='Whether filter fragment by length')
    parser.add_argument('-minl', '--minLen',  type=int,
                        help='Min fragment length.')
    parser.add_argument('-maxl', '--maxLen',  type=int,
                        help='Max fragment length.')
    parser.add_argument('-k', '--k_mer', type=int, default=6,
                        help='The number of motifs to analyze')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Whether to use multithreading')
    return parser.parse_args()



if __name__ == '__main__':
    main()
