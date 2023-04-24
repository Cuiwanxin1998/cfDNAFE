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
    bamPath = args.bamPath
    bamInput = []
    files = os.listdir(bamPath)
    pathTorunIchorCNA = op.join(dpath, 'scripts/runIchorCNA.R')
    for file in files:
        if file.endswith('.bam'):
            bamInput.append(os.path.join(bamPath, file))
    if args.seqinfo == 'hg38':
        if args.gcWig == None:
            gcWig = op.join(dpath, 'data/CNVdependency/gc_hg38_1000kb.wig')
        else:
            gcWig = args.gcWig

        if args.mapWig == None:
            mapWig = op.join(dpath, 'data/CNVdependency/map_hg38_1000kb.wig')
        else:
            mapWig = args.mapWig

        if args.centromere == None:
            centromere = op.join(dpath, 'data/CNVdependency/GRCh38.GCA_000001405.2_centromere_acen.txt')
        else:
            centromere = args.centromere
        seqinfo = op.join(dpath, 'data/CNVdependency/seqinfo_hg38_ucsc.rds')
    if args.seqinfo == 'hg19':
        if args.gcWig == None:
            gcWig = op.join(dpath, 'data/CNVdependency/gc_hg19_1000kb.wig')
        else:
            gcWig = args.gcWig

        if args.mapWig == None:
            mapWig = op.join(dpath, 'data/CNVdependency/map_hg19_1000kb.wig')
        else:
            mapWig = args.mapWig

        if args.centromere == None:
            centromere = op.join(dpath, 'data/CNVdependency/GRCh37.p13_centromere_UCSC-gapTable.txt')
        else:
            centromere = args.centromere
        seqinfo = op.join(dpath, 'data/CNVdependency/seqinfo_hg19_ucsc.rds')
    outputdir = args.outputdir
    quality = args.quality
    ploidy = args.ploidy
    window_size = args.window_size
    chromosome = args.chromosome
    normal = args.normal
    maxCN = args.maxCN
    normalPanel = args.normalPanel
    includeHOMD = args.includeHOMD
    chrs = args.chrs
    chrTrain = args.chrTrain
    estimateNormal = args.estimateNormal
    estimatePloidy = args.estimatePloidy
    estimateScPrevalence = args.estimateScPrevalence
    scStates = args.scStates
    txnE = args.txnE
    txnStrength = args.txnStrength
    threads = args.threads

    runCNV(
        pathTorunIchorCNA=pathTorunIchorCNA,
        bamInput=bamInput,
        outputdir=outputdir,
        gcWig=gcWig,
        mapWig=mapWig,
        centromere=centromere,
        window_size=window_size,
        quality=quality,
        chromosome=chromosome,
        ploidy=ploidy,
        normal=normal,
        maxCN=maxCN,
        normalPanel=normalPanel,
        seqinfo=seqinfo,
        includeHOMD=includeHOMD,
        chrs=chrs,
        chrTrain=chrTrain,
        estimateNormal=estimateNormal,
        estimatePloidy=estimatePloidy,
        estimateScPrevalence=estimateScPrevalence,
        scStates=scStates,
        txnE=txnE,
        txnStrength=txnStrength,
        threads=threads
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bamPath', required=True, type=str,
                        help='this is a file folder, path to bam file')
    parser.add_argument('-c', '--chromosome', type=str,
                        default='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22',
                        help="Chromosomes to be processed")
    parser.add_argument('-o', '--outputdir', type=str,
                        help='output result folder')
    parser.add_argument('-w', '--window_size', type=int, default=1000000,
                        help='input window size')
    parser.add_argument('-q', '--quality', type=int, default=20,
                        help=' remove reads with lower quality')
    parser.add_argument('-ploidy', '--ploidy', type=str, default="'c(2)'",
                        help='Initial tumour ploidy')


    parser.add_argument('-normal', '--normal', type=str, default="'c(0.95, 0.99, 0.995, 0.999)'",
                        help='Initial normal contamination')
    parser.add_argument('-maxCN', '--maxCN', type=str, default=3,
                        help='Total clonal CN states')
    parser.add_argument('-gcWig', '--gcWig', type=str,
                        help='Path to GC-content WIG file')
    parser.add_argument('-mapWig', '--mapWig', type=str,
                        help='Path to mappability score WIG file')
    parser.add_argument('-centromere', '--centromere', type=str,
                        help='Median corrected depth from panel of normals')
    parser.add_argument('-includeHOMD', '--includeHOMD', type=bool, default=False,
                        help='If FALSE, then exclude HOMD state. Useful when using large bins')

    parser.add_argument('-chrs', '--chrs', type=str, default="'c(1:22)'",
                        help='Specify chromosomes to analyze')
    parser.add_argument('-chrTrain', '--chrTrain', type=str, default="'c(1:22)'",
                        help='Estimate normal.')
    parser.add_argument('-normalPanel', '--normalPanel', type=str,
                        help='Median corrected depth from panel of normals')
    parser.add_argument('-estimateNormal', '--estimateNormal', type=bool, default=True,
                        help='Specify chromosomes to estimate params')

    parser.add_argument('-estimatePloidy', '--estimatePloidy', type=bool, default=True,
                        help='Estimate tumour ploidy')

    parser.add_argument('-estimateScPrevalence', '--estimateScPrevalence', type=bool, default=True,
                        help='Estimate subclonal prevalence')

    parser.add_argument('-scStates', '--scStates', type=str, default="'c()'",
                        help='Subclonal states to consider')

    parser.add_argument('-txnE', '--txnE', type=int, default=0.9999999,
                        help='Self-transition probability. Increase to decrease number of segments')

    parser.add_argument('-txnStrength', '--txnStrength', type=int, default=1000000,
                        help='Transition pseudo-counts. Exponent should be the same as the number of decimal places of')

    parser.add_argument('-seqinfo', '--seqinfo', type=str, default="hg38", choices=['hg38', 'hg19'],
                        help='genome build')

    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Whether to use multithreading')
    return parser.parse_args()


if __name__ == '__main__':
    main()
