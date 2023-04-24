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
    bamPath = args.bamPath
    bamInput = []
    files = os.listdir(bamPath)
    for file in files:
        if file.endswith('.bam'):
            bamInput.append(os.path.join(bamPath, file))

            # You can get ocrInput in /cfDNAFE/data/OpenChromatinRegion/
    pathToRunMutation = op.join(dpath, 'scripts/runMutation.R')
    outputdir = args.outputdir
    refGenome = args.refGenome
    genome_reference = args.genome_reference
    baseQ = args.baseQ
    mapQ = args.mapQ
    threads = args.threads
    id = args.id
    runMutation(
        pathToRunMutation=pathToRunMutation,
        bamInput=bamInput,
        genome_reference = genome_reference,
        outputdir=outputdir,
        refGenome=refGenome,
        ID=id,
        baseQ=baseQ,
        mapQ=mapQ,
        threads=threads
    )




def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--bamPath',  required=True, type=str,
            help='this is a file folder, path to bam file')
    parser.add_argument('-g', '--genome_reference', type=str, required=True,
            help="genome reference .fa file")
    parser.add_argument('-o', '--outputdir',  type=str,
            help='output result folder')
    parser.add_argument('-r', '--refGenome',  type=str, choices=['hg19', 'hg38'], default='hg38',
            help='human genome reference hg19 or hg38')
    parser.add_argument('-mq', '--mapQ',  type=int, default=60,
            help='skip bases with baseQ/BAQ smaller than INT')
    parser.add_argument('-bq', '--baseQ',type=int, default=30,
            help='skip alignments with mapQ smaller than INT')
    parser.add_argument('-id', '--id', type=str, default="test",
            help=' group id')
    parser.add_argument('-t', '--threads',  type=int, default=1,
                        help='Whether to use multithreading')
    return parser.parse_args()



if __name__ == '__main__':
    main()
