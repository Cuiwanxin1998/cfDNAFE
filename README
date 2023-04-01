# cfDNAFE

* [Introduction](#introduction)
* [Section 1: Installation Tutorial](#section-1-installation-tutorial)
    * [Section 1.1: System requirement](#section-11-system-requirement)
    * [Section 1.2: Create environment and Install Dependencies](#section-12-create-environment-and-install-dependencies)
* [Section 2: Bam File Data Processing (Function:runBamprocess)](#section-2-bam-file-data-processing)
* [Section 3: Fragment Size Ratrio (FSR) and Fragment Size Coverage (FSC) and Fragment Size Distribution (FSD) (Function: runFSR/runFSC/runFSD) ](#section-3-fragment-size-ratio-and-fsr-fragment-size-coverage-fsc--and-fragment-size-distribution-fsd)
* [Section 4: Windows protection score (WPS) (Function:runWPS) ](#section-4-windows-protection-score-wps)
* [Section 5: Orientation-aware cfDNA fragmentation (OCF) (Function:runOCF)](#section-5-orientation-aware-cfdna-fragmentation-ocf)
* [Section 6: Copy Number variations(CNV)(Function:runCNV) ](#section-6-copy-number-variations-cnv)
* [Section 7: Mutation Signature(Function:runMutation) ](#section-7-mutation-signature)
* [Section 8: UXM fragment-level (Function:runMeth) ](#section-8-uxm-fragment-level)



## Introduction

**cfDNAFE(<u>c</u>ell free DNA <u>F</u>eature <u>E</u>extraction)** is a tool for extracting cfDNA features, it contains **<font color=red>End Motif(EDM)</font>**,  **<font color=red>Breakpoint End(BPM)</font>**, **<font color=red>Motif-Diversity Score(MDS)</font>**, **<font color=red>Fragment Size Ratio (FSR)</font>**,**<font color=red>Fragment Size Distribution (FSD)</font>** ,**<font color=red> Windows protection score(WPS)</font>**,**<font color=red> Orientation-aware cfDNA fragmentation(OCF) value</font>**, **<font color=red>Copy Number variations(CNV)</font>**, **<font color=red>UXM fragment-level</font>** and **<font color=red>mutation signature</font>**.

The main functions are as the following picture.
<br/>
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./pics/workFlow.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">cfDNAFE Function</div>
</center>

<br/>



## Section 1: Installation Tutorial

### Section 1.1: System requirement
Since many whole genome sequencing(WGS) or whole genome bisulfite sequencing(WGBS) analysis toolkits are released on Unix/Linux systems, they are based on different programming languages. Therefore, it is very difficult to rewrite all software in one language. Fortunately, the [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/)  program collects many popular python modules and bioinformatics software, so we can install all dependencies via  [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/).

We recommend using [conda/Anaconda](https://www.anaconda.com/) and create a virtual environment to manage all the dependencies.

### Section 1.2: Create environment and Install Dependencies
First, run the following command. The environment will be created and all the dependencies as well as the latest cfDNAFE will be installed. In order to avoid unexpected errors, we recommend that you use R function **install.packages()** and** BiocManager::install()** to download R packages.

```shell
#download cfDNAFE from github(From cfDNAFE you can get the files necessary for the software to run)
git clone https://github.com/Cuiwanxin1998/cfDNAFE.git

python3 setup.py install
#create a virtual environment and activate environment
conda env create -n cfDNAFE -f environment.yml
conda activate cfDNAFE
```
Second, if you want to extract CNVS and mutational signatures, You must install R packages ichorCNA and Mutationalpattern. You need to get into the R language, You can enter the R environment by typing `R` in the Sheel command.
```
#How to download ichorCNA
#Option1(Recommended), you can use devtools
install.packages("devtools")
library(devtools)
install_github("broadinstitute/ichorCNA")
#Option2
#Checkout the latest release of ichorCNA from GitHub
git clone git@github.com:broadinstitute/ichorCNA.git  
 ##install from CRAN
 install.packages("plyr") 
 ## install packages from
 source("https://bioconductor.org/biocLite.R")
 BiocManager::install("HMMcopy")  
 BiocManager::install("GenomeInfoDb")  
 BiocManager::install("GenomicRanges")  
## from the command line and in the directory where ichorCNA github was cloned.
 R CMD INSTALL ichorCNA   


#How to download MutationalPattern
install.packages("BiocManager")
BiocManager::install('MutationalPatterns')
#install dependencies



 

```


## Section 2: Bam File Data Processing

cfDNAFE mainly processes bam file data, which needs to be indexing by samtools. This function input is the initial step of cfDNAFE, which mainly extracts the bed input files required by the following functions and Motif End, Breakpoint End, MDS.

Human reference genome  can be obtained through [UCSC](https://genome.ucsc.edu/index.html), here we only provide [GRCh37/hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and [GRCh38/hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/).


- Example Usage

```Python
from cfDNAFEE import *
import os
#You can get bamInput by simply entering the path where the bam file is stored, as Lines 3 through 9,
bamPath='/path/to/sotred_bamfile/'
bamInput=[]
files=os.listdir(bamPath)
for file in files:
    if file.endswith('.bam'):
        bamInput.append(os.path.join(bamPath, file)) 

#We provide blacklist file for hg19 and hg38, you can find in /cfDNAFE/data/BlackList/
#If you want to other version of blacklist, you can download from here: https://github.com/Boyle-Lab/Blacklist/
blacklistInput='/path/to/cfDNAFE/data/BlackList/blacklist_region'

#You can go to the UCSC website to download the latest versions of the human reference genome hg38 and hg19 as compressed files and unzip them
genome_reference='/path/to/human_genome_reference'

# You can choose the output folder to save the output file, if there is no output folder, cfDNAFE will choose bamInput path as the output folder
outputdir='/path/to/outputdir'

runBamProcess(
        bamInput=bamInput,
        blacklistInput=blacklistInut,
        outputdir=outputdir,
        genome_reference=genome_reference,
)

```

- Detailed parameters

```
bamInput: list, bam file for WGBS or WGS
blacklistInput: str, regions of blacklist file. download from this site (https://github.com/Boyle-Lab/Blacklist/).
outputdir: str, output result folder, None means the same folder as input files.
genome_reference: str, input .fa file.
CHR: list, Chromosomes to be processed(default ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr18', 'chr19', 'chr20', 'chr21','chr22', 'chrX']).
fragFilter: bool, (default False). Whether filter fragment by length, only for paired data.
minLen: int, Min fragment length.
maxLen: int, Max fragment length.
k_mer: int, (default 6) The number of motifs to analyze.
threads: int (default None) Whether to use multithreading
```

- Output Folder Arrangement

```
output_folders/
├──sample1.bed
├──sample1.bed.gz
├──sample1.bed.gz.tbi
├──sample1_EndMotif.txt
├──sample1_BreakPointMotif.txt
├──sample1_MDS.txt 
```

## Section 3: Fragment Size Ratio (FSR) and Fragment Size Coverage (FSC) and Fragment Size Distribution (FSD)
**FSR**: The fragment sizes were used to construct fragmentation profiles with in-house scripts. The FSR was adapted from the [DELFI method](https://www.nature.com/articles/s41586-019-1272-6) and optimized by introducing an extra fragment size group and using improved cutoff. It was was generated using the short/intermediate/long fragments ratios except using different cutoffs: the short, intermediate and long fragments were defined as 65-150bp, 151-220bp and 221-400bp, according to the overall fragment lengths profile in our cohorts. 

**FSC**: FSC was generated using the coverages of short (65-150bp), intermediate (151-260bp), long (261-400bp), and total (65-400bp) cfDNA fragments. The extended ranges allowed the inclusion of broader size regions than what DELFI has reported. The genome was firstly divided into 100 kB bins. Next, the coverage of the four fragment size groups in each 100 kB bin was calculated and corrected by GC content. We then combined the coverages in every 50 contiguous 100 kB bins to calculate the coverage in the corresponding 5 MB (50 × 100 kB) bin. For each fragmentation size group, the scaled coverage score (z-score) in every 5 MB bin was calculated by comparing the variable value against the overall mean value.

We provide 10kb window files for hg19 and hg38, you can find in **/cfDNAFE/data/ChormosomeBins**, if you want to extract them yourself, Here, we give an example how to get a 10kb bin file.

```Python
import pybedtools
chromsize = '/path/to/genome_size'
a = pybedtools.BedTool(chromsize)
binlen = 100000
bins_init = a.window_maker(w=binlen, g=chromsize)
bins_fin = bins_init.filter(
    lambda x: x.chrom
              in [
                  "chr1",
                  "chr2",
                  "chr3",
                  "chr4",
                  "chr5",
                  "chr6",
                  "chr7",
                  "chr8",
                  "chr9",
                  "chr10",
                  "chr11",
                  "chr12",
                  "chr13",
                  "chr14",
                  "chr15",
                  "chr16",
                  "chr17",
                  "chr18",
                  "chr19",
                  "chr20",
                  "chr21",
                  "chr22",
              ]
)
bins_fin.saveas('10kb_bin.bed')
```



**FSD**: The FSD examined fragment length patterns at a high resolution by grouping cfDNA fragments into length bins of 5bp ranging from 65bp and 400bp and calculating the ratio of fragments in each bin at arm level for each chromosome. A total of 41
chromosome arms were examined.

We provide chromosome arms files for hg19 and hg38, you can find in **/cfDNAFE/data/ChormosomeArms**
- Example Usage

```Python
from cfDNAFE import *
import os
#You can get bedgzInput by simply entering the path where the bam file is stored, as Lines 5 through 10
#You can get the .bedgz file by running the function runBamProcess, the runBamProcess have been introduced in the section2
bedgzPath='/path/to/bedgz_file/'
bedgzInput=[]
files=os.listdir(bedgzPath)
for file in files:
    if file.endswith('.bed.gz'):
        bedgzInput.append(os.path.join(bedgzPath, file))
#the binInput you can find in /cfDNAFE/data/ChormosomeBins/.
binInput='/path/to/cfDNAFE/data/ChormosomeBins/hg38_window_10kb.bed'

#the armsInput you can find in /cfDNAFE/data/ChormosomeArms/.
armsInput='/path/to/cfDNAFE/data/ChormosomeArms/hg38.arms.bed'
outputdir='/path/to/outputdir/'
runFSR(
        bedgzInput=bedgzInput,
        binInput=binInput,
        outputdir=outputdir,
        threads=None
)
runFSC(
        bedgzInput=bedgzInput,
        binInput=binInput,
        outputdir=outputdir,
        threads=None
)
runFSD(
        bedgzInput=bedgzInput,
        armsInput=armsInput,
        outputdir=outputdir,
        threads=None
)

```

- Detailed parameters
**runFSR**:
```
 bedgzInput: list, input bed.gz files.
 binInput: str, regions of chromosome N kb bin file.
 windows:int, the window size.
 outputdir: str, output result folder, None means the same folder as input files.
 threads: int, Whether to use multithreading
```
**runFSC**:
```
 bedgzInput: list, input bed.gz files.
 binInput: str, regions of chromosome N kb bin file.
 windows:int, the window size.
 outputdir: str, output result folder, None means the same folder as input files.
 threads: int, Whether to use multithreading
```
**runFSD**:
```
bedgzInput: list, input bed.gz files.
armsInput: str, regions of chromosome arms file.
outputdir: str, output result folder, None means the same folder as input files.
threads: int, Whether to use multithreading
```

- Output Folder Arrangement

```
output_folders/
├──sample1.FSR.txt
├──sample2.FSR.txt
├──sample1.FSC.txt
├──sample2.FSC.txt
├──sample1.FSD.txt
├──sample2.FSD.txt

```


## Section 4: Windows protection score (WPS)
**WPS**: Both outer alignment coordinates of PE data were extracted for properly paired reads. Both end coordinates of SR alignments were extracted when PE data were collapsed to SR data by adapter trimming. A fragment coverage is defined as all positions rag between the two (inferred) , inclusive of endpoints. We define the [windowed protection score (WPS)](https://www.cell.com/cell/fulltext/S0092-8674(15)01569-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS009286741501569X%3Fshowall%3Dtrue) of a window of size k as the number of molecules spanning the window minus those with an endpoint within the window. 

We provide the gene bodies for hg19 and hg38, you can find in **/cfDNAFE/data/transcriptAnno/**

Moreover, we will illustrate how to get gene bodies from gencode annotation files. Users can download gencode annotation files from [gencode database](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/), the commonly used files are [gencode.v19.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz) for hg19 and [gencode.v37.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz) for hg38. Here, we use hg19 as an example.

```R
library(rtracklayer)
library(dplyr)

anno_raw <- import("gencode.v19.annotation.gtf.gz")

# get all genes
anno_raw <- anno_raw[which(anno_raw$type == "gene"), ]

anno <- data.frame(gene_id = anno_raw$gene_id, chr = seqnames(anno_raw), 
                   start = start(anno_raw), end = end(anno_raw), 
                   strand = strand(anno_raw))

# get genome region downstream 10000bp from TSS
for (i in nrow(anno)) {
    if (anno$strand[i] == "+") {
        anno$start = anno$start - 1
        anno$end = anno$start + 10000
    } else {
        anno$start =anno$end + 1 - 10000
        anno$end = anno$end + 1
    }
}

# remove invalid
anno <- anno[which(anno$chr %in% paste0("chr", c(1:22, "X"))), ]
anno <- anno[which(anno$start > 0), ]

write.table(x = anno, file = "transcriptAnno-v19.tsv", sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

```

- Example Usage

```Python
from cfDNAFE import *
import os
#You can get bedgzInput by simply entering the path where the bam file is stored, as Lines 5 through 10
#You can get the .bedgz file by running the function runBamProcess, the runBamProcess have been introduced in the section2
bedgzPath='/path/to/bedgz_file/'
bedgzInput=[]
files=os.listdir(bedgzPath)
for file in files:
    if file.endswith('.bed.gz'):
        bedgzInput.append(os.path.join(bedgzPath, file)) 

#You can find the tsvInput in the /cfDNAFE/data/transcriptAnno/
tsvInput='/path/to/cfDNAFE/data/transcriptAnno/transcriptAnno-v29-hg38.tsv'
outputdir='/path/to/outputdir'
runWPS(
        bedgzInput=bedgzInput,
        tsvInput=tsvInput,
        outputdir=outputdir
)

```

- Detailed parameters

```
bedgzInput: list, input bed.gz files.
tsvInput: str, regions of transcript file.
outputdir: str, output result folder, None means the same folder as input files.
protectInput: int, base pair protection assumed for elements (default 120).
empty: bool, keep files of empty blocks (default False).
minSize: int, Min fragment length (default 120).
maxSize: int, Max fragment length (default 180).
threads: int (default None) Whether to use multithreading
```

- Output Folder Arrangement

```
output_folders/
├──sample1
│   ├──sample1_ENSG000000003.14.tsv.gz
│   ├──sample1_ENSG000000005.5tsv.gz
│   ├──sample1_ENSG0000000419.12.tsv.gz
│   │.....
├──sample2
│   ├──sample2_ENSG000000003.14.tsv.gz
│   ├──sample2_ENSG000000005.5tsv.gz
│   ├──sample2_ENSG0000000419.12.tsv.gz
│   │.....

```

## Section 5: Orientation-aware cfDNA fragmentation (OCF)
**OCF**: To explore the potential in inferring the relative contributions of various tissues in the plasma DNA pool, [Sun *et al.*](https://genome.cshlp.org/content/29/3/418.full) developed a novel approach to measure the differential phasing of upstream (U) and downstream (D) fragment ends in tissue-specific open chromatin regions. They called this strategy orientation-aware cfDNA fragmentation (OCF) analysis. OCF values are based on the differences in U and D end signals in the center of the relevant open chromatin regions. For tissues that contributed DNA into plasma, one would expect much cfDNA fragmentation to have occurred at the nucleosome-depleted region in the center of the corresponding tissue-specific open chromatin regions. In such a region, U and D ends exhibited the highest read densities (i.e., peaks) at ∼60 bp from the center, whereas the peaks for U and D ends were located on the right- and left-hand sides, respectively. Conversely, this pattern would not be expected for tissue-specific open chromatin regions where the corresponding tissue did not contribute DNA into the plasma. Thus measured the differences of U and D end signals in 20-bp windows around the peaks in the tissue-specific open chromatin regions as the OCF value for the corresponding tissue.

We provide tissue specific open chromatin regions for seven tissues, you can find in **/cfDNAFE/data/OpenChromatinRegion/**


- Example Usage

```Python
from cfDNAFE import *
import os
#You can get bedgzInput by simply entering the path where the bam file is stored, as Lines 5 through 10
#You can get the .bedgz file by running the function runBamProcess, the runBamProcess have been introduced in the section2
bedgzPath='/path/to/bedgz_file/'
bedgzInput=[]
files=os.listdir(bedgzPath)
for file in files:
    if file.endswith('.bed.gz'):
        bedgzInput.append(os.path.join(bedgzPath, file)) 

#You can get ocrInput in /cfDNAFE/data/OpenChromatinRegion/
ocrInput='/path/to /cfDNAFE/data/OpenChromatinRegion/7specificTissue.all.OC.bed'
outputdir='/path/to/outputdir'
runOCF(
        bedgzInput=bedgzInput,
        ocrInput=ocrInput,
        outputdir=outputdir
)

```

- Detailed parameters

```
bedgzInput: list, input bed.gz files.
ocrInput: str, regions of open chromosome file.
outputdir: str, output result folder, None means the same folder as input files.
threads: int (default None) Whether to use multithreading.
```

- Output Folder Arrangement

```
output_folders/
├──sample1
│   ├──Breast.sync.end
│   ├──Intestine.sync.end
│   ├──Liver.sync.end
│   ├──Lung.sync.end
│   ├──Ovary.sync.end
│   ├──Placenta.sync.end
│   ├──Tcell.sync.end
│   ├──all.ocf.csv
├──sample2
│   ├──Breast.sync.end
│   ├──Intestine.sync.end
│   ├──Liver.sync.end
│   ├──Lung.sync.end
│   ├──Ovary.sync.end
│   ├──Placenta.sync.end
│   ├──Tcell.sync.end
│   ├──all.ocf.csv

```


## Section 6: Copy Number variations (CNV)
**CNV**: The Copy Number Variation (CNV) profile was calculated using ichorCNA as reported by [Wan *et al.*](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-6003-8). First, the genome of each sample was divided into 1 MB bins. For each bin, the depth after bin-level GC correction was used by a Hidden Markov Model (HMM) to compare against the software baseline. Then, we calculated the log2 ratio for the CNV score.

There are 2 main steps in this part, generating read count coverage information using readCounter from the HMMcopy suite.
Copy number analysis using ichorCNA, you can find in the **/cfDNAFE/scripts/**
Users can find the [input parameters](https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata). In the output results, we can find the log2 transformed CNV from the fourth column in the **sample.cna.seg** file.

- Example Usage **(Fisrt Step: Generating read count)**

```Python
from cfDNAFEE import *
import os

#You can get bamInput by simply entering the path where the bam file is stored, as Lines 5 through 10
bamPath='/path/to/sotred_bamfile/'
bamInput=[]
files=os.listdir(bamPath)
for file in files:
    if file.endswith('.bam'):
        bamInput.append(os.path.join(bamPath, file)) 
		
#generating read count coverage information
outputdir='/path/to/outputdir'
readCounter(
		bamInput=bamInput,
		outputdir=outputdir
)
```
- Example Usage **(Second Step: Copy Number Analysis)**

```
from cfDNAFE import *
import os

#You can get wigInput and ID by simply entering the path where the bam file is stored, as Lines 6 through 13
#You can get the .wig file by running the function readCounter, the readCounter have been introduced in the section6
wigPath='/path/to/sotred_bamfile/'
wigInput=[]
ID=[]
files=os.listdir(wigPath)
for file in files:
    if file.endswith('.wig'):
        wigInput.append(os.path.join(wigPath, file)) 
		ID.append(os.path.splitext(file)[0] )

#You can find the runIchorCNA.R scripts in the /cfDNAFE/scripts
pathTorunIchorCNA='/path/to/runIchorCNA.R'

#You can find the gcWig and mapWig in /cfDNAFE/data/CNVdependency
gcWig='/path/to/cfDNAFE/data/CNVdependecy/gc_hg38_1000kb.wig'
mapWig='/path/to/cfDNAFE/data/CNVdependecy/map_hg38_1000kb.wig'

#You can find the centromere in /cfDNAFE/data/CNVdependecy
centromere='/path/to/cfDNAFE/data/CNVdependecy/GRCh38.GCA_000001405.2_centromere_acen.txt'
runCNV(
        pathTorunIchorCNA=pathTorunIchorCNA,
		wig=wigInput,
		gcWig=gcWig,
		ID=ID,
		mapWig=mapWig,
		centromere=centromere,
		threads=NONE
)
```

- Detailed parameters

```
#readCounter
pathToHreadCounter: str, /path/to/HMMcopy/bin/readCounter.
bamInput: str, regions of blacklist file. download from this site https://github.com/Boyle-Lab/Blacklist/.
outputdir: str, output result folder, None means the same folder as input files.
window_size: int, (default 1000000). input window size.
quality: int, (default 20). remove reads with lower quality
chromosome: str, Chromosomes to be processed(default 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22').
threads: int (default None) Whether to use multithreading

#runCNV
pathTorunIchorCNA: str, /path/to/ichorCNA/scripts/runIchorCNA.R.
wig: list, Path to tumor WIG file
ploidy: str, Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired.
normal: str, Initial normal contamination; can be more than one value if additional normal initializations are desired.
maxCN: int, (default 3). Total clonal CN states.
gcWig: str, Path to GC-content WIG file.
ID: list, Patient ID.
mapWig: str, Path to mappability score WIG file.
centromere:str, File containing Centromere locations; if not provided then will use hg19 version from ichorCNA package.
normalPanel: str, Median corrected depth from panel of normals.
includeHOMD: bool, If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb).
chrs: str, Specify chromosomes to analyze.
chrTrain: str, Specify chromosomes to estimate params.
estimateNormal: bool, Estimate normal.
estimatePloidy: bool, Estimate tumour ploidy.
estimateScPrevalence, bool, Estimate subclonal prevalence.
scStates, str, Subclonal states to consider.
txnE, int, Self-transition probability. Increase to decrease number of segments
txnStrength, int, Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE.
outputdir: str, output result folder, None means the same folder as input files.
threads: int (default None) Whether to use multithreading.
```

- Output Folder Arrangement

```
output_folders/
├──sample1.wig
├──sample1.CNV
├──sample2.wig
├──sample2.CNV

```

## Section 7: Mutation Signature
**Mutation signature**: Each mutational process is thought to leave its own characteristic mark on the genome. For example, AID/APOBEC activity can specifically cause C > T and C > G substitutions at TpCpA and TpCpT sites (of which the underlined nucleotide is mutated.  Thus, patterns of somatic mutations can serve as readout of the mutational processes that have been active and as proxies for the molecular perturbations in a tumour. These [mutational signatures](https://doi.org/10.1038/nature12477) are characterized by a specific contribution of 96 base substitution types with a certain sequence context.

If you do not have vcf files, cfDNAFE provide a function **BamToVcf** to help user convert bam files to vcf files, you can filter out fragments with lower alignment and base quality.

- Example Usage **(run BamToVcf)**

```Python
from cfDNAFE import *
import os 

#You can get bamInput by simply entering the path where the bam file is stored, as Lines 5 through 10
bamPath='/path/to/sotred_bamfile/'
bamInput=[]
files=os.listdir(bamPath)
for file in files:
    if file.endswith('.bam'):
        bamInput.append(os.path.join(bamPath, file)) 

#You can go to the UCSC website to download the latest versions of the human reference genome hg38 and hg19 as compressed files and unzip them
genome_reference='/path/to/human_genome_reference'

outputdir='path/to/output/'

BamToVcf(
            bamInput=bamInput,
            genome_reference=genome_reference,
            outputdir=outputdir,
            threads=None
			)

```


- Example Usage **(runMutation)**

```Python
from cfDNAEE import *

#You can find the runMutation.R scripts in the /cfDNAFE/scripts/
pathToRunMutation='/path/to/cfDNAFE/scripts/runMutation.R'
vcfPath=['path/to/vcf'1, 'path/to/vcf2']
outputdir='path/to/output/'
ID=['EGA5093','EGA5094']
runMutation(
		   pathToRunMutation=pathToRunMutation,
		   vcfInput=vcfInput,
		   outputdir=outputdir,
		   id=ID
)
```

- Detailed parameters

```
#BamToVcf
 bamInput: list, input .bam files
 outputdir: str, output result folder, None means the same folder as input files.
 genome_reference: str, Human reference genome storage path, input .fa or fa.gz file.
 baseQ: int, (default 30) skip bases with baseQ/BAQ smaller than INT
 mapQ: int, (default 60) skip alignments with mapQ smaller than INT
 threads: int (default None) Whether to use multithreading


#runMutation
pathToRunMutation: str, /path/to/RunMutation.R.
vcfInput: str, Path to .vcf files.
outputdir: str, output result folder, None means the same folder as input files.
id: list, Files ID
cutoff: int, a cosine similarity of more than cutoff with an existing COSMIC signature
refGenome: str, human genome reference 'hg19' or 'hg38'
threads: int (default None) Whether to use multithreading


```

- Output Folder Arrangement

```
output_folders/
├──EGA5093.signatures
├──EGA5094.signatures
```

## Section 8: UXM fragment-level
**UXM fragment-level**: Each fragment was annotated as U (mostly unmethylated), M (mostly methylated) or X (mixed) depending on the number of methylated and unmethylated CpGs64. We then calculated, for each genomic region (marker) and across all cell types, the proportion of [U/X/M](https://www.nature.com/articles/s41586-022-05580-6) fragments with at least k CpGs.

cfDNAFE provide a function **runMeth** to calculate the UXM fragment-level. We collect the top 25 differentially unmethylated regions for each cell type comprise a human cell-type-specific methylation atlas in **/cfDNAFE/data/MethMark/**. And cfDNAFE allows users to use custom disease-specific or tissue-specific marks, which requires that the mark contains the chromosome, the start position of the region, the end position of the region and file names ends with **.bed** file


- Example Usage  

```Python
from cfDNAFE import *
import os 

#You can get bamInput by simply entering the path where the bam file is stored, as Lines 5 through 10
bamPath='/path/to/sotred_bamfile/'
bamInput=[]
files=os.listdir(bamPath)
for file in files:
    if file.endswith('.bam'):
        bamInput.append(os.path.join(bamPath, file)) 

#If you do not have  disease-specific or tissue-specific marks, you can find the cell-type-specific methylation mark in /cfDNAFE/data/MethMark/
markInput='/path/to/cfDNAFE/data/MethMark/Markers.U250.hg38.bed'

outputdir='path/to/output/'

runMeth(
            bamInput=bamInput,
            markInput=markInput,
            outputdir=outputdir,
            threads=None
			)

```




- Detailed parameters

```
bamInput: list, input .bam files.
markInput: str, regions of mark file. including chrom start end.
outputdir: str, output result folder, None means the same folder as input files.
mapQuality: int, Min fragment quality. default:30
minCpG: int, Min fragment CpG number. defalut:3
methyThreshold: float, thresholds of more than or equal to T% methylated CpGs for most methylated reads. defalut:0.75
unmethyThreshold: float, thresholds of less than or equal to T% methylated CpGs for most unmethylated reads. default:0.25
threads: int (default None) Whether to use multithreading


```

- Output Folder Arrangement

```
output_folders/
├── sample1.UXM.tsv
├── sample2.UXM.tsv
```




## Authors

* **wx Cui** - *Thesis code writing work* [Central South University](https://cse.csu.edu.cn/)
* **xq Peng**- * Ideas provided and paper writing* [Central South University](https://life.csu.edu.cn/jsxx.jsp?urltype=news.NewsContentUrl&wbtreeid=1815&wbnewsid=3625)


## Version update

* **cfDNAFE v0.1.0** - * News version 0.1.0, 2023.04.01 The first version.* 