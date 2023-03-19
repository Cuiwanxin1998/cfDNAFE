# cfDNAFFE

* [Introduction](#introduction)
* [Section 1: Installation Tutorial](#section-1-installation-tutorial)
    * [Section 1.1: System requirement](#section-11-system-requirement)
    * [Section 1.2: Create environment and Install Dependencies](#section-12-create-environment-and-install-dependencies)
    * [Section 1.3: Activate Environment and Use cfDNApipe](#section-13-activate-environment-and-use-cfdnapipe)
* [Section 2: WGS Data Processing (Function:Bamprocess)](#section-2-wgs-data-processing-function-bamprocess)
* [Section 3: Fragment Size Coverage (FSC) and Fragment Size Distribution (FSD) (Function: runFSC/runFSD) ](#section-3-fragment-size-coverage-fsc-and-fragment-size-distribution-fsd-function-runfsc-runfsd)
* [Section 4: Windows protection score(WPS) (Function:runWPS) ](#section-4-windows-protection-score-wps-function-runwps)
* [Section 5: Orientation-aware cfDNA fragmentation(OCF) (Function:runOCF)](#section-5-orientation-aware-cfDNA-fragmentation-OCF-Function-runOCF)
* [Section 6: Copy Number variations(CNV)(Function:runCNV) ](#section-6-copy-number-variations-CNV-Function-runCNV)
* [Section 7: Mutation Signature(Function:runMutation) ](#section-7-mutation-signature-function-runmutation)



## Introduction

**cfDNAFFE(<u>c</u>ell <u>f</u>ree <u>DNA</u> <u>F</u>ragment <u>F</u>eature <u>E</u>extraction)** is a tool for extracting cfDNA fragmentation features, tt contains <font color=red>**End Motif(EDM)**</font>,  <font color=red>**Breakpoint End(BPM)**</font>, <font color=red>**Motif-Diversity Score(MDS)**</font>, <font color=red>**Fragment Size Coverage (FSC)**</font>,<font color=red>**Fragment Size Distribution (FSD)**</font> ,<font color=red>** Windows protection score(WPS)**</font>,<font color=red>** Orientation-aware cfDNA fragmentation(OCF) value **</font>, <font color=red>**Copy Number variations(CNV)**</font>  and <font color=red>**mutation signature**</font>.

The main functions are as the following picture.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);"
    src="./pics/workflow.jpg">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">cfDNAFFE Function</div>
</center>


## Section 1: Installation Tutorial

### Section 1.1: System requirement
Since many WGS analysis toolkits are released on Unix/Linux systems, they are based on different programming languages. Therefore, it is very difficult to rewrite all software in one language. Fortunately, the [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/)  program collects many popular python modules and bioinformatics software, so we can install all dependencies via  [conda](https://docs.conda.io/en/latest/)/[bioconda](http://bioconda.github.io/).

We recommend using [conda/Anaconda](https://www.anaconda.com/) and create a virtual environment to manage all the dependencies.

### Section 1.2: Create environment and Install Dependencies
First, run the following command. The environment will be created and all the dependencies as well as the latest cfDNAFFE will be installed. In order to avoid unexpected errors, we recommend that you use R function **install.packages()** and** BiocManager::install()** to download R packages.

```shell
#create a virtual environment
conda env create -n cfDNAFFE -f environment.yml

#install dependencies

#ichorCNA: Users can follow https://github.com/broadinstitute/ichorCNA/wiki/Installation

#MutationalPatterns
BiocManager::install('MutationalPatterns')

```


## Section 2: WGS Data Processing (*Function:<font color=red>Bamprocess</font>*)

cfDNAFFE mainly processes bam file data, which needs to be sorted by samtools. This function input is the initial step of cfDNAFFE, which mainly extracts the bed input files required by the following functions and Motif End, Breakpoint End, MDS.

- Example Usage

```Python
from cfDNAFEE import *

bamInput=['path_to_sample1','path_to_sample2']
blacklistInut=['path_to_blacklist_region']
outputdir='path_to_outputdir'
genome_reference='path_to_human_genome_reference'
CHR = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr18', 'chr19', 'chr20', 'chr21',
               'chr22']

runBamProcess(
        bamInput=bamInput,
        blacklistInput=blacklistInut,
        outputdir=outputdir,
        genome_reference=genome_reference,
        CHR=CHR,
        mapQuality=30,
        fragFilter=False,
        minLen=None,
        maxLen=None,
        k_mer=6,
        threads=None
)

```

- Detailed parameters

```
bamInput: list, input .bam files.
blacklistInput: str, regions of blacklist file. download from this site (https://github.com/Boyle-Lab/Blacklist/).
outputdir: str, output result folder, None means the same folder as input files.
genome_reference: str, input .fa or fa.gz file.
CHR: list, Chromosomes to be processed(default ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr18', 'chr19', 'chr20', 'chr21',
'chr22', 'chrX']).
empty: bool, keep files of empty blocks (default False).
fragFilter: bool, (default False). Whether filter fragment by length, only for paired data.
minLen: int, Min fragment length.
maxLen: int, Max fragment length.
k_mer: int, (default 6) The number of motifs to analyze.
threads: int (default None) Whether to use multithreading
```

- Output Folder Arrangement

```
output_folders/
├──sample1/
│   ├──sample1.bed
│   ├──sample1.bed.gz
│   ├──sample1.bed.gz.tbi
│   ├──sample1_EndMotif.txt
│   ├──sample1_BreakPointMotif.txt
│   ├──sample1_MDS.txt
├──sample2/
│   ├──sample2.bed
│   ├──sample2.bed.gz
│   ├──sample2.bed.gz.tbi
│   ├──sample2_EndMotif.txt
│   ├──sample2_BreakPointMotif.txt
│   ├──sample2_MDS.txt
```

## Section 3: Fragment Size Coverage (FSC) and Fragment Size Distribution (FSD) (*Function:<font color=red>runFSC/runFSD</font>*)
**FSC**: *The fragment sizes were used to construct fragmentation profiles with in-house scripts. The FSC was adapted from the [DELFI method](https://www.nature.com/articles/s41586-019-1272-6) and optimized by introducing an extra fragment size group and using improved cutoff. It was generated using the coverages of short (65-150bp), intermediate (151-260bp), long (261-400bp), and total (65-400bp) cfDNA fragments. The extended ranges allowed the inclusion of broader size regions than what DELFI has reported. The genome was firstly divided into 100 kB bins. Next, the coverage of the four fragment size groups in each 100 kB bin was calculated and corrected by GC content. We then combined the coverages in every 50 contiguous 100 kB bins to calculate the coverage in the corresponding 5 MB (50 × 100 kB) bin. For each fragmentation size group, the scaled coverage score (z-score) in every 5 MB bin was calculated by comparing the variable value against the overall mean value.*

Here, we give an example how to get a 10kb bin file. Human genome chromosome length can be obtained through [UCSC](https://genome.ucsc.edu/index.html), here we only provide [GRCh37/hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and [GRCh38/hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/).

```Python
import pybedtools
chromsize = 'path_to_genome_size'
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



**FSD**: *The FSD feature examined the coverage of cfDNA fragments ranging from 65 bps to 400 bps in 5 bp stepwise (e.g., 65-69 bps, 70-74 bps…) at every chromosome arm. The raw coverage score of FSD was also scaled into the z-score by comparing the variable value against the overall mean value. *
- Example Usage

```Python
from cfDNAFEE import *

bedgzInput=['sample1.bed.gz','sample1.bed.gz']
binInput='path_to_10kb_bin.bed'
armsInput='path_to_arms.bed'
outputdir='path_to_outputdir'
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

```
bedgzInput: list, input bed.gz files.
binInput: str, regions of chromosome 10kb bin file.
armsInput: str, regions of chromosome arms file.
outputdir: str, output result folder, None means the same folder as input files.
threads: int (default None) Whether to use multithreading
```

- Output Folder Arrangement

```
output_folders/
├──sample1.FSC.txt
├──sample2.FSC.txt
├──sample1.FSD.txt
├──sample2.FSD.txt

```


## Section 4: Windows protection score(WPS) (*Function:<font color=red>runWPS</font>*)
**WPS**: *Both outer alignment coordinates of PE data were extracted for properly paired reads. Both end coordinates of SR alignments were extracted when PE data were collapsed to SR data by adapter trimming. A fragment coverage is defined as all positions rag between the two (inferred) , inclusive of endpoints. We define the [windowed protection score (WPS)](https://www.cell.com/cell/fulltext/S0092-8674(15)01569-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS009286741501569X%3Fshowall%3Dtrue) of a window of size k as the number of molecules spanning the window minus those with an endpoint within the window. *

we will illustrate how to get gene bodies from gencode annotation files. Users can download gencode annotation files from [gencode database](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/), the commonly used files are [gencode.v19.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz) for hg19 and [gencode.v37.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz) for hg38. Here, we use hg19 as an example.

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
from cfDNAFEE import *

bedgzInput=['sample1.bed.gz','sample1.bed.gz']
tsvInput='path_to_transcript.tsv'
outputdir='path_to_outputdir'
runWPS(
        bedgzInput=None,
        tsvInput=None,
        outputdir=None,
        protectInput=120,
        empty=False,
        minSize=120,
        maxSize=180,
        threads=None
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

## Section 5: Orientation-aware cfDNA fragmentation(OCF) (*Function:<font color=red>runOCF</font>*)
**OCF**: *To explore the potential in inferring the relative contributions of various tissues in the plasma DNA pool, [Sun *et al.*](https://genome.cshlp.org/content/29/3/418.full) developed a novel approach to measure the differential phasing of upstream (U) and downstream (D) fragment ends in tissue-specific open chromatin regions. They called this strategy orientation-aware cfDNA fragmentation (OCF) analysis. OCF values are based on the differences in U and D end signals in the center of the relevant open chromatin regions. For tissues that contributed DNA into plasma, one would expect much cfDNA fragmentation to have occurred at the nucleosome-depleted region in the center of the corresponding tissue-specific open chromatin regions. In such a region, U and D ends exhibited the highest read densities (i.e., peaks) at ∼60 bp from the center, whereas the peaks for U and D ends were located on the right- and left-hand sides, respectively. Conversely, this pattern would not be expected for tissue-specific open chromatin regions where the corresponding tissue did not contribute DNA into the plasma. Thus measured the differences of U and D end signals in 20-bp windows around the peaks in the tissue-specific open chromatin regions as the OCF value for the corresponding tissue.*

Tissue-specific open chromatin regions can be found in [Supplemental_Code.zip](https://genome.cshlp.org/content/29/3/418/suppl/DC1)

- Example Usage

```Python
from cfDNAFEE import *

bedgzInput=['sample1.bed.gz','sample1.bed.gz']
ocrInput='path_to_open_chromosome_region.bed'
outputdir='path_to_outputdir'
runOCF(
        bedgzInput=None,
        ocrInput=None,
        outputdir=None,
        threads=None,
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


## Section 6: Copy Number variations(CNV)(*Function:<font color=red>runCNV</font>*)
**CNV**: *The Copy Number Variation (CNV) profile was calculated using ichorCNA as reported by [Wan *et al.*](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-6003-8). First, the genome of each sample was divided into 1 MB bins. For each bin, the depth after bin-level GC correction was used by a Hidden Markov Model (HMM) to compare against the software baseline. Then, we calculated the log2 ratio for the CNV score.*

There are 2 main steps in this part, generating read count coverage information using readCounter from the HMMcopy suite.
Copy number analysis using ichorCNA R package. Users can find the [input parameters](https://github.com/broadinstitute/ichorCNA/tree/master/inst/extdata). In the output results, we can find the log2 transformed CNV from the fourth column in the **sample.cna.seg** file.

- Example Usage

```Python
from cfDNAFEE import *

#generating read count coverage information
bamInput=['sample1_sort.bam','sample2_sort.bam']
pathToreadCounter='path_to_HMMcopy_ReadCounter'
outputdir='path_to_outputdir'
readCounter(
        pathToreadCounter=None,
		bamInput=bamInput,
		outputdir=outputdir,
		window_size=1000000
)
#Copy number analysis
pathTorunIchorCNA='path_to_runIchorCNA.R'
wig=['sample1.wig','sample2.wig']
gcWig='path_to_gc_content.wig'
mapWig='path_to_mappability.wig'
ID=['sample1','sample2']
runCNV(
        pathTorunIchorCNA=pathTorunIchorCNA,
		wig=wig,
		ploidy="'c(2)'",
		normal="'c(0.95, 0.99, 0.995, 0.999)'",
		maxCN=3,
		gcWig=gcWig,
		ID=ID,
		mapWig=mapWig,
		centromere=None,
		normalPanel=None,
		includeHOMD=False,
		chrs="'c(1:22)'",
		chrTrain="'c(1:22)'",
		estimateNormal=True,
		estimatePloidy=True,
		estimateScPrevalence=False,
		scStates="'c()'",
		txnE=0.9999999,
		txnStrength=1000000,
		outputdir=None,
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
├──sample1.cna.seg
├──sample2.wig
├──sample2.cna.seg

```

## Section 7:Mutation Signature(*Function:<font color=red>runMutation</font>*)
**Mutation signature**: *Each mutational process is thought to leave its own characteristic mark on the genome. For example, AID/APOBEC activity can specifically cause C > T and C > G substitutions at TpCpA and TpCpT sites (of which the underlined nucleotide is mutated.  Thus, patterns of somatic mutations can serve as readout of the mutational processes that have been active and as proxies for the molecular perturbations in a tumour. These [mutational signatures](https://doi.org/10.1038/nature12477) are characterized by a specific contribution of 96 base substitution types with a certain sequence context.*

The R file Run_mutation in the scripts folder.

- Example Usage

```Python
from cfDNAFEE import *

pathToRunMutation='path_to_RunMutation.R'
vcfInput=['path/to/vcfInput1', 'path/to/vcfInput2']
outputdir='path/to/output/'
ID=['EGA5093','EGA5094']
runMutation(
       pathToRunMutation=pathToRunMutation,
	   vcfInput=vcfInput,
	   outputdir=outputdir,
	   id=ID,
	   threads=None,
)
```

- Detailed parameters

```
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
├──EGA5093_MutationSignature.txt
├──EGA5094_MutationSignature.txt
```
