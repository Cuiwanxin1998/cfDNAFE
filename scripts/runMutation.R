if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("optparse", quietly = TRUE)){
  install.packages("optparse")
}
if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
}
if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
}

library(optparse)

option_list <- list(
  make_option(c("--refGenome"),type = "character", default="hg38", help = "Path to normal WIG file. Default: [%default]"),
  make_option(c("--inputFile"),type = "character", help = "Path to VCF file; Require"), 
  make_option(c("--outputfile1"), type='character', default='test', help='96 mutational profile. Default: [%default]')
  make_option(c("--outputfile2"), type='character', default='test', help='Fit mutation matrix to the COSMIC mutational signatures. Default: [%default]')
  )
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj) 
options(scipen=0, stringsAsFactors=F)
library(MutationalPatterns)
ref_genome = opt$refGenome
if(ref_genome == 'hg19'){
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
  library(ref_genome, character.only = TRUE)
  signatures = get_known_signatures(genome="GRCh37")
}else if(ref_genome == 'hg38'){
  ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
  library(ref_genome, character.only = TRUE)
  signatures = get_known_signatures(genome="GRCh38")
}else{
  stop("The parameter refGenome should be hg19 or hg38.")
}


vcf_files = list.files(opt$inputFile, pattern = ".vcf",full.names=TRUE)
tmp<-strsplit(basename(vcf_files),split=".",fixed=TRUE)
sample_names = unlist(lapply(tmp,head,1)) 


if(length(vcf_files) != length(sample_names)){
  stop("the number of sample Name should be equal the number of vcf files. ")
}

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

  
fit_res <- fit_to_signatures(mut_mat, signatures)
res = fit_res$contribution 
write.table(mut_mat, file=opt$outputfile1, row.names=T, col.names=T, quote=F, sep="\t")
write.table(res, file=opt$outputfile2, row.names=T, col.names=T, quote=F, sep="\t")
