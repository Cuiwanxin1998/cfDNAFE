library(optparse)

option_list <- list(
  make_option(c("--refGenome"),type = "character", default="hg38", help = "Path to normal WIG file. Default: [%default]"),
  make_option(c("--inputFile"),type = "character", help = "Path to VCF file; Require"), 
  make_option(c("--cutoff"), type = "numeric", default=0.85, help = "a cosine similarity of more than cutoff with an existing COSMIC signature. Default: [%default]"),
  make_option(c("--outputdir"),type = 'character', default = NULL, help='Output Dir.'),
  make_option(c("--id"), type='character', default='test', help='Files ID. Default: [%default]')
  )
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
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


vcf_files = list.files(inputFile, pattern = ".vcf",full.names=TRUE)
tmp<-strsplit(basename(vcf_files),split=".",fixed=TRUE)
sample_names = unlist(lapply(tmp,head,1)) 

if(len(vcf_files) != len(sample_names)){
  stop("the number of sample Name should be equal the number of vcf files. ")
}

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
 


nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = opt$cutoff)
colnames(nmf_res$signatures)
fit_res <- fit_to_signatures(mut_mat, signatures)
res = fit_res$contribution
dir.create(paste0(outDir, "/"), recursive = TRUE)
outFile = paste0(outDir, '/', opt$id, "mutation_signatures.txt")
write.table(res, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")