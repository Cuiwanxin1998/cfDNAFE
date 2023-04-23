library("devtools")
install_github("broadinstitute/ichorCNA")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38",force=TRUE)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19",force=TRUE)
BiocManager::install("MutationalPatterns",force=TRUE)

