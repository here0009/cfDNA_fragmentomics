# generate cytosine_ref.rds file
# https://github.com/cancer-genomics/reproduce_lucas_wflow/issues/1
# conda activate r_env
# conda install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg19
# conda install -c bioconda bioconductor-homo.sapiens
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Homo.sapiens)
cytosines <- vmatchPattern("C", Hsapiens)
saveRDS(cytosines, file="cytosine_ref.rds")