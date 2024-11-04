# source: https://xwanglabthu.github.io/cfDNApipe/
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