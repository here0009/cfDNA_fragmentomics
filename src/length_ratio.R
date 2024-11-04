#!/usr/bin/env Rscript
# source: https://github.com/derekwong90/fragmentomics
library(getopt)
args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
bedfile <- options.args[1]
outdir <- options.args[2]
id <- options.args[3]

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Homo.sapiens)
library(biovizBase)
library(data.table)
library(tibble)
library(tidyverse)
library(KernSmooth)

# bedfile <- "filtered_bam_no_blacklist.bed"
# id <- "test"
# outdir <- "."q

# Step 1 GC correction and get 10kb bins
## GC correct function
gc.correct <- function(coverage, bias) {
  i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias, na.action = na.omit)
  coverage.model <- loess(predict(coverage.trend, i) ~ i, na.action = na.omit)
  coverage.pred <- predict(coverage.model, bias)
  coverage.corrected <- coverage - coverage.pred + median(coverage, na.rm=TRUE)
}

## GC prediction function
gc.pred <- function(coverage, bias) {
  i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend <- loess(coverage ~ bias, na.action = na.omit)
  coverage.model <- loess(predict(coverage.trend, i) ~ i, na.action = na.omit)
  coverage.pred <- predict(coverage.model, bias)
}

AB <- read.table("/source_code/delfi/delfi_scripts/data/hic_compartments_100kb_ebv_2014.txt", header=TRUE)
load(file="/source_code/delfi/delfi_scripts/data/filters.hg19.rda")
load(file="/source_code/delfi/delfi_scripts/data/gaps.hg19.rda")

bed <- fread(bedfile)
# setnames(bed, c("chr", "start", "end", "name", "score", "strand", "end_motif_1", "end_motif_2"))
setnames(bed, c('chr', 'start', 'end', 'name', 'score', 'strand', 'breakpoint_F', 'end_motif_F', 'breakpoint_R', 'end_motif_R', 'gc'))
bed <- bed[,start:=start+1]
frags <- makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)
w.all <- width(frags)
frags <- frags[which(w.all >= 90 & w.all <= 220)] # we only analyse fragments with a length between 90 and 220 bp
# frags$gc2 <- GCcontent(Hsapiens, unstrand(frags))
# cor(frags$gc2, frags$gc)
# frags$gc_diff <- abs(frags$gc2 - frags$gc)
# mean(frags$gc_diff)
# frags:
  #             seqnames            ranges strand |                   name
  #                <Rle>         <IRanges>  <Rle> |            <character>
  #         [1]     chr1     750115-750256      - | E100041297L1C034R019..
  #         [2]     chr1     750104-750268      - | E100041297L1C001R034..
  #         [3]     chr1     750114-750273      - | E100041297L1C030R014..
  #         [4]     chr1     750109-750278      - | E100041297L1C006R018..
  #         [5]     chr1     750121-750284      + | E100041297L1C011R024..
  #         ...      ...               ...    ... .                    ...
  # [103932074]    chr21 48119444-48119563      - | E100041297L1C038R016..
  # [103932075]    chr21 48119437-48119597      + | E100041297L1C024R013..
  # [103932076]    chr21 48119450-48119581      - | E100041297L1C038R040..
  # [103932077]    chr21 48119437-48119618      - | E100041297L1C028R028..
  # [103932078]    chr21 48119446-48119635      - | E100041297L1C014R024..
  #                 score end_motif_1 end_motif_2       gc
  #             <integer> <character> <character> <matrix>
  #         [1]        48      CCAATG      CATCTG 0.549296
  #         [2]        48      TGGAGC      GGTTTG 0.563636
  #         [3]        48      GCCAAT      TGATAG 0.550000
  #         [4]        48      CCACAG      CCCTTT 0.552941
  #         [5]        48      TGGGCC      CTGATG 0.548780
  #         ...       ...         ...         ...      ...
  # [103932074]        49      AGGCGG      ACCTTG 0.533333
  # [103932075]        48      ACCGTA      GCTTTC 0.509317
  # [103932076]        48      AGCAGC      CATTCT 0.507576
  # [103932077]        43      ACCGTA      CACAGA 0.510989
  # [103932078]        27      GCGGAG      TCTCGG 0.521053
# load 10kb tiling data
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)
chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))
tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]
# tcmers
# GRanges object with 70 ranges and 1 metadata column:
#        seqnames              ranges strand |       type
#           <Rle>           <IRanges>  <Rle> |   <factor>
#    [1]     chr1 121535434-124535434      * | centromere
#    [2]     chr1             0-10000      * | telomere
#    [3]     chr1 249240621-249250621      * | telomere
#    [4]     chr2   92326171-95326171      * | centromere
#    [5]     chr2             0-10000      * | telomere
#    ...      ...                 ...    ... .        ...
#   [66]    chr21             0-10000      * | telomere
#   [67]    chr21   48119895-48129895      * | telomere
#   [68]    chr22             0-10000      * | telomere
#   [69]    chr22   13000000-16000000      * | centromere
#   [70]    chr22   51294566-51304566      * | telomere
#   -------
#   seqinfo: 24 sequences from hg19 genome
arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
arms$arm <- armlevels
#write.table( x = data.frame(arms), file = "arms.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg19))]
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]
seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(Hsapiens, AB)
## These bins had no coverage
AB <- AB[-c(8780, 13665)]
# GRanges object with 26236 ranges and 4 metadata columns:
#           seqnames            ranges strand |     eigen      domain         arm
#              <Rle>         <IRanges>  <Rle> | <numeric> <character> <character>
#       [1]     chr1     700000-799999      * |  -1.62476        open          1p
#       [2]     chr1     800000-899999      * |  -1.64961        open          1p
#       [3]     chr1     900000-999999      * |  -1.60062        open          1p
#       [4]     chr1   1000000-1099999      * |  -1.57911        open          1p
#       [5]     chr1   1100000-1199999      * |  -1.51213        open          1p
#       ...      ...               ...    ... .       ...         ...         ...
#   [26232]    chr22 50700000-50799999      * | -0.702034        open         22q
#   [26233]    chr22 50800000-50899999      * | -0.954097        open         22q
#   [26234]    chr22 50900000-50999999      * | -0.924447        open         22q
#   [26235]    chr22 51000000-51099999      * | -0.630956        open         22q
#   [26236]    chr22 51100000-51199999      * |  0.101323      closed         22q
#                 gc
#           <matrix>
#       [1]  0.44116
#       [2]  0.56645
#       [3]  0.62135
#       [4]  0.60348
#       [5]  0.61372
#       ...      ...
#   [26232]  0.59216
#   [26233]  0.52411
#   [26234]  0.57211
#   [26235]  0.52015
#   [26236]  0.53297
# Filters
fragments <- frags[-queryHits(findOverlaps(frags, filters.hg19))]
w <- width(fragments) # length of fragments
frag.list <- split(fragments, w)
# frag.list is the fragments grouped by the length of the fragments, so we can count the number of short and long fragments
# GRangesList object of length 6:
# $`90`
# GRanges object with 83257 ranges and 5 metadata columns:
#           seqnames            ranges strand |                   name     score
#              <Rle>         <IRanges>  <Rle> |            <character> <integer>
#       [1]     chr1     779695-779784      + | E100041297L1C009R020..        60
#       [2]     chr1     867714-867803      - | E100041297L1C001R024..        60
#       [3]     chr1     930357-930446      + | E100041297L1C016R021..        60
#       [4]     chr1     949677-949766      + | E100041297L1C039R026..        60
#       [5]     chr1     955846-955935      + | E100041297L1C026R013..        60
#       ...      ...               ...    ... .                    ...       ...
#   [83253]    chr21 47794657-47794746      + | E100041297L1C019R032..        58
#   [83254]    chr21 47836267-47836356      + | E100041297L1C004R028..        60
#   [83255]    chr21 47971668-47971757      + | E100041297L1C011R004..        60
#   [83256]    chr21 48092968-48093057      + | E100041297L1C036R023..        60
#   [83257]    chr21 48100050-48100139      - | E100041297L1C039R030..        60
#           end_motif_1 end_motif_2       gc
#           <character> <character> <matrix>
#       [1]      GGTGCT      CACTGA 0.422222
#       [2]      CAGCGT      GCCCCT 0.677778
#       [3]      CCACCA      ACTGTC 0.477778
#       [4]      ACCTGA      GGAGCT 0.644444
#       [5]      CGCAGC      CTTTCC 0.800000
#       ...         ...         ...      ...
#   [83253]      GGGAGG      TGTTGC 0.600000
#   [83254]      GGCCAT      GGCATT 0.466667
#   [83255]      GTTCCG      CACGTC 0.677778
#   [83256]      CAACAA      GCCACA 0.255556
#   [83257]      GGGAGA      CTCCAA 0.533333
counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
if(min(w) > 90) {
  m0 <- matrix(0, ncol=min(w) - 90, nrow=nrow(counts),
               dimnames=list(rownames(counts), 90:(min(w)-1)))
  counts <- cbind(m0, counts)
}

if(max(w) < 220) {
  m1 <- matrix(0, ncol=220 - max(w), nrow=nrow(counts),
               dimnames=list(rownames(counts), (max(w)+1):220))
  counts <- cbind(counts, m1)
}
# dim(counts) is (26236,131)
# cols of counts are the width of the fragments, rows are the AB segments
olaps <- findOverlaps(fragments, AB)
# Hits object with 97761342 hits and 0 metadata columns:
#              queryHits subjectHits
#              <integer>   <integer>
#          [1]         1           1
#          [2]         2           1
#          [3]         3           1
#          [4]         4           1
#          [5]         5           1
#          ...       ...         ...
#   [97761338] 103772877       25904
#   [97761339] 103772878       25904
#   [97761340] 103772879       25904
#   [97761341] 103772880       25904
#   [97761342] 103772881       25904
#   -------
#   queryLength: 103773158 / subjectLength: 26236
# query is the fragments, subject is the AB segments
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))
# > tail(bin.list)
# GRangesList object of length 6:
# $`26231`
# GRanges object with 2803 ranges and 5 metadata columns:
#          seqnames            ranges strand |                   name     score
#             <Rle>         <IRanges>  <Rle> |            <character> <integer>
#      [1]    chr22 50599827-50600005      - | E100041297L1C024R007..        60
#      [2]    chr22 50599907-50600085      - | E100041297L1C013R023..        60
#      [3]    chr22 50599920-50600092      + | E100041297L1C039R037..        60
#      [4]    chr22 50599963-50600098      - | E100041297L1C029R024..        60
#      [5]    chr22 50599962-50600157      + | E100041297L1C032R026..        60
#      ...      ...               ...    ... .                    ...       ...
#   [2799]    chr22 50699674-50699858      + | E100041297L1C009R003..        60
#   [2800]    chr22 50699854-50699980      + | E100041297L1C042R019..        60
#   [2801]    chr22 50699834-50700007      - | E100041297L1C019R028..        60
#   [2802]    chr22 50699863-50699998      - | E100041297L1C035R017..        60
#   [2803]    chr22 50699882-50700045      - | E100041297L1C003R034..        60
#          end_motif_1 end_motif_2       gc
#          <character> <character> <matrix>
#      [1]      CCCAGT      GCCCCA 0.642458
#      [2]      TGTGCT      GGCGGC 0.653631
#      [3]      GCCACG      CAGGGT 0.653179
#      [4]      CCCTTG      TTATTT 0.625000
#      [5]      GCCCTT      CACCAC 0.622449
#      ...         ...         ...      ...
#   [2799]      ATACAG      CCGGGA 0.794595
#   [2800]      CCCGGT      CCGCCC 0.748031
#   [2801]      CGAGCC      GGAGCC 0.758621
#   [2802]      CGGCGC      TCCCCG 0.757353
#   [2803]      CGTCTT      GACGAG 0.780488
#   -------
#   seqinfo: 23 sequences from an unspecified genome; no seqlengths

# group fragments by their AB segments
# bingc <- rep(NA, length(bin.list))
bingc <- rep(NA, length(AB))
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc)) # the mean GC is not accurate, for length of fragments are not the same in each bin, however, it is a good approximation
# bingc is the mean gc of the fragments in each AB segment

### Get modes
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
modes <- Mode(w)
medians <- median(w)
q25 <- quantile(w, 0.25)
q75 <- quantile(w, 0.75)

short <- rowSums(counts[,1:61])
long <- rowSums(counts[,62:121])
ratio <- short/long
ratio[is.nan(ratio)] <- NA
ratio[is.infinite(ratio)] <- NA
nfrags <- short+long
coverage <- nfrags/sum(nfrags, na.rm=TRUE)

## GC correction and prediction
short.corrected <- gc.correct(short, bingc)
long.corrected <- gc.correct(long, bingc)
short.predicted <- gc.pred(short, bingc)
long.predicted <- gc.pred(long, bingc)
nfrags.predicted <- gc.pred(short+long, bingc)
coverage.predicted <- gc.pred(coverage, bingc)

short.corrected <- ifelse(short.corrected <= 0, 0, short.corrected)
long.corrected <- ifelse(long.corrected <= 0, NA, long.corrected)
nfrags.corrected <- short.corrected+long.corrected
ratio.corrected <- short.corrected/long.corrected
coverage.corrected <- nfrags.corrected/sum(nfrags.corrected, na.rm=TRUE)
combined <- ratio.corrected*coverage.corrected

## Append fragment information
AB$short <- short
AB$long <- long
AB$ratio <- ratio
AB$nfrags <- short+long
AB$coverage <- coverage
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected
AB$coverage.corrected <- coverage.corrected
AB$combined <- combined
AB$short.predicted <- short.predicted
AB$long.predicted <- long.predicted
AB$nfrags.predicted <- nfrags.predicted
AB$ratio.predicted <- short.predicted/long.predicted
AB$coverage.predicted <- coverage.predicted

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc
AB$id <- id

for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]
tib.list <- as_tibble(AB)
tib.list <- tib.list %>% dplyr::select(-matches("X"))
write.table(tib.list, file.path(outdir, paste0(id, "_100kb_bins.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

## Step 2 get 5Mb bins

## Plot GC Correction metrics
pdf(file = file.path(outdir, paste0(id, "_GC_metrics.pdf")))
par(mfrow=c(2,2))
## short
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$short,
              main = "short",
              xlab = "frag_GC", 
              ylab = "short")
## short corrected
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$short.corrected,
              main = "short corrected",
              xlab = "frag_GC", 
              ylab = "short_corrected")
## short vs short predicted
smoothScatter(x = tib.list$short.predicted, 
              y = tib.list$short, 
              main = "short predicted vs actual",
              xlab = "short_predicted", 
              ylab = "short")
## short corrected vs short predicted
smoothScatter(x = tib.list$short.predicted, 
              y = tib.list$short.corrected, 
              main = "short predicted vs corrected",
              xlab = "short_predicted", 
              ylab = "short_corrected")
## long
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$long, 
              main = "long",
              xlab = "frag_GC", 
              ylab = "long")
## long corrected
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$long.corrected, 
              main = "corrected long",
              xlab = "frag_GC", 
              ylab = "long_corrected")
## long vs long predicted
smoothScatter(x = tib.list$long.predicted, 
              y = tib.list$long, 
              main = "long predicted vs actual",
              xlab = "long_predicted", 
              ylab = "long")
## long corrected vs long predicted
smoothScatter(x = tib.list$long.predicted, 
              y = tib.list$long.corrected, 
              main = "long predicted vs corrected",
              xlab = "long_predicted", 
              ylab = "long_corrected")
## total fragments
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$nfrags, 
              main = "nfrags",
              xlab = "frag_GC", 
              ylab = "nfrags")
## corrected total fragments
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$nfrags.corrected, 
              main = "corrected nfrags",
              xlab = "frag_GC", 
              ylab = "nfrags_corrected")
## fragments vs predicted fragments
smoothScatter(x = tib.list$nfrags.predicted, 
              y = tib.list$nfrags, 
              main = "nfrags predicted vs actual",
              xlab = "nfrags_predicted", 
              ylab = "nfrags")
## corrected fragments vs predicted fragments
smoothScatter(x = tib.list$nfrags.predicted, 
              y = tib.list$nfrags.corrected, 
              main = "nfrags predicted vs corrected",
              xlab = "nfrags_predicted", 
              ylab = "nfrags_corrected")
## ratios
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$ratio, 
              main = "ratios",
              xlab = "frag_gc", 
              ylab = "ratio")
## corrected ratios
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$ratio.corrected, 
              main = "corrected ratios",
              xlab = "frag_gc", 
              ylab = "ratio_corrected")
## predicted ratios
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$ratio.predicted, 
              main = "predicted ratios",
              xlab = "frag_gc", 
              ylab = "ratio_predicted")
## bin GC vs frag GC
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$C.G, 
              main = "GC content",
              xlab = "frag_GC", 
              ylab = "bin_GC")
## coverage
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$coverage, 
              main = "coverage",
              xlab = "frag_GC", 
              ylab = "coverage")
## corrected coverage
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$coverage.corrected, 
              main = "corrected coverage",
              xlab = "frag_GC", 
              ylab = "coverage_corrected")
## predicted coverage
smoothScatter(x = tib.list$frag.gc, 
              y = tib.list$coverage.predicted, 
              main = "predicted coverage",
              xlab = "frag_GC", 
              ylab = "coverage_predicted")
dev.off()

## Set arm levels
df.fr2 <- tib.list
rm(tib.list)
armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")


## Combine adjacent 100kb bins to form 5mb bins. We count starting from
## the telomeric end and remove the bin closest to the centromere if it is
## smaller than 5mb.
df.fr2 <- df.fr2 %>% group_by(arm) %>%
  mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                          ceiling((1:length(arm))/50) ))

df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
  dplyr::summarize(short2=sum(short, na.rm=TRUE),
            long2=sum(long, na.rm=TRUE),
            short.corrected2=sum(short.corrected, na.rm=TRUE),
            long.corrected2=sum(long.corrected, na.rm=TRUE),
            gc=mean(C.G, na.rm=TRUE),
            frag.gc2=mean(frag.gc, na.rm=TRUE),
            nfrags2=sum(nfrags, na.rm=TRUE),
            nfrags.corrected2=sum(nfrags.corrected, na.rm=TRUE),
            short.var=var(short.corrected, na.rm=TRUE),
            long.var=var(long.corrected, na.rm=TRUE),
            nfrags.var=var(nfrags.corrected, na.rm=TRUE),
            mode_size=unique(mode, na.rm=TRUE),
            mean_size=unique(mean, na.rm=TRUE),
            median_size=unique(median, na.rm=TRUE),
            q25_size=unique(quantile.25, na.rm=TRUE),
            q75_size=unique(quantile.75, na.rm=TRUE),
            start=start[1],
            end=rev(end)[1],
            binsize = n())

df.fr3$ratio2 <- df.fr3$short2/df.fr3$long2
df.fr3$ratio.corrected2 <- df.fr3$short.corrected2/df.fr3$long.corrected2
df.fr3$coverage2 <- df.fr3$short2/sum(df.fr3$nfrags2)
df.fr3$coverage.corrected2 <- df.fr3$short.corrected2/sum(df.fr3$nfrags.corrected2)
df.fr3$combined2 <- df.fr3$ratio.corrected2*df.fr3$coverage.corrected2

## Assign bins
df.fr3 <- df.fr3 %>% filter(binsize==50)
df.fr3 <- df.fr3 %>% group_by(id) %>% mutate(bin = 1:length(id))

## Center and scale fragment and coverage ratios
df.fr3$ratio.centered <- ((df.fr3$ratio.corrected2 - mean(df.fr3$ratio.corrected2))/sd(df.fr3$ratio.corrected2))*0.01
df.fr3$coverage.centered <- ((df.fr3$coverage.corrected2 - mean(df.fr3$coverage.corrected2))/sd(df.fr3$coverage.corrected2))*0.01
df.fr3$combined.centered <- ((df.fr3$combined2 - mean(df.fr3$combined2))/sd(df.fr3$combined2))*0.01

### Rename and reorder dataframe
df.fr3 <- df.fr3 %>%
  dplyr::rename(
    sample_id = id,
    short = short2,
    long = long2,
    short_corrected = short.corrected2,
    long_corrected = long.corrected2,
    ratio = ratio2,
    ratio_corrected = ratio.corrected2,
    nfrags = nfrags2,
    nfrags_corrected = nfrags.corrected2,
    coverage_corrected = coverage.corrected2,
    short_var = short.var,
    long_var = long.var,
    nfrags_var = nfrags.var,
    ratio_centered = ratio.centered,
    coverage_centered = coverage.centered,
    combined = combined2,
    combined_centered = combined.centered,
    frag_gc = frag.gc2,
    coverage = coverage2
  )

df.fr3 <- df.fr3 %>%
  relocate(sample_id, seqnames, arm, start, end, gc, frag_gc, short, long, nfrags, ratio,
           short_corrected, long_corrected, nfrags_corrected, ratio_corrected, ratio_centered, 
           coverage, coverage_corrected, coverage_centered, combined, combined_centered,
           short_var, long_var, nfrags_var, mode_size,mean_size, median_size, q25_size, q75_size, 
           binsize, bin)

write.table(df.fr3, file.path(outdir, paste0(id, "_5Mb_bins.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

### Step 3 


### Step 4 plot
## Set themes and plot layouts
mytheme <- theme_classic(base_size=12) + theme(
  axis.text.x = element_blank(),
  axis.ticks.x=element_blank(),
  strip.text.x = element_text(size=11),
  strip.text.y = element_text(size=12),
  axis.title.x = element_text(face="bold", size=17),
  axis.title.y = element_text(size=15),
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=15),
  legend.position = "none",
  legend.title = element_text(size=10),
  legend.text = element_text(size=10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background=element_rect(fill="white", color="white"),
  panel.spacing.x=unit(0.1, "lines"))

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr3$arm <- factor(df.fr3$arm, levels=armlevels)
# healthy_median$arm <- factor(healthy_median$arm, levels=armlevels)

arm <- df.fr3 %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

## Plot Fragmentation profile
f1 <- ggplot(df.fr3, aes(x=bin, y=ratio_centered, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.75)
# f1 <- f1 + geom_line(data=healthy_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
f1 <- f1 + labs(x="", y="Fragmentation profile\n", color="")
f1 <- f1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
f1 <- f1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
f1 <- f1 + mytheme
ggsave(file.path(outdir, paste0(id, "_fragment.pdf")), f1, width=15, height=3, units="in")

## Plot short fragment coverage profile
c1 <- ggplot(df.fr3, aes(x=bin, y=coverage_centered, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.75)
# c1 <- c1 + geom_line(data=healthy_median, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
c1 <- c1 + mytheme
ggsave(file.path(outdir, paste0(id, "_coverage.pdf")), c1, width=15, height=3, units="in")

## Plot combined profile
b1 <- ggplot(df.fr3, aes(x=bin, y=combined_centered, group=sample_id, color="red")) + 
  geom_line(size=0.75, alpha=0.75)
# b1 <- b1 + geom_line(data=healthy_median, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
b1 <- b1 + labs(x="", y="Combined profile\n", color="")
b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
b1 <- b1 + mytheme
ggsave(file.path(outdir, paste0(id, "_combined.pdf")), b1, width=15, height=3, units="in")