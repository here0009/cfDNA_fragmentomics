#!/usr/bin/env Rscript
# file: median_table.R
# author: Derek Wong, Ph.D
# date: August 4th, 2021

library(tidyverse)
library(GenomicRanges)
args <- commandArgs(trailingOnly = TRUE)
# output_file <- args[1]
# output_tag <- args[2]
# Read in files and combine into data frame
# control dir
dir_lst <- c('220615_control','220608','220629_control','220610_cancer','220628_cancer','220718_merged_cancer','220803','220729','220805')

outdir <- "results/delfi_220809"

# outdir <- "results/delfi/total"
files <- lapply(dir_lst, function(dir) {
  list.files(dir, pattern = "_5Mb_bins.txt", recursive = TRUE, full.names = TRUE)
})
files <- unlist(files)
print("The input files are:")
print(files)
bins.list <- lapply(files, read.delim)
tib.list <- lapply(bins.list, as_tibble)
rm(bins.list)
# concate 5mb tables together
df.fr <- tib.list %>%
  bind_rows() %>% dplyr::select(everything())

coverage <- df.fr %>% ungroup() %>% group_by(sample_id) %>%
  dplyr::summarize(hqbases_analyzed = 100*sum(nfrags)*2,
                   depth = hqbases_analyzed/(504*5e6)
  )
write.table(df.fr, file.path(outdir, paste0("total_results.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
#    sample_id                                  hqbases_analyzed depth
#    <chr>                                                 <dbl> <dbl>
#  1 S002_SZ20220425057WHB-5_cfdna_genome_72630      20032506600  7.95
#  2 S085_SZ20220425040WHB-c_cfdna_genome_72923      16291067600  6.46
# Generate healthy median
median_table <- df.fr %>%
  group_by(bin) %>% 
  dplyr::summarize(sample_id = "median",
            seqnames = unique(seqnames),
            arm = unique(arm),
            start = unique(start),
            end = unique(end),
            gc = unique(gc),
            median_frag_gc = median(frag_gc, na.rm=TRUE),
            median_ratio=median(short/long, na.rm=TRUE),
            sd_ratio=sd(short/long, na.rm=TRUE),
            median_ratio_corrected=median(short_corrected/long_corrected, na.rm=TRUE),
            sd_ratio_corrected=sd(short_corrected/long_corrected, na.rm=TRUE),
            median_ratio_centered=median(ratio_centered, na.rm=TRUE),
            sd_ratio_centered=sd(ratio_centered, na.rm=TRUE),
            median_coverage=median(coverage, na.rm=TRUE),
            sd_coverage=sd(coverage, na.rm=TRUE),
            median_coverage_corrected=median(coverage_corrected, na.rm=TRUE),
            sd_coverage_corrected=sd(coverage_corrected, na.rm=TRUE),
            median_coverage_centered=median(coverage_centered, na.rm=TRUE),
            sd_coverage_centered=sd(coverage_centered, na.rm=TRUE),
            median_combined=median(combined, na.rm=TRUE),
            sd_combined=sd(combined, na.rm=TRUE),
            median_combined_centered=median(combined_centered, na.rm=TRUE),
            sd_combined_centered=sd(combined_centered, na.rm=TRUE),
            median_mode_size=median(mode_size, na.rm=TRUE),
            median_mean_size=median(mean_size, na.rm=TRUE),
            sd_mean_size=sd(mean_size, na.rm=TRUE),
            median_median_size=median(median_size, na.rm=TRUE),
            sd_median_size=sd(median_size, na.rm=TRUE),
            median_q25_size=median(q25_size, na.rm=TRUE),
            median_q75_size=median(q75_size, na.rm=TRUE))

write.table(median_table, file.path(outdir, paste0("median.hg19.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
save(median_table, file=file.path(outdir, paste0("median.hg19.rda")))

# Generate correlations
summary_df <- df.fr %>% ungroup() %>% group_by(sample_id) %>%
  dplyr::summarize(
            ratio_cor=cor(ratio, median_table$median_ratio, method="pearson", use="complete.obs"),
            ratio_corrected_cor=cor(ratio_corrected, median_table$median_ratio_corrected, method="pearson", use="complete.obs"),
            ratio_centered_cor=cor(ratio_centered, median_table$median_ratio_centered, method="pearson", use="complete.obs"),
            coverage_cor=cor(coverage, median_table$median_coverage, method="pearson", use="complete.obs"),
            coverage_corrected_cor=cor(coverage_corrected, median_table$median_coverage_corrected, method="pearson", use="complete.obs"),
            coverage_centered_cor=cor(coverage_centered, median_table$median_coverage_centered, method="pearson", use="complete.obs"),
            combined_cor=cor(combined, median_table$median_combined, method="pearson", use="complete.obs"),
            combined_centered_cor=cor(combined_centered, median_table$median_combined_centered, method="pearson", use="complete.obs"),
            nfrags = sum(nfrags),
            mode_size=unique(mode_size),
            mean_size=unique(mean_size),
            median_size=unique(median_size),
            q25_size=unique(q25_size),
            q75_size=unique(q75_size),
            hqbases_analyzed = 100*sum(nfrags)*2,
            coverage = hqbases_analyzed/(504*5e6)
  )

write.table(summary_df, file.path(outdir, paste0("summary_correlations.txt")), sep = "\t", row.names = FALSE, quote = FALSE)

# Plot profiles
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
df.fr$arm <- factor(df.fr$arm, levels=armlevels)
median_table$arm <- factor(median_table$arm, levels=armlevels)

arm <- df.fr %>% group_by(arm) %>%
  dplyr::summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "10q", "", "12q", "", "16q",
                         "", "17q", "", "18q",
                         "", "", "", "",
                         "", ""),
                       c("10p", "10q", "12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

# Generate Fragmentation and Coverage plots
g1 <- ggplot(df.fr, aes(x=bin, y=ratio_centered, group=sample_id, color="red")) + 
  geom_line(size=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr, sample_id=="", size=0.75, alpha=1))
g1 <- g1 + geom_line(data=median_table, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
g1 <- g1 + labs(x="", y="Fragmentation profile\n", color="")
g1 <- g1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
g1 <- g1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g1 <- g1 + mytheme
g1
ggsave(file.path(outdir, paste0("fragment.pdf")), g1, width=15, height=3, units="in")

c1 <- ggplot(df.fr, aes(x=bin, y=coverage_centered, group=sample_id, color="red")) + 
  geom_line(size=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr, sample_id==""), size=0.75, alpha=1)
c1 <- c1 + geom_line(data=median_table, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
c1 <- c1 + labs(x="", y="Coverage profile\n", color="")
c1 <- c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
c1 <- c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
c1 <- c1 + mytheme
c1
ggsave(file.path(outdir, paste0("coverage.pdf")), c1, width=15, height=3, units="in")

b1 <- ggplot(df.fr, aes(x=bin, y=combined_centered, group=sample_id, color="red")) + 
  geom_line(size=0.5, alpha=0.5)
  #geom_line(data=subset(df.fr, sample_id==""), size=0.75, alpha=1)
b1 <- b1 + geom_line(data=median_table, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
b1 <- b1 + labs(x="", y="Combined profile\n", color="")
b1 <- b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
b1 <- b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
b1 <- b1 + mytheme
ggsave(file.path(outdir, paste0("combined.pdf")), b1, width=15, height=3, units="in")
