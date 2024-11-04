data <- read.table("SRR2130016.bed_alpha_satellite.tsv.gz")
WPS <- data$V5
x <- seq(34443000, 34446000)
plot(x = x, y = WPS, type = "l", 
     xlab = "genome coordinate", ylab = "WPS", lwd = "1")