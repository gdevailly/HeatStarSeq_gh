###
# shell scripts
### --------------------------

# downloaded from http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/
# 2015-11-04
# file hg19.cage_peak_phase1and2combined_tpm.osc.txt.gz

# file is challenging to load in R. We split headers in two other files, and keep only the matrix
grep -E '^#' hg19.cage_peak_phase1and2combined_tpm.osc.txt | wc -l # 1831
grep -E '^#' hg19.cage_peak_phase1and2combined_tpm.osc.txt > hg19.cage_peak_phase1and2combined_tpm.osc.header_1.txt
wc -l hg19.cage_peak_phase1and2combined_tpm.osc.txt # 203636
# 203636 - 1831
tail -201805 hg19.cage_peak_phase1and2combined_tpm.osc.txt > hg19.cage_peak_phase1and2combined_tpm.osc.reddata.txt
head hg19.cage_peak_phase1and2combined_tpm.osc.reddata.txt | awk '{ print $1, $2,  $3,  $4;}'

head -3 hg19.cage_peak_phase1and2combined_tpm.osc.reddata.txt > hg19.cage_peak_phase1and2combined_tpm.osc.header_2.txt
tail -201802 hg19.cage_peak_phase1and2combined_tpm.osc.reddata.txt > hg19.cage_peak_phase1and2combined_tpm.osc.data.txt
head hg19.cage_peak_phase1and2combined_tpm.osc.data.txt | awk '{ print $1, $2,  $3,  $4;}'

###
# R scripts
### -------------------------
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/CAGE_hg38/")

library(readr)
library(magrittr)
library(GenomicRanges)
library(parallel)
library(ggplot2)
library(vioplot)

load("../../heatmap_shinyApp/data/encode.RData")

a <- Sys.time()
cageData <- list(
    header1 = read_tsv("hg19.cage_peak_phase1and2combined_tpm.osc.header_1.txt", col_names = FALSE),
    header2 = read_tsv("hg19.cage_peak_phase1and2combined_tpm.osc.header_2.txt", col_names = FALSE),
    data = read_tsv("hg19.cage_peak_phase1and2combined_tpm.osc.data.txt", col_names = FALSE)
)
Sys.time() - a # 1 minute

lapply(cageData, dim)
list(
    header1 = head(cageData$header1),
    header2 = cageData$header2[, 1:6],
    data = cageData$data[1:4, 1:6]
    )


regionLocations <- cageData$data[, 1]
length(regionLocations)
a <- Sys.time()
regionLocations %>% strsplit(split = ":|\\.\\.|,", fixed = FALSE) %>% do.call(rbind, .) -> regionLocations
Sys.time() - a # 7 seconds
dim(regionLocations)

colnames(regionLocations) <- c("chr", "start", "end", "strand")
regionLocations %>% as.data.frame(stringsAsFactors = FALSE) -> regionLocations
regionLocations$name <- rep("no_name", length.out = nrow(regionLocations))
regionLocations$score <- rep(0, length.out = nrow(regionLocations))
regionLocations <- regionLocations[, c(1,2,3,5,6,4)]

write.table(regionLocations, file = "all_cage_regions.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# shell, bedtools commands. Checking if all regions are unique. They are
# sortBed -i all_cage_regions.bed > all_cage_regions.s.bed
# bedtools merge -s -c 4 -o count -d -1 -i all_cage_regions.s.bed > all_cage_regions.merged.bed
# wc -l all_cage_regions.bed # 201802
# wc -l all_cage_regions.s.bed # 201802
# wc -l all_cage_regions.merged.bed # 201802 youpi !

regionLocationsGR <- with(regionLocations, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand = strand))

cageMatrix <- as.matrix(cageData$data[, -1])
png(file = "../../Rplots/medianFPM_cage.png")
boxplot(cageMatrix[, sample(1:ncol(cageMatrix), 20)], pch = 20, ylim = c(0, 2), main = "TPM for all regions, 20 samples", ylab = "TPM")
dev.off()

summary(rowSums(cageMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1       82      300     7768     1217 28490000

summary(colSums(cageMatrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 138700  680300  795000  857100  928900 6453000

pdf(file = "../../Rplots/vioPlot_TPMinPeaks_cage.pdf")
vioplot(colSums(cageMatrix), col = "gold")
title(main = "Total TPM in CAGE peaks")
dev.off()

which(colSums(cageMatrix) > 1E6) %>% length # 322

binMatrix_all <- bimMatrix_u1 <- matrix(data = FALSE, nrow = nrow(cageMatrix), ncol = ncol(cageMatrix))
binMatrix_all[cageMatrix > 0] <- TRUE
bimMatrix_u1[cageMatrix >= 1] <- TRUE

cageMatrix[1:8, 1:5]
binMatrix_all[1:8, 1:5]
bimMatrix_u1[1:8, 1:5]

a <- Sys.time()
myCorMatrics <- list(
    all = cor(binMatrix_all),
    u1 = cor(bimMatrix_u1)
)
Sys.time() - a # 40 minutes ...

lapply(myCorMatrics, dim)
save(myCorMatrics, file = "myCorMatrics_cage.RData")


a <- Sys.time()
pdf(file = "../../Rplots/heatmap_all_cage.pdf", width = 13.85, height = 13.85)
heatmap(myCorMatrics$all,
        scale = "none",
        margins = rep(1, 2),
        cexRow = 1,
        cexCol = 1,
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()
Sys.time() - a # 3 minutes


a <- Sys.time()
pdf(file = "../../Rplots/heatmap_all_u1.pdf", width = 13.85, height = 13.85)
heatmap(myCorMatrics$u1,
        scale = "none",
        margins = rep(1, 2),
        cexRow = 1,
        cexCol = 1,
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()
Sys.time() - a # 3 minutes

a <- Sys.time()
myCorMatrics$score <- cor(cageMatrix)
Sys.time() - a # 20 min

lapply(myCorMatrics, dim)
save(myCorMatrics, file = "myCorMatrics_cage.RData")

a <- Sys.time()
pdf(file = "../../Rplots/heatmap_all_socre.pdf", width = 13.85, height = 13.85)
heatmap(myCorMatrics$score,
        scale = "none",
        margins = rep(1, 2),
        cexRow = 1,
        cexCol = 1,
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()
Sys.time() - a # 3 minutes


summary(as.vector(myCorMatrics$all))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0291  0.4238  0.4828  0.4783  0.5372  1.0000

summary(as.vector(myCorMatrics$u1))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04649 0.48740 0.54790 0.54310 0.60250 1.00000

summary(as.vector(myCorMatrics$score))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0004882 0.3098000 0.4534000 0.4538000 0.5984000 1.0000000

pdf(file = "../../Rplots/distrib_of_correlations_cage.pdf")
boxplot(list(
    all_non_null = as.vector(myCorMatrics$all),
    greater_than_1_fpm = as.vector(myCorMatrics$u1),
    with_scores = myCorMatrics$score
), pch = 20, main = "distribution of correlations", ylim = c(0,1))
dev.off()

cageData$header2[, 1:5]

expNames <- cageData$header2[1, -1] %>% as.matrix %>% as.vector()
expNames <- sub("tpm.", "", expNames, fixed = TRUE)
expNames <- gsub("%20", "_", expNames, fixed = TRUE)
expNames <- gsub("%2c", "", expNames, fixed = TRUE)
head(expNames)
tail(expNames)

length(expNames)
grep("_rep", expNames, fixed = TRUE) %>% length # 661 with replicates












