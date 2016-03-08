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
library(dplyr)
library(GenomicRanges)
library(parallel)
library(ggplot2)
library(vioplot)

#load("../../heatmap_shinyApp/data/encode.RData")

a <- Sys.time()
cageData <- list(
    header1 = read_tsv("hg19.cage_peak_phase1and2combined_tpm.osc.header_1.txt", col_names = FALSE),
    header2 = read_tsv("hg19.cage_peak_phase1and2combined_tpm.osc.header_2.txt", col_names = FALSE, col_types = do.call(paste0, as.list(rep("c", 1830)))),
    data = read_tsv("hg19.cage_peak_phase1and2combined_tpm.osc.data.txt", col_names = FALSE)
)
Sys.time() - a # 1 minute
cageData <- lapply(cageData, as.data.frame)
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
cageMatrix <- as.matrix(cageData$data[, -1])
regionLocationsGR <- with(regionLocations, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand = strand))

# write.table(regionLocations, file = "all_cage_regions.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# shell, bedtools commands. Checking if all regions are unique. They are
# sortBed -i all_cage_regions.bed > all_cage_regions.s.bed
# bedtools merge -s -c 4 -o count -d -1 -i all_cage_regions.s.bed > all_cage_regions.merged.bed
# wc -l all_cage_regions.bed # 201802
# wc -l all_cage_regions.s.bed # 201802
# wc -l all_cage_regions.merged.bed # 201802 youpi !


# summary(rowSums(cageMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 1       82      300     7768     1217 28490000

# summary(colSums(cageMatrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 138700  680300  795000  857100  928900 6453000

# t0 <- Sys.time()
# myCor <- cor(cageMatrix)
# Sys.time() - t0 # 20 minutes ...
#
# save(myCor, file = "myCorMatrics_cage.RData")
load("myCorMatrics_cage.RData")


# a <- Sys.time()
# pdf(file = "../../Rplots/heatmap_all_socre.pdf", width = 13.85, height = 13.85)
# heatmap(myCorMatrics$score,
#         scale = "none",
#         margins = rep(1, 2),
#         cexRow = 1,
#         cexCol = 1,
#         breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
#         col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
#         useRaster = TRUE
# )
# dev.off()
# Sys.time() - a # 3 minutes

summary(as.vector(myCor))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0004882 0.3098000 0.4534000 0.4538000 0.5984000 1.0000000

cageData$header2[, 1:5]

expNames <- cageData$header2[1, -1] %>% as.matrix %>% as.vector()
expNames <- sub("tpm.", "", expNames, fixed = TRUE)
expNames <- gsub("%20", "_", expNames, fixed = TRUE)
expNames <- gsub("%2c", "", expNames, fixed = TRUE)
expNames <- gsub("%2e", "", expNames, fixed = TRUE)
expNames <- gsub("%2f", "", expNames, fixed = TRUE)
expNames <- gsub("%3a", "_", expNames, fixed = TRUE)
expNames <- gsub("%28", "", expNames, fixed = TRUE)
expNames <- gsub("%29", "", expNames, fixed = TRUE)
expNames <- gsub("%27", "", expNames, fixed = TRUE)
expNames <- gsub("%2b", "+", expNames, fixed = TRUE)
expNames <- gsub("%5e", "", expNames, fixed = TRUE)
expNames <- gsub("_-_", "_", expNames, fixed = TRUE)
expNames <- gsub("___", "_", expNames, fixed = TRUE)
expNames <- gsub("__", "_", expNames, fixed = TRUE)

head(expNames)
tail(expNames)

length(expNames)
# write.table(expNames, file = "cage_exp_name.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
grep("_rep", expNames, fixed = TRUE) %>% length # 661 with replicates

# dataset is to big, we remove time course data
# expNames[grep("[0-9]hr", expNames, fixed = FALSE)]
# expNames[grep("day[0-9][0-9]", expNames, fixed = FALSE)]
# expNames[grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
timeCourse <- expNames[grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
notTimeCourse <- expNames[-grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
length(timeCourse) # 771
length(notTimeCourse) # 1058
# visual check of our RegExp result
# write.table(timeCourse, file = "cage_exp_name_timeCourse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(notTimeCourse, file = "cage_exp_name_notTimeCourse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

removeThose <- grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)
notTimeCourse <- expNames[-grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
dim(cageMatrix)
cageMatrix <- cageMatrix[, -removeThose]
dim(cageMatrix)

which(rowSums(cageMatrix) == 0) %>% length # 219

regionLocations <- regionLocations[which(rowSums(cageMatrix) != 0), ]
regionLocationsGR <- with(regionLocations, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand = strand))
cageMatrix <- cageMatrix[which(rowSums(cageMatrix) != 0), ]

dim(cageMatrix)
regionLocationsGR
length(notTimeCourse)

cellType <- notTimeCourse
cellType <- gsub("_donor[0-9]*", "", cellType)
cellType <- gsub("_biol_rep[0.9]*", "", cellType)
cellType <- gsub("_tech_rep[0.9]*", "", cellType)
cellType <- gsub("_MC-[0.9]*", "", cellType)
cellType <- gsub("_donation[0.9]*", "", cellType)
cellType <- gsub("_ENCODE", "", cellType)
cellType <- gsub("_lineage[0.9]", "", cellType)
cellType <- gsub("_derived[0.9]", "", cellType)
cellType <- gsub("_pluriselect[0.9]", "", cellType)
cellType <- gsub("_control[0.9]", "", cellType)

cellType <- gsub("\\.[[:alnum:]]*\\.[[:alnum:]]*-[[:alnum:]]*$", "", cellType)
cellType <- gsub("SOC-57-02", "", cellType)
cellType <- gsub("SOC-57-02-G", "", cellType)


head(cellType, 20)

anno <- data.frame(
    "name" = notTimeCourse,
    "tissue" = cellType,
    stringsAsFactors = FALSE
)
anno$isCellType <- rep(FALSE, nrow(anno))
anno[grep("_cell_line", anno$tissue), "isCellType"] <- TRUE

write.table(anno,  file = "cage_exp_anno.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") # regExp result checking
write.table(unique(anno$tissue),  file = "cage_exp_anno_small.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") # regExp result checking

t0 <- Sys.time()
correlationMatrix <- cor(cageMatrix)
Sys.time() - t0 # 7 min

fantom5_human_cage <- list(
    "dataMatrix" = cageMatrix,
    "regionMetaData" = regionLocationsGR,
    "correlationMatrix" =  correlationMatrix,
    "annotation" = anno
)

fantom5_human_cage$anno$url <- "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/"
colnames(fantom5_human_cage$dataMatrix) <- NULL
colnames(fantom5_human_cage$correlationMatrix) <- NULL
rownames(fantom5_human_cage$dataMatrix) <- NULL
rownames(fantom5_human_cage$correlationMatrix) <- NULL
rownames(fantom5_human_cage$annotation) <- NULL


save(fantom5_human_cage, file = "../../heatcageseq/data/fantom5_human_cage.RData")







