###
# shell scripts
### --------------------------

# downloaded from http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/
# 2015-11-04
# file 	mm9.cage_peak_phase1and2combined_tpm.osc.txt.gz

# file is challenging to load in R. We split headers in two other files, and keep only the matrix
grep -E '^#' mm9.cage_peak_phase1and2combined_tpm.osc.txt | wc -l # 1075
grep -E '^#' mm9.cage_peak_phase1and2combined_tpm.osc.txt > mm9.cage_peak_phase1and2combined_tpm.osc.txt.header_1.txt
wc -l mm9.cage_peak_phase1and2combined_tpm.osc.txt # 160044
# 160044 - 1075
tail -158969 mm9.cage_peak_phase1and2combined_tpm.osc.txt > mm9.cage_peak_phase1and2combined_tpm.osc.reddata.txt
head mm9.cage_peak_phase1and2combined_tpm.osc.reddata.txt | awk '{ print $1, $2,  $3,  $4;}'

head -3 mm9.cage_peak_phase1and2combined_tpm.osc.reddata.txt > mm9.cage_peak_phase1and2combined_tpm.osc.header_2.txt
tail -158966 mm9.cage_peak_phase1and2combined_tpm.osc.reddata.txt > mm9.cage_peak_phase1and2combined_tpm.osc.data.txt
head mm9.cage_peak_phase1and2combined_tpm.osc.data.txt | awk '{ print $1, $2,  $3,  $4;}'

###
# R
###
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/CAGE_mm9/")

library(readr)
library(dplyr)
library(GenomicRanges)
library(parallel)
library(ggplot2)

a <- Sys.time()
cageData <- list(
    header1 = read_tsv("mm9.cage_peak_phase1and2combined_tpm.osc.txt.header_1.txt", col_names = FALSE),
    header2 = read_tsv("mm9.cage_peak_phase1and2combined_tpm.osc.header_2.txt", col_names = FALSE, col_types = do.call(paste0, as.list(rep("c", 1074)))),
    data = read_tsv("mm9.cage_peak_phase1and2combined_tpm.osc.data.txt", col_names = FALSE)
)
Sys.time() - a # 1 minute


# cageData <- lapply(cageData, as.data.frame)
lapply(cageData, dim)
list(
    header1 = head(cageData$header1),
    header2 = cageData$header2[, 1:6],
    data = cageData$data[1:4, 1:6]
)

regionLocations <- cageData$data[, 1] %>% unlist
length(regionLocations)
a <- Sys.time()
regionLocations %>% strsplit(split = ":|\\.\\.|,", fixed = FALSE) %>% do.call(rbind, .) -> regionLocations
Sys.time() - a # 6 seconds
dim(regionLocations)
colnames(regionLocations) <- c("chr", "start", "end", "strand")
regionLocations <- data.frame(regionLocations, stringsAsFactors = FALSE)
regionLocations$name <- rep("no_name", length.out = nrow(regionLocations))
regionLocations$score <- rep(0, length.out = nrow(regionLocations))
regionLocations <- regionLocations[, c(1,2,3,5,6,4)]

# write.table(regionLocations, file = "all_cage_regions.bed", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# shell, bedtools commands. Checking if all regions are unique. They are.
# sortBed -i all_cage_regions.bed > all_cage_regions.s.bed
# bedtools merge -s -c 4 -o count -d -1 -i all_cage_regions.s.bed > all_cage_regions.merged.bed
# wc -l all_cage_regions.bed # 158966
# wc -l all_cage_regions.s.bed # 158966
# wc -l all_cage_regions.merged.bed # 158966 youpi !

cageMatrix <- as.matrix(cageData$data[, -1])

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

timeCourse <- expNames[grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
notTimeCourse <- expNames[-grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
length(timeCourse) # 583
length(notTimeCourse) # 490
# visual check of our RegExp result
# write.table(timeCourse, file = "cage_exp_name_timeCourse_mouse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(notTimeCourse, file = "cage_exp_name_notTimeCourse_mouse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
# meh... that will do for now. But A bit 2 many are in timeCourse.

removeThose <- grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)
notTimeCourse <- expNames[-grep("[0-9]hr|day[0-9][0-9]", expNames, fixed = FALSE)]
dim(cageMatrix)
cageMatrix <- cageMatrix[, -removeThose]
dim(cageMatrix)

which(rowSums(cageMatrix) == 0) %>% length # 569

regionLocations <- regionLocations[which(rowSums(cageMatrix) != 0), ]
regionLocationsGR <- with(regionLocations, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand = strand))
cageMatrix <- cageMatrix[which(rowSums(cageMatrix) != 0), ]

dim(cageMatrix)
regionLocationsGR
length(notTimeCourse)

cellType <- notTimeCourse
cellType <- gsub("_donor[0-9]*", "", cellType)
cellType <- gsub("_biol_rep[0-9]*", "", cellType)
cellType <- gsub("_tech_rep[0-9]*", "", cellType)
cellType <- gsub("_MC-[0-9]*", "", cellType)
cellType <- gsub("_donation[0-9]*", "", cellType)
cellType <- gsub("_ENCODE", "", cellType)
cellType <- gsub("_lineage[0-9]", "", cellType)
cellType <- gsub("_derived[0-9]", "", cellType)
cellType <- gsub("_pluriselect[0-9]", "", cellType)
cellType <- gsub("_control[0-9]", "", cellType)
cellType <- gsub("_pool[0-9]*", "", cellType)

cellType <- gsub("\\.[[:alnum:]]*\\.[[:alnum:]]*-[[:alnum:]]*$", "", cellType)
cellType <- gsub("SOC-57-02", "", cellType)
cellType <- gsub("SOC-57-02-G", "", cellType)
cellType <- gsub("Mouse_", "", cellType)

anno <- data.frame(
    "name" = notTimeCourse,
    "tissue" = cellType,
    stringsAsFactors = FALSE
)

head(anno, 20)

# write.table(anno,  file = "cage_exp_anno_mouse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") # regExp result checking
# write.table(unique(anno$tissue),  file = "cage_exp_anno_small_mouse.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t") # regExp result checking

t0 <- Sys.time()
correlationMatrix <- cor(cageMatrix)
Sys.time() - t0 # 30 s

fantom5_mouse_cage <- list(
    "dataMatrix" = cageMatrix,
    "regionMetaData" = regionLocationsGR,
    "correlationMatrix" =  correlationMatrix,
    "annotation" = anno
)

fantom5_mouse_cage$annotation$url <- "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/"
colnames(fantom5_mouse_cage$dataMatrix) <- NULL
colnames(fantom5_mouse_cage$correlationMatrix) <- NULL
rownames(fantom5_mouse_cage$dataMatrix) <- NULL
rownames(fantom5_mouse_cage$correlationMatrix) <- NULL
rownames(fantom5_mouse_cage$annotation) <- NULL

save(fantom5_mouse_cage, file = "../../heatcageseq/data/fantom5_mouse_cage.RData")
