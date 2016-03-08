setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
library(GenomicRanges)
library(magrittr)
library(readr)

load("heatcageseq/data/fantom5_human_cage.RData")

userCageFile <- read_tsv("data/CAGE_hg38/Hep-2%20cells%20mock%20treated%2c%20biol_rep1.CNhs13479.11898-125E8.hg19.ctss.bed", col_names = FALSE)
head(userCageFile, 25)
summary(userCageFile$X5)
summary(fantom5_human_cage$dataMatrix[, 1][fantom5_human_cage$dataMatrix[, 1] != 0])
colnames(userCageFile) <- c("chr", "start", "end", "name", "score", "strand")

sum(userCageFile$score)
userCageFile$score <- userCageFile$score * 1E6 /sum(userCageFile$score)
userCageFile <- userCageFile[userCageFile$score >= 1,]

userCageFileGR <- with(userCageFile, GRanges(chr, IRanges(start, end), strand = strand, score = score))

t0 <- Sys.time()
myO <- findOverlaps(fantom5_human_cage$regionMetaData, userCageFileGR, select = "arbitrary")
Sys.time() - t0 # fast

system.time(userCuratedValues <- numeric(nrow(fantom5_human_cage$dataMatrix)))
t0 <- Sys.time()
userCuratedValues[!is.na(myO)] <- userCageFile[myO[!is.na(myO)], "score"]
Sys.time() - t0 # fast

summary(userCageFile$score)
summary(userCuratedValues[userCuratedValues!=0])
summary(userCuratedValues)

t0 <- Sys.time()
userCor <- cor(userCuratedValues, fantom5_human_cage$dataMatrix) %>% as.vector
Sys.time() - t0 # 3 sec

summary(userCor)

head(fantom5_human_cage$annotation[order(userCor, decreasing = T),])


# exemple file:

myData <- read_tsv("data/CAGE_hg38/hg19.cage_peak_phase1and2combined_tpm.osc.data.txt", col_names = FALSE)

regionLocations <- myData[, 1]
length(regionLocations)
a <- Sys.time()
regionLocations %>% strsplit(split = ":|\\.\\.|,", fixed = FALSE) %>% do.call(rbind, .) -> regionLocations
Sys.time() - a # 7 seconds

exemple1 <- data.frame(
    chr = regionLocations[, 1],
    start = as.character(regionLocations[, 2]),
    end = as.character(regionLocations[, 3]),
    name = "ex1",
    score = myData[, 2],
    strand = regionLocations[, 4],
    stringsAsFactors = FALSE
)
exemple2 <- data.frame(
    chr = regionLocations[, 1],
    start = as.character(regionLocations[, 2]),
    end = as.character(regionLocations[, 3]),
    name = "ex2",
    score = myData[, 78],
    strand = regionLocations[, 4],
    stringsAsFactors = FALSE
)
exemple3 <- data.frame(
    chr = regionLocations[, 1],
    start = as.character(regionLocations[, 2]),
    end = as.character(regionLocations[, 3]),
    name = "ex1",
    score = myData[, 150],
    strand = regionLocations[, 4],
    stringsAsFactors = FALSE
)
exemple1 <- subset(exemple1, exemple1$score != 0)
exemple2 <- subset(exemple2, exemple2$score != 0)
exemple3 <- subset(exemple3, exemple2$score != 0)

write.table(exemple1, file = "exemple1.cage.hg19.bed",  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(exemple2, file = "exemple2.cage.hg19.bed",  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(exemple3, file = "exemple3.cage.hg19.bed",  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


exemple4 <- data.frame(
    chr = regionLocations[, 1],
    start = as.character(regionLocations[, 2]),
    end = as.character(regionLocations[, 3]),
    name = "ex4",
    score = myData[, 320],
    strand = regionLocations[, 4],
    stringsAsFactors = FALSE
)
exemple4 <- subset(exemple4, exemple4$score != 0)
write.table(exemple4, file = "exemple4.cage.hg19.bed",  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

