# 2016/03/21
# download all file from: ftp://data.modencode.org/D.melanogaster/Transcriptional-Factor/ChIP-seq/computed-peaks_gff3/
# we merged gff3 files
# gunzip *.gz
# grep -v '^#' --no-filename *.gff3  > ../mergedPeaks.gff3

# negative start values and 0 start values. :(
# removing them manually, then:

# sortBed -i mergedPeaks.gff3 > mergedPeaks_sorted.gff3
# bedtools merge -i mergedPeaks_sorted.gff3 > mergedPeaks_merged.gff3

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/modEncodeDM/computed-peaks_gff3")

library(dplyr)
library(readr)
library(parallel)
library(GenomicRanges)

mergedPeaks <- read_tsv("../mergedPeaks_merged.gff3", col_names = FALSE, col_types = "cii")
names(mergedPeaks) <- c("chr", "start", "end")
table(mergedPeaks$chr)
validChr <- c("X", "2L", "2R", "3L", "3R", "4", "Y", "M", "2CEN", "3CEN", "XY", "XHet", "YHet", "2LHet", "2RHet", "3LHet", "3RHet", "2", "3")
dim(mergedPeaks)
mergedPeaks <- filter(mergedPeaks, chr %in% validChr)
dim(mergedPeaks)
mergedPeaksGR <- with(mergedPeaks, GRanges(chr, IRanges(start, end)))

myFiles <- list.files()

t0 <- Sys.time()
modData <- mclapply(myFiles, function(x) read_tsv(x, col_names = FALSE, skip = 2, col_types = "c????????"), mc.cores = 32)
Sys.time() - t0 # 1 sec

t0 <- Sys.time()
modDataGR <- lapply(
    modData,
    function(x) with(x, GRanges(X1, IRanges(X4, X5)))
)
Sys.time() - t0 # 1 s

t0 <- Sys.time()
listOverlap <- lapply(
    modDataGR,
    function(x) overlapsAny(mergedPeaksGR, x)
)
Sys.time() - t0 # 3 sec

table(sapply(listOverlap, length)) # ok

dataMtrix <- do.call(cbind, listOverlap)

corMtrix <- cor(dataMtrix)

metadata <- sub(".gff3", "", myFiles)
strsplit(metadata, split = ":", fixed = TRUE) %>% sapply(. , length) %>% table
strsplit(metadata, split = ":", fixed = TRUE) %>% do.call(rbind, .) %>% data.frame(stringsAsFactors = FALSE) -> metadata
metadata <- metadata[, c(1, 2, 7)]
colnames(metadata) <- c("antibody", "developmentalStage", "id")

strsplit(metadata$developmentalStage, split = "#", fixed = TRUE) %>% sapply(. , length) %>% table

tempChr <- strsplit(metadata$developmentalStage, split = "#", fixed = TRUE)

library(stringr)

myMD <- data.frame(
    cellLine = str_match(metadata$developmentalStage, ".*Cell-Line=(.*)")[, 2],
    devStage = str_match(metadata$developmentalStage, ".*Developmental-Stage=(.*)")[, 2],
    strain = str_match(metadata$developmentalStage, ".*Strain=(.*)")[, 2],
    stringsAsFactors = FALSE
)
myMD <- lapply(myMD, function(x) gsub("#.*", "", x)) %>% data.frame(stringsAsFactors = FALSE)

metadata <- cbind(metadata[, c("antibody", "id")], myMD)

filter(metadata, antibody == "eGFP")

metadata[!is.na(metadata$cellLine), "devStage"] <- metadata[!is.na(metadata$cellLine), "cellLine"]
# http://stackoverflow.com/questions/13673894/suppress-nas-in-paste
paste3 <- function(... ,sep=" ") {
    L <- list(...)
    L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
    ret <- gsub(paste0("(^",sep,"|",sep,"$)"),"",
                gsub(paste0(sep,sep),sep,
                    do.call(paste,c(L,list(sep=sep)))))
    is.na(ret) <- ret==""
    ret
}

metadata$name <- paste3(metadata$antibody, metadata$devStage, metadata$strain, metadata$id)
length(metadata$name)
length(unique(metadata$name))
metadata$url <- paste0(
    "http://data.modencode.org/cgi-bin/findFiles.cgi?download=",
    sub("modENCODE_", "", metadata$id, fixed = TRUE)
)

metadata <- metadata[, c(1:2, 4:7)]
rownames(metadata) <- NULL

modEncodeD_ChIPseq <- list(
    dataMatrix = dataMtrix,
    regionMetaData = mergedPeaksGR,
    correlationMatrix = cor(dataMtrix),
    annotation = metadata
)

modEncodeD_ChIPseq$annotation$strain[is.na(modEncodeD_ChIPseq$annotation$strain)] <- "unspecify"

save(modEncodeD_ChIPseq, file = "../../../heatchipseq/data/modEncodeD_ChIPseq.RData")

