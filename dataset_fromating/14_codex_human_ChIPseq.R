# 2016/02/04
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/Codex_hgxx/")
library(magrittr)
library(readr)
metadata <- read_tsv("CODEX_hgxx_metadata.txt", col_names = FALSE) # it is hg19
metadata <- metadata[-nrow(metadata),] # error with excel, last line is empty
command <- paste(
    "wget",
    metadata[, 6]
)
# write.table(command, file = "download_codex_hgxx.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)
colnames(metadata) <- c("cell_type", "cell_subtype", "tf", "gse", "gsm", "url")
metadata$filename <- strsplit(metadata$url, split = "/") %>% lapply(function(x) tail(x, 1)) %>% unlist

filesList <- list.files()
summary(filesList %in% metadata$filename)
filesList[!(filesList %in% metadata$filename)] # \o/
metadata$filename[!(metadata$filename %in% filesList)] # /o\
# "GSM1291197_TAL1.bed" "GSM1291199_TAL1.bed" "GSM1085400_TAL1.bed"
command[grep("TAL1", command, fixed = TRUE)]
command[grep("GSM1291197_TAL1.bed", command, fixed = TRUE)] # file not found
command[grep("GSM1291199_TAL1.bed", command, fixed = TRUE)] # file not found
command[grep("GSM1085400_TAL1.bed", command, fixed = TRUE)] # file not found
# to: bg200@cam.ac.uk.
# Dear HSCL team,
#
# I am a Postdoc at the Rolsin Institute (Edinburgh) in Anagha Joshi's group.
# I am trying to download several peak files from the CODEX human database.
# I have been mostly successful, but three files seems to escape my reach: "Files not found"
# http://codex.stemcells.cam.ac.uk/data/bed/hg19/GSM1291197_TAL1.bed
# http://codex.stemcells.cam.ac.uk/data/bed/hg19/GSM1291199_TAL1.bed
# http://codex.stemcells.cam.ac.uk/data/bed/hg19/GSM1085400_TAL1.bed
#
# Could you please provide me the correct link for those file? And maybe correct the links on the codex table?
#
# Kind regards,
#
# Guillaume
#
# reply: this were remanant of re-process files, ie dupplicate entries
#
# in the mean time:
metadata <- subset(metadata, !(metadata$filename %in% c("GSM1291197_TAL1.bed", "GSM1291199_TAL1.bed", "GSM1085400_TAL1.bed")))

t0 <- Sys.time()
peakList <- lapply(metadata$filename, function(x) read_tsv(x, col_names = FALSE))
Sys.time() - t0 # 9s

metadata$name <- sub(".bed", "", metadata$filename, fixed = TRUE)
names(peakList) <- metadata$name
peakList <- lapply(peakList, function(x) x[, 1:3])
peakList <- lapply(1:length(peakList), function(x) {
    peakList[[x]]$exp <- names(peakList)[x]
    return(peakList[[x]])
})

a <- Sys.time()
bigBed <- do.call(rbind, peakList)
Sys.time() - a # 12 s
dim(bigBed)
#[1] 3591101       4

# write.table(bigBed, file = "allCodexHg19.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
# then bedtools :
# sortBed -i allCodexHg19.bed > allCodexHg19_sorted.bed
# bedtools merge -c 4 -o distinct -i allCodexHg19_sorted.bed > allCodexHg19_merged.bed

#########
# new R session
#########
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/Codex_hgxx/")
library(magrittr)
library(readr)
library(parallel)
library(GenomicRanges)

codexMerged <- read_tsv("allCodexHg19_merged.bed", col_names = FALSE)
tail(codexMerged[, 1]) # hrM, seriously CODEX? -> might have been corrected 2016/02/04 after email contact
myChr <- c(paste0("chr", 1:22), "chrX", "chrY") # we discard chrM, should we discard chrY to? and chrX?
codexMerged <- subset(codexMerged, codexMerged[, 1] %in% myChr)
unique(codexMerged[, 1])
dim(codexMerged)
# [1] 702011      4
metadata <- read_tsv("CODEX_hgxx_metadata.txt", col_names = FALSE)
metadata <- metadata[-nrow(metadata),]
colnames(metadata) <- c("cellType", "cellSubtype", "tf", "gse", "gsm", "url")
metadata$filename <- strsplit(metadata$url, split = "/") %>% lapply(function(x) tail(x, 1)) %>% unlist
metadata$name <- sub(".bed", "", metadata$filename, fixed = TRUE)
metadata <- subset(metadata, !(metadata$filename %in% c("GSM1291197_TAL1.bed", "GSM1291199_TAL1.bed", "GSM1085400_TAL1.bed")))

templateVector <- vector(mode = "logical", length = nrow(metadata))
names(templateVector) <- metadata$name

createBoolVectorForLine <- function(line) {
    templateVector[strsplit(codexMerged[line, "X4"], split=",")[[1]]] <- TRUE
    return(templateVector)
}
t0 <- Sys.time()
dataMatrix <- do.call(rbind, mclapply(1:nrow(codexMerged), createBoolVectorForLine, mc.cores = 32))
Sys.time() - t0 # 1 minutes, mostly slow rbind

dim(dataMatrix)
#[1]   702011    238
summary(as.vector(dataMatrix))
# Mode     FALSE      TRUE      NA's
# logical 163589570   3489048         0
which(!apply(dataMatrix, 2, any))
# named integer(0)
regionMetaData <-  GRanges(codexMerged[,1], IRanges(codexMerged[,2], codexMerged[,3]))
regionMetaData

# precomputation off correlations
t0 <- Sys.time()
correlationMatrix <- cor(dataMatrix)
Sys.time() - t0 # 1.2 minutes.

dim(correlationMatrix)
# [1] 238 238
summary(as.vector(correlationMatrix))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.06998  0.01790  0.06116  0.11810  0.15800  1.00000

codex_human_chip <- list("dataMatrix" = dataMatrix,
              "regionMetaData" = regionMetaData,
              "correlationMatrix" = correlationMatrix,
              "annotation" = metadata)
codex_human_chip$metadata <- codex_human_chip$metadata[, c(8, 1:6)]

unique(codex_human_chip$annotation$cellType)
# [1] "Acute Myeloid Leukemia"
# [2] "acute myeloid leukemia"
# [3] "Acute myeloid leukemia (AML)"
# [4] "B-cells"
# [5] "B-Cells"
# [6] "Cervical cancer cells"
# [7] "Chronic myelogenous leukemia cells"
# [8] "Dendritic"
# [9] "Embryonic kidney cells"
# [10] "Embryonic Stem Cell"
# [11] "Embryonic stem cell"
# [12] "Endoderm"
# [13] "Endometrial Stromal Cell"
# [14] "endothelial progenitor cells"
# [15] "Erythroblast"
# [16] "Erythrocytic leukaemia cells"
# [17] "Erythroid Leukaemia"
# [18] "Erythroid Progenitor"
# [19] "Hematopoietic Stem and Progenitor Cells"
# [20] "Jurkat"
# [21] "Jurkat cell line"
# [22] "Leukaemia"
# [23] "Leukaemia cell"
# [24] "Leukemia"
# [25] "Leukemia cell"
# [26] "Lymphoma cell"
# [27] "Lymphoma cell line"
# [28] "macrophage"
# [29] "Macrophage"
# [30] "Megakaryocyte"
# [31] "Mesenchymal stem cell"
# [32] "Mesendoderm cell"
# [33] "Monocyte"
# [34] "Multiple myeloma cell"
# [35] "Myeloid Leukemia"
# [36] "Neuronal Progenitor cell"
# [37] "Raji cells"
# [38] "T-Cell"
# [39] "T-cell acute lymphoblastic leukemia cells"
# [40] "T-lymphoblastic leukemia (T-ALL)"
# [41] "T_ALL"
# [42] "T_ALL primagraft"
# [43] "Trophectoderm cell"
# [44] "Umbilical Cord Blood Stem and Progenitor Cells"

ct <- factor(codex_human_chip$annotation$cellType)
lct <- levels(ct)

lct[lct == "acute myeloid leukemia"] <- "Acute Myeloid Leukemia"
lct[lct == "Acute myeloid leukemia (AML)"] <- "Acute Myeloid Leukemia"
lct[lct == "B-Cells"] <- "B-cells"
lct[lct == "Embryonic stem cell"] <- "Embryonic Stem Cell"
lct[lct == "Jurkat"] <- "Jurkat cell line"
lct[lct == "Leukemia"] <- "Leukaemia"
lct[lct == "Leukemia cell"] <- "Leukaemia cell"
lct[lct == "Macrophage"] <- "macrophage"
lct[lct == "T_ALL"] <- "T-lymphoblastic leukemia (T-ALL)"

levels(ct) <- lct
codex_human_chip$annotation$cellType <- as.character(ct)
unique(codex_human_chip$annotation$cellType)

unique(codex_human_chip$annotation$cellSubtype)
ct <- factor(codex_human_chip$annotation$cellSubtype)
lct <- levels(ct)
lct[lct == "[CL] Kasumi-1 acute myeloid leukaemia +mismatch siRNA."] <- "[CL] Kasumi-1 acute myeloid leukaemia +mismatch siRNA"
lct[lct == "[PC] Peripherial blood derived proerythroblast"] <- "[PC] Peripheral blood derived proerythroblast"
levels(ct) <- lct
codex_human_chip$annotation$cellSubtype <- as.character(ct)
# maybe jurkat also
codex_human_chip$annotation <- codex_human_chip$annotation[, c(1:6, 8)]
codex_human_chip$annotation$name <- paste(codex_human_chip$annotation$name, codex_human_chip$annotation$cellSubtype)
object.size(codex_human_chip) # 67Mo
save(codex_human_chip, file = "../../heatchipseq/data/codex_human_chip.RData")


