library(readr)
library(magrittr)
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/bgee_human")
bgee_human <- list(
        GSE30352 = read_tsv("Homo_sapiens_RNA-Seq_read_counts_RPKM_GSE30352.tsv"),
        GSE30611 = read_tsv("Homo_sapiens_RNA-Seq_read_counts_RPKM_GSE30611.tsv"),
        GSE43520 = read_tsv("Homo_sapiens_RNA-Seq_read_counts_RPKM_GSE43520.tsv")
)
bgee_human_big <- do.call(rbind, bgee_human)

nPerGene <- table(bgee_human_big[,"Gene ID"])
as.vector(nPerGene) %>% summary # 77

apply(bgee_human_big, 2, function(x) length(unique(x)))
# Experiment ID             Library ID           Library type
# 3                     77                      2
# Gene ID   Anatomical entity ID Anatomical entity name
# 63677                     22                     22
# Stage ID             Stage name             Read count
# 26                     26                  28711
# RPKM         Detection flag      Detection quality
# 1068217                      2                      1
# State in Bgee
# 4

bgeeExp <- unique(bgee_human_big[, "Library ID"])
names(bgeeExp) <- bgeeExp
bgee_human_bigl <- lapply(bgeeExp, function(x) subset(bgee_human_big, bgee_human_big[, "Library ID"] == x))
lapply(bgee_human_bigl, function(x) tail(x[, c("Gene ID", "RPKM")]))
lapply(bgee_human_bigl, dim)

load("../../heatrnaseq/data/encode_rnaseq.RData")
encode_rnaseq$dataMatrix[1:10, 1:10]
dim(encode_rnaseq$dataMatrix)

bgee_human <- list()
bgee_human$dataMatrix <- matrix(0, nrow = 63677, ncol = 77)
for(i in 1:77) {
    bgee_human$dataMatrix[, i] <- bgee_human_bigl[[i]]$RPKM
}
bgee_human$geneName <- bgee_human_bigl[[1]][, "Gene ID"]
bgee_human$correlationMatrix <- cor(bgee_human$dataMatrix)

bgee_human_annotation <- do.call(rbind, lapply(bgee_human_bigl, function(x) x[1, ]))
bgee_human$annotation <- bgee_human_annotation[, c("Library ID", "Anatomical entity name", "Stage name", "Library type")]
names(bgee_human$annotation) <- c("File accession", "Biosample term name", "Stage name", "Library type")
bgee_human$annotation$fileURL <- paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", bgee_human$annotation[, "File accession"])
bgee_human$annotation[which(bgee_human$annotation[, "Biosample term name"] == "multi-cellular organism"), "Biosample term name"] <- "16 tissues mixture"

names(bgee_human)
lapply(bgee_human, dim)

nonZero <- which(rowSums(bgee_human$dataMatrix) != 0)
bgee_human$dataMatrix <- bgee_human$dataMatrix[nonZero, ]
bgee_human$geneName <- bgee_human$geneName[nonZero]
bgee_human$correlationMatrix <- cor(bgee_human$dataMatrix)
bgee_human$annotation$name <- paste(bgee_human$annotation[, "Biosample term name"], bgee_human$annotation[, "File accession"])
bgee_human$annotation <- bgee_human$annotation[, c(1, 2, 6, 3, 4, 5)]
lapply(bgee_human, head)
row.names(bgee_human$annotation) <- NULL
lapply(bgee_human, row.names)
colnames(bgee_human$annotation) <- c("geoAccession", "tissue", "name", "stage", "libraryType", "url")

object.size(bgee_human) # 37Mo

save(bgee_human, file = "/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatrnaseq/data/bgee_human.RData")


