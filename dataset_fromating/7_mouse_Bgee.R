# shell script
# 2016/01/06
wget ftp://ftp.bgee.org/current/download/processed_expr_values/rna_seq/Mus_musculus/Mus_musculus_RNA-Seq_read_counts_RPKM.zip
unzip Mus_musculus_RNA-Seq_read_counts_RPKM.zip
rm Mus_musculus_RNA-Seq_read_counts_RPKM.zip
#unzip all zip
rm *.zip

# R script
library(readr)
library(magrittr)
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/bgee_mouse")
bgee_mouse_files <- list.files()
bgee_mouse <- lapply(bgee_mouse_files, read_tsv)
length(bgee_mouse)
lapply(bgee_mouse, head)
lapply(bgee_mouse, dim)

bgee_mouse_big <- do.call(rbind, bgee_mouse)

nPerGene <- table(bgee_mouse_big[, "Gene ID"])
as.vector(nPerGene) %>% summary # 109

apply(bgee_mouse_big, 2, function(x) length(unique(x)))
# Experiment ID             Library ID           Library type
# 7                    109                      2
# Gene ID   Anatomical entity ID Anatomical entity name
# 39179                     20                     20
# Stage ID             Stage name             Read count
# 6                      6                  33611
# RPKM         Detection flag      Detection quality
# 1232604                      2                      1
# State in Bgee
# 4

bgeeExp <- unique(bgee_mouse_big[, "Library ID"])
names(bgeeExp) <- bgeeExp
bgee_mouse_bigl <- lapply(bgeeExp, function(x) subset(bgee_mouse_big, bgee_mouse_big[, "Library ID"] == x))
lapply(bgee_mouse_bigl, function(x) tail(x[, c("Gene ID", "RPKM")]))
lapply(bgee_mouse_bigl, dim)

load("../../heatrnaseq/data/bgee_human.RData")

bgee_mouse <- list()
bgee_mouse$dataMatrix <- matrix(0, nrow = 39179, ncol = 109)
for(i in 1:109) {
    bgee_mouse$dataMatrix[, i] <- bgee_mouse_bigl[[i]]$RPKM
}
bgee_mouse$geneName <- bgee_mouse_bigl[[1]][, "Gene ID"]
bgee_mouse$correlationMatrix <- cor(bgee_mouse$dataMatrix)

bgee_mouse_annotation <- do.call(rbind, lapply(bgee_mouse_bigl, function(x) x[1, ]))
bgee_mouse$annotation <- bgee_mouse_annotation[, c("Library ID", "Anatomical entity name", "Stage name", "Library type")]
names(bgee_mouse$annotation) <- c("File accession", "Biosample term name", "Stage name", "Library type")
bgee_mouse$annotation$fileURL <- paste0("http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", bgee_mouse$annotation[, "File accession"])

names(bgee_mouse)
lapply(bgee_mouse, dim)

nonZero <- which(rowSums(bgee_mouse$dataMatrix) != 0)
bgee_mouse$dataMatrix <- bgee_mouse$dataMatrix[nonZero, ]
bgee_mouse$geneName <- bgee_mouse$geneName[nonZero]
bgee_mouse$correlationMatrix <- cor(bgee_mouse$dataMatrix)
bgee_mouse$annotation$name <- paste(bgee_mouse$annotation[, "Biosample term name"], bgee_mouse$annotation[, "File accession"])
bgee_mouse$annotation <- bgee_mouse$annotation[, c(1, 2, 6, 3, 4, 5)]
lapply(bgee_mouse, head)
row.names(bgee_mouse$annotation) <- NULL
lapply(bgee_mouse, row.names)
colnames(bgee_mouse$annotation) <-  c("geoAccession", "tissue", "name", "stage", "libraryType", "url")

object.size(bgee_mouse) # 36Mo

# log transform
# log10(FPKM) + 1
t0 <- Sys.time()
bgee_mouse$dataMatrix <- log10(bgee_mouse$dataMatrix + 1)
bgee_mouse$correlationMatrix <- cor(bgee_mouse$dataMatrix)
Sys.time() - t0 #

save(bgee_mouse, file = "/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatrnaseq/data/bgee_mouse.RData")

# various old tests
save(bgee_mouse, file = "/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/application4/data/bgee_mouse.RData")
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/")

# spearman of FPKM
bgee_mouse_spearman <- bgee_mouse
t0 <- Sys.time()
bgee_mouse_spearman$correlationMatrix <- cor(bgee_mouse_spearman$dataMatrix, method = "spearman")
Sys.time() - t0 # 1 sec
save(bgee_mouse_spearman, file = "application4/data/bgee_mouse_spearman.RData")

# log10(FPKM) + 1
bgee_mouse_log10 <- bgee_mouse
t0 <- Sys.time()
bgee_mouse_log10$dataMatrix <- log10(bgee_mouse$dataMatrix + 1)
bgee_mouse_log10$correlationMatrix <- cor(bgee_mouse_log10$dataMatrix)
Sys.time() - t0 #
save(bgee_mouse_log10, file = "application4/data/bgee_mouse_log10.RData")

# asinh
bgee_mouse_asinh <- bgee_mouse
t0 <- Sys.time()
bgee_mouse_asinh$dataMatrix <- asinh(bgee_mouse$dataMatrix)
bgee_mouse_asinh$correlationMatrix <- cor(bgee_mouse_asinh$dataMatrix)
Sys.time() - t0 #
save(bgee_mouse_asinh, file = "application4/data/bgee_mouse_asinh.RData")

# invnorm
bgee_mouse_invnorm <- bgee_mouse
t0 <- Sys.time()
bgee_mouse_invnorm$dataMatrix <- 1/(bgee_mouse$dataMatrix + 1)
bgee_mouse_invnorm$correlationMatrix <- cor(bgee_mouse_invnorm$dataMatrix)
Sys.time() - t0 #
save(bgee_mouse_invnorm, file = "application4/data/bgee_mouse_invnorm.RData")

# log2(FPKM + 1)
bgee_mouse_log2 <- bgee_mouse
t0 <- Sys.time()
bgee_mouse_log2$dataMatrix <- log2(bgee_mouse$dataMatrix + 1)
bgee_mouse_log2$correlationMatrix <- cor(bgee_mouse_log2$dataMatrix)
Sys.time() - t0 #
save(bgee_mouse_log2, file = "application4/data/bgee_mouse_log2.RData")

