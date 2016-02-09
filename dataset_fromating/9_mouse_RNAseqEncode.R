#  mouse_encode_rnaseq_tsv.txt
# file generated trhough:
# https://www.encodeproject.org/search/?type=Experiment&assay_term_name=RNA-seq&assembly=mm10&files.file_type=tsv
# 2016-01-22

setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/mEncodeRnaseq")

library(readr)
library(magrittr)
library(parallel)
library(biomaRt)

metadata <- read_tsv("metadata.tsv",
                     col_types = paste0(
                         paste0(rep("?", 11), collapse = ""),
                         "c",
                         paste0(rep("?", 41 - 12), collapse = ""),
                         collapse = "")
)
dim(metadata) # 416
length(list.files()) # 418 / 2 metadata

metadata <- subset(metadata, metadata[, "Output type"] == "gene quantifications")
dim(metadata) # 208
apply(metadata[, c(12:19)],
      2,
      unique
)

metadata[, "Biosample term name"] %>% unique
metadata[, "Biosample type"] %>% unique
metadata[, "File format"] %>% unique

t0 <- Sys.time()
allFiles <- mclapply(metadata[, "File accession"], function(x) read_tsv(paste0(x, ".tsv"), col_types = "cc?????????????"), mc.preschedule = FALSE, mc.cores = 9)
Sys.time() - t0

table(sapply(allFiles, nrow))
# 65268 67472 69690
# 8     8   192

table(sapply(allFiles, ncol))
# 15
# 208

sapply(allFiles, nrow)
head(allFiles[[1]])
allFiles[[1]][40000:40005,]
head(allFiles[[72]]) # 65268 does have ENSM + FPKM
allFiles[[72]][40000:40005,]
head(allFiles[[73]]) # idem for 67472
allFiles[[73]][40000:40005,]

ENSM_case1 <- allFiles[[1]]$gene_id
ENSM_case2 <- allFiles[[72]]$gene_id
ENSM_case3 <- allFiles[[73]]$gene_id

ENSM_case1 <- ENSM_case1[grep("ENSMUS", ENSM_case1, fixed = TRUE)]
ENSM_case2 <- ENSM_case2[grep("ENSMUS", ENSM_case2, fixed = TRUE)]
ENSM_case3 <- ENSM_case3[grep("ENSMUS", ENSM_case3, fixed = TRUE)]

length(ENSM_case1) # 43346
length(ENSM_case2) # 38924
length(ENSM_case3) # 41128

# we will keep only the 192 files with 69690 rows

case1_keep <- which(sapply(allFiles, nrow) == 69690)

case1 <- list(
    metadata = metadata[case1_keep,],
    geneName = allFiles[case1_keep][[1]]$gene_id,
    files = allFiles[case1_keep]
)

case1$geneName %>% length
case1$files[[1]] %>% dim

keepENSMUS <- grep("ENSMUS", case1$geneName, fixed = TRUE)
ENSMUS <- case1$geneName[keepENSMUS]
ENSMUS <- sub("\\.[0-9]*", "", ENSMUS)
case1$geneName <- ENSMUS
case1$files <- lapply(case1$files, function(x) x[keepENSMUS, ])

case1$geneName %>% length
case1$files[[1]] %>% dim

lapply(case1$files[1:5], tail)

case1$fpkmMatrix <- matrix(0, nrow = nrow(case1$files[[1]]), ncol = length(case1$files))
t0 <- Sys.time()
for (i in 1:length(case1$files)) {
    case1$fpkmMatrix[, i] <- case1$files[[i]][, "FPKM"]
}
Sys.time() - t0 # 1 sec !

t0 <- Sys.time()
case1$corMatrix <- cor(case1$fpkmMatrix)
Sys.time() - t0 # 2 sec

dim(case1$corMatrix)

# pdf(file = "../../Rplots/heatmap_case1_encode_rnaseq_muouse.pdf", width = 13.85, height = 13.85)
# heatmap(case1$corMatrix,
#         scale = "none",
#         margins = rep(1, 2),
#         cexRow = 1,
#         cexCol = 1,
#         breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
#         col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
#         useRaster = FALSE
# )
# dev.off()
# system("firefox ../../Rplots/heatmap_case1_encode_rnaseq_muouse.pdf")

unique(case1$metadata[, "Biosample subcellular fraction term name"])
unique(case1$metadata[, "Library made from"])
table(case1$metadata[, "Library made from"])
unique(case1$metadata[, "Library depleted in"])
table(case1$metadata[, "Library depleted in"])

encode_mouse_rnaseq <- list(
    "dataMatrix" = case1$fpkmMatrix,
    "geneName" = case1$geneName,
    "correlationMatrix" = case1$corMatrix,
    "annotation" = case1$metadata[, c("File accession", "Biosample term name", "Biosample type",
                                      "Library depleted in",
                                      "File download URL")]
)

encode_mouse_rnaseq$annotation$name <- paste(encode_mouse_rnaseq$annotation[, "Biosample term name"], encode_mouse_rnaseq$annotation[, "File accession"])
object.size(encode_mouse_rnaseq) # 70 Mo

apply(encode_mouse_rnaseq$annotation[, 2:4],
      2,
      unique
)

nonZero <- which(rowSums(encode_mouse_rnaseq$dataMatrix) != 0)
encode_mouse_rnaseq$dataMatrix <- encode_mouse_rnaseq$dataMatrix[nonZero, ]
encode_mouse_rnaseq$geneName <- encode_mouse_rnaseq$geneName[nonZero]
encode_mouse_rnaseq$correlationMatrix <- cor(encode_mouse_rnaseq$dataMatrix)
encode_mouse_rnaseq$annotation <- encode_mouse_rnaseq$annotation[, c(1, 2, 6, 3, 4, 5)]

lapply(encode_mouse_rnaseq, head)
row.names(encode_mouse_rnaseq$annotation) <- NULL
lapply(encode_mouse_rnaseq, row.names)

colnames(encode_mouse_rnaseq$annotation) <- c("encodeAccession", "tissue", "name", "sampleType", "rnaFraction", "url")

object.size(encode_mouse_rnaseq) # 58 Mo

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
save(encode_mouse_rnaseq, file = "heatrnaseq/data/encode_mouse_rnaseq.RData")



