setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/EncodeRNAseq")

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
dim(metadata) # 710
length(list.files()) # 712 / 2 metadata

metadata <- subset(metadata, metadata[, "Output type"] == "gene quantifications")
dim(metadata) # 390
apply(metadata[, c(12:19)],
      2,
      unique
      )

metadata[, "Biosample term name"] %>% unique
metadata[, "Biosample type"] %>% unique
metadata[, "File format"] %>% unique

t0 <- Sys.time()
allFiles <- mclapply(metadata[, "File accession"], function(x) read_tsv(paste0(x, ".tsv")), mc.preschedule = FALSE, mc.cores = 32)
Sys.time() - t0 # 43 seconds

table(sapply(allFiles, nrow))
# 57904 58447 58540
# 50    20   320

table(sapply(allFiles, ncol))
# 4  15
# 70 320

###
# case 1
###
case1_keep <- which(sapply(allFiles, nrow) == 58540)
case1 <- list(
    metadata = metadata[case1_keep,],
    geneName = allFiles[case1_keep][[1]]$gene_id,
    files = allFiles[case1_keep]
)

length(unique(case1$files[[1]]$gene_id)) # 58540
length(grep("ENSG", case1$files[[1]]$gene_id, fixed = TRUE)) # 57820

case1$fpkmMatrix <- matrix(0, nrow = nrow(case1$files[[1]]), ncol = length(case1$files))

t0 <- Sys.time()
for (i in 1:length(case1$files)) {
    case1$fpkmMatrix[, i] <- case1$files[[i]][, "FPKM"]
}
Sys.time() - t0 # 1 sec !

ENSG <- grep("ENSG", case1$files[[1]]$gene_id, fixed = TRUE)
case1 <- list(
    metadata = case1$metadata,
    geneName = case1$geneName[ENSG],
    fpkmMatrix = case1$fpkmMatrix[ENSG, ]
)

t0 <- Sys.time()
case1$corMatrix <- cor(case1$fpkmMatrix)
Sys.time() - t0 # 5 sec

dim(case1$corMatrix)

# pdf(file = "../../Rplots/heatmap_case1_encode_rnaseq.pdf", width = 13.85, height = 13.85)
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
# Sys.time() - a # 3 minutes

case1$geneName <- sub("\\.[0-9]*", "", case1$geneName)
    
# mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
# t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

###
# case 2
###
case2_keep <- which(sapply(allFiles, nrow) == 57904)
allFiles[case2_keep][[1]] %>% head(50) # no FPKM values...

###
# case 3
###
case3_keep <- which(sapply(allFiles, nrow) == 58447)
allFiles[case3_keep][[1]] %>% head(50) # no FPKM values...

# case 1 metadata
unique(case1$metadata[, "Biosample subcellular fraction term name"])
unique(case1$metadata[, "Library made from"])
unique(case1$metadata[, "Library depleted in"])

case1$metadata[, "Biosample subcellular fraction term name"][is.na(case1$metadata[, "Biosample subcellular fraction term name"])] <- "all"
case1$metadata[, "Library depleted in"][is.na(case1$metadata[, "Library depleted in"])] <- ""
case1$metadata[, "Library depleted in"][case1$metadata[, "Library depleted in"] == "rRNA"] <- ", rRNA_depletion"
case1$metadata[, "Library depleted in"][case1$metadata[, "Library depleted in"] == "rRNA"] <- ", polyA_mRNA_and_rRNA_depletion"

unique(case1$metadata[, "Biosample subcellular fraction term name"])
unique(case1$metadata[, "Library made from"])
unique(case1$metadata[, "Library depleted in"])

case1$metadata <- as.data.frame(case1$metadata, stringsAsFactors = FALSE)

case1$metadata$rnaFraction <- paste0(
    case1$metadata[, "Biosample subcellular fraction term name"],
    ", ",
    case1$metadata[, "Library made from"],
    case1$metadata[, "Library depleted in"]
)
unique(case1$metadata$rnaFraction)

encode_rnaseq <- list(
    "dataMatrix" = case1$fpkmMatrix,
    "geneName" = case1$geneName,
    "correlationMatrix" = case1$corMatrix,
    "annotation" = case1$metadata[, c("File accession", "Biosample term name", "Biosample type",
                                      "rnaFraction", 
                                      "File download URL")]
)

apply(encode_rnaseq$annotation[, 2:5],
      2,
      unique
)

nonZero <- which(rowSums(encode_rnaseq$dataMatrix) != 0)
encode_rnaseq$dataMatrix <- encode_rnaseq$dataMatrix[nonZero, ]
encode_rnaseq$geneName <- encode_rnaseq$geneName[nonZero]
encode_rnaseq$correlationMatrix <- cor(encode_rnaseq$dataMatrix)
encode_rnaseq$annotation$name <- paste(encode_rnaseq$annotation[, "Biosample term name"], encode_rnaseq$annotation[, "File accession"])
encode_rnaseq$annotation <- encode_rnaseq$annotation[, c(1, 2, 6, 3, 4, 5)] 

lapply(encode_rnaseq, head)
row.names(encode_rnaseq$annotation) <- NULL
lapply(encode_rnaseq, row.names)

object.size(encode_rnaseq) # 125 Mo

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
save(encode_rnaseq, file = "heatrnaseq/data/encode_rnaseq.RData")


