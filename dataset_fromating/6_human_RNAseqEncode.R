# 2016-05-09
# from encode DCC website:
# https://www.encodeproject.org/batch_download/type%3DExperiment%26assay_slims%3DTranscription%26assay_title%3DpolyA%2BmRNA%2BRNA-seq%26assembly%3Dhg19%26assay_title%3DRNA-seq%26assay_title%3DshRNA%2BRNA-seq%26assay_title%3DpolyA%2Bdepleted%2BRNA-seq%26files.file_type%3Dtsv
# clic on link, no wget (I think?)
# then xargs -n 1 curl -O -L < files.txt
# as ENCODE suggest
# 4496 files!? 2 files / experiment + biological repplicates...
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/EncodeRNAseq/new/")

library(readr)
library(dplyr)
library(parallel)
library(biomaRt)

metadata <- read_tsv("metadata.tsv",
                     col_types = paste0(
                         paste0(rep("?", 11), collapse = ""),
                         "c",
                         paste0(rep("?", 41 - 12), collapse = ""),
                         collapse = "")
                     )
dim(metadata) #  4496   41
length(list.files()) # 4498

colnames(metadata) <- gsub(" ", "", colnames(metadata))
unique(metadata$Fileformat) # cool
unique(metadata$Outputtype)
fmetadata <- filter(metadata, Outputtype == "gene quantifications")
unique(fmetadata$Assay) # cool
table(fmetadata$Assay) # divide dataset? probably...
# RNA-seq shRNA knockdown followed by RNA-seq
# 606                                1641
fmetadata <- filter(fmetadata, Assay == "RNA-seq")
unique(fmetadata$Biosampletype) # ok
unique(fmetadata$Biosampleorganism) # ok
table(fmetadata$Biosampletreatments) # 4 sampls treated with dimethyl sulfoxide
table(fmetadata$Biosamplesubcellularfractiontermname)
unique(fmetadata$Biosamplephase) # good
unique(fmetadata$Biosamplesynchronizationstage) # good
table(fmetadata$Librarymadefrom)
# polyadenylated mRNA                 RNA
# 106                 500
table(fmetadata$Librarydepletedin)
# rRNA rRNA, polyadenylated mRNA
# 526                        40
table(fmetadata$Project) # ok
table(fmetadata$Assembly)
# GRCh38   hg19 # sight...
# 302    304

hg19 <- filter(fmetadata, Assembly == "hg19") %>% arrange(Experimentaccession)
hg38 <- filter(fmetadata, Assembly == "GRCh38") %>% arrange(Experimentaccession)

union(hg19$Experimentaccession, hg38$Experimentaccession) %>% length # 155
intersect(hg19$Experimentaccession, hg38$Experimentaccession) %>% length # 154
setdiff(hg19$Experimentaccession, hg38$Experimentaccession) # ENCSR000AAM

filter(fmetadata, Experimentaccession == "ENCSR000AAM") %>% glimpse # why don't you like it ENCODE?

dim(fmetadata)
paste0(fmetadata$Fileaccession, ".tsv") %in% list.files() %>% which %>% length # 606 \o/

t0 <- Sys.time()
allFiles <- mclapply(hg38$Fileaccession, function(x) read_tsv(paste0(x, ".tsv")), mc.preschedule = FALSE, mc.cores = 32)
Sys.time() - t0 # 1 min
names(allFiles) <- hg38$Fileaccession

table(sapply(allFiles, nrow))
# 61471
# 302

table(sapply(allFiles, ncol))
# 15
# 302

fpkmMatrix <- matrix(0, nrow = nrow(allFiles[[1]]), ncol = length(allFiles))

t0 <- Sys.time()
for (i in seq_along(allFiles)) {
    fpkmMatrix[, i] <- allFiles[[i]]$"FPKM"
}
Sys.time() - t0 # 1 sec !

geneName <- allFiles[[1]]$gene_id
dim(fpkmMatrix)
fpkmMatrix <- fpkmMatrix[grep("ENSG", geneName, fixed = TRUE), ]
geneName <- geneName[grep("ENSG", geneName, fixed = TRUE)]
dim(fpkmMatrix)

summary(as.vector(fpkmMatrix))
keep <- rowSums(fpkmMatrix) > 0
fpkmMatrix <- fpkmMatrix[keep, ]
geneName <- geneName[keep]

t0 <- Sys.time()
corMat <- cor(fpkmMatrix)
Sys.time() - t0 # 14 sec

hg38md <- dplyr::select(hg38, Fileaccession, Experimentaccession, Biosampletermname,
                 Biosampletype, Biosamplesubcellularfractiontermname, Librarymadefrom,
                 Librarydepletedin, Derivedfrom, FiledownloadURL)

hg38md$name <- paste(hg38md$Fileaccession, hg38md$Biosampletermname)
hg38md$url <- paste0("https://www.encodeproject.org/experiments/", hg38md$Experimentaccession, "/")

unique(hg38md$Biosamplesubcellularfractiontermname)
unique(hg38md$Librarymadefrom)
unique(hg38md$Librarydepletedin)

hg38md$Librarydepletedin2 <- ifelse(is.na(hg38md$Librarydepletedin), "", paste(", depleted in", hg38md$Librarydepletedin))
hg38md$Biosamplesubcellularfractiontermname2 <- ifelse(is.na(hg38md$Biosamplesubcellularfractiontermname), "", paste(",", hg38md$Biosamplesubcellularfractiontermname))

hg38md$rnaFraction <- paste0(
    hg38md$Librarymadefrom,
    hg38md$Librarydepletedin2,
    hg38md$Biosamplesubcellularfractiontermname2
)
table(hg38md$rnaFraction)

geneName <-

metadata <- data.frame(
    name = hg38md$name,
    fileAccession = hg38md$Fileaccession,
    experimentAccession = hg38md$Experimentaccession,
    biosampleTermname = hg38md$Biosampletermname,
    biosampletype = hg38md$Biosampletype,
    rnaFraction = hg38md$rnaFraction,
    derivedFrom = hg38md$Derivedfrom,
    url = hg38md$url,
    stringsAsFactors = FALSE
)


encode_rnaseq <- list(
    "dataMatrix" = fpkmMatrix,
    "geneName" = geneName,
    "correlationMatrix" = corMat,
    "annotation" = metadata
)

lapply(encode_rnaseq, dim)
length(geneName)

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
save(encode_rnaseq, file = "heatrnaseq/data/encode_rnaseq.RData")

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

colnames(encode_rnaseq$annotation) <- c("encodeAccession", "tissue", "name", "sampleType", "rnaFraction", "url")
object.size(encode_rnaseq) # 125 Mo

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
save(encode_rnaseq, file = "heatrnaseq/data/encode_rnaseq.RData")


