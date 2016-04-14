
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/gtex")
library(dplyr)
library(readr)
library(cba)

t0 <- Sys.time()
gtexdata <- read_tsv("GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", skip = 2)
Sys.time() - t0 # 3 min

metadata <- read_tsv("GTEx_Data_V6_Annotations_SampleAttributesDS.txt")

gtexNumData <- as.matrix(gtexdata[, c(-1,-2)])

length(table(metadata$SMTSD))
table(metadata$SAMPID)

experiences <- colnames(gtexdata)[c(-1,-2)]
dim(metadata)
metadata <- filter(metadata, SAMPID %in% experiences)
dim(metadata)
metadata <- arrange(metadata, SAMPID)
dim(gtexNumData)
gtexNumData <- gtexNumData[, metadata$SAMPID]
dim(gtexNumData)


t0 <- Sys.time()
lmat <- lapply(unique(metadata$SMTSD), function(x) {
    gtexNumData[, metadata$SAMPID[metadata$SMTSD == x]]
})
Sys.time() - t0

t0 <- Sys.time()
lmat_median <- lapply(lmat, function(x) apply(x, 1, median))
Sys.time() - t0 # 5 mins

lapply(lmat_median, length)
gtex_small_data <- do.call(cbind, lmat_median)
colnames(gtex_small_data) <- unique(metadata$SMTSD)

keep <- which(rowSums(gtex_small_data) > 0)
nrow(gtex_small_data)
length(keep)

geneNames <- sub("\\.[0-9]*", "", gtexdata$Name)

geneNames <- geneNames[keep]
gtex_small_data <- gtex_small_data[keep, ]

myCor <- cor(gtex_small_data)

myMetadata <- data.frame(
    name = colnames(gtex_small_data),
    url = "http://gtexportal.org/home/"
)

gtex_small <- list(
    dataMatrix = gtex_small_data,
    geneName = geneNames,
    correlationMatrix = myCor,
    annotation = as.data.frame(myMetadata)
)

save(gtex_small, file = "../../heatrnaseq/data/gtex_small.RData")

