# GTEX data, 2016-04-14
# wget http://gtexportal.org/static/datasets/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz
# wget http://gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/gtex")
library(dplyr)
library(readr)
library(cba)

t0 <- Sys.time()
gtexdata <- read_tsv("GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", skip = 2)
Sys.time() - t0 # 3 min

metadata <- read_tsv("GTEx_Data_V6_Annotations_SampleAttributesDS.txt")

gtexNumData <- as.matrix(gtexdata[, c(-1,-2)])

dim(gtexdata)
dim(gtexNumData)

keep <- which(rowSums(gtexNumData) > 0)
nrow(gtexNumData)
length(keep)

geneNames <- sub("\\.[0-9]*", "", gtexdata$Name)
table(nchar(geneNames))

geneNames <- geneNames[keep]
gtexNumData <- gtexNumData[keep, ]

experiences <- colnames(gtexdata)[c(-1,-2)]

length(experiences)
intersect(experiences, metadata$SAMPID) %>% length # \o/

dim(metadata)
metadata <- filter(metadata, SAMPID %in% experiences)
dim(metadata)
metadata <- arrange(metadata, SAMPID)
dim(gtexNumData)
gtexNumData <- gtexNumData[, metadata$SAMPID]
dim(gtexNumData)

t0 <- Sys.time()
corMatrix <- cor(gtexNumData)
Sys.time() - t0 # 2h

# t0 <- Sys.time()
# d <- dist(corMatrix) # to long, did not run
# Sys.time() - t0
# myClust <- hclust(d)
# Sys.time() - t0
# co <- order.optimal(d, myClust$merge)
# Sys.time() - t0
# myClust$merge <- co$merge
# myClust$order <- co$order
# dendro <- as.dendrogram(myClust)
# Sys.time() - t0

metadata <- select(metadata, SAMPID, SMTS, SMTSD, SMNABTCHT)
metadata$name <- paste(metadata$SMTSD, gsub("GTEX-", "", metadata$SAMPID))
metadata$url <- paste0("http://www.ncbi.nlm.nih.gov/sra/?term=", metadata$SAMPID)

gtex_large <- list(
    dataMatrix = gtexNumData,
    geneName = geneNames,
    correlationMatrix = corMatrix,
    annotation = as.data.frame(metadata)
)
gtex_large$annotation$name <- as.character(gtex_large$annotation$name)
gtex_large$annotation$SAMPID <- as.character(gtex_large$annotation$SAMPID)
gtex_large$annotation$SMTS <- as.character(gtex_large$annotation$SMTS)
gtex_large$annotation$SMTSD <- as.character(gtex_large$annotation$SMTSD)
gtex_large$annotation$SMNABTCHT <- as.character(gtex_large$annotation$SMNABTCHT)
gtex_large$annotation$url <- as.character(gtex_large$annotation$url)

save(gtex_large, file = "gtex_large.RData")





###
t0 <- Sys.time()
a <- dist(gtex_large$correlationMatrix[, 1:10])
Sys.time() - t0 # 5 s

t0 <- Sys.time()
a <- dist(gtex_large$correlationMatrix[, 1:100])
Sys.time() - t0 # 40 s

t0 <- Sys.time()
a <- dist(gtex_large$correlationMatrix[, 1:1000])
Sys.time() - t0 # 30 minutes...

t0 <- Sys.time()
a <- as.dist(1 - gtex_large$correlationMatrix)
Sys.time() - t0 # 3 sec

t0 <- Sys.time()
myClust <- hclust(a)
Sys.time() - t0 # 5sec
co <- order.optimal(a, myClust$merge)
Sys.time() - t0 # 5 min
myClust$merge <- co$merge
myClust$order <- co$order
dendro <- as.dendrogram(myClust)
Sys.time() - t0 # fast



