# 2016/02/08
# lftp ftp.ebi.ac.uk
# cd pub/databases/blueprint/data/homo_sapiens/
# mget -d GRCh38/*/*/*/RNA-Seq/*/*.gene_quantification.*.results
# mget -d GRCh37/*/*/*/RNA-Seq/*/*.gene_quantification.*.gff

setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/blueprint")
library(readr)
library(magrittr)

metadata <- data.frame(
    genome = "dummystring",
    level1 = "dummystring",
    level2 = "dummystring",
    level3 = "dummystring",
    level4 = "dummystring",
    fileName = "dummystring",
    stringsAsFactors = FALSE
)

hg38 <- "GRCh38"
t0<- Sys.time()
for (level1 in list.files(hg38)) {
    for(level2 in list.files(paste(hg38, level1, sep = "/"))) {
        for(level3 in list.files(paste(hg38, level1, level2, sep = "/"))) {
            for(level4 in list.files(paste(hg38, level1, level2, level3, "RNA-Seq",  sep = "/"))) {
                for(myFile in list.files(paste(hg38, level1, level2, level3, "RNA-Seq", level4, sep = "/"))) {
                    metadata <- rbind(
                        metadata,
                        c(hg38,level1, level2, level3, level4, myFile)
                    )
                }
            }
        }
    }
}
Sys.time() - t0 # <1s
metadata <- metadata[-1, ]

dim(metadata) # [1] 163   6

with(metadata, paste(genome, level1, level2, level3, level4)) %>% subset(., table(.) == 2)
subset(metadata, metadata$level2 == "N00031407013121") # should be allright

metadata$name <- sapply(metadata$fileName, function(x) strsplit(x, split = ".", fixed = TRUE)[[1]][1])
length(unique(metadata$name)) # cool

t0 <- Sys.time()
myListOfFile <- lapply(1:nrow(metadata), function(x) read_tsv(
    paste(metadata$genome[x], metadata$level1[x], metadata$level2[x], metadata$level3[x], "RNA-Seq", metadata$level4[x], metadata$fileName[x], sep = "/")
))
Sys.time() - t0 # 50s

names(myListOfFile) <- metadata$name
sapply(myListOfFile, dim) %>% table # \o/
# 15 60483
# 163   163

head(myListOfFile[[1]])
head(myListOfFile[[75]])

t0 <- Sys.time()
dataMatrix <- matrix(0, nrow = nrow(myListOfFile[[1]]), ncol = nrow(metadata))
for(i in 1:nrow(metadata)) {
    dataMatrix[, i] <- myListOfFile[[i]]$FPKM
}
Sys.time() - t0 # 0.787 sec

geneName <- sub("\\.[0-9]*", "", myListOfFile[[1]]$gene_id)
nonZero <- which(rowSums(dataMatrix) != 0) # 53667
dataMatrix <- dataMatrix[nonZero, ]
geneName <- geneName[nonZero]
corMatrix <- cor(dataMatrix)

# alt_corMatrix1 <- cor(dataMatrix, method = "kendall") # far to long, didn't complete
alt_corMatrix2 <- cor(dataMatrix, method = "spearman")

summary(as.vector(corMatrix))
summary(as.vector(alt_corMatrix2))
# > summary(as.vector(corMatrix))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0001763 0.0099270 0.1300000 0.3600000 0.7561000 1.0000000
# > summary(as.vector(alt_corMatrix2))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4825  0.6918  0.7404  0.7465  0.8125  1.0000 # interesting...
#
# pdf("testCor.pdf")
# heatmap(corMatrix)
# heatmap(alt_corMatrix2)
# dev.off()
# system("firefox testCor.pdf &")

unique(metadata$level1)
unique(metadata$level2)
unique(metadata$level3)
unique(metadata$level4)

metadata$url <- with(metadata, paste0(
    "ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/",
    level1, "/",
    level2, "/",
    level3,
    "/RNA-Seq/",
    level4, "/",
    fileName
))

lessMetadata <- metadata[, c("name", "level1", "level3", "url")]

blueprint_rnaseq <- list(
    "dataMatrix" = dataMatrix,
    "geneName" = geneName,
    "correlationMatrix" = corMatrix,
    "annotation" = lessMetadata
)

blueprint_rnaseq$annotation$name <- paste(blueprint_rnaseq$annotation$name, blueprint_rnaseq$annotation$level3)
colnames(blueprint_rnaseq$annotation) <- c("name", "sampleSource", "cellType", "url")
# log transform
blueprint_rnaseq$dataMatrix <- log10(blueprint_rnaseq$dataMatrix + 1)
blueprint_rnaseq$correlationMatrix <- cor(blueprint_rnaseq$dataMatrix)

save(blueprint_rnaseq, file = "../../heatrnaseq/data/blueprint_rnaseq.RData")





