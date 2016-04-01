setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
library(readr)
library(parallel)
library(GenomicRanges)

# data loading and formating
a <- Sys.time()
merged_data <- read_tsv("data/merged_wgEncodeRegTfbsClusteredV3.bed", col_names = FALSE)
Sys.time() - a # 5 seconds
annotation <- read_tsv("data/wgEncodeRegTfbsClusteredInputsV3.tab", col_names = FALSE)
Sys.time() - a # instant

colnames(merged_data) <- c("chr", "start", "end", "experiment")
head(merged_data)
colnames(annotation) <- c("peakFile", "name", "TF", "antibody", "cell_line", "treatment", "lab")
annotation <- cbind("index" = rownames(annotation), annotation)
head(annotation)

# creation of a binary data matrix
templateVector <- vector(mode = "logical", length = nrow(annotation))
names(templateVector) <- 1:nrow(annotation)

createBoolVectorForLine <- function(line) {
    templateVector[as.numeric(strsplit(merged_data[line, "experiment"], split=",")[[1]]) + 1] <- TRUE # they start at 0, hence the "+1", thanks R...
    return(templateVector)
}

a <- Sys.time()
dataMatrix <- do.call(rbind, mclapply(1:nrow(merged_data), createBoolVectorForLine, mc.cores = 32))
Sys.time() - a # 2 minutes, mostly slow rbind

#summary(as.vector(dataMatrix))
# long
#  Mode     FALSE      TRUE      NA's
# logical 503140147  12020753         0

#which(!apply(dataMatrix, 2, any))
# named integer(0)

# bed as GRanges for future overlap with user peak list
regionMetaData <- with(merged_data, GRanges(chr, IRanges(start, end)))
regionMetaData

# precomputation off correlations
a <- Sys.time()
correlationMatrix <- cor(dataMatrix)
Sys.time() - a # ridiculously long, like 5 minutes.

dim(correlationMatrix)
summary(as.vector(correlationMatrix))

# > dim(correlationMatrix)
# [1] 691 691
# > summary(as.vector(correlationMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.05517  0.04175  0.09748  0.15110  0.19740  1.00000

# plot.ly don't like duplicate names...
expNames <- table(annotation$name)
expNames[which(expNames>1)]
for (n in names(expNames[which(expNames>1)])) {
    annotation[which(annotation$name == n)[1], "name"] <- paste0(n,"_(1)")
    annotation[which(annotation$name == n)[1], "name"] <- paste0(n,"_(2)") # this is NOT a typo
    message(n)
}
expNames <- table(annotation$name)
expNames[which(expNames>1)]

encode <- list("dataMatrix" = dataMatrix,
               "regionMetaData" = regionMetaData,
               "correlationMatrix" = correlationMatrix,
               "annotation" = annotation)
encode$annotation <- encode$annotation[, 2:8]
encode$annotation$url <- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/"
names(encode$annotation) <- c("peakFile", "name", "tf", "antibody", "cellLine", "treatment", "laboratory", "url")

object.size(encode) # 2Go :-S
save(encode, file = "heatchipseq/data/encode.RData")
