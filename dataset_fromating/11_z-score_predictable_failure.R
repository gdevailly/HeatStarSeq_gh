# z scores doesn't modify pearson correlation value.
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/")
library(shiny)
library(GenomicRanges)
library(magrittr)
library(readr)
library(preprocessCore)

###
load("heatrnaseq/data/encode_mouse_rnaseq.RData")
dataset <- encode_mouse_rnaseq
myFile <- read_tsv("mK3_1.txt")
###

userExpressionFile_temp <- myFile
colnames(userExpressionFile_temp) <- c("ensembl_id","exp_value")
userExpressionFile_temp_v <- userExpressionFile_temp$exp_value
names(userExpressionFile_temp_v) <- userExpressionFile_temp$ensembl_id

userExpressionFile <- userExpressionFile_temp_v[dataset$geneName]
userExpressionFile[is.na(userExpressionFile)] <- 0
userCorrelations <- cor(userExpressionFile, dataset$dataMatrix) %>% as.vector

dim(dataset$dataMatrix)
t0 <- Sys.time()
normData <- scale(dataset$dataMatrix)
Sys.time() - t0 # 2 sec

cornormData <- cor(normData)
normUserData <- scale(userExpressionFile)
normUserCorrelations <- cor(normUserData, normData) %>% as.vector

summary(as.vector(dataset$correlationMatrix))
summary(as.vector(cornormData))
summary(userCorrelations)
summary(normUserCorrelations)

# bgee human + Blood
# > summary(as.vector(dataset$correlationMatrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1339  0.3721  0.5689  0.5646  0.7538  1.0000
# > summary(as.vector(cornormData))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1339  0.3721  0.5689  0.5646  0.7538  1.0000
# > summary(userCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1730  0.3287  0.3921  0.3631  0.4117  0.4891
# > summary(normUserCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1730  0.3287  0.3921  0.3631  0.4117  0.4891

# encode mouse + mk3
# > summary(as.vector(dataset$correlationMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.002985 0.230300 0.442600 0.476000 0.737200 1.000000
# > summary(as.vector(cornormData))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.002985 0.230300 0.442600 0.476000 0.737200 1.000000
# > summary(userCorrelations)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.001383 0.073440 0.125500 0.114200 0.152800 0.333200
# > summary(normUserCorrelations)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.001383 0.073440 0.125500 0.114200 0.152800 0.333200

