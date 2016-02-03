# NOISeq TMM
# tmm() doesn't change correlation value
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/")
library(shiny)
library(GenomicRanges)
library(magrittr)
library(readr)
library(NOISeq)

###
load("heatrnaseq/data/encode_rnaseq.RData")
dataset <- encode_rnaseq
myFile <- read_tsv("Blood_microaaray_1.txt")
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
# normData <- tmm(cbind(dataset$dataMatrix, userExpressionFile))
normData <- uqua(cbind(dataset$dataMatrix, userExpressionFile))
Sys.time() - t0 # 8 sec for tmm, 3 sec for uqua

cornormData <- cor(normData)
normUserCorrelations <- cornormData[1:(nrow(cornormData) - 1), ncol(cornormData)]
cornormData <- cornormData[-nrow(cornormData), -nrow(cornormData)]

summary(as.vector(dataset$correlationMatrix))
summary(as.vector(cornormData))
summary(userCorrelations)
summary(normUserCorrelations)

# tmm, encode mouse rnaseq, mk3
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

# uqua, encode mouse rnaseq, mk3
# > summary(as.vector(dataset$correlationMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.002985 0.230300 0.442600 0.476000 0.737200 1.000000
# > summary(as.vector(cornormData))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# -0.0008  0.2519  0.5353  0.5219  0.7823  1.0000    1140
# > summary(userCorrelations)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.001383 0.073440 0.125500 0.114200 0.152800 0.333200
# > summary(normUserCorrelations)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's
# 0.000322 0.184400 0.291300 0.273600 0.356600 0.765500        3

# \o/

which(normUserCorrelations == max(normUserCorrelations, na.rm = TRUE)) # 70

# uqua, gbee human, blood microarray
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

# uqua, encode human, blood microarray
# > summary(as.vector(dataset$correlationMatrix))
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0001169  0.2031000  0.7595000  0.5957000  0.9561000  1.0000000
# > summary(as.vector(cornormData))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# -0.0041  0.2893  0.8703  0.6549  0.9769  1.0000    1908
# > summary(userCorrelations)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0014490 -0.0005227  0.0029590  0.0581200  0.0338300  0.4004000
# > summary(normUserCorrelations)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's
# -0.009219 -0.008462 -0.006169  0.057810  0.024150  0.473500         3

# uqua, bgee mouse, mk3
# > summary(as.vector(dataset$correlationMatrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02765 0.18880 0.32290 0.39750 0.56890 1.00000
# > summary(as.vector(cornormData))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02765 0.18880 0.32290 0.39750 0.56890 1.00000
# > summary(userCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01796 0.09100 0.15920 0.16220 0.23670 0.37030
# > summary(normUserCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01796 0.09100 0.15920 0.16220 0.23670 0.37030

# why bugs ??? -> due to 3rd quantile being 0, messing all values. UQUA do not change correlation values.

myData <- list(
    noNorm = cbind(dataset$dataMatrix, userExpressionFile),
    noK = uqua(cbind(dataset$dataMatrix, userExpressionFile)),
    kequal1 = uqua(cbind(dataset$dataMatrix, userExpressionFile), k = 1)
)
myCors <- lapply(myData, cor)
lapply(myCors, function(x) summary(as.vector(x)))
lapply(myData, function(x) x[1:10, 1:4])

head(myData$noNorm)
head(myData$noK)
dim(myData$noNorm) # [1] 57820   321
summary(myData$noK[, 320])
summary(myData$noK[, 319])
summary(myData$noNorm[, 320])
summary(myData$noNorm[, 319])

length(which(rowSums(dataset$dataMatrix) == 0))

min(dataset$dataMatrix[dataset$dataMatrix != 0])
dataset$dataMatrix[dataset$dataMatrix == 0] <- 0.01
