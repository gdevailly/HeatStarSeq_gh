setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/")
library(magrittr)
library(readr)
library(preprocessCore)

###
load("heatrnaseq/data/bgee_mouse.RData")
dataset <- bgee_mouse
myFile <- read_tsv("mK3_1.txt")
###

userExpressionFile_temp <- myFile
colnames(userExpressionFile_temp) <- c("ensembl_id","exp_value")
userExpressionFile_temp_v <- userExpressionFile_temp$exp_value
names(userExpressionFile_temp_v) <- userExpressionFile_temp$ensembl_id

userExpressionFile <- userExpressionFile_temp_v[dataset$geneName]
userExpressionFile[is.na(userExpressionFile)] <- 0
userCorrelations <- cor(userExpressionFile, dataset$dataMatrix) %>% as.vector

t0 <- Sys.time()
normData <- normalize.quantiles(dataset$dataMatrix)
Sys.time() - t0 # 3 sec

cornormData <- cor(normData)

normUserData <- normalize.quantiles.use.target(
    as.matrix(userExpressionFile),
    normalize.quantiles.determine.target(dataset$dataMatrix)
) %>% as.vector

normUserCorrelations <- cor(normUserData, normData) %>% as.vector

summary(as.vector(dataset$correlationMatrix))
summary(as.vector(cornormData))
summary(userCorrelations)
summary(normUserCorrelations)

# encode mouse + mk3
# > summary(as.vector(dataset$correlationMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.002884 0.230000 0.442300 0.475600 0.736600 1.000000
# > summary(as.vector(cornormData))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.005868 0.115200 0.221100 0.329700 0.498600 1.000000
# > summary(userCorrelations)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.001331 0.073290 0.125300 0.114000 0.152600 0.333200
# > summary(normUserCorrelations)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.008046 0.054060 0.068970 0.088000 0.103300 0.438500
dataset$annotation[which(userCorrelations == max(userCorrelations)), ]
dataset$annotation[which(normUserCorrelations == max(normUserCorrelations)), ]

# Bgee mouse + mK3
# > summary(as.vector(dataset$correlationMatrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.02753 0.18860 0.32270 0.39740 0.56900 1.00000
# > summary(as.vector(cornormData))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01005 0.21920 0.35810 0.41710 0.58850 1.00000
# > summary(userCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01782 0.09094 0.15920 0.16210 0.23660 0.37030
# > summary(normUserCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01616 0.17990 0.29010 0.29820 0.42720 0.59460
dataset$annotation[which(userCorrelations == max(userCorrelations)), ]
dataset$annotation[which(normUserCorrelations == max(normUserCorrelations)), ]


# encode human + blood
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0001849  0.2029000  0.7595000  0.5957000  0.9561000  1.0000000
# > summary(as.vector(cornormData))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000047 0.0161500 0.6738000 0.5547000 0.9633000 1.0000000
# > summary(userCorrelations)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0018070 -0.0009043  0.0023980  0.0573000  0.0328300  0.3997000
# > summary(normUserCorrelations)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -4.936e-05  1.471e-05  1.362e-04  6.227e-03  6.254e-04  1.470e-01


# bgee human + blood
# summary(as.vector(dataset$correlationMatrix))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1335  0.3718  0.5685  0.5642  0.7537  1.0000
# > summary(as.vector(cornormData))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.08426 0.33120 0.58540 0.54340 0.75630 1.00000
# > summary(userCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1726  0.3290  0.3914  0.3620  0.4116  0.4889
# > summary(normUserCorrelations)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04943 0.36370 0.46670 0.41360 0.52950 0.64040

    








