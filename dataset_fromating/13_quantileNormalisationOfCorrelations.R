setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/")
library(shiny)
library(magrittr)
library(readr)
library(NOISeq)
library(preprocessCore)
load("heatrnaseq/data/encode_rnaseq.RData")
dataset <- encode_rnaseq
myFile <- read_tsv("Blood_microaaray_1.txt")


userExpressionFile_temp <- myFile
colnames(userExpressionFile_temp) <- c("ensembl_id","exp_value")
userExpressionFile_temp_v <- userExpressionFile_temp$exp_value
names(userExpressionFile_temp_v) <- userExpressionFile_temp$ensembl_id

userExpressionFile <- userExpressionFile_temp_v[dataset$geneName]
userExpressionFile[is.na(userExpressionFile)] <- 0
userCorrelations <- cor(userExpressionFile, dataset$dataMatrix) %>% as.vector

t0 <- Sys.time()
normalize.quantiles.determine.target(dataset$correlationMatrix)
Sys.time() - t0

t0 <- Sys.time()
normUserCorrelations <- normalize.quantiles.use.target(
    as.matrix(c(userCorrelations, 1)),
    normalize.quantiles.determine.target(dataset$correlationMatrix)
) %>% as.vector
Sys.time() - t0

summary(userCorrelations)
summary(normUserCorrelations)


normUserCorrelations2 <- normalize.quantiles.use.target(
    as.matrix(userCorrelations),
    normalize.quantiles.determine.target(dataset$correlationMatrix)
) %>% as.vector

pdf(file = "Rplots/quantNormOfCorValue_encode_rnaseq_blood.pdf")
plot(c(userCorrelations, 1), normUserCorrelations, pch = 19, xlim = c(-0.1, 1), ylim = c(-0.1, 1),
     xlab = "Raw correlation", ylab = "Adjusted correlation", main = "bgee Mouse vs mK3")
points(c(userCorrelations, 1), c(userCorrelations*0.95/max(userCorrelations), 1), pch = 18, col = "blue")
abline(a = 0, b = 1)
legend("bottomright", legend = c("Quant. norm. of cor. values", "Linear adj. 0.95"), pch = c(19, 18), col = c("black", "blue"))
dev.off()

system("firefox Rplots/quantNormOfCorValue_encode_rnaseq_blood.pdf &")
