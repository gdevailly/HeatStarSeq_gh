# flybase rnaseq
# file downloaded from ftp://ftp.flybase.net/releases/current/precomputed_files/genes/ 2016-03-15
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/flybase_rnaseq")
library(dplyr)
library(readr)
library(parallel)

fbdata <- read_tsv("gene_rpkm_report_fb_2016_01.tsv", skip = 5)
# end of file is messy. Seriously flybase?

names(fbdata) <- c("Release_ID", "FBgn", "GeneSymbol", "Parent_library_FBlc", "Parent_library_name", "RNASource_FBlc", "RNASource_name",
                 "RPKM_value", "Bin_value", "Unique_exon_base_count", "Total_exon_base_count", "Count_used")
fbdata <- filter(fbdata, !is.na(FBgn))
fbdata <- filter(fbdata, Count_used == "Unique")

length(unique(fbdata$RNASource_name)) # 124
table(fbdata$RNASource_name) # messy
table(fbdata$RNASource_name) %>% range # 17154 17680
length(unique(fbdata$GeneSymbol)) # 17571
length(unique(fbdata$FBgn)) # 17681, this will be fun

uniqueGenes <- unique(fbdata$FBgn)[order(unique(fbdata$FBgn))]

t0 <- Sys.time()
fbdatalist <- mclapply(unique(fbdata$RNASource_name), function(x) {
    subset(fbdata, fbdata$RNASource_name == x)
}, mc.cores = 12)
Sys.time() - t0 # 7 sec
names(fbdatalist) <- unique(fbdata$RNASource_name)

setdiff(fbdatalist$mE_mRNA_P8_CNS$FBgn, fbdatalist$mE_mRNA_P8$FBgn)

dataMatrix <- matrix(0, nrow = length(uniqueGenes), ncol = length(unique(fbdata$RNASource_name)))

t0 <- Sys.time()
for (i in 1:length(fbdatalist)) {
    exprVector <- fbdatalist[[i]]$RPKM_value
    names(exprVector) <- fbdatalist[[i]]$FBgn
    dataMatrix[, i] <- exprVector[uniqueGenes]
}
Sys.time() - t0 # 1.38 sec
dataMatrix[is.na(dataMatrix)] <- 0

corMatrix <- cor(dataMatrix)
summary(as.vector(corMatrix))
# heatmap(corMatrix)

infoOfInterest <- function(myTable) {
    return(c(
        myTable$Parent_library_name[1],
        myTable$Parent_library_FBlc[1],
        myTable$RNASource_name[1],
        myTable$RNASource_FBlc[1]
    ))
}

metadata <- do.call(rbind, lapply(fbdatalist, infoOfInterest)) %>% data.frame(stringsAsFactors = FALSE)
rownames(metadata) <- NULL
colnames(metadata) <- c("parentLibrary", "parentLibraryCode", "name", "nameCode")
metadata$url <- "ftp://ftp.flybase.net/releases/current/precomputed_files/genes/"

flybase_rnaseq <- list(
    dataMatrix = dataMatrix,
    geneName = uniqueGenes,
    correlationMatrix = corMatrix,
    annotation = metadata
)

# log transform
flybase_rnaseq$dataMatrix <- log10(flybase_rnaseq$dataMatrix + 1)
flybase_rnaseq$correlationMatrix <- cor(flybase_rnaseq$dataMatrix)

save(flybase_rnaseq, file = "../../heatrnaseq/data/flybase_rnaseq.RData")
