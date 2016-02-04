library(magrittr)
library(readr)
library(parallel)
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
# extract from CODEX 2015-09-04
codexMD <- read_tsv("data/CODEX_mm10/CODEX_mm10.txt")

# download file
codexMD$getURL <- paste0("wget ", codexMD$peakURL) 
# write.table(codexMD$getURL, file = "data/CODEX_mm10/peakFiles/getPeakFile.sh", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "/t")
codexMD$getURL <- NULL

# format file for bedmerge
codexMD$filename <- strsplit(codexMD$peakURL, split = "/") %>% lapply(function(x) tail(x, 1)) %>% unlist

# file with parenthesis and url ...
filesList <- list.files("data/CODEX_mm10/peakFiles/")
summary(filesList %in% codexMD$filename)
filesList[!(filesList %in% codexMD$filename)]
subset(codexMD, codexMD$Sample == "GSM1350476")
filesList[(filesList == "GSM1350476_RORgamma%28t%29.bed")]
codexMD[631, "filename"] <- "GSM1350476_RORgamma(t).bed"
filesList[!(filesList %in% codexMD$filename)]
codexMD$filename[!(codexMD$filename %in% filesList)]
codexMD$name <- sub(".bed", "", codexMD$filename, fixed = TRUE) 
#
a <- Sys.time()
peakList <- lapply(codexMD$filename, function(x) read_tsv(paste0("data/CODEX_mm10/peakFiles/", x), col_names = FALSE))
Sys.time() - a
names(peakList) <- codexMD$name

peakList <- lapply(peakList, function(x) x[,1:3])
peakList <- lapply(1:length(peakList), function(x) {
                        peakList[[x]]$exp <- names(peakList)[x]
                        return(peakList[[x]])        
                    })

a <- Sys.time()
bigBed <- do.call(rbind, peakList)
Sys.time() - a # < 2 min

# > dim(bigBed)
# [1] 9507832       4

# write.table(bigBed, file = "data/CODEX_mm10/allCodexMM10.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
# then bedtools :
# sortBed -i allCodexMM10.bed > sorted_allCodexMM10.bed
# bedtools merge -c 4 -o distinct -i sorted_allCodexMM10.bed > merged_allCodexMM10.bed

#########
# new R session
#########
library(magrittr)
library(readr)
library(parallel)
library(GenomicRanges)

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
# extract from CODEX 2015-09-04
codexMerged <- read_tsv("data/CODEX_mm10/merged_allCodexMM10.bed", col_names = FALSE)
dim(codexMerged)
# [1] 697160      4
codexMD <- read_tsv("data/CODEX_mm10/CODEX_mm10.txt")
codexMD$filename <- strsplit(codexMD$peakURL, split = "/") %>% lapply(function(x) tail(x, 1)) %>% unlist
codexMD[631, "filename"] <- "GSM1350476_RORgamma(t).bed"
codexMD$name <- sub(".bed", "", codexMD$filename, fixed = TRUE) 
filesList <- list.files("data/CODEX_mm10/peakFiles/")
filesList[!(filesList %in% codexMD$filename)]
codexMD$filename[!(codexMD$filename %in% filesList)]

# creation of a binary data matrix
templateVector <- vector(mode = "logical", length = nrow(codexMD))
names(templateVector) <- codexMD$name

createBoolVectorForLine <- function(line) {
    templateVector[strsplit(codexMerged[line, "X4"], split=",")[[1]]] <- TRUE 
    return(templateVector)
}

a <- Sys.time()
dataMatrix <- do.call(rbind, mclapply(1:nrow(codexMerged), createBoolVectorForLine, mc.cores = 32))
Sys.time() - a # 3 minutes, mostly slow rbind

dim(dataMatrix)
#[1] 697160    651
summary(as.vector(dataMatrix))
# Mode     FALSE      TRUE      NA's
#   logical 444775544   9075616         0
which(!apply(dataMatrix, 2, any))
# named integer(0)

# bed as GRanges for future overlap with user peak list
regionMetaData <-  GRanges(codexMerged[,1], IRanges(codexMerged[,2], codexMerged[,3]))
regionMetaData

# precomputation off correlations
a <- Sys.time()
correlationMatrix <- cor(dataMatrix)
Sys.time() - a # ridiculously long, like 4 minutes.

dim(correlationMatrix)
# [1] 651 651
summary(as.vector(correlationMatrix))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.05623  0.01721  0.05037  0.09696  0.12780  1.00000

codex <- list("dataMatrix" = dataMatrix,
               "regionMetaData" = regionMetaData,
               "correlationMatrix" = correlationMatrix,
               "annotation" = codexMD[,c("Cell type", "Cell subtype", "TF", "name")])
object.size(codex) # 1.8Go :-S
save(codex, file = "heatmap_shinyApp/data/codex.RData")
