# 2016-05-09
# from encode DCC website:
# https://www.encodeproject.org/batch_download/type%3DExperiment%26assay_slims%3DTranscription%26assay_title%3DpolyA%2BmRNA%2BRNA-seq%26assembly%3Dhg19%26assay_title%3DRNA-seq%26assay_title%3DshRNA%2BRNA-seq%26assay_title%3DpolyA%2Bdepleted%2BRNA-seq%26files.file_type%3Dtsv
# clic on link, no wget (I think?)
# then xargs -n 1 curl -O -L < files.txt
# as ENCODE suggest
# 4496 files!? 2 files / experiment + biological repplicates...
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/EncodeRNAseq/new/")

library(readr)
library(dplyr)
library(parallel)
library(biomaRt)

metadata <- read_tsv("metadata.tsv",
                     col_types = paste0(
                         paste0(rep("?", 11), collapse = ""),
                         "c",
                         paste0(rep("?", 41 - 12), collapse = ""),
                         collapse = "")
                     )
dim(metadata) #  4496   41
length(list.files()) # 4498

colnames(metadata) <- gsub(" ", "", colnames(metadata))
unique(metadata$Fileformat) # cool
unique(metadata$Outputtype)
fmetadata <- filter(metadata, Outputtype == "gene quantifications")
unique(fmetadata$Assay) # cool
table(fmetadata$Assay) # divide dataset? probably...
# RNA-seq shRNA knockdown followed by RNA-seq
# 606                                1641
fmetadata <- filter(fmetadata, Assay == "RNA-seq")
unique(fmetadata$Biosampletype) # ok
unique(fmetadata$Biosampleorganism) # ok
table(fmetadata$Biosampletreatments) # 4 sampls treated with dimethyl sulfoxide
table(fmetadata$Biosamplesubcellularfractiontermname)
unique(fmetadata$Biosamplephase) # good
unique(fmetadata$Biosamplesynchronizationstage) # good
table(fmetadata$Librarymadefrom)
# polyadenylated mRNA                 RNA
# 106                 500
table(fmetadata$Librarydepletedin)
# rRNA rRNA, polyadenylated mRNA
# 526                        40
table(fmetadata$Project) # ok
table(fmetadata$Assembly)
# GRCh38   hg19 # sight...
# 302    304

hg19 <- filter(fmetadata, Assembly == "hg19") %>% arrange(Experimentaccession)
hg38 <- filter(fmetadata, Assembly == "GRCh38") %>% arrange(Experimentaccession)

union(hg19$Experimentaccession, hg38$Experimentaccession) %>% length # 155
intersect(hg19$Experimentaccession, hg38$Experimentaccession) %>% length # 154
setdiff(hg19$Experimentaccession, hg38$Experimentaccession) # ENCSR000AAM

filter(fmetadata, Experimentaccession == "ENCSR000AAM") %>% glimpse # why don't you like it ENCODE?

dim(fmetadata)
paste0(fmetadata$Fileaccession, ".tsv") %in% list.files() %>% which %>% length # 606 \o/

t0 <- Sys.time()
allFiles <- mclapply(hg38$Fileaccession, function(x) read_tsv(paste0(x, ".tsv")), mc.preschedule = FALSE, mc.cores = 32)
Sys.time() - t0 # 1 min
names(allFiles) <- hg38$Fileaccession

table(sapply(allFiles, nrow))
# 61471
# 302

table(sapply(allFiles, ncol))
# 15
# 302

fpkmMatrix <- matrix(0, nrow = nrow(allFiles[[1]]), ncol = length(allFiles))

t0 <- Sys.time()
for (i in seq_along(allFiles)) {
    fpkmMatrix[, i] <- allFiles[[i]]$"FPKM"
}
Sys.time() - t0 # 1 sec !

geneName <- allFiles[[1]]$gene_id
dim(fpkmMatrix)
fpkmMatrix <- fpkmMatrix[grep("ENSG", geneName, fixed = TRUE), ]
geneName <- geneName[grep("ENSG", geneName, fixed = TRUE)]
dim(fpkmMatrix)

summary(as.vector(fpkmMatrix))
keep <- rowSums(fpkmMatrix) > 0
fpkmMatrix <- fpkmMatrix[keep, ]
geneName <- geneName[keep]

t0 <- Sys.time()
corMat <- cor(fpkmMatrix)
Sys.time() - t0 # 14 sec

hg38md <- dplyr::select(hg38, Fileaccession, Experimentaccession, Biosampletermname,
                 Biosampletype, Biosamplesubcellularfractiontermname, Librarymadefrom,
                 Librarydepletedin, Derivedfrom, FiledownloadURL)

hg38md$name <- paste(hg38md$Fileaccession, hg38md$Biosampletermname)
hg38md$url <- paste0("https://www.encodeproject.org/experiments/", hg38md$Experimentaccession, "/")

unique(hg38md$Biosamplesubcellularfractiontermname)
unique(hg38md$Librarymadefrom)
unique(hg38md$Librarydepletedin)

hg38md$Librarydepletedin2 <- ifelse(is.na(hg38md$Librarydepletedin), "", paste(", depleted in", hg38md$Librarydepletedin))
hg38md$Biosamplesubcellularfractiontermname2 <- ifelse(is.na(hg38md$Biosamplesubcellularfractiontermname), "", paste(",", hg38md$Biosamplesubcellularfractiontermname))

hg38md$rnaFraction <- paste0(
    hg38md$Librarymadefrom,
    hg38md$Librarydepletedin2,
    hg38md$Biosamplesubcellularfractiontermname2
)
table(hg38md$rnaFraction)

geneName <- sub("\\.[0-9]*", "", geneName)

metadata <- data.frame(
    name = hg38md$name,
    fileAccession = hg38md$Fileaccession,
    experimentAccession = hg38md$Experimentaccession,
    biosampleTermname = hg38md$Biosampletermname,
    biosampletype = hg38md$Biosampletype,
    rnaFraction = hg38md$rnaFraction,
    derivedFrom = hg38md$Derivedfrom,
    url = hg38md$url,
    stringsAsFactors = FALSE
)


encode_rnaseq <- list(
    "dataMatrix" = fpkmMatrix,
    "geneName" = geneName,
    "correlationMatrix" = corMat,
    "annotation" = metadata
)

lapply(encode_rnaseq, dim)
length(geneName)

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
save(encode_rnaseq, file = "heatrnaseq/data/encode_rnaseq.RData")
