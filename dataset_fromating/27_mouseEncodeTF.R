# file download ------------
# from https://www.encodeproject.org/search/?type=Experiment&assay_title=ChIP-seq&replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&target.investigated_as=transcription+factor&target.investigated_as=chromatin+remodeller&files.file_type=bed+narrowPeak
# 2017-06-07
# 179 files
# filters: ChIP-seq, Mus musculus, transcription factor & chromatin remodeller, bed narrowPeak

library(dplyr)
library(readr)
library(purrr)
library(future); plan(multiprocess(workers = 16))

setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/mouseEncodeTF")

metadata <- read_tsv("metadata.tsv", progress = FALSE)
glimpse(metadata)

table(metadata[, "Output type"] %>% unlist)

unique(metadata[, "File format"] %>% unlist) # ok
unique(metadata[, "Output type"] %>% unlist) # not ok
unique(metadata[, "Assay"] %>% unlist) # ok
unique(metadata[, "Biosample term name"] %>% unlist) # to keep
unique(metadata[, "Biosample type"] %>% unlist) #
unique(metadata[, "Biosample organism"] %>% unlist) # ok
unique(metadata[, "Biosample treatments"] %>% unlist) # erf
unique(metadata[, "Biosample subcellular fraction term name"] %>% unlist) # ok
unique(metadata[, "Biosample phase"] %>% unlist) # ok
unique(metadata[, "Biosample synchronization stage"] %>% unlist) # ok
unique(metadata[, "Library made from"] %>% unlist) # ok
unique(metadata[, "Antibody accession"] %>% unlist) #
unique(metadata[, "Library depleted in"] %>% unlist) #
table(metadata[, "Assembly"] %>% unlist) # f*
unique(metadata[, "Technical replicate"] %>% unlist) # ok
unique(metadata[, "Biological replicate(s)"] %>% unlist) #

smallMD <- metadata[, c(
    "File accession",
    "Output type",
    "Experiment accession",
    "Biosample term name",
    "Biosample type",
    "Biosample sex",
    "Biosample treatments",
    "Experiment target",
    "Antibody accession",
    "File download URL",
    "Biological replicate(s)",
    "Assembly"
)]
colnames(smallMD) <- gsub(" ", "_", colnames(smallMD))
colnames(smallMD)[11] <- "Biological_replicate"
table(smallMD$Assembly, smallMD$Output_type)
smallMD <- filter(smallMD, Assembly == "mm10", Output_type == "optimal idr thresholded peaks") %>%
    select(
        File_accession,
        Experiment_accession,
        Biosample_term_name,
        Biosample_type,
        Biosample_sex,
        Biosample_treatments,
        Experiment_target,
        Antibody_accession,
        File_download_URL,
        Biological_replicate
) %>% mutate(
    Experiment_target = gsub("-mouse$", "", Experiment_target),
    name = paste0(
        Experiment_target, "_",
        Biosample_term_name, "_",
        ifelse(is.na(Biosample_treatments), "", paste0(Biosample_treatments, "_")),
        File_accession
    )
)

nrow(smallMD) # 156
length(unique(smallMD$name)) # 156 \o/
smallMD$name[duplicated(smallMD$name)]


myPeakTmp <- read_tsv("ENCFF128TGM.bed.gz", progress = FALSE, col_names = FALSE, col_types = "ciicicdddi")

t0 <- Sys.time()
myPeaks <- map_df(seq_len(nrow(smallMD)), function(i) {
    return(
        read_tsv(
            paste0(smallMD[i, ]$File_accession, ".bed.gz"),
            col_names = FALSE,
            progress = FALSE,
            col_types = "ciicicdddi"
        ) %>%
            select(X1, X2, X3) %>%
            rename(chr = X1, start = X2, end = X3) %>%
            mutate(name = smallMD[i, ]$name)
    )
})
Sys.time() - t0

nrow(myPeaks) # 3800966
length(unique(myPeaks$name)) # 156

write.table(myPeaks, file = "mouse_encode_tf_all_peaks.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# bedtools merge ---------------------
# module load apps/gcc/BEDTools/2.25.0
# sortBed -i mouse_encode_tf_all_peaks.bed > mouse_encode_tf_all_peaks.sorted.bed
# bedtools merge -c 4 -o distinct -i mouse_encode_tf_all_peaks.sorted.bed > mouse_encode_tf_all_peaks.merged.bed
#

myMergedPeaks <- read_tsv(
    "mouse_encode_tf_all_peaks.merged.bed",
    col_names = FALSE,
    progress = FALSE,
    col_types = "ciic"
)
colnames(myMergedPeaks) <- c("chr", "start", "end", "experiments")

table(myMergedPeaks$chr)
# removal on sexual chormosomes and unconventional ones
myMergedPeaks <- dplyr::filter(
    myMergedPeaks,
    chr %in% paste0("chr", 1:19)
)
table(myMergedPeaks$chr)

summary(myMergedPeaks$end - myMergedPeaks$start) # seems ok

template <- matrix(FALSE, nrow = 1, ncol = nrow(smallMD), dimnames = list(NULL, smallMD$name)) %>%
    as_data_frame

t0 <- Sys.time() # 8 minutes
dataMatrix2 <- map_df(
    seq_len(nrow(myMergedPeaks)),
    function(i) {
        template[, strsplit(myMergedPeaks[i, ]$experiments, split = ",")[[1]]] <- TRUE
        return(template)
    }
) %>% as.matrix
Sys.time() - t0

colnames(dataMatrix2) <- NULL
length(which(rowSums(dataMatrix2) == 0))

t0 <- Sys.time() # 30 sec
correlationMatrix <- cor(dataMatrix2)
Sys.time() - t0

library(GenomicRanges)
regionMetaData <- with(myMergedPeaks, GRanges(chr, IRanges(start, end)))

annotation <- smallMD
annotation$url <- paste0("https://www.encodeproject.org/files/", annotation$File_accession)
annotation <- dplyr::select(
    annotation,
    name, Experiment_target, Biosample_term_name, Biosample_treatments,
    File_accession, Experiment_accession, Biosample_type, Biosample_sex,
    Antibody_accession, url
) %>%
    dplyr::rename(
        target = Experiment_target, cell_type = Biosample_term_name,
        treatment = Biosample_treatments
    )

encode_mouse <- list(
    dataMatrix = dataMatrix2,
    regionMetaData = regionMetaData,
    correlationMatrix = correlationMatrix,
    annotation = data.frame(annotation, stringsAsFactors = FALSE)
)

map(encode_mouse, dim)
save(encode_mouse, file = "../../heatchipseq/data/encode_mouse.RData")

nullifyDataForFasterPreloading <- function(myList) {
    myList[["dataMatrix"]] <- NULL
    myList[["regionMetaData"]] <- NULL
    myList[["correlationMatrix"]] <- NULL
    return(myList)
}
object.size(encode_mouse)
encode_mouse <- nullifyDataForFasterPreloading(encode_mouse)
object.size(encode_mouse)
save(encode_mouse, file = "../../heatchipseq/data/encode_mouse_preload.RData")


