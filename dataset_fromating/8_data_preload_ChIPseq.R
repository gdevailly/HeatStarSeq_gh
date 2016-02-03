### ChIP seq -----------
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatchipseq/")
library(GenomicRanges)
library(magrittr)

load("data/encode.RData")
load("data/codex_currated.RData")

names(encode)
names(codex)

nullifyDataForFasterPreloading <- function(myList) {
    myList[["dataMatrix"]] <- NULL
    myList[["regionMetaData"]] <- NULL
    myList[["correlationMatrix"]] <- NULL
    return(myList)
}

object.size(encode)
encode <- nullifyDataForFasterPreloading(encode)
object.size(encode)
object.size(codex)
codex <- nullifyDataForFasterPreloading(codex)
object.size(codex)

save(encode, file = "data/encode_preload.RData")
save(codex, file = "data/codex_currated_preload.RData")

load("data/encode.RData")
object.size(encode)
load("data/encode_preload.RData")
object.size(encode)

### RNA seq -----------
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatrnaseq/")
library(GenomicRanges)
library(magrittr)
load("data/encode_rnaseq.RData")
load("data/bgee_human.RData")
load("data/bgee_mouse.RData")
load("data/encode_mouse_rnaseq.RData")

names(encode_rnaseq)
names(bgee_human)
names(bgee_mouse)
names(encode_mouse_rnaseq)

nullifyDataForFasterPreloading <- function(myList) {
    myList[["dataMatrix"]] <- NULL
    myList[["geneName"]] <- NULL
    myList[["correlationMatrix"]] <- NULL
    return(myList)
}

object.size(encode_rnaseq)
encode_rnaseq <- nullifyDataForFasterPreloading(encode_rnaseq)
object.size(encode_rnaseq)
object.size(bgee_human)
bgee_human <- nullifyDataForFasterPreloading(bgee_human)
object.size(bgee_human)
object.size(bgee_mouse)
bgee_mouse <- nullifyDataForFasterPreloading(bgee_mouse)
object.size(bgee_mouse)
object.size(encode_mouse_rnaseq)
encode_mouse_rnaseq <- nullifyDataForFasterPreloading(encode_mouse_rnaseq)
object.size(encode_mouse_rnaseq)

save(encode_rnaseq, file = "data/encode_rnaseq_preload.RData")
save(bgee_human, file = "data/bgee_human_preload.RData")
save(bgee_mouse, file = "data/bgee_mouse_preload.RData")
save(encode_mouse_rnaseq, file = "data/encode_mouse_rnaseq_preload.RData")

