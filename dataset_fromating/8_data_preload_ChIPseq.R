### ChIP seq -----------
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatchipseq/")
library(GenomicRanges)
library(magrittr)

load("data/encode.RData")
load("data/codex.RData")
load("data/codex_human_chip.RData")

names(encode)
names(codex)
names(codex_human_chip)

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
object.size(codex_human_chip)
codex_human_chip <- nullifyDataForFasterPreloading(codex_human_chip)
object.size(codex_human_chip)

save(encode, file = "data/encode_preload.RData")
save(codex, file = "data/codex_preload.RData")
save(codex_human_chip, file = "data/codex_human_chip_preload.RData")

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
load("data/blueprint_rnaseq.RData")
load("data/roadmap_rnaseq.RData")

names(encode_rnaseq)
names(bgee_human)
names(bgee_mouse)
names(encode_mouse_rnaseq)
names(blueprint_rnaseq)
names(roadmap_rnaseq)

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
object.size(blueprint_rnaseq)
blueprint_rnaseq <- nullifyDataForFasterPreloading(blueprint_rnaseq)
object.size(blueprint_rnaseq)
object.size(roadmap_rnaseq)
roadmap_rnaseq <- nullifyDataForFasterPreloading(roadmap_rnaseq)
object.size(roadmap_rnaseq)

save(encode_rnaseq, file = "data/encode_rnaseq_preload.RData")
save(bgee_human, file = "data/bgee_human_preload.RData")
save(bgee_mouse, file = "data/bgee_mouse_preload.RData")
save(encode_mouse_rnaseq, file = "data/encode_mouse_rnaseq_preload.RData")
save(blueprint_rnaseq, file = "data/blueprint_rnaseq_preload.RData")
save(roadmap_rnaseq, file = "data/roadmap_rnaseq_preload.RData")

# CAGE seq
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatcageseq/")
library(GenomicRanges)
library(magrittr)

load("data/fantom5_human_cage.RData")
load("data/fantom5_mouse_cage.RData")
names(fantom5_human_cage)
names(fantom5_mouse_cage)

nullifyDataForFasterPreloading <- function(myList) {
    myList[["dataMatrix"]] <- NULL
    myList[["regionMetaData"]] <- NULL
    myList[["correlationMatrix"]] <- NULL
    return(myList)
}

object.size(fantom5_human_cage)
fantom5_human_cage <- nullifyDataForFasterPreloading(fantom5_human_cage)
object.size(fantom5_human_cage)
object.size(fantom5_mouse_cage)
fantom5_mouse_cage <- nullifyDataForFasterPreloading(fantom5_mouse_cage)
object.size(fantom5_mouse_cage)

save(fantom5_human_cage, file = "data/fantom5_human_cage_preload.RData")
save(fantom5_mouse_cage, file = "data/fantom5_mouse_cage_preload.RData")


