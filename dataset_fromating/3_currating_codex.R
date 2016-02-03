#########
# manual curation of CODEX mouse data...
# So much FUN!!!
#########
library(magrittr)
library(parallel)
library(GenomicRanges)
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap")
load("heatmap_shinyApp/data/codex.RData")

ct <- factor(codex$annotation[,"Cell type"])
lct <- levels(ct)
table(ct)
lct[lct == "embryonic stem cell"] <- "Embryonic Stem Cell"
lct[lct == "Embryonic stem cell"] <- "Embryonic Stem Cell"
lct[lct == "Embryonic Stem cell"] <- "Embryonic Stem Cell"
lct[lct == "embryonic stem cells"] <- "Embryonic Stem Cell"
lct[lct == "Embryonic stem cells"] <- "Embryonic Stem Cell"
lct[lct == "Embryonic stem cells"] <- "Embryonic Stem Cell"
lct[lct == "Embryonic Stem cells"] <- "Embryonic Stem Cell"
lct[lct == "Embryonic Stem Cells"] <- "Embryonic Stem Cell"
lct[lct == "Embyonic stem cell"] <- "Embryonic Stem Cell"
lct[lct == "ES cells"] <- "Embryonic Stem Cell"

lct[lct == "Erythroid progenitor"] <- "Erythroid Progenitors"

lct[lct == "Haematopoietic progenitors"] <- "Haematopoietic Progenitor"
lct[lct == "Haematopoietic progenitor"] <- "Haematopoietic Progenitor"
lct[lct == "Hematopoietic Progenitor"] <- "Haematopoietic Progenitor"
lct[lct == "Hematopoietic Progenitor Cells"] <- "Haematopoietic Progenitor"

lct[lct == "Macrophages"] <- "Macrophage"

lct[lct == "Embryonic fibroblast cells"] <- "Mouse Embryonic fibroblasts"

lct[lct == "Mouse ErythroLeukaemic "] <- "Mouse ErythroLeukaemic"

lct[lct == "Myeloid progenitor cells"] <- "Myeloid Progenitors"

lct[lct == "neural progenitor cells"] <- "Neural progenitor cells"

lct[lct == "Bone Marrow "] <- "Bone Marrow"

levels(ct) <- lct
table(ct)

codex$annotation[,"Cell type"] <- as.character(ct)
save(codex, file = "heatmap_shinyApp/data/codex_currated.RData")

