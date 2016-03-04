###
# 2016-03-04
# files: 57epigenomes.RPKM.nc.gz  57epigenomes.RPKM.pc.gz  57epigenomes.RPKM.rb.gz  EG.name.txt  Ensembl_v65.Gencode_v10.ENSG.gene_info
# downloaded from: http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/
###
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/roadmap")
library(dplyr)
library(readr)

roadmap <- list(
    pc = read_tsv("57epigenomes.RPKM.pc"),
    nc = read_tsv("57epigenomes.RPKM.nc"),
    rb = read_tsv("57epigenomes.RPKM.rb")
)
# pbs seems to be due to tab at the end of lines

lapply(roadmap, dim)
lapply(roadmap, head)

roadmap$all <- do.call(rbind, roadmap)
roadmap$pc_nc <- rbind(roadmap$pc, roadmap$pcnc)

roadmap$all_mat <- roadmap$all[, -1] %>% as.matrix
roadmap$pc_nc_mat <-  roadmap$pc_nc[, -1] %>% as.matrix

t0 <- Sys.time()
roadmap$cor_all <- cor(roadmap$all_mat)
roadmap$cor_pc_nc <- cor(roadmap$pc_nc_mat)
Sys.time() - t0 # 0.5 sec

pdf(file = "heatmaps_roadmap.pdf")
heatmap(roadmap$cor_all,
        scale = "none",
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255)
)
heatmap(roadmap$cor_pc_nc,
        scale = "none",
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255)
)
dev.off()
system("firefox heatmaps_roadmap.pdf &")

geneName <- c(roadmap$pc$gene_id, roadmap$nc$gene_id, roadmap$rb$gene_id)
metadata <- read_tsv("EG.name.txt", col_names = FALSE)
metadata <- dplyr::rename(metadata, codeName = X1, name = X2)

length(metadata$name)
length(unique(metadata$name)) # \o/

metadata$url <- "http://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/"

roadmap_rna <- list(
    dataMatrix = roadmap$all_mat,
    geneName = geneName,
    correlationMatrix = roadmap$cor_all,
    annotation = metadata %>% as.data.frame(stringsAsFacor = FALSE)
)

colnames(roadmap_rna$dataMatrix) <- NULL
colnames(roadmap_rna$correlationMatrix) <- NULL
rownames(roadmap_rna$annotation) <- NULL
roadmap_rnaseq <- roadmap_rna
roadmap_rnaseq$annotation <- roadmap_rnaseq$annotation[-nrow(roadmap_rnaseq$annotation),]

save(roadmap_rnaseq, file ="../../heatrnaseq/data/roadmap_rnaseq.RData")




