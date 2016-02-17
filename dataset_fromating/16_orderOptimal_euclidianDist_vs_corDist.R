setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/heatrnaseq/")

library(dplyr)
library(gplots)
library(svglite)
library(cba)

load("data/bgee_mouse.RData")

colnames(bgee_mouse$correlationMatrix) <- rownames(bgee_mouse$correlationMatrix) <- bgee_mouse$geneName

oldClust <- hclust(dist(bgee_mouse$correlationMatrix)) %>% as.dendrogram
newClust <- hclust(as.dist(1 - bgee_mouse$correlationMatrix)) %>% as.dendrogram

d <- as.dist(1 - bgee_mouse$correlationMatrix)
hc <- hclust(d)
co <- order.optimal(d, hc$merge)
ho <- hc
ho$merge <- co$merge
ho$order <- co$order
newOoClust <- as.dendrogram(ho)

d <- dist(bgee_mouse$correlationMatrix)
hc <- hclust(d)
co <- order.optimal(d, hc$merge)
ho <- hc
ho$merge <- co$merge
ho$order <- co$order
oldOoClust <- as.dendrogram(ho)


svglite("../Rplots/compHM_1.svg")
heatmap(bgee_mouse$correlationMatrix,
        Rowv = oldClust,
        Colv = "Rowv",
        scale = "none",
        margins = rep(5, 2),
        cexRow = 0.4,
        cexCol = 0.4,
        breaks = seq(-0.5, 3, length.out = 256),
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()
svglite("../Rplots/compHM_2.svg")
heatmap(bgee_mouse$correlationMatrix,
        Rowv = newClust,
        Colv = "Rowv",
        scale = "none",
        margins = rep(5, 2),
        cexRow = 0.4,
        cexCol = 0.4,
        breaks = seq(-0.5, 3, length.out = 256),
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()
svglite("../Rplots/compHM_3.svg")
heatmap(bgee_mouse$correlationMatrix,
        Rowv = newOoClust,
        Colv = "Rowv",
        scale = "none",
        margins = rep(5, 2),
        cexRow = 0.4,
        cexCol = 0.4,
        breaks = seq(-0.5, 3, length.out = 256),
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()
svglite("../Rplots/compHM_4.svg")
heatmap(bgee_mouse$correlationMatrix,
        Rowv = oldOoClust,
        Colv = "Rowv",
        scale = "none",
        margins = rep(5, 2),
        cexRow = 0.4,
        cexCol = 0.4,
        breaks = seq(-0.5, 3, length.out = 256),
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)
dev.off()






