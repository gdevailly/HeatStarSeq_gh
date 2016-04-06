library(gplots)
library(cba)
myMat <- cor(cbind(runif(5), runif(5), runif(5)+5, runif(5)+5, runif(5)))
myLabs <- 1:5

d <- dist(myMat)

myClust <- hclust(d)
co <- order.optimal(d, myClust$merge)
myClust$merge <- co$merge
myClust$order <- co$order
dendro <- as.dendrogram(myClust)
clusterDat <- list("dend" = dendro, "order" = myClust$order)

heatmap(myMat,
        Rowv = clusterDat$dend,
        Colv = "Rowv",
        scale = "none",
        labRow = myLabs,
        labCol = myLabs,
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)

getParam <- function(x) {return(clusterDat$dend)}

heatmap(myMat[, clusterDat$order],
        Rowv = clusterDat$dend,
        Colv = getParam("both"),
        scale = "none",
        labRow = myLabs,
        labCol = myLabs,
        breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
        col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
        useRaster = TRUE
)

heatmap(myMat, scale = "none")
heatmap(myMat, Rowv = c(4,1,2,3,5), scale = "none", keep.dendro = TRUE)


heatmap.2(myMat, scale = "none", trace = "none", key = FALSE)
heatmap.2(myMat, scale = "none", trace = "none", key = FALSE)
