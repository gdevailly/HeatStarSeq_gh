library(dplyr)
library(cowplot)
theme_set(theme_cowplot() + theme(axis.line.x = element_line(), axis.line.y = element_line()))

peakSpace <- rep(FALSE, 500000)
refPeaks <- sample(length(peakSpace), 40000)
proportion <- seq(0.01, 1, by = 0.01) %>% rev
vectorRef <- peakSpace
vectorRef[refPeaks] <- TRUE

t0 <- Sys.time()
results <- lapply(proportion, function(x) {
    subsetPeaks <- sample(
        refPeaks,
        round(length(refPeaks) * x)
    )
    vectorSubset <- peakSpace
    vectorSubset[subsetPeaks] <- TRUE
    cc <- cor(vectorRef, vectorSubset)
    jaccard <- length(which(vectorSubset & vectorRef)) / length(which(vectorSubset | vectorRef))
    return(data.frame(
        prop = rep(x, 2),
        type = c("correlation", "jaccard"),
        value = c(cc, jaccard)
    ))
}) %>% bind_rows
Sys.time() - t0

p <- ggplot(results, aes(x = prop, y = value, col = type)) + geom_line() + scale_x_reverse() + labs(x = "relative subset size")
ggsave("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/Rplots/jaccard_vs_cor.png", p, width = 5, height = 4)
