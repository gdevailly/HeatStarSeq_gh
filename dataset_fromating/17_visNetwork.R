library(dplyr)
library(networkD3)
library(reshape2)

setwd("otherProject/ChIP_heatmap/heatchipseq/")
load("data/encode.RData")

myThreshold <- 0.8

myEdges <- encode$correlationMatrix %>% melt
myEdges[(myEdges$Var2 > myEdges$Var1) & myEdges$value >= myThreshold, ]

t0 <- Sys.time()

simpleNetwork(
    myEdges
)

Sys.time() - t0

# the execution time is fast, but the display takes ages... DO NOT USE networkD3 !!!



#####
library(dplyr)
library(visNetwork)
library(reshape2)

setwd("otherProject/ChIP_heatmap/heatchipseq/")
load("data/encode.RData")

myThreshold <- 0.8

myEdges <- encode$correlationMatrix %>% melt
myEdges[(myEdges$Var2 > myEdges$Var1) & myEdges$value >= myThreshold, ]

t0 <- Sys.time()
visNetwork(
    data.frame(
        id = 1:nrow(encode$correlationMatrix)
    ),
    data.frame(
        from = myEdges$Var1,
        to = myEdges$Var2
    ),
    width = "100%"
)
Sys.time() - t0

# the execution time is fast, but the display takes ages... DO NOT USE visNetwork !!!

