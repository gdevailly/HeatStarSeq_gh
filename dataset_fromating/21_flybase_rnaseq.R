# flybase rnaseq
# file downloaded from ftp://ftp.flybase.net/releases/current/precomputed_files/genes/ 2016-03-15
setwd("/groups2/joshi_grp/guillaume/otherProject/ChIP_heatmap/data/flybase_rnaseq")
library(dplyr)
library(readr)


fbdata <- read_tsv("gene_rpkm_report_fb_2016_01.tsv", skip = 5)
# end of file is messy. Seriously flybase?

names(fbdata) <- c("Release_ID", "FBgn", "GeneSymbol", "Parent_library_FBlc", "Parent_library_name", "RNASource_FBlc", "RNASource_name",
                 "RPKM_value", "Bin_value", "Unique_exon_base_count", "Total_exon_base_count", "Count_used")
fbdata <- filter(fbdata, !is.na(FBgn))

length(unique(fbdata$RNASource_name)) # 124
table(fbdata$RNASource_name) # messy
table(fbdata$RNASource_name) %>% range # 17527 18056
length(unique(fbdata$GeneSymbol)) # 17943
length(unique(fbdata$FBgn)) # 18057, this will be fun

uniqueGenes <- unique(fbdata$FBgn)[order(unique(fbdata$FBgn))]

library(tidyr)
tidyrMagix <- spread(fbdata, RNASource_name, RPKM_value) # 1 minute ?

dataMatrix <- matrix(0, nrow = length(uniqueGenes), ncol = length(unique(fbdata$RNASource_name)))

