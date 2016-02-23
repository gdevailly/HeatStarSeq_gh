# Heat\*Seq: A R/shiny web application to browse and compare public RNA-seq / ChIP-seq / CAGE datasets
by Guillaume Devailly, Anna Mantsoki, and Anagha Joshi

Contact: [@G_Devailly](https://twitter.com/G_Devailly) / guillaume.devailly _at_ roslin.ed.ac.uk

## Summary

###### How to use Heat\*Seq
- The public server
- Running Heat\*Seq locally
- Creating your Heat\*Seq server

###### How to add new datasets
- Gene expression data
- ChIP-seq data

## How to use Heat\*Seq

### The public server
Heat\*Seq is available for now [at this address](http://www.chipcompare.roslin.ed.ac.uk/).
**The Application is still in early development.**
**The url of the public server WILL change in the near future.**

### Running Heat\*Seq locally
Download the Github folder (for example, from [here](https://github.com/gdevailly/HeatStarSeq_gh/archive/master.zip)). Extract the  .zip archive. You will need R (at least 3.2), and need to install several R packages from CRAN and Bioconductor (package list to come soon). Launch R, and goes to one of the two following directories `application/heatrnaseq/` or `application/heatchipseq/`, using for example the `setwd()` R command. Finally, execute the following:
```
library(shiny)
runApp()
```
Some datasets are heavy, so you may need more than 3 Gb of available memory.
### Creating your Heat\*Seq server
Heat\*Seq is a [Shiny application](http://shiny.rstudio.com/), so you will need a **Shiny server**.
You can install the [Open source edition of Shiny server](https://github.com/rstudio/shiny-server) on any compatible web-server, or use [shinyapps.io](http://www.shinyapps.io/).

Heat\*Seq is composed of three independent Shiny Apps, [HeatStarSeq](https://github.com/gdevailly/HeatStarSeq_gh/tree/master/application/heatstarseq), [HeatRNAseq](https://github.com/gdevailly/HeatStarSeq_gh/tree/master/application/heatrnaseq) and [HeatChIPseq](https://github.com/gdevailly/HeatStarSeq_gh/tree/master/application/heatchipseq). You will need to edit the [server_adresses.R](https://github.com/gdevailly/HeatStarSeq_gh/blob/master/application/heatstarseq/data/server_adresses.R) script with the proper web url:
```
URL_HEATSTARSEQ <- "http://www.chipcompare.roslin.ed.ac.uk"
URL_HEATRNASEQ <- "http://www.chipcompare.roslin.ed.ac.uk/heatrnaseq"
URL_HEATCHIPSEQ <- "http://www.chipcompare.roslin.ed.ac.uk/heatchipseq"
```
and copy it in the `data` folder of **each** application.

If you created a mirror of Heat\*Seq, I will be very pleased if you contact me so that I can advertised it.

## How to add new datasets

### Gene expression data
