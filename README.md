# Heat\*Seq: A R/shiny web application to browse and compare public RNA-seq / ChIP-seq / CAGE datasets
by Guillaume Devailly, Anna Mantsoki, and Anagha Joshi

Contact: [@G_Devailly](https://twitter.com/G_Devailly) / guillaume.devailly _at_ roslin.ed.ac.uk

![HeatRNAseq screenshot](heatstarseq_screenshot.png)

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

**1) Formatting the dataset**

One needs to create a R list object, hereafter named `newDataset`, which contain the following elements (the element name matters, not the order of them in the list):
- `newDataset$dataMatrix`, a numeric matrix of one row per [GENCODE gene](http://www.gencodegenes.org/) and one column per sample. Each value should be a measure of gene expression (usually FPKM or TPM) for that gene in that sample. I would strongly advice to remove genes with no expression in the dataset `which(rowSums(newDataset$dataMatrix) != 0)`, to not name the rows and columns, and to replace all NAs by 0s.
- `newDataset$geneName`, a character vector of as many elements as there is rows in dataMatrix, containing the [GENCODE gene](http://www.gencodegenes.org/) name of each gene, in the same order as in the dataMatrix. Use GENCODE name without the number of transcripts (ENSG00000134046 and not ENSG00000134046.5).
- `newDataset$correlationMatrix`, the output of `cor(newDataset$dataMatrix)`:
```
newDataset$correlationMatrix <- cor(newDataset$dataMatrix)
```
- `newDataset$annotation`, a data.frame (with string as character, not factor) of one line per experiments in the datatset. The number of line of the annotation table must equal the number of columns of the dataMatrix. The annotation table MUST contain a `name` column containing UNIQUE character strings describing the experiment. It is strongly advised to add a `url` column storing a link to the original experiment, and to add one or more column for subseting the dataset, such as a cell type column or a library type column. Please, nullify the rownames of the annotation table before saving the object:
```
row.names(newDataset$annotation) <- NULL
```
- The, the `newDataset` object should be saved in the relevant directory:
```
save(newDataset, file = "heatrnaseq/data/newDataset.RData")
```

**2) Generating a newDataset_preload.RData**

Depending on the number of experiments, those dataset can be quite heavy, so we will load only one at a time on the shiny server. However, each dataset need to be preloaded. One can generate a `newDataset_perload.RData` by doing the following:
```
load("heatrnaseq/data/newDataset.RData")
nullifyDataForFasterPreloading <- function(myList) {
    myList[["dataMatrix"]] <- NULL
    myList[["geneName"]] <- NULL
    myList[["correlationMatrix"]] <- NULL
    return(myList)
}
object.size(newDataset)
newDataset <- nullifyDataForFasterPreloading(newDataset)
object.size(newDataset)
save(newDataset, file = "heatrnaseq/data/newDataset_preload.RData")
```
Which is what is done in [this script](dataset_fromating/8_data_preload_ChIPseq.R).

**3) Edit the shiny app to accept the new dataset**
- In the [gobal.R](application/heatrnaseq/global.R), add the following line at the end:
```
[...]
load("data/blueprint_rnaseq_preload.RData")
load("data/newDataset_preload.RData") # <- this line
```
- In the [ui.R](application/heatrnaseq/ui.R), modify the `selectInput("dataset",` function:
```
                    selectInput("dataset", label = NULL, choices = c(
                        "Bgee RNA-seq (human)",
                        "Blueprint RNA-seq (human)",
                        "ENCODE RNA-seq (human)",
                        "Bgee RNA-seq (mouse)",
                        "ENCODE RNA-seq (mouse)",
                        "the new dataset (relevant species)" # <- this line
                    )),
```
We will now setup the subseting widget for the new dataset. Subseting is done according to relevant column(s) of the `annotation` table. One can add subseting option for cell line/tissue, library preparation method, laboratory of origin, etc. Look for the line `# we addapt filtering widgets to the various datasets` and insert a new `div` function below the other filtering widgets for other datasets:
```
                    div(id = "widgetForNewDataset",
                        selectInput("tissue_newDataset", "Tissue of origin:",
                                    choices = unique(newDataset$annotation$tissue)[order(unique(newDataset$annotation$tissue))],
                                    selected = NULL, multiple = TRUE)
                    ),
```
Here we only add one subseting field, but one can include as many as she/he wants. Look for `div(id = "widgetForBgeeHuman",` to see how to include multiple filtering fields.
- In the [server.R](application/heatrnaseq/server.R), modify the second `observe({` function so that the filtering widget is displayed only when this dataset was selected:
```
observe({
        shinyjs::hide("widgetForBgeeHuman")
        shinyjs::hide("widgetForBlueprintHuman")
        shinyjs::hide("widgetForEncodeHuman")
        shinyjs::hide("widgetForEncodeMouse")
        shinyjs::hide("widgetForBgeeMouse")
        shinyjs::hide("widgetForNewDataset") # <- insert this line, same string as in the div() defined in the ui.R

        if (input$dataset == "ENCODE RNA-seq (human)") {
            shinyjs::show("widgetForEncodeHuman")
        } else if (input$dataset == "Bgee RNA-seq (human)") {
            shinyjs::show("widgetForBgeeHuman")
        } else if (input$dataset == "Blueprint RNA-seq (human)") {
            shinyjs::show("widgetForBlueprintHuman")
        } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
            shinyjs::show("widgetForEncodeMouse")
        } else if (input$dataset == "Bgee RNA-seq (mouse)") {
            shinyjs::show("widgetForBgeeMouse")
        } else if (input$dataset == "the new dataset (relevant species)") {  # <- insert this line, same sting as the dataset string defined in the ui.R
            shinyjs::show("widgetForNewDataset") # <- insert this line
        }
    })

```
Modify the `getSelectedDataset` reactive function, so that the dataset is loaded when the user select it:
```
    getSelectedDataset <- reactive({
        withProgress(value = 1, message = "Loading dataset: ", detail = "removing old dataset", {
            load("data/encode_rnaseq_preload.RData")
            load("data/bgee_human_preload.RData")
            load("data/encode_mouse_rnaseq_preload.RData")
            load("data/bgee_mouse_preload.RData")
            load("data/newDataset_preload.RData") # <- insert this line
            setProgress(value = 1, detail = "loading new dataset")
            if (input$dataset == "ENCODE RNA-seq (human)") {
                load("data/encode_rnaseq.RData")
                dataset <- encode_rnaseq
            } else if (input$dataset == "Bgee RNA-seq (human)") {
                load("data/bgee_human.RData")
                dataset <- bgee_human
            } else if (input$dataset == "Blueprint RNA-seq (human)") {
                [...]
            } else if (input$dataset == "the new dataset (relevant species)") {  # <- insert this line
                load("data/newDataset.RData")  # <- insert this line
                dataset <- newDataset  # <- insert this line
            }
            setProgress(value = 1, detail = "done!")
        })
        return(dataset)
    })
```
Finally, modify the tedious `subsetMatrix` reactive function, inserting the following at the bottom of the big chunk of `if else`:
```
        } else if (input$dataset == "the new dataset (relevant species)") {
            if (is.null(input$tissue_newDataset)) { # same input line as defined in the ui.R
                temp_tissue_newDataset <- unique(newDataset$annotation$tissue)
            } else {
                temp_tissue_newDataset <- input$cells_encode_m
            }
            keep <- which(
                dataset$annotation$tissue %in% temp_tissue_newDataset
            )
        }
```
Please, look into the `} else if (input$dataset == "Bgee RNA-seq (human)") {` section to see how to include multiple filtering conditions.

After debugging for typos, missing commas, parenthesis and brackets, it should work! Please, do push a merge request if you implemented a new dataset, or updated an old one!
