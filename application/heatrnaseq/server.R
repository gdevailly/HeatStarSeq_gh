library(readr)
library(svglite)
library(cba)

options(shiny.maxRequestSize = 10*1024^2) # max file size, 10Mb

shinyServer(function(input, output, session) {

    # UI elements activations -------------------
    observe({
        junkVar <- input$fileFormatInstructions
        shinyjs::toggle("div_fileFormatInstructions")
    })

    observe({
        if(input$fileToUse == "Upload your expression file") {
            shinyjs::show("div_fileupload")
            shinyjs::hide("div_exampleInUse")
        } else {
            shinyjs::hide("div_fileupload")
            shinyjs::show("div_exampleInUse")
        }
    })

    observe({
        junkVar <- input$advClustOptions
        shinyjs::toggle("widgetForClustOptions")
    })

    observe({
        junkVar <- input$coloursOptions
        shinyjs::toggle("div_colourOptions")
    })

    observe({
        if(input$correlationCorrection == "Linear scaling") {
            shinyjs::show("div_maxCorrelation")
        } else {
            shinyjs::hide("div_maxCorrelation")
        }
    })

    observe({
        if (is.null(input$expressionFile) & input$fileToUse != "Use the example file") {
            shinyjs::hide("downloadUserExpressionFile")
            shinyjs::hide("downloadUserCorrelationTable")
        } else {
            shinyjs::show("downloadUserExpressionFile")
            shinyjs::show("downloadUserCorrelationTable")
        }
    })

    observe({
        shinyjs::hide("widgetForBgeeHuman")
        shinyjs::hide("widgetForBlueprintHuman")
        shinyjs::hide("widgetForRoadmapHuman")
        shinyjs::hide("widgetForGtexSmall")
        shinyjs::hide("widgetForGtexLarge")
        shinyjs::hide("widgetForEncodeHuman")
        shinyjs::hide("widgetForEncodeMouse")
        shinyjs::hide("widgetForBgeeMouse")
        shinyjs::hide("widgetForFlybase")

        if (input$dataset == "ENCODE RNA-seq (human)") {
            shinyjs::show("widgetForEncodeHuman")
        } else if (input$dataset == "Bgee RNA-seq (human)") {
            shinyjs::show("widgetForBgeeHuman")
        } else if (input$dataset == "Blueprint RNA-seq (human)") {
            shinyjs::show("widgetForBlueprintHuman")
        } else if (input$dataset == "Roadmap Epigenomics RNA-seq (human)") {
            shinyjs::show("widgetForRoadmapHuman")
        } else if (input$dataset == "GTEx summary (human)") {
            shinyjs::show("widgetForGtexSmall")
        } else if (input$dataset == "GTEx - all samples (human)") {
            shinyjs::show("widgetForGtexLarge")
        } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
            shinyjs::show("widgetForEncodeMouse")
        } else if (input$dataset == "Bgee RNA-seq (mouse)") {
            shinyjs::show("widgetForBgeeMouse")
        } else if (input$dataset == "Flybase RNA-seq (drosophila)") {
            shinyjs::show("widgetForFlybase")
        } else if (input$dataset == "other (soon)") {
            # fill this
        }
    })

    observe({
        if(input$myPanels == "Static heatmap" | input$myPanels == "Tree") {
            shinyjs::show("widgetForLabels")
        } else {
            shinyjs::hide("widgetForLabels")
        }
        if(input$myPanels == "Static heatmap") {
            shinyjs::show("div_widgetHMoptions")
        } else {
            shinyjs::hide("div_widgetHMoptions")
        }
    })

    observe({
        if(input$labelOption == "Automatic") {
            shinyjs::hide("widgetForLabelsManual")
        } else {
            shinyjs::show("widgetForLabelsManual")
        }
    })

    observe({
        # range of colour slider should be precisely defined to keep ordering:
        updateSliderInput(session, "col_val1", max = min(max(input$col_val2 - 0.05, -0.95), 0.7))
        updateSliderInput(session, "col_val2", min = max(input$col_val1 + 0.05, -0.9), max = min(input$col_val3 - 0.05, 0.8))
        updateSliderInput(session, "col_val3", min = max(input$col_val2 + 0.05, -0.8), max = min(input$col_val4 - 0.05, 0.9))
        updateSliderInput(session, "col_val4", min = max(input$col_val3 + 0.05, -0.7))
    })

    observe({
        if(input$myPanels == "Pairwise plot") {
            shinyjs::show("widgetForPairwisePlots")
            shinyjs::hide("widgetsForHeatmap")
        } else {
            shinyjs::hide("widgetForPairwisePlots")
            shinyjs::show("widgetsForHeatmap")
        }
    })

    # output computations ----------------
    getSelectedDataset <- reactive({
        withProgress(value = 1, message = "Loading dataset: ", detail = "removing old dataset", {
            load("data/encode_rnaseq_preload.RData")
            load("data/bgee_human_preload.RData")
            load("data/encode_mouse_rnaseq_preload.RData")
            load("data/bgee_mouse_preload.RData")
            load("data/blueprint_rnaseq_preload.RData")
            load("data/roadmap_rnaseq_preload.RData")
            load("data/gtex_small_preload.RData")
            load("data/gtex_large_preload.RData")
            load("data/flybase_rnaseq_preload.RData")
            setProgress(value = 1, detail = "loading new dataset")
            if (input$dataset == "ENCODE RNA-seq (human)") {
                load("data/encode_rnaseq.RData")
                dataset <- encode_rnaseq
            } else if (input$dataset == "Bgee RNA-seq (human)") {
                load("data/bgee_human.RData")
                dataset <- bgee_human
            } else if (input$dataset == "Blueprint RNA-seq (human)") {
                load("data/blueprint_rnaseq.RData")
                dataset <- blueprint_rnaseq
            } else if (input$dataset == "Roadmap Epigenomics RNA-seq (human)") {
                load("data/roadmap_rnaseq.RData")
                dataset <- roadmap_rnaseq
            } else if (input$dataset == "GTEx summary (human)") {
                load("data/gtex_small.RData")
                dataset <- gtex_small
            } else if (input$dataset == "GTEx - all samples (human)") {
                load("data/gtex_large.RData")
                dataset <- gtex_large
            } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
                load("data/encode_mouse_rnaseq.RData")
                dataset <- encode_mouse_rnaseq
            } else if (input$dataset == "Bgee RNA-seq (mouse)") {
                load("data/bgee_mouse.RData")
                dataset <- bgee_mouse
            } else if (input$dataset == "Flybase RNA-seq (drosophila)") {
                load("data/flybase_rnaseq.RData")
                dataset <- flybase_rnaseq
            } else if (input$dataset == "other (soon)") {
                # fill this
            }
            setProgress(value = 1, detail = "done!")
        })
        isolate(
            if (input$fileToUse == "Use the example file") {
                if (!input$dataset %in% c("Bgee RNA-seq (mouse)", "ENCODE RNA-seq (mouse)")) {
                    updateRadioButtons(session = session, "fileToUse", label = NULL, choices = NULL, selected = "Upload your expression file")
                }
            }
        )
        return(dataset)
    })

    userExpressionFileAnalysis_1 <- reactive({

        if (input$fileToUse == "Use the example file") {
            withProgress(value = 1, message = "Loading example: ", detail = "reading file", {
                userExpressionFile_temp <- read_tsv("www/rnaseq_mouse_GSE70732_NP1071_FPKM_all_genes_with_header.txt", col_names = TRUE)[, 1:2]
                colnames(userExpressionFile_temp) <- c("ensembl_id","exp_value")
                setProgress(value = 1, detail = "intersecting gene names")
                userExpressionFile_temp_v <- userExpressionFile_temp$exp_value
                names(userExpressionFile_temp_v) <- userExpressionFile_temp$ensembl_id
                dataset <- getSelectedDataset()
                validate(need(
                    input$dataset %in% c("Bgee RNA-seq (mouse)", "ENCODE RNA-seq (mouse)"),
                    "The example file is a mouse RNA-seq experiment. Please choose a mouse dataset or unload the example."
                ))
                userExpressionFile <- userExpressionFile_temp_v[dataset$geneName]
                userExpressionFile[is.na(userExpressionFile)] <- 0
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(log10(userExpressionFile + 1), dataset$dataMatrix) %>% as.vector
                userExpressionFile <- data.frame(
                    geneName = dataset$geneName,
                    value = userExpressionFile,
                    stringsAsFactors = FALSE
                )
                setProgress(value = 1, detail = "done!")
            })
            return(list(
                "expression" = userExpressionFile,
                "correlations" = userCorrelations,
                "linearNormCorrelations" = NULL # to be fill below
            ))
        }

        userExpressionFileName <- input$expressionFile
        if (is.null(userExpressionFileName)){
            userExpressionFile <- NULL
            userCorrelations <- NULL
        } else {
            withProgress(value = 1, message = "User expression file: ", detail = "reading file", {
                userExpressionFile_temp <- read_tsv(userExpressionFileName$datapath, col_names = input$header)[, 1:2]
                colnames(userExpressionFile_temp) <- c("ensembl_id","exp_value")
                setProgress(value = 1, detail = "intersecting gene names")
                userExpressionFile_temp_v <- userExpressionFile_temp$exp_value
                names(userExpressionFile_temp_v) <- userExpressionFile_temp$ensembl_id
                dataset <- getSelectedDataset()
                userExpressionFile <- userExpressionFile_temp_v[dataset$geneName]
                userExpressionFile[is.na(userExpressionFile)] <- 0
                validate(need(
                    length(unique(userExpressionFile)) > 1,
                    "Something went wrong, sorry. :(
                    \nUnable to compute correlations.
                    \nThis is likely because:
                    \n-the wrong dataset is selected (ie you may have uploaded a mouse expression file, but the selected dataset is for human).
                    \n-the file formatting is not recognized. It should be a two columns, tab delimited text file. First columns should contain the ensembl gene name id (ie ENSG00000000003) or Flybase id (FBgn0000003). Second column should contains expression values.
                    \nPlease contact us if you think you have found a bug."
                ))
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(log10(userExpressionFile + 1), dataset$dataMatrix) %>% as.vector
                userExpressionFile <- data.frame(
                    geneName = dataset$geneName,
                    value = userExpressionFile,
                    stringsAsFactors = FALSE
                )
                setProgress(value = 1, detail = "done!")
            })
        }
        return(list(
            "expression" = userExpressionFile,
            "correlations" = userCorrelations,
            "linearNormCorrelations" = NULL # to be fill below
        ))
    })

    userExpressionFileAnalysis <- reactive({
        userExpressionFileData <- userExpressionFileAnalysis_1()
        userCorrelations <- userExpressionFileData$correlations
        if (!is.null(input$expressionFile) | input$fileToUse == "Use the example file") {
            if (input$maxCorrelation != 0) {
                userExpressionFileData$linearNormCorrelations <- userCorrelations * input$maxCorrelation / max(userCorrelations)
            }
        }
        return(userExpressionFileData)
    })

    getExperimentList <- reactive({
        datasetSampleNames <- sort(getSelectedDataset()$annotation$name)
        if (!is.null(input$expressionFile) | input$fileToUse == "Use the example file") {
            datasetSampleNames <- c(input$nameOfExpressionFile, datasetSampleNames)
        }
        return(datasetSampleNames)
    })

    output$scatterPlotSample1 <- renderUI({
        withProgress(value = 1, message = "Creating list of experiments 1/2", {
            selectInput("scatterPlotSample1Defined", "Experiment 1:", getExperimentList(), selected = getExperimentList()[1])
        })
    })

    output$scatterPlotSample2 <- renderUI({
        withProgress(value = 1, message = "Creating list of experiments 2/2", {
            selectInput("scatterPlotSample2Defined", "Experiment 2:", getExperimentList(), selected = getExperimentList()[2])
        })
    })

    output$tabUserExpressionFile <- renderDataTable({
        validate(
            need(!is.null(input$expressionFile) | input$fileToUse == "Use the example file", "Upload an expression file, or click on a 'heatmap' tab to explore the dataset.")
        )
        userExpressionFileAnalysis()$expression
    })

    output$downloadUserExpressionFile <- downloadHandler("uploaded_expression_file.txt",
                                                         content = function(file) {
                                                             write.table(
                                                                 userExpressionFileAnalysis()$expression,
                                                                 file = file, row.names = FALSE, quote = FALSE, sep = "\t"
                                                             )
                                                         },
                                                         contentType = "text/tsv")

    getCorrelationTable <- reactive({
        validate(
            need(!is.null(input$expressionFile) | input$fileToUse == "Use the example file", "Upload an expression file, or click on a 'heatmap' tab to explore the dataset.")
        )
        if(input$correlationCorrection == "Linear scaling") {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userExpressionFileAnalysis()$correlations,
                "scaledCorrelation" = userExpressionFileAnalysis()$linearNormCorrelations,
                stringsAsFactors = FALSE
            )[order(userExpressionFileAnalysis()$correlations, decreasing = TRUE), ])
        } else {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userExpressionFileAnalysis()$correlations,
                stringsAsFactors = FALSE
            )[order(userExpressionFileAnalysis()$correlations, decreasing = TRUE), ])
        }
    })

    output$tabUserCorrelationTable <- renderDataTable(getCorrelationTable())

    output$downloadUserCorrelationTable <- downloadHandler("heatRNAseq_correlations.txt",
                                                           content = function(file) {
                                                              write.table(getCorrelationTable(), file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                                           },
                                                           contentType = "text/tsv")
    # functions and reactives ------------------
    subsetMatrix <- reactive({
        dataset <- getSelectedDataset()
        workingMatrix <- dataset$correlationMatrix
        keep<- 1:nrow(workingMatrix)
        # filtering is dataset-dependent...
        if (input$dataset == "ENCODE RNA-seq (human)") {
            if (is.null(input$cells)) {
                temp_cells <- unique(encode_rnaseq$annotation$biosampleTermname)
            } else {
                temp_cells <- input$cells
            }
            if (is.null(input$sampleType)) {
                temp_sampleType <- unique(encode_rnaseq$annotation$biosampletype)
            } else {
                temp_sampleType <- input$sampleType
            }
            if (is.null(input$rnaExtract)) {
                temp_rnaExtract <- unique(encode_rnaseq$annotation$rnaFraction)
            } else {
                temp_rnaExtract <- input$rnaExtract
            }
            keep <- which(
                dataset$annotation$biosampleTermname %in% temp_cells &
                dataset$annotation$biosampletype %in% temp_sampleType &
                dataset$annotation$rnaFraction %in% temp_rnaExtract
            )
        } else if (input$dataset == "Bgee RNA-seq (human)") {
            if (is.null(input$tissus_bgee_h)) {
                temp_tissus_bgee_h <- unique(bgee_human$annotation$tissue)
            } else {
                temp_tissus_bgee_h <- input$tissus_bgee_h
            }
            if (is.null(input$dvp_bgee_h)) {
                temp_dvp_bgee_h <- unique(bgee_human$annotation$stage)
            } else {
                temp_dvp_bgee_h <- input$dvp_bgee_h
            }
            if (is.null(input$library_bgee_h)) {
                temp_library_bgee_h <- unique(bgee_human$annotation$libraryType)
            } else {
                temp_library_bgee_h <- input$library_bgee_h
            }
            keep <- which(
                dataset$annotation$tissue %in% temp_tissus_bgee_h &
                dataset$annotation$stage %in% temp_dvp_bgee_h &
                dataset$annotation$libraryType %in% temp_library_bgee_h
            )
        } else if (input$dataset == "Blueprint RNA-seq (human)") {
            if (is.null(input$tissue_blueprint_h)) {
                temp_tissue_blueprint_h <- unique(blueprint_rnaseq$annotation$sampleSource)
            } else {
                temp_tissue_blueprint_h <- input$tissue_blueprint_h
            }
            if (is.null(input$celltype_blueprint_h)) {
                temp_celltype_blueprint_h <- unique(blueprint_rnaseq$annotation$cellType)
            } else {
                temp_celltype_blueprint_h <- input$celltype_blueprint_h
            }
            keep <- which(
                dataset$annotation$sampleSource %in% temp_tissue_blueprint_h &
                dataset$annotation$cellType %in% temp_celltype_blueprint_h
            )
        } else if (input$dataset == "Roadmap Epigenomics RNA-seq (human)") {
            if (is.null(input$celltype_roadmap_h)) {
                temp_celltype_roadmap_h <- unique(roadmap_rnaseq$annotation$name)
            } else {
                temp_celltype_roadmap_h <- input$celltype_roadmap_h
            }
            keep <- which(
                dataset$annotation$name %in% temp_celltype_roadmap_h
            )
        } else if (input$dataset == "GTEx summary (human)") {
            if (is.null(input$celltype_gtex_small_h)) {
                temp_celltype_gtex_small_h <- unique(gtex_small$annotation$name)
            } else {
                temp_celltype_gtex_small_h <- input$celltype_gtex_small_h
            }
            keep <- which(
                dataset$annotation$name %in% temp_celltype_gtex_small_h
            )
        } else if (input$dataset == "GTEx - all samples (human)") {
            if (is.null(input$celltype_gtex_large_h)) {
                temp_celltype_gtex_large_h <- unique(gtex_large$annotation$SMTSD)
            } else {
                temp_celltype_gtex_large_h <- input$celltype_gtex_large_h
            }
            keep <- which(
                dataset$annotation$SMTSD %in% temp_celltype_gtex_large_h
            )
        } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
            if (is.null(input$cells_encode_m)) {
                temp_cells_encode_m <- unique(encode_mouse_rnaseq$annotation$tissue)
            } else {
                temp_cells_encode_m <- input$cells_encode_m
            }
            if (is.null(input$sampleType_encode_m)) {
                temp_sampleType_encode_m <- unique(encode_mouse_rnaseq$annotation$sampleType)
            } else {
                temp_sampleType_encode_m <- input$sampleType_encode_m
            }
            keep <- which(
                dataset$annotation$tissue %in% temp_cells_encode_m &
                dataset$annotation$sampleType %in% temp_sampleType_encode_m
            )
        } else if (input$dataset == "Bgee RNA-seq (mouse)") {
            if (is.null(input$tissus_bgee_m)) {
                temp_tissus_bgee_m <- unique(bgee_mouse$annotation$tissue)
            } else {
                temp_tissus_bgee_m <- input$tissus_bgee_m
            }
            if (is.null(input$dvp_bgee_m)) {
                temp_dvp_bgee_m <- unique(bgee_mouse$annotation$stage)
            } else {
                temp_dvp_bgee_m <- input$dvp_bgee_m
            }
            if (is.null(input$library_bgee_m)) {
                temp_library_bgee_m <- unique(bgee_mouse$annotation$libraryType)
            } else {
                temp_library_bgee_m <- input$library_bgee_m
            }
            keep <- which(
                dataset$annotation$tissue %in% temp_tissus_bgee_m &
                dataset$annotation$stage %in% temp_dvp_bgee_m &
                dataset$annotation$libraryType %in% temp_library_bgee_m
            )
        } else if (input$dataset == "Flybase RNA-seq (drosophila)") {
            if (is.null(input$sample_flybase)) {
                temp_sample_flybase <- unique(flybase_rnaseq$annotation$name)
            } else {
                temp_sample_flybase <- input$sample_flybase
            }
            if (is.null(input$library_flybase)) {
                temp_library_flybase <- unique(flybase_rnaseq$annotation$parentLibrary)
            } else {
                temp_library_flybase <- input$library_flybase
            }
            keep <- which(
                dataset$annotation$name %in% temp_sample_flybase &
                dataset$annotation$parentLibrary %in% temp_library_flybase
            )
        } else { # modify accordingly for new datasets!

        }

        myLabels <- dataset$annotation$name
        validate(
            need(length(keep) >= 3, "Less than 3 experiments match your criteria. Please select more experiments.")
        )
        validate(
            need(length(keep) <= 1100, "More than 1100 experiments match your criteria. Please select fewer experiments.")
        )
        if(length(keep) >= 3) {
            workingMatrix <- workingMatrix[keep, keep]
            myLabels <- myLabels[keep]
        }

        # merging user data in correlation matrix
        if(!is.null(userExpressionFileAnalysis()$correlations)) {
            userCorrelations <- userExpressionFileAnalysis()$correlations
            if (input$correlationCorrection == "Linear scaling") {
                userCorrelations <- userExpressionFileAnalysis()$linearNormCorrelations
            }
            workingMatrix <- rbind(workingMatrix, userCorrelations[keep])
            workingMatrix <- cbind(workingMatrix, c(userCorrelations[keep], 1))
            myLabels <- c(myLabels, input$nameOfExpressionFile)
        }
        return(list("mat" = workingMatrix, "myLabels" = myLabels))
    })

    doTheClustering <- reactive({
        withProgress(value = 1, message = "Clustering: ", detail = "hierarchical clustering", {
            myData <- subsetMatrix()
            myMat <- myData$mat
            colnames(myMat) <- rownames(myMat) <- myData$myLabels
            if (input$distOption == "1 - Pearson's correlation coefficient") {
                d <- as.dist(1 - myMat)
            } else {
                d <- dist(myMat, method = input$distOption)
            }
            myClust <- hclust(d, method = input$hclustMethod)
            co <- order.optimal(d, myClust$merge)
            myClust$merge <- co$merge
            myClust$order <- co$order
            dendro <- as.dendrogram(myClust)
            setProgress(value = 1, detail = "done!")
        })
        return(list("dend" = dendro, "order" = myClust$order))
    })

    matAfterHighlight <- reactive({
        matAA <- subsetMatrix()
        if (input$highlight & (!is.null(input$expressionFile) | input$fileToUse == "Use the example file")) {
            # we increase user cor value by 2 to do the color trick
            matAA$mat[, ncol(matAA$mat)] <- matAA$mat[, ncol(matAA$mat)] + 2
            matAA$mat[nrow(matAA$mat), 1:(ncol(matAA$mat) - 1)] <- matAA$mat[nrow(matAA$mat), 1:(ncol(matAA$mat) - 1)] + 2
        }
        return(matAA)
    })

    getColourPalette <- reactive({
        junkVar <- input$applyColoursOptions
        isolate({
            myCols <- c(input$col_col1, input$col_col2, input$col_col3, input$col_col4)
            myBreaks <- c(input$col_val1, input$col_val2, input$col_val3, input$col_val4)
        })
        return(list("myCols" = myCols, "myBreaks" = myBreaks))
    })

    # static heatmap ------------------------
    myRenderPlot <- function() {
        matData <- matAfterHighlight()
        clusterDat <- doTheClustering()

        nSample <- length(matData$myLabels)

        newLabSize <- round((6/nSample^0.6)*10)/10
        if(newLabSize > 3) newLabSize <- 3
        if(newLabSize < 0.1) newLabSize <- 0.1
        if (input$labelOption == "Automatic") updateSliderInput(session, "labCex", label = NULL, value = newLabSize, min = NULL, max = NULL, step = NULL)

        newMargin <- round((66/nSample^0.25)*10)/10
        if(newMargin > 50) newMargin <- 50
        if(newMargin < 1) newMargin <- 1
        if (input$labelOption == "Automatic") updateSliderInput(session, "margin", label = NULL, value = newMargin, min = NULL, max = NULL, step = NULL)

        if (input$showDend == "both") {
            myMat <- matData$mat
            Rowv <- clusterDat$dend
            Colv <- "Rowv"
            labRow <- if (input$showLabels %in% c("both", "row")) matData$myLabels else NA
            labCol <-  if (input$showLabels %in% c("both", "column")) matData$myLabels else NA
        } else if (input$showDend == "row") {
            myMat <- matData$mat[, clusterDat$order]
            Rowv <- rev(clusterDat$dend)
            Colv <- NA
            labRow <- if (input$showLabels %in% c("both", "row")) matData$myLabels else NA
            labCol <-  if (input$showLabels %in% c("both", "column")) matData$myLabels[clusterDat$order] else NA
        } else if (input$showDend == "column") {
            myMat <- matData$mat[rev(clusterDat$order), ]
            Rowv <- NA
            Colv <- clusterDat$dend
            labRow <- if (input$showLabels %in% c("both", "row")) matData$myLabels[rev(clusterDat$order)] else NA
            labCol <-  if (input$showLabels %in% c("both", "column")) matData$myLabels else NA
        } else if (input$showDend == "none") {
            myMat <- matData$mat[rev(clusterDat$order), clusterDat$order]
            Rowv <- NA
            Colv <- NA
            labRow <- if (input$showLabels %in% c("both", "row")) matData$myLabels[rev(clusterDat$order)] else NA
            labCol <-  if (input$showLabels %in% c("both", "column")) matData$myLabels[clusterDat$order] else NA
        }

        myPalette <- getColourPalette()

        heatmap(myMat,
                Rowv = Rowv,
                Colv = Colv,
                scale = "none",
                labRow = labRow,
                labCol = labCol,
                margins = if(input$labelOption == "Automatic") rep(newMargin, 2) else  rep(input$margin, 2),
                cexRow = if(input$labelOption == "Automatic") newLabSize else input$labCex,
                cexCol = if(input$labelOption == "Automatic") newLabSize else input$labCex,
                breaks = unique(c(
                    -1,
                    seq(myPalette$myBreaks[1], myPalette$myBreaks[2], length.out = 40),
                    seq(myPalette$myBreaks[2], myPalette$myBreaks[3], length.out = 41),
                    seq(myPalette$myBreaks[3], myPalette$myBreaks[4], length.out = 41),
                    seq(myPalette$myBreaks[4], 1.01, length.out = 2),
                    seq(1.01, 3, length.out = 121)
                )),
                col = c(myPalette$myCols[1], colorRampPalette(myPalette$myCols)(119), myPalette$myCols[4], colorRampPalette(c("blue", "yellow", "green", "black"))(120)),
                useRaster = TRUE
        )
    }

    output$downloadHMpng <- downloadHandler("heatmap.png",
                                            content = function(file) {
                                                png(file, width = 950, height = 950)
                                                myRenderPlot()
                                                dev.off()
                                            },
                                            contentType = "image/png")

    output$downloadHMpdf <- downloadHandler("heatmap.pdf",
                                            content = function(file) {
                                                pdf(file, width = 13.85, height = 13.85)
                                                myRenderPlot()
                                                dev.off()
                                            },
                                            contentType = "image/pdf")

    output$downloadHMsvg <- downloadHandler("heatmap.svg",
                                            content = function(file) {
                                                svglite(file, width = 13.85, height = 13.85)
                                                myRenderPlot()
                                                dev.off()
                                            },
                                            contentType = "image/svg")

    output$downloadHMdata <- downloadHandler("clusteredMatrix.txt",
                                            content = function(file) {
                                                matData <- subsetMatrix()
                                                clusterDat <- doTheClustering()
                                                clusteredMatrix <- matData$mat[clusterDat$order, clusterDat$order]
                                                colnames(clusteredMatrix) <- rownames(clusteredMatrix) <- matData$myLabels[clusterDat$order]
                                                write.table(clusteredMatrix, file = file, quote = FALSE, sep = "\t")
                                            },
                                            contentType = "text/tsv")

    output$myHeatmap <- renderPlot({
        withProgress(value = 1, message = "Ploting...", {
            myRenderPlot()
            setProgress(value = 1, detail = "done!")
        })
    })

    renderColourKey <- function() {
        myPalette <- getColourPalette()
        breaks <- unique(c(
            seq(myPalette$myBreaks[1], myPalette$myBreaks[2], length.out = 40),
            seq(myPalette$myBreaks[2], myPalette$myBreaks[3], length.out = 41),
            seq(myPalette$myBreaks[3], myPalette$myBreaks[4], length.out = 41)
        ))
        col <- colorRampPalette(myPalette$myCols)(119)
        oldPar <- list(mar = par()$mar, cex = par()$cex, mgp = par()$mgp)
        par(mar = c(1.5, 0.8, 1.5, 0.8), cex = 1.4, mgp = c(3, 0.4, 0))
        image(matrix(breaks), col = col, axes = FALSE)
        axis(1, at = c(0, 1/3, 2/3, 1), labels = myPalette$myBreaks)
        title(main = "Correlation coefficient", line = 0.2)
        box()
        par(oldPar)
    }

    output$colourKey1 <- renderPlot(renderColourKey())

    # responsive heatmap -----------------------
    output$myPlotlyHeatmap <- renderPlotly({
        matData <- matAfterHighlight()
        clusterDat <- doTheClustering()
        myPalette <- getColourPalette()
        breaksForPlotly <- (c(myPalette$myBreaks, 1.01, 1.5, 2, 2.5, 3) - myPalette$myBreaks[1])/(3 - myPalette$myBreaks[1])
        colsForPlotly <- paste0(
            "rgb(", col2rgb(myPalette$myCols)[1,], ",",
            col2rgb(myPalette$myCols)[2,], ",",
            col2rgb(myPalette$myCols)[3,], ")"
        )
        p <- plot_ly(z = matData$mat[rev(clusterDat$order), clusterDat$order],
                     x = matData$myLabels[clusterDat$order],
                     y = matData$myLabels[rev(clusterDat$order)],
                     colorscale = list(
                         c(breaksForPlotly[1], colsForPlotly[1]),
                         c(breaksForPlotly[2], colsForPlotly[2]),
                         c(breaksForPlotly[3], colsForPlotly[3]),
                         c(breaksForPlotly[4], colsForPlotly[4]),
                         c(breaksForPlotly[5], colsForPlotly[4]),
                         c(breaksForPlotly[6], "rgb(0,0,255)"),
                         c(breaksForPlotly[7], "rgb(255,255,0)"),
                         c(breaksForPlotly[8], "rgb(0,255,0)"),
                         c(breaksForPlotly[9], "rgb(0,0,0)")
                     ),
                     zmin = myPalette$myBreaks[1],
                     zmax = 3, # color trick for the highlight
                     type = "heatmap",
                     showscale = FALSE
        )
        # style
        layout(p,
               title = "Pairwise correlations",
               autosize = FALSE,
               width = 980,
               height = 980,
               margin = list(l = 250, r = 50, b = 250, t = 50, pad = 4),
               font = list(size = 12),
               xaxis = list(title = ""),
               yaxis = list(title = "")
        )
    })

    output$colourKey2 <- renderPlot(renderColourKey())

    # dendrogram plot --------------------
    myRenderTreePlot <- function() {
        matData <- matAfterHighlight()
        clusterDat <- doTheClustering()

        nSample <- length(matData$myLabels)

        newLabSize <- round((6/nSample^0.6)*10)/10
        if(newLabSize > 1.4) newLabSize <- 1.4
        if(newLabSize < 0.1) newLabSize <- 0.1
        if (input$labelOption == "Automatic") updateSliderInput(session, "labCex", label = NULL, value = newLabSize, min = NULL, max = NULL, step = NULL)

        newMargin <- round((66/nSample^0.25)*10)/10
        if(newMargin > 20) newMargin <- 20
        if(newMargin < 1) newMargin <- 1
        if (input$labelOption == "Automatic") updateSliderInput(session, "margin", label = NULL, value = newMargin, min = NULL, max = NULL, step = NULL)

        oldPar <- list(mar = par()$mar, cex = par()$cex)
        par(mar = c(5, 4, 4, if (input$labelOption == "Automatic") newMargin else input$margin), cex = if (input$labelOption == "Automatic") newLabSize else input$labCex)
        plot(
            clusterDat$dend,
            horiz = TRUE
        )
        par(oldPar)
    }

    output$myTree <- renderPlot(myRenderTreePlot())

    output$downloadTreePng <- downloadHandler("dendrogram.png",
                                            content = function(file) {
                                                png(file, width = 500, height = 950)
                                                myRenderTreePlot()
                                                dev.off()
                                            },
                                            contentType = "image/png")

    output$downloadTreePdf <- downloadHandler("dendrogram.pdf",
                                            content = function(file) {
                                                pdf(file, width = 7.29, height = 13.85)
                                                myRenderTreePlot()
                                                dev.off()
                                            },
                                            contentType = "image/pdf")

    output$downloadTreeSvg <- downloadHandler("dendrogram.svg",
                                              content = function(file) {
                                                  svglite(file, width = 7.29, height = 13.85)
                                                  myRenderTreePlot()
                                                  dev.off()
                                              },
                                              contentType = "image/svg")

    # scatter plot -----------------
    getPairWiseData <- reactive({
        withProgress(value = 1, message = "Extracting data", {
            dataset <- getSelectedDataset()
            if (input$scatterPlotSample1Defined == input$nameOfExpressionFile) {
                expression_file1 <- log10(userExpressionFileAnalysis()$expression$value + 1)
            } else {
                expression_file1 <- dataset$dataMatrix[, which(dataset$annotation$name == input$scatterPlotSample1Defined)]
            }
            if (input$scatterPlotSample2Defined == input$nameOfExpressionFile) {
                expression_file2 <- log10(userExpressionFileAnalysis()$expression$value + 1)
            } else {
                expression_file2 <- dataset$dataMatrix[, which(dataset$annotation$name == input$scatterPlotSample2Defined)]
            }
            if (input$scatterPlotDataScaling == "none") {
                expression_file1 <- 10^expression_file1 - 1
                expression_file2 <- 10^expression_file2 - 1
            } else if (input$scatterPlotDataScaling == "log(e + 1)") {
                expression_file1 <- expression_file1 * log(10)
                expression_file2 <- expression_file2 * log(10)
            } else if (input$scatterPlotDataScaling == "log2(e + 1)") {
                expression_file1 <- expression_file1 * log(10)/log(2)
                expression_file2 <- expression_file2 * log(10)/log(2)
            } else if (input$scatterPlotDataScaling == "asinh(e)") {
                expression_file1 <- asinh(10^expression_file1 - 1)
                expression_file2 <- asinh(10^expression_file2 - 1)
            } else if (input$scatterPlotDataScaling == "1/(1 + e)") {
                expression_file1 <- 1/(10^expression_file1)
                expression_file2 <- 1/(10^expression_file2)
            }
            validate(need(length(expression_file1) == length(expression_file2), "Loading..."))
        })
        return(data.frame(
            geneID = dataset$geneName,
            exp1 = expression_file1,
            exp2 = expression_file2,
            stringsAsFactors = FALSE
        ))
    })

    renderMyScatterPlot <- reactive({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 3, "Loading..."))
        if (input$scatterPlotType == "XY") {
            myScatterPlot <- ggplot(myDF, aes(x = exp1, y = exp2)) + geom_point(size = 1) +
                labs(x = input$scatterPlotSample1Defined, y = input$scatterPlotSample2Defined)
            if (input$scatterPlotGuide) myScatterPlot <- myScatterPlot + geom_abline(intercept = 0, slope = 1, colour = "red")
        } else if (input$scatterPlotType == "MA") {
            myDF$mean <- rowMeans(myDF[, c("exp1", "exp2")])
            myDF$delta <- myDF$exp1 - myDF$exp2
            myScatterPlot <- ggplot(myDF, aes(x = mean, y = delta)) + geom_point(size = 1) +
                labs(x = "mean", y = "difference", title = paste(input$scatterPlotSample1Defined, "\n-", input$scatterPlotSample2Defined))
            if (input$scatterPlotGuide) myScatterPlot <- myScatterPlot + geom_hline(yintercept = 0, colour = "red")
        }
        if (input$scatterPlotRegression) myScatterPlot <- myScatterPlot + geom_smooth(method = "lm")
        myScatterPlot <- myScatterPlot + theme_bw(base_size = 20)
        return(myScatterPlot)
    })

    output$myScatterPlot <- renderPlot({
        withProgress(value = 1, message = "Making scatter plot...",
                     print(renderMyScatterPlot())
        )
    })

    output$scatterPlotMetricsPearson <- renderText({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 3, ""))
        pcc <- cor.test(myDF$exp1, myDF$exp2, method = "pearson")
        paste0(
            "Pearson correlation coefficient: ",
            signif(pcc$estimate, digits = 4),
            " p-value:",
            signif(pcc$p.value, digits = 4)
        )
    })

    output$scatterPlotMetricsSpearman <- renderText({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 3, ""))
        scc <- cor.test(myDF$exp1, myDF$exp2, method = "spearman")
        paste0(
            "Spearman correlation coefficient: ",
            signif(scc$estimate, digits = 4),
            " p-value:",
            signif(scc$p.value, digits = 4)
        )
    })

    output$downloadScatterPlotPng <- downloadHandler("scatterplot.png",
                                            content = function(file) {
                                                png(file, width = 500, height = 500)
                                                print(renderMyScatterPlot())
                                                dev.off()
                                            },
                                            contentType = "image/png")

    output$downloadScatterPlotPdf <- downloadHandler("scatterplot.pdf",
                                            content = function(file) {
                                                pdf(file, width = 7.29, height = 7.29)
                                                print(renderMyScatterPlot())
                                                dev.off()
                                            },
                                            contentType = "image/pdf")

    output$downloadScatterPlotSvg <- downloadHandler("scatterplot.svg",
                                            content = function(file) {
                                                svglite(file, width = 7.29, height = 7.29)
                                                print(renderMyScatterPlot())
                                                dev.off()
                                            },
                                            contentType = "image/svg")

    output$downloadScatterPlotData <- downloadHandler("clusteredMatrix.txt",
                                             content = function(file) {
                                                 myDF <- getPairWiseData()
                                                 colnames(myDF) <- c("geneID", input$scatterPlotSample1Defined, input$scatterPlotSample2Defined)
                                                 write.table(myDF, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
                                             },
                                             contentType = "text/tsv")

    # Metadata table -----------------
    output$tabSampleList <- renderDataTable({
        myTable <- getSelectedDataset()$annotation
        if(!is.null(myTable$url)) {
            myTable$url <- paste0('<a href="', myTable$url, '">link</a>')
        }
        return(myTable)
    }, escape = FALSE)

    output$downloadDatasetTable <- downloadHandler("dataset_table.txt",
                                                         content = function(file) {
                                                             write.table(
                                                                 getSelectedDataset()$annotation,
                                                                 file = file, row.names = FALSE, quote = FALSE, sep = "\t"
                                                             )
                                                         },
                                                         contentType = "text/tsv")

    output$downloadExempleFile <- downloadHandler("rnaseq_mouse_GSE70732_NP1071_FPKM_all_genes_with_header.txt",
                                                   content = function(file) {
                                                       file.copy("www/rnaseq_mouse_GSE70732_NP1071_FPKM_all_genes_with_header.txt", file)
                                                   },
                                                   contentType = "text/tsv")

    # end of server -------------
})