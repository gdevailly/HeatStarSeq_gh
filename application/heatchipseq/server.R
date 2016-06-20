library(readr)
library(svglite)
library(cba)

options(shiny.maxRequestSize = 10*1024^2) # max file size, 10Mb

shinyServer(function(input, output, session) {

    # UI elements activations
    observe({
        junkVar <- input$fileFormatInstructions
        shinyjs::toggle("div_fileFormatInstructions")
    })

    observe({
        if(input$fileToUse == "Upload your peak file") {
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
        if (is.null(input$peakFile) & input$fileToUse != "Use the example file") {
            shinyjs::hide("downloadUserPeaks")
            shinyjs::hide("downloadUserCorrelationTable")
        } else {
            shinyjs::show("downloadUserPeaks")
            shinyjs::show("downloadUserCorrelationTable")
        }
    })

    observe({
        shinyjs::hide("widgetForEncodeHuman")
        shinyjs::hide("widgetForCodexHuman")
        shinyjs::hide("widgetForCodexMouse")
        shinyjs::hide("widgetForModEncodeD")
        if (input$selectedDataset == "ENCODE TFBS ChIP-seq (human, hg19)") {
            shinyjs::show("widgetForEncodeHuman")
        } else if (input$selectedDataset == "CODEX ChIP-seq (human, hg19)") {
            shinyjs::show("widgetForCodexHuman")
        } else if (input$selectedDataset == "CODEX ChIP-seq (mouse, mm10)") {
            shinyjs::show("widgetForCodexMouse")
        } else if (input$selectedDataset == "modEncode TF ChIP-seq (drosophila, r5)") {
            shinyjs::show("widgetForModEncodeD")
        } else if (input$selectedDataset == "other (soon)") {
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

    # output computations
    getSelectedDataset <- reactive({
        withProgress(value = 1, message = "Loading dataset: ", detail = "removing old dataset", {
            # this remove from memmory non-selected datasets
            load("data/encode_preload.RData")
            load("data/codex_preload.RData")
            load("data/codex_human_chip_preload.RData")
            load("data/modEncodeD_ChIPseq_preload.RData")
            setProgress(value = 1, detail = "loading new dataset")
            if(input$selectedDataset == "ENCODE TFBS ChIP-seq (human, hg19)") {
                load("data/encode.RData")
                dataset <- encode
            } else if(input$selectedDataset == "CODEX ChIP-seq (mouse, mm10)") {
                load("data/codex.RData")
                dataset <- codex
            } else if(input$selectedDataset == "CODEX ChIP-seq (human, hg19)") {
                load("data/codex_human_chip.RData")
                dataset <- codex_human_chip
            } else if(input$selectedDataset == "modEncode TF ChIP-seq (drosophila, r5)") {
                load("data/modEncodeD_ChIPseq.RData")
                dataset <- modEncodeD_ChIPseq
            }
            setProgress(value = 1, detail = "done!")
        })
        isolate(
            if (input$fileToUse == "Use the example file") {
                if (!input$selectedDataset %in% c("ENCODE TFBS ChIP-seq (human, hg19)", "CODEX ChIP-seq (human, hg19)")) {
                    updateRadioButtons(session = session, "fileToUse", label = NULL, choices = NULL, selected = "Upload your peak file")
                }
            }
        )
        return(dataset)
    })

    userPeakFileAnalysis_1 <- reactive({

        if (input$fileToUse == "Use the example file") {
            withProgress(value = 1, message = "User peak file: ", detail = "reading file", {
                userPeakFile <- read_tsv("www/chipseq_human_hg19_GSM1890761_ERa_peaks_no_header.bed", col_names = FALSE)[, 1:3]
                colnames(userPeakFile) <- c("chr","start","end")
                setProgress(value = 1, detail = "intersecting peaks")
                userPeakFileGR <- with(userPeakFile, GRanges(chr, IRanges(start, end)))
                dataset <- getSelectedDataset()
                validate(need(
                    input$selectedDataset %in% c("ENCODE TFBS ChIP-seq (human, hg19)", "CODEX ChIP-seq (human, hg19)"),
                    "The example file is a  human ChIP-seq experiment. Please choose a human dataset or unload the example."
                ))
                # warnings about missing chromosomes in one of the 2 sets. Let's assume user now what it is uploading...
                # Which encode region overlaps?
                # we ignore peaks presents only in user peak lists for some reasons.
                suppressWarnings(userOverlap <- overlapsAny(dataset$regionMetaData, userPeakFileGR))
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(userOverlap,
                                        dataset$dataMatrix
                ) %>% as.vector
                setProgress(value = 1, detail = "done!")
            })
            return(list(
                "peaks" = userPeakFile,
                "correlations" = userCorrelations,
                "linearNormCorrelations" = NULL, # to be fill below
                "overlaps" = userOverlap
            ))
        }

        userPeakFileName <- input$peakFile
        if (is.null(userPeakFileName)){
            userPeakFile <- NULL
            userCorrelations <- NULL
            userOverlap <- NULL
        } else {
            withProgress(value = 1, message = "User peak file: ", detail = "reading file", {
                userPeakFile <- read_tsv(userPeakFileName$datapath, col_names = input$header)[, 1:3]
                colnames(userPeakFile) <- c("chr","start","end")
                setProgress(value = 1, detail = "intersecting peaks")
                userPeakFileGR <- with(userPeakFile, GRanges(chr, IRanges(start, end)))
                dataset <- getSelectedDataset()
                # warnings about missing chromosomes in one of the 2 sets. Let's assume user now what it is uploading...
                # Which encode region overlaps?
                # we ignore peaks presents only in user peak lists for some reasons.
                suppressWarnings(userOverlap <- overlapsAny(dataset$regionMetaData, userPeakFileGR))
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(userOverlap,
                                        dataset$dataMatrix
                                        ) %>% as.vector
                setProgress(value = 1, detail = "done!")
            })
        }
        return(list(
            "peaks" = userPeakFile,
            "correlations" = userCorrelations,
            "linearNormCorrelations" = NULL, # to be fill below
            "overlaps" = userOverlap
        ))
    })

    userPeakFileAnalysis <- reactive({
        userPeakFileData <- userPeakFileAnalysis_1()
        userCorrelations <- userPeakFileData$correlations
        if (!is.null(input$peakFile) | input$fileToUse == "Use the example file"){
            if (input$maxCorrelation != 0) {
                userPeakFileData$linearNormCorrelations <- userCorrelations * input$maxCorrelation / max(userCorrelations)
            }
        }
        return(userPeakFileData)
    })

    getExperimentList <- reactive({
        datasetSampleNames <- sort(getSelectedDataset()$annotation$name)
        if (!is.null(input$peakFile) | input$fileToUse == "Use the example file") {
            datasetSampleNames <- c(input$nameOfPeakFile, datasetSampleNames)
        }
        return(datasetSampleNames)
    })

    output$barPlotSample1 <- renderUI({
        selectInput("barPlotSample1Defined", "Experiment 1:", getExperimentList(), selected = getExperimentList()[1])
    })

    output$barPlotSample2 <- renderUI({
        selectInput("barPlotSample2Defined", "Experiment 2:", getExperimentList(), selected = getExperimentList()[2])

    })

    output$tabUserPeaks <- renderDataTable({
        validate(
            need(!is.null(input$peakFile) | input$fileToUse == "Use the example file", "Upload an peak file, or click on a 'heatmap' tab to explore the dataset.")
        )
        userPeakFileAnalysis()$peaks
    })

    output$downloadUserPeaks <- downloadHandler("uploaded_peak_list.txt",
                                                 content = function(file) {
                                                     write.table(userPeakFileAnalysis()$peaks, file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                                 },
                                                 contentType = "text/tsv")

    getCorrelationTable <- reactive({
        validate(
            need(!is.null(input$peakFile) | input$fileToUse == "Use the example file", "Upload an peak file, or click on a 'heatmap' tab to explore the dataset.")
        )
        if(input$correlationCorrection == "Linear scaling") {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userPeakFileAnalysis()$correlations,
                "scaledCorrelation" = userPeakFileAnalysis()$linearNormCorrelations,
                stringsAsFactors = FALSE
            )[order(userPeakFileAnalysis()$correlations, decreasing = TRUE), ])
        } else {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userPeakFileAnalysis()$correlations,
                stringsAsFactors = FALSE
            )[order(userPeakFileAnalysis()$correlations, decreasing = TRUE), ])
        }
    })

    output$tabUserCorrelationTable <- renderDataTable(getCorrelationTable())

    output$downloadUserCorrelationTable <- downloadHandler("heatChIPseq_correlations.txt",
                                                           content = function(file) {
                                                              write.table(getCorrelationTable(), file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                                           },
                                                           contentType = "text/tsv")

    subsetMatrix <- reactive({
        dataset <- getSelectedDataset()
        workingMatrix <- dataset$correlationMatrix
        keep<- 1:nrow(workingMatrix)
        # filtering is dataset-depedent...
        if (input$selectedDataset == "ENCODE TFBS ChIP-seq (human, hg19)") {
            if (is.null(input$cells)) {
                temp_cells <- unique(encode$annotation$cellLine)
            } else {
                temp_cells <- input$cells
            }
            if (is.null(input$TF)) {
                temp_TF <- unique(encode$annotation$tf)
            } else {
                temp_TF <- input$TF
            }
            keep <- which(
                dataset$annotation$cellLine %in% temp_cells &
                dataset$annotation$tf %in% temp_TF
            )

        } else if (input$selectedDataset == "CODEX ChIP-seq (human, hg19)") {
            if (is.null(input$TF_ch)) {
                temp_TF_ch <- unique(codex_human_chip$annotation$tf)
            } else {
                temp_TF_ch <- input$TF_ch
            }
            if (input$filterCellsBy_ch == "Cell type") {
                if (is.null(input$cell_types_ch)) {
                    temp_cell_types_ch <- unique(codex_human_chip$annotation$cellType)
                } else {
                    temp_cell_types_ch <- input$cell_types_ch
                }
                keep <- which(
                    dataset$annotation$cellType %in% temp_cell_types_ch &
                    dataset$annotation$tf %in% temp_TF_ch
                )
            } else if (input$filterCellsBy_ch == "Cell subtype") {
                if(is.null(input$cell_subtypes_ch)) {
                    temp_cell_subtypes_ch <- codex_human_chip$annotation$cellSubtype
                } else {
                    temp_cell_subtypes_ch <- input$cell_subtypes_ch
                }
                keep <- which(
                    dataset$annotation$cellSubtype %in% temp_cell_subtypes_ch &
                    dataset$annotation$tf %in% temp_TF_ch
                )
            }

        } else if (input$selectedDataset == "CODEX ChIP-seq (mouse, mm10)") {
            if (is.null(input$TF_m)) {
                temp_TF_m <- unique(codex$annotation$tf)
            } else {
                temp_TF_m <- input$TF_m
            }
            if (input$filterCellsBy == "Cell type") {
                if (is.null(input$cell_types_m)) {
                    temp_cell_types_m <- unique(codex$annotation$cellType)
                } else {
                    temp_cell_types_m <- input$cell_types_m
                }
                keep <- which(
                    dataset$annotation$cellType %in% temp_cell_types_m &
                    dataset$annotation$tf %in% temp_TF_m
                )
            } else if (input$filterCellsBy == "Cell subtype") {
                if(is.null(input$cell_subtypes_m)) {
                    temp_cell_subtypes_m <- codex$annotation$cellSubtype
                } else {
                    temp_cell_subtypes_m <- input$cell_subtypes_m
                }
                keep <- which(
                    dataset$annotation$cellSubtype %in% temp_cell_subtypes_m &
                    dataset$annotation$tf %in% temp_TF_m
                )
            }
        } else if (input$selectedDataset == "modEncode TF ChIP-seq (drosophila, r5)") {
            if (is.null(input$tf_med)) {
                temp_tf_med <- unique(modEncodeD_ChIPseq$annotation$antibody)
            } else {
                temp_tf_med <- input$tf_med
            }
            if (is.null(input$stage_med)) {
                temp_stage_med <- unique(modEncodeD_ChIPseq$annotation$devStage)
            } else {
                temp_stage_med <- input$stage_med
            }
            if (is.null(input$strain_med)) {
                temp_strain_med <- unique(modEncodeD_ChIPseq$annotation$strain)
            } else {
                temp_strain_med <- input$strain_med
            }
            keep <- which(
                dataset$annotation$antibody %in% temp_tf_med &
                dataset$annotation$devStage %in% temp_stage_med &
                dataset$annotation$strain %in% temp_strain_med
            )
        }

        myLabels <- dataset$annotation$name # annotation must have a name column, with _unique_ elements
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
        if(!is.null(userPeakFileAnalysis()$correlations)) {
            userCorrelations <- userPeakFileAnalysis()$correlations
            if (input$correlationCorrection == "Linear scaling") {
                userCorrelations <- userPeakFileAnalysis()$linearNormCorrelations
            }
            workingMatrix <- rbind(workingMatrix, userCorrelations[keep])
            workingMatrix <- cbind(workingMatrix, c(userCorrelations[keep], 1))
            myLabels <- c(myLabels, input$nameOfPeakFile)
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
        if (input$highlight & (!is.null(input$peakFile) | input$fileToUse == "Use the example file")) {
            # we increase user cor value by 2 to do the color trick
            matAA$mat[, ncol(matAA$mat)] <- matAA$mat[, ncol(matAA$mat)] + 2
            matAA$mat[nrow(matAA$mat), 1:(ncol(matAA$mat) - 1)] <- matAA$mat[nrow(matAA$mat), 1:(ncol(matAA$mat) - 1)] + 2
        }
        return(matAA)
    })

    # heatmap
    getColourPalette <- reactive({
        junkVar <- input$applyColoursOptions
        isolate({
            myCols <- c(input$col_col1, input$col_col2, input$col_col3, input$col_col4)
            myBreaks <- c(input$col_val1, input$col_val2, input$col_val3, input$col_val4)
        })
        return(list("myCols" = myCols, "myBreaks" = myBreaks))
    })

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

    output$myHeatmap <- renderPlot(
        withProgress(value = 1, message = "Ploting...", {
            myRenderPlot()
            setProgress(value = 1, detail = "done!")
        })
    )

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

    # plotly heatmap
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

    # dendrogram
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

    output$myTree <- renderPlot(myRenderTreePlot())

    # pairwise barPlot
    getPairWiseData <- reactive({
        dataset <- getSelectedDataset()
        if (input$barPlotSample1Defined == input$nameOfPeakFile) {
            peaks1 <- userPeakFileAnalysis()$overlaps
        } else {
            peaks1 <- dataset$dataMatrix[, which(dataset$annotation$name == input$barPlotSample1Defined)]
        }
        if (input$barPlotSample2Defined == input$nameOfPeakFile) {
            peaks2 <- userPeakFileAnalysis()$overlaps
        } else {
            peaks2 <- dataset$dataMatrix[, which(dataset$annotation$name == input$barPlotSample2Defined)]
        }
        validate(need(length(peaks1) == length(peaks2), "Loading..."))
        return(data.frame(
            as.data.frame(dataset$regionMetaData)[, 1:3],
            peaks1 = peaks1,
            peaks2 = peaks2,
            stringsAsFactors = FALSE
        ))
    })

    renderMyBarPlot <- reactive({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 5, "Loading..."))
        name1 <- input$barPlotSample1Defined %>%
            gsub("+", "\n", ., fixed = TRUE) %>%
            gsub(" ", "\n", ., fixed = TRUE) %>%
            gsub("[\n]+", "\n", ., fixed = FALSE)
        name2 <- input$barPlotSample2Defined %>%
            gsub("+", "\n", ., fixed = TRUE) %>%
            gsub(" ", "\n", ., fixed = TRUE) %>%
            gsub("[\n]+", "\n", ., fixed = FALSE)
        dfForBarplot <- data.frame(
            experiment = factor(
                c(rep(name1, 2), rep(name2, 2)),
                levels = c(name1, name2)
            ),
            status = rep(c("common", "unique"), 2),
            nPeaks = c(
                length(which(myDF$peaks1 &  myDF$peaks2)),
                length(which(myDF$peaks1 & !myDF$peaks2)),
                length(which(myDF$peaks2 &  myDF$peaks1)),
                length(which(myDF$peaks2 & !myDF$peaks1))
            ),
            stringsAsFactors = FALSE
        )
        myBarPlot <- ggplot(dfForBarplot, aes(x = experiment, y = nPeaks, fill = status)) + geom_bar(stat = "identity") +
            labs(y = "number of peaks")
        if (input$barPlotGuide) myBarPlot <- myBarPlot + geom_hline(yintercept = dfForBarplot$nPeaks[1])
        myBarPlot <- myBarPlot + theme_bw(base_size = 20) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "top")
        return(myBarPlot)
    })

    output$myBarPlot <- renderPlot({
        withProgress(value = 1, message = "Making bar plot...",
                     print(renderMyBarPlot())
        )
    })

    dataForBarPlotTests <- reactive({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 5, "Loading..."))
        myTable <- data.frame(
            category = c(
                "common peaks",
                paste(input$barPlotSample1Defined, "only"),
                paste(input$barPlotSample2Defined, "only"),
                "peak space"
            ),
            nPeaks = c(
                length(which(myDF$peaks1 &  myDF$peaks2)),
                length(which(myDF$peaks1 & !myDF$peaks2)),
                length(which(myDF$peaks2 & !myDF$peaks1)),
                nrow(myDF)
            )
        )
        return(myTable)
    })

    output$tabBarPlot <- renderTable(dataForBarPlotTests(), include.rownames = FALSE)

    output$barPlotCorrelation <- renderText({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 5, "Loading..."))
        cc <- cor.test(as.numeric(myDF$peaks1), as.numeric(myDF$peaks2))
        return(paste("Correlation coefficient:", round(cc$estimate, digits = 4), "  p-value:", signif(cc$p.value, digits = 4)))
    })

    output$barPlotJaccard <- renderText({
        myDF <- dataForBarPlotTests()
        myJaccard <- myDF$nPeaks[1] / (myDF$nPeaks[1] + myDF$nPeaks[2] + myDF$nPeaks[3])
        return(paste("Jaccard index:", round(myJaccard, digits = 4)))
    })

    output$downloadBarPlotPng <- downloadHandler("barplot.png",
                                                     content = function(file) {
                                                         png(file, width = 500, height = 500)
                                                         print(renderMyBarPlot())
                                                         dev.off()
                                                     },
                                                     contentType = "image/png")

    output$downloadBarPlotPdf <- downloadHandler("barplot.pdf",
                                                     content = function(file) {
                                                         pdf(file, width = 7.29, height = 7.29)
                                                         print(renderMyBarPlot())
                                                         dev.off()
                                                     },
                                                     contentType = "image/pdf")

    output$downloadBarPlotSvg <- downloadHandler("barplot.svg",
                                                     content = function(file) {
                                                         svglite(file, width = 7.29, height = 7.29)
                                                         print(renderMyBarPlot())
                                                         dev.off()
                                                     },
                                                     contentType = "image/svg")

    output$downloadBarPlotData <- downloadHandler("barplot.txt",
                                                      content = function(file) {
                                                          myDF <- getPairWiseData()
                                                          colnames(myDF) <- c("chr", "start", "end", input$barPlotSample1Defined, input$barPlotSample2Defined)
                                                          write.table(myDF, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
                                                      },
                                                      contentType = "text/tsv")

    # metadata table
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

    output$downloadExempleFile <- downloadHandler("chipseq_human_hg19_GSM1890761_ERa_peaks_no_header.bed",
                                                  content = function(file) {
                                                      file.copy("www/chipseq_human_hg19_GSM1890761_ERa_peaks_no_header.bed", file)
                                                  },
                                                  contentType = "text/tsv")

})
