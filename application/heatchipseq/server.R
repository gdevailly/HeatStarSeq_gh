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
                    "Something is wrong, sorry. :(
                    \nThe example file is a human ChIP-seq experiment. You should choose a human dataset, I am sure it will work better."
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
                "linearNormCorrelations" = NULL # to be fill below
            ))
        }

        userPeakFileName <- input$peakFile
        if (is.null(userPeakFileName)){
            userPeakFile <- NULL
            userCorrelations <- NULL
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
            "linearNormCorrelations" = NULL # to be fill below
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
            need(length(keep) >= 3, "Fewer than 3 experiments match your criteria. Please select more experiments.")
        )
        validate(
            need(length(keep) <= 1100, "More than 1100 experiments match your criteria. Please select less experiments.")
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
            if (input$distOption == "1 - Pearson correlation coefficient") {
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

        heatmap(myMat,
                Rowv = Rowv,
                Colv = Colv,
                scale = "none",
                labRow = labRow,
                labCol = labCol,
                margins = if(input$labelOption == "Automatic") rep(newMargin, 2) else  rep(input$margin, 2),
                cexRow = if(input$labelOption == "Automatic") newLabSize else input$labCex,
                cexCol = if(input$labelOption == "Automatic") newLabSize else input$labCex,
                breaks = seq(-0.5, 3, length.out = 256), # color trick for the highlight
                col = colorRampPalette(c("blue", "white", "red", "black", "blue", "yellow", "green", "black"))(255),
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

    output$myHeatmap <- renderPlot(myRenderPlot())

    output$myPlotlyHeatmap <- renderPlotly({
        matData <- matAfterHighlight()
        clusterDat <- doTheClustering()
        p <- plot_ly(z = matData$mat[rev(clusterDat$order), clusterDat$order],
                     x = matData$myLabels[clusterDat$order],
                     y = matData$myLabels[rev(clusterDat$order)],
                     colorscale = list(
                                       c(0.0, "rgb(0,0,255)"),
                                       c(0.1428571, "rgb(255,255,255)"),
                                       c(0.2857143, "rgb(255,0,0)"),
                                       c(0.4285714, "rgb(0,0,0)"),
                                       c(0.5714286, "rgb(0,0,255)"),
                                       c(0.7142857, "rgb(255,255,0)"),
                                       c(0.8571429, "rgb(0,255,0)"),
                                       c(1, "rgb(0,0,0)")
                                       ),
                     zmin = -0.5,
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