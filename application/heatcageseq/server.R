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
        if(input$fileToUse == "Upload your result file") {
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
            shinyjs::hide("downloadUserCageFile")
            shinyjs::hide("downloadUserCorrelationTable")
        } else {
            shinyjs::show("downloadUserCageFile")
            shinyjs::show("downloadUserCorrelationTable")
        }
    })

    observe({
        shinyjs::hide("widgetForFantom5Human")
        shinyjs::hide("widgetForFantom5Mouse")
        if (input$selectedDataset == "FANTOM5 (human, hg19)") {
            shinyjs::show("widgetForFantom5Human")
        } else if (input$selectedDataset == "FANTOM5 (mouse, mm9)") {
                shinyjs::show("widgetForFantom5Mouse")
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
            load("data/fantom5_human_cage_preload.RData")
            load("data/fantom5_mouse_cage_preload.RData")
            setProgress(value = 1, detail = "loading new dataset")
            if(input$selectedDataset == "FANTOM5 (human, hg19)") {
                load("data/fantom5_human_cage.RData")
                dataset <- fantom5_human_cage
            } else if(input$selectedDataset == "FANTOM5 (mouse, mm9)") {
                load("data/fantom5_mouse_cage.RData")
                dataset <- fantom5_mouse_cage
            }
            setProgress(value = 1, detail = "done!")
        })
        isolate(
            if (input$fileToUse == "Use the example file") {
                if (!input$selectedDataset %in% c("FANTOM5 (mouse, mm9)")) {
                    updateRadioButtons(session = session, "fileToUse", label = NULL, choices = NULL, selected = "Upload your result file")
                }
            }
        )
        return(dataset)
    })

    userCageFileAnalysis_1 <- reactive({

        dataset <- getSelectedDataset()

        if (input$fileToUse == "Use the example file") {
            withProgress(value = 1, message = "Usercage file: ", detail = "reading file", {
                userCageFile <- read_tsv("www/Expression.mm9.HCC.49096peaks_column1.bed", col_names = FALSE)[, 1:6]
                setProgress(value = 1, detail = "intersecting CAGE regions")
                colnames(userCageFile) <- c("chr", "start", "end", "name", "score", "strand")
                userCageFileGR <- with(userCageFile, GRanges(chr, IRanges(start, end), strand = strand, score = score))
                validate(need(
                    input$selectedDataset %in% c("FANTOM5 (mouse, mm9)"),
                    "The example file is a mouse RNA-seq experiment. Please choose a mouse dataset or unload the example."
                ))
                # warnings about missing chromosomes in one of the 2 sets. Let's assume user now what he/her is uploading...
                # Which region overlaps?
                # we ignore peaks presents only in user peak lists.
                suppressWarnings(userOverlap <- findOverlaps(dataset$regionMetaData, userCageFileGR, select = "arbitrary"))
                userCuratedValues <- numeric(nrow(dataset$dataMatrix))
                userCuratedValues[!is.na(userOverlap)] <- as.data.frame(userCageFile)[userOverlap[!is.na(userOverlap)], "score"]
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(log10(userCuratedValues + 1),
                                        dataset$dataMatrix
                ) %>% as.vector
                setProgress(value = 1, detail = "done!")
            })
            return(list(
                "peaks" =  data.frame(as.data.frame(dataset$regionMetaData)[1:3], score = userCuratedValues),
                "correlations" = userCorrelations,
                "linearNormCorrelations" = NULL # to be fill below
            ))
        }

        userCageFileName <- input$cageFile
        if (is.null(userCageFileName)){
            return(list(
                "peaks" = NULL,
                "correlations" = NULL,
                "linearNormCorrelations" = NULL # to be fill below
            ))
        } else {
            withProgress(value = 1, message = "Usercage file: ", detail = "reading file", {
                userCageFile <- read_tsv(userCageFileName$datapath, col_names = input$header)[, 1:6]
                setProgress(value = 1, detail = "intersecting CAGE regions")
                colnames(userCageFile) <- c("chr", "start", "end", "name", "score", "strand")
                userCageFileGR <- with(userCageFile, GRanges(chr, IRanges(start, end), strand = strand, score = score))
                # warnings about missing chromosomes in one of the 2 sets. Let's assume user now what he/her is uploading...
                # Which region overlaps?
                # we ignore peaks presents only in user peak lists.
                suppressWarnings(userOverlap <- findOverlaps(dataset$regionMetaData, userCageFileGR, select = "arbitrary"))
                userCuratedValues <- numeric(nrow(dataset$dataMatrix))
                userCuratedValues[!is.na(userOverlap)] <- as.data.frame(userCageFile)[userOverlap[!is.na(userOverlap)], "score"]
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(log10(userCuratedValues + 1),
                                        dataset$dataMatrix
                                        ) %>% as.vector
                setProgress(value = 1, detail = "done!")
                return(list(
                    "peaks" = data.frame(as.data.frame(dataset$regionMetaData)[1:3], score = userCuratedValues),
                    "correlations" = userCorrelations,
                    "linearNormCorrelations" = NULL # to be fill below
                ))
            })
        }
    })

    userCageFileAnalysis <- reactive({
        userCageFileData <- userCageFileAnalysis_1()
        userCorrelations <- userCageFileData$correlations
        if (!is.null(input$cageFile)  | input$fileToUse == "Use the example file"){
            if (input$maxCorrelation != 0) {
                userCageFileData$linearNormCorrelations <- userCorrelations * input$maxCorrelation / max(userCorrelations)
            }
        }
        return(userCageFileData)
    })

    getExperimentList <- reactive({
        datasetSampleNames <- getSelectedDataset()$annotation$name
        if (!is.null(input$cageFile) | input$fileToUse == "Use the example file") {
            datasetSampleNames <- c(input$nameOfCageFile, datasetSampleNames)
        }
        return(datasetSampleNames)
    })

    output$scatterPlotSample1 <- renderUI({
        selectInput("scatterPlotSample1Defined", "Experiment 1:", getExperimentList(), selected = getExperimentList()[1])
    })

    output$scatterPlotSample2 <- renderUI({
        selectInput("scatterPlotSample2Defined", "Experiment 2:", getExperimentList(), selected = getExperimentList()[2])
    })

    output$tabUserCageFile <- renderDataTable({
        validate(
            need(!is.null(input$cageFile) | input$fileToUse == "Use the example file", "Upload an bed file, or click on a 'heatmap' tab to explore the dataset.")
        )
        userCageFileAnalysis()$peaks
    })

    output$downloadUserCageFile <- downloadHandler("uploaded_cage_file.txt",
                                                 content = function(file) {
                                                     write.table(userCageFileAnalysis()$peaks, file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                                 },
                                                 contentType = "text/tsv")

    getCorrelationTable <- reactive({
        validate(
            need(!is.null(input$cageFile) | input$fileToUse == "Use the example file", "Upload an bed file, or click on a 'heatmap' tab to explore the dataset.")
        )
        if(input$correlationCorrection == "Linear scaling") {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userCageFileAnalysis()$correlations,
                "scaledCorrelation" = userCageFileAnalysis()$linearNormCorrelations,
                stringsAsFactors = FALSE
            )[order(userCageFileAnalysis()$correlations, decreasing = TRUE),])
        } else {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userCageFileAnalysis()$correlations,
                stringsAsFactors = FALSE
            )[order(userCageFileAnalysis()$correlations, decreasing = TRUE),])
        }
    })

    output$tabUserCorrelationTable <- renderDataTable(getCorrelationTable())

    output$downloadUserCorrelationTable <- downloadHandler("heatCAGEseq_correlations.txt",
                                                           content = function(file) {
                                                              write.table(getCorrelationTable(), file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                                           },
                                                           contentType = "text/tsv")

    subsetMatrix <- reactive({
        dataset <- getSelectedDataset()
        workingMatrix <- dataset$correlationMatrix
        keep<- 1:nrow(workingMatrix)
        # filtering is dataset-depedent...

        if (input$selectedDataset == "FANTOM5 (human, hg19)") {
            if (is.null(input$f5h_cells)) {
                temp_f5h_cells <- unique(fantom5_human_cage$annotation$tissue)
            } else {
                temp_f5h_cells <- input$f5h_cells
            }
            if (input$f5h_isCellLine == "all") {
                temp_f5h_isCellLine <- rep(TRUE, nrow(fantom5_human_cage$annotation))
            } else if (input$f5h_isCellLine == "only cell line") {
                temp_f5h_isCellLine <- fantom5_human_cage$annotation$isCellType
            } else if (input$f5h_isCellLine == "only non cell line") {
                temp_f5h_isCellLine <- !fantom5_human_cage$annotation$isCellType
            }
            keep <- which(
                dataset$annotation$tissue %in% temp_f5h_cells &
                temp_f5h_isCellLine
            )
        } else if (input$selectedDataset == "FANTOM5 (mouse, mm9)") {
            if (is.null(input$f5m_cells)) {
                temp_f5m_cells <- unique(fantom5_mouse_cage$annotation$tissue)
            } else {
                temp_f5m_cells <- input$f5m_cells
            }
            keep <- which(
                dataset$annotation$tissue %in% temp_f5m_cells
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
        if(!is.null(userCageFileAnalysis()$correlations)) {
            userCorrelations <- userCageFileAnalysis()$correlations
            if (input$correlationCorrection == "Linear scaling") {
                userCorrelations <- userCageFileAnalysis()$linearNormCorrelations
            }
            workingMatrix <- rbind(workingMatrix, userCorrelations[keep])
            workingMatrix <- cbind(workingMatrix, c(userCorrelations[keep], 1))
            myLabels <- c(myLabels, input$nameOfCageFile)
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
        if (input$highlight & (!is.null(input$cageFile) | input$fileToUse == "Use the example file")) {
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

        heatmap(matData$mat,
                Rowv = clusterDat$dend,
                Colv = "Rowv",
                scale = "none",
                labRow = matData$myLabels,
                labCol = matData$myLabels,
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

    output$colourKey2 <- renderPlot(renderColourKey())

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

    # scatter plot -----------------
    getPairWiseData <- reactive({
        dataset <- getSelectedDataset()
        if (input$scatterPlotSample1Defined == input$nameOfCageFile) {
            expression_file1 <- log10(userCageFileAnalysis()$peaks$score + 1)
        } else {
            expression_file1 <- dataset$dataMatrix[, which(dataset$annotation$name == input$scatterPlotSample1Defined)]
        }
        if (input$scatterPlotSample2Defined == input$nameOfCageFile) {
            expression_file2 <- log10(userCageFileAnalysis()$peaks$score + 1)
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
        return(data.frame(
            as.data.frame(dataset$regionMetaData)[,1:3],
            exp1 = expression_file1,
            exp2 = expression_file2,
            stringsAsFactors = FALSE
        ))
    })

    renderMyScatterPlot <- reactive({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 5, "Loading..."))
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
        myScatterPlot <- myScatterPlot + theme_bw(base_size = 18)
        return(myScatterPlot)
    })

    output$myScatterPlot <- renderPlot({
        withProgress(value = 1, message = "Making scatter plot...",
            print(renderMyScatterPlot())
        )
    })

    output$scatterPlotMetricsPearson <- renderText({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 5, ""))
        pcc <- cor(myDF$exp1, myDF$exp2, method = "pearson")
        paste0(
            "Pearson correlation coefficient: ",
            round(pcc, digits = 4)
        )
    })

    output$scatterPlotMetricsSpearman <- renderText({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 5, ""))
        scc <- cor(myDF$exp1, myDF$exp2, method = "spearman")
        paste0(
            "Spearman correlation coefficient: ",
            round(scc, digits = 4)
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
                                                          colnames(myDF) <- c("chr", "start", "end", input$scatterPlotSample1Defined, input$scatterPlotSample2Defined)
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

    output$downloadExempleFile <- downloadHandler("Expression.mm9.HCC.49096peaks_column1.bed",
                                                  content = function(file) {
                                                      file.copy("www/Expression.mm9.HCC.49096peaks_column1.bed", file)
                                                  },
                                                  contentType = "text/tsv")

})