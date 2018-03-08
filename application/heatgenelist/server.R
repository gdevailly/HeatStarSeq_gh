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
        if(input$fileToUse == "Upload your gene list") {
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
        shinyjs::hide("widgetForHumanKo")
        shinyjs::hide("widgetForMouseKo")

        if (input$dataset == "Knock-out (human)") {
            shinyjs::show("widgetForHumanKo")
        } else if (input$dataset == "Knock-out (mouse)") {
            shinyjs::show("widgetForMouseKo")
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
            load("data/human_ko_preload.Rdata")
            load("data/mouse_ko_preload.Rdata")
            setProgress(value = 1, detail = "loading new dataset")
            if (input$dataset == "Knock-out (human)") {
                load("data/human_ko.Rdata")
                dataset <- human_ko
            } else if (input$dataset == "Knock-out (mouse)") {
                load("data/mouse_ko.Rdata")
                dataset <- mouse_ko
            } else if (input$dataset == "other (soon)") {
                # fill this
            }
            setProgress(value = 1, detail = "done!")
        })
        isolate(
            if (input$fileToUse == "Use the example file") {
                if (!input$dataset %in% c("Knock-out (human)")) {
                    updateRadioButtons(session = session, "fileToUse", label = NULL, choices = NULL, selected = "Upload your gene list")
                }
            }
        )
        return(dataset)
    })

    userGeneListAnalysis_1 <- reactive({

        if (input$fileToUse == "Use the example file") {
            withProgress(value = 1, message = "Loading example: ", detail = "reading file", {
                userGeneList_temp <- read_tsv("www/ESR1_breastcancercells_PMID15023342.txt", col_names = FALSE)[, 1, drop = TRUE]
                setProgress(value = 1, detail = "intersecting gene names")
                userGeneList <- rep(TRUE, time = length(userGeneList_temp))
                names(userGeneList) <- tolower(userGeneList_temp)
                dataset <- getSelectedDataset()
                validate(need(
                    input$dataset %in% c("Knock-out (human)"),
                    "The example file is human gene list. Please choose a human dataset or unload the example."
                ))
                userGeneList <- userGeneList[tolower(dataset$geneName)]
                userGeneList[is.na(userGeneList)] <- FALSE
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(userGeneList, dataset$dataMatrix) %>% as.vector
                userGeneList <- data.frame(
                    geneName = dataset$geneName,
                    isPresent = userGeneList,
                    stringsAsFactors = FALSE
                )
                setProgress(value = 1, detail = "done!")
            })
            return(list(
                "expression" = userGeneList,
                "correlations" = userCorrelations,
                "linearNormCorrelations" = NULL # to be fill below
            ))
        }

        userExpressionFileName <- input$expressionFile
        if (is.null(userExpressionFileName)){
            userGeneList <- NULL
            userCorrelations <- NULL
        } else {
            withProgress(value = 1, message = "User gene list: ", detail = "reading file", {
                userGeneList_temp <- read_tsv(userExpressionFileName$datapath, col_names = input$header)[, 1, drop = TRUE]
                setProgress(value = 1, detail = "intersecting gene names")
                userGeneList <- rep(TRUE, time = length(userGeneList_temp))
                names(userGeneList) <- tolower(userGeneList_temp)
                dataset <- getSelectedDataset()
                userGeneList <- userGeneList[tolower(dataset$geneName)]
                userGeneList[is.na(userGeneList)] <- FALSE
                validate(need(
                    length(userGeneList) > 1,
                    "Something went wrong, sorry. :(
                    \nUnable to compute correlations.
                    \nThis is likely because:
                    \n-the wrong dataset is selected (ie you may have uploaded a mouse gene list, but the selected dataset is for human).
                    \n-the file formatting is not recognized. It should be a single column text file containing gene symbol (i.e. MBD2 (human) or Mbd2 (Mouse)). Case insensitive.
                    \nPlease contact us if you think you have found a bug."
                ))
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(userGeneList, dataset$dataMatrix) %>% as.vector
                userGeneList <- data.frame(
                    geneName = dataset$geneName,
                    isPresent = userGeneList,
                    stringsAsFactors = FALSE
                )
                setProgress(value = 1, detail = "done!")
            })
        }
        return(list(
            "expression" = userGeneList,
            "correlations" = userCorrelations
        ))
    })

    userGeneListAnalysis <- reactive({
        userGeneListData <- userGeneListAnalysis_1()
        userCorrelations <- userGeneListData$correlations
        if (!is.null(input$expressionFile) | input$fileToUse == "Use the example file") {
            if (input$maxCorrelation != 0 & input$correlationCorrection != "None") {
                userGeneListData$linearNormCorrelations <- userCorrelations * input$maxCorrelation / max(userCorrelations)
            }
        }
        return(userGeneListData)
    })

    getExperimentList <- reactive({
        datasetSampleNames <- sort(getSelectedDataset()$annotation$name)
        if (!is.null(input$expressionFile) | input$fileToUse == "Use the example file") {
            datasetSampleNames <- c(input$nameOfExpressionFile, datasetSampleNames)
        }
        return(datasetSampleNames)
    })

    output$barPlotSample1 <- renderUI({
        withProgress(value = 1, message = "Creating list of experiments 1/2", {
            selectInput("barPlotSample1Defined", "Experiment 1:", getExperimentList(), selected = getExperimentList()[1])
        })
    })

    output$barPlotSample2 <- renderUI({
        withProgress(value = 1, message = "Creating list of experiments 2/2", {
            selectInput("barPlotSample2Defined", "Experiment 2:", getExperimentList(), selected = getExperimentList()[2])
        })
    })

    output$tabUserExpressionFile <- renderDataTable({
        validate(
            need(!is.null(input$expressionFile) | input$fileToUse == "Use the example file", "Upload an expression file, or click on a 'heatmap' tab to explore the dataset.")
        )
        userGeneListAnalysis()$expression
    })

    output$downloadUserGeneList <- downloadHandler("uploaded_expression_file.txt",
                                                         content = function(file) {
                                                             write.table(
                                                                 userGeneListAnalysis()$expression,
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
                "correlation" = userGeneListAnalysis()$correlations,
                "scaledCorrelation" = userGeneListAnalysis()$linearNormCorrelations,
                stringsAsFactors = FALSE
            )[order(userGeneListAnalysis()$correlations, decreasing = TRUE), ])
        } else {
            return(data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userGeneListAnalysis()$correlations,
                stringsAsFactors = FALSE
            )[order(userGeneListAnalysis()$correlations, decreasing = TRUE), ])
        }
    })

    output$tabUserCorrelationTable <- renderDataTable(getCorrelationTable())

    output$downloadUserCorrelationTable <- downloadHandler("heatGeneList_correlations.txt",
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
        if (input$dataset == "Knock-out (human)") {
            if (is.null(input$cell_human_ko)) {
                temp_cell_human_ko <- unique(human_ko$annotation$cell_type)
            } else {
                temp_cell_human_ko <- input$cell_human_ko
            }
            if (is.null(input$tf_human_ko)) {
                temp_tf_human_ko <- unique(human_ko$annotation$ko_gene)
            } else {
                temp_tf_human_ko <- input$tf_human_ko
            }
            keep <- which(
                dataset$annotation$cell_type %in% temp_cell_human_ko &
                    dataset$annotation$ko_gene %in% temp_tf_human_ko
            )
        } else if (input$dataset == "Knock-out (mouse)") {
            if (is.null(input$cell_mouse_ko)) {
                temp_cell_mouse_ko <- unique(mouse_ko$annotation$cell_type)
            } else {
                temp_cell_mouse_ko <- input$cell_mouse_ko
            }
            if (is.null(input$tf_mouse_ko)) {
                temp_tf_mouse_ko <- unique(mouse_ko$annotation$ko_gene)
            } else {
                temp_tf_mouse_ko <- input$tf_mouse_ko
            }
            keep <- which(
                dataset$annotation$cell_type %in% temp_cell_mouse_ko &
                    dataset$annotation$ko_gene %in% temp_tf_mouse_ko
            )
        } else { # modify accordingly for new datasets!

        }

        myLabels <- dataset$annotation$name
        validate(
            need(length(keep) >= 3, "Less than 3 experiments match your criteria. Please select more experiments.")
        )
        validate(
            need(length(keep) <= 2000, "More than 2000 experiments match your criteria. Please select fewer experiments.")
        )
        if(length(keep) >= 3) {
            workingMatrix <- workingMatrix[keep, keep]
            myLabels <- myLabels[keep]
        }

        # merging user data in correlation matrix
        if(!is.null(userGeneListAnalysis()$correlations)) {
            userCorrelations <- userGeneListAnalysis()$correlations
            if (input$correlationCorrection == "Linear scaling") {
                userCorrelations <- userGeneListAnalysis()$linearNormCorrelations
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
        if (newLabSize > 3) newLabSize <- 3
        if (newLabSize < 0.1) newLabSize <- 0.1
        if (input$labelOption == "Automatic") updateSliderInput(session, "labCex", label = NULL, value = newLabSize, min = NULL, max = NULL, step = NULL)

        newMargin <- round((66/nSample^0.25)*10)/10
        if (newMargin > 50) newMargin <- 50
        if (newMargin < 1) newMargin <- 1
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

    # bar plot -----------------
    getPairWiseData <- reactive({
        withProgress(value = 1, message = "Extracting data", {
            dataset <- getSelectedDataset()
            if (input$barPlotSample1Defined == input$nameOfExpressionFile) {
                genelist_file1 <- userGeneListAnalysis()$expression$isPresent
            } else {
                genelist_file1 <- dataset$dataMatrix[, which(dataset$annotation$name == input$barPlotSample1Defined)]
            }
            if (input$barPlotSample2Defined == input$nameOfExpressionFile) {
                genelist_file2 <- userGeneListAnalysis()$expression$isPresent
            } else {
                genelist_file2 <- dataset$dataMatrix[, which(dataset$annotation$name == input$barPlotSample2Defined)]
            }
            validate(need(length(genelist_file1) == length(genelist_file2), "Loading..."))
        })
        return(data.frame(
            geneID = dataset$geneName,
            exp1 = genelist_file1,
            exp2 = genelist_file2,
            stringsAsFactors = FALSE
        ))
    })

    renderMyBarPlot <- reactive({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 3, "Loading..."))
        name1 <- input$barPlotSample1Defined
        name2 <- input$barPlotSample2Defined
        dfForBarplot <- data.frame(
            experiment = factor(
                c(rep(name1, 2), rep(name2, 2)),
                levels = c(name1, name2)
            ),
            status = factor(rep(c("common", "unique"), 2), levels = c("unique", "common")),
            nGenes = c(
                length(which(myDF$exp1 &  myDF$exp2)),
                length(which(myDF$exp1 & !myDF$exp2)),
                length(which(myDF$exp2 &  myDF$exp1)),
                length(which(myDF$exp2 & !myDF$exp1))
            ),
            stringsAsFactors = FALSE
        )
        myBarPlot <- ggplot(dfForBarplot, aes(x = experiment, y = nGenes, fill = status)) + geom_bar(stat = "identity") +
            labs(y = "number of genes")
        if (input$barPlotGuide) myBarPlot <- myBarPlot + geom_hline(yintercept = dfForBarplot$nGenes[1])
        myBarPlot <- myBarPlot + theme_bw(base_size = 20) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "top")
        return(myBarPlot)
    })

    dataForBarPlotTests <- reactive({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 3, "Loading..."))
        myTable <- data.frame(
            category = c(
                "common peaks",
                paste(input$barPlotSample1Defined, "only"),
                paste(input$barPlotSample2Defined, "only"),
                "gene space"
            ),
            nGenes = c(
                length(which(myDF$exp1 &  myDF$exp2)),
                length(which(myDF$exp1 & !myDF$exp2)),
                length(which(myDF$exp2 & !myDF$exp1)),
                nrow(myDF)
            )
        )
        return(myTable)
    })

    output$myBarPlot <- renderPlot({
        withProgress(value = 1, message = "Making bar plot...",
                     print(renderMyBarPlot())
        )
    })

    output$tabBarPlot <- renderTable(dataForBarPlotTests(), include.rownames = FALSE)

    output$barPlotCorrelation <- renderText({
        myDF <- getPairWiseData()
        validate(need(ncol(myDF) == 3, "Loading..."))
        cc <- cor.test(as.numeric(myDF$exp1), as.numeric(myDF$exp2))
        return(paste("Correlation coefficient:", round(cc$estimate, digits = 4), "  p-value:", signif(cc$p.value, digits = 4)))
    })

    output$barPlotJaccard <- renderText({
        myDF <- dataForBarPlotTests()
        myJaccard <- myDF$nGenes[1] / (myDF$nGenes[1] + myDF$nGenes[2] + myDF$nGenes[3])
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

    # to edit
    output$downloadExempleFile <- downloadHandler("ESR1_breastcancercells_PMID15023342.txt",
                                                   content = function(file) {
                                                       file.copy("www/ESR1_breastcancercells_PMID15023342.txt", file)
                                                   },
                                                   contentType = "text/tsv")

    # end of server -------------
})
