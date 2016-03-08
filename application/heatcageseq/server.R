library(readr)
library(svglite)
library(cba)

options(shiny.maxRequestSize = 10*1024^2) # max file size, 10Mb

shinyServer(function(input, output) {

    # UI elements activations
    observe({
        if (is.null(input$peakFile)) {
            shinyjs::hide("downloadUserCageFile")
            shinyjs::hide("downloadUserCorrelationTable")
        } else {
            shinyjs::show("downloadUserCageFile")
            shinyjs::show("downloadUserCorrelationTable")
        }
    })

    observe({
        shinyjs::hide("widgetForFantom5Human")
        if (input$selectedDataset == "FANTOM5 (human, hg19)") {
            shinyjs::show("widgetForFantom5Human")
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
    })

    observe({
        junkVar <- input$advClustOptions
        shinyjs::toggle("widgetForClustOptions")
    })

    # output computations
    getSelectedDataset <- reactive({
        withProgress(value = 1, message = "Loading dataset: ", detail = "removing old dataset", {
            # this remove from memmory non-selected datasets
            load("data/fantom5_human_cage_preload.RData")
            setProgress(value = 1, detail = "loading new dataset")
            if(input$selectedDataset == "FANTOM5 (human, hg19)") {
                load("data/fantom5_human_cage.RData")
                dataset <- fantom5_human_cage
            }
            setProgress(value = 1, detail = "done!")
        })
        return(dataset)
    })

    userCageFileAnalysis_1 <- reactive({
        userCageFileName <- input$cageFile
        if (is.null(userCageFileName)){
            userCageFile <- NULL
            userCorrelations <- NULL
        } else {
            withProgress(value = 1, message = "Usercage file: ", detail = "reading file", {
                userCageFile <- read_tsv(userCageFileName$datapath, col_names = input$header)[, 1:6]
                setProgress(value = 1, detail = "intersecting CAGE regions")
                colnames(userCageFile) <- c("chr", "start", "end", "name", "score", "strand")
                userCageFileGR <- with(userCageFile, GRanges(chr, IRanges(start, end), strand = strand, score = score))
                dataset <- getSelectedDataset()
                # warnings about missing chromosomes in one of the 2 sets. Let's assume user now what he/her is uploading...
                # Which region overlaps?
                # we ignore peaks presents only in user peak lists.
                suppressWarnings(userOverlap <- findOverlaps(dataset$regionMetaData, userCageFileGR, select = "arbitrary"))
                userCuratedValues <- numeric(nrow(dataset$dataMatrix))
                userCuratedValues[!is.na(userOverlap)] <- userCageFile[userOverlap[!is.na(userOverlap)], "score"]
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(userCuratedValues,
                                         dataset$dataMatrix
                                        ) %>% as.vector
                setProgress(value = 1, detail = "done!")
            })
        }
        return(list(
            "peaks" = userCageFile,
            "correlations" = userCorrelations,
            "linearNormCorrelations" = NULL # to be fill below
        ))
    })

    userCageFileAnalysis <- reactive({
        userCageFileData <- userCageFileAnalysis_1()
        userCorrelations <- userCageFileData$correlations
        if (!is.null(input$cageFile)){
            if (input$maxCorrelation != 0) {
                userCageFileData$linearNormCorrelations <- userCorrelations * input$maxCorrelation / max(userCorrelations)
            }
        }
        return(userCageFileData)
    })

    output$tabUserCageFile <- renderDataTable({
        validate(
            need(!is.null(input$cageFile), "Upload an bed file, or click on a 'heatmap' tab to explore the dataset.")
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
            need(!is.null(input$cageFile), "Upload an bed file, or click on a 'heatmap' tab to explore the dataset.")
        )
            data.frame(
                "experiment" = getSelectedDataset()$annotation$name,
                "correlation" = userCageFileAnalysis()$correlations,
                "scaledCorrelation" = userCageFileAnalysis()$linearNormCorrelations,
                stringsAsFactors = FALSE
            )[order(userCageFileAnalysis()$correlations, decreasing = TRUE),]
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
        }

        myLabels <- dataset$annotation$name # annotation must have a name column, with _unique_ elements
        validate(
            need(length(keep) >= 3, "Less than 3 experiments match your criteria. Please selecet more experiments.")
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
            if (input$distOption == "1 - correlations") {
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
        if (input$highlight & !is.null(input$cageFile)) {
            # we increase user cor value by 2 to do the color trick
            matAA$mat[, ncol(matAA$mat)] <- matAA$mat[, ncol(matAA$mat)] + 2
            matAA$mat[nrow(matAA$mat), 1:(ncol(matAA$mat) - 1)] <- matAA$mat[nrow(matAA$mat), 1:(ncol(matAA$mat) - 1)] + 2
        }
        return(matAA)
    })

    myRenderPlot <- function() {
        matData <- matAfterHighlight()
        clusterDat <- doTheClustering()
        heatmap(matData$mat,
                Rowv = clusterDat$dend,
                Colv = "Rowv",
                scale = "none",
                labRow = matData$myLabels,
                labCol = matData$myLabels,
                margins = rep(input$margin, 2),
                cexRow = input$labCex,
                cexCol = input$labCex,
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
        clusterDat <- doTheClustering()
        oldPar <- list(mar = par()$mar, cex = par()$cex)
        par(mar = c(5, 4, 4, input$margin), cex = input$labCex)
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

    # output$downloadExempleFile <- downloadHandler("chipseq_human_hg19_GSM1890761_ERa_peaks_no_header.bed",
    #                                               content = function(file) {
    #                                                   file.copy("www/chipseq_human_hg19_GSM1890761_ERa_peaks_no_header.bed", file)
    #                                               },
    #                                               contentType = "text/tsv")

})
