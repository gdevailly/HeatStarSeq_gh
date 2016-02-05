library(gplots)
library(readr)
library(preprocessCore)
options(shiny.maxRequestSize = 10*1024^2) # max file size, 10Mb

shinyServer(function(input, output) {
    
    # UI elements activations
    observe({
        if (is.null(input$expressionFile)) {
            shinyjs::hide("downloadUserExpressionFile")
            shinyjs::hide("downloadUserCorrelationTable")
        } else {
            shinyjs::show("downloadUserExpressionFile")
            shinyjs::show("downloadUserCorrelationTable")
        }
    })
    
    observe({
        shinyjs::hide("widgetForEncodeHuman")
        shinyjs::hide("widgetForBgeeHuman")
        shinyjs::hide("widgetForEncodeMouse")
        shinyjs::hide("widgetForBgeeMouse")
        if (input$dataset == "ENCODE RNA-seq (human)") {
            shinyjs::show("widgetForEncodeHuman")
        } else if (input$dataset == "Bgee RNA-seq (human)") {
            shinyjs::show("widgetForBgeeHuman")
        } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
            shinyjs::show("widgetForEncodeMouse")
        } else if (input$dataset == "Bgee RNA-seq (mouse)") {
            shinyjs::show("widgetForBgeeMouse")
        } else if (input$dataset == "other (soon)") {
            # fill this
        }
    })
    
    observe({
        if(input$myPanels == "Static Heatmap" | input$myPanels == "Tree") {
            shinyjs::show("widgetForLabels")
        } else {
            shinyjs::hide("widgetForLabels")
        }
    })
    
    # output computations
    getSelectedDataset <- reactive({
        withProgress(value = 1, message = "Loading dataset: ", detail = "removing old dataset", {
            load("data/encode_rnaseq_preload.RData")
            load("data/bgee_human_preload.RData")
            load("data/encode_mouse_rnaseq_preload.RData")
            load("data/bgee_mouse_preload.RData")
            setProgress(value = 1, detail = "loading new dataset")
            if (input$dataset == "ENCODE RNA-seq (human)") {
                load("data/encode_rnaseq.RData")
                dataset <- encode_rnaseq
            } else if (input$dataset == "Bgee RNA-seq (human)") {
                load("data/bgee_human.RData")
                dataset <- bgee_human
            } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
                load("data/encode_mouse_rnaseq.RData")
                dataset <- encode_mouse_rnaseq
            } else if (input$dataset == "Bgee RNA-seq (mouse)") {
                load("data/bgee_mouse.RData")
                dataset <- bgee_mouse
            } else if (input$dataset == "other (soon)") {
                # fill this
            }
            setProgress(value = 1, detail = "done!")
        })
        return(dataset)    
    })

    userExpressionFileAnalysis <- reactive({
        userExpressionFileName <- input$expressionFile
        if (is.null(userExpressionFileName)){
            userExpressionFile <- NULL
            userCorrelations <- NULL
            normUserCorrelations <- NULL
            normUserCorrelations_linear <- NULL
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
                    "Something went wrong: all genes have the same expression value (probably 0). Cannot calculate correlations.
                    \nThis is likely because:
                    \n-the wrong dataset is selected (ie you may have uploaded a mouse expression file, but the selected dataset is for human).
                    \n-the file formating is not recognized. It should be a two columns, tab delimited text file. First columns must contains the ensembl gene name (ie ENSG00000000003). Second column must contains expression values."
                ))
                # correlation calculation
                setProgress(value = 1, detail = "correlations calculation")
                userCorrelations <- cor(userExpressionFile, dataset$dataMatrix) %>% as.vector
                normUserCorrelations <- normalize.quantiles.use.target(
                    as.matrix(c(userCorrelations, 1)),
                    normalize.quantiles.determine.target(dataset$correlationMatrix)
                ) %>% as.vector
                if (max(userCorrelations) != 0) {
                    normUserCorrelations_linear <- userCorrelations * max(normUserCorrelations[-length(normUserCorrelations)]) / max(userCorrelations)
                }
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
            "quantNormCorrelations" = normUserCorrelations[-length(normUserCorrelations)],
            "linearNormCorrelations" = normUserCorrelations_linear
        ))
    })
    
    output$tabUserExpressionFile <- renderDataTable({
        validate(
            need(!is.null(input$expressionFile), "Upload an expression file, or click on a 'heatmap' tab to explore the dataset.")
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
            need(!is.null(input$expressionFile), "Upload an expression file, or click on a 'heatmap' tab to explore the dataset.")
        )
            data.frame(
                "Experiment" = getSelectedDataset()$annotation$name,
                #"Experiment" = paste(getSelectedDataset()$annotation[, "File accession"], getSelectedDataset()$annotation[, "Biosample term name"]),
                "Correlation" = userExpressionFileAnalysis()$correlations,
                "QuantNorm Correlation" = userExpressionFileAnalysis()$quantNormCorrelations,
                "LinearNorm Correlation" = userExpressionFileAnalysis()$linearNormCorrelations,
                stringsAsFactors = FALSE
            )[order(userExpressionFileAnalysis()$correlations, decreasing = TRUE), ]
    })
    
    output$tabUserCorrelationTable <- renderDataTable(getCorrelationTable())
    
    output$downloadUserCorrelationTable <- downloadHandler("heatRNAseq_correlations.txt",
                                                           content = function(file) {
                                                              write.table(getCorrelationTable(), file = file, row.names = FALSE, quote = FALSE, sep = "\t")
                                                           },
                                                           contentType = "text/tsv")
    
    subsetMatrix <- reactive({
        dataset <- getSelectedDataset()
        workingMatrix <- dataset$correlationMatrix
        keep<- 1:nrow(workingMatrix)
        # filtering is dataset-dependent...
        if (input$dataset == "ENCODE RNA-seq (human)") {
            if (is.null(input$cells)) {
                temp_cells <- unique(encode_rnaseq$annotation[, "Biosample term name"])
            } else {
                temp_cells <- input$cells
            }
            if (is.null(input$sampleType)) {
                temp_sampleType <- unique(encode_rnaseq$annotation[, "Biosample type"])
            } else {
                temp_sampleType <- input$sampleType
            }
            if (is.null(input$rnaExtract)) {
                temp_rnaExtract <- unique(encode_rnaseq$annotation[, "rnaFraction"])
            } else {
                temp_rnaExtract <- input$rnaExtract
            }
            keep <- which(
                dataset$annotation[, "Biosample term name"] %in% temp_cells &
                dataset$annotation[, "Biosample type"] %in% temp_sampleType &
                dataset$annotation[, "rnaFraction"] %in% temp_rnaExtract   
            )
        } else if (input$dataset == "Bgee RNA-seq (human)") {
            if (is.null(input$tissus_bgee_h)) {
                temp_tissus_bgee_h <- unique(bgee_human$annotation[, "Biosample term name"])
            } else {
                temp_tissus_bgee_h <- input$tissus_bgee_h
            }
            if (is.null(input$dvp_bgee_h)) {
                temp_dvp_bgee_h <- unique(bgee_human$annotation[, "Stage name"])
            } else {
                temp_dvp_bgee_h <- input$dvp_bgee_h
            }
            if (is.null(input$library_bgee_h)) {
                temp_library_bgee_h <- unique(bgee_human$annotation[, "Library type"])
            } else {
                temp_library_bgee_h <- input$library_bgee_h
            }
            keep <- which(
                dataset$annotation[, "Biosample term name"] %in% temp_tissus_bgee_h &
                dataset$annotation[, "Stage name"] %in% temp_dvp_bgee_h &
                dataset$annotation[, "Library type"] %in% temp_library_bgee_h   
            )
        } else if (input$dataset == "ENCODE RNA-seq (mouse)") {
            if (is.null(input$cells_encode_m)) {
                temp_cells_encode_m <- unique(encode_mouse_rnaseq$annotation[, "Biosample term name"])
            } else {
                temp_cells_encode_m <- input$cells_encode_m
            }
            if (is.null(input$sampleType_encode_m)) {
                temp_sampleType_encode_m <- unique(encode_mouse_rnaseq$annotation[, "Biosample type"])
            } else {
                temp_sampleType_encode_m <- input$sampleType_encode_m
            }
            keep <- which(
                dataset$annotation[, "Biosample term name"] %in% temp_cells_encode_m &
                dataset$annotation[, "Biosample type"] %in% temp_sampleType_encode_m
            )
        } else if (input$dataset == "Bgee RNA-seq (mouse)") {
            if (is.null(input$tissus_bgee_m)) {
                temp_tissus_bgee_m <- unique(bgee_mouse$annotation[, "Biosample term name"])
            } else {
                temp_tissus_bgee_m <- input$tissus_bgee_m
            }
            if (is.null(input$dvp_bgee_m)) {
                temp_dvp_bgee_m <- unique(bgee_mouse$annotation[, "Stage name"])
            } else {
                temp_dvp_bgee_m <- input$dvp_bgee_m
            }
            if (is.null(input$library_bgee_m)) {
                temp_library_bgee_m <- unique(bgee_mouse$annotation[, "Library type"])
            } else {
                temp_library_bgee_m <- input$library_bgee_m
            }
            keep <- which(
                dataset$annotation[, "Biosample term name"] %in% temp_tissus_bgee_m &
                dataset$annotation[, "Stage name"] %in% temp_dvp_bgee_m &
                dataset$annotation[, "Library type"] %in% temp_library_bgee_m   
            )
        } else { # modify accordingly for new datasets!

        }

        myLabels <- dataset$annotation$name
        if(length(keep) >= 2) {
            workingMatrix <- workingMatrix[keep, keep]
            myLabels <- myLabels[keep]
        }
        
        # merging user data in correlation matrix
        if(!is.null(userExpressionFileAnalysis()$correlations)) {
            userCorrelations <- userExpressionFileAnalysis()$correlations
            if (input$correlationCorrection == "Quantile normalisation") {
                userCorrelations <- userExpressionFileAnalysis()$quantNormCorrelations
            } else if (input$correlationCorrection == "Linear scaling") {
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
            myClust <- hclust(dist(myMat), method = input$hclustMethod)
            dendro <- as.dendrogram(myClust)
            setProgress(value = 1, detail = "done!")
        })
        return(list("dend" = dendro, "order" = myClust$order))
    })
    
    matAfterHighlight <- reactive({
        matAA <- subsetMatrix()
        if (input$highlight & !is.null(input$expressionFile)) {
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
        par(oldPar) # not all par() can be set
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
    
    output$myTree <- renderPlot(myRenderTreePlot())
    
    output$tabSampleList <- renderDataTable(getSelectedDataset()$annotation)
    
})
