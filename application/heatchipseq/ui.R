# UI
source("data/server_adresses.R")
shinyUI(navbarPage("HeatChIPseq",

    tabPanel("Instructions",
             icon("home"),
             a("Back to main page", href = URL_HEATSTARSEQ),
             h2("Welcome !"),
             h3("Intsructions:"),
             p("1 - load a peak file."),
             p("2 - wait."),
             p("3 - I will write some interesting stuff here in the future."),
             h3("FAQ"),
             p("..."),
             h3("Info"),
             p("Source code available on GitHub (soon)."),
             p("Made by the Josh group (add link)."),
             p("How to cite:"),
             p("... unpublished")
             ),
    
    tabPanel("Use application",
            sidebarLayout(
                
                sidebarPanel(
                    icon("home"),
                    a("Back to main page", href = URL_HEATSTARSEQ),
                    h3("1 - Select a dataset"),
                    selectInput("dataset", label = NULL, choices = c(
                        "ENCODE TFBS ChIP-seq (human, hg19)",
                        "CODEX ChIP-seq (mouse, mm10)",
                        "CODEX ChIP-seq (human, hg19)"
                    )),
                    h3("2 - Load your data"),
                    p("Upload a bed-like peak file. Tab delimited, first three columns must be chromsome, peak start and peak end.
                      Maximum size: 10MB. Please, use the same reference genome version than the selected dataset."),
                    textInput("nameOfPeakFile", "Name of your experiment", value="my ChIP experiment"),
                    checkboxInput("header", "My peak file contains a header.", FALSE),
                    fileInput("peakFile", "Upload your peak file:", accept = "text/tab-separated-values"),
                    img(src = "legend_small.png"),
                    h3("3 - Heatmap customisation"),
                    checkboxInput("highlight", "Highlight my experiment in the heatmap", FALSE),
                    conditionalPanel(condition = "input.dataset == 'ENCODE TFBS ChIP-seq (human, hg19)'",
                        selectInput("TF", "Subset for TF(s) (empty to select all):",
                                    choices = unique(encode$annotation$TF)[order(unique(encode$annotation$TF))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("cells", "Subset for cell line(s) (empty to select all):",
                                    choices = levels(factor(encode$annotation$cell_line))[order(levels(factor(encode$annotation$cell_line)))],
                                    selected = NULL, multiple = TRUE)
                    ),
                    conditionalPanel(condition = "input.dataset == 'CODEX ChIP-seq (mouse, mm10)'",
                         selectInput("TF_m", "Subset for TF(s) (empty to select all):",
                                     choices = unique(codex$annotation$TF)[order(unique(codex$annotation$TF))],
                                     selected = NULL, multiple = TRUE),
                         radioButtons("filterCellsBy", label = "Filter by", choices = list("Cell type", "Cell subtype")),
                         conditionalPanel(condition = "input.filterCellsBy == 'Cell type'",
                                          selectInput("cell_types_m", "Subset for cell type(s) (empty to select all):",
                                                      choices = unique(codex$annotation[,"Cell type"])[order(unique(codex$annotation[,"Cell type"]))],
                                                      selected = NULL, multiple = TRUE)
                         ),
                         conditionalPanel(condition = "input.filterCellsBy == 'Cell subtype'",
                                          selectInput("cell_subtypes_m", "Subset for cell subtypes(s) (empty to select all):",
                                                      choices = unique(codex$annotation[,"Cell subtype"])[order(unique(codex$annotation[,"Cell subtype"]))],
                                                      selected = NULL, multiple = TRUE)
                         )             
                    ),
                    selectInput("hclustMethod", 
                                label = "Clusterisation method",
                                choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                selected = "complete"),
                    conditionalPanel(condition = "input.myPanels == 'Static Heatmap'",    
                                     numericInput("margin", label = "Label margin [1 - 50] (static heatmap only)", value = 20, min = 1, max = 50),
                                     numericInput("labCex", label = "Label size [0.1 - 3] (static heatmap only)", value = 1.2, min = 0.1, max = 3)
                    )
                ),
                
                mainPanel(
                    tabsetPanel(
                        tabPanel("My peaks", 
                                 dataTableOutput("tabUserPeaks"),
                                 downloadButton("downloadUserPeaks", label = "Save as tab delimited .txt")
                                 ),
                        tabPanel("Correlation table",
                                 dataTableOutput("tabUserCorrelationTable"),
                                 downloadButton("downloadUserCorrelationTable", label = "Save as tab delimited .txt")
                                 ),
                        tabPanel("Static Heatmap",
                                 plotOutput("myHeatmap", width = "950px", height = "950px"),
                                 downloadButton("downloadHMpng", label = "Save as png"),
                                 downloadButton("downloadHMpdf", label = "Save as pdf")
                                 ), 
                        tabPanel("Responsive Heatmap", plotlyOutput("myPlotlyHeatmap", width = "1000px", height = "1000px")),
                        tabPanel("TF", dataTableOutput("tabTF")),
                        tabPanel("Cells", dataTableOutput("tabCells"))
                        , id = "myPanels"
                    )
                )
                
            )
    ),
    
    theme = "bootstrap.css"
    
))

