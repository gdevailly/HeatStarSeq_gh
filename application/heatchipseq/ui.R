# UI
source("data/server_adresses.R")
shinyUI(tagList(useShinyjs(), navbarPage("HeatChIPseq",

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
                    selectInput("selectedDataset", label = NULL, choices = c(
                        "ENCODE TFBS ChIP-seq (human, hg19)",
                        "CODEX ChIP-seq (human, hg19)",
                        "CODEX ChIP-seq (mouse, mm10)"
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
                    div(id = "widgetForEncodeHuman",
                        selectInput("TF", "Subset for TF(s) (empty to select all):",
                                    choices = unique(encode$annotation$tf)[order(unique(encode$annotation$tf))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("cells", "Subset for cell line(s) (empty to select all):",
                                    choices = levels(factor(encode$annotation$cellLine))[order(levels(factor(encode$annotation$cellLine)))],
                                    selected = NULL, multiple = TRUE)
                    ),
                    div(id = "widgetForCodexHuman",
                        selectInput("TF_ch", "Subset for TF(s) (empty to select all):",
                                    choices = unique(codex_human_chip$annotation$tf)[order(unique(codex_human_chip$annotation$tf))],
                                    selected = NULL, multiple = TRUE),
                        radioButtons("filterCellsBy_ch", label = "Filter by", choices = list("Cell type", "Cell subtype")),
                        conditionalPanel(condition = "input.filterCellsBy_ch == 'Cell type'",
                                         selectInput("cell_types_ch", "Subset for cell type(s) (empty to select all):",
                                                     choices = unique(codex_human_chip$annotation$cellType)[order(unique(codex_human_chip$annotation$cellType))],
                                                     selected = NULL, multiple = TRUE)
                        ),
                        conditionalPanel(condition = "input.filterCellsBy_ch == 'Cell subtype'",
                                         selectInput("cell_subtypes_ch", "Subset for cell subtypes(s) (empty to select all):",
                                                     choices = unique(codex_human_chip$annotation$cellSubtype)[order(unique(codex_human_chip$annotation$cellSubtype))],
                                                     selected = NULL, multiple = TRUE)
                        )
                    ),
                    div(id = "widgetForCodexMouse",
                         selectInput("TF_m", "Subset for TF(s) (empty to select all):",
                                     choices = unique(codex$annotation$tf)[order(unique(codex$annotation$tf))],
                                     selected = NULL, multiple = TRUE),
                         radioButtons("filterCellsBy", label = "Filter by", choices = list("Cell type", "Cell subtype")),
                         conditionalPanel(condition = "input.filterCellsBy == 'Cell type'",
                                          selectInput("cell_types_m", "Subset for cell type(s) (empty to select all):",
                                                      choices = unique(codex$annotation$cellType)[order(unique(codex$annotation$cellType))],
                                                      selected = NULL, multiple = TRUE)
                         ),
                         conditionalPanel(condition = "input.filterCellsBy == 'Cell subtype'",
                                          selectInput("cell_subtypes_m", "Subset for cell subtypes(s) (empty to select all):",
                                                      choices = unique(codex$annotation$cellSubtype)[order(unique(codex$annotation$cellSubtype))],
                                                      selected = NULL, multiple = TRUE)
                         )
                    ),
                    selectInput("hclustMethod",
                                label = "Clusterisation method",
                                choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                selected = "complete"),
                    div(id = "widgetForLabels",
                        sliderInput("margin", label = "Sample name margin", value = 20, min = 1, max = 50, step = 1),
                        sliderInput("labCex", label = "Sample name size", value = 1.2, min = 0.1, max = 3, step = 0.1)
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
                        tabPanel("Tree",
                                 plotOutput("myTree", width = "500px", height = "950px"),
                                 downloadButton("downloadTreePng", label = "Save as png"),
                                 downloadButton("downloadTreePdf", label = "Save as pdf"),
                                 downloadButton("downloadTreeSvg", label = "Save as svg")
                        ),
                        tabPanel("Samples metadata",
                                 dataTableOutput("tabSampleList"),
                                 downloadButton("downloadDatasetTable", label = "Save as tab delimited .txt")
                        )
                        , id = "myPanels"
                    )
                )

            )
    ),

    theme = "bootstrap.css"

)))

