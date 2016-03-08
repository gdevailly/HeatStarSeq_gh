source("data/server_adresses.R")

shinyUI(tagList(useShinyjs(), navbarPage("HeatCAGEseq",

    tabPanel("Instructions",
             icon("home"),
             a("Back to main page", href = URL_HEATSTARSEQ),
             includeHTML("www/Instructions_heatcageseq.html")
    ),

    tabPanel("Use application",
             sidebarLayout(

                sidebarPanel(
                    icon("home"),
                    a("Back to main page", href = URL_HEATSTARSEQ),
                    h3("1 - Select a dataset"),
                    selectInput("selectedDataset", label = NULL, choices = c(
                        "FANTOM5 (human, hg19)",
                        "FANTOM5 (mouse, mm10) (soon)"
                    )),
                    h3("2 - Load your data"),
                    p("Upload a 6 column bed file, with the score column scoring CAGE expression.
                      Maximum size: 10MB. Please, use the same reference genome version than the selected dataset."),
                    # p("You can download ",
                    # downloadLink("downloadExempleFile", label = " an example file"),
                    #   ". It is a human ESR1 (ERalpha) ChIP-seq experiment in MCF7 cells, kindly provided by Dr. Stromblad.
                    #   It corresponds to the MCF7_ERa_E2_ChIP sample from ",
                    #   a("this GEO dataset", href = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73320"), "."),
                    checkboxInput("header", "My file contains a header.", TRUE),
                    fileInput("cageFile", strong("Upload your bed file:"), accept = "text/tab-separated-values"),
                    textInput("nameOfCageFile", "Name of your experiment", value="my CAGE experiment"),
                    h3("3 - Plot customisation"),
                    checkboxInput("highlight", strong("Highlight my experiment in the heatmap"), FALSE),
                    # we addapt filtering widgets to the various datasets
                    div(id = "widgetForFantom5Human",
                        selectInput("f5h_cells", "Subset for cell type(s) (empty to select all):",
                                    choices = unique(fantom5_human_cage$annotation$tissue)[order(unique(fantom5_human_cage$annotation$tissue))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("f5h_isCellLine", "Cell lines or not? (empty to select all)",
                                    choices = c("all", "only cell line", "only non cell line"),
                                    selected = NULL, multiple = FALSE)
                    ),
                    selectInput("correlationCorrection",
                                label = "Uploaded experiment correlation correction:",
                                choices = list("None", "Linear scaling"),
                                selected = "None"),
                    sliderInput("maxCorrelation",
                                label = "Maximum expected correlation value for linear scaling correction",
                                min = 0.1, max = 1, value = 0.95, step = 0.01),
                    actionButton("advClustOptions", label = "Advance clustering options"),
                    div(id = "widgetForClustOptions",
                        selectInput("distOption", label = "Distance calculation:",
                                    choices = list("euclidean", "1 - correlations", "maximum", "manhattan", "canberra"),
                                    selected = 1),
                        selectInput("hclustMethod",
                                    label = "Clusterisation method",
                                    choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                    selected = "complete")
                    ),
                    div(id = "widgetForLabels",
                        sliderInput("labCex", label = "Sample name size", value = 1.2, min = 0.1, max = 3, step = 0.1),
                        sliderInput("margin", label = "Sample name margin", value = 20, min = 1, max = 50, step = 1)
                    )
                ),

                mainPanel(
                    tabsetPanel(
                        tabPanel("My bed file",
                                 dataTableOutput("tabUserCageFile"),
                                 downloadButton("downloadUserCageFile", label = "Save as tab delimited .txt")
                                 ),
                        tabPanel("Correlation table",
                                 dataTableOutput("tabUserCorrelationTable"),
                                 downloadButton("downloadUserCorrelationTable", label = "Save as tab delimited .txt")
                                 ),
                        tabPanel("Static heatmap",
                                 plotOutput("myHeatmap", width = "950px", height = "950px"),
                                 img(src = "legend_small.png"),
                                 downloadButton("downloadHMpng", label = "Save as png"),
                                 downloadButton("downloadHMpdf", label = "Save as pdf"),
                                 downloadButton("downloadHMsvg", label = "Save as svg")
                                 ),
                        tabPanel("Responsive heatmap",
                                 plotlyOutput("myPlotlyHeatmap", width = "1000px", height = "1000px"),
                                 img(src = "legend_small.png")
                        ),
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

