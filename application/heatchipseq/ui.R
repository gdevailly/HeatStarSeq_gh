source("data/server_adresses.R")

shinyUI(tagList(useShinyjs(), navbarPage("HeatChIPseq",

    tabPanel("Instructions",
             icon("home"),
             a("Back to main page", href = URL_HEATSTARSEQ),
             includeHTML("www/Instructions_heatchipseq.html")
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
                        "CODEX ChIP-seq (mouse, mm10)",
                        "modEncode TF ChIP-seq (drosophila, r5)"
                    )),
                    h3("2 - Load your data"),
                    p("Upload a bed-like peak file. Tab delimited, first three columns must be chromsome, peak start and peak end.
                      Maximum size: 10MB. Please, use the same reference genome version than the selected dataset."),
                    p("You can download ",
                    downloadLink("downloadExempleFile", label = " an example file"),
                      ". It is a human ESR1 (ERalpha) ChIP-seq experiment in MCF7 cells, kindly provided by Dr. Stromblad.
                      It corresponds to the MCF7_ERa_E2_ChIP sample from ",
                      a("this GEO dataset", href = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73320"), "."),
                    checkboxInput("header", "My peak file contains a header.", TRUE),
                    fileInput("peakFile", strong("Upload your peak file:"), accept = "text/tab-separated-values"),
                    textInput("nameOfPeakFile", "Name of your experiment", value="my ChIP experiment"),
                    h3("3 - Plot customization"),
                    checkboxInput("highlight", strong("Highlight my experiment in the heatmap"), FALSE),
                    # we adapt filtering widgets to the various datasets
                    div(id = "widgetForEncodeHuman",
                        selectInput("TF", "Subset for TF(s) (empty to select all):",
                                    choices = unique(encode$annotation$tf)[order(unique(encode$annotation$tf))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("cells", "Subset for cell line(s) (empty to select all):",
                                    choices = unique(encode$annotation$cellLine)[order(unique(encode$annotation$cellLine))],
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
                    div(id = "widgetForModEncodeD",
                        selectInput("tf_med", "Subset for antibody (empty to select all):",
                                    choices = unique(modEncodeD_ChIPseq$annotation$antibody)[order(unique(modEncodeD_ChIPseq$annotation$antibody))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("stage_med", "Subset for developmental stage (empty to select all):",
                                    choices = unique(modEncodeD_ChIPseq$annotation$devStage)[order(unique(modEncodeD_ChIPseq$annotation$devStage))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("strain_med", "Subset for strain (empty to select all):",
                                    choices = unique(modEncodeD_ChIPseq$annotation$strain)[order(unique(modEncodeD_ChIPseq$annotation$strain))],
                                    selected = NULL, multiple = TRUE)
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
                                    label = "Clustering method",
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
                        tabPanel("My peaks",
                                 dataTableOutput("tabUserPeaks"),
                                 downloadButton("downloadUserPeaks", label = "Save as tab delimited .txt")
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
                                 downloadButton("downloadHMsvg", label = "Save as svg"),
                                 downloadButton("downloadHMdata", label = "Export data as tab delimited .txt")
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

