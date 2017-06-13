source("data/server_adresses.R")
shinyUI(tagList(
    useShinyjs(),
    tags$head(list(
        includeScript("www/google_analytics.js"),
        includeScript("www/events.js")
    )),
    navbarPage(
        a(div(icon("home"), "Heat*seq"), href = URL_HEATSTARSEQ),

        tabPanel("Use application",
                 sidebarLayout(

                    sidebarPanel(
                        h2("HeatChIPseq"),
                        h3("1 - Select a dataset"),
                        selectInput("selectedDataset", label = NULL, choices = c(
                            "ENCODE TFBS ChIP-seq (human, hg19)",
                            "CODEX ChIP-seq (human, hg19)",
                            "ENCODE TFBS ChIP-seq (mouse, mm10)",
                            "CODEX ChIP-seq (mouse, mm10)",
                            "modEncode TF ChIP-seq (drosophila, r5)"
                        ), selected = "ENCODE TFBS ChIP-seq (human, hg19)"),
                        h3("2 - Load your data (optional)"),
                        actionButton("fileFormatInstructions", label = "File formating instructions"),
                        div(id = "div_fileFormatInstructions",
                            p("Upload a bed-like peak file. Tab delimited, first three columns should be chromosome, peak start and peak end.
                              Maximum size: 10MB. Please, use the same reference genome version as the selected dataset."),
                            HTML("<p>
    	                         First lines of your file should look like this:
                                 <table style=\"width:50%\">
                                 <tr>
                                 <th>chr</th>
                                 <th>start</th>
                                 <th>end</th>
                                 </tr>
                                 <tr>
                                 <td>chr1</td>
                                 <td>125423</td>
                                 <td>125891</td>
                                 </tr>
                                 <tr>
                                 <td>chr1</td>
                                 <td>8545032</td>
                                 <td>8546254</td>
                                 </tr>
                                 <tr>
                                 <td>chr4</td>
                                 <td>4523698</td>
                                 <td>4524785</td>
                                 </tr>
                                 <tr>
                                 <td>chr12</td>
                                 <td>854120</td>
                                 <td>854870</td>
                                 </tr>
                                 <tr>
                                 <td>chrX</td>
                                 <td>2458750</td>
                                 <td>2459872</td>
                                 </tr>
                                 </table>
                                 </p>
                                 "),
                            p("You can download ",
                            downloadLink("downloadExempleFile", label = " an example file"),
                              ". It is a human ESR1 (ERalpha) ChIP-seq experiment in MCF7 cells, without a header.")
                        ),
                        radioButtons("fileToUse", label = NULL, choices = c("Upload your peak file", "Use the example file")),
                        div(id = "div_fileupload",
                            fileInput("peakFile", strong("Choose a file:"), accept = "text/tab-separated-values"),
                            checkboxInput("header", "My peak file contains a header.", FALSE)
                        ),
                        div(id = "div_exampleInUse", "The example file is a ESR1 (ERalpha) ChIP-seq in MCF7 cells. Please select only human datasets."),
                        textInput("nameOfPeakFile", "Name of your experiment", value="my ChIP experiment"),
                        h3("3 - Plot customization"),
                        div(id = "widgetsForHeatmap",
                            checkboxInput("highlight", strong("Highlight my experiment in the heatmap"), TRUE),
                            # we adapt filtering widgets to the various datasets
                            div(id = "widgetForEncodeHuman",
                                selectInput("TF", "Transcription factor(s) (empty to select all):",
                                            choices = unique(encode$annotation$tf)[order(unique(encode$annotation$tf))],
                                            selected = NULL, multiple = TRUE),
                                selectInput("cells", "Cell line(s) (empty to select all):",
                                            choices = unique(encode$annotation$cellLine)[order(unique(encode$annotation$cellLine))],
                                            selected = NULL, multiple = TRUE)
                            ),
                            div(id = "widgetForCodexHuman",
                                selectInput("TF_ch", "Transcription factor(s) (empty to select all):",
                                            choices = unique(codex_human_chip$annotation$tf)[order(unique(codex_human_chip$annotation$tf))],
                                            selected = NULL, multiple = TRUE),
                                radioButtons("filterCellsBy_ch", label = "Filter by", choices = list("Cell type", "Cell subtype")),
                                conditionalPanel(condition = "input.filterCellsBy_ch == 'Cell type'",
                                                 selectInput("cell_types_ch", "Cell type(s) (empty to select all):",
                                                             choices = sort(unique(codex_human_chip$annotation$cellType)),
                                                             selected = NULL, multiple = TRUE)
                                ),
                                conditionalPanel(condition = "input.filterCellsBy_ch == 'Cell subtype'",
                                                 selectInput("cell_subtypes_ch", "Cell subtypes(s) (empty to select all):",
                                                             choices = sort(unique(codex_human_chip$annotation$cellSubtype)),
                                                             selected = NULL, multiple = TRUE)
                                )
                            ),
                            div(id = "widgetForEncodeMouse",
                                selectInput("TF_m_enc", "Transcription factor(s) (empty to select all):",
                                            choices = unique(encode_mouse$annotation$target)[order(unique(encode_mouse$annotation$target))],
                                            selected = NULL, multiple = TRUE),
                                selectInput("cells_m_enc", "Cell line(s) (empty to select all):",
                                            choices = unique(encode_mouse$annotation$cell_type)[order(unique(encode_mouse$annotation$cell_type))],
                                            selected = NULL, multiple = TRUE)
                            ),
                            div(id = "widgetForCodexMouse",
                                 selectInput("TF_m", "Transcription factor(s) (empty to select all):",
                                             choices = unique(codex$annotation$tf)[order(unique(codex$annotation$tf))],
                                             selected = NULL, multiple = TRUE),
                                 radioButtons("filterCellsBy", label = "Filter by", choices = list("Cell type", "Cell subtype")),
                                 conditionalPanel(condition = "input.filterCellsBy == 'Cell type'",
                                                  selectInput("cell_types_m", "SCell type(s) (empty to select all):",
                                                              choices = sort(unique(codex$annotation$cellType)),
                                                              selected = NULL, multiple = TRUE)
                                 ),
                                 conditionalPanel(condition = "input.filterCellsBy == 'Cell subtype'",
                                                  selectInput("cell_subtypes_m", "Cell subtypes(s) (empty to select all):",
                                                              choices = sort(unique(codex$annotation$cellSubtype)),
                                                              selected = NULL, multiple = TRUE)
                                 )
                            ),
                            div(id = "widgetForModEncodeD",
                                selectInput("tf_med", "Antibody (empty to select all):",
                                            choices = sort(unique(modEncodeD_ChIPseq$annotation$antibody)),
                                            selected = NULL, multiple = TRUE),
                                selectInput("stage_med", "Developmental stage (empty to select all):",
                                            choices = sort(unique(modEncodeD_ChIPseq$annotation$devStage)),
                                            selected = NULL, multiple = TRUE),
                                selectInput("strain_med", "Strain (empty to select all):",
                                            choices = sort(unique(modEncodeD_ChIPseq$annotation$strain)),
                                            selected = NULL, multiple = TRUE)
                            ),
                            selectInput("correlationCorrection",
                                        label = "Uploaded experiment correlation correction:",
                                        choices = list("None", "Linear scaling"),
                                        selected = "None"),
                            div(id = "div_maxCorrelation", sliderInput("maxCorrelation",
                                        label = "Maximum expected correlation value for linear scaling correction",
                                        min = 0.1, max = 1, value = 0.95, step = 0.01)),
                            actionButton("advClustOptions", label = "Advanced clustering options"),
                            div(id = "widgetForClustOptions",
                                selectInput("distOption", label = "Distance calculation:",
                                            choices = list("euclidean", "1 - Pearson's correlation coefficient", "maximum", "manhattan", "canberra"),
                                            selected = 1),
                                selectInput("hclustMethod",
                                            label = "Clustering method",
                                            choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                            selected = "complete")
                            ),
                            div(id = "widgetForLabels",
                                radioButtons("labelOption", label = "Label size:", choices = list("Automatic", "Adjust manually")),
                                div(id = "widgetForLabelsManual",
                                    sliderInput("labCex", label = "Sample name size", value = 0.1, min = 0.1, max = 3, step = 0.1),
                                    sliderInput("margin", label = "Sample name margin", value = 20, min = 1, max = 50, step = 1)
                                )
                            ),

                            div(id = "div_widgetHMoptions",
                                selectInput("showDend", "Show dendrogram(s)?",
                                            choices = c("both", "row", "column", "none"),
                                            selected = "both", multiple = FALSE),
                                selectInput("showLabels", "Show labels?",
                                            choices = c("both", "row", "column", "none"),
                                            selected = "both", multiple = FALSE)
                            ),

                            actionButton("coloursOptions", label = "Customise colours"),
                            div(id = "div_colourOptions",
                                colourInput("col_col1", label = "Colour 1:", value = "blue"),
                                sliderInput("col_val1", label = "Value 1:", value = -0.5, min = -1, max = 1, step = 0.05),
                                colourInput("col_col2", label = "Colour 2:", value = "white"),
                                sliderInput("col_val2", label = "Value 2:", value = 0, min = -1, max = 1, step = 0.05),
                                colourInput("col_col3", label = "Colour 3:", value = "red"),
                                sliderInput("col_val3", label = "Value 3:", value = 0.5, min = -1, max = 1, step = 0.05),
                                colourInput("col_col4", label = "Colour 4:", value = "black"),
                                sliderInput("col_val4", label = "Value 4:", value = 1, min = -1, max = 1, step = 0.05),
                                actionButton("applyColoursOptions", label = "Apply colour changes")
                            )
                        ),

                        div(id = "widgetForPairwisePlots",
                            uiOutput("barPlotSample1"),
                            uiOutput("barPlotSample2"),
                            checkboxInput("barPlotGuide", label = "Add guide line", value = TRUE)
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
                                     plotOutput("colourKey1", width = "250px", height = "100px"),
                                     downloadButton("downloadHMpng", label = "Save as png"),
                                     downloadButton("downloadHMpdf", label = "Save as pdf"),
                                     downloadButton("downloadHMsvg", label = "Save as svg"),
                                     downloadButton("downloadHMdata", label = "Export data as tab delimited .txt")
                            ),
                            tabPanel("Responsive heatmap",
                                     plotlyOutput("myPlotlyHeatmap", width = "1000px", height = "1000px"),
                                     plotOutput("colourKey2", width = "250px", height = "100px")
                            ),
                            tabPanel("Tree",
                                     plotOutput("myTree", width = "500px", height = "950px"),
                                     downloadButton("downloadTreePng", label = "Save as png"),
                                     downloadButton("downloadTreePdf", label = "Save as pdf"),
                                     downloadButton("downloadTreeSvg", label = "Save as svg")
                            ),
                            tabPanel("Pairwise plot",
                                     plotOutput("myBarPlot", width = "400px", height = "600px"),
                                     textOutput("barPlotCorrelation"),
                                     textOutput("barPlotJaccard"),
                                     tableOutput("tabBarPlot"),
                                     downloadButton("downloadBarPlotPng", label = "Save as png"),
                                     downloadButton("downloadBarPlotPdf", label = "Save as pdf"),
                                     downloadButton("downloadBarPlotSvg", label = "Save as svg"),
                                     downloadButton("downloadBarPlotData", label = "Export data as tab delimited .txt")
                            ),
                            tabPanel("Samples metadata",
                                     dataTableOutput("tabSampleList"),
                                     downloadButton("downloadDatasetTable", label = "Save as tab delimited .txt")
                            ), id = "myPanels", selected = "Static heatmap"
                        )
                    )

                )
        ),

        tabPanel("Instructions",
                 h1("HeatChIPseq", align = "center"),
                 includeHTML("www/Instructions_heatchipseq.html"),
                 p(a(img(src = "The_Roslin_Institute_logo.gif"), href = "http://www.roslin.ed.ac.uk/"),
                   a(img(src = "BBSRC_logo.gif"), href = "http://www.bbsrc.ac.uk/"), align = "center")
        ),

        theme = "bootstrap.css",

        windowTitle = "HeatChIPseq"

)))