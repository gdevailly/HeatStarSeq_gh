source("data/server_adresses.R")
shinyUI(tagList(
    useShinyjs(),
    tags$head(
        includeScript("www/google_analytics.js"),
        includeScript("www/events.js")
    ),
    navbarPage(
        a(div(icon("home"), "Heat*seq"), href = URL_HEATSTARSEQ),

        tabPanel("Use application",
                 sidebarLayout(

                    sidebarPanel(
                        h2("HeatCAGEseq"),
                        h3("1 - Select a dataset"),
                        selectInput("selectedDataset", label = NULL, choices = c(
                            "FANTOM5 (human, hg19)",
                            "FANTOM5 (mouse, mm9)"
                        ), selected = "FANTOM5 (mouse, mm9)"),
                        h3("2 - Load your data  (optional)"),
                        actionButton("fileFormatInstructions", label = "File formating instructions"),
                        div(id = "div_fileFormatInstructions",
                            p("Upload a 6 column bed file, with the score column scoring CAGE expression.
                              Maximum size: 10MB. Please, use the same reference genome version than the selected dataset."),
                            HTML("<p>
    	                         First lines of your file should look like this:
                                 <table style=\"width:100%\">
                                 <tr>
                                 <th>chr</th>
                                 <th>start</th>
                                 <th>end</th>
                                 <th>name</th>
                                 <th>rpm</th>
                                 <th>strand</th>
                                 </tr>
                                 <tr>
                                 <td>chr2</td>
                                 <td>25486325</td>
                                 <td>25486487</td>
                                 <td>CAGEpeak_1</td>
                                 <td>458.12</td>
                                 <td>+</td>
                                 </tr>
                                 <tr>
                                 <td>chr6</td>
                                 <td>5896321</td>
                                 <td>5896380</td>
                                 <td>CAGEpeak_2</td>
                                 <td>25.03</td>
                                 <td>+</td>
                                 </tr>
                                 <tr>
                                 <td>chr6</td>
                                 <td>223541</td>
                                 <td>223602</td>
                                 <td>CAGEpeak_3</td>
                                 <td>1.23</td>
                                 <td>-</td>
                                 </tr>
                                 <tr>
                                 <td>chr17</td>
                                 <td>5012035</td>
                                 <td>5012100</td>
                                 <td>CAGEpeak_4</td>
                                 <td>45.3</td>
                                 <td>+</td>
                                 </tr>
                                 <tr>
                                 <td>chr21</td>
                                 <td>960032</td>
                                 <td>960098</td>
                                 <td>CAGEpeak_5</td>
                                 <td>8.70</td>
                                 <td>-</td>
                                 </tr>
                                 </table>
                                 </p>"),
                            p("You can download ",
                              downloadLink("downloadExempleFile", label = " an example file"),
                              ". It is a mouse CAGE experiment from the liver (mm9), without a header.")
                        ),
                        radioButtons("fileToUse", label = NULL, choices = c("Upload your result file", "Use the example file")),
                        div(id = "div_fileupload",
                            fileInput("cageFile", strong("Upload your bed file:"), accept = "text/tab-separated-values"),
                            checkboxInput("header", "My file contains a header.", TRUE)
                        ),
                        div(id = "div_exampleInUse", "The example file is a mouse CAGE experiment from the liver. Please select only the mouse dataset."),
                        textInput("nameOfCageFile", "Name of your experiment", value="my CAGE experiment"),
                        h3("3 - Plot customization"),
                        div(id = "widgetsForHeatmap",
                            checkboxInput("highlight", strong("Highlight my experiment in the heatmap"), TRUE),
                            # we adapt filtering widgets to the various datasets
                            div(id = "widgetForFantom5Human",
                                selectInput("f5h_cells", "Cell type(s) (empty to select all):",
                                            choices = sort(unique(fantom5_human_cage$annotation$tissue)),
                                            selected = NULL, multiple = TRUE),
                                selectInput("f5h_isCellLine", "Sample origin:",
                                            choices = c("all", "only cell line", "only non cell line"),
                                            selected = NULL, multiple = FALSE)
                            ),
                            div(id = "widgetForFantom5Mouse",
                                selectInput("f5m_cells", "Cell type(s) (empty to select all):",
                                            choices = sort(unique(fantom5_mouse_cage$annotation$tissue)),
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
                                    sliderInput("labCex", label = "Sample name size:", value = 1.2, min = 0.1, max = 3, step = 0.1),
                                    sliderInput("margin", label = "Sample name margin:", value = 20, min = 1, max = 50, step = 1)
                                )                    ),
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
                                sliderInput("col_val1", label = "Value 1:", value = 0.25, min = -1, max = 1, step = 0.05),
                                colourInput("col_col2", label = "Colour 2:", value = "white"),
                                sliderInput("col_val2", label = "Value 2:", value = 0.5, min = -1, max = 1, step = 0.05),
                                colourInput("col_col3", label = "Colour 3:", value = "red"),
                                sliderInput("col_val3", label = "Value 3:", value = 0.75, min = -1, max = 1, step = 0.05),
                                colourInput("col_col4", label = "Colour 4:", value = "black"),
                                sliderInput("col_val4", label = "Value 4:", value = 1, min = -1, max = 1, step = 0.05),
                                actionButton("applyColoursOptions", label = "Apply colour changes")
                            )
                        ),

                        div(id = "widgetForPairwisePlots",
                            uiOutput("scatterPlotSample1"),
                            uiOutput("scatterPlotSample2"),
                            selectInput("scatterPlotType", label = "Plot type:", choices = list("XY", "MA")),
                            selectInput("scatterPlotDataScaling", label = "Data scaling:", choices = list(
                                "none",
                                "log10(e + 1)",
                                "log(e + 1)",
                                "log2(e + 1)",
                                "asinh(e)",
                                "1/(1 + e)"
                            ), selected = "log10(e + 1)"),
                            checkboxInput("scatterPlotRegression", label = "Add regression line (blue)", value = TRUE),
                            checkboxInput("scatterPlotGuide", label = "Add guide line (red)", value = TRUE)
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
                                     plotOutput("myScatterPlot", width = "500px", height = "500px"),
                                     textOutput("scatterPlotMetricsPearson"),
                                     textOutput("scatterPlotMetricsSpearman"),
                                     downloadButton("downloadScatterPlotPng", label = "Save as png"),
                                     downloadButton("downloadScatterPlotPdf", label = "Save as pdf"),
                                     downloadButton("downloadScatterPlotSvg", label = "Save as svg"),
                                     downloadButton("downloadScatterPlotData", label = "Export data as tab delimited .txt")
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
                 h1("HeatCAGEseq", align = "center"),
                 includeHTML("www/Instructions_heatcageseq.html"),
                 p(a(img(src = "The_Roslin_Institute_logo.gif"), href = "http://www.roslin.ed.ac.uk/"),
                   a(img(src = "BBSRC_logo.gif"), href = "http://www.bbsrc.ac.uk/"), align = "center")
        ),

        theme = "bootstrap.css",

        windowTitle = "HeatCAGEseq"

)))
