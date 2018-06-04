source("data/server_adresses.R")

shinyUI(tagList(
    useShinyjs(),
    tags$head(list(
        includeScript("www/events.js")
    )),
    navbarPage(
        a(div(icon("home"), "Heat*seq"), href = URL_HEATSTARSEQ),

        tabPanel("Use application",
                 sidebarLayout(

                    sidebarPanel(
                        h2("HeatGeneList"),
                        h3("1 - Select a dataset"),
                        selectInput("dataset", label = NULL, choices = c(
                            "Knock-out (human)",
                            "Knock-out (mouse)"
                        ), selected = "Knock-out (human)"),
                        h3("2 - Load your data (optional)"),
                        actionButton("fileFormatInstructions", label = "File formating instructions"),
                        div(id = "div_fileFormatInstructions",
                            p("Upload a tab delimited text file of at least one columns. First column should contain gene name (i.e. HUGO symbols),
                               Other columns will be ignore. Case should not matter. Maximum size: 10MB."),
                            HTML("<p>
    			                First lines of your file should look like this:
                                <table style=\"width:80%\">
                                <tr>
                                <th>geneID</th>
                                </tr>
                                <tr>
                                <td>MBD2</td>
                                </tr>
                                <tr>
                                <td>CLDN6</td>
                                </tr>
                                <tr>
                                <td>CDKN1A</td>
                                </tr>
                                <tr>
                                <td>NTN1</td>
                                </tr>
                                <tr>
                                <td>UNC5B</td>
                                </tr>
                                </table>
                                </p>"),
                            p("You can download ",
                              downloadLink("downloadExempleFile", label = " an example file"),
                              ". It is a somewhat random human gene list.")
                        ),
                        radioButtons("fileToUse", label = NULL, choices = c("Upload your gene list", "Use the example file")),
                        div(id = "div_fileupload",
                            fileInput("expressionFile", label = "Choose a file:", accept = "text/tab-separated-values"),
                            checkboxInput("header", strong("The expression file contains a header."), TRUE)
                        ),
                        div(id = "div_exampleInUse", "The example file is a Human gene list. Please select only human datasets."),
                        textInput("nameOfExpressionFile", "Name of your experiment:", value="my gene list"),
                        h3("3 - Plot customization"),
                        div(id = "widgetsForHeatmap",
                            checkboxInput("highlight", strong("Highlight my gene list in the heatmap."), TRUE),
                            # we adapt filtering widgets to the various datasets
                            div(id = "widgetForHumanKo",
                                 selectInput("cell_human_ko", "Subset for cell type (empty to select all):",
                                             choices = sort(unique(human_ko$annotation$cell_type)),
                                             selected = NULL, multiple = TRUE),
                                 selectInput("tf_human_ko", "Subset for knock-out gene (empty to select all):",
                                             choices = sort(unique(human_ko$annotation$ko_gene)),
                                             selected = NULL, multiple = TRUE)
                            ),
                            div(id = "widgetForMouseKo",
                                selectInput("cell_mouse_ko", "Subset for cell type (empty to select all):",
                                            choices = sort(unique(mouse_ko$annotation$cell_type)),
                                            selected = NULL, multiple = TRUE),
                                selectInput("tf_mouse_ko", "Subset for knock-out gene (empty to select all):",
                                            choices = sort(unique(mouse_ko$annotation$ko_gene)),
                                            selected = NULL, multiple = TRUE)
                            ),
                            selectInput("correlationCorrection",
                                        label = "Uploaded experiment correlation correction:",
                                        choices = list("None", "Linear scaling"),
                                        selected = "None"),
                            div(id = "div_maxCorrelation", sliderInput("maxCorrelation",
                                        label = "Maximum expected correlation value for linear scaling correction:",
                                        min = 0.1, max = 1, value = 0.95, step = 0.01)),
                            actionButton("advClustOptions", label = "Advanced clustering options"),
                            div(id = "widgetForClustOptions",
                                selectInput("distOption", label = "Distance calculation:",
                                            choices = list("euclidean", "1 - Pearson's correlation coefficient", "maximum", "manhattan", "canberra"),
                                            selected = 1),
                                selectInput("hclustMethod",
                                            label = "Clustering method:",
                                            choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                            selected = "complete")
                            ),
                            div(id = "widgetForLabels",
                                radioButtons("labelOption", label = "Label size:", choices = list("Automatic", "Adjust manually")),
                                div(id = "widgetForLabelsManual",
                                    sliderInput("labCex", label = "Sample name size:", value = 1.2, min = 0.1, max = 3, step = 0.1),
                                    sliderInput("margin", label = "Sample name margin:", value = 20, min = 1, max = 50, step = 1)
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
                                colourpicker::colourInput("col_col1", label = "Colour 1:", value = "blue"),
                                sliderInput("col_val1", label = "Value 1:", value = -0.5, min = -1, max = 1, step = 0.05),
                                colourpicker::colourInput("col_col2", label = "Colour 2:", value = "white"),
                                sliderInput("col_val2", label = "Value 2:", value = 0, min = -1, max = 1, step = 0.05),
                                colourpicker::colourInput("col_col3", label = "Colour 3:", value = "red"),
                                sliderInput("col_val3", label = "Value 3:", value = 0.5, min = -1, max = 1, step = 0.05),
                                colourpicker::colourInput("col_col4", label = "Colour 4:", value = "black"),
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
                            tabPanel("My gene list",
                                     dataTableOutput("tabUserExpressionFile"),
                                     downloadButton("downloadUserGeneList", label = "Save as tab delimited .txt")
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
                 h1("HeatGeneLits", align = "center"),
                 includeHTML("www/Instructions_heatgenelist.html"),
                 p(a(img(src = "The_Roslin_Institute_logo.gif"), href = "http://www.roslin.ed.ac.uk/"),
                   a(img(src = "BBSRC_logo.gif"), href = "http://www.bbsrc.ac.uk/"), align = "center")
        ),

        theme = "bootstrap.css",

        windowTitle = "HeatGeneList"

)))
