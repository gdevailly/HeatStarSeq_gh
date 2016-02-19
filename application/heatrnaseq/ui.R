source("data/server_adresses.R")
shinyUI(tagList(useShinyjs(), navbarPage("HeatRNAseq",

    tabPanel("Instructions",
             icon("home"),
             a("Back to main page", href = URL_HEATSTARSEQ),
             h2("Welcome !"),
             h3("Intsructions:"),
             p("1 - load a exprssion file."),
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
                        "Bgee RNA-seq (human)",
                        "Blueprint RNA-seq (human)",
                        "ENCODE RNA-seq (human)",
                        "Bgee RNA-seq (mouse)",
                        "ENCODE RNA-seq (mouse)"
                    )),
                    h3("2 - Load your data (optional)"),
                    p("Upload a tab delimited text file of at least two columns. First column must contains ensemble gene id, second column must contains normalised expression value
                      (ie FPKM or TPM)"),
                    textInput("nameOfExpressionFile", "Name of your experiment", value="my RNA-seq"),
                    checkboxInput("header", strong("My expression file contains a header."), FALSE),
                    fileInput("expressionFile", "Upload your expression file:", accept = "text/tab-separated-values"),
                    h3("3 - Plot customisation"),
                    checkboxInput("highlight", strong("Highlight my experiment in the heatmap"), FALSE),
                    # we addapt filtering widgets to the various datasets
                    div(id = "widgetForBgeeHuman",
                         selectInput("tissus_bgee_h", "Subset for tissue (empty to select all):",
                                     choices = unique(bgee_human$annotation$tissue)[order(unique(bgee_human$annotation$tissue))],
                                     selected = NULL, multiple = TRUE),
                         selectInput("dvp_bgee_h", "Subset for devlopemntal stage (empty to select all):",
                                     choices = unique(bgee_human$annotation$stage)[order(unique(bgee_human$annotation$stage))],
                                     selected = NULL, multiple = TRUE),
                         selectInput("library_bgee_h", "Subset for library type (empty to select all):",
                                     choices = unique(bgee_human$annotation$libraryType)[order(unique(bgee_human$annotation$libraryType))],
                                     selected = NULL, multiple = TRUE)
                    ),
                    div(id = "widgetForBlueprintHuman",
                        selectInput("tissue_blueprint_h", "Extract from (empty to select all):",
                                    choices = unique(blueprint_rnaseq$annotation$sampleSource)[order(unique(blueprint_rnaseq$annotation$sampleSource))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("celltype_blueprint_h", "Cell type (empty to select all):",
                                    choices = unique(blueprint_rnaseq$annotation$cellType)[order(unique(blueprint_rnaseq$annotation$cellType))],
                                    selected = NULL, multiple = TRUE)
                    ),
                    div(id = "widgetForEncodeHuman",
                        selectInput("cells", "Subset for tissue/cell line (empty to select all):",
                                    choices = unique(encode_rnaseq$annotation$tissue)[order(unique(encode_rnaseq$annotation$tissue))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("sampleType", "Subset for sample type (empty to select all):",
                                    choices = unique(encode_rnaseq$annotation$sampleType)[order(unique(encode_rnaseq$annotation$sampleType))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("rnaExtract", "Subset for RNA fraction (empty to select all):",
                                    choices = unique(encode_rnaseq$annotation$rnaFraction)[order(unique(encode_rnaseq$annotation$rnaFraction))],
                                    selected = NULL, multiple = TRUE)
                    ),
                    div(id = "widgetForEncodeMouse",
                                     selectInput("cells_encode_m", "Subset for tissue/cell line (empty to select all):",
                                                 choices = unique(encode_mouse_rnaseq$annotation$tissue)[order(unique(encode_mouse_rnaseq$annotation$tissue))],
                                                 selected = NULL, multiple = TRUE),
                                     selectInput("sampleType_encode_m", "Subset for sample type (empty to select all):",
                                                 choices = unique(encode_mouse_rnaseq$annotation$sampleType)[order(unique(encode_mouse_rnaseq$annotation$sampleType))],
                                                 selected = NULL, multiple = TRUE)
                    ),
                    div(id = "widgetForBgeeMouse",
                                     selectInput("tissus_bgee_m", "Subset for tissue (empty to select all):",
                                                 choices = unique(bgee_mouse$annotation$tissue)[order(unique(bgee_mouse$annotation$tissue))],
                                                 selected = NULL, multiple = TRUE),
                                     selectInput("dvp_bgee_m", "Subset for devlopemntal stage (empty to select all):",
                                                 choices = unique(bgee_mouse$annotation$stage)[order(unique(bgee_mouse$annotation$stage))],
                                                 selected = NULL, multiple = TRUE),
                                     selectInput("library_bgee_m", "Subset for library type (empty to select all):",
                                                 choices = unique(bgee_mouse$annotation$libraryType)[order(unique(bgee_mouse$annotation$libraryType))],
                                                 selected = NULL, multiple = TRUE)
                    ),
                    selectInput("correlationCorrection",
                                label = "Uploaded experiment correlation correction:",
                                choices = list("None", "Linear scaling"),
                                selected = "None"),
                    sliderInput("maxCorrelation",
                                label = "Maximum expected correlation value for Linear scaling correction",
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
                        sliderInput("margin", label = "Sample name margin", value = 20, min = 1, max = 50, step = 1),
                        sliderInput("labCex", label = "Sample name size", value = 1.2, min = 0.1, max = 3, step = 0.1)
                    )
                ),

                mainPanel(
                    tabsetPanel(
                        tabPanel("My expression file",
                                 dataTableOutput("tabUserExpressionFile"),
                                 downloadButton("downloadUserExpressionFile", label = "Save as tab delimited .txt")
                                 ),
                        tabPanel("Correlation table",
                                 dataTableOutput("tabUserCorrelationTable"),
                                 downloadButton("downloadUserCorrelationTable", label = "Save as tab delimited .txt")
                                 ),
                        tabPanel("Static Heatmap",
                                 plotOutput("myHeatmap", width = "950px", height = "950px"),
                                 img(src = "legend_small.png"),
                                 downloadButton("downloadHMpng", label = "Save as png"),
                                 downloadButton("downloadHMpdf", label = "Save as pdf"),
                                 downloadButton("downloadHMsvg", label = "Save as svg")
                                 ),
                        tabPanel("Responsive Heatmap",
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