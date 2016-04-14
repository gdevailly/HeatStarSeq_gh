source("data/server_adresses.R")
shinyUI(tagList(useShinyjs(), navbarPage(a(div(icon("home"), "Heat*seq"), href = URL_HEATSTARSEQ),

    tabPanel("Use application",
             sidebarLayout(

                sidebarPanel(
                    h2("HeatRNAseq"),
                    h3("1 - Select a dataset"),
                    selectInput("dataset", label = NULL, choices = c(
                        "Bgee RNA-seq (human)",
                        "Blueprint RNA-seq (human)",
                        "Roadmap Epigenomics RNA-seq (human)",
                        "GTEX summary (human)",
                        "ENCODE RNA-seq (human)",
                        "Bgee RNA-seq (mouse)",
                        "ENCODE RNA-seq (mouse)",
                        "Flybase RNA-seq (drosophila)"
                    ), selected = "Bgee RNA-seq (mouse)"),
                    h3("2 - Load your data (optional)"),
                    actionButton("fileFormatInstructions", label = "File formating instructions"),
                    div(id = "div_fileFormatInstructions",
                        p("Upload a tab delimited text file of at least two columns. First column should contain gene id (Ensembl or Flybase),
                       second column should contain normalised expression value (i.e. FPKM or TPM). Maximum size: 10MB."),
                        HTML("<p>
			                First lines of your file should look like this:
                            <table style=\"width:80%\">
                            <tr>
                            <th>geneID</th>
                            <th>tpm</th>
                            </tr>
                            <tr>
                            <td>ENSG00000134046</td>
                            <td>120.12</td>
                            </tr>
                            <tr>
                            <td>ENSG00000141644</td>
                            <td>0</td>
                            </tr>
                            <tr>
                            <td>ENSG00000169057</td>
                            <td>85.24</td>
                            </tr>
                            <tr>
                            <td>ENSG00000174282</td>
                            <td>0.54</td>
                            </tr>
                            <tr>
                            <td>ENSG00000187098</td>
                            <td>42</td>
                            </tr>
                            </table>
                            </p>"),
                        p("You can download ",
                          downloadLink("downloadExempleFile", label = " an example file"),
                          ". It is a mouse RNA-seq experiment from the brain, with a header.")
                    ),
                    radioButtons("fileToUse", label = NULL, choices = c("Upload your expression file", "Use the example file")),
                    div(id = "div_fileupload", fileInput("expressionFile", label = "Choose a file:", accept = "text/tab-separated-values")),
                    div(id = "div_exampleInUse", "The example file is a mouse brain RNA-seq, with a header. Please select only mouse datasets."),
                    checkboxInput("header", strong("The expression file contains a header."), TRUE),
                    textInput("nameOfExpressionFile", "Name of your experiment:", value="my RNA-seq"),
                    h3("3 - Plot customization"),
                    checkboxInput("highlight", strong("Highlight my experiment in the heatmap."), TRUE),
                    # we adapt filtering widgets to the various datasets
                    div(id = "widgetForBgeeHuman",
                         selectInput("tissus_bgee_h", "Subset for tissue (empty to select all):",
                                     choices = unique(bgee_human$annotation$tissue)[order(unique(bgee_human$annotation$tissue))],
                                     selected = NULL, multiple = TRUE),
                         selectInput("dvp_bgee_h", "Subset for developmental stage (empty to select all):",
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
                    div(id = "widgetForRoadmapHuman",
                        selectInput("celltype_roadmap_h", "Cell type (empty to select all):",
                                    choices = unique(roadmap_rnaseq$annotation$name)[order(unique(roadmap_rnaseq$annotation$name))],
                                    selected = NULL, multiple = TRUE)
                    ),
                    div(id = "widgetForGtexSmall",
                        selectInput("celltype_gtex_small_h", "Tissue (empty to select all):",
                                    choices = unique(gtex_small$annotation$name)[order(unique(gtex_small$annotation$name))],
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
                    div(id = "widgetForFlybase",
                        selectInput("sample_flybase", "Subset for sample (empty to select all):",
                                    choices = unique(flybase_rnaseq$annotation$name)[order(unique(flybase_rnaseq$annotation$name))],
                                    selected = NULL, multiple = TRUE),
                        selectInput("library_flybase", "Subset for parent library (empty to select all):",
                                    choices = unique(flybase_rnaseq$annotation$parentLibrary)[order(unique(flybase_rnaseq$annotation$parentLibrary))],
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
                                    choices = list("euclidean", "1 - correlations", "maximum", "manhattan", "canberra"),
                                    selected = 1),
                        selectInput("hclustMethod",
                                    label = "Clustering method:",
                                    choices = list("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                                    selected = "complete")
                    ),
                    div(id = "widgetForLabels",
                        radioButtons("labelOption", label = "Label size:", choices = list("Automatic", "Adjust manualy")),
                        div(id = "widgetForLabelsManual",
                            sliderInput("labCex", label = "Sample name size:", value = 1.2, min = 0.1, max = 3, step = 0.1),
                            sliderInput("margin", label = "Sample name margin:", value = 20, min = 1, max = 50, step = 1)
                        )
                    ),
                    div(id = "div_widgetHMoptions",
                        selectInput("showDend", "Show which dendrogram(s)?",
                                    choices = c("both", "row", "column", "none"),
                                    selected = "both", multiple = FALSE),
                        selectInput("showLabels", "Show which labels?",
                                    choices = c("both", "row", "column", "none"),
                                    selected = "both", multiple = FALSE)
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
                        , id = "myPanels", selected = "Responsive heatmap"
                    )
                )

            )
    ),

    tabPanel("Instructions",
             h1("HeatRNAseq", align = "center"),
             includeHTML("www/Instructions_heatrnaseq.html"),
             p(a(img(src = "The_Roslin_Institute_logo.gif"), href = "http://www.roslin.ed.ac.uk/"),
               a(img(src = "BBSRC_logo.gif"), href = "http://www.bbsrc.ac.uk/"), align = "center")
    ),

    theme = "bootstrap.css",

    windowTitle = "HeatRNAseq"

)))