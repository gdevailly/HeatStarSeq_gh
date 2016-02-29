library(shiny)
source("data/server_adresses.R")
shinyUI(
    navbarPage(
        "Heat*seq",
        tabPanel(
            "Home",
            h1("Heat*seq", align = "center"),
            h3("Correlation heatmaps for hight-througput sequencing experiments", align = "center"),
            fluidRow(
               column(4,
                      h2("Gene Expression data", align = "center"),
                      p(a(
                          actionButton("rnaButton", "HeatRNAseq", class = "btn-success", style = "font-size:150%"),
                          href = URL_HEATRNASEQ
                      ), align = "center"),
                      "Data file format you could upload:",
                      p("Tab-delimited file of expression values (2 columns)."),
                      dataTableOutput("RNAseqExempleTable")
               ),
               column(4,
                      h2("ChIP-seq data", align = "center"),
                      p(a(
                          actionButton("chipButton", "HeatChIPseq", class = "btn-info", style = "font-size:150%"),
                          href = URL_HEATCHIPSEQ
                      ), align = "center"),
                      "Data file format you could upload:",
                      p("Tab-delimited bed-like file of peak coordinates (minimum 3 columns)."),
                      dataTableOutput("ChIPseqExempleTable")
               ),
               column(4,
                      h2("CAGE-seq data", align = "center"),
                      p(a(
                          actionButton("cageButton", "HeatCAGEseq", class = "btn-warning", style = "font-size:150%"),
                          href = URL_HEATCAGESEQ
                      ), align = "center"),
                      "Data file format you could upload:",
                      p("Tab-delimited bed-like file of peak coordinates, strand and expression value (minimum 5 columns)."),
                      dataTableOutput("CAGEseqExempleTable")
               )
            )

        ),

        tabPanel(
            "About",
            includeHTML("www/Instructions_heatstarseq.html"),
            a(img(src = "The_Roslin_Institute_logo.gif"), href = "http://www.roslin.ed.ac.uk/"),
            a(img(src = "BBSRC_logo.gif"), href = "http://www.bbsrc.ac.uk/")

        ),

        theme = "bootstrap.css"
    )
)