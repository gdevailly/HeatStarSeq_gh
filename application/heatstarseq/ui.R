library(shiny)
source("data/server_adresses.R")
shinyUI(tagList(
    navbarPage(
        "",
        tabPanel(
            "Heat*seq Home",
            h1("Heat*seq", align = "center"),
            h3("An interactive web tool for high-throughput sequencing experiment comparison with public data", align = "center"),
            h3(" "),
            fluidRow(
               column(2),
               column(2,
                      h2("RNA-seq", align = "center"),
                      p(a(
                          actionButton("rnaButton", "HeatRNAseq", class = "btn-success", style = "font-size:150%"),
                          href = URL_HEATRNASEQ
                      ), align = "center")

               ),
               column(1),
               column(2,
                      h2("ChIP-seq", align = "center"),
                      p(a(
                          actionButton("chipButton", "HeatChIPseq", class = "btn-info", style = "font-size:150%"),
                          href = URL_HEATCHIPSEQ
                      ), align = "center")
               ),
               column(1),
               column(2,
                      h2("CAGE", align = "center"),
                      p(a(
                          actionButton("cageButton", "HeatCAGEseq", class = "btn-warning", style = "font-size:150%"),
                          href = URL_HEATCAGESEQ
                      ), align = "center")
               ),
               column(2)
            ),
            fluidRow(
                column(2),
                column(2,
                       h2("Gene lists", align = "center"),
                       p(a(
                           actionButton("geneButton", "HeatGeneList", class = "btn-danger", style = "font-size:150%"),
                           href = URL_HEATGENELIST
                       ), align = "center")
                ),
                column(8)
            ),
            h1(img(src = "heatmap.svg", width = 700), align = "center")


        ),

        tabPanel(
            "About",
            includeHTML("www/Instructions_heatstarseq.html"),
            p(a(img(src = "The_Roslin_Institute_logo.gif"), href = "http://www.roslin.ed.ac.uk/"),
            a(img(src = "BBSRC_logo.gif"), href = "http://www.bbsrc.ac.uk/"), align = "center")

        ),

        theme = "bootstrap.css", windowTitle = "Heat*seq"
    )
))
