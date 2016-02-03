library(shiny)
source("data/server_adresses.R")
shinyUI(
    navbarPage(
        "Heat*seq",
        tabPanel(
            "Home",
            h1("Heat*seq", align = "center"),
            h3("Correlation heatmaps for hight-througput sequencing expreriments", align = "center"),
            fluidRow(
               column(4,
                      h2("RNA-seq data", align = "center"),
                      p(a( 
                          actionButton("chipButton", "HeatRNAseq", class = "btn-success", style = "font-size:150%"),
                          href = URL_HEATRNASEQ
                      ), align = "center"),
                      t("Data file format you could upload:"),
                      p("Tab-delimited file of expression values (2 columns)."),
                      dataTableOutput("RNAseqExempleTable")
               ),
               column(4,
                      h2("ChIP-seq data", align = "center"),
                      p(a( 
                          actionButton("chipButton", "HeatChIPseq", class = "btn-info", style = "font-size:150%"),
                          href = URL_HEATCHIPSEQ
                      ), align = "center"),
                      t("Data file format you could upload:"),
                      p("Tab-delimited bed-like file of peak coordinates (minimum 3 columns)."),
                      dataTableOutput("ChIPseqExempleTable")
               ),
               column(4,
                      h2("CAGE-seq data", align = "center"),
                      p(a( 
                          actionButton("cageButton", "HeatCAGEseq", class = "btn-warning", style = "font-size:150%"),
                          href = URL_HEATCAGESEQ
                      ), align = "center"),
                      t("Data file format you could upload:"),
                      p("Tab-delimited bed-like file of peak coordinates, strand and expression value (minimum 5 columns)."),
                      dataTableOutput("CAGEseqExempleTable")
               )
            )
           
        ),
        
        tabPanel(
            "About",
            t("Write meaningfull things here.")
            
        ),
        
        theme = "bootstrap.css"
    )
)