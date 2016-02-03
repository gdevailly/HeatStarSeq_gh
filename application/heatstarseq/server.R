library(shiny)

shinyServer(function(input, output) {
    
    output$RNAseqExempleTable <- renderDataTable(
        data.frame(
            geneID = c("ENSG00000134046", "ENSG00000141644", "ENSG00000169057", "ENSG00000174282", "ENSG00000187098"),
            tpm = c("120.12", "0", "85.24", "0.54", "42")
        ),
        options = list(paging = FALSE, searching = FALSE)
    )
    
    output$ChIPseqExempleTable <- renderDataTable(
        data.frame(
            chr = c("chr1", "chr1", "chr4", "chr12", "chrX"),
            start = c("125423", "8545032", "4523698", "854120", "2458750"),
            end =   c("125891", "8546254", "4524785", "854870", "2459872")
        ),
        options = list(paging = FALSE, searching = FALSE)
    )
    
    output$CAGEseqExempleTable <- renderDataTable(
        data.frame(
            chr = c("chr2", "chr6", "chr6", "chr17", "chr21"),
            start = c("25486325", "5896321", "223541", "5012035", "960032"),
            end =   c("25486487", "5896380", "223602", "5012100", "960098"),
            name = c("CAGEpeak_1", "CAGEpeak_2", "CAGEpeak_3", "CAGEpeak_4" , "CAGEpeak_5"),
            rpm = c(458.12, 25.03, 1.23, 45.3, 8.70),
            strand = c("+", "+", "-", "+", "-")
        ),
        options = list(paging = FALSE, searching = FALSE)
    )
    
})