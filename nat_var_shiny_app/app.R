library(shiny)
library(biomaRt)

ui <- fluidPage(
  headerPanel("Arabidopsis Natural Variation Webtool"),
  "This app will provide an inteface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset",
  tags$h5('style'="color:red", "this app is a work in progress"),
  tabsetPanel(

    tabPanel("SNP Stats",
        ## Tab 1
      textAreaInput(inputId = "gene_ids", label = "Type a list of gene loci in the box below, separated by commas ",
                    width = 600, height = 250, value = "AT3G62980, AT4G03190, AT3G26810, AT1G12820, AT4G24390, AT5G49980, AT2G39940" ),
      actionButton(inputId="STATS_submit", label = "Submit"),
      tags$hr(),
      tableOutput("table"),
      downloadButton("downloadData","Download"),
      tableOutput("SNPStats_Table")
    ),

    tabPanel("Diversity Plot",
        ## Tab 2
      textInput(inputId = "plotGene", label = "Type a single gene locus in the box below",
                value = "AT1G80490"),
      actionButton(inputId="tab2Submit", label = "Submit"),
      tags$hr(),
      tableOutput("tab2GeneInfo"),
      h4("Nucleotide Diversity Statistic by Codon"),
      plotOutput("diversityPlot", brush="plot_brush", click="plot_click", height = 400),
      verbatimTextOutput("info"),
      downloadButton("downloadSNPData","Download"),
      tableOutput("Diversity_Table")
    )


  ),

  "THIS IS THE FOOTER"
)
#=================================================================

source("VCF_Utils.R")
source("Strains_and_Gene_Families.R")

parseInput <- function (textIn) {
  names <- str_extract_all(textIn, "AT[1-5]G[0-9]{5}")
  return (names[[1]])
}

load_tab_2_Data <- function (geneInfo){
  tab2VCF <- VCFByTranscript(geneInfo[1, ], strains)
  tab2data <- tab2VCF$dat
  tab2data <- parseEFF(tab2data, geneInfo[1, "transcript_ID"])
  tab2data <- Nucleotide_diversity(tab2data)
  
  coding_variants <- coding_Diversity_Plot(tab2data)
  
  return(coding_variants)
  
}

#writeData <- function (table, file){
#  write.csv(tableData, file)
#}

plotPi <- function(unique_coding_variants) {
  plot <- ggplot(unique_coding_variants, aes(x=Codon_Number,y=Diversity, colour=Effect)) +
    geom_point() +
    scale_y_log10(breaks=c(0.0001, 0.001, 0.01, 0.1),limits=c(0.0001, 1)) +
    #scale_colour_manual(values=c(synonymous_diversity="blue", missense_diversity="red")) +
    ylab("nucleotide diversity, log scale")
  return(plot)
  
}


server <- function(input, output){

  ##   _________
  ##  /  tab1   \
  ##             --------------------------------------------------
  ## Tab 1 stuff:

  observeEvent(input$STATS_submit,{
    names <- parseInput(input$gene_ids)
    genes <- getGeneInfo(names)
    req(genes != FALSE)
    output$table <- renderTable(genes)
    output$downloadData <- downloadHandler(filename="yourData.csv", content = function(file) {
      write.csv(genes, file)
    })
    SNPStats <- polymorphTable(genes, strains)
    output$SNPStats_Table <- renderTable(SNPStats, rownames=TRUE, digits=5)
  })

  ##                 _________
  ##                /  tab2   \
  ## ---------------           -------------------------------------
  ## Tab 2 stuff:

  tab2Genes <- eventReactive(input$tab2Submit, {
      #gene Info for gene on tab 2, updates on 'submit' button press
    names <- parseInput(input$plotGene)
    genes <- getGeneInfo(names[1])
    return(genes)
  })

  output$tab2GeneInfo <- renderTable(tab2Genes())
    #rendered table of Gene info

  tab2tableData <- reactive({load_tab_2_Data(tab2Genes())})
    #SNP reactive data

  output$Diversity_Table <- renderTable(tab2tableData(), rownames=TRUE)
    #render table of diversity data

  output$downloadSNPData <- downloadHandler(filename="yourData.csv", content = function(file) {
      #download for diversity data
      write.csv(tab2tableData(), file)
  })

  output$diversityPlot <- renderPlot(plotPi(tab2tableData()))
    #plot output

  output$info <- renderPrint({
    brushedPoints(tab2tableData(), input$plot_brush)
  })

}


shinyApp(ui = ui, server = server)

