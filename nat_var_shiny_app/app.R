library(shiny)
library(biomaRt)
library(leaflet)
library(RColorBrewer)



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
      actionButton(inputId="tab2.Submit", label = "Submit"),
      tags$hr(),
      tableOutput("tab2.GeneInfo"),
      h4("Nucleotide Diversity Statistic by Codon"),
      plotOutput("diversityPlot", brush="plot_brush", click="plot_click", height = 400),
      verbatimTextOutput("info"),
      downloadButton("downloadSNPData","Download"),
      tableOutput("Diversity_Table")
    ),
    
    tabPanel("SNP Mapping",
             ## Tab 3
             textInput(inputId = "tab3.Gene", label = "Type a single gene locus in the box below",
                       value = "AT1G80490"),
             actionButton(inputId="tab3.Submit", label = "Submit"),
             tags$hr(),
             leafletOutput("tab3.map"),
             tableOutput("tab3Table")
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

  tab2.Genes <- eventReactive(input$tab2.Submit, {
      #gene Info for gene on tab 2, updates on 'submit' button press
    names <- parseInput(input$plotGene)
    genes <- getGeneInfo(names[1])
    return(genes)
  })

  output$tab2.GeneInfo <- renderTable(tab2.Genes())
    #rendered table of Gene info

  tab2.tableData <- reactive({load_tab_2_Data(tab2.Genes())})
    #SNP reactive data

  output$Diversity_Table <- renderTable(tab2.tableData(), rownames=TRUE)
    #render table of diversity data

  output$downloadSNPData <- downloadHandler(filename="yourData.csv", content = function(file) {
      #download for diversity data
      write.csv(tab2.tableData(), file)
  })

  output$diversityPlot <- renderPlot(plotPi(tab2.tableData()))
    #plot output

  output$info <- renderPrint({
    brushedPoints(tab2.tableData(), input$plot_brush)
  })

  
  ##                            _________
  ##                           /  tab3   \
  ## --------------------------           ----------------------------
  ## Tab 3 stuff:
  
  tab3.Genes <- eventReactive(input$tab3.Submit, {
    #gene Info for gene on tab 2, updates on 'submit' button press
    names <- parseInput(input$tab3.Gene)
    genes <- getGeneInfo(names[1])
    return(genes)
  })
  
  
  tab3.tidyData <- reactive({
    tidyVCF <- VCFByTranscript(tab3.Genes()[1, ], strains)
    data <- tidyVCF$dat[tidyVCF$dat$gt_GT != "0|0",]
    # Parse the EFF field
    data <- parseEFF(tidyVCF = data, Transcript_ID = tab3.Genes()$transcript_ID[1])
    
    # calculate diversity
    data <- Nucleotide_diversity(data)
    
    return(data)
  })
  
  
  tab3.labeledSNPs <- reactive({
    
    data <- tab3.tidyData()
    keyPOS <- unique(data[which(data$Diversity >= 0.5*max(data$Diversity)), "POS"])
    
    keydata <- data[data$POS %in% keyPOS, ]
    keydata_labeled <- label_bySNPs(keydata)
    return(keydata_labeled)
    
  })
  
  output$tab3.map <- renderLeaflet({
    
    mapdata <- tab3.labeledSNPs()
    
    # Reorganize to plot NA's underneath non NA's
    mapdata <- rbind(mapdata[is.na(mapdata$SNPs), ], mapdata[!is.na(mapdata$SNPs), ])
    
    # make a field with text to be displayed when clicking on a marker
    mapdata$popup <- paste("EcoID:",  mapdata$Indiv,"Name:", mapdata$Name, " SNPs:", mapdata$SNPs)
    
    pallet <- colorFactor(palette="Set1", domain=mapdata$SNPs )
    
    map <- leaflet()
    
    map <- addProviderTiles(map, providers$Stamen.TonerLite,
                     options = providerTileOptions(noWrap = TRUE))
    
    map <- addCircleMarkers(map, data=mapdata, color= ~pallet(SNPs), 
                            radius=6, popup= ~popup, stroke=FALSE, fillOpacity=0.85)
    
    map <- addLegend(map, position="bottomright", pal=pallet, 
                     values=mapdata$SNPs, title="Marker Colors", opacity=1)

    
    return(map)

  })
  
  output$tab3Table <- renderTable(tab3.labeledSNPs()[1:350, ], digits=5)

  
}


shinyApp(ui = ui, server = server)

