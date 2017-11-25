library(shiny)
library(biomaRt)
library(leaflet)
library(RColorBrewer)



ui <- fluidPage(
  headerPanel("Arabidopsis Natural Variation Webtool"),
  "This app will provide an interface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset",
  tags$h5('style'="color:red", "this app is a work in progress"),
  tabsetPanel(

    tabPanel("SNP Stats",
        ## Tab 1
      textAreaInput(inputId = "gene_ids", label = "Type a list of gene loci in the box below, separated by commas ",
                    width = 600, height = 75, value = "AT3G62980, AT4G03190, AT3G26810, AT1G12820, AT4G24390, AT5G49980, AT2G39940" ),
      actionButton(inputId="STATS_submit", label = "Submit"),
      tags$hr(),
      tags$br(),
      h3("Selected Gene Information"),
      downloadButton("tab1.downloadGeneInfo","Download Content of Table Below"),
      tableOutput("table"),
      tags$hr(),
      tags$br(),
      tags$h3("SNP Type and Diversity Statistics"),
      downloadButton("tab1.downloadStats","Download Content of Table Below"),
      tableOutput("SNPStats_Table")
    ),

    tabPanel("Diversity Plot",
        ## Tab 2
      textInput(inputId = "plotGene", label = "Type a single gene locus in the box below",
                value = "AT1G80490"),
      actionButton(inputId="tab2.Submit", label = "Submit"),
      tags$hr(),
      tags$br(),
      tableOutput("tab2.GeneInfo"),
      h3("Plot of Nucleotide Diversity Statistic by Codon"),
      h5("click and drag a box accross the plot below to see details on specific points"),
      plotOutput("diversityPlot", brush="plot_brush", click="plot_click", height = 400),
      verbatimTextOutput("info"),
      tags$hr(),
      tags$br(),
      tags$h3("Diversity Plot Data"),
      downloadButton("tab2.downloadSNPData","Download Content of Table Below"),
      tableOutput("Diversity_Table")
    ),
    
    tabPanel("SNP Mapping",
             ## Tab 3
             textInput(inputId="tab3.Gene", label="Type a single gene locus in the box below",
                       value="AT1G80490"),
             actionButton(inputId="tab3.Submit", label="Submit"),
             sliderInput(inputId="tab3.filter_value", label="Nucleotide diversity filter limit",
                         min=0, max=.7, value=0.01),
             radioButtons("tab3.SNPtype", "Type of SNP to mark", 
                          choices=c("All", "Coding", "Missense")),
             verbatimTextOutput("tab3.debug"),
             tags$hr(),
             tags$br(),
             tags$h3("Accession Map"),
             tags$h5("Zoom with scroll wheel, click and drag to pan, click on individual point to see details"),
             leafletOutput("tab3.map"),
             tags$hr(),
             tags$br(),
             tags$h3("Map Data"),
             downloadButton("tab3.downloadMapData","Download Content of Table Below"),
             tableOutput("tab3Table")
    ),
    
    tabPanel("Accessions and Mutations",
             ## Tab 4
             textInput(inputId="tab4.Gene", label="Type list of gene ID's",
                       value="AT1G80490"),
             actionButton(inputId="tab4.Submit", label="Get Data"),
             
             selectInput(inputId="tab4.filter_type", label="select an option",
                         choices=c("list accessions by mutation", "list mutations by accession")),
             
             
             fluidRow(
               column(6, 
                      h3("Accession Filter"),
                      textAreaInput(inputId = "tab4.ecoIDs", label = "Ecotype ID(s)",
                                    width = 400, height = 250, value = "" )
                      ),
               column(6,
                      h3("Mutation Filter"),
                      textAreaInput(inputId = "tab4.SNPs", label = "SNP(s)",
                                    width = 400, height = 250, value = "" )
                      )
             ),
             
             actionButton(inputId="tab4.reFilter", label="Submit"),
             
             #tableOutput("tab4Table")
             DT::dataTableOutput("tab4.dataTable")
             
             
             
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

parsebysep <- function (textIn, sep) {
  text <- gsub(" ", "", ttext, fixed = TRUE)
  output <- strsplit(text, sep)
  
  
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


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



server <- function(input, output){

  ##   _________
  ##  /  tab1   \
  ##             --------------------------------------------------
  ## Tab 1 stuff:

  tab1.Genes <- eventReactive(input$STATS_submit,{
    names <- parseInput(input$gene_ids)
    genes <- getGeneInfo(names)
    req(genes != FALSE)
    return(genes)
  })
  
  output$table <- renderTable(tab1.Genes())
    
  SNPStats <- reactive({polymorphTable(tab1.Genes(), strains)})
  
  output$SNPStats_Table <- renderTable(SNPStats(), rownames=TRUE, digits=5)
  
  output$tab1.downloadStats <- downloadHandler(
    filename=function(){
      paste("SNPStats-", Sys.Date(), ".csv", sep="")
    }, 
    content = function(file) {
      write.csv(SNPStats(), file, row.names=FALSE)
    }
  )
  
  output$tab1.downloadGeneInfo <- downloadHandler(
    filename=function(){
      paste("GeneInfo-", Sys.Date(), ".csv", sep="")
    }, 
    content = function(file) {
      write.csv(tab1.Genes(), file, row.names=FALSE)
    }
  )
  

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
  
  output$tab2.downloadSNPData <- downloadHandler(
    filename=function(){
      paste("SNPData-", Sys.Date(), ".csv", sep="")
    }, 
    content = function(file) {
      write.csv(tab2.tableData(), file, row.names=FALSE)
    }
  )
  
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
  
  tab3.EffectValues <- reactive({    
    effects <- c("5_prime_UTR_variant",
                 "intron_variant",
                 "3_prime_UTR_variant",
                 "synonymous_variant",
                 "missense_variant",
                 "upstream_gene_variant",
                 "downstream_gene_variant")
    switch(input$tab3.SNPtype, 
           "All"=effects,
           "Missense"="missense_variant",
           "Coding"= c("missense_variant", "synonymous_variant"))
    #return(input$tab3.filter_value)
  })
  
  output$tab3.debug <- renderPrint({
    tab3.EffectValues()
  })
    
  tab3.labeledSNPs <- reactive({

    data <- tab3.tidyData()
    #keyPOS <- unique(data[which(data$Diversity >= 0.5*max(data$Diversity)), "POS"])
    
    data2 <- data[data$Effect %in% tab3.EffectValues(), ]

    keyPOS <- unique(data2[which(data2$Diversity >= input$tab3.filter_value), "POS"])

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
    
    pal <- c(brewer.pal(8, "Set1"), brewer.pal(7, "Dark2"))
    pallet <- colorFactor(palette=pal, domain=mapdata$SNPs )
    
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

  output$tab3.downloadMapData <- downloadHandler(
    filename=function(){
      paste("MapData-", Sys.Date(), ".csv", sep="")
    }, 
    content = function(file) {
      write.csv(tab3.labeledSNPs(), file, row.names=FALSE)
    }
  )
  
  
  
  ##                                        _________
  ##                                       /  tab4   \ 
  ## --------------------------------------           ----------------
  ## Tab 4 stuff:
  
  tab4.Genes <- eventReactive(input$tab4.Submit, {
    #gene Info for gene on tab 2, updates on 'submit' button press
    names <- parseInput(input$tab4.Gene)
    genes <- getGeneInfo(names[1])
    return(genes)
  })
  
  tab4.tidyData <- reactive({
    tidyVCF <- VCFByTranscript(tab4.Genes()[1, ], strains)
    data <- tidyVCF$dat[tidyVCF$dat$gt_GT != "0|0",]
    # Parse the EFF field
    data <- parseEFF(tidyVCF = data, Transcript_ID = tab4.Genes()$transcript_ID[1])
    
    # calculate diversity
    data <- Nucleotide_diversity(data)
    data <- add_ecotype_details(data)
    data <- subset(data, select=-c(EFF))
    
    
    return(data)
  })
  

  tab4.filteredData <- eventReactive(input$tab4.reFilter, {
    data <- tab4.tidyData
    
    if (input$tab4.ecoIDs != ""){
      
    }
    
    
  })
  
  
  output$tab4.dataTable <- DT::renderDataTable(DT::datatable(tab4.tidyData()))
  
  
}


shinyApp(ui = ui, server = server)

