library(shiny)
library(biomaRt)
library(leaflet)
library(RColorBrewer)

CSSCode <- tags$head(tags$style(
   HTML("
      .checkbox-format {
         -webkit-column-width: 350px;
         -moz-column-width: 350px;
         column-width: 350px;
      }
                                    
      .input-format {
         background-color: #dddddd;
         border: 1px solid #dddddd;
         border-radius: 12px;
         padding:1px 15px 10px 10px;
      }

      .output-format {
         border: 5px solid #dddddd;
         border-radius: 12px;
         padding:1px 15px 15px 20px;
      }

      .btn-default{ 
         color: #333;
         background-color: #eeeeee;
         border-color: #ccc;
      }  

      .form-control{ 
         color: #333;
         background-color: #eeeeee;
         border-color: #ccc;
      }  

      h1 {
         font-family: Helvetica;
         font-weight: 500;
         line-height: 1.1;
         color: #48ca3b;
         background-color: #dce4f2;
         border: 10px solid #dce4f2;
         border-radius: 12px;
      }
                                     
   ")
                                  
))


ui <- fluidPage(
  CSSCode,
  headerPanel("Arabidopsis Natural Variation Webtool"),
  "This app will provide an interface to examine the natural variation of specified genes of interest in the 1001 Genomes project dataset",
  tags$h5('style'="color:red", "this app is a work in progress"),
  tabsetPanel(
    tabPanel("SNP Stats",
        ## Tab 1 ###############################################################
      tags$br(),
      tags$div(class="input-format",
               tags$h3("Gene Selection"),
               tags$h5("Type a list of gene loci in the box below, separated by commas "),
               textAreaInput(inputId = "gene_ids", label = NULL,
                             width = 600, height = 75, value = "AT3G62980, AT4G03190, AT3G26810, AT1G12820, AT4G24390, AT5G49980, AT2G39940" ),
               actionButton(inputId="STATS_submit", label = "Submit"),
               tags$br()
      ),
       #tags$hr(),
       #uiOutput("tab1.gene_table_ui"),
      tags$hr(),
      
      tags$div(class="output-format",
               tags$h3("Selected Gene Information"),
               tags$h5("this table provides details on the gene(s) input above, including transcript IDs, and chromosome position information on the start and end of the transcript"),
               downloadButton("tab1.downloadGeneInfo","Download Content of Table Below"),
               DT::dataTableOutput("tab1.genes_table")
        
      ),
      tags$br(),
      tags$div(class="output-format",
      tags$h3("SNP Type and Diversity Statistics"),
          HTML("<h5>
                   This table provides basic statistics on the polymorphisms present in the given gene.
                   <br/>the columns \"5_prime_UTR_variant\" through \"coding_total\" are total, non unique numbers of variants, \"coding_total\" is the sum of missense and synonymous variants.
                   <br/>the final four columns are Nucleotide Diversity values, for different sections and SNP types
                   <br/> <strong>NOTE:</strong> download button downloads content of both tables to a single file.
                </h5>"),
          downloadButton("tab1.downloadStats","Download Content of Tables Below"),
          #tableOutput("SNPStats_Table")
          tags$h4("SNP Counts"),
          DT::dataTableOutput("tab1.SNPcounts"),
          tags$h4("Nucleotide Diversity Statistic"),
          DT::dataTableOutput("tab1.Diversity_table")
      )
    ),

    tabPanel("Diversity Plot",
        ## Tab 2 ###############################################################
      tags$br(),
      tags$div(class="input-format", 
               tags$h3("Gene Select"),
               tags$h5("Select a transcript ID in the box below"),
               uiOutput("tab2.selectGene")
               # textInput(inputId = "plotGene", label =NULL,
               #           value = "AT1G80490"),
               #actionButton(inputId="tab2.Submit", label = "Submit")
      ),
      
      tags$hr(),
      tags$div(class="output-format", 
          tags$h3("Selected Gene Information"),
          tableOutput("tab2.GeneInfo")
      ),
      tags$br(),
      
      tags$div(class="output-format", 
          tags$h3("Plot of Nucleotide Diversity Statistic by Codon"),
          tags$h5("click and drag a box accross the plot below to see details on specific points"),
          plotOutput("diversityPlot", brush="plot_brush", click="plot_click", height = 400),
          verbatimTextOutput("info")
      ),
      tags$br(),    
      tags$div(class="output-format", 

          tags$h3("Diversity Plot Data"),
          tags$h5("This table provides the raw data from the plot. \"POS\" is the chromosomal position of the SNP, 
                  the Amino_Acid_Change field provides both the amino acid as well as the base change"),
          downloadButton("tab2.downloadSNPData","Download Content of Table Below"),
          DT::dataTableOutput("Diversity_Table")
      )
    ),
    
    tabPanel("SNP Mapping",
             ## Tab 3 ##########################################################
             tags$br(),
             tags$div(class="input-format", 
                 tags$h3("Gene Select and Filter Parameters"),
                 tags$h5("Type a single gene locus in the box below"),
                 # textInput(inputId="tab3.Gene", label=NULL,
                 #           value="AT1G80490"),
                 uiOutput("tab3.selectGene"),
                 # actionButton(inputId="tab3.Submit", label="Submit"),
                 sliderInput(inputId="tab3.filter_value", label="Log Nucleotide diversity filter limit",
                             min=-4, max=0, value=-2, step=0.05),
                 radioButtons("tab3.SNPtype", "Type of SNP to mark", 
                              choices=c("All", "Coding", "Missense")),
                 verbatimTextOutput("tab3.debug")
             ),
             
             tags$br(),

             uiOutput("tab3.mutation_checkbox"), 

             tags$hr(),
             tags$div(class="output-format",
             tags$h3("Accession Map"),
             tags$h5("Zoom with scroll wheel, click and drag to pan. click on individual point to see details. 
                     use the layers pop out to the lower left of the map to hide or show sets of markers with similar mutations"),
             leafletOutput("tab3.map", height="650")
             ),
             tags$br(),
             tags$div(class="output-format",
                 tags$h3("Map Data"),
                 downloadButton("tab3.downloadMapData","Download Content of Table Below"),
                 DT::dataTableOutput("tab3.dataTable")
             )
    ),
    
    tabPanel("About",
             ## Tab 4 ##########################################################
             tags$br(),
             column(6, 
                    tags$div(class="output-format",
                             includeHTML("Glossary.html")
                    )
             ),
             column(6, 
                    tags$div(class="output-format",
                             includeMarkdown("Bibliography.Rmd")
                    )                    
             )
             
             
    )
    
    # tabPanel("Accessions and Mutations",
    #          ## Tab 4
    #          textInput(inputId="tab4.Gene", label="Type list of gene ID's",
    #                    value="AT1G80490"),
    #          actionButton(inputId="tab4.Submit", label="Get Data"),
    #          
    #          selectInput(inputId="tab4.filter_type", label="select an option",
    #                      choices=c("list accessions by mutation", "list mutations by accession")),
    #          
    #          uiOutput("tab4.ui2"),
    #          
    #          fluidRow(
    #            column(6, 
    #                   h3("Accession Filter"),
    #                   textAreaInput(inputId = "tab4.ecoIDs", label = "Ecotype ID(s)",
    #                                 width = 400, height = 250, value = "" )
    #                   ),
    #            column(6,
    #                   h3("Mutation Filter"),
    #                   textAreaInput(inputId = "tab4.SNPs", label = "SNP(s)",
    #                                 width = 400, height = 250, value = "" )
    #                   )
    #          ),
    #          
    #          actionButton(inputId="tab4.reFilter", label="Submit"),
    #          
    #          #tableOutput("tab4Table")
    #          DT::dataTableOutput("tab4.dataTable")
    #)
    
    
  )

  # "THIS IS THE FOOTER"
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
  tab2data <- parseEFF(tab2data)
  tab2data <- Nucleotide_diversity(tab2data)
  
  coding_variants <- coding_Diversity_Plot(tab2data)
  
  return(coding_variants)
}


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

  all.Genes <- eventReactive(input$STATS_submit,{
    # list of genes for tab 1, updated on pressing submit button
    names <- parseInput(input$gene_ids)
    genes <- getGeneInfo(names)
    req(genes != FALSE)
    return(genes)
  })
  
  output$tab1.genes_table <- DT::renderDataTable(all.Genes()[, -c(5,6)], options=list(paging=FALSE, searching=FALSE))
    
  #SNPStats <- reactive({polymorphTable(tab1.Genes(), strains)})
  
  all.VCFList <- reactive({
    output <- VCFList(all.Genes())
    output <- llply(output, parseEFF)
    output <- llply(output, Nucleotide_diversity)
    return(output)
  })
  
  SNPStats <- reactive({ ldply(all.VCFList(), polymorphRow, geneInfo=all.Genes(), .id="transcript_ID") })
  
  output$SNPStats_Table <- renderTable(SNPStats())
  
  output$tab1.SNPcounts <- DT::renderDataTable(SNPStats()[,1:8], options=list(paging=FALSE, searching=FALSE))
  output$tab1.Diversity_table <- DT::renderDataTable(SNPStats()[, c(1,9:13)], options=list(paging=FALSE, searching=FALSE))
  
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
      write.csv(all.Genes(), file, row.names=FALSE)
    }
  )
  

    output$tab1.gene_table_ui <- renderUI({
      if (input$STATS_submit==0){
        return()
      }
      
          
      tagList(
        tags$div(class="input-format",
                 tags$h3("Mutation select"),
                 tags$h5("Select the SNPs you want to see on the map by clicking the checkboxes"),
                 tags$div(class="checkbox-format", 
                          checkboxGroupInput("tab3.mutation_select", "select_mutations to display", choices=tab3.mutationList())
                 ),
                 
                 actionButton(inputId="tab3.update_map", label = "Update Map")
        )
      )
      
      # tags$div(class="output-format",
                # tags$h3("Selected Gene Information"),
                # tags$h5("this table provides details on the gene(s) input above, including transcript IDs, and chromosome position information on the start and end of the transcript"),
                # downloadButton("tab1.downloadGeneInfo","Download Content of Table Below")
                # DT::dataTableOutput("tab1.genes_table")
      # )
    })
 
  ##                 _________
  ##                /  tab2   \
  ## ---------------           -------------------------------------
  ## Tab 2 stuff:

  output$tab2.selectGene <- renderUI({
    tagList(
      selectInput("tab2.transcript_ID", label=NULL, choices=all.Genes()$transcript_ID),
      actionButton(inputId="tab2.Submit", label = "Submit")
    )
  })
    
  tab2.Genes <- eventReactive(input$tab2.Submit, {
      #gene Info for gene on tab 2, updates on 'submit' button press
    # names <- parseInput(input$plotGene)
    # genes <- getGeneInfo(names[1])
    # return(genes)
    return(all.Genes()[ all.Genes()$transcript_ID == input$tab2.transcript_ID,])
  })

  output$tab2.GeneInfo <- renderTable(tab2.Genes())
    #rendered table of Gene info

  #tab2.tableData <- reactive({load_tab_2_Data(tab2.Genes())})
    #SNP reactive data
  tab2.tableData <- eventReactive(input$tab2.Submit, {
    tab2data <- all.VCFList()[[input$tab2.transcript_ID]]
    coding_variants <- coding_Diversity_Plot(tab2data)
    return(coding_variants)
  })
  

  output$Diversity_Table <- DT::renderDataTable(tab2.tableData())
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
  
  output$tab3.selectGene <- renderUI({
    tagList(
      selectInput("tab3.transcript_ID", label=NULL, choices=all.Genes()$transcript_ID),
      actionButton(inputId="tab3.Submit", label = "Submit")
    )
  })
  
  tab3.Genes <- eventReactive(input$tab3.Submit, {
    #gene Info for gene on tab 3, updates on 'submit' button press
    return(all.Genes()[ all.Genes()$transcript_ID == input$tab3.transcript_ID,])
  })
  
  
  tab3.tidyData <- eventReactive(input$tab3.Submit, {
    
    # Get the data
    # tidyVCF <- VCFByTranscript(tab3.Genes()[1, ], strains)
    # data <- tidyVCF$dat
    # # Parse the EFF field
    # data <- parseEFF(tidyVCF = data)
    # 
    # # calculate diversity
    # data <- Nucleotide_diversity(data)
    
    data <- all.VCFList()[[input$tab3.transcript_ID]]
    
    # remove 0|0 genotypes
    data <- data[data$gt_GT != "0|0",]
    
    return(data)
  })
  
  tab3.EffectValues <- reactive({    
    # effects <- c("5_prime_UTR_variant",
    #              "intron_variant",
    #              "3_prime_UTR_variant",
    #              "synonymous_variant",
    #              "missense_variant",
    #              "upstream_gene_variant",
    #              "downstream_gene_variant")
    
    effects <- unique(tab3.tidyData()$Effect)

    return( switch(input$tab3.SNPtype, 
            "All"=effects,
            "Missense"="missense_variant",
            "Coding"= c("missense_variant", "synonymous_variant"))
    )
  })
  
  output$tab3.debug <- renderPrint({
    # temporary debug output
    input$tab3.mutation_select
    input$tab3.Submit
  })
  
  tab3.filteredByDiv <- reactive({
    # filter by diversity slider and SNP type radio button then add SNPs column
    
    data <- tab3.tidyData()
    
    # filter by effect type (all, coding, or missense)
    data2 <- data[data$Effect %in% tab3.EffectValues(), ]
    
    # filter on positions with diversity greater than or equal to the 10^slider value
    keyPOS <- unique(data2[which(data2$Diversity >= 10**input$tab3.filter_value), "POS"])
    keydata <- data[data$POS %in% keyPOS, ]
    
    return(keydata)
  })
  
  tab3.mutationList <- reactive({
    mutList <- label_bySNPs(tab3.filteredByDiv(), collapse=FALSE)$SNPs
    mutList <- unique(mutList[!is.na(mutList)])
    return(mutList)
  })
  
  output$tab3.mutation_checkbox <- renderUI({
    tagList(
      tags$div(class="input-format",
          tags$h3("Mutation select"),
          tags$h5("Select the SNPs you want to see on the map by clicking the checkboxes"),
          tags$div(class="checkbox-format", 
                   checkboxGroupInput("tab3.mutation_select", "select_mutations to display", choices=tab3.mutationList())
          ),
          actionButton(inputId="tab3.update_map", label = "Update Map")
      )
    )
    
  })
    
  tab3.labeled <- eventReactive(input$tab3.update_map, {
    # a dataframe with a single row per accession, containing accession info, 
    
    # start with the data filtered by the diversity slider and type buttons
    data <- tab3.filteredByDiv()
    # label by SNPs creates column SNPs with text strings formatted [transcriptID|AA_Change]
    data <- label_bySNPs(data, collapse=FALSE)
    # filter on selected SNPs
    data <- data[data$SNPs %in% input$tab3.mutation_select, ]
    # combine mutations to single row (this is slow)
    data <- ddply(data, "Indiv", summarise, SNPs=paste(SNPs, collapse=","))
    # add back ecotype details
    data <- add_ecotype_details(data)
    
    return(data)
  })
  
  output$tab3.map <- renderLeaflet({
    
    mapdata <- tab3.labeled()
    
    # Reorganize to plot NA's underneath non NA's
    mapdata <- rbind(mapdata[is.na(mapdata$SNPs), ], mapdata[!is.na(mapdata$SNPs), ])
    
    # make a field with text to be displayed when clicking on a marker
    mapdata$popup <- paste("EcoID:",  mapdata$Indiv,"Name:", mapdata$Name, " SNPs:", mapdata$SNPs)
    
    # create the color pallet for the map points
    pal <- brewer.pal(8, "Set1")
    pallet <- colorFactor(palette=pal, domain=mapdata$SNPs)
    
    # create a new leaflet map
    map <- leaflet()
    map <- addProviderTiles(map, providers$Stamen.TonerLite,
                     options = providerTileOptions(noWrap = TRUE))
    
    # groupnames to be used by draw groups of points as separate layers below
    groupnames <- unique(mapdata$SNPs)
    groupnames <- groupnames[!is.na(groupnames)]
    
    # add markers for NA points first so they are furthest back layer
    map <- addCircleMarkers(map, data=mapdata[is.na(mapdata$SNPs), ], color= "#9b9b9b", group="NA", 
                            radius=6, popup= ~popup, stroke=FALSE, fillOpacity=0.6)
   
    # for each of the group names, add a set of markers  
    for (SNP in groupnames){
          map <- addCircleMarkers(map, data=mapdata[mapdata$SNPs == SNP, ], color= ~pallet(SNPs), group= SNP, 
                            radius=6, popup= ~popup, stroke=FALSE, fillOpacity=0.85)
    }

    # add the legend to the map
    map <- addLegend(map, position="bottomright", pal=pallet, 
                     values=mapdata$SNPs, title="Marker Colors", opacity=1)
    
    # add layer control to map to turn on or off groups of points
    map <- addLayersControl(map, overlayGroups=c(groupnames, "NA"), 
                            options = layersControlOptions(collapsed = TRUE),
                            position="bottomleft")
    
    return(map)
  })
  
  output$tab3.dataTable <- DT::renderDataTable(tab3.labeled())

  output$tab3.downloadMapData <- downloadHandler(
    filename=function(){
      paste("MapData-", Sys.Date(), ".csv", sep="")
    }, 
    content = function(file) {
      write.csv(tab3.labeled(), file, row.names=FALSE)
    }
  )
  
  
  
  ##                                        _________
  ##                                       /  tab4   \ 
  ## --------------------------------------           ----------------
  ## Tab 4 stuff:
  
  # tab4.Genes <- eventReactive(input$tab4.Submit, {
  #   #gene Info for gene on tab 2, updates on 'submit' button press
  #   names <- parseInput(input$tab4.Gene)
  #   genes <- getGeneInfo(names[1])
  #   return(genes)
  # })
  # 
  # tab4.tidyData <- reactive({
  #   tidyVCF <- VCFByTranscript(tab4.Genes()[1, ], strains)
  #   data <- tidyVCF$dat[tidyVCF$dat$gt_GT != "0|0",]
  #   # Parse the EFF field
  #   data <- parseEFF(tidyVCF = data, Transcript_ID = tab4.Genes()$transcript_ID[1])
  #   
  #   # calculate diversity
  #   data <- Nucleotide_diversity(data)
  #   data <- add_ecotype_details(data)
  #   data <- subset(data, select=-c(EFF))
  #   
  #   
  #   return(data)
  # })
  # 
  # output$tab4.ui2 <- renderUI({
  #   
  #   if (input$tab4.submit==0) {
  #     return(tags$h4("mot submitted"))
  #   }
  #   
  #   return(tags$h4("submitted"))
  #   
  #   
  #   
  # })
  # 
  # tab4.filteredData <- eventReactive(input$tab4.reFilter, {
  #   data <- tab4.tidyData
  #   
  #   if (input$tab4.ecoIDs != ""){
  #   }
  #   
  #   
  # })
  # 
  # output$tab4.dataTable <- DT::renderDataTable(DT::datatable(tab4.tidyData()))
  
  
}


shinyApp(ui = ui, server = server)

