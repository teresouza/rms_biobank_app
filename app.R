##Source file with paths to data and packages
source("global.R")

server <- function(input, output, session){
  #############################################################
  ###################### DEFINE FUNCTIONS  #################### 
  #############################################################
  ## User defined functions
  ## Boxplot gene expression ##
  exp.boxplot.fun <- function(gene.name){
    df <- rna[rna$Gene.name == gene.name,]
    validate(need(input$geneExp %in% df$Gene.name, "Please input official gene symbol and press 'Submit' to generate plot"))
    ggplot(df, aes(source, log2TPM, fill = type)) + geom_boxplot() + theme_few() +
      scale_fill_npg() + ylab("Log2 TPM") + xlab("Source") + 
      ggtitle(paste(input$geneExp, "expression in RMS tumor and tumoroids")) + 
      theme(axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14),
            axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            strip.text.x = element_text(size = 14)) +
      guides(fill=guide_legend(title="RMS subtype"))
  }
  
  ## Dotplot gene expression ##
  exp.dotplot.fun <- function(gene.name){
    df <- rna[rna$Gene.name == gene.name,]
    validate(need(input$geneExp %in% df$Gene.name, "Please input official gene name and press 'Submit' to generate plot"))
    ggplot(df, aes(x = id2, y = log2TPM, colour = passage, size = 3)) + geom_point() + 
      scale_color_npg() + ylab("Log2 TPM") + xlab("Sample") + 
      ggtitle(paste(input$geneExp, "expression in RMS tumor and tumoroids")) +
      facet_grid(type~source, scales = "free_y", space = "free_y") + coord_flip() + theme_few() +
      theme(panel.grid.major = element_line(colour="lightgrey", size=0.5),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 14),
            axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            strip.text.x = element_text(size = 14)) +
      guides(colour = guide_legend(title = "Passage", override.aes = list(size=10)),
             size=F)
  }
  
  ## Get name of pdf files 
  cna.plot.fun <- function(sample.name){
    pngNames <- list.files("./www", pattern = sample.name)
  }
  
  #############################################################
  ###################### TESTS AND CHECKS #################### 
  #############################################################
  ## Validate inputs and capitalize them if not the case
  observe({
    updateTextInput(session, "geneExp", value = toupper(as.character(input$geneExp)))
  })
  observe({
    updateTextInput(session, "geneMut", value = toupper(as.character(input$geneMut)))
  })
  observe({
    updateTextInput(session, "geneDrug", value = toupper(as.character(input$geneDrug)))
  })
  observe({
    updateTextInput(session, "geneCNATable", value = toupper(as.character(input$geneCNATable)))
  })
  
  #############################################################
  ###################### GENE EXPRESSION ###################### 
  #############################################################
  
  ## Generates boxplot, with groupwise expression for a gene defined by the user
  output$exp.boxplot <- renderPlot({
    #validate(need(input$geneExp %in% data()$Gene.name, "Please input a valid gene name and press 'Submit' to generate plot"))
    input$geneExpButton
    isolate(exp.boxplot.fun(input$geneExp))
  })
  
  ## Generates dotplots for individual samples for a gene defined by the user
  output$exp.sample <- renderPlot({
    input$geneExpButton
    isolate(exp.dotplot.fun(input$geneExp))
  })
  
  #############################################################
  ######################## GENE FUSION ######################## 
  #############################################################
  ## Generates a table with the fusions
  df.fusion <- eventReactive(input$fusionButton, {
    df1 <- fusion %>% dplyr::select(Sample, Source, Passage, Fusion, FFPM,  `Gene 1`, `Gene 2`, `Protein fusion type`) %>%
      add_column(`Left Gene` = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$`Gene 1`), "' target='_blank'>", .$`Gene 1`, "</a>"),
                 `Right Gene` = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$`Gene 2`), "' target='_blank'>", .$`Gene 2`, "</a>")) 
  })

  ## Renders a data table with the fusions specified by the user
  output$fusionTable1 <- DT::renderDataTable({
    df.fusion() %>% .[.$Fusion == input$geneFusion,] %>% 
      dplyr::select(Sample, Source, Passage, Fusion, FFPM,  `Left Gene`, `Right Gene`, `Protein fusion type`)
  }, escape = F)
  ## Renders a data table with all the fusions detected across models and tumors
  output$fusionTable2 = DT::renderDataTable({
    fusion %>% dplyr::select(Sample, Source, Passage, Fusion, FFPM,  `Gene 1`, `Gene 2`, `Protein fusion type`) %>%
      add_column(`Left Gene` = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$`Gene 1`), "' target='_blank'>", .$`Gene 1`, "</a>"),
                 `Right Gene` = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$`Gene 2`), "' target='_blank'>", .$`Gene 2`, "</a>")) %>% 
      dplyr::select(Sample, Source, Passage, Fusion, FFPM, `Left Gene`, `Right Gene`, `Protein fusion type`)
  }, escape = F)
  
   
  #############################################################
  ########################### SNVS ########################### 
  #############################################################
  ## Generates an interactive bubble plot that shows only cancer driver genes (categorized as TSG or oncogenes in COSMIC). A threshold for 
  ## AF is set t 10%.
  output$mutPlot <- renderPlotly({
    dataMutBubble <- snvs[snvs$`Oncogene or TSG` != "" & snvs$AF >= 0.1,] %>%
      arrange(`Sample ID`, Source, Passage, Gene, `Type mutation`, AF)
    dataMutBubble$AF <- as.numeric(dataMutBubble$AF)
    h1 <- ggplot(dataMutBubble, aes(x = interaction(`Sample ID`, Source, Passage, lex.order = TRUE,sep = " "),
                                    y = Gene, fill = `Type mutation`, 
                                    text = gsub("NA", "", paste(paste0("AF: ", AF), `Genomic position`, 
                                                                `AA mutation`, `COSMIC`, sep = "\n")))) +
      geom_point(aes(size = AF)) +
      geom_hline(yintercept = 1:length(dataMutBubble$Gene)+0.5, color = "grey") +
      scale_fill_npg() + theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) + coord_flip() +
      scale_x_discrete(position = 'top') +
      #scale_y_discrete(position = 'right') +
      scale_size_continuous(breaks=c(0.1, 0.25, 0.75, 0.9)) +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.05, face = "bold"),
            #axis.text.x = element_blank(),
            axis.text.y = element_text(face="bold", hjust = 0)) + 
      labs(fill = "SNV type")
    ggplotly(h1,tooltip = "text")
  })
  
  ## Generates mutation tables depending on user input 
  dfMut <- eventReactive(input$mutButton, {
    req(input$queryMut, input$AFvalue)
    mut.df <- snvs %>% select(Sample, Source, Passage, Gene, `Sample ID`, `Genomic position`, AF,`Type mutation`,`COSMIC`, 
                              `AA mutation`, SIFT, PolyPhen2, `Oncogene or TSG`, `CGC Tier`, Legacy)
    if(input$queryMut == "geneMut"){
      validate(need(input$geneMut %in% snvs$Gene, "Please input official gene symbol and press 'Submit' to generate plot"))
      df2 <- mut.df[mut.df$Gene == input$geneMut,] %>% select(-`Sample ID`) %>%
        dplyr::mutate(Gene = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$Gene),"' target='_blank'>", .$Gene,"</a>"),
                      COSMIC= paste0("<a href='", paste0("https://cancer.sanger.ac.uk/cosmic/search?q=", .$Legacy),"' target='_blank'>", .$COSMIC, "</a>")) %>% 
        dplyr::mutate(COSMIC = gsub("NA", "", .$COSMIC)) %>% select(-Legacy)
    } else if(input$queryMut == "sampleMut"){
      df2 <- mut.df[mut.df$`Sample ID` == input$sampleMut,] %>% select(-`Sample ID`)  %>%
        .[.$AF <= input$AFvalue[2] & .$AF >= input$AFvalue[1],] %>% 
        dplyr::mutate(Gene = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$Gene),"' target='_blank'>", .$Gene,"</a>"),
                      COSMIC= paste0("<a href='", paste0("https://cancer.sanger.ac.uk/cosmic/search?q=", .$Legacy),"' target='_blank'>", .$COSMIC, "</a>")) %>% 
        dplyr::mutate(COSMIC = gsub("NA", "", .$COSMIC)) %>% select(-Legacy)
    }
    return(df2)
  })
  output$mutTable <-  DT::renderDataTable({
    dfMut()
  }, escape = F)
  
  
  #############################################################
  ########################### CNAs ########################### 
  #############################################################
  ## Generates a CN table if the user selects "gene" as input
  dfCNATable <- eventReactive(input$cnaButton, {
    req(input$geneCNATable)
    validate(need(input$geneCNATable, "Please input official gene symbol and press 'Submit' to generate table"))
    cna.df <- cna %>% dplyr::select(Sample, Source, Passage, Gene, `Sample ID`, `CR meanLog2`, CR, Call) %>% 
      dplyr::filter (Call != "Neutral") %>%
      filter(Gene == input$geneCNATable) %>% select(-`Sample ID`, -Call) %>%
      dplyr::mutate(Gene = paste0("<a href='", paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",.$Gene),"'target='_blank'>", .$Gene, "</a>"))
  })
  output$CNtable <-  DT::renderDataTable({
    dfCNATable()
  }, escape = F)
  
  
  ## Retrieve png file names that will be pulled from www/ folder
  ## dependiing on user's choice on selectInput drop down menu 
  ## It also binds to cnaButton being pressed
  pngNames <- eventReactive(input$cnaButton, {
    pngNames <- list.files(path = "www")[grepl(input$sampleCN, 
                                               list.files(path = "www",
                                                          pattern = "png", 
                                                          full.names = T))]
    isolate(return(pngNames))
  })
  
  output$CNplotOut <- renderUI({
    input$cnaButton
    validate(need(pngNames(), "Please choose a sample"))
    if(input$sampleCN %in% c("RMS006", "RMS007", "RMS335")){
      tags$div(img(src=pngNames()[1], style="width:75%"), 
               img(src=pngNames()[2], style="width:75%")) 
    } else if(input$sampleCN %in% c("RMS012", "RMS444")){
      tags$div(img(src=pngNames()[1], style="width:75%"), 
               img(src=pngNames()[2], style="width:75%"), 
               img(src=pngNames()[3], style="width:75%"))
    } else {
    isolate(img(src=pngNames(), style="width:75%"))
    }
    # isolate(tags$iframe(#style = 'height: 680px; width: 960px;',
    #   src = paste0(input$sampleCN, "*.pdf")))
    # tags$iframe(style="height:50%; width:50%; scrolling=no", 
    #                     src=pngNames)
    #list(src = pngNames, contentType = "image/png")
    })
  
  

  
  
  #############################################################
  ####################### DRUG SCREENING ##################### 
  #############################################################
  ##Generates the drug info data tables 
  output$drugInfo = DT::renderDataTable({
    drugs.cor[["class"]][,c(2,5:11)]
  })
  
  ## Generates the input for correlation plots based on user choice ("drug" or "drug class")
  corPlotInput <- eventReactive(input$drugButton, {
    if(input$corType == "drug"){
      gene.df <- rna[rna$Gene.name == input$geneDrug & rna$source == "Tumoroid" & rna$passage == "Standard",] %>%
        select(log2TPM, id2) %>% data.table
      drug.df <- drugs.cor[["AUC"]][drugs.cor[["AUC"]]$Drug == input$drug,]
      df4 <- drug.df[gene.df, on=c(sample ="id2"), nomatch=NULL]
    } else if(input$corType == "drugClass"){
      df4 <- drugs.cor[["cor"]][drugs.cor[["cor"]]$gene == input$geneDrug & drugs.cor[["cor"]]$TargetType == input$drugClass,] %>%
        mutate(col = ifelse(r > 0.5, "darkred", ifelse(r < -0.5, "darkgreen", "gray")))
    }
    isolate(return(df4))
  })
  
  ## Generates a linear correlation plot (AUC versus log2TPM expression in tumoroids)
  output$corPlot <- renderPlot({
    validate(need(corPlotInput()$log2TPM, "Please choose a drug"))
    input$drugButton
    isolate(ggplot(corPlotInput(), aes(x = log2TPM, y = AUC)) + 
              ylab("AUC") + xlab("log2 TPM") + 
              ggtitle(paste0(input$drug, " correlation to ", input$geneDrug, " expression in RMS tumoroids (r =", 
                             round(cor(corPlotInput()$AUC, corPlotInput()$log2TPM),digits = 2), ")")) +
              scale_color_npg() + 
              geom_smooth(method = "lm", se = F, color = "gray") + theme_classic() + 
              geom_text(aes(label = sample)) + 
              theme(axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=14),
                    strip.text.x = element_text(size = 14)))
  })
  
  ## Generates a dotplot with the correlations for a specific drug class
  output$corPlot2 <- renderPlot({
    validate(need(corPlotInput()$r, "Please choose a drug class"))
    input$drugButton
    isolate(ggplot(corPlotInput(), aes(x =r, y = reorder(drug, -r), color = col, size = 3)) + 
              geom_point() +
              scale_color_identity() +
              ylab("Drug") + xlab("Pearson's r") + 
              ggtitle(paste(input$drugClass, "drug class correlation to", input$geneDrug, "expression in RMS tumoroids")) +
              geom_vline(xintercept = c(-0.5, 0.5), linetype="dashed", 
                         color = "lightgray") + xlim(c(-1, 1)) +
              theme_classic() +
              theme(axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size = 14),
                    axis.text=element_text(size=14),
                    axis.title=element_text(size=14),
                    strip.text.x = element_text(size = 14)) + 
              guides(size=F, fill =F))
  })
  
  #############################################################
  ####################### METADATA ########################### 
  #############################################################
  ## Renders data tables containing the pseudonymized metadata
  output$metadata = DT::renderDataTable({
    metadata
  })
  
  #############################################################
  #################### ABBREVIATIONS ########################## 
  #############################################################
  ## Renders a staticc table, alphabetically ordered, containing the definitions of abbreviations used throughout the app
  output$abbreviations = renderTable({
    abvs
  })
}


ui <- navbarPage("PMC RMS biobank",
                 tabPanel("About",
                          tags$head(tags$style(
                            HTML('#sidebar {background-color: white;}'))),
                          fluidPage(
                            titlePanel("About the PMC RMS biobank app"),
                            fluidRow(
                              column(3, br(), align = "right",
                                     br(),
                                     img(src = "pm_logo2.png", height = 150, width = 300)),
                              br(),
                              br(),
                              column(8, offset =1, align = "right",
                                     p("This app was developed at the Holstege group at the Prinses Máxima Center for Pediatric Oncology."),
                                     p("This app contains molecular data from rhabdomyosarcoma (RMS) samples, including ", 
                                       strong("gene expression"), " and ", strong("predicted gene fusion "), "from bulk mRNA sequencing data,", 
                                       strong(" copy number alterations")," against a panel of normals (PoN) and ", 
                                       strong("somatic mutations"), " against paired germline samples."),
                                     p("For more information about data collection and analysis, please refer to the original paper.")
                              )
                            ))),
                 tabPanel("Gene Expression",
                          fluidRow(
                            headerPanel(h3("Gene Expression")),
                            column(12, h6("Check the expression (log2 transcripts per million, TPM) of genes in RMS samples.")),
                            column(3,
                                   textInput("geneExp", h4("Gene symbol")),
                                   actionButton("geneExpButton", "Submit", icon("play"), width = "99%"), 
                                   h4("Guide"),
                                   p(strong("eRMS"),": Embryonal RMS"),
                                   p(strong("P3F aRMS"), ": PAX3-FOXO1 alveolar RMS"),
                                   p(strong("P7F aRMS"), ": PAX7-FOXO1 alveolar RMS"),
                                   p(strong("FN aRMS"),": Fusion negative alveolar RMS"),
                                   p(strong("Standard"),": for tumoroids, the standrad model generated from the primary tumor"),
                                   p(strong("Second"),": for tumoroids, a secondary model generated from the primary tumor"),
                                   p(strong("Late"),": for tumoroids, the standard model at a later passage")
                                   ),
                            mainPanel(  
                              tabsetPanel(type = "tabs",
                                          tabPanel("Expression per subtype", plotOutput("exp.boxplot")),
                                          tabPanel("Expression per sample", plotOutput("exp.sample", height = 750))
                              )
                            ),
                            tags$head(
                              tags$style(HTML('#geneExpButton{background-color:orange}')))
                          )
                 ),
                 tabPanel("Gene fusion",
                          fluidRow(
                            headerPanel(h3("Gene fusion")),
                            column(12, h6("Check fusions detected in the RNA sequencing data in RMS samples.")),
                            column(3,
                                   selectInput("geneFusion", h4("Select a fusion"),
                                               unique(fusion$Fusion), selected = NULL),
                                   actionButton("fusionButton", "Submit", icon("play"), width = "99%")),
                            column(9,  
                                   tabsetPanel(type = "tabs",
                                               tabPanel("Query result", DT::dataTableOutput("fusionTable1")),
                                               tabPanel("All samples", DT::dataTableOutput("fusionTable2"))
                                   )
                            ),
                            tags$head(
                              tags$style(HTML('#fusionButton{background-color:orange}')))
                          )
                 ),
                 navbarMenu("Genomic alterations",
                            tabPanel("Somatic mutations",
                                     fluidRow(
                                       headerPanel(h3("Somatic mutations")),
                                       column(12, h6("Check non-synonymous somatic mutations in cancer genes (1,914) detected by whole-genome sequencing in RMS samples.")),
                                       column(2,
                                              radioButtons("queryMut", h4("Search by"),
                                                           c("Gene" = "geneMut","Sample" = "sampleMut"))),
                                       column(3,conditionalPanel(condition = "input.queryMut == 'geneMut'",
                                                                 textInput("geneMut", h4("Gene symbol"))),
                                              conditionalPanel(condition = "input.queryMut == 'sampleMut'",
                                                               selectInput("sampleMut", h4("Select a sample"),
                                                                           menu.data$models), selected = NULL)),
                                       column(4,conditionalPanel(condition = "input.queryMut == 'sampleMut'",
                                                                 sliderInput("AFvalue", h4("Allele frequency (AF)"), 
                                                                             min = 0, max =1, 
                                                                             value = c(0.4,1)))),
                                       column(3 ,conditionalPanel(condition = "input.queryMut == 'sampleMut'",
                                                                        h4("Use the slider to specify the ranges for variant allele frequency (AF)")))),
                                     fluidRow(column(12,
                                                     actionButton("mutButton", "Submit", icon("play"), width = "99%"))),
                                     tags$head(
                                       tags$style(HTML('#mutButton{background-color:orange}'))
                                     ),
                                     br(),
                                     fluidRow(
                                       column(12, align = 'center', tabsetPanel(type = "tabs",
                                                                                tabPanel("Query result",
                                                                                         DT::dataTableOutput("mutTable")),
                                                                                tabPanel("Cancer driver genes", 
                                                                                         plotlyOutput("mutPlot", height = 750))
                                       )))
                            ),
                            tabPanel("Copy number alterations",
                                     fluidRow(
                                       headerPanel(h3("Copy number alterations (CNAs)")),
                                       column(12, h6("Check copy number alterations in cancer genes (1,1914) detected by whole-genome sequencing in RMS samples.", tags$br(),
                                                     "When selecting genes only significant alterations (called as gains or losses 
                                 compared to a Panel of Normals [PoN]) will be displayed.", tags$br(),
                                 "By selecting samples, circos plots containing genome-wide alterations from tumors (outer circle) and tumoroids (inner circle) will be displayed.", tags$br(),
                                 "No cutoffs for gains or losses are applied in the circos plots; the y-axis is represented in log2 scale.")),
                                       column(3, 
                                              radioButtons("queryCN", h4("Query CNAs"),
                                                           c("Search gene" = "geneCN","Generate CN plots" = "CNplot"))),
                                       column(3,
                                              conditionalPanel(condition = "input.queryCN == 'geneCN'",
                                                               textInput("geneCNATable", h4("Gene symbol"))),
                                              conditionalPanel(condition = "input.queryCN == 'CNplot'",
                                                               selectInput("sampleCN", h4("Select a sample"),
                                                                           menu.data$models, selected = NULL)))
                                     ),
                                     fluidRow(column(12,
                                                     actionButton("cnaButton", "Submit", icon("play"), width = "99%"))),
                                     tags$head(
                                       tags$style(HTML('#cnaButton{background-color:orange}'))
                                     ),
                                     br(),
                                     tabsetPanel(type = "tabs", 
                                                 tabPanel(align = 'center', "Query result",
                                                          conditionalPanel(condition = "input.queryCN == 'geneCN'",
                                                                           DT::dataTableOutput("CNtable")),
                                                          conditionalPanel(condition = "input.queryCN == 'CNplot'",
                                                                           uiOutput("CNplotOut")))
                                     )
                            )
                 ),
                 tabPanel("Drug screening",
                          fluidRow(
                            headerPanel(h3("Gene expression and drug response")),
                            column(12, h6("Find correlations between drug area under the curve (AUC) and baseline gene expression (log2TPM) in RMS tumoroid samples. Only protein coding genes were used to perform correlations.")),
                            column(3,
                                   textInput("geneDrug", h4("Gene symbol"))),
                            column(3,
                                   radioButtons("corType", h4("Correlate expression to"),
                                                c("Drug" = "drug","Drug class" = "drugClass"))),
                            column(3,
                                   conditionalPanel(condition = "input.corType == 'drug'",
                                                    selectInput("drug", h4("Select a drug"),
                                                                menu.data$drug.names)),
                                   conditionalPanel(condition = "input.corType == 'drugClass'",
                                                    selectInput("drugClass", h4("Select a drug class"),
                                                                menu.data$drug.classes)))
                          ), 
                          fluidRow(column(12,
                                          actionButton("drugButton", "Submit", icon("play"), width = "99%"))),
                          tags$head(
                            tags$style(HTML('#drugButton{background-color:orange}'))
                          ),
                          br(),
                          fluidRow(column(12, align = 'center', 
                                          tabsetPanel(type = "tabs", 
                                                      tabPanel("Query result", 
                                                               conditionalPanel(condition = "input.corType == 'drug'",
                                                                                plotOutput("corPlot", width = "75%")),
                                                               conditionalPanel(condition = "input.corType == 'drugClass'",
                                                                                plotOutput("corPlot2", height = 750, width ="75%") %>% 
                                                                                  withSpinner(color = "orange", size = 2, type = 6))),
                                                      tabPanel("Drug info", DT::dataTableOutput("drugInfo")))
                          )
                          )
                 ),
                 # tabPanel("Generate report",
                 #          fluidRow(
                 #            headerPanel(h3("Generate a a molecular report")),
                 #            column(12, h6("UNDER CONSTRUCTION")),
                 #          )
                 # ),
                 navbarMenu("More",
                            tabPanel("List of abbreviations",
                                     fluidRow(
                                       column(10, align = "center", tabPanel("List of abbreviations and features", tableOutput("abbreviations"))
                                     )
                                     )
                                     ),
                            tabPanel("Metadata",
                                     fluidRow(tabPanel("Metadata", DT::dataTableOutput("metadata")
                                     )
                                     )
                            ),
                            tabPanel("Contact",
                                     fluidRow(
                                       column(12, h4("Our team"), align = "center"),
                                       br(),
                                       column(12, align = "center",
                                              p(strong("Michael Meister (Coordinator of the Soft Tissue Sarcoma division)"), ": M.T.Meister@prinsesmaximacentrum.nl"),
                                              p(strong("Marian Groot Koerkamp (Head technician and Data steward)"), ": M.GrootKoerkamp@prinsesmaximacentrum.nl"),
                                              p(strong("Terezinha Souza (Data analyst and App developer)"), ": T.M.deSouza@prinsesmaximacentrum.nl"),
                                              p(strong("Frank Holstege (Principal Investigator)"), ": F.C.P.Holstege@prinsesmaximacentrum.nl"))
                                     )
                            )
                 )
)

# Run the application 
shinyApp(ui = ui, server = server)