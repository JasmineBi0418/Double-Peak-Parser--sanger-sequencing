library(shiny)
library(bslib)
library(sangerseqR)
library(Biostrings)
library(ggplotify)
library(plotly)
library(data.table)
library(stringr)
library(stringi)
library(DT)
library(shinyBS)
library(DECIPHER)
library(BiocManager)
library(ggplot2)
library(ggbreak)
library(shinyWidgets)
library(shinyalert)
library(tidyr)
library(msa)
library(shinyjs)
library(dplyr)
library(shinydashboard)
library(shinyWidgets)

options(repos = BiocManager::repositories())


#resources for genetic code
select_choices<-c("Standard"="1",
                  "Vertebrate Mitochondrial"="2",
                  "Yeast Mitochondrial"="3",
                  "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"="4",
                  "Invertebrate Mitochondrial"="5",
                  "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"="6",
                  "Echinoderm Mitochondrial; Flatworm Mitochondrial"="9",
                  "Euplotid Nuclear"="10",
                  "Bacterial, Archaeal and Plant Plastid"="11",
                  "Alternative Yeast Nuclear"="12",
                  "Ascidian Mitochondrial"="13",
                  "Alternative Flatworm Mitochondrial"="14",
                  "Blepharisma Macronuclear"="15",
                  "Chlorophycean Mitochondrial"="16",
                  "Trematode Mitochondrial"="21",
                  "Scenedesmus obliquus Mitochondrial"="22",
                  "Thraustochytrium Mitochondrial"="23",
                  "Pterobranchia Mitochondrial"="24",
                  "Candidate Division SR1 and Gracilibacteria"="25",
                  "Pachysolen tannophilus Nuclear"="26")
# Define UI ----
ui <- fluidPage(theme = bs_theme(version = 4, bootswatch = "flatly"),
                tags$head(
                  tags$style(HTML("
            .footer {
                position: fixed;
                right: 10px;
                bottom: 10px;
                color: #555;
            }
        "))
                ),
                useShinyjs(),  # Initialize shinyjs
                tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
                
                hidden(div(id = "loader", class = "loader")),  # Loader div
                navbarPage("Double Peak Parser",
                           tabPanel("How to Use with Example Data",
                                    fluidPage(
                                      h4("Purpose of this app:"),
                                      p("This Shiny app is designed for analyzing sanger sequencing double peak data and determining the most likely wild type. It provides insights into the analysis results."),
                                      br(),
                                      br(),
                                      h4("Using the App with Example Data"),
                                      p("Follow these steps to navigate and utilize the Double Peak Parser app with the provided example data:"),
                                      tags$head(
                                        # Include custom CSS to style the paragraph
                                        tags$style(HTML(".highlight {
                                        color: maroon; # Change color as needed
                                        font-weight: bold;}"))
                                      ),tags$head(
                                        tags$link(rel = "stylesheet", 
                                                  type = "text/css", 
                                                  href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.3/css/all.min.css")
                                      ),tags$head(
                                        tags$style(HTML("
                                        .help-tip {
                                        position: relative;
                                        cursor: pointer;
                                        }
                                        .help-tip .fa-question-circle {
                                        color: #029be5;
                                        }
                                        .help-tip:hover .tooltiptext {
                                        visibility: visible;}
                                        .help-tip .tooltiptext {
                                        visibility: hidden;
                                        width: 140px;
                                        background-color: black;
                                        color: #fff;
                                        text-align: center;
                                        border-radius: 6px;
                                        padding: 5px 0;
                                        position: absolute;
                                        z-index: 1;
                                        bottom: 150%;
                                        left: 50%;
                                        margin-left: -75px;
                                        }
                                                        "))
                                      )
                                      ,
                                      tags$ol(
                                        tags$li("Toggle the 'Load Example Data' switch on the 'Chromatogram & Quality plot' tab to load the pre-loaded example ABI file."),
                                        tags$li("Explore different sections of the app to see how the example data is processed and visualized."),
                                        tags$li("Utilize the 'Chromatogram & Quality plot' tab to examine the sequencing data's quality and details."),
                                        tags$li("Analyze the 'Parsing Results' and 'Amino Acid Translate' tabs to interpret the analysis outcomes."),
                                        tags$li("Find detailed information about genetic codes and relevant data in the 'References' tab."),
                                      ),
                                      
                                      p("This example demonstrates the app's capabilities and features. You can apply similar steps to your data by uploading it in the respective sections."),
                                      
                                      br(),
                                      
                                      
                                      h5("Note:"),
                                      p("The example data is preloaded and doesn't require any file uploads; simply activate the switch to view the data."),
                                      p("The example data will be automatically unloaded when you upload your own interested ABI file."),
                                      p("Please allow a few minutes for data processing since the app is hosted on a free server. We appreciate your patience!", class = "highlight"),
                                      p("Notifications and error messages will appear in the bottom right corner. Please read them as needed."),
                                      
                                      h5("Important:"),
                                      p("To generate parsing results, it's essential to either paste your RefSeq sequence or upload it to the corresponding area."),
                                      p("Please ensure that you provide the reference sequence to enable accurate parsing results.")
                                    ),
                                    div(class = "footer", "Designed by Huaien.Wang & Jas.B")
                           )
                           ,
                           tabPanel("Chromatogram & Quality plot",
                                    dashboardBody(
                                      fluidRow(
                                        box(
                                          title = "", width = 2, solidHeader = TRUE, status = "primary",
                                          
                                          fileInput("file1", h6("Please Upload your ABI file here!"),accept = c('.ab1'))
                                        ),
                                        
                                        box(title = "", width = 2, solidHeader = TRUE, status = "primary",
                                            br(),
                                            #br(),
                                            radioButtons("Direction", "Your Read Direction Feature",  choices=list("Forward" = 1, "Reverse" = 2), selected = 1)
                                        ),
                                        box(title = "", width = 2, solidHeader = FALSE, status = "primary",
                                            numericInput("Ratio2",h6("Signal Ratio Cutoff"), value = '0.33',min = "0",step = "0.01"
                                        ),p('Adjust the Ratio to highlight regions with double peaks!')),
                                        box(title = "", width = 2, solidHeader = TRUE, status = "primary",
                                            br(),
                                            br(),
                                            #br(),
                                            #br(),
                                            materialSwitch(inputId = "id_check", label = "Load Example!",value = TRUE,status = "success")),
                                        uiOutput("dynamic_sliders"),
                                        box(title = "", width = 12, solidHeader = TRUE, status = "primary",
                                            plotlyOutput("quality", height = "1600px"))
                                      ))
                           ),
                           tabPanel("Parsing Results",
                                      dashboardBody(
                                        fluidRow(
                                          box(textAreaInput(inputId = "ref", 
                                                            label = "Paste your Reference DNA Sequence here", 
                                                            rows=6,
                                                            width = "auto",
                                                            value = "ACTCTTTCCCTACACGACGCTCTTCCGATCTAAAGGGACTTGGGCTTTGGTGTTGGGCTGGTAGGCTGAGAACACAGTCCTGAGGGTTCTTCTTTGTCCCATCCACAGCCATCACCTCACTGCATGGACGATCTGTTGCTGCCCCAGGATGTTGAGGAGTTTTTTGAAGGCCCAAGTGAAGCCCTCCGAGTGTCAGGAGCTCCTGCAGCACAGGACCCTGTCACCGAGACCCCTGGGCCAGTGGCCCCTGCCCCAGCCACTCCATGGCCCCTGTCATCTTTTGTCCCTTCTCAAAAAACTTACCAGGGCAACTATGGCTTCCACCTGGGCTTCCTGCAGTCTGGGACAGCCAAGTCTGTTATGTGCACGGTGAGAGATCGGAAGAGCACACGTCTGAACTC",
                                                            resize = "both")
                                          ),
                                          box(title = "", width = 3, solidHeader = TRUE, status = "success",
                                              fileInput("file2", h6("or Upload a Reference sequence(accept .fasta format)"), accept = c(".fasta",".fa",".fna",".dna"))
                                          ),
                                          box(title = "", width = 3,
                                              solidHeader =FALSE, status = "success",style = "margin-top: 90px;",
                                              actionButton(inputId = "update_results", 
                                                           label = "Update Inputs and Show Alignments",
                                                           style = "color: white; background-color: #337ab7; border-color: #2e6da4;" )
                                          )
                                        ),
                                        fluidRow(
                                          box(title = "", width = 4, solidHeader = FALSE, status = "primary",
                                              numericInput("Ratio",h6("Signal Ratio Cutoff"), value = '0.33',min = "0",step = "0.01")
                                          ),
                                          box(title = "", width = 4, solidHeader = TRUE, status = "primary",
                                              numericInput("Trim_5", h6("Trim 5' end"), value = '100', min = "0",step = "1")
                                          ),
                                          box(title = "", width = 4, solidHeader = TRUE, status = "primary",
                                              numericInput("Trim_3", h6("Trim 3' end"), value = '60', min = "0",step = "1")
                                          )
                                        ) ,
                                        fluidRow(
                                          box(title = h5("Processing Status"), width = 12, solidHeader = TRUE, status = "warning",
                                              progressBar(id = "pb", value = 0, striped = TRUE,status = "success",
                                                          display_pct = TRUE)
                                          )
                                        ),
                                        box(title = h5("Primary sequence mapped to the Reference"), width = 12, solidHeader = TRUE, status = "primary",
                                            # Button to copy primary sequence
                                            actionButton("copy_primary_seq", "Copy Primary Sequence"),
                                            htmlOutput('alignment1'),
                                            downloadButton("downloadAlignment1", "Download Alignment"),
                                            htmlOutput("alignmentOutput1")
                                        ),
                                        box(title = h5("Secondary sequence mapped to the Reference"), width = 12, solidHeader = TRUE, status = "primary",
                                            # Button to copy secondary sequence
                                            actionButton("copy_secondary_seq", "Copy Secondary Sequence"),
                                            htmlOutput('alignment2'),
                                            downloadButton("downloadAlignment2", "Download Alignment"),
                                            htmlOutput("alignmentOutput2")
                                        ),
                                        #box(title = h5("Results"), width = 12, solidHeader = TRUE, status = "info",
                                        #      uiOutput("summaryReport")
                                        #),
                                        br(),
                                        box(materialSwitch(inputId = "id_summary", label = "Show full Summary Table!",value = FALSE,status = "success")
                                            ),
                                        box(width =12,DT::dataTableOutput("summary_table"))
                                      )
                           )
                           ,
                           tags$script(HTML("shinyjs.clickUpdateButton = function() {
                           $('#update_results').click();
                           }")),
                           tabPanel("Translation",
                                    fluidRow(
                                      column(4,
                                             textAreaInput(inputId = "translate",
                                                           label = "DNA Sequence to be translated paste here!",
                                                           rows=4,
                                                           width = 1600,
                                                           value = "CTATGGAAGGTGCACACGCATGCCATGCCATTTGTGCAGAAGACATATTTCGGGAGTTGCTTC",
                                                           resize = "both")
                                      ),
                                      column(3,
                                             selectInput("genetic_source", "Genetic Code Table:",select_choices,selected = "Standard")
                                      ),
                                      br(),
                                      column(3,
                                             selectInput("Frame", "Frame shift",c("Frame shift 1"="1",
                                                                                  "Frame shift 2"="2",
                                                                                  "Frame shift 3"="3"),selected = "Frame shift 1")

                                      ),
                                      br(),
                                      br(),
                                      column(2,
                                             materialSwitch(inputId = "id_genetic", label = "Show Genetic Code Table!",value = FALSE,status = "success")
                                      ),
                                      column(12,
                                             DT::dataTableOutput("AA_out")
                                      ),
                                      column(12,
                                             DT::dataTableOutput("source_table"))
                                      
                                    )
                           ),
                           tabPanel("References",
                                    fluidPage(#title = "IUPAC_CODE_MAP",
                                      #dataTableOutput("IUPAC_CODE_MAP")
                                      h5("Citations for R packages"),
                                      # Output for the citations
                                      verbatimTextOutput("citationOutput")
                                    )
                                    
                           )
                           
                           
                           
                )
)
readABIFile <- function(file_path) {
  tryCatch({
    readsangerseq(file_path)
  }, error = function(e) {
    stop("Error reading ABI file: ", e$message)
  })
}
count_ATCGN<-function(seq){
  num_G <- str_count(seq, "G")
  num_A <- str_count(seq, "A")
  num_T <- str_count(seq, "T")
  num_C <- str_count(seq, "C")
  num_N <- str_count(seq, "N")
  seq_length<-str_length(seq)
  GC_content<-round((num_G + num_C) / str_length(seq) * 100,2)
  count_data<-cbind(seq_length,num_G,num_A,num_T,num_C,num_N,GC_content)
  count_data <-data.frame(count_data)
  return(count_data)
}
getDiff <- function(sa, sb){
  n = 0
  for (i in 1:(nchar(sa))) {
    if (substring(sa, i, i) != substring(sb, i, i)) {
      n = n + 1
    }
  }
  return (n)
}
getTwoSeq <- function(seq1, seq2, sl, nd){
  for (i in 1:(nchar(seq1)-sl)) {
    twoSeq <-c()
    if (substring(seq1, i, i) != substring(seq2, i, i)) {
      temp1 <- substring(seq1, i, i + sl - 1)
      temp2 <- substring(seq2, i, i + sl - 1)
      if (getDiff(temp1, temp2) >= nd){
        twoSeq <- c(temp1, temp2, i)
        break
      }
    }
  }
  return (twoSeq)
}
getList <- function(SEQ1, SEQ2){
  s1 <- unlist(strsplit(SEQ1, split = ""))
  s2 <- unlist(strsplit(SEQ2, split = ""))
  strList <- c()
  inputLength <- length(s1)
  for (i in 1:inputLength){
    if (i == 1){
      for (j in 1:2**inputLength){
        if (j <= 2**(inputLength - 1)){
          strList <- c(strList, s1[1])
        }else{
          strList <- c(strList, s2[1])
        }
      }
    }else{
      for (k in 1:2**inputLength){
        if (k %% 2**(inputLength-i+1) != 0 & (k %% 2**(inputLength-i+1)) <= 2**(inputLength-i)){
          strList[k] <- paste(strList[k],s1[i], sep = "") 
        }else{
          strList[k] <- paste(strList[k],s2[i], sep = "") 
        }
      }
    }
  }
  return (strList[!duplicated(strList)])
}
getOutputSequences <- function(refSeq, seq1, seq2){
  out1 <- "" # the first sequence
  out2 <- "" # the second sequence
  twoSeq <- getTwoSeq(seq1, seq2, sl, nd)
  refStart <- -1
  if (! is.null(twoSeq)){
    theList <- getList(twoSeq[1], twoSeq[2])
    for (s in theList){
      if (! is.na(str_locate(refSeq, s)[1])){
        refStart <- str_locate(refSeq, s)[1]
        refSeqFinal <- substring(refSeq, refStart - strtoi(twoSeq[3]) + 2, refStart + nchar(seq1) - strtoi(twoSeq[3]) + 1)
        break
      }
    }
    for (i in twoSeq[3]:nchar(seq1)) {
      
      if (substring(seq1, i, i) == substring(refSeqFinal, i-1, i-1)) {
        out1 <- paste(out1, substring(seq1, i, i), sep="") 
        if (substring(seq1, i, i) == substring(seq2, i, i)){
          out2 <- paste(out2, substring(seq1, i, i), sep="")
        }else{
          out2 <- paste(out2, substring(seq2, i, i), sep="")
        }
      }else{
        out1 <- paste(out1, substring(seq2, i, i), sep="")
        out2 <- paste(out2, substring(seq1, i, i), sep="")
      }
    }
    
  }
  out1 <- paste(substring(seq1, 1, strtoi(twoSeq[3]) - 1), out1, sep = "") 
  out2 <- paste(substring(seq2, 1, strtoi(twoSeq[3]) - 1), out2, sep = "") 
  return (c(out1, out2, refSeqFinal))
}
sl = 10 # variable 6
nd = 3 # variable 7

# get_citation_text <- function(package) {
#   tryCatch({
#     citation_text <- capture.output(citation(package))
#     paste(citation_text, collapse = "\n")
#   }, error = function(e) {
#     paste("Error retrieving citation for", package, ":", e$message)
#   })
# }
# packages <- c("shiny", "bslib", "sangerseqR", "Biostrings", "ggplotify", "plotly",
#               "data.table", "stringr", "stringi", "DT", "shinyBS", "DECIPHER",
#               "BiocManager", "ggplot2", "ggbreak", "shinyWidgets", "shinyalert",
#               "tidyr", "msa", "shinyjs", "shinydashboard", "shinyWidgets")
# packages <- unique(packages)
# all_citations_text <- sapply(packages, get_citation_text, USE.NAMES = FALSE)
all_citations_text <- "1.Chang, W. et al. Shiny: Web Application Framework for R. R package version 1.8.0 (2023). Available at: https://CRAN.R-project.org/package=shiny.

2.Sievert, C. et al. Bslib: Custom 'Bootstrap' 'Sass' Themes for 'Shiny' and 'Rmarkdown'. R package version 0.6.0 (2023). Available at: https://CRAN.R-project.org/package=bslib.

3.Hill, J.T. et al. SangerseqR: Poly Peak Parser: Method and Software for Identification of Unknown Indels Using Sanger Sequencing of Polymerase Chain Reaction Products. Developmental Dynamics 243:1632-1636 (2014).

4.Pagès, H. et al. Biostrings: Efficient Manipulation of Biological Strings. R package version 2.68.1 (2023). Available at: https://bioconductor.org/packages/Biostrings.

5.Yu, G. Ggplotify: Convert Plot to 'Grob' or 'Ggplot' Object. R package version 0.1.2 (2023). Available at: https://CRAN.R-project.org/package=ggplotify.

6.Sievert, C. Interactive Web-Based Data Visualization with R, Plotly, and Shiny. Chapman and Hall/CRC Florida (2020).

7.Dowle, M. et al. Data.table: Extension of Data.frame. R package version 1.14.8 (2023). Available at: https://CRAN.R-project.org/package=data.table.

8.Wickham, H. Stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.5.1 (2023). Available at: https://CRAN.R-project.org/package=stringr.

9.Gagolewski, M. Stringi: Fast and Portable Character String Processing in R. Journal of Statistical Software, 103(2), 1-59 (2022). DOI: 10.18637/jss.v103.i02.

10.Xie, Y. et al. DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.30 (2023). Available at: https://CRAN.R-project.org/package=DT.

11.Bailey, E. ShinyBS: Twitter Bootstrap Components for Shiny. R package version 0.61.1 (2022). Available at: https://CRAN.R-project.org/package=shinyBS.

12.Wright, E.S. Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R. The R Journal, 8(1), 352-359 (2016).

13.Morgan, M. et al. BiocManager: Access the Bioconductor Project Package Repository. R package version 1.30.22 (2023). Available at: https://CRAN.R-project.org/package=BiocManager.

14.Wickham, H. Ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York (2016).

15.Xu, S. et al. Use Ggbreak to Effectively Utilize Plotting Space to Deal with Large Datasets and Outliers. Frontiers in Genetics. 2021, 12:774846. DOI: 10.3389/fgene.2021.774846.

16.Perrier, V. et al. ShinyWidgets: Custom Inputs Widgets for Shiny. R package version 0.8.0 (2023). Available at: https://CRAN.R-project.org/package=shinyWidgets.

17.Attali, D. et al. Shinyalert: Easily Create Pretty Popup Messages (Modals) in 'Shiny'. R package version 3.0.0 (2021). Available at: https://CRAN.R-project.org/package=shinyalert.

18.Wickham, H. et al. Tidyr: Tidy Messy Data. R package version 1.3.0 (2023). Available at: https://CRAN.R-project.org/package=tidyr.

19.Bodenhofer, U. et al. Msa: An R Package for Multiple Sequence Alignment. Bioinformatics 31(24):3997-9999 (2015). DOI: 10.1093/bioinformatics/btv176.

20.Attali, D. Shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 2.1.0 (2021). Available at: https://CRAN.R-project.org/package=shinyjs.

21.Chang, W. et al. Shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.2 (2021). Available at: https://CRAN.R-project.org/package=shinydashboard.

22.D.R. Lide, Handbook of Chemistry and Physics, 72nd Edition, CRC Press, Boca Raton, FL, 1991.
"

# Define server logic ----
server <- function(input, output,session) {
  fileUploaded <- reactiveVal(FALSE)
  observe({
    if (!is.null(input$file1)) {
      fileUploaded(TRUE)  # Update the file uploaded status
      updateMaterialSwitch(session, "id_check", value = FALSE)  # Turn off the switch
      showNotification("Your file has been uploaded. Unloading the example data.", type = "message")
    }
  })
  observe({
    shinyjs::runjs("shinyjs.clickUpdateButton()")
  })
  # Output the citations text
  output$citationOutput <- renderText({ all_citations_text })
  #define function counting ATCGNs from ABI file
  abiData <- reactive({
    req(input$file1)
    readABIFile(input$file1$datapath)
  })
  #stats of sequences
  output$info <- renderTable({
    req(abiData())
    data_ATCG <- countATCGN(abiData()$primaryseq)
  })
  # Code to process file1
  observeEvent(input$file1, {
    tryCatch({
      # Code to process file1
    }, error = function(e) {
      showNotification("Error: ", e$message, type = "error")
    })
  })
  #alignment errors
  output$alignment1 <- renderUI({
    tryCatch({
      req(input$file1, input$file2)
      # Alignment logic
    }, error = function(e) {
      showNotification("trigger the update button to show alignments", e$message, type = "error")
    })
  })
  output$alignment2 <- renderUI({
    tryCatch({
      req(input$file1, input$file2)
      # Alignment logic
    }, error = function(e) {
      showNotification("trigger the update button to show alignments", e$message, type = "error")
    })
  })
  #show chromatogram plot and quality plot errors
  output$quality <- renderPlotly({
    tryCatch({
      req(abiData())
      # Plotting logic
    }, error = function(e) {
      showNotification("Plotting Error: ", e$message, type = "error")
    })
  })
  #calculate the counts for each nucleotides
  ATCGN_abi<-function(abifile){
    obj<-readsangerseq(abifile)
    primaryseq<-substring(primarySeq(obj,string = TRUE),input$Trim_5+1,nchar(primarySeq(obj,string = TRUE))-input$Trim_3)
    secondaryseq<-substring(secondarySeq(obj,string = TRUE),input$Trim_5+1,nchar(secondarySeq(obj,string = TRUE))-input$Trim_3)
    count_data<-rbind(count_ATCGN(primaryseq),count_ATCGN(secondaryseq))
    rownames(count_data)<-c("PrimarySeq","SecondarySeq")
    return(count_data)
  }
  #output basic sequence stats from inputs
  observeEvent(input$id_check, {
    
    if(input$id_check=="TRUE"){
      output$info <-renderTable({
        data_ATCG<-ATCGN_abi("www/AP009_6_Trp53_Long_F.ab1")
        rownames(data_ATCG)<-c("Primary Sequence","Secondary Sequence")
        colnames(data_ATCG)<-c("Sequence Length","Count of G","Count of A","Count of T","Count of C","Count of N","GC content(%)")
        data_ATCG
      }, rownames = TRUE,width=1400,hover = TRUE)
    }else{
      output$info <-renderTable({
        inFile1 <- input$file1
        if (is.null(inFile1))
          return(NULL)
        data_ATCG<-ATCGN_abi(inFile1$datapath)
        rownames(data_ATCG)<-c("Primary Sequence","Secondary Sequence")
        colnames(data_ATCG)<-c("Sequence Length","Count of G","Count of A","Count of T","Count of C","Count of N","GC content(%)")
        data_ATCG
      }, rownames = TRUE,width=1400,hover = TRUE)
    }
  })
  # Render only when input$id_genetic is TRUE
  observeEvent(input$id_genetic, {
    if(input$id_genetic) {
      print("observeEvent genetic code table triggered")
      output$source_table <- renderDT({
        data1<-as.data.frame(GENETIC_CODE_TABLE)
        datatable(
          data1, extensions = c('Buttons','FixedColumns'), 
          options = list(
            searchHighlight = TRUE,
            scrollY = 1000,
            scroller = TRUE,
            scrollX= TRUE,
            fixedColumns = TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
            lengthMenu = list(c(50,100, -1) # declare values
                              , c(50,100, "All") # declare titles
            ) # end of lengthMenu customization
            , pageLength = -1
          )
        )
      }) 
    }else{
      output$source_table <- renderDT({
        NULL
      })
    }
     #load user input data
  })
  ab1File<-reactive({
    inFile1 <- input$file1
    if (is.null(inFile1))
      return(NULL)
    ab1File <- readsangerseq(inFile1$datapath)
  })
  #peak data extracted from abi file
  peak_data<- reactive({
    if(input$id_check==TRUE){
      file_name<-"www/AP009_6_Trp53_Long_F.ab1"
      obj<-readsangerseq(file_name)
      seq.abif<-read.abif(file_name)
      primarypeaks <- obj@peakPosMatrix[,1]
      Apeaks <- obj@traceMatrix[,1][primarypeaks]
      Cpeaks <- obj@traceMatrix[,2][primarypeaks]
      Gpeaks <- obj@traceMatrix[,3][primarypeaks]
      Tpeaks <- obj@traceMatrix[,4][primarypeaks]
      trace<-as.data.frame(cbind(Apeaks,Cpeaks,Gpeaks,Tpeaks))
      colnames(trace)<-c("A","C","G","T")
      qualityscore<-seq.abif@data$PCON.1
      #-find peak channel-------------------------------------------------------------
      peak_data<-trace
      #find second highest signal
      for (i in 1:nrow(peak_data)){
        #print(i)
        row<-peak_data[i,1:4]
        new_row<-apply(row,1,sort)
        ratio2_1<-round(new_row[3,]/new_row[4,],3)
        ratio3_1<-round(new_row[2,]/new_row[4,],3)
        ratio4_1<-round(new_row[1,]/new_row[4,],3)
        peak_data[i,5:8]<-rownames(new_row)
        peak_data[i,9:11]<-c(ratio2_1,ratio3_1,ratio4_1)
      }
      colnames(peak_data)[5:11]<-c("peak_4th","peak_3rd","peak_2nd","peak_1st","ratio_peak2_peak1","ratio_peak3_peak1","ratio_peak4_peak1")
      peak_data$pos<-as.numeric(rownames(peak_data))
      peak_data$quality<-seq.abif@data$PCON.1
      peak_data
    }else{
      inFile1 <- input$file1
      if (is.null(inFile1))
        return(NULL)
      obj<-readsangerseq(inFile1$datapath)
      seq.abif<-read.abif(inFile1$datapath)
      primarypeaks <- obj@peakPosMatrix[,1]
      Apeaks <- obj@traceMatrix[,1][primarypeaks]
      Cpeaks <- obj@traceMatrix[,2][primarypeaks]
      Gpeaks <- obj@traceMatrix[,3][primarypeaks]
      Tpeaks <- obj@traceMatrix[,4][primarypeaks]
      trace<-as.data.frame(cbind(Apeaks,Cpeaks,Gpeaks,Tpeaks))
      colnames(trace)<-c("A","C","G","T")
      qualityscore<-seq.abif@data$PCON.1
      #-find peak channel-------------------------------------------------------------
      peak_data<-trace
      #find second highest signal
      for (i in 1:nrow(peak_data)){
        #print(i)
        row<-peak_data[i,1:4]
        new_row<-apply(row,1,sort)
        ratio2_1<-round(new_row[3,]/new_row[4,],3)
        ratio3_1<-round(new_row[2,]/new_row[4,],3)
        ratio4_1<-round(new_row[1,]/new_row[4,],3)
        peak_data[i,5:8]<-rownames(new_row)
        peak_data[i,9:11]<-c(ratio2_1,ratio3_1,ratio4_1)
      }
      colnames(peak_data)[5:11]<-c("peak_4th","peak_3rd","peak_2nd","peak_1st","ratio_peak2_peak1","ratio_peak3_peak1","ratio_peak4_peak1")
      peak_data$pos<-as.numeric(rownames(peak_data))
      peak_data$quality<-seq.abif@data$PCON.1
      peak_data
      
    }
  })
  #update the trimming for chromatogram and quality plot
  output$dynamic_sliders <- renderUI({
    req(peak_data())  # Ensure peak_data is available before rendering sliders
    max_base_pairs <- nrow(peak_data()) - 1  # Assuming peak_data rows correspond to base pairs
    
    tagList(
      sliderInput(
        "cut_5", 
        label = HTML("Trim 5' End: <span title=\"Slide to remove bases from the left end (5's end) of the sequence.\">?</span>"),
        value = 0, 
        min = 0, 
        max = max_base_pairs, 
        step = 1
      ),
      
      sliderInput(
        "cut_3", 
        label = HTML("Trim 3' End: <span title=\"Slide to remove bases from the right end (3's end) of the sequence.\">?</span>"),
        value = 0, 
        min = 0, 
        max = max_base_pairs, 
        step = 1
      )
    )
  })
  #output quality plot
  output$quality<-renderPlotly({
    peak_data <- req(peak_data())
    
    if(!is.null(peak_data)){
      peak_data$QC<-"unknown"
      for (i in 1:nrow(peak_data)){
        #print(i)
        if (peak_data[i,"quality"]<20){
          peak_data[i,"QC"]<-"Bad(score < 20)"
        }else if(peak_data[i,"quality"]>=20){
          peak_data[i,"QC"]<-"Good(score >= 20)"
        }
      }
      from <- max(input$cut_5 + 1, 1)  # Ensure from is at least 1
      to <- min(nrow(peak_data) - input$cut_3, nrow(peak_data))  # Ensure to is not greater than nrow(peak_data)
      if (from <= to) {
        peak_data <- peak_data[from:to,]
      } else {
        # Handle the case when the range is invalid, perhaps by returning NULL or a message
        return(NULL)
      }
      
      req(nrow(peak_data) > 0)
      peak_data_pivot<-peak_data[,c(1:4,8,9,12,13,14)]%>%
        pivot_longer(cols = c(1:4),names_to = "Channel",values_to = "Signal")
      data_1<-peak_data_pivot
      max_signal<-max(data_1$Signal)
      yrange<-max_signal
      expand_dataset_with_zeros <- function(data) {
        # Expand the dataset by duplicating each row three times
        expanded_data <- data %>%
          rowwise() %>%
          do(data.frame(pos = rep(.$pos, 2),
                        quality = c(NA, .$quality),
                        Channel = rep(.$Channel, 2),
                        QC = rep(.$QC, 2),
                        ratio = rep(.$ratio, 2),
                        Signal = c(0, .$Signal))) %>%
          ungroup()
        
        return(expanded_data)
      }
      expanded_data <- expand_dataset_with_zeros(data_1)
      combined_expanded_data <- expanded_data %>%
        group_by(Channel) %>%
        mutate(expanded_x = row_number()) %>%
        ungroup()
      unique_mapping <- combined_expanded_data %>%
        distinct(expanded_x, pos) %>%
        arrange(expanded_x)
      combined_expanded_data$hover_text <- ifelse(combined_expanded_data$expanded_x %% 2 == 0, 
                                                  paste("Position: ", combined_expanded_data$pos, 
                                                        "<br>Signal: ", combined_expanded_data$Signal, 
                                                        "<br>Channel: ", combined_expanded_data$Channel),
                                                  "")
      cutoff <- input$Ratio2
      
      # Filter positions with ratio above cutoff
      highlight_positions <- combined_expanded_data %>%
        filter(ratio > cutoff) %>%
        distinct(expanded_x) %>%
        arrange(expanded_x)
      
      # Create shapes for highlighting
      shapes_list <- lapply(1:nrow(highlight_positions), function(i) {
        list(
          type = 'rect',
          x0 = highlight_positions$expanded_x[i] - 0,  # Adjust as necessary for the width of the highlight
          x1 = highlight_positions$expanded_x[i] + 1,
          y0 = 0,
          y1 = Inf,
          fillcolor = 'blue',
          opacity = 0.2,
          line = list(width = 0)
        )
      })
      # Now, you can plot this combined dataset directly
      fig <- plot_ly(data = combined_expanded_data, x = ~expanded_x, y = ~Signal, type = 'scatter', mode = 'lines',
                     hoverinfo = 'text',  text = ~hover_text, 
                     line = list(shape = "spline"), color = ~Channel) %>%
        layout(title = '',
               xaxis = list(title = '', 
                            tickvals = NULL,  # Do not specify tick values
                            ticktext = NULL,  # Do not specify tick labels
                            showticklabels = FALSE,  # Hide tick labels
                            rangeslider = list(visible = TRUE)),
               yaxis = list(title = 'Signal Intensity',range = list(0, yrange)),
               shapes = shapes_list)
      combined_expanded_data$custom_hover_text <- ifelse(combined_expanded_data$expanded_x %% 2 == 0, 
                                                         paste("Position: ", combined_expanded_data$pos, 
                                                               "\nQuality: ", combined_expanded_data$quality, 
                                                               "\nQC: ", combined_expanded_data$QC),
                                                         "")
      p0 <- ggplot(data = combined_expanded_data, aes(x = expanded_x, y = quality, color = QC, text = custom_hover_text)) +
        geom_point() +
        geom_rug() +
        ylab("Quality Score(Phred 33)") +
        geom_hline(yintercept = 20, color = "red", linetype = "dashed") +
        xlab("") +
        theme_minimal() +
        theme(legend.position = "top") +
        #xlim(from, to) + # Set x-axis limits
        scale_color_manual(values = c("Good(score >= 20)" = "darkgreen", "Bad(score < 20)" = "maroon"))
      
      fig1 <- ggplotly(p0, tooltip = "text") %>% layout(xaxis = list(title = '', 
                                                   tickvals = NULL,  # Do not specify tick values
                                                   ticktext = NULL,  # Do not specify tick labels
                                                   showticklabels = FALSE,
                                                   rangeslider = list(visible = TRUE)))
      fig_all <- subplot( fig, fig1,nrows = 2, shareX = TRUE) %>%
        layout(title = list(text = "Chromatogram Plot with Quality Scores(Phred 33)", 
                            font = list(size = 16),
                            family = "Arial, bold",
                            pad = list(t = 400) ),
               #margin = list(t = 100), 
               legend = list(title = list(text = ''), orientation = "h", xanchor = "center", x = 0.5, y = -0.6),
               hovermode = "x"
        )
      
      # Return the combined plot
      fig_all%>% layout(height = 800)%>%config(displayModeBar = FALSE)
      #fig_all
    }else{
      return(NULL) 
    }
    
  })
  # validate the user input
  valid_nucleotides <- names(Biostrings::IUPAC_CODE_MAP)
  # processing refseq
  refSeq <- reactive({
    if (!is.null(input$ref) && input$ref != "") {
      unique_chars <- unique(unlist(strsplit(as.character(toupper(input$ref)), "")))
      invalid_chars <- setdiff(unique_chars, valid_nucleotides)
      if (length(invalid_chars) > 0) {
        stop(paste("Invalid DNA sequence. The following characters are not valid nucleotides:", paste(invalid_chars, collapse=", ")))
      } else {
        refSeq<-toupper(input$ref)
      }
      
    }else{
      refSeq<-readDNAStringSet(inFile2$datapath)
      unique_chars <- unique(unlist(strsplit(as.character(refSeq), "")))
      invalid_chars <- setdiff(unique_chars, valid_nucleotides)
      if (length(invalid_chars) > 0) {
        stop(paste("Invalid DNA sequence. The following characters are not valid nucleotides:", paste(invalid_chars, collapse=", ")))
      } else {
        refSeq
      }
    }
  })
  primaryseq <- reactiveVal()
  secondaryseq <- reactiveVal()
  mismatch_table<-reactive({
    primary_alignments <- pairwiseAlignment(primaryseq(), refSeq(),type="global-local")
    secondary_alignments <- pairwiseAlignment(secondaryseq(), refSeq(),type="global-local")
    mistable1<-mismatchTable(primary_alignments)
    if (nrow(mistable1)>0){
      mistable1<-cbind(mistable1,"primary")
      colnames(mistable1)[10]<-"SequenceType"
    }
    mistable2<-mismatchTable(secondary_alignments)
    if (nrow(mistable2)>0){
      mistable2<-cbind(mistable2,"secondary")
      colnames(mistable2)[10]<-"SequenceType"
    }
    mistable<-rbind(mistable1,mistable2)
  })
  observeEvent(input$update_results, {
    reactiveValuesToUpdate$Ratio <- input$Ratio
    reactiveValuesToUpdate$Trim_5 <- input$Trim_5
    reactiveValuesToUpdate$Trim_3 <- input$Trim_3
    reactiveValuesToUpdate$id_check <- input$id_check
    reactiveValuesToUpdate$refSeq <- refSeq()
    if(reactiveValuesToUpdate$id_check=="TRUE"){
      data<-readsangerseq("www/AP009_6_Trp53_Long_F.ab1")
      refSeq<-reactiveValuesToUpdate$refSeq
      ratio <- reactiveValuesToUpdate$Ratio
      ab1BaseCall <- makeBaseCalls(data, ratio=as.numeric(ratio))
      ab1FirstCall <- primarySeq(ab1BaseCall, string = TRUE)
      ab1SecondCall <- secondarySeq(ab1BaseCall, string = TRUE)
      trim5 <- reactiveValuesToUpdate$Trim_5
      trim3 <- reactiveValuesToUpdate$Trim_3
      seq1_1<- substring(toString(ab1FirstCall), as.numeric(trim5+1), nchar(ab1FirstCall)-as.numeric(trim3))
      seq2_1<- substring(toString(ab1SecondCall), as.numeric(trim5+1), nchar(ab1SecondCall)-as.numeric(trim3))
      twoSequences <- getOutputSequences(refSeq, seq1_1, seq2_1)
      seq1 <-twoSequences[1]
      seq2 <-twoSequences[2]
      refSeq <-twoSequences[3]
      primaryseq(seq1)
      secondaryseq(seq2)
    }else{
      print("unloaded example data!")
      inFile1 <- input$file1
      inFile2 <- input$file2
      if (is.null(inFile1))
        return(NULL)
      data<-ab1File()
      ratio <- reactiveValuesToUpdate$Ratio
      ab1BaseCall <- makeBaseCalls(data, ratio=as.numeric(ratio))
      ab1FirstCall <- primarySeq(ab1BaseCall, string = TRUE)
      ab1SecondCall <- secondarySeq(ab1BaseCall, string = TRUE)
      trim5 <- reactiveValuesToUpdate$Trim_5
      trim3 <- reactiveValuesToUpdate$Trim_3
      primaryseq_seq1 <- substring(toString(ab1FirstCall), as.numeric(trim5+1), nchar(ab1FirstCall)-as.numeric(trim3))
      secondaryseq_seq1 <- substring(toString(ab1SecondCall), as.numeric(trim5+1), nchar(ab1SecondCall)-as.numeric(trim3))
      twoSequences <- getOutputSequences(refSeq, primaryseq_seq1, secondaryseq_seq1)
      primaryseq_seq <-twoSequences[1]
      secondaryseq_seq<-twoSequences[2]
      if (input$Direction == 1) {
        primaryseq(primaryseq_seq)
        secondaryseq(secondaryseq_seq)
      } else {
        primaryseq(stri_reverse(rev_complement(primaryseq_seq)))
        secondaryseq(stri_reverse(rev_complement(secondaryseq_seq)))
      }
    }
  })
  # JavaScript for copying primary sequence
  observeEvent(input$copy_primary_seq, {
    shinyjs::runjs(sprintf("navigator.clipboard.writeText('%s');", primaryseq()))
  })
  # JavaScript for copying secondary sequence
  observeEvent(input$copy_secondary_seq, {
    shinyjs::runjs(sprintf("navigator.clipboard.writeText('%s');", secondaryseq()))
  })
  observeEvent(input$update_results, {
    # Generate the primary and secondary alignments
    primary_alignments <- pairwiseAlignment(pattern = primaryseq(), subject = refSeq(), type = "global-local")
    secondary_alignments <- pairwiseAlignment(pattern = secondaryseq(), subject = refSeq(), type = "global-local")
    
    # Write the primary alignments to a file
    writePairwiseAlignments(primary_alignments,block.width = 100, file = "primary_alignments.txt")
    # Write the secondary alignments to a file
    writePairwiseAlignments(secondary_alignments,block.width = 100, file = "secondary_alignments.txt")
  })
  # Read the primary alignment file and output it
  observeEvent(input$update_results, {
    output$alignmentOutput1 <- renderUI({
    if(file.exists("primary_alignments.txt")) {
      contents <- readLines("primary_alignments.txt")
      contents <- paste(contents, collapse = "\n")
      tags$pre(style = "white-space: pre-wrap; font-family: monospace;", contents)
    } else {
      "No primary alignment available. Please generate it."
    }
  })
  })
  # Read the secondary alignment file and output it
  observeEvent(input$update_results, {
    output$alignmentOutput2 <- renderUI({
    if(file.exists("secondary_alignments.txt")) {
      contents <- readLines("secondary_alignments.txt")
      contents <- paste(contents, collapse = "\n")
      tags$pre(style = "white-space: pre-wrap; font-family: monospace;", contents)
    } else {
      "No secondary alignment available. Please generate it."
    }
  })
  })
  # Provide download handlers for the alignments
  output$downloadAlignment1 <- downloadHandler(
    filename = function() { "primary_alignments.txt" },
    content = function(file) { file.copy("primary_alignments.txt", file) }
  )
  
  output$downloadAlignment2 <- downloadHandler(
    filename = function() { "secondary_alignments.txt" },
    content = function(file) { file.copy("secondary_alignments.txt", file) }
  )
  # functions of generating summary reports
  createReport <- function(mistable, sequence_name) {
    if (nrow(mistable) == 0) {
      return(paste(sequence_name, ": 100% mapped to the reference, reckoned as wild type."))
    } else {
      mutations <- paste(mistable$PatternStart, "->", mistable$PatternSubstring, collapse = ", ")
      return(paste(sequence_name, ": Mutations at positions:", mutations))
    }
  }
  #reverse complement function
  rev_complement<-function(seq){
    rev_seq<-NULL
    for (i in 1:nchar(seq)){
      if (substring(seq,i,i)=="A"){
        base<-"T"
        rev_seq<-c(rev_seq,base)  
      }else if(substring(seq,i,i)=="T"){
        base<-"A"
        rev_seq<-c(rev_seq,base) 
      }else if(substring(seq,i,i)=="G"){
        base<-"C"
        rev_seq<-c(rev_seq,base) 
      }else if(substring(seq,i,i)=="C"){
        base<-"G"
        rev_seq<-c(rev_seq,base) 
      }
    }
    rev_seq<-paste(rev_seq,collapse = "")
    return(rev_seq)
  }
  reactiveValuesToUpdate <- reactiveValues()
  observeEvent(input$update_results, {
    # Initialize progress bar
    showNotification("Processing... Please wait...", type = "message", duration = 5)
    
    updateProgressBar(session, "pb", value = 0)
    Sys.sleep(0.5)
    reactiveValuesToUpdate$Ratio <- input$Ratio
    reactiveValuesToUpdate$Trim_5 <- input$Trim_5
    reactiveValuesToUpdate$Trim_3 <- input$Trim_3
    reactiveValuesToUpdate$id_check <- input$id_check
    reactiveValuesToUpdate$refSeq <- refSeq()
    Sys.sleep(0.5)
    updateProgressBar(session, "pb", value = 10)
    Sys.sleep(0.5)
    updateProgressBar(session, "pb", value = 30)
    #automatically unload example data when user upload the data, processing alignment results
    observeEvent(input$id_check, {
      if(reactiveValuesToUpdate$id_check=="TRUE"){
        data<-readsangerseq("www/AP009_6_Trp53_Long_F.ab1")
        refSeq<-reactiveValuesToUpdate$refSeq
        ab1BaseCall <- makeBaseCalls(data, ratio=as.numeric(input$Ratio))
        ab1FirstCall <- primarySeq(ab1BaseCall, string = TRUE)
        ab1SecondCall <- secondarySeq(ab1BaseCall, string = TRUE)
        trim5 <- reactiveValuesToUpdate$Trim_5
        trim3 <- reactiveValuesToUpdate$Trim_3
        seq1_1 <- substring(toString(ab1FirstCall),as.numeric(trim5+1),nchar(ab1FirstCall)-as.numeric(trim3))
        seq2_1 <- substring(toString(ab1SecondCall),as.numeric(trim5+1),nchar(ab1SecondCall)-as.numeric(trim3))
        twoSequences <- getOutputSequences(refSeq, seq1_1, seq2_1)
        seq1 <-twoSequences[1]
        seq2 <-twoSequences[2]
        refSeq <-twoSequences[3]
        output$alignment1 <-renderUI({
          mySequences<-DNAStringSet(c(seq1,refSeq),use.names = TRUE)
          
          aligned<-AlignSeqs(mySequences)
          myhtml<-print(BrowseSeqs(openURL = FALSE,DNAStringSet(aligned),highlight = 0,patterns=c("A", "C", "G", "T","N"),colWidth = Inf,
                                   colors=c("#ff7f0e", "#9467bd", "#1f77b4", "#17becf","#e377c2"))
          )
          includeHTML(myhtml)
        }) 
        output$alignment2 <-renderUI({
           mySequences<-DNAStringSet(c(seq2,refSeq),use.names = TRUE)
          
          aligned<-AlignSeqs(mySequences)
          myhtml<-print(BrowseSeqs(openURL = FALSE,DNAStringSet(aligned),highlight = 0,patterns=c("A", "C", "G", "T","N"),colWidth = Inf,
                                   colors=c("#ff7f0e", "#9467bd", "#1f77b4", "#17becf","#e377c2"))
          )
          includeHTML(myhtml)
        })
      }else{
        print("unloaded example data!")
        inFile1 <- input$file1
        inFile2 <- input$file2
        if (is.null(inFile1))
          return(NULL)
        data<-ab1File()
        ratio <- reactiveValuesToUpdate$Ratio
        ab1BaseCall <- makeBaseCalls(data, ratio=as.numeric(ratio))
        ab1FirstCall <- primarySeq(ab1BaseCall, string = TRUE)
        ab1SecondCall <- secondarySeq(ab1BaseCall, string = TRUE)
        trim5 <- reactiveValuesToUpdate$Trim_5
        trim3 <- reactiveValuesToUpdate$Trim_3
        seq1_1 <- substring(toString(ab1FirstCall),as.numeric(trim5+1),nchar(ab1FirstCall)-as.numeric(trim3))
        seq2_1 <- substring(toString(ab1SecondCall),as.numeric(trim5+1),nchar(ab1SecondCall)-as.numeric(trim3))
        twoSequences <- getOutputSequences(refSeq(), seq1_1, seq2_1)
        seq1 <-twoSequences[1]
        seq2 <-twoSequences[2]
        refseq <-twoSequences[3]
        output$alignment1 <-renderUI({
                   if (input$Direction==1){
            refSeq<-reactiveValuesToUpdate$refSeq
            mySequences<-DNAStringSet(c(seq1,refSeq),use.names = TRUE)
            aligned<-AlignSeqs(mySequences)
          }else{
            rev_seq1<-stri_reverse(rev_complement(seq1))
            #rev_seq2<-stri_reverse(rev_complement(seq2))
            refSeq<-reactiveValuesToUpdate$refSeq
            mySequences<-DNAStringSet(c(rev_seq1,refSeq),use.names = TRUE)
            aligned<-AlignSeqs(mySequences)
            
          }
          myhtml<-print(BrowseSeqs(openURL = FALSE,DNAStringSet(aligned),highlight = 0,patterns=c("A", "C", "G", "T","N"),colWidth = Inf,
                                   colors=c("#ff7f0e", "#9467bd", "#1f77b4", "#17becf","#e377c2")
          ))
          includeHTML(myhtml)
        }) 
        output$alignment2 <-renderUI({
           if (input$Direction==1){
            refSeq<-reactiveValuesToUpdate$refSeq
            mySequences<-DNAStringSet(c(seq2,refSeq),use.names = TRUE)
            aligned<-AlignSeqs(mySequences)
            
          }else{
            #rev_seq1<-stri_reverse(rev_complement(seq1))
            rev_seq2<-stri_reverse(rev_complement(seq2))
            refSeq<-reactiveValuesToUpdate$refSeq
            mySequences<-DNAStringSet(c(rev_seq2,refSeq),use.names = TRUE)
            aligned<-AlignSeqs(mySequences)
            
          }
          myhtml<-print(BrowseSeqs(openURL = FALSE,DNAStringSet(aligned),highlight = 0,patterns=c("A", "C", "G", "T","N"),colWidth = Inf,
                                   colors=c("#ff7f0e", "#9467bd", "#1f77b4", "#17becf","#e377c2")
          ))
          includeHTML(myhtml)
        }) 
      }
    })
    Sys.sleep(1)
    updateProgressBar(session, "pb", value = 50)
    Sys.sleep(1)
    updateProgressBar(session, "pb", value = 70)
    #summary reports
    observeEvent(input$id_check, {
      if(reactiveValuesToUpdate$id_check) {  # If TRUE, use example data
        output$summaryReport <- renderUI({
          data<-readsangerseq("www/AP009_6_Trp53_Long_F.ab1")
          ratio <- reactiveValuesToUpdate$Ratio
          ab1BaseCall <- makeBaseCalls(data, ratio=as.numeric(ratio))
          ab1FirstCall <- primarySeq(ab1BaseCall, string = TRUE)
          ab1SecondCall <- secondarySeq(ab1BaseCall, string = TRUE)
          refSeq<-reactiveValuesToUpdate$refSeq
          trim5 <- reactiveValuesToUpdate$Trim_5
          trim3 <- reactiveValuesToUpdate$Trim_3
          seq1_1 <- substring(toString(ab1FirstCall),as.numeric(trim5+1),nchar(ab1FirstCall)-as.numeric(trim3))
          seq2_1 <- substring(toString(ab1SecondCall),as.numeric(trim5+1),nchar(ab1SecondCall)-as.numeric(trim3))
          twoSequences <- getOutputSequences(refSeq, seq1_1, seq2_1)
          seq1 <-twoSequences[1]
          seq2 <-twoSequences[2]
          refSeq <-twoSequences[3]
          primary_alignments<-pairwiseAlignment(pattern = seq1, subject = refSeq,type="global-local")
          secondary_alignments<-pairwiseAlignment(pattern = seq2, subject = refSeq,type="global-local")
          # Determine wild type based on mismatches
          mistable1 <- mismatchTable(primary_alignments)
          mistable2 <- mismatchTable(secondary_alignments)
          mistable1<- as.data.frame(mistable1)
          mistable2 <-as.data.frame(mistable2)
          # Assuming mistable1 and mistable2 are defined elsewhere in your app
          primary_report <- createReport(mistable1, "Primary sequence")
          secondary_report <- createReport(mistable2, "Secondary sequence")
          
          HTML(paste(
            "Summary report:",
            "<br>- ", primary_report,
            "<br>- ", secondary_report
          ))
          
          
        })
      } else {  # If FALSE, use user input data
        output$summaryReport <- renderUI({
          inFile1 <- input$file1
          if (is.null(inFile1)){
            return(NULL)
          } else {
            # Initialize the report variable
            report <- ""
            data<-data<-ab1File()
            ratio <- reactiveValuesToUpdate$Ratio
            ab1BaseCall <- makeBaseCalls(data, ratio=as.numeric(ratio))
            ab1FirstCall <- primarySeq(ab1BaseCall, string = TRUE)
            ab1SecondCall <- secondarySeq(ab1BaseCall, string = TRUE)
            refSeq<-reactiveValuesToUpdate$refSeq
            trim5 <- reactiveValuesToUpdate$Trim_5
            trim3 <- reactiveValuesToUpdate$Trim_3
            seq1_1 <- substring(toString(ab1FirstCall),as.numeric(trim5+1),nchar(ab1FirstCall)-as.numeric(trim3))
            seq2_1 <- substring(toString(ab1SecondCall),as.numeric(trim5+1),nchar(ab1SecondCall)-as.numeric(trim3))
            twoSequences <- getOutputSequences(refSeq, seq1_1, seq2_1)
            seq1 <-twoSequences[1]
            seq2 <-twoSequences[2]
            refseq <-twoSequences[3]
            primary_alignments<-pairwiseAlignment(pattern = seq1, subject = refSeq,type="global-local")
            secondary_alignments<-pairwiseAlignment(pattern = seq2, subject = refSeq,type="global-local")
            # Determine wild type based on mismatches
            mistable1 <- mismatchTable(primary_alignments)
            mistable2 <- mismatchTable(secondary_alignments)
            
            # Assuming mistable1 and mistable2 are defined elsewhere in your app
            primary_report <- createReport(mistable1, "Primary sequence")
            secondary_report <- createReport(mistable2, "Secondary sequence")
            
            HTML(paste(
              "Summary report:",
              "<br>- ", primary_report,
              "<br>- ", secondary_report
            ))
            
          }
          
          
        })
      }
    })
    # Update progress bar during processing
    # for(i in 1:10) {
    #   # Simulate some processing
    #   Sys.sleep(0.5)
    #   updateProgressBar(session, "pb", value = i * 10)
    # }
    #Sys.sleep(100)
    #showNotification("Processing complete!", type = "message")
    Sys.sleep(1.5)
    updateProgressBar(session, "pb", value = 90)
    Sys.sleep(1.5)
    updateProgressBar(session, "pb", value = 100)
  })
  #output mutation summary data full or concise.
  observeEvent(input$id_summary, {
    if(input$id_summary) {
      output$summary_table <- renderDT({
        data0<-mismatch_table()
        
        #---------
        datatable(
          data0, extensions = c('Buttons','Scroller'), 
          options = list(
            searchHighlight = TRUE,
            scrollY = 1000,
            scroller = TRUE,
            scrollX= TRUE,
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
            lengthMenu = list(c(50,100, -1) # declare values
                              , c(50,100, "All") # declare titles
            ) # end of lengthMenu customization
            , pageLength = -1
          )
        )
      }) 
    }
    else{
      output$summary_table <- renderDT({
        NULL
      })
      
    }
  })
  #define translation
  AA_translate<-function(genetic_source,seq,start,seq_name){
    genetic_code<-getGeneticCode(id_or_name2 = genetic_source,as.data.frame = TRUE)
    seq_frame<-strsplit(substring(seq,start,nchar(seq)), "(?<=.{3})", perl = TRUE)[[1]]
    seq_frame<-data.frame(codon = seq_frame)
    for (i in 1:nrow(seq_frame)){
      codon<-seq_frame[i,"codon"]
      #print(codon)
      if (nchar(codon)==3){
        AA<-genetic_code[codon,"AA"]
        seq_frame[i,"AA"]<-AA
      }else{
        seq_frame[i,"AA"]<-"NA"
      }
      translate<-data.frame(t(seq_frame))
      colnames(translate)<-as.numeric(rownames(seq_frame))
      names_used<-paste(seq_name,c("Codon","1-letter code"),sep=" ")
      rownames(translate)<-names_used
    }
    return(translate)
  }
  #display the amino acid translation
  output$AA_out<-renderDT({
    data_AA<-AA_translate(input$genetic_source,toupper(input$translate),input$Frame,"seq")
    data_AA_t<-as.data.frame(t(data_AA))
    AA<-readxl::read_excel("www/aa_table.xlsx")
    AA<- as.data.frame(AA)
    AA[c("pKa", "pKb", "pKx", "pl")] <- lapply(AA[c("pKa", "pKb", "pKx", "pl")], function(x) as.numeric(gsub("–", NA, x)))
    AA$pKa<-round(as.numeric(AA$pKa),2)
    AA$pKb<-round(as.numeric(AA$pKb),2)
    AA$pKx<-round(as.numeric(AA$pKx),2)
    AA$pl<-round(as.numeric(AA$pl),2)
    #AA[pKa, pKb, pKx, pl] <- lapply(AA[pKa, pKb, pKx, pl], as.numeric)
    colnames(data_AA_t)[2]<-"1-letter code"
    data_full<-left_join(data_AA_t,AA,by="1-letter code")
    rownames(data_full)<-1:nrow(data_full)
    
    datatable(data_full,
              extensions = 'Buttons', 
              rownames=TRUE,
              caption = htmltools::tags$caption(
                style = 'caption-side: bottom; text-align: left;',
                htmltools::em("pKa is the negative of the logarithm of the dissociation constant for the -COOH group."),br(),
                htmltools::em("pKb is the negative of the logarithm of the dissociation constant for the -NH3 group."),br(),
                htmltools::em("pKx is the negative of the logarithm of the dissociation constant for any other group in the molecule."),br(),
                htmltools::em( "pl is the pH at the isoelectric point.")
              ),
              options = list(
                searchHighlight = TRUE,
                scrollY = 1000,
                scroller = TRUE,
                scrollX= TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                lengthMenu = list(c(25,50,100, -1) # declare values
                                  , c(25,50,100, "All") # declare titles
                ) # end of lengthMenu customization
                , pageLength = -1
              ))
  })
  #ambiguous code
  output$IUPAC_CODE_MAP<-renderDataTable({
    IUPAC_code<-as.data.frame(IUPAC_CODE_MAP)
    IUPAC_code$code<-rownames(IUPAC_code)
    IUPAC_code<-as.data.frame(t(IUPAC_code))
    datatable(IUPAC_code,extensions = 'Buttons', rownames=FALSE,
              options = list(
                searchHighlight = TRUE,
                scrollY = 1000,
                scroller = TRUE,
                scrollX= TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                lengthMenu = list(c(50,100, -1) # declare values
                                  , c(50,100, "All") # declare titles
                ) # end of lengthMenu customization
                , pageLength = -1
              ))
  })
}
# Run the app ----
shinyApp(ui = ui, server = server)

