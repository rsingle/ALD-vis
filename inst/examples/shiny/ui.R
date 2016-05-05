library(shiny)

shinyUI(navbarPage("Asymmetric LD (ALD)",
                   
  tabPanel("Upload Data", 
    sidebarLayout(
      sidebarPanel(
        fileInput('file1', 'Upload a Haplotype Frequency File',
                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
        #checkboxInput('header', 'Header', TRUE),      
        tags$a(href="http://www.uvm.edu/~rsingle/software/ALD/file_formats.html", target="_blank", "Click here to see the file formats"),

        radioButtons('sep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
        numericInput('tol', label = "tolerance (for sum of haplo.freqs)", value = 0.01, min = 0.01, max = 0.1, step = 0.01),
       #tags$hr(),
        tags$p("__________________________"),
        tags$h4("ALD plot options"),
        checkboxInput('values', 'Show ALD values', TRUE)
      ),
      mainPanel(
       #titlePanel("Asymmetric Linkage Disequilibrium (ALD)"),
        tags$h2(tags$img(src="ALD_small.png", width="18%"),"Asymmetric Linkage Disequilibrium"),
        tags$hr(),
        helpText("NOTE: the sum of hapltoype frequencies needs to be within +/- tolerance of 1.0."),
        helpText("If the sum is not within this range, you must increase the tolerance value."),
        textOutput("text_rawdata1"),
        textOutput("text_rawdata2"),
        tags$head(tags$style("#text_rawdata2{color: red;}") ),
       #tags$img(src="ALD.png", width="50%"),
        tags$hr(),
        plotOutput('heatmap')
      )
    )
  ),

  tabPanel("ALD Table", 
           htmlOutput("text_ALDdata1"),
           dataTableOutput('plot_data')
  ),
  
  tabPanel("Allele Specific Homozygosity", 
    sidebarLayout( 
      sidebarPanel(
        tags$h4("Allele Specific Homozygosity"),
        uiOutput('choose_locus_pair'), #input$ var: selected_pair
        uiOutput('choose_locus')       #input$ var: selected_locus
      ),
      mainPanel(
        tags$head(
          tags$script(src="d3.js"),
          tags$script(src="lodash.js"),
          tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
        ),
        uiOutput("asf_display"),
        htmlOutput("maxVal_f"),
        htmlOutput("minVal_f"),
        htmlOutput("maxVal_h"),
        htmlOutput("minVal_h")        
      )
    )
  ),

  tabPanel("Raw Data", 
           dataTableOutput('raw_data')
  ),
  
  tabPanel("Citation", 
           includeMarkdown('cite.md')
  )

))

  