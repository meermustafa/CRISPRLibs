#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/


# setwd('~/Documents/NYU/CRISPRLibs/')


library(shiny)

# # SERVER
# ##############################################################################
# # Define server logic required to draw a histogram
# server <- function(input, output) {
#   
#   output$distPlot <- renderPlot({
#     # generate bins based on input$bins from ui.R
#     x    <- iris$Sepal.Length 
#     bins <- seq(min(x), max(x), length.out = input$bins + 1)
#     
#     # draw the histogram with the specified number of bins
#     hist(x, breaks = bins, col = 'darkgray', border = 'white')
#   })
# }
# 
# 
# 
# 
# USER INTERFACE
##############################################################################
# Define UI for application that draws a histogram
# ui <- fluidPage(
#   
#   # Application title
#   titlePanel("Old Faithful Geyser Data"),
#   
#   # Sidebar with a slider input for number of bins 
#   sidebarLayout(
#     sidebarPanel(
#       sliderInput("bins",
#                   "Number of bins:",
#                   min = 1,
#                   max = 50,
#                   value = 30)
#     ),
#     
#     # Show a plot of the generated distribution
#     mainPanel(
#       plotOutput("distPlot")
#     )
#   )
# )
# 
# 
# 
# 
# # Run the application 
# shinyApp(ui = ui, server = server)








############### MORE COMPLEX EXAMPLE WITH DOWNLOAD BUTTON

# initially run 
source(file = 'create_CRISPR_Cas9_library.R')
create_Cas9_library(PAM = 'NGG', chr = 'chr8', chrstart = 125000000, chrend = 125000500)
head(combined_top_and_bottom_strand_library)

#### DOWNLOAD BUTTON -----
#server.R

server = shinyServer(function(input, output) {
  datasetInput <- reactive({
    
    # Fetch the appropriate data object, depending on the value
    # of input$dataset.
    switch(input$dataset,
           "SpCas9" = combined_top_and_bottom_strand_library, 
           "LbCpf1" = pressure,
           "Cars" = cars)
  })
  
  # render the output, which is a table
  output$table <- renderTable({
    datasetInput()
  })
  
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$dataset,  # e.g. SpCas9
            input$filetype,  # e.g. csv 
            sep = ".")  # puts the previous two names together = SpCas9.csv
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype, 
                    "tsv" = "\t",
                    "csv" = ","
                    )
      
      # Write to a file specified by the 'file' argument
      write.table(datasetInput(), 
                  file, 
                  sep = sep,
                  row.names = FALSE)
    }
  )
  
  
  # TRY IMPLEMENTING PLOT e.g. input$plotOfSpacingBetweenGuideCuts
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    x = iris$Sepal.Length
    bins = seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
  })
  
  
  
})


#ui.R
ui = shinyUI(fluidPage(
  titlePanel('Design your own CRISPR libraries.'),
  sidebarLayout(
    sidebarPanel(
      
      # show what selectors the user has
      selectInput("dataset", "Choose a CRISPR enzyme:", 
                  choices = c("SpCas9", "LbCpf1", "Cars")),
      
      # choose file type to download
      radioButtons("filetype", "File type:",
                   choices = c("csv", "tsv")),
      
      # give option to download the table
      downloadButton('downloadData', 'Download')
    ),
    
    # mainPanel(
    #   tableOutput('table')
    # )
    
    
    #multiple tab page
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("distPlot")),
        tabPanel("Summary", verbatimTextOutput("summary")),
        tabPanel("Table", tableOutput("table"))
      )
    )
    
    
    
  )
))


# Run the application 
shinyApp(ui = ui, server = server)






# deploy app to the internet

# 1. enable a connection:
# rsconnect::setAccountInfo(name='meermustafa',
#                           token='E79294993B9449087DB96548EA8D99A3',
#                           secret='PQQlLwvEs7l8DnSu6KPH+0ByBW1cPjBFzytn8all')
# 
# # 2. connect
# rsconnect::deployApp('~/Documents/NYU/CRISPRLibs/')
