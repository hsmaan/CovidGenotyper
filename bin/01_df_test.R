library(data.table)
library(umap)
library(shiny)
library(reticulate)

# Load dataframe 


# RShiny 

ui <- fluidPage(
    
  fileInput(inputId = "input_fasta", label = "Fasta files"),
  
  sliderInput(inputId = "test", label = "test", min = 0, max = 100, value = 5),
  
  plotOutput("hist")
    # Input
    # Output
                  
)

server <- function(input, output) {
  
  output$hist <- renderPlot({ hist(rnorm(input$test)) })
  
}

shinyApp(ui = ui, server = server)

# Server side

u = umap.UMAP(metric = "precomputed")


# UI

