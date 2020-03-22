library(data.table)
library(shiny)
library(reticulate)
library(ggplot2)
library(ggthemes)

# Load testing data 

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("data")

pre_aligned <- loadRData("02_covid_93_aligned.RData")
pre_meta <- loadRData("02_covid_dataframe.RData")
countries <- pre_meta[,c("Geo_Location")]
rownames(countries) <- pre_meta$Accession

setwd("..")

# Load color palette 

kev_palette <- c(
  "dodgerblue2", "#E31A1C",
  "green4",
  "#6A3D9A", 
  "#FF7F00", 
  "black", "gold1",
  "skyblue2", "#FB9A99", 
  "palegreen2",
  "#CAB2D6", 
  "#FDBF6F", 
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# Source scripts

source_python("bin/umap_get.py")
source("bin/clustal_dist.R")

# RShiny 

ui <- fluidPage(
    
  fileInput(inputId = "input_fasta", label = "Add fasta files"),
  
  plotOutput("umap")
                  
)

server <- function(input, output) {
  
  dat <- reactive({
    if (is.null(input$input_fasta)) {
      dist <- dist_get_null(pre_aligned)
      umap_res <- umap_process(dist)
      rownames(umap_res) <- pre_meta$Accession
      umap_res <- cbind(umap_res, countries)
      colnames(umap_res) <- c("UMAP_1", "UMAP_2", "Country")
      return(umap_res)
    } else {
      dist <- dist_get(input$input_fasta$datapath, pre_aligned, pre_meta$Accession)
      umap_res <- umap_process(dist)
      novel <- (length(umap_res[,1]) - length(countries))
      countries <- c(countries, rep("Novel", novel))
      umap_res <- cbind(umap_res, countries)
      colnames(umap_res) <- c("UMAP_1", "UMAP_2", "Country")
      return(umap_res)
    }
  })
  
  output$umap <- renderPlot({ 
    
    ggplot(data = dat(), aes(x = UMAP_1, y = UMAP_2)) +
      theme_few() +
      geom_point(aes(color = Country), size = 2) +
      scale_color_manual(name = "Country", values = kev_palette[1:length(unique(dat()[,3]))]) +
      geom_point(data = dat()[dat()[,3] %in% c("Novel"), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 3) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.ticks.y = element_blank()) +
      theme(axis.text.y = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.y = element_text(size = 14, face = "bold")) +
      theme(axis.title.x = element_text(size = 14, face = "bold")) +
      theme(legend.title = element_text(size = 14, face = "bold")) +
      theme(legend.text = element_text(size = 12))
    
    })
  
}

shinyApp(ui = ui, server = server)




