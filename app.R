library(data.table)
library(shiny)
library(shinythemes)
library(reticulate)
library(ggplot2)
library(ggthemes)
library(dendextend)
library(plyr)
library(colorspace)

# Load testing data 

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("data")

pre_aligned <- loadRData("dec_aligned_fastas.RData")
pre_meta <- loadRData("02_covid_dataframe.RData")
pre_dist <- loadRData("dec_fasta_dist.RData")
countries <- as.data.frame(pre_meta[,c("Accession", "Geo_Location")])

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
source("bin/cmds.R")

# RShiny 

ui <- fluidPage(theme = shinytheme("flatly"),
    
  titlePanel("Covid-19 Genotyping Tool (CGT)"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      h4(strong("Instructions")), 
      
      p(
        "Add covid-19 sequences in fasta format, in one file.",
        br(),
        "Please ensure fasta sequences have headers."
      ),
      
      fileInput(inputId = "input_fasta", label = "Upload fasta"),
      
      p(
        "Developed by Hassaan Maan",
        br(),
        "Bo Wang Lab",
        br(),
        "University Health Network",
        br(),
        "Toronto, Canada"
      )
      
    ),
  
    mainPanel(
    
      tabsetPanel(type = "tabs", tabPanel("MDS", plotOutput("cmds")), tabPanel("UMAP", plotOutput("umap")), tabPanel("HClust", plotOutput("dend")))
    
  )
  )
                  
)

server <- function(input, output) {
  
  dist_reac <- reactive ({
    if (is.null(input$input_fasta)) {
      dist <- pre_dist
      return(dist)
    } else {
      dist <- dist_get(input$input_fasta$datapath, pre_aligned)
      return(dist)
    }
  })
  
  new_countries <- reactive ({
    if (is.null(input$input_fasta)) {
      return(countries)
    } else {
      new_accessions <- rownames(dist_reac())[(nrow(countries)+1):nrow(dist_reac())]
      countries_new <- data.frame("Accession" = new_accessions, "Geo_Location" = paste("Novel", seq(1, length(new_accessions), 1)))
      new_countries <- rbind(countries, countries_new)
      return(new_countries)
    }
  })
  
  mds <- reactive ({
      cmds <- multi_dim_scale(dist_reac())
      cmds <- merge(cmds, new_countries())
      return(cmds)
  })
  
  umap <- reactive ({
    umap_res <- umap_process(dist_reac())
    acc_names = rownames(dist_reac())
    umap_df <- data.frame("Accession" = acc_names, "UMAP_1" = umap_res[,1], "UMAP_2" = umap_res[,2])
    umap_df <- merge(umap_df, new_countries())
    return(umap_df)
  })
  
  dend <- reactive ({
    h_clust <- hclust(as.dist(dist_reac()))
    h_dend <- as.dendrogram(h_clust)
    h_dend <- raise.dendrogram(h_dend, median(get_branches_heights(h_dend)))
    return(h_dend)
  })
  
  dend_cols <- reactive ({
    cols_unique <- kev_palette[1:length(unique(new_countries()[,2]))]
    cols_assigned <- cols_unique[factor(new_countries()[,2])]
    return(cols_assigned)
  })

  output$cmds <- renderPlot ({
    ggplot(data = mds(), aes(x = MDS_1, y = MDS_2)) +
      theme_few() +
      geom_point(aes(color = Geo_Location), size = 2) +
      scale_color_manual(name = "Country", values = kev_palette[1:length(unique(mds()[,4]))]) +
      geom_point(data = mds()[grep("Novel", mds()[,4]), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
      labs(x = "MDS 1", y = "MDS 2") +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.ticks.y = element_blank()) +
      theme(axis.text.y = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.y = element_text(size = 14, face = "bold")) +
      theme(axis.title.x = element_text(size = 14, face = "bold")) +
      theme(legend.title = element_text(size = 14, face = "bold")) +
      theme(legend.text = element_text(size = 12)) +
      theme(aspect.ratio = 0.6)
    })
  
  output$umap <- renderPlot ({
    ggplot(data = umap(), aes(x = UMAP_1, y = UMAP_2)) +
      theme_few() +
      geom_point(aes(color = Geo_Location), size = 2) +
      scale_color_manual(name = "Country", values = kev_palette[1:length(unique(umap()[,4]))]) +
      geom_point(data = umap()[grep("Novel", umap()[,4]), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.ticks.y = element_blank()) +
      theme(axis.text.y = element_blank()) +
      theme(axis.text.x = element_blank()) +
      theme(axis.title.y = element_text(size = 14, face = "bold")) +
      theme(axis.title.x = element_text(size = 14, face = "bold")) +
      theme(legend.title = element_text(size = 14, face = "bold")) +
      theme(legend.text = element_text(size = 12)) +
      theme(aspect.ratio = 0.6)
  })
  
  output$dend <- renderPlot ({
    par(mar = c(6,6,1,1))
    dend() %>% set("labels_cex", NA) %>% plot()
    legend("topright", legend = levels(factor(new_countries()[,2])), fill = kev_palette[1:length(unique(new_countries()[,2]))], pt.cex = 1, cex = 1, text.font = 2)
    colored_bars(dend_cols(), dend = dend(), rowLabels = c("Country"), cex.rowLabels = 1.25)
  })
  

}

shinyApp(ui = ui, server = server)




