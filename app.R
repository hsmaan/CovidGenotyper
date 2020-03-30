library(data.table)
library(shiny)
library(shinythemes)
library(ggplot2)
library(ggthemes)
library(plyr)
library(RColorBrewer)
library(shinycssloaders)
library(shinyWidgets)

# Load testing data 

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

setwd("data")

file_list <- list.files() 

pre_aligned_filtered <- loadRData(grep("dec_aligned_filtered", file_list, value = TRUE))
pre_meta <- loadRData(grep("covid_filtered_meta", file_list, value = TRUE))
pre_dist <- loadRData(grep("dec_fasta_dist", file_list, value = TRUE))
countries <- as.data.frame(pre_meta[,c("Accession", "Region")])

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

qual_palettes = brewer.pal.info[brewer.pal.info$category == "qual", ]
qual_vector = unlist(mapply(brewer.pal, qual_palettes$maxcolors, rownames(qual_palettes)))

# Source scripts

source("bin/align.R")
source("bin/clustal_dist.R")
source("bin/umap_get.R")
source("bin/mst_graph.R")
source("bin/snps.R")

# RShiny 

ui <- fluidPage(theme = shinytheme("flatly"),
                
  titlePanel(title = strong("COVID-19 GENOTYPING TOOL (Alpha Testing)")),

  sidebarLayout(
    
    sidebarPanel(
      
      h4(strong("Please note that the tools and UI are under development, but analyses are functional.")), 
      
      h4(strong("Instructions")), 
      
      p(
        "Add complete Covid-19 sequences in fasta format, in one file. Please ensure fasta sequences have headers. Files are not stored in any way."
      ),
      
      fileInput(inputId = "input_fasta", label = "Upload fasta"),
      
      h4(strong("Details")),
      
      p(
        "Input fasta sequences are aligned to pre-aligned Covid-19 sequences uploaded to", a("GISAID.", href= "https://www.gisaid.org/"), "Each fasta sequence is assigned a name", strong("(Novel, number)"), "and is presented with respect to the public sequencing data. The following visualizations are done to determine sequence variation:", .noWS = c("after-begin", "before-end")
      ),
      
      p(
        strong("UMAP"), "- Uniform manifold approximation and projection",
        br(),
        strong("MST"), "- Minimum spanning tree of sequence network",
        br(),
        strong("SNP"), "- Prevalent single-nucleotide polymorphisms",
      ),
      
      p(
        "All methods approximate genomic differences using DNA distance determined by the Kimura-80 model of DNA evolution."
      ),
        
      p(
        a("Bo Wang Lab", href="https://wanglab.ml/"),
        br(),
        "University Health Network",
        br(),
        "Toronto, Canada"
      )
      
    ),
  
    mainPanel(
      
      navbarPage(title = NULL, tabPanel("UMAP", withSpinner(plotOutput("umap", height = "550px", width = "950px"), color = getOption("spinner.color", default = "#18BC9C"))), tabPanel("MST", withSpinner(plotOutput("mst", height = "550px", width = "950px"), color = getOption("spinner.color", default = "#18BC9C"))), tabPanel("SNP", withSpinner(plotOutput("snps", height = "550px", width = "950px"), color = getOption("spinner.color", default = "#18BC9C"))))

    
  )
  )
                  
)
  
server <- function(input, output) {
  
  align <- reactive ({
    if (is.null(input$input_fasta)) {
      align <- pre_aligned_filtered
      return(align)
    } else {
      align <- align_get(input$input_fasta$datapath, pre_aligned_filtered)
      return(align)
    }
  })
  
  dist_reac <- reactive ({
    if (is.null(input$input_fasta)) {
      dist <- pre_dist
      return(dist)
    } else {
      dist <- dist_get(align())
      return(dist)
    }
  })
  
  new_countries <- reactive ({
    if (is.null(input$input_fasta)) {
      return(countries)
    } else {
      new_accessions <- rownames(dist_reac())[(nrow(countries)+1):nrow(dist_reac())]
      countries_new <- data.frame("Accession" = new_accessions, "Region" = paste("Novel", seq(1, length(new_accessions), 1)))
      new_countries <- rbind(countries, countries_new)
      return(new_countries)
    }
  })
  
  snps <- reactive ({
    snp_df <- snps_get(align(), new_countries())
    return(snp_df)
  })
  
  umap <- reactive ({
    umap_res <- umap_process(dist_reac())
    acc_names = rownames(dist_reac())
    umap_df <- data.frame("Accession" = acc_names, "UMAP_1" = umap_res[,1], "UMAP_2" = umap_res[,2])
    umap_df <- merge(umap_df, new_countries())
    return(umap_df)
  })
  
  graph_m <- reactive ({
    g <- mst_graph(dist_reac(), new_countries(), kev_palette)
    return(g)
  })
  
  observe ({
  
    output$umap <- renderPlot ({
      ggplot(data = umap(), aes(x = UMAP_1, y = UMAP_2)) +
        theme_few() +
        geom_jitter(aes(fill = Region), size = 3, position = "jitter", colour = "black", pch = 21, stroke = 0.25) +
        scale_fill_manual(name = "Region", values = kev_palette[1:length(unique(umap()[,4]))]) +
        geom_point(data = umap()[grep("Novel", umap()[,4]), ], pch = 21, fill = NA, size = 4, colour = "firebrick1", stroke = 4) +
        labs(x = "UMAP 1", y = "UMAP 2") +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.ticks.y = element_blank()) +
        theme(axis.text.y = element_blank()) +
        theme(axis.text.x = element_blank()) +
        theme(axis.title.y = element_text(size = 16, face = "bold")) +
        theme(axis.title.x = element_text(size = 16, face = "bold")) +
        theme(legend.title = element_text(size = 15, face = "bold")) +
        theme(legend.text = element_text(size = 14)) +
        theme(aspect.ratio = 0.6)
    })
  
    output$mst <- renderPlot ({
      lay <- layout_with_graphopt(graph_m(), niter = 1000)
      plot.igraph(graph_m(), vertex.label = NA, vertex.size = 4, edge.width = 1, layout = lay, edge.color = "gray25")
      legend("topleft", legend = levels(factor(new_countries()[,2])), fill = kev_palette[1:length(unique(new_countries()[,2]))], pt.cex = 1, cex = 1, text.font = 2)
    })
    
    output$snps <- renderPlot ({
      ggplot(data = snps(), aes(x = Allele, y = Freq)) +
        theme_few () +
        geom_bar(stat = "identity", position = "dodge2", aes(fill = Meta), color = "black") +
        scale_fill_manual(name = "Region", values = kev_palette[1:length(unique(new_countries()[,2]))]) +
        facet_wrap(~Position, scales = "free") +
        labs(x = "Allele", y = "Frequency") + 
        theme(axis.text.y = element_text(size = 12)) +
        theme(axis.text.x = element_text(size = 12)) +
        theme(axis.title.y = element_text(size = 16, face = "bold")) +
        theme(axis.title.x = element_text(size = 16, face = "bold")) +
        theme(legend.title = element_text(size = 15, face = "bold")) +
        theme(legend.text = element_text(size = 14)) +
        theme(strip.text = element_text(size = 14, face = "bold"))
    })
    
  })
  
}
  
shinyApp(ui = ui, server = server)




