library(data.table)
library(shiny)
library(shinythemes)
library(ggplot2)
library(ggthemes)
library(plyr)
library(RColorBrewer)
library(shinycssloaders)
library(shinyWidgets)
library(Cairo)

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
meta <- as.data.frame(pre_meta[,c("Accession", "Region", "Geo_Location", "Datetime")])
pre_umap <- loadRData(grep("umap_preloaded", file_list, value = TRUE))
pre_mst <- loadRData(grep("mst_preloaded", file_list, value = TRUE))
pre_snp <- loadRData(grep("snps_preloaded", file_list, value = TRUE))

setwd("..")

# Load color palettes    

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

source("global.R")

# RShiny 

ui <- fluidPage(theme = shinytheme("flatly"),
                
  titlePanel(title = strong("COVID-19 GENOTYPING TOOL (Alpha Testing)")),
      
  sidebarLayout(
    
    tabsetPanel(
      
      tabPanel(strong("Input"), 
      
        sidebarPanel(
        
          h4(strong(strong("Please note that the tools and UI are under development, but analyses are functional."))), 
          
          h4(strong("Instructions")), 
          
          p(
            "Add complete Covid-19 sequences in fasta format, in one file. Please ensure fasta sequences have headers.", strong("Allow up to 5 minutes for processing.")
            ),
          
          fileInput(inputId = "input_fasta", label = h4(strong("Upload fasta"))),
          
          radioButtons("metatype", label = h4(strong("Metadata")), choices = list("Region" = 1, "Country" = 2, "Collection date" = 3), selected = 1),
          
          p(
            a("Bo Wang Lab", href="https://wanglab.ml/"),
            br(),
            "University Health Network",
            br(),
            "Toronto, Canada"
          )
        )
      ),
    
    tabPanel(strong("Details"),
                 
       sidebarPanel(
                  
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
            "Further information and detailed documentation are available at", a("hsmaan/CovidGenotyper.", href="https://github.com/hsmaan/CovidGenotyper")
          ),
          
          p(
            strong("Author:"), "Hassaan Maan",
            a("(Github)", href = "https://github.com/hsmaan")
          ),
          p(
            a("Bo Wang Lab", href="https://wanglab.ml/"),
            br(),
            "University Health Network",
            br(),
            "Toronto, Canada"
          )
       )
    ),
    
    tabPanel(strong("Data"),
         
       sidebarPanel(
           
           h2(strong("To be completed")),
           
           p(
             a("Bo Wang Lab", href="https://wanglab.ml/"),
             br(),
             "University Health Network",
             br(),
             "Toronto, Canada"
           )
         
         )
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
  
  meta_reac <- reactive ({
    if (is.null(input$input_fasta)) {
      new_meta <- meta
      return(new_meta)
    } else {
      new_accessions <- rownames(dist_reac())[(nrow(meta)+1):nrow(dist_reac())]
      meta_new <- data.frame("Accession" = new_accessions, "Region" = paste("Novel", seq(1, length(new_accessions), 1)), "Geo_Location" = paste("Novel", seq(1, length(new_accessions), 1)), "Datetime" = rep((unclass(Sys.Date()) - unclass(as.Date("2019-12-01", format = "%Y-%m-%d"))), length(new_accessions)))
      new_meta <- rbind(meta, meta_new)
      return(new_meta)
    }
  })
  
  umap <- reactive ({
    if (is.null(input$input_fasta)) {
      umap_df <- pre_umap
      return(umap_df)
    } else {
      umap_df <- umap_process(dist_reac(), meta_reac())
      return(umap_df)
    }
  })
  
  graph_m <- reactive ({
    if (is.null(input$input_fasta)) { 
      g <- pre_mst
      return(g)
    } else {
      g <- mst_graph(dist_reac(), meta_reac())
      return(g)
    }
  })
  
  snps <- reactive ({
    if (is.null(input$input_fasta)) { 
      snp_dfs <- pre_snp
      return(snp_dfs)
    } else {
      snp_dfs <- snps_get(align(), meta_reac())
      return(snp_dfs)
    }
  })
  
  umap_plots <- reactive ({
    plots <- umap_plotter(umap())
    return(plots)
  })
  
  mst_plots <- reactive ({
    plots <- mst_plotter(graph_m(), meta_reac())
    return(plots)
  })
  
  snp_plots <- reactive ({
    plots <- snp_plotter(snps(), meta_reac())
    return(plots)
  })
  
  observe ({
    
    output$umap <- renderPlot ({
      umap_plots()[as.numeric(input$metatype)]
    })
  
    output$mst <- renderPlot ({
      mst_plots()[as.numeric(input$metatype)]
    })
    
    output$snps <- renderPlot ({
      snp_plots()[as.numeric(input$metatype)]
    })
    
  })
  
}
  
shinyApp(ui = ui, server = server)




