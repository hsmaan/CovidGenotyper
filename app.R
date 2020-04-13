library(shiny)
library(shinythemes)
library(plyr)
library(shinycssloaders)
library(shinyWidgets)
library(Cairo)
library(intergraph)
library(ggnetwork)

options(repos = BiocManager::repositories())

# Source scripts

source("R/global.R")

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
vars_freq <- loadRData(grep("var_freq_sub*", file_list, value = TRUE))
date_file <- stringr::str_split_fixed((str_split_fixed(grep("dec_aligned_filtered", file_list, value = TRUE), "_", 4)[,4][1]), stringr::fixed("."), 2)[,1][1]

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

# RShiny 

ui <- fluidPage(theme = shinytheme("flatly"),
                
  titlePanel(title = "COVID-19 GENOTYPING TOOL (Alpha Testing)"),
      
  sidebarLayout(fluid = TRUE,
    
    tabsetPanel(
      
      tabPanel(strong("Input"), 
      
        sidebarPanel(
        
          h4(strong(strong("Please note that the application is under continuous development, but analyses are functional."))), 
          
          h4(strong("Instructions")), 
          
          p(
            "Add", strong("up to 10"), "complete Covid-19 sequences in fasta format, in one file. Please ensure fasta sequences have headers.", strong("Allow up to 5 minutes for processing."), "Hover over visualizations and use the plotly interface for navigation. Plots can be saved as a png using plotly, but for best quality we encourage users to save the webpage as a pdf and crop accordingly. Metadata can be toggled from below."
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
            "CGT is a tool to visualize public COVID-19 viral genome sequence information with respect to user uploaded sequences. Input fasta sequences are aligned to pre-aligned Covid-19 sequences uploaded to", a("GISAID.", href= "https://www.gisaid.org/"), "Each fasta sequence is assigned a name", strong("(Novel, number)"), "and is presented with respect to the public sequencing data. The following visualizations are done to determine sequence variation:", .noWS = c("after-begin", "before-end")
          ),
          
          p(
            strong("UMAP"), "- Uniform manifold approximation and projection",
            br(),
            strong("MST"), "- Minimum spanning tree of sequence network",
            br(),
            strong("SNP"), "- Prevalent non-synonymous coding single-nucleotide polymorphisms (Position - Protein - Variant Type)",
          ),
          
          p(
            "All methods approximate genomic differences using DNA distance determined by the Kimura-80 model of DNA evolution.",
          ),
          
          p(
            "Details on metadata:"
          ),
          
          p(
            strong("Region"), "- Major geographic region where sample was processed",
            br(),
            strong("Country"), "- Country where sample was processed",
            br(),
            strong("Collection date"), "- Sample collection date in terms of days since first case (Dec 1, 2019)",
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
           
           h4(
             strong("GISAID data:"),
           ),
           
           p(
             "COVID-19 viral genome sequences from GISAID are downloaded daily and processed using the CGT computational pipeline. No sequence information is published on the website, as per the", 
             a("GISAID data usage policy.", href = "https://www.gisaid.org/registration/terms-of-use"),
             "We thank all of the GISAID contributers for sharing their data. Full acknowledgements of COVID-19 sequence resources are available",
             a("here.", href = "404")
           ),
           
           h4(
             strong("Nextstrain metadata:")
           ),
           
           p(
             "Metadata for GISAID sequences from",
             a("Nextstrain's ncov repository", href = "https://github.com/nextstrain/ncov/tree/master/data"),
             "is utilized in conjunction with the up-to-date GISAID data. All credit goes to the Nextstrain team for curating this data."
           ),
           
           h4(
             strong("User privacy:")
           ),
           
           p(
             "CGT does not perform any server-side storage of user uploaded sequence data. Visualizations of user analyzed data are downloadable as png images. CGT simply processes user sequence data to create the visualizations using the R-Shiny reactive framework."
           ),
           
           p(
             strong(paste("Last data update:", date_file, sep = " "))
           ),

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
      
      navbarPage(title = NULL, tabPanel("UMAP", withSpinner(plotlyOutput("umap", height = "550px", width = "950px"), color = getOption("spinner.color", default = "#18BC9C"))), tabPanel("MST", withSpinner(plotlyOutput("mst", height = "550px", width = "950px"), color = getOption("spinner.color", default = "#18BC9C"))), tabPanel("SNP", withSpinner(plotlyOutput("snps", height = "550px", width = "950px"), color = getOption("spinner.color", default = "#18BC9C"))))

    
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
      snp_dfs <- snps_get(align(), meta_reac(), vars_freq)
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
    
    output$umap <- renderPlotly ({
      ggplotly(umap_plots()[[as.numeric(input$metatype)]]) %>% 
        layout(legend = list(font = list(size = 15)))
    })
  
    output$mst <- renderPlotly ({
      ggplotly(mst_plots()[[as.numeric(input$metatype)]]) %>%
        layout(legend = list(font = list(size = 15)))
    })
    
    output$snps <- renderPlotly ({
      ggplotly(snp_plots()[[as.numeric(input$metatype)]]) %>%
        layout(legend = list(font = list(size = 15)))
    })
    
  })
  
}
  
shinyApp(ui = ui, server = server)




