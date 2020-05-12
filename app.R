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
                
  tags$head(tags$link(rel = "shortcut icon", type = "image/png", href = "cgt_icon.png"),
            tags$title("Covid-19 Genotyping Tool")),
                
  titlePanel(title = img(height = 80, width = 80, src = "cgt_logo.png", align = "left")), h3(strong("COVID-19 GENOTYPING TOOL", align = "left"), img(src = "sunnybrook.png", height = 35, align = "right"), img(src = "mcmaster_logo.png", height = 60, align = "right"), img(src = "vector_logo.png", height = 60, align = "right"), img(src = "pmcc.jpg", height = 40, align = "right")),
  
  sidebarLayout(fluid = TRUE,
    
    tabsetPanel(
      
      tabPanel(strong("Input"), 
      
        sidebarPanel(
          
         
        
          p("The COVID-19 Genotyping Tool (CGT) is a visualization toolbox for SARS-CoV-2 whole genome sequencing data. Public sequences from GISAID is preloaded for inspection, and users have the option of uploading in-house SARS-CoV-2 sequencing data for concurrent analysis with public data."),
          
          h4(strong("Instructions")), 
          
          p(
            "Add", strong("up to 10"), "complete (>29000 bp) COVID-19 sequences in fasta format, in one file. Ensure fasta sequences have headers.", strong("Processing all visualizations should take up to 15 minutes, but may be longer if uploaded sequences are sufficiently dissimilar to public domain data."), "Please", strong("DO NOT REFRESH"), "the page after uploading. Uploaded data will be lost on refresh as we do not store or cache user data in any way.", "Hover over visualizations and use the plotly interface for navigation. Plots can be saved as a png using plotly. Metadata can be toggled from below."
            ),
          
          fileInput(inputId = "input_fasta", label = h4(strong("Upload fasta"))),
          
          radioButtons("metatype", label = h4(strong("Metadata")), choices = list("Region" = 1, "Country" = 2, "Collection date" = 3), selected = 1),
          
          p(
            a("Bo Wang Lab", href="https://wanglab.ml/"),
            br(),
            "University Health Network",
            br(),
            "Vector Institute",
            br(),
            "Toronto, Canada"
          )
        )
      ),
    
    tabPanel(strong("Details"),
                 
       sidebarPanel(
         
                  
          p(
            "CGT is a tool to visualize public COVID-19 viral genome sequence information with respect to user uploaded sequences. Input fasta sequences are aligned to pre-aligned COVID-19 sequences uploaded to", a("GISAID.", href= "https://www.gisaid.org/"), "Each fasta sequence is assigned a name", strong("(Novel, number)"), "and is presented with respect to the public sequencing data. The following visualizations are done to determine sequence variation:", .noWS = c("after-begin", "before-end"),
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
            "Further information and detailed documentation are available at", a("hsmaan/CovidGenotyper.", href="https://github.com/hsmaan/CovidGenotyper"), "Recommendations and problems with the application can be opened up as an issue on GitHub. General feedback can be forwarded to", a("hurmaan99@gmail.com", href = "mailto:hurmaan99@gmail.com")
          ),
          
          p(
            strong("Author/Maintainer:"), "Hassaan Maan",
            a("(Github)", href = "https://github.com/hsmaan")
          ),
          p(
            a("Bo Wang Lab", href="https://wanglab.ml/"),
            br(),
            "University Health Network",
            br(),
            "Vector Institute",
            br(),
            "Toronto, Canada"
          )
       )
    ),
    
    tabPanel(strong("Data"),
         
       sidebarPanel(
              
           h4(
             strong("GISAID data"),
           ),
           
           p(
             "COVID-19 viral genome sequences and sequence metadata from GISAID are downloaded and processed using the CGT computational pipeline, updated on a weekly basis. No sequence information is explicitly published on the website, as per the", 
             a("GISAID data usage policy.", href = "https://www.gisaid.org/registration/terms-of-use"),
             "We thank all of the GISAID contributers for sharing their data. Full up-to-date acknowledgements of COVID-19 sequence resources available", 
             a("here.", href = "https://github.com/hsmaan/CovidGenotyper/ack")
           ),
           
           h4(
             strong("User privacy")
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
             "Vector Institute",
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
      gc()
      return(align)
    }
  })
  
  dist_reac <- reactive ({
    if (is.null(input$input_fasta)) {
      dist <- pre_dist
      return(dist)
    } else {
      dist <- dist_get_heur(align(), input$input_fasta$datapath, pre_dist)
      gc()
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
      gc()
      return(new_meta)
    }
  })
  
  umap <- reactive ({
    if (is.null(input$input_fasta)) {
      umap_df <- pre_umap
      return(umap_df)
    } else {
      umap_df <- umap_process_heur(align(), input$input_fasta$datapath, dist_reac(), meta_reac(), pre_umap)
      gc()
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
      gc()
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




