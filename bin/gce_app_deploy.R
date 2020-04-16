library(googleComputeEngineR)

# Load previous vm instance

setwd("../data")

files_list <- list.files()

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

vm_file <- grep("shiny_vm_instance*", files_list, value = TRUE)

vm <- loadRData(vm_file)

# Close previous instance and delete files

try(gce_vm_delete(vm), silent = TRUE) # Throws an error but instance deletes regardless

file.remove(vm_file)

rm(vm)

# Wait for instance to be deleted

Sys.sleep(300)

# Instantiate new vm instance of shiny app
      
vm <- gce_vm("cv", 
             template = "shiny",
             predefined_type = "n1-standard-16",
             externalIP = "34.74.33.29",
             dynamic_image = gce_tag_container("cgt:latest", "cgt-100"))

save(vm, file = "shiny_vm_instance.RData")
  
