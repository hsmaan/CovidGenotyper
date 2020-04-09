library(googleComputeEngineR)

# Test vm instance of shiny app

vm <- gce_vm("cv", 
             template = "shiny",
             predefined_type = "n1-standard-2",
             dynamic_image = gce_tag_container("cgt", "mcarthurlabpoc"))


# Close 

gce_vm_stop(vm)

