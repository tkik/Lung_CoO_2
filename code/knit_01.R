###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="Z:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="Y:/"} else {
    PATH_Y="/C010-datasets/"}

DATA = paste0(PATH, "Reka/33_CoO_lung/data_02_10/")
RESULT = paste0(PATH, "Reka/33_CoO_lung/results_02_10/")
CALLS = paste0(PATH_Y, "External/2018-10-Sotillo/data/methylDackel/")
SNPS = paste0(PATH_Y, "External/2018-10-Sotillo/data/snps/")


rmarkdown::render(input = "01_read_in_MethylDackel.Rmd", output_dir = RESULT, 
                  knit_root_dir = RESULT)
