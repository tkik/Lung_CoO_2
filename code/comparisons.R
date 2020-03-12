###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}

DATA = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/")
RESULT = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/output/")
CALLS = paste0(PATH_Y, "External/2018-10-Sotillo/data/methylDackel/")


comparisons <- data.frame(name=as.character(c("MsCCSP_tumor_control_comparison",
                                              "MsCCSP_tumor_rfp_control_comparison",
                                              "MsCCSP_tumor_rfp_tumor_comparison",
                                              "MsCCSP_MsSPC_comparison",
                                              "MsCCSP_MsSPC_control_comparison",
                                              "MsSPC_tumor_control_comparison",
                                              "MsSPC_tumor_rfp_control_comparison",
                                              "MsSPC_tumor_rfp_tumor_comparison",
                                              "MsCCSP_MsSPC_tumor_comparison")))

comparisons$cell_type1 <- c("MsCCSP",
                            "MsCCSP",
                            "MsCCSP",
                            "MsCCSP",
                            "MsCCSP",
                            "MsSPC",
                            "MsSPC",
                            "MsSPC",
                            "MsCCSP")

comparisons$cell_type2 <- c("MsCCSP",
                            "MsCCSP",
                            "MsCCSP",
                            "MsSPC",
                            "MsSPC",
                            "MsSPC",
                            "MsSPC",
                            "MsSPC",
                            "MsSPC")
comparisons$sample_type1 <- c("tumor",
                              "tumor-rfp",
                              "tumor-rfp",
                              NA,
                              "control",
                              "tumor",
                              "tumor-rfp",
                              "tumor-rfp",
                              "tumor")
comparisons$sample_type2 <- c("control",
                              "control",
                              "tumor",
                              NA,
                              "control",
                              "control",
                              "control",
                              "tumor",
                              "tumor")

saveRDS(comparisons, file = file.path(DATA, "comparisons.RDS"))
