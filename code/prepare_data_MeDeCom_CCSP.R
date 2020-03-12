library(methrix)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(ggplot2)
library(rtracklayer)
library(limma)
library(knitr)
library(pheatmap)
library(MeDeCom)
library(ggsci)

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
DOCS = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/docs/")

res <- readRDS(paste0(DATA, "no_snps_methrix.RDS"))

res_SSPC <- methrix::subset_methrix(res,
                                    samples = attr(res@colData, "rownames")[grep( "MsSPC_control", attr(res@colData, "listData")$full_name)])

mat_SSPC <- as.data.frame(get_matrix(res_SSPC, add_loci = T))
rownames(mat_SSPC) <- paste0(mat_SSPC$chr, "_", mat_SSPC$start)
mat_SSPC <- mat_SSPC[,-(1:3)]

mat_CCSP <-readRDS( file = file.path(PATH_Y, "External/Sotillo_mouse_lung_cell/data/MeDeCom_MCCSP_prepared_filt.rds"))

mat_SSPC <- mat_SSPC[rownames(mat_CCSP),]
mat <- as.matrix(cbind(mat_CCSP, mat_SSPC))
mat <- mat[complete.cases(mat),]

saveRDS(mat, file = file.path(DATA, "MeDeCom_MCCSP_SPC_prepared_filt.rds"))


