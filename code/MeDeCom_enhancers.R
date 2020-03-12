
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

mat  <- readRDS(paste0(DATA, "MsCCPS_MsSPC_enhancers.rds"))

medecom.result <- MeDeCom::runMeDeCom(as.matrix(mat), 2:5, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
saveRDS(medecom.result, file = paste0(DATA, "medecom_MsCCSP_MsSPC_enhancers.RDS"))

MeDeCom::plotParameters(MeDeComSet = medecom.result)

proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.1)
colnames(proportions) <- colnames(mat)

p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12)


mat  <- readRDS(paste0(DATA, "MsCCPS_MsSPC_promoters.rds"))


medecom.result <- MeDeCom::runMeDeCom(as.matrix(mat), 2:5, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
saveRDS(medecom.result, file = paste0(DATA, "medecom_MsCCSP_MsSPC_promoters.RDS"))

MeDeCom::plotParameters(MeDeComSet = medecom.result)

proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.1)
colnames(proportions) <- colnames(mat)

p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12)



mat  <- readRDS(paste0(DATA, "most_variable_in_controls.rds"))
mat <- mat [complete.cases(mat),]

medecom.result <- MeDeCom::runMeDeCom(as.matrix(mat), 2:5, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
saveRDS(medecom.result, file = paste0(DATA, "medecom_most_variable_in_controls.RDS"))

MeDeCom::plotParameters(MeDeComSet = medecom.result)

proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.1)
colnames(proportions) <- colnames(mat)

p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12)

mat  <- readRDS(paste0(DATA, "random_sites_all.rds"))
mat <- mat [complete.cases(mat),]

medecom.result <- MeDeCom::runMeDeCom(as.matrix(mat), 2:6, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
saveRDS(medecom.result, file = paste0(DATA, "medecom_random_sites_all.RDS"))

MeDeCom::plotParameters(MeDeComSet = medecom.result)

proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.01)
colnames(proportions) <- colnames(mat)

annotation <- data.frame(staining=gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\1", colnames(mat)),
                         type= gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\2", colnames(mat)),
                         origin = ifelse(grepl("rfp", colnames(mat)), "RFP", "GFP"))
rownames(annotation) <- colnames(mat)


p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12)

mat  <- readRDS(paste0(DATA, "random_sites.rds"))
mat <- mat [complete.cases(mat),]

medecom.result <- MeDeCom::runMeDeCom(as.matrix(mat), 2:8, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
saveRDS(medecom.result, file = paste0(DATA, "medecom_random_sites.RDS"))
#medecom.result <- readRDS( file = paste0(DATA, "medecom_random_sites.RDS"))

MeDeCom::plotParameters(MeDeComSet = medecom.result)

proportions <- MeDeCom::getProportions(medecom.result, K=5, lambda=0.01)
colnames(proportions) <- colnames(mat)

annotation <- data.frame(staining=gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\1", colnames(mat)),
                         type= gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\2", colnames(mat)),
                         origin = ifelse(grepl("rfp", colnames(mat)), "RFP", "GFP"))
rownames(annotation) <- colnames(mat)


p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12)

