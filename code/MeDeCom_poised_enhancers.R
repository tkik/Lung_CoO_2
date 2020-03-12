
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
library(TCA)

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
#
# mat  <- readRDS(paste0(DATA, "medecom/smoothed_bivalent_enhancers.RDS"))
# rowsds <- rowSds(mat, na.rm=T)
# mat <- mat[order(rowsds, decreasing = T)[1:100000],]
#
#
# #medecom.result <- MeDeCom::runMeDeCom(as.matrix(mat), 2:5, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
# #saveRDS(medecom.result, file = paste0(DATA, "medecom/smoothed_poised_enhancers.RDS"))
# medecom.result <- readRDS( file = paste0(DATA, "medecom/smoothed_poised_enhancers.RDS"))


mat  <- as.data.frame(readRDS(paste0(DATA, "medecom/bivalent_enhancers.RDS")))
rowsds <- rowSds(as.matrix(mat[,-(1:3)]), na.rm=T)
mat <- mat[order(rowsds, decreasing = T)[1:100000],]
mat <- mat[complete.cases(mat),]
mat_control <- mat[,grep("control", colnames(mat))]

#medecom.result2 <- MeDeCom::runMeDeCom(as.matrix(mat_control), 2:6, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
#saveRDS(medecom.result2, file = paste0(DATA, "medecom/poised_enhancers_controls.RDS"))

medecom.result <-  readRDS( file = paste0(DATA, "medecom/poised_enhancers_controls.RDS"))
MeDeCom::plotParameters(MeDeComSet = medecom.result)


K_sel <- 4
lambda_sel <- 0.0001
proportions <- MeDeCom::getProportions(medecom.result, K=K_sel, lambda=lambda_sel)
colnames(proportions) <- colnames(mat_control)

LMCs <- MeDeCom::getLMCs(medecom.result, K=K_sel, lambda=lambda_sel)
colnames(LMCs) <- paste0("LMC", 1:K_sel)
LMCs <- cbind(mat[,1:3], LMCs)



annotation <- data.frame(staining=gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\1", colnames(mat_control)),
                         type= gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\2", colnames(mat_control)),
                         origin = ifelse(grepl("rfp", colnames(mat_control)), "RFP", "GFP"))
rownames(annotation) <- colnames(mat_control)


p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12, annotation_col = annotation,
             file= file.path(RESULT, "poised_enhancers_controls_heatmap.pdf"))




which_mat <- which(medecom.result@parameters$Ks==K_sel)+(which(medecom.result@parameters$lambdas==lambda_sel)-1)*length(medecom.result@parameters$Ks)

all_Ts <- MeDeCom::factorize.regr(D= as.matrix(mat[,-(1:3)]), Tt = medecom.result@outputs$`1`$T[[which_mat]])

proportions_all <- all_Ts$A
colnames(proportions_all) <- colnames(mat[,-(1:3)])
rownames(proportions_all) <- rownames(proportions)

annotation <- data.frame(staining=gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\1", colnames(mat[,-(1:3)])),
                         type= gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\2", colnames(mat[,-(1:3)])),
                         origin = ifelse(grepl("rfp", colnames(mat[,-(1:3)])), "RFP", "GFP"))
rownames(annotation) <- colnames(mat[,-(1:3)])

annotation$ID <- colnames(mat[,-(1:3)])



p <- pheatmap(proportions_all,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12, annotation_col = annotation,
              filename = file.path(RESULT, "poised_enhancers_all_heatmap.pdf"))




sequencing_annotation_table <- read.csv(file.path(DATA, "sequencing_annotation_table.csv"))
sequencing_annotation_table$SampleID <-gsub("-", "", paste0(gsub("B220_", "", sequencing_annotation_table$Individual), "_",
                                              gsub("-", ".", tolower(sequencing_annotation_table$Sample.Type), fixed=T)))

annotation <- merge(annotation, sequencing_annotation_table[,c("SampleID", "File.exists")], by.x="ID", by.y="SampleID", sort=F, all=F)
annotation <- annotation[!duplicated(annotation$ID),]
annotation$type <- as.numeric(annotation$type)
annotation$File.exists <- as.numeric(annotation$File.exists)
rownames(annotation) <- annotation$ID
#annotation$uniform <- 1
results <- tca(as.matrix(mat[,-(1:3)]), t(as.matrix(proportions_all)), C1=as.matrix(annotation[,"type", drop=F]),
               C2=as.matrix(annotation[,"File.exists", drop=F]), refit_W = T)
#saveRDS(results, file.path(RESULT, "tca_res.RDS"))
results <- readRDS(file.path(RESULT, "tca_res.RDS"))

res_cal <- results$mus_hat
colnames(res_cal) <- paste0("calc_", colnames(res_cal))
combined <- cbind(LMCs[,-(1:3)], res_cal)
pheatmap(combined[1:10000,])



##Promoters
candidate_promoters <- readRDS(file.path(DATA, "candidate_promoter_methylation.RDS"))
candidate_promoters <- candidate_promoters[complete.cases(candidate_promoters),]
rownames(candidate_promoters) <- 1:nrow(candidate_promoters)


results <- tca(as.matrix(candidate_promoters), t(as.matrix(proportions_all)), C1=as.matrix(annotation[,"type", drop=F]),
              C2=as.matrix(annotation[,"File.exists", drop=F]), refit_W = T)

saveRDS(results, file.path(RESULT, "tca_res_candidate_promoters.RDS"))
results <- readRDS(file.path(RESULT, "tca_res_candidate_promoters.RDS"))

res_cal <- results$mus_hat
colnames(res_cal) <- paste0("calc_", colnames(res_cal))
combined <- cbind(LMCs[,-(1:3)], res_cal)
pheatmap(combined[1:10000,])
