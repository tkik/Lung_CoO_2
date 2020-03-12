
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
library(annotatr)
library(ggfortify)


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


remove_correlated <- function(mat, threshold=0.95, iter=100){
  mat_removed <- list()

  #lscomb = list()
  #rmpt = list()
  #size = 30*10^3

  rmpi = NULL

  pb = txtProgressBar(
    min = 1,
    max = iter,
    initial = 1,
    style = 3
  )

  for (j in 1:iter) {
    setTxtProgressBar(pb, j)
    set.seed(j)
    random.site = sample(1:nrow(mat), size = 1)

    retain.sites = unlist(parallel::mclapply(seq_along(1:nrow(mat)), function(i) {
      #  lmat <- as.data.frame(t(mat[c(random.site, i)]))

      return(cor.test(as.numeric(mat[random.site, ]), as.numeric(mat[i, ]))$estimate)

    }, mc.cores = 1))

    remove.sites = which(retain.sites > threshold)
    remove.sites = remove.sites[-which(remove.sites==random.site)]
    #browser()
    if (length(remove.sites)!=0){
      if (!rownames(mat)[random.site] %in% names(mat_removed)){
        mat_removed[[rownames(mat)[random.site]]] <- mat[remove.sites,]
      } else {
        mat_removed[[rownames(mat)[random.site]]] <- rbind(mat_removed[[rownames(mat)[random.site]]], mat[remove.sites,])
      }

      mat = mat[-remove.sites, ]}
    rmpi = c(rmpi, length(remove.sites))

    if (all(tail(rmpi)==0))
      stop("finished with the selection")
  }
  return(mat)
}

res <- readRDS(paste0(DATA, "no_snps_methrix.RDS"))
res <- methrix::subset_methrix(res, samples = attr(res@colData, "rownames")[-which(attr(res@colData, "listData")$full_name=="MsCCSP_control01")])

annotation <- read.delim(file=file.path(DATA, "annotation_table_TAGWGBS_02_10.txt"),stringsAsFactors = F)

########different filtering techniques for MeDeCom

# 1. select the most variable sites in controls (check with tumors)

m_2 <- subset_methrix(res, samples = rownames(res@colData)[res@colData$sample_name %in% c("control")])

m_2 <- coverage_filter(m_2, cov_thr = 8, min_samples = 5)

beta <- as.data.frame(get_matrix(m_2, add_loci = T))
rownames(beta) <- paste0(beta$chr, "_", beta$start)
beta <- beta[,-(1:3)]

mat <- beta[order(matrixStats::rowSds(as.matrix(beta), na.rm=T), decreasing = T)[1:500000],]

mat <- remove_correlated(mat)

saveRDS(mat, file= file.path(DATA, "most_variable_in_controls.rds"))
# 2. select the most variable sites in tumors

# random selection of 100000 sites.
beta <- get_matrix(m_2, add_loci = T, in_granges = T)
m_2_2 <- subset_methrix(res, beta)

beta <- as.data.frame(get_matrix(m_2_2, add_loci = T))
rownames(beta) <- paste0(beta$chr, "_", beta$start)
beta <- beta[,-(1:3)]
beta <- beta[complete.cases(beta),]

for(i in c(1,5363,87,3,83,546387, 534)){
set.seed(i)
mat <- beta[sample(1:nrow(beta), 100000, replace=F),]

pr <- prcomp(t(mat))
p <- autoplot(pr, data = as.data.frame(m_2_2@colData), colour="sample_name", shape="cell_type")
print(p)
}

set.seed(83)
mat <- beta[sample(1:nrow(beta), 100000, replace=F),]

saveRDS(mat, file= file.path(DATA, "random_sites_all.rds"))


beta <- get_matrix(m_2, add_loci = T)
rownames(beta) <- paste0(beta$chr, "_", beta$start)
beta <- beta[,-(1:3)]
beta <- beta[complete.cases(beta),]

for(i in c(1,5363,87,3,83,546387, 534)){
  set.seed(i)
  mat <- beta[sample(1:nrow(beta), 100000, replace=F),]

  pr <- prcomp(t(mat))
  p <- autoplot(pr, data = as.data.frame(m_2@colData), colour="sample_name", shape="cell_type")
  print(p)
}

set.seed(83)
mat <- beta[sample(1:nrow(beta), 100000, replace=F),]

saveRDS(mat, file= file.path(DATA, "random_sites.rds"))
