library(methrix)


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

DMRs <- readRDS(file.path(DATA, "DMRs_noMsCCSP_control01_smoothed.RDS"))
comparisons_all <- names(DMRs)

labels <- gsub("_comparison", "", comparisons_all)
labels <- gsub("tumor_rfp", "tumor-rfp", labels)
labels <- do.call(rbind.data.frame, strsplit(labels, "_", fixed=T))
colnames(labels) <- c("cell", "first", "sec")
labels <- apply(labels, 2, as.character)

##correct the labels for 2 comparisons
labels[4,] <- c("", labels[4,1], labels[4,2])
labels[5,] <- c("", paste0(labels[5,1], "_control"), paste0(labels[5,2], "_control"))
labels[9,] <- c("", paste0(labels[9,1], "_tumor"), paste0(labels[9,2], "_tumor"))

for (i in 1:nrow(labels)){
  labels[i,c(2,3)] <- sort(c(labels[i,2], labels[i,3]), decreasing = F)
}
labels <- as.data.frame(labels)
labels$comparisons <- comparisons_all

saveRDS(labels, file.path(DATA, "labels.RDS"))
