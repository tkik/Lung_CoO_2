
###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}

#get segmentation in R
library(BSgenome.Mmusculus.UCSC.mm10)
seq <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
gr <- GRanges(
  seqnames =names(seq),
  ranges = IRanges(start =  1,
                   end =  seq
  )
)
gr<- gr[seqnames(gr) %in% paste0( "chr",1:19)]
gr_tiled <- unlist(tile(gr, width=1000))
export.bed(gr_tiled, file.path(PATH_Y, "External/2018-10-Sotillo/cnv_calling/mm10_100bp_binning.bed"))

