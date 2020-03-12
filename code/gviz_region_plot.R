gviz_region_plot <- function(beta, anno = anno, location, genome="mm10", gene_name){

itrack <- IdeogramTrack(genome=genome, chromosome=as.character(seqnames(location)))

plot_data <- beta
plot_data$rowmeans_MsCCSP_control <- rowMeans(as.data.frame(mcols(plot_data[,which(colnames(mcols(plot_data)) %in% rownames(anno)[anno$sample_name=="control" & anno$cell_type=="MsCCSP"])])), na.rm=T)
plot_data$rowmeans_MsCCSP_tumor <- rowMeans(as.data.frame(mcols(plot_data[,which(colnames(mcols(plot_data)) %in% rownames(anno)[anno$sample_name=="tumor" & anno$cell_type=="MsCCSP"])])), na.rm=T)
plot_data$rowmeans_MsSPC_control <- rowMeans(as.data.frame(mcols(plot_data[,which(colnames(mcols(plot_data)) %in% rownames(anno)[anno$sample_name=="control" & anno$cell_type=="MsSPC"])])), na.rm=T)
plot_data$rowmeans_MsSPC_tumor <- rowMeans(as.data.frame(mcols(plot_data[,which(colnames(mcols(plot_data)) %in% rownames(anno)[anno$sample_name=="tumor" & anno$cell_type=="MsSPC"])])), na.rm=T)



dTrack <- DataTrack(plot_data[,grep("rowmeans", colnames(mcols(plot_data)))], name="DNA \nmethylation")
dTrack_MsCCSP_control <- DataTrack( plot_data[,"rowmeans_MsCCSP_control"], name="MsCCSP \ncontrol")
dTrack_MsCCSP_tumor <- DataTrack( plot_data[,"rowmeans_MsCCSP_tumor"], name="MsCCSP \ntumor")
dTrack_MsSPC_control <- DataTrack( plot_data[,"rowmeans_MsSPC_control"], name="MsSPC \ncontrol")
dTrack_MsSPC_tumor <- DataTrack( plot_data[,"rowmeans_MsSPC_tumor"], name="MsSPC \ntumor")

displayPars(dTrack_MsCCSP_control) <- list(col.histogram=mypal[1], fill.histogram=mypal[1])
displayPars(dTrack_MsCCSP_tumor) <- list(col.histogram=mypal[2], fill.histogram=mypal[2])
displayPars(dTrack_MsSPC_control) <- list(col.histogram=mypal[3], fill.histogram=mypal[3])
displayPars(dTrack_MsSPC_tumor) <- list(col.histogram=mypal[4], fill.histogram=mypal[4])

ucscGenes <- UcscTrack(genome="mm10", table="ncbiRefSeq", track = 'NCBI RefSeq',
                       trackType="GeneRegionTrack", chromosome=as.character(seqnames(location)),
                       rstarts = "exonStarts", rends = "exonEnds", gene = "name",
                       symbol = 'name', transcript = "name", name="UCSC genes",
                       nstrand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Mm.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z
#table="cpgIslandExt", track = 'CpG Islands',
ucscCpG <- UcscTrack(genome="mm10", chromosome=as.character(seqnames(location)),
                     trackType="AnnotationTrack",track="cpgIslandExt", from=start(location),
                     to=end(location), start="chromStart", end="chromEnd", id="name",
                     shape="box", fill="grey", name="CpG islands")



cat("\n")
res1 <- try(plotTracks(list(itrack, ucscCpG, ucscGenes2, dTrack_MsCCSP_control,
                    dTrack_MsCCSP_tumor, dTrack_MsSPC_control, dTrack_MsSPC_tumor),
               sizes = c(0.5,0.75, 1.5, 1,1,1,1),
               #groups=paste(anno$cell_type, anno$sample_name, sep="_"),
               type=c("histogram"),
               legend=TRUE,
               collapseTracks=T,
               aggregateGroups=TRUE,
               aggregation="mean",
               na.rm=T,
               col=mypal[1:7],
               cex=0.7,
               from = start(location),
               to = end(location),
               transcriptAnnotation="symbol",
               main = paste0("Gene plot, ", gene_name)), silent = T)

groups <- gsub("rowmeans_", "", colnames(mcols(plot_data[,grep("rowmeans", colnames(mcols(plot_data)))])))
res2 <- try(plotTracks(list(itrack, ucscCpG, ucscGenes2, dTrack),
               sizes = c(0.5,0.75, 1.5, 4),
               groups=groups,
               type=c("p", "smooth"),
               legend=TRUE,
               na.rm=TRUE,
               col=mypal[1:7],
               cex=0.1,
               from = start(location),
               to = end(location),
               transcriptAnnotation="symbol",
               main = paste0("Gene plot, ", gene_name)), silent = T)
#dev.off()
return(list(res1, res2))

}
