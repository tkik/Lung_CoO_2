conda activate gatk

# GATK segmentation 
#get segmentation in R
#library(BSgenome.Mmusculus.UCSC.mm10)
#seq <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
#gr <- GRanges(
#  seqnames =names(seq),
#  ranges = IRanges(start =  1,
#                   end =  seq
#  )
#)
#gr<- gr[seqnames(gr) %in% paste0( "chr",1:19)]
#gr_tiled <- unlist(tile(gr, width=1000))
#export.bed(gr_tiled, "c010-datasets/External/2018-10-Sotillo/cnv_calling/mm10_100bp_binning.bed")
#bed to interval list
picard BedToIntervalList \
      I=c010-datasets/External/2018-10-Sotillo/cnv_calling/mm10_100bp_binning.bed \
      O=c010-datasets/External/2018-10-Sotillo/cnv_calling/mm10list.interval_list \
      SD=ngs_share/bcbio/genomes/Mmusculus/mm10/seq/mm10.dict

## Directories
OUTPUT_DIR=/C010-datasets/External/2018-10-Sotillo/cnv_calling/gatk/readCounts/
BAM_DIR=/C010-datasets/External/2018-10-Sotillo/bam_folder/
## Generate dict file (only if you dont have it already). 
## It creates a sequence dictionary file from a fasta file
#samtools dict -a GRCh37 -s Homo_sapiens -o $REFERENCE_DIR/hs37d5.dict $REFERENCE_DIR/hs37d5.fa
#Üsamtools dict -a GRCh37 -s Homo_sapiens -o $REFERENCE_DIR/hs37d5.chr7.dict $REFERENCE_DIR/hs37d5/hs37d5.chr.7.fa
#ngs_share/bcbio/genomes/Mmusculus/mm10/seq/mm10.dict

## Remove unwanted contigs from the genome build
#cat $REFERENCE_DIR/hs37d5_PhiX_1kb.interval_list | grep -v "^GL" | grep -v "^NC" | grep -v "^hs" | grep -v "^phiX" | grep -v "^MT" > $REFERENCE_DIR/hs37d5_PhiX_1kb.interval_list1 
#mv $REFERENCE_DIR/hs37d5_PhiX_1kb.interval_list1 $REFERENCE_DIR/hs37d5_PhiX_1kb.interval_list


## Loop over samples and the folders
PID=`ls -l $BAM_DIR | grep -v '.bai' | awk '{print $9}'`
echo $PID


#control01_B220_MsSPC_merged.mdup.bam 
#control02_B220_MsCCSP_merged.mdup.bam 
#control02_B220_MsSPC_merged.mdup.bam 
#control03_B220_MsCCSP_merged.mdup.bam 
#control03_B220_MsSPC_merged.mdup.bam 
#control04_B220_MsCCSP_merged.mdup.bam 
#tumor01_B220_MsCCSP_merged.mdup.bam 
#tumor01_B220_MsSPC_merged.mdup.bam 
#tumor02_B220_MsCCSP_merged.mdup.bam 
#tumor02_B220_MsSPC_merged.mdup.bam 
#tumor03_B220_MsCCSP_merged.mdup.bam 
#tumor03_B220_MsSPC_merged.mdup.bam



#for patient in $PID;
#do
#samtools index $BAM_DIR/${patient}
#done


for patient in $PID;
do
echo $patient;
## collect read counts for all samples
gatk CollectReadCounts -I $BAM_DIR/${patient} \
-L /C010-datasets/External/2018-10-Sotillo/cnv_calling/mm10list.interval_list \
--interval-merging-rule OVERLAPPING_ONLY -O $OUTPUT_DIR/${patient}.counts.h5
done

## Create a panel for control sample. One can also put all normals together to create a panel!
gatk CreateReadCountPanelOfNormals \
-I $OUTPUT_DIR/control01_B220_MsSPC_merged.mdup.bam.counts.h5 \
-I $OUTPUT_DIR/control02_B220_MsSPC_merged.mdup.bam.counts.h5 \
-I $OUTPUT_DIR/control03_B220_MsSPC_merged.mdup.bam.counts.h5 \
-I $OUTPUT_DIR/control02_B220_MsCCSP_merged.mdup.bam.counts.h5 \
-I $OUTPUT_DIR/control03_B220_MsCCSP_merged.mdup.bam.counts.h5 \
--minimum-interval-median-percentile 5.0 --maximum-zeros-in-sample-percentage 20 \
-O $OUTPUT_DIR/NMF_panel.h5 

mkdir $OUTPUT_DIR/normal/
mv $OUTPUT_DIR/*NMF* $OUTPUT_DIR/normal/

## Standardize and denoise case read counts against the PoN with DenoiseReadCounts

mkdir $OUTPUT_DIR/denoised_copy/
mkdir $OUTPUT_DIR/stand_copy/

PID=`ls -l $OUTPUT_DIR | grep  '.h5' | awk '{print $9}'`
echo $PID

for patient in $PID;
do
echo $patient
gatk DenoiseReadCounts -I $OUTPUT_DIR/${patient} \
--count-panel-of-normals $OUTPUT_DIR/normal/NMF_panel.h5 \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy/${patient}_denoised.tsv 
done

## Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios.
mkdir $OUTPUT_DIR/plots/

for patient in $PID;
do
echo $patient
gatk PlotDenoisedCopyRatios \
--standardized-copy-ratios $OUTPUT_DIR/stand_copy/${patient}_stand.tsv \
--denoised-copy-ratios $OUTPUT_DIR/denoised_copy/${patient}_denoised.tsv \
--sequence-dictionary /ngs_share/bcbio/genomes/Mmusculus/mm10/seq/mm10.dict \
-O $OUTPUT_DIR/plots/ --output-prefix ${patient}
done

R_DIR=/C010-datasets/External/2018-10-Sotillo/cnv_calling/gatk/
#anlyse segmented geneomes and their cnv

for patient in $PID;
do
echo $patient
Rscript $R_DIR/segment_denoiseCR.R \
$OUTPUT_DIR/denoised_copy/$patient"_denoised.tsv"
done
#This R script currently works with an old version of optparse. One has to install it from the archives

#packageurl <- "http://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.3.2.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")