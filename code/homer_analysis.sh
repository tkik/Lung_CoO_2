module load anaconda3/2019.07
source activate homer

cd /icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer/bin/
bed_dir=/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/beds
output_dir=/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer_res/
beds=`ls -l $bed_dir | awk '{print $9}'`

perl ../configureHomer.pl -install mm10

for bed in $beds
do
comp="${bed%.*}"
echo $comp
#mkdir $output_dir/$comp
/icgc/dkfzlsdf/analysis/C010/lung_mouse_cells/methylation_analysis/homer_analysis/homer/bin/findMotifsGenome.pl $bed_dir/$bed mm10  $output_dir/$comp  -len 8,10,12 -size 100 -S 8 -p 8 -cache 6921 -fdr 0
done

