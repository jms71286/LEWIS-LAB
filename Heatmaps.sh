After FiberTools_Pt3 - make gtf extract of genes of interest

#NOT SCRIPT - run line by line in terminal.

make a new .bed file from a list of genes in .txt using command line in terminal aka not in a script

First CD into the folder where your genome files are

grep -Ff my_gene_list.txt reference.gtf > extracted_genes.gtf
awk '$3=="gene" {print $1 "\t" $4-1 "\t" $5 "\t" $9}' extracted_genes.gtf > extracted_genes.bed

line1: grep -Ff AFUB_conidia.txt GCA_051225625.1_liftoff.gff3 > extractedgenes_conidia.gtf
Line2: awk '$3=="gene" {print $1 "\t" $4-1 "\t" $5 "\t" $9}' extractedgenes_conidia.gtf > extractedgenes_conidia.bed


Change my_gene_list to the name of your txt file; the name of the reference should be .gff/gtf
Change extracted_genes.gtf to my_gene_list.gff/gtf
Change extracted_genes.bed to my_gene_list.bed
Nano my_gene_list.bed  to see that it made the 6 columns with the genes of interest

DO LINE BY LINE:

ml deepTools/3.5.5-gfbf-2023a

#example:
#computeMatrix reference-point -p 12 -R "/home/jms71286/research/genomes/CEA10/GCA_051225625.1_TSS_CEA10.bed" -S ../BigWigs/157-43_ChIP_tetRsuz12_36hr_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw ../BigWigs/157-47_ChIP_tetRsuz12_36hr_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw ../BigWigs/157-53_ChIP_epr1_tetRsuz12_36hr_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw ../BigWigs/157-57_ChIP_epr1_tetRsuz12_36hr_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw	-o EPR1_K27Genes_K27me3_K36me3_V1.mat 	--sortRegions keep --

computeMatrix reference-point -p 12 -R "/home/jms71286/research/genomes/CEA10/GCA_051225625.1_TSS_CEA10.bed" -S /lustre2/scratch/jms71286/FIBER/ONTRun13/fibertools_results/ONTRun13/ONTRun13_CEA10_conidia3_barcode09_merged/ONTRun13_CEA10_conidia3_barcode09_merged.m6Apileup.bw /lustre2/scratch/jms71286/FIBER/ONTRun13/fibertools_results/ONTRun13/ONTRun13_CEA10_conidia2_barcode08_merged/ONTRun13_CEA10_conidia2_barcode08_merged.m6Apileup.bw /lustre2/scratch/jms71286/FIBER/ONTRun13/fibertools_results/ONTRun13/ONTRun13_CEA10_conidia1_barcode07_merged/ONTRun13_CEA10_conidia1_barcode07_merged.m6Apileup.bw	-o /lustre2/scratch/jms71286/FIBER/ONTRun13/HeatmapWGenome/CEA10m6A_HM.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 1500 -b 1000

plotHeatmap -m /lustre2/scratch/jms71286/FIBER/ONTRun13/HeatmapWGenome/CEA10m6A_HM.mat -o /lustre2/scratch/jms71286/FIBER/ONTRun13/HeatmapWGenome/CEA10m6A_HM.png --sortRegions descend  --sortUsingSamples 1 2  --outFileSortedRegions /lustre2/scratch/jms71286/FIBER/ONTRun13/HeatmapWGenome/CEA10m6A_HM.bed --startLabel "5'" -max 15 8 15 8 --endLabel "3'"  --boxAroundHeatmaps no  --samplesLabel "m6A_R1" "m6A_R2" "m6A_R3" --colorMap 'Greens' 'Greens' 'Greens' --sortUsing mean --heatmapHeight 16  --heatmapWidth 3 --yMin 2 --yMax 8

#TO VIEW IN IGV:
Download GCA_051225625.1_TSS_CEA10.bed to desktop
Download heatmap.bed to desktop (from /scratch/jms71286/FIBER/ONTRun13/HeatmapWGenome/CEA10m6A_HM.bed)
Download _merged.m6Apileup.bw (from analysis) to desktop
Download GCA_051225625.1_liftoff.gff3 to desktop
upload files in IGV
