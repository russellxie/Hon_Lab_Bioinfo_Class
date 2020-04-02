#!/bin/tcsh

#SBATCH --job-name=Bin_generator                          # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=10:00:00                                   # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                         # standard output file name
#SBATCH --error=serialJob.%j.time                         # standard error output file name
#SBATCH --mail-user=username@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

module load subread
module load bedtools
module load samtools
module load python/3.6.1-2-anaconda

set DATA_DIR = /project/GCRB/Hon_lab/s160875/03.analysis/Mosaic-seq/CROP-DE-analysis_10X-66K_no_downsampling-CPM.hg38/meta_analysis/ENCODE_data_hg38w

cat enhancer_regions.hg38.bed | ./bin_generator.py > enhancer_bins.hg38.saf

#count the total reads of each file
rm bam_read_count.txt
foreach BAM_FILE (`cut -f1 $DATA_DIR/ENCODE_annot.txt`)
    samtools view $DATA_DIR/$BAM_FILE | wc -l >> bam_read_count.txt
end

paste $DATA_DIR/ENCODE_annot.txt bam_read_count.txt > ENCODE_annot_counts.txt

foreach marker(`cut -f2 $DATA_DIR/ENCODE_annot.txt | sort | uniq`)

    set BAM_1 = `grep $marker $DATA_DIR/ENCODE_annot.txt | cut -f1 | head -1`
    set BAM_2 = `grep $marker $DATA_DIR/ENCODE_annot.txt | cut -f1 | head -2 | tail -1`

    featureCounts -a enhancer_bins.hg38.saf\
		  -F SAF\
		  -o enhancer_bins.hg38.$marker.counts\
		  -O\
		  $DATA_DIR/$BAM_1\
		  $DATA_DIR/$BAM_2
	      
    tail -n+3 enhancer_bins.hg38.$marker.counts\
    | ./heatmap_mat_generator.py `grep $marker ./ENCODE_annot_counts.txt | cut -f3 | sed -e 'N;s/\n/\t/g' | perl -ne 'chomp; @_ = split(/\t/,$_); (print($_[0] + $_[1]).'\n');'`\
    > enhancer_bins.hg38.$marker.cpm_matrix.txt
end



