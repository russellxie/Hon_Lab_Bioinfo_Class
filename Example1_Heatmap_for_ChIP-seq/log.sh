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

set DATA_DIR = /project/GCRB/Hon_lab/s160875/03.analysis/Mosaic-seq/CROP-DE-analysis_10X-66K_no_downsampling-CPM.hg38/meta_analysis/ENCODE_data_hg38

foreach SAMPLE (`cut -f2 ENCODE_annot.txt`)
    echo $SAMPLE
    set BAM_FILE = `grep $SAMPLE ENCODE_annot.txt | cut -f1`
    echo $BAM_FILE
    ./heatmap_generator.py --region enhancer_regions.hg38.bed\
                        --out_dir ENCODE_heatmap\
                        --sample_name $SAMPLE\
                        --bam $DATA_DIR/$BAM_FILE\
                        --bin_num 200\
                        --bin_size 50\

end



