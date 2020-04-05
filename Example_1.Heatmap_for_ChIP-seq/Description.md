# Example 1. Heatmap plotting for ChIP-seq / ATAC-seq / DNase-seq Samples

## Introduction

Genomic assays such as ChIP-seq, ATAC-seq and DNase-seq can robustly capture the epigenetic features across the whole genome. However, it is not that intuitive to illustrate data from large number of genomic loci. In this example, we will cover a common strategy to achieve this goal: plotting heatmap.

Think about this question:

> Give a list of regions, how to illustrate multiple epigenetic features of these regions in the same figure?

## The Input File You Have  
- A Bed File containing the region information, for example:  
> chr10 101364638 101365038  
> chr10 16996176 16996576  
> chr10 17003096 17003496  

- The bam files for each sequencing libraries, in which each line represents one read.

## Hints  
- Ignore the programming part and focus on what steps you need to take.
- This is what we expect to have:  

![drawing](Heatmap_allmarkers_200X50.png)

## File Description  
`log.sh` Main log file for job submission.  
`heatmap_generator.py` Main python script to generate the matrix required for heatmap plotting.  
`Heatmap_plotting.ipynb` Jupyter notebook file for plotting.  
`enhancer_regions.hg38.bed` Region info.
