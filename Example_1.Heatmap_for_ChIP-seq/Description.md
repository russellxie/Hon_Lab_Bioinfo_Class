# Example 1. Heatmap plotting for ChIP-seq / ATAC-seq / DNase-seq Samples

## Description
Genomic assays such as ChIP-seq, ATAC-seq and DNase-seq can robustly capture the epigenetic features across the whole genome. However, 

Think about this question:
> Give a list of regions, how to illustrate multiple epigenetic features of these regions in the same figure? 

## The Input File You Have:

- A Bed File containing the region information, for example:

```
	chr10	101364638	101365038
	chr10	16996176	16996576
	chr10	17003096	17003496
```
- The bam files for each sequencing libraries, in which each line represents one read.

## Hints:

- Ignore the programming part and think about how would you do it as a human being.
- This is what we expect to have:

<img src="Heatmap_allmarkers_200X50.png" alt="test image size" height="60%" width="60%">
