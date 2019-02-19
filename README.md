Imputing IBD Immuno Chip Genetic Data
================
Translational Genomics Group
2/22/2019

Organization
------------

Each workflow step has its own folder and README. Often each workflow step is repeated for each of our identified cohorts. Genetic data and derivatives are not updated to this repository, but are available.

Original Data
-------------

1.  Genotype data from IBDichip1to7eBBCTOP\_unfilteredHG19b.\*
2.  Updated annotation file 25\_jul\_2018
3.  List of batches each sample was run under for cohort splits (e.g. ichip 1 vs 7, from PCA\_Investigation Project)

Workflow
--------

1.  Split genetic data into relevant cohorts [split\_cohorts](workflow/1.split_cohorts/)
2.  Estimate ancestry and split cohorts based on ancestry [estimate\_ancestry](workflow/2.estimate_ancestry/)
3.  Perform ancestry specific QC [qc](wokflow/3.qc/)
4.  Prepare genotype data for Imputation QC [pre\_impute\_prep](workflow/4.pre_impute_prep/)
5.  Post imputation QC [post\_impute\_qc](workflow/5.post_impute_qc/)

Summary
-------

All samples were gentoyped using the Illumina ImmunoChip platform and markers were assigned hg19 coordinates based on the manufacturers alignment files. Markers mapping to insertions and deletions, non-automosomal chromosomes, or markers that were were unable to be mapped to an hg19 coordinate were removed. Samples were split into three cohorts based on the version of the platform they were gentoyped on (ichip1-6, ichip7, bbc). For each cohort, ancestry estimation was performed using the ADMIXTURE tool. Allele frequencies learned in a supervised manner from labeled data provided by 1000 Genomes samples belonging to the EUR, AFR and EAS populations. Then ancestry proportions were estimated and samples with &gt; 75% estimated european ancestry were split into a european cohort. A standardized QC was performed for each of the 5 ancestry defined cohorts and then vcf files per chromosome were uploaded to the Michigan Imputation Server. The hrc.r1.1.2016 reference panel was used for european ancestry samples, phasing was performed by Eagle and imputation was performed by Minimac. Post imputation markers were filtered for R2 &gt; 0.3 and stored per chromosome and all chromosomes combined.
