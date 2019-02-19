Imputing IBD Immuno Chip Genetic Data
================
Translational Genomics Group
8/31/2018

Organization
------------

Each workflow step has its own folder and README. Often each workflow step is repeated for each of our identified cohorts. Genetic data and derivatives are not updated to this repository, but are avaialble.

Original Data
-------------

1.  Genotype data from IBDichip1to7eBBCTOP\_unfilteredHG19b.\*
2.  Updated annotation file 25\_jul\_2018
3.  List of batches each sample was run under for cohort splits (e.g. ichip 1 vs 7, from PCA\_Investigation Project)

Workflow
--------

1.  Split genetic data into relevant cohorts [split\_cohorts](workflow/1.split_cohorts/)
2.  Estimate ancestry and split cohorts based on ancestry [estimate\_ancestry](workflow/2.estimate_ancestry/)
3.  Perform ancestry specific QC [qc](3.qc/)
4.  Prepare genotype data for Imputation QC [pre\_impute\_prep](workflow/4.pre_impute_prep/)
5.  Post imputation QC [post\_impute\_qc](workflow/5.post_impute_qc/)

Summary
-------

All samples were gentoyped using the Illumina ImmunoChip platform and markers were assigned hg19 coordinates based on the manufacturers alignment files. Markers mapping to insertions and deletions, non-automosomal chromosomes, or were unable to be mapped to an hg19 coordiante were removed.
