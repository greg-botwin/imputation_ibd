Imputing IBD Immuno Chip Genetic Data
================
Translational Genomics Group
8/31/2018

Organization
------------

Each workflow step has its own folder and README. Often each workflow step is repeated for each of our identified cohorts. Genetic data and derivatives are not updated to this repository, but are avaialble.

Workflow
--------

1.  Split genetic data into relevant cohorts [split\_cohorts](split_cohorts/)
2.  Estimate ancestry and split cohorts based on ancestry [estimate\_ancestry](estimate_ancestry/)
3.  Perform ancestry specific QC [qc](qc/)
4.  Prepare genotype data for Imputation QC [pre\_impute\_prep](pre_impute_prep/)
5.  Post imputation QC [post\_impute\_qc](post_impute_qc/)
