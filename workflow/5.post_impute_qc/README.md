Prepare Genotype Data for Imputation - README
================
Prepare Genotype Data for Imputation
2/22/2019

Tasks
-----

1.  Post Imputation QC Plots were generated for all SNPs per cohort and per chromosome
2.  Variants were filtered with `bcftools filter -i 'R2>.3'`
3.  Each chromosome VCF was indexed `tabix -p vcf`
4.  An avaialble concatenated vcf across all chromosomes is avaialble
5.  The data is stored in: `/mnt/YanX/imputation_results/ichip1t6_eur/extracted_vcfs/filtered` `/mnt/YanX/imputation_results/ichip7_eur/extracted_vcfs/filtered`
