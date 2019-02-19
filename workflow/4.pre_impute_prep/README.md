Prepare Genotype Data for Imputation - README
================
Prepare Genotype Data for Imputation
2/22/2019

Tasks
-----

European cohorts will be imputed against the Haplotype\_Reference\_Consortium reference

1.  iChip data was aligned to HRC data using will rayner strand file and check vcf tools
    -   Flip SNPs to Convert TOP to Forward/Ref
    -   Set Ref/Alt Allele
    -   Remove SNPs with different positions
    -   Convert to VCF using Plink
    -   Create Sorted and Compressed VCF using VCFtools and tabix
2.  Files were submitted to Micigan Imputation Server for Imputation
    -   European Reference Panel: hrc.r1.1.2016
    -   Phasing was performed by Eagle
    -   Imputation was perfoemed by Minimac
