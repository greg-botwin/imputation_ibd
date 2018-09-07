Split Cohorts and Perform Global Operations - Workflow
================
Translational Genomics Group
06 September, 2018

1. Update Bim Positions
-----------------------

### A. Read in Annotation File and Bim File

``` r
# read anno file
annotation <- read_tsv("../../original_data/ichip_v2_v1_grch37_no_indels_annovar_annotation_25_jul_2018.tsv", 
                       col_types = cols(.default = "c"))

# select only necessary columns and change class
annotation <- annotation %>%
  select(Name, IlmnID, Chr, SNP, MapInfo_GRCh37, deCODEcM_GRCh37, RsID_dbsnp138, 
         For_Rev_Allele1_GRCh37, For_Rev_Allele2_GRCh37, Top_AlleleA, Top_AlleleB,
         MappingComment, insertion_deleton, platform, Ref_GRCh37, Alt_GRCh37,
         bim_A1, bim_A2) %>%
  mutate(Chr = as.integer(Chr), MapInfo_GRCh37 = as.integer(MapInfo_GRCh37))

# read bim
bim <- read_table2("../../original_data/IBDichip1to7eBBCTOP_unfilteredHG19b.bim",
                  col_names = c("CHR", "Name", "POS", "BP", "A1", "A2"))
```

    ## Parsed with column specification:
    ## cols(
    ##   CHR = col_integer(),
    ##   Name = col_character(),
    ##   POS = col_integer(),
    ##   BP = col_integer(),
    ##   A1 = col_character(),
    ##   A2 = col_character()
    ## )

### B. Combine and Write File with Name, New BP

``` r
# combine bim and annotation
bim_annotated <- left_join(bim, annotation, by = "Name")

# write update map file
bim_annotated %>%
  select(Name, MapInfo_GRCh37) %>%
  mutate(MapInfo_GRCh37 = if_else(is.na(MapInfo_GRCh37), as.integer(0), MapInfo_GRCh37)) %>%
  write_tsv("temp_update_map.tsv", col_names = FALSE)
```

### C. Update Map in Plink

``` bash
plink \
--bfile ../../original_data/IBDichip1to7eBBCTOP_unfilteredHG19b \
--update-map temp_update_map.tsv 2 1 \
--make-bed \
--out temp_1
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_1.log.
    ## Options in effect:
    ##   --bfile ../../original_data/IBDichip1to7eBBCTOP_unfilteredHG19b
    ##   --make-bed
    ##   --out temp_1
    ##   --update-map temp_update_map.tsv 2 1
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 274093 variants loaded from .bim file.
    ## 14856 people (7208 males, 7609 females, 39 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_1.nosex .
    ## 4224 phenotype values loaded from .fam.
    ## --update-map: 274093 values updated.
    ## Warning: Base-pair positions are now unsorted!
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 14856 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 19438 het. haploid genotypes present (see temp_1.hh ); many commands
    ## treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.547589.
    ## 274093 variants and 14856 people pass filters and QC.
    ## Among remaining phenotypes, 0 are cases and 4224 are controls.  (10632
    ## phenotypes are missing.)
    ## --make-bed to temp_1.bed + temp_1.bim + temp_1.fam ... 0%0%1%1%2%2%3%3%4%4%5%5%6%6%7%7%8%8%9%9%10%10%11%11%12%12%13%13%14%14%15%15%16%16%17%17%18%18%19%19%20%20%21%21%22%22%23%23%24%24%25%25%26%26%27%27%28%28%29%29%30%30%31%31%32%32%33%33%34%34%35%35%36%36%37%37%38%38%39%39%40%40%41%41%42%42%43%43%44%44%45%45%46%46%47%47%48%48%49%49%50%50%51%51%52%52%53%53%54%54%55%55%56%56%57%57%58%58%59%59%60%60%61%61%62%62%63%63%64%64%65%65%66%66%67%67%68%68%69%69%70%70%71%71%72%72%73%73%74%74%75%75%76%76%77%77%78%78%79%79%80%80%81%81%82%82%83%83%84%84%85%85%86%86%87%87%88%88%89%89%90%90%91%91%92%92%93%93%94%94%95%95%96%96%97%97%98%98%99%done.

2. Filter Markers for Issues
----------------------------

Remove SNPs with:
1. Has a mapping issue
2. No BP Position (doesn't add anything)
3. Insertion/Deletion

### A. Create List of Passing SNPs

``` r
bim_annotated %>%
  filter(MappingComment == "None" | 
           MappingComment == "updated map position differs from position in ichip1to7top_unfiltered_e_bim" |
           MappingComment == "Previously mapped to Chr 5, but probe correctly matched to Chr 15.") %>%
  filter(!is.na(MapInfo_GRCh37)) %>%
  filter(insertion_deleton == FALSE) %>%
  select(Name) %>%
  write_tsv("temp_passing_markers.tsv", col_names = FALSE)
```

### B. Filter Genotype Data

``` bash
plink \
--bfile temp_1 \
--extract temp_passing_markers.tsv \
--make-bed \
--out temp_2
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_2.log.
    ## Options in effect:
    ##   --bfile temp_1
    ##   --extract temp_passing_markers.tsv
    ##   --make-bed
    ##   --out temp_2
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 274093 variants loaded from .bim file.
    ## 14856 people (7208 males, 7609 females, 39 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_2.nosex .
    ## 4224 phenotype values loaded from .fam.
    ## --extract: 270813 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 14856 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 16041 het. haploid genotypes present (see temp_2.hh ); many commands
    ## treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.550839.
    ## 270813 variants and 14856 people pass filters and QC.
    ## Among remaining phenotypes, 0 are cases and 4224 are controls.  (10632
    ## phenotypes are missing.)
    ## --make-bed to temp_2.bed + temp_2.bim + temp_2.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

3. Seperate into the 3 cohorts
------------------------------

I am going to seperate the samples into 3 cohorts. - iChip 1 - 6 - iChip 7 - BBC

This is necessary because if I include all the samples in one batch, the missingness is high for specific SNPs and those SNPs fail pre-imputation QC on the server.

### A. Create Lists of of Samples Identifying Which Cohort Each Sample Belongs to

There are 83 samples that I could not assign batch on. These samples generally have updated genetic ids or are replicated. I assigned all of these samples to the ichip 1 to 6 cohort.

``` r
fam <- read_table2("../../original_data/IBDichip1to7eBBCTOP_unfilteredHG19b.fam", 
                   col_names = c("FID", "IID", "WFID_F", "WFID_M", "SEX", "PHENO"))
```

    ## Parsed with column specification:
    ## cols(
    ##   FID = col_character(),
    ##   IID = col_integer(),
    ##   WFID_F = col_integer(),
    ##   WFID_M = col_integer(),
    ##   SEX = col_integer(),
    ##   PHENO = col_integer()
    ## )

``` r
batches <- read_tsv("../../original_data/batch_tidy_all.tsv", col_types = cols(LABID2 = col_character(),
                                                                               Run = col_character()))


bbc <- batches %>%
  filter(Run == "BBC") %>%
  select(LABID2) %>%
  write_tsv("temp_bbc_samples.tsv", col_names = FALSE)

ichip7 <- batches %>%
  filter(Run == "Immunochip7") %>%
  select(LABID2) %>%
  write_tsv("temp_ichip7_samples.tsv", col_names = FALSE)

fam %>%
  filter(!FID %in% c(bbc$LABID2, ichip7$LABID2)) %>%
  select(FID) %>%
  write_tsv("temp_ichip1t6_samples.tsv", col_names = FALSE)
```

### B. Create ichip1t6 Cohort

``` bash
plink \
--bfile temp_2 \
--geno 0.99 \
--keep-fam temp_ichip1t6_samples.tsv \
--make-bed \
--out cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp_2
    ##   --geno 0.99
    ##   --keep-fam temp_ichip1t6_samples.tsv
    ##   --make-bed
    ##   --out cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 270813 variants loaded from .bim file.
    ## 14856 people (7208 males, 7609 females, 39 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to cohort_split_ichip1t6.nosex .
    ## 4224 phenotype values loaded from .fam.
    ## --keep-fam: 9971 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9971 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 9917 het. haploid genotypes present (see cohort_split_ichip1t6.hh );
    ## many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.559539.
    ## 91819 variants removed due to missing genotype data (--geno).
    ## 178994 variants and 9971 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to cohort_split_ichip1t6.bed + cohort_split_ichip1t6.bim +
    ## cohort_split_ichip1t6.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### C. Create ichip7 Cohort

``` bash
plink \
--bfile temp_2 \
--geno 0.99 \
--keep-fam temp_ichip7_samples.tsv \
--make-bed \
--out cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile temp_2
    ##   --geno 0.99
    ##   --keep-fam temp_ichip7_samples.tsv
    ##   --make-bed
    ##   --out cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 270813 variants loaded from .bim file.
    ## 14856 people (7208 males, 7609 females, 39 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to cohort_split_ichip7.nosex .
    ## 4224 phenotype values loaded from .fam.
    ## --keep-fam: 661 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 661 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see cohort_split_ichip7.hh ); many
    ## commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.915066.
    ## 22910 variants removed due to missing genotype data (--geno).
    ## 247903 variants and 661 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to cohort_split_ichip7.bed + cohort_split_ichip7.bim +
    ## cohort_split_ichip7.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### D. Create BBC Cohort

``` bash
plink \
--bfile temp_2 \
--geno 0.99 \
--keep-fam temp_bbc_samples.tsv \
--make-bed \
--out cohort_split_bbc
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to cohort_split_bbc.log.
    ## Options in effect:
    ##   --bfile temp_2
    ##   --geno 0.99
    ##   --keep-fam temp_bbc_samples.tsv
    ##   --make-bed
    ##   --out cohort_split_bbc
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 270813 variants loaded from .bim file.
    ## 14856 people (7208 males, 7609 females, 39 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to cohort_split_bbc.nosex .
    ## 4224 phenotype values loaded from .fam.
    ## --keep-fam: 4224 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 4224 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 6114 het. haploid genotypes present (see cohort_split_bbc.hh ); many
    ## commands treat these as missing.
    ## Total genotyping rate in remaining samples is 0.473304.
    ## 142621 variants removed due to missing genotype data (--geno).
    ## 128192 variants and 4224 people pass filters and QC.
    ## Among remaining phenotypes, 0 are cases and 4224 are controls.
    ## --make-bed to cohort_split_bbc.bed + cohort_split_bbc.bim +
    ## cohort_split_bbc.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

Clean Up
--------

``` r
file.remove(list.files(pattern = "^temp", full.names = TRUE))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE
