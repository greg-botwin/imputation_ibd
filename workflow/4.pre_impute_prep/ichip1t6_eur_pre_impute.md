iChip 1 - 6 European Samples Prepare for Imputation Submission
================
Translational Genomics Group
29 October, 2018

Prepare VCF for Imputation Submission
-------------------------------------

### Download HRC Reference

This file has already been downloaded to a shared directory.

``` bash
wget -p /mnt/share6/SHARED_DATASETS/Haplotype_Reference_Consortium ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz &>/dev/null
```

``` bash
gunzip -k /mnt/share6/SHARED_DATASETS/Haplotype_Reference_Consortium/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz &>/dev/null
```

### Download HRC Check Tool

``` bash
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.9.zip -O temp.zip; unzip temp.zip; rm temp.zip
```

### Download Will Rayner Strand FIle to Align to Reference

Download Immuno\_BeadChip\_11419691\_B Strand Files for iChip 1 to 6 and Ref/Alt File

``` bash
wget http://www.well.ox.ac.uk/~wrayner/strand/Immuno_BeadChip_11419691_B-b37-strand.zip -O temp.zip; unzip temp.zip; rm temp.zip &>/dev/null

wget http://www.well.ox.ac.uk/~wrayner/strand/RefAlt/Immuno_BeadChip_11419691_B-b37.strand.RefAlt.zip -O temp.zip; unzip temp.zip; rm temp.zip &>/dev/null
```

### B. Create List of SNPs to Flip

``` r
strand <- read_tsv("immuno_beadchip_11419691_b-b37.strand", 
                   col_names = c("SNP", "Chr", "BP", "Match", "Strand", "Allele"))

strand %>%
  filter(Strand == "-") %>%
  select(SNP) %>%
  write_tsv("temp_flip_ichip1t6.tsv", col_names = FALSE)
```

### C. Flip SNPs to Convert TOP to Forward/Ref in Plink

``` bash
plink \
--bfile ../3.qc/qc_eur_cohort_split_ichip1t6 \
--flip temp_flip_ichip1t6.tsv \
--make-bed \
--out temp_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_ichip1t6.log.
    ## Options in effect:
    ##   --bfile ../3.qc/qc_eur_cohort_split_ichip1t6
    ##   --flip temp_flip_ichip1t6.tsv
    ##   --make-bed
    ##   --out temp_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 114148 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## --flip: 57179 SNPs flipped, 40514 SNP IDs not present.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 42 het. haploid genotypes present (see temp_ichip1t6.hh ); many
    ## commands treat these as missing.
    ## Total genotyping rate is 0.999835.
    ## 114148 variants and 8150 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp_ichip1t6.bed + temp_ichip1t6.bim + temp_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### D. Set Ref/Alt Allele

``` bash
plink \
--bfile temp_ichip1t6 \
--reference-allele Immuno_BeadChip_11419691_B-b37.strand.RefAlt 2 \
--chr 1-22 \
--make-bed \
--out temp1_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp1_ichip1t6.log.
    ## Options in effect:
    ##   --a1-allele Immuno_BeadChip_11419691_B-b37.strand.RefAlt 2
    ##   --bfile temp_ichip1t6
    ##   --chr 1-22
    ##   --make-bed
    ##   --out temp1_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 113489 out of 114148 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.999836.
    ## --a1-allele: 112524 assignments made.
    ## 113489 variants and 8150 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp1_ichip1t6.bed + temp1_ichip1t6.bim + temp1_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### E. Calcualte Allele Frequency

``` bash
plink \
--bfile temp1_ichip1t6 \
--keep-allele-order \
--freq \
--out temp1_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp1_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp1_ichip1t6
    ##   --freq
    ##   --keep-allele-order
    ##   --out temp1_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 113489 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.999836.
    ## --freq: Allele frequencies (founders only) written to temp1_ichip1t6.frq .

### F. Final Checks

``` bash
perl HRC-1000G-check-bim.pl \
-b temp1_ichip1t6.bim \
-f temp1_ichip1t6.frq \
-r /mnt/share6/SHARED_DATASETS/Haplotype_Reference_Consortium/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-h
```

-   4 markers, match by name to HRC but have a different position (same CHR), will exclude. rs1377587 rs1319548 rs9258651 rs2085508
-   66 markers are on the wrong strand when compared to HRC reference after strand flips from Will's files. Will re-flip possible difference references used.
-   854 with wrong ref assignemnt

### Remove 4 SNPs with different positions

``` bash
plink \
--bfile temp1_ichip1t6 \
--exclude Position-temp1_ichip1t6-HRC.txt \
--make-bed \
--out temp2_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp2_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp1_ichip1t6
    ##   --exclude Position-temp1_ichip1t6-HRC.txt
    ##   --make-bed
    ##   --out temp2_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 113489 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## --exclude: 113485 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.999836.
    ## 113485 variants and 8150 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp2_ichip1t6.bed + temp2_ichip1t6.bim + temp2_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### Exclude SNPS

``` bash
plink \
--bfile temp2_ichip1t6 \
--exclude Exclude-temp1_ichip1t6-HRC.txt \
--make-bed \
--out temp3_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp3_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp2_ichip1t6
    ##   --exclude Exclude-temp1_ichip1t6-HRC.txt
    ##   --make-bed
    ##   --out temp3_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 113485 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## --exclude: 109633 variants remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.999838.
    ## 109633 variants and 8150 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp3_ichip1t6.bed + temp3_ichip1t6.bim + temp3_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### Flip

``` bash
plink \
--bfile temp3_ichip1t6 \
--flip Strand-Flip-temp1_ichip1t6-HRC.txt \
--make-bed \
--out temp4_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp4_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp3_ichip1t6
    ##   --flip Strand-Flip-temp1_ichip1t6-HRC.txt
    ##   --make-bed
    ##   --out temp4_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 109633 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## --flip: 66 SNPs flipped.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.999838.
    ## 109633 variants and 8150 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp4_ichip1t6.bed + temp4_ichip1t6.bim + temp4_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

### Set Ref/Alt

``` bash
plink \
--bfile temp4_ichip1t6 \
--reference-allele Force-Allele1-temp1_ichip1t6-HRC.txt \
--make-bed \
--out ichip1t6_eur_qc_pre_impute
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to ichip1t6_eur_qc_pre_impute.log.
    ## Options in effect:
    ##   --a1-allele Force-Allele1-temp1_ichip1t6-HRC.txt
    ##   --bfile temp4_ichip1t6
    ##   --make-bed
    ##   --out ichip1t6_eur_qc_pre_impute
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 109633 variants loaded from .bim file.
    ## 8150 people (4125 males, 4025 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 8150 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.999838.
    ## --a1-allele: 109633 assignments made.
    ## 109633 variants and 8150 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to ichip1t6_eur_qc_pre_impute.bed + ichip1t6_eur_qc_pre_impute.bim +
    ## ichip1t6_eur_qc_pre_impute.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

``` r
freq_plot <- read_tsv("FreqPlot-temp1_ichip1t6-HRC.txt", col_names = FALSE)
```

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_character(),
    ##   X2 = col_double(),
    ##   X3 = col_double(),
    ##   X4 = col_double(),
    ##   X5 = col_integer()
    ## )

``` r
freq_plot %>%
  ggplot(aes(x = X2, y = X3)) +
  geom_point(alpha = 1/500)
```

![](ichip1t6_eur_pre_impute_files/figure-markdown_github/unnamed-chunk-14-1.png)

Convert to VCF using Plink
--------------------------

Using Plink1.9 becuase plink2 outputs VCFv4.3 which is not yet supported by the server. Plink1.9 on vcf recode defualts to A2 allele as reference, so forcing a2 allele to curent a1.

``` bash
for i in {1..22}; do
    plink --bfile ichip1t6_eur_qc_pre_impute \
    --chr ${i} \
    --real-ref-alleles \
    --recode vcf \
    --out ichip1t6_eur_vcfs/ichip1t6_eur_chr${i}
done
```

Create Sorted and Compressed VCF using VCFtools and tabix (including bgzip)
---------------------------------------------------------------------------

``` bash
for i in {1..22}; do
    vcf-sort ichip1t6_eur_vcfs/ichip1t6_eur_chr${i}.vcf | bgzip -c > ichip1t6_eur_vcfs/ichip1t6_eur_chr${i}.vcf.gz
done
```

``` r
file.remove(list.files(pattern = "^temp", full.names = TRUE))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

Submit to Michigan Imputation Server.
