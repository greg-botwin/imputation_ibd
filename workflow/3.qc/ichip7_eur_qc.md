iChip 7 European QC
================
Translational Genomics Group
07 September, 2018

Set-up
------

This document is broken down into three parts. - Part 1 Identifies QC issues - Part 2 Fixes and applies QC parameters - Part 3 Re-Checks genotype data post QC Changes

Source Data
-----------

European IBD IChip, Run 7, Illumina's TOP Allele Called, aligned to HG19.

Part 1 - Identify QC Issues
===========================

Identification of individuals with discordant sex information WITH --geno and --maf filters
-------------------------------------------------------------------------------------------

-   Plink uses chrX data to determine sex (based on heterozygosity rates).
-   Default PLINK thresholds of .2 for Females and .8 for Males when assessing homozygosity rates, but use .4 and .8.
-   When the homozygosity rate is more than 0.4 but less than 0.8, the genotype data are inconclusive regarding the sex of an individual and these are marked in column 4 with a 0.
-   maf and geno filter added with check-sex run
-   PEDSEX=sex as recorded in pedfile (1=male, 2=female)
-   SNPSEX=sex as predicted based on genetic data (1=male, 2=female, 0=unknown)
-   If needed, compare with GenomeStudio gender estimates to ID true problems
-   We suggest running --check-sex once without parameters, eyeballing the distribution of F estimates (there should be a clear gap between a very tight male clump at the right side of the distribution and the females everywhere else), and then rerunning with parameters corresponding to the empirical gap.

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--geno 0.03 \
--maf 0.05 \
--check-sex .4 .8 \
--out temp_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile ../2.estimate_ancestry/eur_cohort_split_ichip7
    ##   --check-sex .4 .8
    ##   --geno 0.03
    ##   --maf 0.05
    ##   --out temp_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 580 people (273 males, 273 females, 34 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_eur_cohort_split_ichip7.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 580 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see temp_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999625.
    ## 18 variants removed due to missing genotype data (--geno).
    ## 66680 variants removed due to minor allele threshold(s)
    ## (--maf/--max-maf/--mac/--max-mac).
    ## 181205 variants and 580 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --check-sex: 788 Xchr and 0 Ychr variant(s) scanned, 38 problems detected.
    ## Report written to temp_eur_cohort_split_ichip7.sexcheck .

Using 788 Xchr variants to check sex; 38 problems detected

``` r
sexcheck <- read_table2(file = "temp_eur_cohort_split_ichip7.sexcheck")
sexcheck %>%
  filter(STATUS == "PROBLEM") %>%
  kable(caption = "List of Individuals with Sex Check Problems")
```

| FID        |  IID|  PEDSEX|  SNPSEX| STATUS  |         F|
|:-----------|----:|-------:|-------:|:--------|---------:|
| 0500928    |    1|       2|       0| PROBLEM |  0.418900|
| 1201279    |    1|       2|       0| PROBLEM |  0.400700|
| 1400212    |    1|       2|       0| PROBLEM |  0.433400|
| 1500141    |    1|       2|       0| PROBLEM |  0.426100|
| Control1   |    1|       0|       2| PROBLEM |  0.008438|
| Control2   |    1|       0|       2| PROBLEM |  0.008438|
| Control3   |    1|       0|       2| PROBLEM |  0.008438|
| Control4   |    1|       0|       2| PROBLEM |  0.008438|
| Control5   |    1|       0|       2| PROBLEM |  0.008438|
| Control6   |    1|       0|       2| PROBLEM |  0.008438|
| Control7   |    1|       0|       2| PROBLEM |  0.007384|
| Control8   |    1|       0|       2| PROBLEM |  0.010750|
| Control9   |    1|       0|       2| PROBLEM |  0.008438|
| Control10  |    1|       0|       2| PROBLEM |  0.010300|
| Control11  |    1|       0|       2| PROBLEM |  0.008438|
| Control12  |    1|       0|       2| PROBLEM |  0.008438|
| Control13  |    1|       0|       2| PROBLEM |  0.008438|
| Control14  |    1|       0|       2| PROBLEM |  0.008438|
| Control15  |    1|       0|       2| PROBLEM |  0.008438|
| Control16  |    1|       0|       2| PROBLEM |  0.008438|
| Control17  |    1|       0|       2| PROBLEM |  0.008438|
| Control18  |    1|       0|       2| PROBLEM |  0.008438|
| Control19  |    1|       0|       2| PROBLEM |  0.008438|
| Control20  |    1|       0|       2| PROBLEM |  0.008438|
| Control21  |    1|       0|       2| PROBLEM |  0.008438|
| Control22  |    1|       0|       2| PROBLEM |  0.008438|
| Control23  |    1|       0|       2| PROBLEM |  0.008438|
| Control24  |    1|       0|       2| PROBLEM |  0.008438|
| Control25  |    1|       0|       2| PROBLEM |  0.008438|
| Control26  |    1|       0|       2| PROBLEM |  0.008438|
| Control27  |    1|       0|       2| PROBLEM |  0.008438|
| Control29  |    1|       0|       2| PROBLEM |  0.008438|
| Control29R |    1|       0|       2| PROBLEM |  0.007384|
| Control30  |    1|       0|       2| PROBLEM |  0.008438|
| Control30R |    1|       0|       2| PROBLEM |  0.008438|
| Control31  |    1|       0|       2| PROBLEM |  0.008438|
| Control32  |    1|       0|       2| PROBLEM |  0.008438|
| Control33  |    1|       0|       2| PROBLEM |  0.008438|

38 problems detected, but this is beause many control samples are blank for Sex. Not real problems. 4 actual sex problems.

``` r
sexcheck %>%
  filter(STATUS == "PROBLEM") %>%
  filter(PEDSEX != 0) %>%
  kable(caption = "Sex Problems with a Reported Sex in Ped")
```

| FID     |  IID|  PEDSEX|  SNPSEX| STATUS  |       F|
|:--------|----:|-------:|-------:|:--------|-------:|
| 0500928 |    1|       2|       0| PROBLEM |  0.4189|
| 1201279 |    1|       2|       0| PROBLEM |  0.4007|
| 1400212 |    1|       2|       0| PROBLEM |  0.4334|
| 1500141 |    1|       2|       0| PROBLEM |  0.4261|

4 of 4 sex problems come from a reported PEDSEX of Female but an F &gt; 0.4.

``` r
sexcheck %>%
  ggplot(aes(x = F)) +
  geom_density() +
  labs(title = "Heterozygoisty Rates for Samples Filtered with geno/maf included") 
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-4-1.png)

Plotting F from check-sex based on chrX heterozygosity rates (F&lt;0.4, M&gt;0.8 expected). We can see that they nicely seperate as expected.

Identification of individuals with elevated missing data rates or outlying heterozygosity rate
----------------------------------------------------------------------------------------------

### Number of missing SNPS and Proportion of Missings SNPs per Individual

-   no filters included here
-   {output.imiss} for individuals (F\_MISS will give proportion of missing SNPs per individual)
-   {output.lmiss} for snps (F\_MISS will give proportion of samples missing per SNP)
-   stricter missingness should apply for low MAF snps
-   evaluate SNP missingness rates per cohort, genotyping batch, case-control status

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--missing \
--out temp_eur_cohort_split_ichip7 
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile ../2.estimate_ancestry/eur_cohort_split_ichip7
    ##   --missing
    ##   --out temp_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 580 people (273 males, 273 females, 34 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_eur_cohort_split_ichip7.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 580 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see temp_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999625.
    ## --missing: Sample missing data report written to
    ## temp_eur_cohort_split_ichip7.imiss, and variant-based missing data report
    ## written to temp_eur_cohort_split_ichip7.lmiss.

Total genotyping rate is 0.999625.

``` r
lmiss <- read_table("temp_eur_cohort_split_ichip7.lmiss")
lmiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_density() +
  labs(title = "SNP Level Missingess")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-6-1.png)

Low missingness across iChip7

``` r
imiss <- read_table("temp_eur_cohort_split_ichip7.imiss")
imiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_histogram(bins = 100) +
  labs(title = "Sample Level Missingess")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-7-1.png)

Small amount of missingness

``` r
imiss %>%
  filter(F_MISS >= 0.03) 
```

    ## # A tibble: 0 x 6
    ## # ... with 6 variables: FID <chr>, IID <int>, MISS_PHENO <chr>,
    ## #   N_MISS <int>, N_GENO <int>, F_MISS <dbl>

No samples with &gt; 3% missingness.

### Samples with outlying heterozygosity rates

-   To calculate individual inbreeding F / heterozygosity With whole genome data, can be applied to pruned subset plink output includes:
-   O(HOM) Observed number of homozygotes
-   E(HOM) Expected number of homozygotes
-   N(NM) Number of non-missing genotypes
-   F F inbreeding coefficient estimate
-   approx norm range: -0.2 to 0.2, mode around zero
-   Low heterozygosity (high F; positive value) may indicate inbreeding;
-   high heterozygosity (low F; negative value) may indicate contamination.

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--het \
--out temp_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile ../2.estimate_ancestry/eur_cohort_split_ichip7
    ##   --het
    ##   --out temp_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 580 people (273 males, 273 females, 34 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_eur_cohort_split_ichip7.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 580 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see temp_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999625.
    ## 247903 variants and 580 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --het: 228885 variants scanned, report written to
    ## temp_eur_cohort_split_ichip7.het .

``` r
het <- read_table("temp_eur_cohort_split_ichip7.het") 
het <- het %>%
  mutate(obs_het_rate = (`N(NM)` - `O(HOM)`)/`E(HOM)`)
imiss_het <- left_join(imiss, het, by = "FID")

ggplot(imiss_het, aes(x = F_MISS, y = obs_het_rate)) +
  geom_point(color = densCols(log10(imiss_het$F_MISS), imiss_het$obs_het_rate)) +
  labs(x = "Proportion of missing genotypes", y = "Heterozygosity rate") +
  scale_x_log10(limits = c(0.001, 1), minor_breaks = c(0.01, 0.1)) +
  scale_y_continuous(limits = c(0, .5)) +
  geom_vline(xintercept = 0.03, color = "red") +
  geom_hline(yintercept = (mean(imiss_het$obs_het_rate)+(3*sd(imiss_het$obs_het_rate))), color = "red") +
  geom_hline(yintercept = (mean(imiss_het$obs_het_rate)-(3*sd(imiss_het$obs_het_rate))), color = "red")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-10-1.png)

### Identify het outliers

``` r
imiss_het %>%
  filter(obs_het_rate >= (mean(imiss_het$obs_het_rate)+(3*sd(imiss_het$obs_het_rate))) |
           obs_het_rate <= (mean(imiss_het$obs_het_rate)-(3*sd(imiss_het$obs_het_rate)))) %>%
  kable(caption = "Samples That Fail Het Check")
```

| FID      |  IID.x| MISS\_PHENO |  N\_MISS|  N\_GENO|    F\_MISS|  IID.y|  O(HOM)|  E(HOM)|   N(NM)|         F|  obs\_het\_rate|
|:---------|------:|:------------|--------:|--------:|----------:|------:|-------:|-------:|-------:|---------:|---------------:|
| 1400090  |      1| Y           |      682|   247903|  0.0027510|      1|  170585|  165400|  228338|   0.08205|       0.3491717|
| 1400152  |      1| Y           |       32|   247903|  0.0001291|      1|  170130|  165800|  228854|   0.06881|       0.3541858|
| 1400186  |      1| Y           |       33|   247903|  0.0001331|      1|  161621|  165800|  228854|  -0.06607|       0.4055066|
| 1500019  |      1| Y           |       42|   247903|  0.0001694|      1|  171164|  165800|  228847|   0.08530|       0.3479071|
| 1500052  |      1| Y           |      164|   247903|  0.0006615|      1|  161324|  165700|  228734|  -0.06944|       0.4068196|
| 1600092S |      1| Y           |       57|   246816|  0.0002309|      1|  172252|  165800|  228836|   0.10270|       0.3412786|

``` r
imiss_het %>%
  filter(obs_het_rate >= (mean(imiss_het$obs_het_rate)+(3*sd(imiss_het$obs_het_rate))) |
           obs_het_rate <= (mean(imiss_het$obs_het_rate)-(3*sd(imiss_het$obs_het_rate)))) %>%
  dplyr::select(FID, IID.x) %>%
  write_tsv("fail-het-outlier.txt", col_names = FALSE)
```

1 samples failed het check at +/- 4SD, 6 +/1 3SD. Only use 4SD with &gt; 1000 Samples.

Calculate Overall project MAF before Filtering
----------------------------------------------

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--freq \
--out temp_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile ../2.estimate_ancestry/eur_cohort_split_ichip7
    ##   --freq
    ##   --out temp_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 580 people (273 males, 273 females, 34 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_eur_cohort_split_ichip7.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 580 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see temp_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999625.
    ## --freq: Allele frequencies (founders only) written to
    ## temp_eur_cohort_split_ichip7.frq .

``` r
maffreq <- read_table2("temp_eur_cohort_split_ichip7.frq")
maffreq %>%
  ggplot(aes(x = MAF)) +
  geom_histogram(aes(y =..density..)) +
  geom_density(col=2) +
  labs(title = "Overall MAF Prior to Filtering")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-13-1.png)

Identification of Duplicated or Related Individuals
---------------------------------------------------

-   Prune dataset for temporary use of calculating cryptic relatedness and PCA as both work best under assumption of no LD among SNPs
-   Prior to calcuating identity by state, IBS, prune SNPs to only independent SNPs and remove regions with extended linkage disequiblibirum such as HLA region. The below removes snps within a 50kb window, with an r2 &gt; .2 and variant count to shift the window at the end of each step of 5kb.

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--exclude range ../../original_data/highLDregions.txt \
--indep 50 5 1.8 \
--out temp_eur_cohort_split_ichip7 &>/dev/null
```

Pruning complete. 12368 variants excluded removed in high ld. 77285 of 121314 variants removed. I will keep only the prune.in snps for the subsequent analysis.

-   Can add --min 0.12 to identify minimum pihat for genome output to manage size of output dataset (will only output pihat &gt;0.12)
-   PIHAT 1.0 = monozygotic twins or known replicates
-   PIHAT 0.5 = 1st degree relatives: P-C, sibs
-   PIHAT 0.25= 2nd degree relatives: half-sib, grandparents
-   PIHAT 0.125= 3rd degree relatives: full cousins

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--extract temp_eur_cohort_split_ichip7.prune.in \
--genome \
--min 0.12 \
--out temp_eur_cohort_split_ichip7 &>/dev/null
```

``` r
genome <- read_table2("temp_eur_cohort_split_ichip7.genome")

genome <- genome %>%
  mutate(PI_HAT = as.double(PI_HAT)) %>%
  mutate(color = if_else(PI_HAT <.15, "~3rd degree",
                         if_else(PI_HAT >=.15 & PI_HAT <.35, "~2nd degree",
                                 if_else(PI_HAT >= .35 & PI_HAT < .65, "~1st degree",
                                         if_else(PI_HAT > .65, "~Replicates or twins", "???")))))

genome %>%
  ggplot(aes(x = as.double(Z0), y = as.double(Z1), color = color)) +
  geom_point(alpha = 1/20) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x = "Z0 the proportion of loci where the pair shares zero alleles", 
       y = "Z1 the proportion of loci where the pair shares one allele") 
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-16-1.png)

The above plot is restricted to samples with a PI\_HAT atleast 0.12. Samples can be represented more than once if multiple relations are found.

``` r
genome %>%
  ggplot(aes(x = as.double(PI_HAT), fill = color)) +
  geom_histogram(bins = 100) +
  labs(title = "Distribution of PI_HAT for Related Individuals >0.12", x = "PI_HAT", y = "Count (non-unique)")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
genome %>%
  filter(PI_HAT > 0.8)
```

    ## # A tibble: 567 x 15
    ##    FID1   IID1 FID2  IID2 RT    EZ       Z0    Z1    Z2 PI_HAT   PHE   DST
    ##    <chr> <int> <ch> <int> <chr> <chr> <dbl> <dbl> <dbl>  <dbl> <int> <dbl>
    ##  1 1400…     1 140…     1 UN    <NA>      0     0     1      1    -1     1
    ##  2 1400…     1 140…     1 UN    <NA>      0     0     1      1    -1     1
    ##  3 1500…     1 150…     1 UN    <NA>      0     0     1      1    -1     1
    ##  4 1500…     1 150…     1 UN    <NA>      0     0     1      1    -1     1
    ##  5 1500…     1 150…     1 UN    <NA>      0     0     1      1    -1     1
    ##  6 1500…     1 150…     1 UN    <NA>      0     0     1      1    -1     1
    ##  7 Cont…     1 Con…     1 UN    <NA>      0     0     1      1    -1     1
    ##  8 Cont…     1 Con…     1 UN    <NA>      0     0     1      1    -1     1
    ##  9 Cont…     1 Con…     1 UN    <NA>      0     0     1      1    -1     1
    ## 10 Cont…     1 Con…     1 UN    <NA>      0     0     1      1    -1     1
    ## # ... with 557 more rows, and 3 more variables: PPC <dbl>, RATIO <dbl>,
    ## #   color <chr>

``` r
genome %>%
  filter(PI_HAT > 0.8) %>%
  write_tsv("possible_duplicates.tsv")
```

Since I am only interested in removing duplicates, will set a PI-HAT threshold of 0.8. IN other analysis, we might want to restrict to non-relateds.

Part 2 Apply QC Filters and Fix Errors
======================================

Clean Sexes and remove samples and SNPs that fail heterozygoisty, genotyping, sex-checks and IBD &gt;= 0.8/
-----------------------------------------------------------------------------------------------------------

### Clean Sexs

#### First update sex for subjects with PEDSEX == 0 and approriate F stat

If PED sex 0 and F &lt; 0.4 updated PED to Female 2 (not needed in this set) If PED sex 0 and F &gt; 0.8 update PED to male 1 All samples with previously un-assigned PED SEX assigned

``` r
sexcheck %>%
  filter(PEDSEX == 0) %>%
  mutate(PEDSEX = ifelse(F < 0.4, 2,
                          ifelse(F > 0.8, 1, PEDSEX))) %>% 
  filter(PEDSEX != 0) %>%
  dplyr::select(FID, IID, PEDSEX)  %>%
  write_tsv("update-missing-sex.txt", col_names = FALSE)
  
sexcheck %>%
  filter(PEDSEX == 0) %>%
  mutate(PEDSEX = ifelse(F < 0.4, 2,
                          ifelse(F > 0.8, 1, PEDSEX))) %>% 
  filter(PEDSEX != 0) %>%
  dplyr::select(FID, IID, PEDSEX) %>%
  kable(caption = "Subjects to Update Sex")
```

| FID        |  IID|  PEDSEX|
|:-----------|----:|-------:|
| Control1   |    1|       2|
| Control2   |    1|       2|
| Control3   |    1|       2|
| Control4   |    1|       2|
| Control5   |    1|       2|
| Control6   |    1|       2|
| Control7   |    1|       2|
| Control8   |    1|       2|
| Control9   |    1|       2|
| Control10  |    1|       2|
| Control11  |    1|       2|
| Control12  |    1|       2|
| Control13  |    1|       2|
| Control14  |    1|       2|
| Control15  |    1|       2|
| Control16  |    1|       2|
| Control17  |    1|       2|
| Control18  |    1|       2|
| Control19  |    1|       2|
| Control20  |    1|       2|
| Control21  |    1|       2|
| Control22  |    1|       2|
| Control23  |    1|       2|
| Control24  |    1|       2|
| Control25  |    1|       2|
| Control26  |    1|       2|
| Control27  |    1|       2|
| Control29  |    1|       2|
| Control29R |    1|       2|
| Control30  |    1|       2|
| Control30R |    1|       2|
| Control31  |    1|       2|
| Control32  |    1|       2|
| Control33  |    1|       2|

#### Create List of Subjects that Fail Sex Check After Update

``` r
sexcheck %>%
  filter(PEDSEX == 1 & SNPSEX == 2 |
           PEDSEX == 2 & SNPSEX == 1 |
           PEDSEX == 2 & SNPSEX == 0 & F > 0.4 |
           PEDSEX == 1 & SNPSEX == 0 & F < 0.8) %>%
  kable(caption = "Subjects that Fail Sex Check")
```

| FID     |  IID|  PEDSEX|  SNPSEX| STATUS  |       F|
|:--------|----:|-------:|-------:|:--------|-------:|
| 0500928 |    1|       2|       0| PROBLEM |  0.4189|
| 1201279 |    1|       2|       0| PROBLEM |  0.4007|
| 1400212 |    1|       2|       0| PROBLEM |  0.4334|
| 1500141 |    1|       2|       0| PROBLEM |  0.4261|

``` r
sexcheck %>%
  filter(PEDSEX == 1 & SNPSEX == 2 |
           PEDSEX == 2 & SNPSEX == 1 |
           PEDSEX == 2 & SNPSEX == 0 & F > 0.4 |
           PEDSEX == 1 & SNPSEX == 0 & F < 0.8) %>%
  dplyr::select(FID, IID) %>%
  write_tsv(path = "fail-updated-sex-check.txt", col_names = FALSE)
```

These are 3 female samples, that have out of range but quasi-reasonable F. I will remove. Disucssed with Talin.

#### Update Sex in PLINK

``` bash
plink \
--bfile ../2.estimate_ancestry/eur_cohort_split_ichip7 \
--update-sex update-missing-sex.txt \
--remove fail-updated-sex-check.txt \
--make-bed \
--out temp1_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp1_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile ../2.estimate_ancestry/eur_cohort_split_ichip7
    ##   --make-bed
    ##   --out temp1_eur_cohort_split_ichip7
    ##   --remove fail-updated-sex-check.txt
    ##   --update-sex update-missing-sex.txt
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 580 people (273 males, 273 females, 34 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp1_eur_cohort_split_ichip7.nosex .
    ## --update-sex: 34 people updated.
    ## --remove: 576 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 576 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see
    ## temp1_eur_cohort_split_ichip7.hh ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999623.
    ## 247903 variants and 576 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp1_eur_cohort_split_ichip7.bed +
    ## temp1_eur_cohort_split_ichip7.bim + temp1_eur_cohort_split_ichip7.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

--update-sex: 34 people updated. --Removed 4 samples. Went from 580 to 876 subjects. All F

### Exclude Samples with IBD &gt; 0.8

For the purposes of this QC, I want to identify all pairs of individuals with an IBD &gt;=0.8 and then remove the individual in the pair with the lower genotyping rate. IN some cases, you will want to remove the individual with the lower genotyping rate dependent on the phenotype. The ID of the individuals with the lower genotype rate will be stored in ‘fail-IBD-QC.txt’ for subsequent removal.

``` r
relateds <- genome %>%
  filter(PI_HAT > 0.8)

compare_missingness <- function(FID1, IID1, FID2, IID2, imiss){
  output <- data.frame(fid1_miss = numeric(length(FID1)),
                       fid2_miss = numeric(length(FID1)))
  for(i in seq_along(FID1)) {
    fid1_miss <- imiss %>%
      filter(FID == FID1[[i]] & IID == IID1[[i]]) %>%
      select(F_MISS) %>%
      unlist() 
    
    fid2_miss <- imiss %>%
      filter(FID == FID2[[i]] & IID == IID2[[i]]) %>%
      select(F_MISS) %>%
      unlist()
    
    output[i,] <- c(fid1_miss, fid2_miss)
  }
  return(output)
}

output <- compare_missingness(relateds$FID1, relateds$IID1, relateds$FID2, relateds$IID2, imiss)

relateds <- relateds %>%
  bind_cols(., output) %>%
  mutate(better_sample = if_else(fid1_miss <= fid2_miss, paste("FID1"), paste("FID2")))

duplicates_to_remove <- relateds %>%
  mutate(remove_fid = if_else(better_sample == "FID1", FID2,
                              if_else(better_sample == "FID2", FID1, "Wilson"))) %>%
  mutate(remove_iid = ifelse(better_sample == "FID1", IID2,
                              ifelse(better_sample == "FID2", IID1, "Wilson"))) %>%
  select(remove_fid, remove_iid) %>%
  distinct() %>%
  write_tsv("fail-IBD-QC.tsv", col_names = FALSE)
```

``` bash
plink \
--bfile temp1_eur_cohort_split_ichip7 \
--remove fail-IBD-QC.tsv \
--make-bed \
--out temp2_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp2_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile temp1_eur_cohort_split_ichip7
    ##   --make-bed
    ##   --out temp2_eur_cohort_split_ichip7
    ##   --remove fail-IBD-QC.tsv
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 576 people (273 males, 303 females) loaded from .fam.
    ## --remove: 537 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 537 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see
    ## temp2_eur_cohort_split_ichip7.hh ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999606.
    ## 247903 variants and 537 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp2_eur_cohort_split_ichip7.bed +
    ## temp2_eur_cohort_split_ichip7.bim + temp2_eur_cohort_split_ichip7.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

8292 people to 8159 people. 133 duplicates removed.

### Failed SNPs missingness &gt;97% and MAF 1%, and HWE 10^-6, Sample missingnes 3%

``` bash
plink \
--bfile temp2_eur_cohort_split_ichip7 \
--geno 0.03 \
--maf 0.01 \
--hwe 0.000001 \
--mind 0.03 \
--make-bed \
--out temp3_eur_cohort_split_ichip7 
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp3_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile temp2_eur_cohort_split_ichip7
    ##   --geno 0.03
    ##   --hwe 0.000001
    ##   --maf 0.01
    ##   --make-bed
    ##   --mind 0.03
    ##   --out temp3_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 247903 variants loaded from .bim file.
    ## 537 people (271 males, 266 females) loaded from .fam.
    ## 0 people removed due to missing genotype data (--mind).
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 537 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see
    ## temp3_eur_cohort_split_ichip7.hh ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999606.
    ## 19 variants removed due to missing genotype data (--geno).
    ## Warning: --hwe observation counts vary by more than 10%, due to the X
    ## chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    ## chromosome variants.
    ## --hwe: 12 variants removed due to Hardy-Weinberg exact test.
    ## 35747 variants removed due to minor allele threshold(s)
    ## (--maf/--max-maf/--mac/--max-mac).
    ## 212125 variants and 537 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp3_eur_cohort_split_ichip7.bed +
    ## temp3_eur_cohort_split_ichip7.bim + temp3_eur_cohort_split_ichip7.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

19 variants removed due to missing genotype data (--geno). 12 variants removed due to Hardy-Weinberg exact test (--hwe) 35747 variants removed due to minor allele threshold(s)(--maf) 0 people removed due to missing genotype data (--mind)

### Remove samples that fail het check

``` bash
plink \
--bfile temp3_eur_cohort_split_ichip7 \
--remove fail-het-outlier.txt \
--make-bed \
--out qc_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile temp3_eur_cohort_split_ichip7
    ##   --make-bed
    ##   --out qc_eur_cohort_split_ichip7
    ##   --remove fail-het-outlier.txt
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 212125 variants loaded from .bim file.
    ## 537 people (271 males, 266 females) loaded from .fam.
    ## --remove: 531 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 531 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see qc_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999627.
    ## 212125 variants and 531 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to qc_eur_cohort_split_ichip7.bed + qc_eur_cohort_split_ichip7.bim +
    ## qc_eur_cohort_split_ichip7.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

537 people to 531 people. 6 samples as ientified.

Part 3 Re-Check File Post QC
============================

Recheck Missingess and Sex with Cleaned File
--------------------------------------------

``` bash
plink \
--bfile qc_eur_cohort_split_ichip7 \
--missing \
--out qc_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile qc_eur_cohort_split_ichip7
    ##   --missing
    ##   --out qc_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 212125 variants loaded from .bim file.
    ## 531 people (266 males, 265 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 531 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see qc_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999627.
    ## --missing: Sample missing data report written to
    ## qc_eur_cohort_split_ichip7.imiss, and variant-based missing data report written
    ## to qc_eur_cohort_split_ichip7.lmiss.

``` r
lmiss <- read_table("qc_eur_cohort_split_ichip7.lmiss")
lmiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_histogram(bins = 10) +
  labs(title = "SNP Level Missingess after SNP level Filtering")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-27-1.png)

``` r
imiss <- read_table("qc_eur_cohort_split_ichip7.imiss")
imiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_histogram(bins = 100) +
  labs(title = "Sample Level Missingess after Sample Level Filtering")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-28-1.png) Batch effect resolved.

Re-Check Hets
-------------

``` bash
plink \
--bfile qc_eur_cohort_split_ichip7 \
--het \
--out qc_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile qc_eur_cohort_split_ichip7
    ##   --het
    ##   --out qc_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 212125 variants loaded from .bim file.
    ## 531 people (266 males, 265 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 531 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see qc_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999627.
    ## 212125 variants and 531 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --het: 210901 variants scanned, report written to
    ## qc_eur_cohort_split_ichip7.het .

``` r
imiss <- read_table("qc_eur_cohort_split_ichip7.imiss")
het <- read_table("qc_eur_cohort_split_ichip7.het") %>%
  mutate(obs_het_rate = (`N(NM)` - `O(HOM)`)/`E(HOM)`)
imiss_het <- left_join(imiss, het, by = "FID")

ggplot(imiss_het, aes(x = F_MISS, y = obs_het_rate)) +
  geom_point(color = densCols(log10(imiss_het$F_MISS), imiss_het$obs_het_rate)) +
  labs(x = "Proportion of missing genotypes", y = "Heterozygosity rate") +
  scale_x_log10(limits = c(0.001, 1), minor_breaks = c(0.01, 0.1)) +
  scale_y_continuous(limits = c(0, .5)) +
  geom_vline(xintercept = 0.03, color = "red") +
  geom_hline(yintercept = (mean(imiss_het$obs_het_rate)+(3*sd(imiss_het$obs_het_rate))), color = "red") +
  geom_hline(yintercept = (mean(imiss_het$obs_het_rate)-(3*sd(imiss_het$obs_het_rate))), color = "red")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-30-1.png) Looks great.

Calculate Overall project MAF Adter Filtering
---------------------------------------------

``` bash
plink \
--bfile qc_eur_cohort_split_ichip7 \
--freq \
--out qc_eur_cohort_split_ichip7
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_eur_cohort_split_ichip7.log.
    ## Options in effect:
    ##   --bfile qc_eur_cohort_split_ichip7
    ##   --freq
    ##   --out qc_eur_cohort_split_ichip7
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 212125 variants loaded from .bim file.
    ## 531 people (266 males, 265 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 531 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 10 het. haploid genotypes present (see qc_eur_cohort_split_ichip7.hh
    ## ); many commands treat these as missing.
    ## Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
    ## treat these as missing.
    ## Total genotyping rate is 0.999627.
    ## --freq: Allele frequencies (founders only) written to
    ## qc_eur_cohort_split_ichip7.frq .

``` r
maffreq <- read_table2("qc_eur_cohort_split_ichip7.frq")
maffreq %>%
  ggplot(aes(x = MAF)) +
  geom_histogram(aes(y =..density..)) +
  geom_density(col=2) +
  labs(title = "Overall MAF After Filtering")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-32-1.png)

Identification of Duplicated or Related Individuals
---------------------------------------------------

-   Prune dataset for temporary use of calculating cryptic relatedness and PCA as both work best under assumption of no LD among SNPs
-   Prior to calcuating identity by state, IBS, prune SNPs to only independent SNPs and remove regions with extended linkage disequiblibirum such as HLA region. The below removes snps within a 50kb window, with an r2 &gt; .2 and variant count to shift the window at the end of each step of 5kb.

``` bash
plink \
--bfile qc_eur_cohort_split_ichip7 \
--exclude ../../original_data/highLDregions.txt \
--range \
--indep 50 5 1.8 \
--out qc_eur_cohort_split_ichip7 &>/dev/null
```

``` bash
plink \
--bfile qc_eur_cohort_split_ichip7 \
--extract qc_eur_cohort_split_ichip7.prune.in \
--genome \
--min 0.12 \
--out qc_eur_cohort_split_ichip7 &>/dev/null
```

``` r
genome <- read_table2("qc_eur_cohort_split_ichip7.genome")

genome <- genome %>%
  mutate(PI_HAT = as.double(PI_HAT)) %>%
  mutate(color = if_else(PI_HAT <.15, "~3rd degree",
                         if_else(PI_HAT >=.15 & PI_HAT <.35, "~2nd degree",
                                 if_else(PI_HAT >= .35 & PI_HAT < .65, "~1st degree",
                                         if_else(PI_HAT > .65, "~Replicates or twins", "???")))))

genome %>%
  ggplot(aes(x = as.double(Z0), y = as.double(Z1), color = color)) +
  geom_point(alpha =1 ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  labs(x = "Z0 the proportion of loci where the pair shares zero alleles", 
       y = "Z1 the proportion of loci where the pair shares one allele") 
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-35-1.png)

Replicates (purple) removed

``` r
genome %>%
  ggplot(aes(x = as.double(PI_HAT), fill = color)) +
  geom_histogram(bins = 100) +
  labs(title = "Distribution of PI_HAT for Related Individuals >0.12", x = "PI_HAT", y = "Count (non-unique)")
```

![](ichip7_eur_qc_files/figure-markdown_github/unnamed-chunk-36-1.png)

``` r
file.remove(list.files(pattern = "^temp", full.names = TRUE))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
