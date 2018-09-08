iChip 1 through 6 All Ancestry QC
================
Translational Genomics Group
08 September, 2018

Set-up
------

This document is broken down into three parts. - Part 1 Identifies QC issues - Part 2 Fixes and applies QC parameters - Part 3 Re-Checks genotype data post QC Changes

Source Data
-----------

All Ancestry IChip, Runs 1-6, Illumina's TOP Allele Called, including GOLDR removed BBC, filtered to SNPS in common, and aligned to HG19.

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
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--geno 0.03 \
--maf 0.05 \
--check-sex .4 .8 \
--out temp_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile ../1.split_cohorts/cohort_split_ichip1t6
    ##   --check-sex .4 .8
    ##   --geno 0.03
    ##   --maf 0.05
    ##   --out temp_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9971 people (4912 males, 5054 females, 5 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_all_cohort_split_ichip1t6.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9971 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate is 0.998891.
    ## 514 variants removed due to missing genotype data (--geno).
    ## 36439 variants removed due to minor allele threshold(s)
    ## (--maf/--max-maf/--mac/--max-mac).
    ## 96297 variants and 9971 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --check-sex: 599 Xchr and 0 Ychr variant(s) scanned, 18 problems detected.
    ## Report written to temp_all_cohort_split_ichip1t6.sexcheck .

Using 599 Xchr variants to check sex; 18 problems detected

``` r
sexcheck <- read_table2(file = "temp_all_cohort_split_ichip1t6.sexcheck")
sexcheck %>%
  filter(STATUS == "PROBLEM") %>%
  kable(caption = "List of Individuals with Sex Check Problems")
```

| FID      |  IID|  PEDSEX|  SNPSEX| STATUS  |        F|
|:---------|----:|-------:|-------:|:--------|--------:|
| 0201354  |   11|       2|       0| PROBLEM |   0.4379|
| 0201478  |   11|       2|       0| PROBLEM |   0.4276|
| 0602622  |    1|       2|       0| PROBLEM |   0.4101|
| 0700060  |    1|       2|       0| PROBLEM |   0.4101|
| 0700795  |    1|       2|       0| PROBLEM |   0.4008|
| 0801111  |    1|       2|       0| PROBLEM |   0.4147|
| 0901475  |    1|       2|       0| PROBLEM |   0.4008|
| 1101433  |    1|       0|       2| PROBLEM |  -0.3006|
| 1101582  |    1|       0|       1| PROBLEM |   1.0000|
| 1101587  |   11|       0|       2| PROBLEM |   0.1091|
| 1101600  |    1|       0|       2| PROBLEM |  -0.1613|
| 1101662  |   11|       0|       1| PROBLEM |   1.0000|
| 1300076  |   11|       2|       0| PROBLEM |   0.4658|
| 1300792  |   11|       2|       0| PROBLEM |   0.4472|
| 1300792R |   11|       2|       0| PROBLEM |   0.4472|
| 1400078  |   11|       2|       0| PROBLEM |   0.4797|
| 1400079  |   11|       2|       0| PROBLEM |   0.4472|
| 1400397  |    1|       2|       0| PROBLEM |   0.4426|

``` r
sexcheck %>%
  filter(STATUS == "PROBLEM") %>%
  filter(PEDSEX != 0) %>%
  kable(caption = "Sex Problems with a Reported Sex in Ped")
```

| FID      |  IID|  PEDSEX|  SNPSEX| STATUS  |       F|
|:---------|----:|-------:|-------:|:--------|-------:|
| 0201354  |   11|       2|       0| PROBLEM |  0.4379|
| 0201478  |   11|       2|       0| PROBLEM |  0.4276|
| 0602622  |    1|       2|       0| PROBLEM |  0.4101|
| 0700060  |    1|       2|       0| PROBLEM |  0.4101|
| 0700795  |    1|       2|       0| PROBLEM |  0.4008|
| 0801111  |    1|       2|       0| PROBLEM |  0.4147|
| 0901475  |    1|       2|       0| PROBLEM |  0.4008|
| 1300076  |   11|       2|       0| PROBLEM |  0.4658|
| 1300792  |   11|       2|       0| PROBLEM |  0.4472|
| 1300792R |   11|       2|       0| PROBLEM |  0.4472|
| 1400078  |   11|       2|       0| PROBLEM |  0.4797|
| 1400079  |   11|       2|       0| PROBLEM |  0.4472|
| 1400397  |    1|       2|       0| PROBLEM |  0.4426|

13 of 18 sex problems come from a reported PEDSEX of Female but an F &gt; 0.4.

``` r
sexcheck %>%
  filter(STATUS == "PROBLEM") %>%
  filter(PEDSEX == 0) %>%
  kable(caption = "Sex Probelms without Sex Recorded in Ped")
```

| FID     |  IID|  PEDSEX|  SNPSEX| STATUS  |        F|
|:--------|----:|-------:|-------:|:--------|--------:|
| 1101433 |    1|       0|       2| PROBLEM |  -0.3006|
| 1101582 |    1|       0|       1| PROBLEM |   1.0000|
| 1101587 |   11|       0|       2| PROBLEM |   0.1091|
| 1101600 |    1|       0|       2| PROBLEM |  -0.1613|
| 1101662 |   11|       0|       1| PROBLEM |   1.0000|

The remaining 5 problem comes from a sample without a reported PEDSEX, but are easily classified as Male or Female.

``` r
sexcheck %>%
  ggplot(aes(x = F)) +
  geom_density() +
  labs(title = "Heterozygoisty Rates for Samples Filtered with geno/maf included") 
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-5-1.png)

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
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--missing \
--out temp_all_cohort_split_ichip1t6 
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile ../1.split_cohorts/cohort_split_ichip1t6
    ##   --missing
    ##   --out temp_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9971 people (4912 males, 5054 females, 5 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_all_cohort_split_ichip1t6.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9971 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate is 0.998891.
    ## --missing: Sample missing data report written to
    ## temp_all_cohort_split_ichip1t6.imiss, and variant-based missing data report
    ## written to temp_all_cohort_split_ichip1t6.lmiss.

Total genotyping rate is 0.998891. The batches have aready been filtered for SNPs in common that passed QC at each run.

``` r
lmiss <- read_table("temp_all_cohort_split_ichip1t6.lmiss")
lmiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_density() +
  labs(title = "SNP Level Missingess")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-7-1.png)

Low missingness across iChip1-6

``` r
imiss <- read_table("temp_all_cohort_split_ichip1t6.imiss")
imiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_histogram(bins = 100) +
  labs(title = "Sample Level Missingess")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-8-1.png)

Small amount of missingness, likely batch related due to unique peaks, but still quite low.

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
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--het \
--out temp_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile ../1.split_cohorts/cohort_split_ichip1t6
    ##   --het
    ##   --out temp_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9971 people (4912 males, 5054 females, 5 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_all_cohort_split_ichip1t6.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9971 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate is 0.998891.
    ## 133250 variants and 9971 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --het: 132405 variants scanned, report written to
    ## temp_all_cohort_split_ichip1t6.het .

``` r
het <- read_table("temp_all_cohort_split_ichip1t6.het") 
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

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-11-1.png)

### Identify het outliers

``` r
imiss_het %>%
  filter(obs_het_rate >= (mean(imiss_het$obs_het_rate)+(4*sd(imiss_het$obs_het_rate))) |
           obs_het_rate <= (mean(imiss_het$obs_het_rate)-(4*sd(imiss_het$obs_het_rate)))) %>%
  kable(caption = "Samples That Fail Het Check")
```

| FID       |  IID.x| MISS\_PHENO |  N\_MISS|  N\_GENO|    F\_MISS|  IID.y|  O(HOM)|  E(HOM)|   N(NM)|       F|  obs\_het\_rate|
|:----------|------:|:------------|--------:|--------:|----------:|------:|-------:|-------:|-------:|-------:|---------------:|
| 06029     |      1| Y           |       16|   133250|  0.0001201|      1|  105767|   98900|  132389|  0.2051|       0.2691810|
| 07022     |      1| Y           |       35|   133250|  0.0002627|      1|  105394|   98880|  132370|  0.1944|       0.2728155|
| 08001     |      1| Y           |       27|   133250|  0.0002026|      1|  105216|   98890|  132378|  0.1889|       0.2746688|
| 0400202   |      1| Y           |      520|   133250|  0.0039020|      1|  104994|   98540|  131892|  0.1936|       0.2729653|
| 0400203   |      1| Y           |       23|   133250|  0.0001726|      1|  106783|   98890|  132382|  0.2356|       0.2588634|
| 0401725   |      1| Y           |       18|   133250|  0.0001351|      1|  105626|   98900|  132387|  0.2009|       0.2705865|
| 0500316   |      1| Y           |      525|   133250|  0.0039400|      1|  105133|   98530|  131887|  0.1979|       0.2715315|
| 0600489   |      1| Y           |       15|   133250|  0.0001126|      1|  105747|   98900|  132391|  0.2044|       0.2694034|
| 0602142   |      1| Y           |       31|   133250|  0.0002326|      1|  105286|   98890|  132374|  0.1911|       0.2739205|
| 0602246   |      1| Y           |      520|   133250|  0.0039020|      1|  104569|   98540|  131892|  0.1809|       0.2772783|
| 0700085   |      1| Y           |       12|   133250|  0.0000901|      1|  105513|   98900|  132393|  0.1974|       0.2717897|
| 0700152   |      1| Y           |       18|   133250|  0.0001351|      1|  105752|   98900|  132387|  0.2047|       0.2693124|
| 0700154   |      1| Y           |        9|   133250|  0.0000675|      1|  105135|   98900|  132396|  0.1861|       0.2756421|
| 0700772   |      1| Y           |       30|   133250|  0.0002251|      1|  108375|   98890|  132375|  0.2833|       0.2426939|
| 0701351   |      1| Y           |      538|   133250|  0.0040380|      1|  105616|   98530|  131875|  0.2126|       0.2665077|
| 0701352   |      1| Y           |      529|   133250|  0.0039700|      1|  105553|   98530|  131883|  0.2106|       0.2672283|
| 0701430   |      1| Y           |      615|   133250|  0.0046150|      1|  105777|   98470|  131797|  0.2192|       0.2642429|
| 0701440   |      1| Y           |       20|   133250|  0.0001501|      1|  104944|   98900|  132387|  0.1805|       0.2774823|
| 0800169   |      1| Y           |      526|   133250|  0.0039470|      1|  104803|   98530|  131887|  0.1880|       0.2748807|
| 0800329   |      1| Y           |      532|   133250|  0.0039920|      1|  105657|   98530|  131880|  0.2138|       0.2661423|
| 0800487   |      1| Y           |      558|   133250|  0.0041880|      1|  104464|   98510|  131854|  0.1786|       0.2780428|
| 0800918   |      1| Y           |      562|   133250|  0.0042180|      1|  104981|   98510|  131850|  0.1942|       0.2727540|
| 0801111   |      1| Y           |       35|   133250|  0.0002627|      1|  105260|   98890|  132370|  0.1904|       0.2741430|
| 0801197   |      1| Y           |      532|   133250|  0.0039920|      1|  104947|   98530|  131880|  0.1924|       0.2733482|
| 0801332   |      1| Y           |      522|   133250|  0.0039170|      1|  107152|   98540|  131890|  0.2583|       0.2510453|
| 0801349   |      1| Y           |      583|   133250|  0.0043750|      1|  104695|   98490|  131830|  0.1861|       0.2755102|
| 0801400   |      1| Y           |      523|   133250|  0.0039250|      1|  105424|   98530|  131889|  0.2065|       0.2685984|
| 0900158   |      1| Y           |       33|   133250|  0.0002477|      1|  104907|   98890|  132372|  0.1798|       0.2777328|
| 0900227   |      1| Y           |       36|   133250|  0.0002702|      1|  106314|   98880|  132369|  0.2219|       0.2635012|
| 0900303   |      1| Y           |       11|   133250|  0.0000826|      1|  104960|   98900|  132394|  0.1809|       0.2773913|
| 0900519   |      1| Y           |      522|   133250|  0.0039170|      1|  105903|   98540|  131890|  0.2209|       0.2637203|
| 0900580   |      1| Y           |      525|   133250|  0.0039400|      1|  104508|   98530|  131887|  0.1791|       0.2778748|
| 0900585   |      1| Y           |       14|   133250|  0.0001051|      1|  106414|   98900|  132391|  0.2244|       0.2626593|
| 0900645   |      1| Y           |      549|   133250|  0.0041200|      1|  104507|   98510|  131863|  0.1797|       0.2776977|
| 0900925   |      1| Y           |      520|   133250|  0.0039020|      1|  104893|   98540|  131892|  0.1905|       0.2739903|
| 0900927   |      1| Y           |      528|   133250|  0.0039620|      1|  105716|   98530|  131884|  0.2154|       0.2655841|
| 0901335   |      1| Y           |       21|   133250|  0.0001576|      1|  104916|   98890|  132384|  0.1798|       0.2777632|
| 0901558   |      1| Y           |      520|   133250|  0.0039020|      1|  104751|   98540|  131892|  0.1863|       0.2754313|
| 0901560   |      1| Y           |      523|   133250|  0.0039250|      1|  105005|   98530|  131889|  0.1940|       0.2728509|
| 0901714   |      1| Y           |      518|   133250|  0.0038870|      1|  105048|   98540|  131894|  0.1951|       0.2724376|
| 0901715   |      1| Y           |      520|   133250|  0.0039020|      1|  104555|   98540|  131892|  0.1804|       0.2774203|
| 0901716   |      1| Y           |      532|   133250|  0.0039920|      1|  104513|   98530|  131880|  0.1794|       0.2777530|
| 1000074   |      1| Y           |      521|   133250|  0.0039100|      1|  105540|   98540|  131891|  0.2100|       0.2674142|
| 1000075   |      1| Y           |      518|   133250|  0.0038870|      1|  106144|   98540|  131894|  0.2280|       0.2613152|
| 1000206   |      1| Y           |      518|   133250|  0.0038870|      1|  108952|   98540|  131894|  0.3122|       0.2328192|
| 1000207   |      1| Y           |      520|   133250|  0.0039020|      1|  109027|   98540|  131892|  0.3145|       0.2320378|
| 1000208   |      1| Y           |      515|   133250|  0.0038650|      1|  108107|   98540|  131897|  0.2868|       0.2414248|
| 1000818   |      1| Y           |       55|   133250|  0.0004128|      1|  104870|   98870|  132351|  0.1792|       0.2779508|
| 1000820   |      1| Y           |       33|   133250|  0.0002477|      1|  105475|   98880|  132372|  0.1968|       0.2720166|
| 1000935   |      1| Y           |       17|   133250|  0.0001276|      1|  105233|   98900|  132388|  0.1892|       0.2745703|
| 1001514   |      1| Y           |        9|   133250|  0.0000675|      1|  105786|   98910|  132397|  0.2054|       0.2690426|
| 1001891   |      1| Y           |      528|   133250|  0.0039620|      1|  107246|   98530|  131884|  0.2613|       0.2500558|
| 1001892   |      1| Y           |      525|   133250|  0.0039400|      1|  105874|   98530|  131887|  0.2201|       0.2640110|
| 1100169   |      1| Y           |      551|   133250|  0.0041350|      1|  106036|   98510|  131861|  0.2255|       0.2621561|
| 1100970   |      1| Y           |      521|   133250|  0.0039100|      1|  104626|   98540|  131891|  0.1825|       0.2766897|
| 1101021   |      1| Y           |      530|   133250|  0.0039770|      1|  104563|   98530|  131882|  0.1809|       0.2772658|
| 1101021R  |      1| Y           |      521|   133250|  0.0039100|      1|  104572|   98540|  131891|  0.1810|       0.2772377|
| 1101275   |      1| Y           |      530|   133250|  0.0039770|      1|  105437|   98530|  131882|  0.2071|       0.2683954|
| 1101516   |      1| Y           |      524|   133250|  0.0039320|      1|  105474|   98540|  131888|  0.2080|       0.2680536|
| 1200124   |      1| Y           |      522|   133250|  0.0039170|      1|  105300|   98540|  131890|  0.2028|       0.2698397|
| 1200247   |      1| Y           |      545|   133250|  0.0040900|      1|  105356|   98520|  131867|  0.2050|       0.2690926|
| 1200319   |      1| Y           |      523|   133250|  0.0039250|      1|  105045|   98540|  131889|  0.1952|       0.2724173|
| 1200636   |      1| Y           |      524|   133250|  0.0039320|      1|  105124|   98530|  131888|  0.1976|       0.2716330|
| 1300050   |      1| Y           |       27|   133250|  0.0002026|      1|  106258|   98890|  132378|  0.2200|       0.2641319|
| 1300076   |     11| Y           |       27|   133250|  0.0002026|     11|  104976|   98890|  132378|  0.1817|       0.2770958|
| 1300716   |      1| Y           |       20|   133250|  0.0001501|      1|  104977|   98900|  132385|  0.1816|       0.2771284|
| 1400078   |     11| Y           |       28|   133250|  0.0002101|     11|  104887|   98890|  132377|  0.1791|       0.2779856|
| 1400147   |      1| Y           |       14|   133250|  0.0001051|      1|  106316|   98900|  132391|  0.2214|       0.2636502|
| 8900681   |      1| Y           |       17|   133250|  0.0001276|      1|  106328|   98900|  132388|  0.2218|       0.2634985|
| 9600127   |      1| Y           |       21|   133250|  0.0001576|      1|  105290|   98890|  132385|  0.1910|       0.2739913|
| 9901000   |     11| Y           |       20|   133250|  0.0001501|     11|  104912|   98900|  132385|  0.1796|       0.2777856|
| 9901005   |      1| Y           |      533|   133250|  0.0040000|      1|  105744|   98530|  131879|  0.2164|       0.2652492|
| 9901324   |      1| Y           |      577|   133250|  0.0043300|      1|  106079|   98500|  131836|  0.2274|       0.2614924|
| 9901458   |      1| Y           |       40|   133250|  0.0003002|      1|  104916|   98880|  132365|  0.1803|       0.2775991|
| 100020902 |      1| Y           |      527|   133250|  0.0039550|      1|  110025|   98530|  131885|  0.3446|       0.2218614|

``` r
imiss_het %>%
  filter(obs_het_rate >= (mean(imiss_het$obs_het_rate)+(4*sd(imiss_het$obs_het_rate))) |
           obs_het_rate <= (mean(imiss_het$obs_het_rate)-(4*sd(imiss_het$obs_het_rate)))) %>%
  dplyr::select(FID, IID.x) %>%
  write_tsv("fail-het-outlier.txt", col_names = FALSE)
```

75 samples failed het check. They are overalll pretty good. Using 4 SD due to large sample.

Calculate Overall project MAF before Filtering
----------------------------------------------

``` bash
plink \
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--freq \
--out temp_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile ../1.split_cohorts/cohort_split_ichip1t6
    ##   --freq
    ##   --out temp_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9971 people (4912 males, 5054 females, 5 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp_all_cohort_split_ichip1t6.nosex .
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9971 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate is 0.998891.
    ## --freq: Allele frequencies (founders only) written to
    ## temp_all_cohort_split_ichip1t6.frq .

``` r
maffreq <- read_table2("temp_all_cohort_split_ichip1t6.frq")
maffreq %>%
  ggplot(aes(x = MAF)) +
  geom_histogram(aes(y =..density..)) +
  geom_density(col=2) +
  labs(title = "Overall MAF Prior to Filtering")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-14-1.png)

Identification of Duplicated or Related Individuals
---------------------------------------------------

-   Prune dataset for temporary use of calculating cryptic relatedness and PCA as both work best under assumption of no LD among SNPs
-   Prior to calcuating identity by state, IBS, prune SNPs to only independent SNPs and remove regions with extended linkage disequiblibirum such as HLA region. The below removes snps within a 50kb window, with an r2 &gt; .2 and variant count to shift the window at the end of each step of 5kb.

``` bash
plink \
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--exclude range ../../original_data/highLDregions.txt \
--indep 50 5 1.8 \
--out temp_all_cohort_split_ichip1t6 &>/dev/null
```

Pruning complete. 12358 variants excluded removed in high ld. 77007 of 120892 variants removed variants removed. I will keep only the prune.in snps for the subsequent analysis.

-   Can add --min 0.12 to identify minimum pihat for genome output to manage size of output dataset (will only output pihat &gt;0.12)
-   PIHAT 1.0 = monozygotic twins or known replicates
-   PIHAT 0.5 = 1st degree relatives: P-C, sibs
-   PIHAT 0.25= 2nd degree relatives: half-sib, grandparents
-   PIHAT 0.125= 3rd degree relatives: full cousins

``` bash
plink \
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--extract temp_all_cohort_split_ichip1t6.prune.in \
--genome \
--min 0.12 \
--out temp_all_cohort_split_ichip1t6 &>/dev/null
```

``` r
genome <- read_table2("temp_all_cohort_split_ichip1t6.genome")

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

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-17-1.png)

We have a related cohort. The above plot is restricted to samples with a PI\_HAT atleast 0.12. Samples can be represented more than once if multiple relations are found.

``` r
genome %>%
  ggplot(aes(x = as.double(PI_HAT), fill = color)) +
  geom_histogram(bins = 100) +
  labs(title = "Distribution of PI_HAT for Related Individuals >0.12", x = "PI_HAT", y = "Count (non-unique)")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
genome %>%
  filter(PI_HAT > 0.8)
```

    ## # A tibble: 155 x 15
    ##    FID1   IID1 FID2  IID2 RT    EZ       Z0    Z1    Z2 PI_HAT   PHE   DST
    ##    <chr> <int> <ch> <int> <chr> <chr> <dbl> <dbl> <dbl>  <dbl> <int> <dbl>
    ##  1 0000…     1 000…     1 UN    <NA>      0     0     1      1    -1     1
    ##  2 0000…     1 000…     1 UN    <NA>      0     0     1      1    -1     1
    ##  3 0000…     1 000…     1 UN    <NA>      0     0     1      1    -1     1
    ##  4 06003     1 070…     1 UN    <NA>      0     0     1      1    -1     1
    ##  5 06006     1 110…     1 UN    <NA>      0     0     1      1    -1     1
    ##  6 06010     1 110…     1 UN    <NA>      0     0     1      1    -1     1
    ##  7 06013     1 110…     1 UN    <NA>      0     0     1      1    -1     1
    ##  8 06015     1 110…     1 UN    <NA>      0     0     1      1    -1     1
    ##  9 06017     1 110…     1 UN    <NA>      0     0     1      1    -1     1
    ## 10 06032     1 110…     1 UN    <NA>      0     0     1      1    -1     1
    ## # ... with 145 more rows, and 3 more variables: PPC <dbl>, RATIO <dbl>,
    ## #   color <chr>

``` r
genome %>%
  filter(PI_HAT > 0.8) %>%
  write_tsv("possible_duplicates.tsv")
```

Since I am only interested in removing duplicates, will set a PI-HAT threshold of 0.8. IN alter analysis, we might want to restrict to non-relateds.

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

| FID        |    IID|        PEDSEX|
|:-----------|------:|-------------:|
| 1101433    |      1|             2|
| 1101582    |      1|             1|
| 1101587    |     11|             2|
| 1101600    |      1|             2|
| 1101662    |     11|             1|
| 5 subjects |  to up|  date PED sex|

#### Create List of Subjects that Fail Sex Check After Update

``` r
sexcheck %>%
  filter(PEDSEX == 1 & SNPSEX == 2 |
           PEDSEX == 2 & SNPSEX == 1 |
           PEDSEX == 2 & SNPSEX == 0 & F > 0.4 |
           PEDSEX == 1 & SNPSEX == 0 & F < 0.8) %>%
  kable(caption = "Subjects that Fail Sex Check")
```

| FID      |  IID|  PEDSEX|  SNPSEX| STATUS  |       F|
|:---------|----:|-------:|-------:|:--------|-------:|
| 0201354  |   11|       2|       0| PROBLEM |  0.4379|
| 0201478  |   11|       2|       0| PROBLEM |  0.4276|
| 0602622  |    1|       2|       0| PROBLEM |  0.4101|
| 0700060  |    1|       2|       0| PROBLEM |  0.4101|
| 0700795  |    1|       2|       0| PROBLEM |  0.4008|
| 0801111  |    1|       2|       0| PROBLEM |  0.4147|
| 0901475  |    1|       2|       0| PROBLEM |  0.4008|
| 1300076  |   11|       2|       0| PROBLEM |  0.4658|
| 1300792  |   11|       2|       0| PROBLEM |  0.4472|
| 1300792R |   11|       2|       0| PROBLEM |  0.4472|
| 1400078  |   11|       2|       0| PROBLEM |  0.4797|
| 1400079  |   11|       2|       0| PROBLEM |  0.4472|
| 1400397  |    1|       2|       0| PROBLEM |  0.4426|

``` r
sexcheck %>%
  filter(PEDSEX == 1 & SNPSEX == 2 |
           PEDSEX == 2 & SNPSEX == 1 |
           PEDSEX == 2 & SNPSEX == 0 & F > 0.4 |
           PEDSEX == 1 & SNPSEX == 0 & F < 0.8) %>%
  dplyr::select(FID, IID) %>%
  write_tsv(path = "fail-updated-sex-check.txt", col_names = FALSE)
```

These are 13 female samples, that have out of range but quasi-reasonable F. I will remove. Disucssed with Talin.

#### Update Sex in PLINK

``` bash
plink \
--bfile ../1.split_cohorts/cohort_split_ichip1t6 \
--update-sex update-missing-sex.txt \
--remove fail-updated-sex-check.txt \
--make-bed \
--out temp1_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp1_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile ../1.split_cohorts/cohort_split_ichip1t6
    ##   --make-bed
    ##   --out temp1_all_cohort_split_ichip1t6
    ##   --remove fail-updated-sex-check.txt
    ##   --update-sex update-missing-sex.txt
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9971 people (4912 males, 5054 females, 5 ambiguous) loaded from .fam.
    ## Ambiguous sex IDs written to temp1_all_cohort_split_ichip1t6.nosex .
    ## --update-sex: 5 people updated.
    ## --remove: 9958 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9958 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp1_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate in remaining samples is 0.99889.
    ## 133250 variants and 9958 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp1_all_cohort_split_ichip1t6.bed +
    ## temp1_all_cohort_split_ichip1t6.bim + temp1_all_cohort_split_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

--update-sex: 5 people updated. --Removed 13 samples. Went from 9971 to 9958 subjects. All F

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
  write_tsv("fail-IBD-QC.tsv", col_names = FALSE)
```

``` bash
plink \
--bfile temp1_all_cohort_split_ichip1t6 \
--remove fail-IBD-QC.tsv \
--make-bed \
--out temp2_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp2_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp1_all_cohort_split_ichip1t6
    ##   --make-bed
    ##   --out temp2_all_cohort_split_ichip1t6
    ##   --remove fail-IBD-QC.tsv
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9958 people (4914 males, 5044 females) loaded from .fam.
    ## --remove: 9804 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9804 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp2_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate in remaining samples is 0.998898.
    ## 133250 variants and 9804 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp2_all_cohort_split_ichip1t6.bed +
    ## temp2_all_cohort_split_ichip1t6.bim + temp2_all_cohort_split_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

9958 people to 9804 people. 154 duplicates removed.

### Failed SNPs missingness &gt;97% and MAF 1%, and HWE 10^-6

``` bash
plink \
--bfile temp2_all_cohort_split_ichip1t6 \
--geno 0.03 \
--maf 0.01 \
--hwe 0.000001 \
--make-bed \
--out temp3_all_cohort_split_ichip1t6 
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp3_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp2_all_cohort_split_ichip1t6
    ##   --geno 0.03
    ##   --hwe 0.000001
    ##   --maf 0.01
    ##   --make-bed
    ##   --out temp3_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 133250 variants loaded from .bim file.
    ## 9804 people (4856 males, 4948 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9804 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 48 het. haploid genotypes present (see
    ## temp3_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate is 0.998898.
    ## 514 variants removed due to missing genotype data (--geno).
    ## Warning: --hwe observation counts vary by more than 10%, due to the X
    ## chromosome.  You may want to use a less stringent --hwe p-value threshold for X
    ## chromosome variants.
    ## --hwe: 10135 variants removed due to Hardy-Weinberg exact test.
    ## 15144 variants removed due to minor allele threshold(s)
    ## (--maf/--max-maf/--mac/--max-mac).
    ## 107457 variants and 9804 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp3_all_cohort_split_ichip1t6.bed +
    ## temp3_all_cohort_split_ichip1t6.bim + temp3_all_cohort_split_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

514 variants removed due to missing genotype data (--geno). 10135 variants removed due to Hardy-Weinberg exact test (--hwe) 15144 variants removed due to minor allele threshold(s)(--maf)

I am not applying a sample level genotype filter at this stage, becuase there may be some batch differences that I frist want to resolve at the SNP level rather than at the sample level. See below

### Failed Sample missingness &gt;97%

``` bash
plink \
--bfile temp3_all_cohort_split_ichip1t6 \
--mind 0.03 \
--make-bed \
--out temp4_all_cohort_split_ichip1t6 
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to temp4_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp3_all_cohort_split_ichip1t6
    ##   --make-bed
    ##   --mind 0.03
    ##   --out temp4_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 107457 variants loaded from .bim file.
    ## 9804 people (4856 males, 4948 females) loaded from .fam.
    ## 0 people removed due to missing genotype data (--mind).
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9804 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 45 het. haploid genotypes present (see
    ## temp4_all_cohort_split_ichip1t6.hh ); many commands treat these as missing.
    ## Total genotyping rate is 0.999824.
    ## 107457 variants and 9804 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to temp4_all_cohort_split_ichip1t6.bed +
    ## temp4_all_cohort_split_ichip1t6.bim + temp4_all_cohort_split_ichip1t6.fam ...
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

No samples are removed at this stage and overall genotyping rate is quite good. Continuing on with this sub-set of SNPs.

### Remove samples that fail het check

``` bash
plink \
--bfile temp4_all_cohort_split_ichip1t6 \
--remove fail-het-outlier.txt \
--make-bed \
--out qc_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile temp4_all_cohort_split_ichip1t6
    ##   --make-bed
    ##   --out qc_all_cohort_split_ichip1t6
    ##   --remove fail-het-outlier.txt
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 107457 variants loaded from .bim file.
    ## 9804 people (4856 males, 4948 females) loaded from .fam.
    ## --remove: 9733 people remaining.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9733 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 45 het. haploid genotypes present (see qc_all_cohort_split_ichip1t6.hh
    ## ); many commands treat these as missing.
    ## Total genotyping rate in remaining samples is 0.999824.
    ## 107457 variants and 9733 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --make-bed to qc_all_cohort_split_ichip1t6.bed +
    ## qc_all_cohort_split_ichip1t6.bim + qc_all_cohort_split_ichip1t6.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

9804 people to 9733 people. 71 samples out of 75 identified (3 previosuly removed).

Part 3 Re-Check File Post QC
============================

Recheck Missingess and Sex with Cleaned File
--------------------------------------------

``` bash
plink \
--bfile qc_all_cohort_split_ichip1t6 \
--missing \
--out qc_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile qc_all_cohort_split_ichip1t6
    ##   --missing
    ##   --out qc_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 107457 variants loaded from .bim file.
    ## 9733 people (4830 males, 4903 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9733 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 45 het. haploid genotypes present (see qc_all_cohort_split_ichip1t6.hh
    ## ); many commands treat these as missing.
    ## Total genotyping rate is 0.999824.
    ## --missing: Sample missing data report written to
    ## qc_all_cohort_split_ichip1t6.imiss, and variant-based missing data report
    ## written to qc_all_cohort_split_ichip1t6.lmiss.

``` r
lmiss <- read_table("qc_all_cohort_split_ichip1t6.lmiss")
lmiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_histogram(bins = 10) +
  labs(title = "SNP Level Missingess after SNP level Filtering")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-29-1.png)

``` r
imiss <- read_table("qc_all_cohort_split_ichip1t6.imiss")
imiss %>%
  ggplot(aes(x = F_MISS)) +
  geom_histogram(bins = 100) +
  labs(title = "Sample Level Missingess after Sample Level Filtering")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-30-1.png) Batch effect resolved.

Re-Check Hets
-------------

``` bash
plink \
--bfile qc_all_cohort_split_ichip1t6 \
--het \
--out qc_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile qc_all_cohort_split_ichip1t6
    ##   --het
    ##   --out qc_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 107457 variants loaded from .bim file.
    ## 9733 people (4830 males, 4903 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9733 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 45 het. haploid genotypes present (see qc_all_cohort_split_ichip1t6.hh
    ## ); many commands treat these as missing.
    ## Total genotyping rate is 0.999824.
    ## 107457 variants and 9733 people pass filters and QC.
    ## Note: No phenotypes present.
    ## --het: 106880 variants scanned, report written to
    ## qc_all_cohort_split_ichip1t6.het .

``` r
imiss <- read_table("qc_all_cohort_split_ichip1t6.imiss")
het <- read_table("qc_all_cohort_split_ichip1t6.het") %>%
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

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-32-1.png) Looks great.

Calculate Overall project MAF Adter Filtering
---------------------------------------------

``` bash
plink \
--bfile qc_all_cohort_split_ichip1t6 \
--freq \
--out qc_all_cohort_split_ichip1t6
```

    ## PLINK v1.90b5.4 64-bit (10 Apr 2018)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to qc_all_cohort_split_ichip1t6.log.
    ## Options in effect:
    ##   --bfile qc_all_cohort_split_ichip1t6
    ##   --freq
    ##   --out qc_all_cohort_split_ichip1t6
    ## 
    ## 128908 MB RAM detected; reserving 64454 MB for main workspace.
    ## 107457 variants loaded from .bim file.
    ## 9733 people (4830 males, 4903 females) loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 9733 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Warning: 45 het. haploid genotypes present (see qc_all_cohort_split_ichip1t6.hh
    ## ); many commands treat these as missing.
    ## Total genotyping rate is 0.999824.
    ## --freq: Allele frequencies (founders only) written to
    ## qc_all_cohort_split_ichip1t6.frq .

``` r
maffreq <- read_table2("qc_all_cohort_split_ichip1t6.frq")
maffreq %>%
  ggplot(aes(x = MAF)) +
  geom_histogram(aes(y =..density..)) +
  geom_density(col=2) +
  labs(title = "Overall MAF After Filtering")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-34-1.png)

Identification of Duplicated or Related Individuals
---------------------------------------------------

-   Prune dataset for temporary use of calculating cryptic relatedness and PCA as both work best under assumption of no LD among SNPs
-   Prior to calcuating identity by state, IBS, prune SNPs to only independent SNPs and remove regions with extended linkage disequiblibirum such as HLA region. The below removes snps within a 50kb window, with an r2 &gt; .2 and variant count to shift the window at the end of each step of 5kb.

``` bash
plink \
--bfile qc_all_cohort_split_ichip1t6 \
--exclude ../../original_data/highLDregions.txt \
--range \
--indep 50 5 1.8 \
--out qc_all_cohort_split_ichip1t6 &>/dev/null
```

``` bash
plink \
--bfile qc_all_cohort_split_ichip1t6 \
--extract qc_all_cohort_split_ichip1t6.prune.in \
--genome \
--min 0.12 \
--out qc_all_cohort_split_ichip1t6 &>/dev/null
```

``` r
genome <- read_table2("qc_all_cohort_split_ichip1t6.genome")

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

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-37-1.png)

Replicates (purple) removed

``` r
genome %>%
  ggplot(aes(x = as.double(PI_HAT), fill = color)) +
  geom_histogram(bins = 100) +
  labs(title = "Distribution of PI_HAT for Related Individuals >0.12", x = "PI_HAT", y = "Count (non-unique)")
```

![](ichip1t6_all_qc_files/figure-markdown_github/unnamed-chunk-38-1.png)

``` r
file.remove(list.files(pattern = "^temp", full.names = TRUE))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [29] TRUE TRUE TRUE TRUE
