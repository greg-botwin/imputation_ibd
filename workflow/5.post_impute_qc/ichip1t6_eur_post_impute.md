iChip 1 - 6 European Samples Prepare for Imputation Submission
================
Translational Genomics Group
30 October, 2018

``` r
df %>%
  group_by(Chr, Genotyped) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = as.factor(Chr), y = n , fill = Genotyped)) +
  geom_col() +
  facet_grid(Genotyped ~., scales="free_y") +
  labs(title = "Number of SNPs Genotyped and Imputed per Chromosome")
```

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# total number of snps genotyped and imputed
df %>%
  group_by(Genotyped) %>%
  summarise(n = n()) %>%
  kable(caption = "Total Number of SNPs Genotyped and Imputed")
```

| Genotyped |        n|
|:----------|--------:|
| Genotyped |    23713|
| Imputed   |  5506329|

``` r
df %>%
  mutate(bin = cut(Rsq, breaks = 10)) %>%
  filter(Genotyped == "Imputed") %>%
  group_by(bin) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = bin, y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -1) +
  labs(x = "Rsq Bins", y = "Number of Markers", title = "Total Number of Imputed Markers per Rsq Bin")
```

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
chrs <- c(1,6)
x2 <- list()

for(i in chrs) {
  x <- post_impute_qc(df, i)
  x2 <- list(x2, x)
}
```

``` r
x2
```

    ## [[1]]
    ## [[1]][[1]]
    ## list()
    ## 
    ## [[1]][[2]]
    ## [[1]][[2]][[1]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-1.png)

    ## 
    ## [[1]][[2]][[2]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-2.png)

    ## 
    ## [[1]][[2]][[3]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-3.png)

    ## 
    ## [[1]][[2]][[4]]
    ## 
    ## 
    ## Table: Chr 1 Number of Imputed and Genotyped Markers
    ## 
    ## Genotyped          n
    ## ----------  --------
    ## Genotyped      11740
    ## Imputed      3058191
    ## 
    ## [[1]][[2]][[5]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-4.png)

    ## 
    ## [[1]][[2]][[6]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-5.png)

    ## 
    ## [[1]][[2]][[7]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-6.png)

    ## 
    ## [[1]][[2]][[8]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-7.png)

    ## 
    ## 
    ## 
    ## [[2]]
    ## [[2]][[1]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-8.png)

    ## 
    ## [[2]][[2]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-9.png)

    ## 
    ## [[2]][[3]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-10.png)

    ## 
    ## [[2]][[4]]
    ## 
    ## 
    ## Table: Chr 6 Number of Imputed and Genotyped Markers
    ## 
    ## Genotyped          n
    ## ----------  --------
    ## Genotyped      11973
    ## Imputed      2448138
    ## 
    ## [[2]][[5]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-11.png)

    ## 
    ## [[2]][[6]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-12.png)

    ## 
    ## [[2]][[7]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-13.png)

    ## 
    ## [[2]][[8]]

![](ichip1t6_eur_post_impute_files/figure-markdown_github/unnamed-chunk-7-14.png)
