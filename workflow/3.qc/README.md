Perform Ancestry Specific QC - README
================
Perform Ancestry Specific QC
2/22/2019

Tasks
-----

A consistent QC is performed across all 5 of our cohorts (e.g. ichip123456\_eur, ichip123456\_all, ichip7\_eur, ichip7\_all, bbc). The QC is broken down into three steps. - Part 1 Identifies QC issues - Part 2 Fixes and applies QC parameters - Part 3 Re-Checks genotype data post QC Changes

The QC thresholds applied are as follows:

1.  Sex Checks were perfomed
    -   If PED sex 0 and F &lt; 0.4 updated PED to Female 2
    -   If PED sex 0 and F &gt; 0.8 update PED to male 1 e
    -   samples that have out of range F were removed using F &lt; 0.4 = female and F &gt; 0.8 male.
2.  Identical samples were removed
    -   all pairs of individuals with an IBD &gt;=0.8 identified
    -   individual in the pair with the lower genotyping rate removed
3.  Markers were filtered for missigness, maf and HWE
    -   --geno 0.03
    -   --maf 0.01
    -   --hwe 0.000001
4.  Samples were filtered for missigness
    -   --mind 0.03
5.  Samples with outlying heterozygoisty were filterd
    -   if cohort &gt; 1000 samples uses `mean(imiss_het$obs_het_rate)+(4*sd(imiss_het$obs_het_rate)`
    -   if cohort &lt; 1000 samples uses `mean(imiss_het$obs_het_rate)+(3*sd(imiss_het$obs_het_rate)`
