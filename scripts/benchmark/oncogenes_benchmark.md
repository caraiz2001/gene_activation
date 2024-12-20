Gene Activation Analysis: Benchmark of NB-act using oncogenes
================

# Introduction:

This analysis explores the performance of nb-act method for detecting
activation outliers, that is, genes that are usually inactive in most of
the samples (very low expression) and become active in a few samples
(high expression). For benchmarking purposes, the method will be applied
to the activation of oncogenes in patients with cancer.

In the normal state, oncogenes are typically involved in promoting cell
growth, division, and survival. When mutated or overexpressed, can
potentially cause normal cells become cancerous. Therefore, we can
expect to observe a high number of activated oncogenes in cancer
patients.

Specifically, we will use the Munich Leukemia Laboratory containing 3760
blood samples of male and female individuals with leukemia.

## Load libraries, data and functions:

``` r
library(dplyr)
library(SummarizedExperiment)
library(tibble)
library(tidyr)
library(ggplot2)
library(patchwork)
library(BiocGenerics)
```

``` r
source("/data/nasif12/home_if12/l_araiz/workspace/gene_activation/scripts/function/plots_functions.R")
source("/data/nasif12/home_if12/l_araiz/workspace/gene_activation/scripts/function/help_functions.R")
source("/data/nasif12/home_if12/l_araiz/workspace/gene_activation/scripts/function/nb_act.R")
```

## Load data

``` r
MLL <- readRDS("/s/project/mll/drop_2023feb/processed_results/aberrant_expression/v33b/outrider/leukemia_14group/ods_unfitted.Rds")
MLL <- OUTRIDER::estimateSizeFactors(MLL)


oncogenes_df <- read.table("/data/nasif12/home_if12/l_araiz/workspace/gene_act/oncogenes.tsv", sep = "\t", row.names = 1)
```

## Pre-processing

Select those genes that are rarely expressed (\< 1 RPKM in \> 95% of the
samples). We would expect oncogenes to be rarely expressed, as there are
usually inactive but become active in some cancer patients.

``` r
#Filter for rarely expressed genes (6401 genes)
rare_MLL <- filterRarely(MLL)
```

    ## count matrix has 20084 genes and 3760 samples

    ## rpkm matrix has 20084 genes and 3760 samples

    ## There are 6401 RE genes

    ##                            MLL_32448 MLL_13008 MLL_28126 MLL_10807 MLL_18412
    ## ENSG00000000003.15_4              44         3         4        18        38
    ## ENSG00000000005.6_3                0         0         0         0         0
    ## ENSG00000001626.16_8               0         0         0         0         0
    ## ENSG00000002586.20_6_PAR_Y         0         0         0         0         0
    ## ENSG00000002746.15_4               4         0         0         2         4

``` r
# There are 318 oncogenes
nrow(oncogenes_df)
```

    ## [1] 318

``` r
# 314 of them are present in the MLL dataset
length(intersect(sub("\\..*", "", rownames(oncogenes_df)), sub("\\..*", "", rownames(rare_MLL$rpkm))))
```

    ## [1] 314

``` r
# 62 of them are rarely expressed. 
length(intersect(sub("\\..*", "", rownames(oncogenes_df)), sub("\\..*", "", rownames(rare_MLL$rarely_exp_rpkm))))
```

    ## [1] 62

For benchmarking purposes, let’s see if we can find a subset of samples
with more rarely expressed oncogenes.

We can have a subset with 67 instead of 62 (not much better).

``` r
# Load your RPKM matrix and oncogene list
# rpkm: rows are genes, columns are samples
# oncogenes: vector of oncogene IDs (e.g., ENSG000xxxx)

# Define the number of samples you want to select
subset_size <- 500 # Adjust as needed

# Run the optimization
result <- findOptimalSubset(rare_MLL$rpkm, oncogenes_df, subset_size)

# Extract the best subset of samples
best_samples <- result$best_subset
max_rare_genes <- result$max_rare_genes

cat("Best subset includes", max_rare_genes, "rarely expressed oncogenes.\n")
```

    ## Best subset includes 66 rarely expressed oncogenes.

``` r
# Just to check
counts_MLL <- assays(MLL)$counts
subset_counts_MLL <- counts_MLL[, best_samples]

rare_MLL_subset <- filterRarely(MLL, subset_counts = subset_counts_MLL)
```

    ## count matrix has 20084 genes and 500 samples

    ## rpkm matrix has 20084 genes and 500 samples

    ## There are 6476 RE genes

    ##                            MLL_19288 MLL_54965 MLL_16704 MLL_12004 MLL_11284
    ## ENSG00000000003.15_4               1        16         4         1        43
    ## ENSG00000000005.6_3                0         0         0         0         0
    ## ENSG00000001626.16_8               0         0         0         0         0
    ## ENSG00000002586.20_6_PAR_Y         0         0         0         0         0
    ## ENSG00000002746.15_4               0         0         0         0         3

``` r
length(intersect(sub("\\..*", "", rownames(oncogenes_df)), sub("\\..*", "", rownames(rare_MLL_subset$rarely_exp_rpkm))))
```

    ## [1] 66

## Parameter exploration

``` r
# Parameter values to test 
threshold_vals <- c(0.1, 1, 5, 10, 50, 100)
theta_vals <- c(0.1, 1, 10, 100)

results_long_all_MLL <- exploreParamsOnco(MLL, rare_MLL_subset$rarely_exp_rpkm, theta_vals = theta_vals, threshold_vals = threshold_vals, interest_list = rownames(oncogenes_df))
```

``` r
results_long_all_MLL$results_long_all_mini <- results_long_all_MLL$results_long_all[which(results_long_all_MLL$results_long_all$interest == TRUE),]
curvePlot(results_long_all_MLL$results_long_all_mini, title = "Cumulative sum of oncogenes")
```

![](oncogenes_benchmark_files/figure-gfm/curve%20MLL%20onco-1.svg)<!-- -->

``` r
n50Plot(results_long_all_MLL, plot_type = "line")
```

![](oncogenes_benchmark_files/figure-gfm/N50%20MLL_onco-1.svg)<!-- -->

In this case, low values for threshold and/or theta work pretty bad (=
0.1, 1). Best performance is achieved with very high values for both
threshold (making it more conservative) and theta (making it less
conservative). This behaviour differs from the pattern that we were
observing in the Y-chr benchmark.

Reasons might be that activated oncogenes have a higher expression (in
absolute terms) than Y-chr genes, especially if the Y-chr gene is not
highly expressed in the tissue. That is, a better signal to noise ratio.
As a consequence, a higher threshold for describing “inactive” genes
might remove part of this “noise”.

# CONCLUSION:

Y-chr and oncogenes benchmark present different patterns in terms of
optimal parameters for modeling the distribution of inactive genes. In
the case of Y-chr (both in GTEx and MLL), lower values of threshold
(0.1, 1) and theta (\<= 10) achieve a better performance. In the case of
Oncogenes, higher values of threshold and theta are better.

In general, oncogenes do not necessarily appear in the first positions
when ranking outliers according to pvalue. Moreover, defining all
oncogenes in all samples as positives might be overstimating the number
of positive genes, as not every oncogene is expected to be activated in
every patient. For this reasons, the oncogene approach is less suited
for benchmarking purposes.

1)  Optimal parameters may depend on the type of gene that is being
    taken as a reference, rather than the dataset (Y-chr in GTEx and MLL
    is more similar than y-chr MLL and oncogenes MLL)

2)  Y-chr works better than oncogenes.
