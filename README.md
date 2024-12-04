# Gene Activation Analysis
This project is focused on the development and benchmark of NB-act for detecting activation outliers in RNA-seq data. 

Previously, we developed OUTRIDER, a denoising-autoencoder based method for detecting gene 
expression outliers in RNA-seq data (https://doi.org/10.1016/j.ajhg.2018.10.025). 
OUTRIDER has been successfully applied to cohorts of rare-disease patients and helped to pinpoint causal genes. 
However, OUTRIDER can currently only fit genes that are generally expressed across samples (> 1 FPKM in > 5% of the samples). 

In cancer, genes that are usually not expressed but become activated in cancer cells can 
be potential oncogenes and are therefore of interest. Not expressed or activated genes can be pathogenic (as oncogenes developmental genes). 
Since then, we started exploring activation itself and we published a prototype method of NB-act in our leukemia paper. 
The goal is to benchmark and compare it against other emerging method by Vanderstichele et al and test whether we can improve it. 
The final goal is to have a working version of the method in our Detection of RNA Outliers Pipeline (https://github.com/gagneurlab/drop).


## NB-act: 
NB-act (Negative Binomial activation) focuses on the analysis of rarely expressed genes, which can not be analysed with OUTRIDER due
to statistical limitations.

This method calculates p-values for observed counts for each gene in each sample 
under the null hypothesis that the gene is **not expressed** in the sample. 
Specifically, NB-act computes the probability of observing a certain number of fragments or more for a gene in a particular sample, 
assuming a negative binomial distribution with a fixed expected baseline expression and a fixed dispersion parameter of 10. 

The original prototype was published with expected baseline expression = 1 FPKM, 
where 1 corresponds to the threshold separating expressed from non-expressed genes in OUTRIDER. 
The dispersion parameter was originally 0.02, corresponds to the empirically observed 
lowest dispersion values estimated by OUTRIDER on expressed genes. 
As low dispersion corresponds to high variance, we chose a low dispersion value for 
NB-act to be conservative.

The optimal values for the baseline expression (= expression threshold) and dispersion (= theta) have been modified
during this benchmark (1 FPKM and 10). 

Specifically, NB-act works as follows: 

1. Filter for rarely expressed genes using OUTRIDER criteria (> 1 RPKM in less than 5% of the samples).

2. Calculate expected counts corresponding to 1 RPKM (baseline expression of **inactive** genes) considering gene length, library size and size factors. 
This results in different expected counts for each gene in each sample and consequently, **each gene in each sample has a different distribution**. 

The expected counts (mu_ij), in the scale of raw counts, can be calculated using the following ecuation:


$$
\mu_{ij} (\text{reads}) = \frac{1 \, \text{RPKM} \left( \frac{\text{reads} \cdot 10^6}{\text{kb}} \right) \cdot \text{gene length (kb)} \cdot \text{median library size (reads)} \cdot S_j}{1000000 \, (\text{reads})}
$$

3. Calculate the probability of the observed counts (each gene in each sample) assuming a Negative Binomial distribution
with parameters mu_ij, theta.

4. Correct pvalue for multiple testing. 


Conclusion: we assume that the gene is inactive (<= 1 RPKM, **expression threshold**), we recover the absolute counts corresponding 1 RPKM 
for that gene in that particular sample (**expected counts**) and we test if the **observed counts** (k_ij) are significantly higher (= activation). 
Significance is assesed based on the p-value corresponding to a Negative Binomial distribution (k_ij ~ NB(mu_ij, theta)). 


## Project structure: 
The repository of this project is organized as follows: 

1. **Analysis**: main analysis (NB-act function). 
  - Input: the OUTRIDER dataset containing count matrix and gene/sample information (gene length, Sex, etc.). 
  
  - analysis: filtering for rarely expressed genes and running nb-act
  
  - Output: output is a list containing pvalues, padj, expected values (all 3 are genes x samples matrix ) and a data frame containing 
  relevant statistics (median of outliers per gene/sample, etc.)

2. **Benchmark**: Several markdown documents explaining different benchmark approaches: 
 
  - Y-chromosome activation: select a subset with a 95% females and 5% males (so that genes in the Y-chr are rarely expressed)
  and see if NB-act identifies Y-chr genes as outliers.
  
  - Oncogene activation: using a dataset of cancer patients (MLL), test if activation outliers are enriched for oncogenes. 
  
3. **Functions**: reusable functions for analysis and benchmark. 

  - nb_act: main NB-act function
  
  - help_functions: accessory functions for fitlering, interpreting the results, etc.
  
  - plot_functions: functions for specific plots. 
