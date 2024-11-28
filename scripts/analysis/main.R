## Negative Binomial Activation (NB-act):

#' This functions calculates the p-value of the observing k counts given a Negative 
#' Binomial distribution with expected value mu (assuming that the 
#' gene is inactive) and fixed theta (dispersion) parameter. 
#' One independent test is done for each gene in each sample. 


# Load outrider object containing counts, sample and gene information
out_obj <- readRDS("/data/.../ods.Rds") # Specify input path


# Filter for rarely expressed genes: < 1 RPKM in > 95% of the samples
rare_exp <- filterRarely(out_obj)


# Run NB-act on the subset of rarely expressed genes and adjust pval for multiple testing
results <- nb_act(out_obj, rare_exp$rarely_exp_rpkm, adj = "BH")

# Output of the nb_act method is a list with: 
# - gene x samples matrix containing p-values
# - gene x samples matrix containing p-adjusted
# - gene x samples matrix containing expected values (mu)
# - data frame containing relevant statistics (median of outliers per gene/sample, etc.)


