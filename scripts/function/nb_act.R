#' Run negative binomial activation 
#'
#' This functions calculates the pvalue of the observing k counts given a Negative 
#' Binomial distribution with calculated mu (expected value assuming that the 
#' gene is inactive) and fixed theta (dispersion) parameter. 
#' (for each gene in each sample)  
#'
#' @param out_obj Outrider object containing counts, sample and gene information
#' @param rarely_exp_genes Count/RPKM matrix filtered for rarely expressed genes (>1 RPKM in <5% of the samples)
#' @param threshold Inactive expression threshold (RPKM = 5 per default) 
#' @param theta Dispersion parameter, corresponds to size in the NB distribution
#' @param adj Optionally, you can adjust for multiple testing (FDR) with adj = "BH"
#' @param sig_threshold outliers with a pval/padj lower than this value will be consider significant
#' @return list containing pval matrix, padj matrix, mu matrix, threshold, theta.
nb_act <- function(out_obj, rarely_exp_rpkm, threshold = 5, theta = 10, adj = "NO", sig_threshold = 0.05){
  # define the subset of genes and samples
  rarely_exp_genes <- rownames(rarely_exp_rpkm)
  subset_samples <- colnames(rarely_exp_rpkm)
  
  # Calculate expected values mu_ij: 
  size_factors <- as.data.frame(colData(out_obj)$sizeFactor) # Extract size factors
  
  gene_length <- as.data.frame(elementMetadata(out_obj)$basepairs)  # Extract gene length
  subset_gene_length <- gene_length[rarely_exp_genes, ,drop= FALSE] # Subset gene length
  
  # Extract raw counts: 
  counts <- as.data.frame(assays(out_obj)$counts)
  
  if (expected_calc == "org"){
    # Calculate mu
    mu_ij <- outer((threshold*subset_gene_length[,1]/1000), size_factors[, 1], "*")
    
  } else if (expected_calc == "mod"){
    #Calculate mu: 
    mu_ij <- outer((threshold*subset_gene_length[,1]/1000), (median(colSums(counts))*size_factors[,1])/1000000, "*")
    
  }
  
  rownames(mu_ij) <- rarely_exp_genes
  colnames(mu_ij) <- colnames(counts)
  mu_ij <- mu_ij[rarely_exp_genes, subset_samples, drop=FALSE]
  
  subset_counts <- counts[rownames(counts) %in% rarely_exp_genes, subset_samples , drop=FALSE]
  
  # Calculate pval
  nb_probs <- pnbinom(as.matrix(subset_counts), mu = mu_ij, size = theta, lower.tail = FALSE)
  
  colnames(nb_probs) <- colnames(subset_counts)
  rownames(nb_probs) <- rownames(subset_counts)
  
  if (adj == "BH"){
    # take all pvalues
    p_values_vector <- as.vector(nb_probs)
    
    # Adjust the p-values (e.g., using Benjamini-Hochberg)
    adjusted_p_values <- p.adjust(p_values_vector, method = "BH")
    
    # Reshape it back into the original matrix shape
    adjusted_matrix <- matrix(adjusted_p_values, nrow = nrow(nb_probs), ncol = ncol(nb_probs))
    
    colnames(adjusted_matrix) <-colnames(subset_counts)
    rownames(adjusted_matrix) <- rownames(subset_counts)
    
  } else {
    adjusted_matrix <- FALSE
  }
  
  # Summary statistics
  if (adj == "BH") {
    number_outliers <- sum(adjusted_matrix < sig_threshold)
    observed_higher_than_expected <- sum(subset_counts > mu_ij)
    activated_genes <- sum(rowSums(adjusted_matrix < sig_threshold) > 0)
    activated_samples <- sum(colSums(adjusted_matrix < sig_threshold) > 0)
    median_outliers_per_gene <- median(rowSums(adjusted_matrix < sig_threshold))
    median_outliers_per_sample <- median(colSums(adjusted_matrix < sig_threshold))
  } else {
    number_outliers <- sum(nb_probs < sig_threshold)
    observed_higher_than_expected <- sum(subset_counts > mu_ij)
    activated_genes <- sum(rowSums(nb_probs < sig_threshold) > 0)
    activated_samples <- sum(colSums(nb_probs < sig_threshold) > 0)
    median_outliers_per_gene <- median(rowSums(nb_probs < sig_threshold))
    median_outliers_per_sample <- median(colSums(nb_probs < sig_threshold))
  }
  
  
  overview <- data.frame(expected_calc = expected_calc, threshold = threshold, theta = theta,
                         number_outliers = number_outliers, observed_higher_than_expected = observed_higher_than_expected,
                         activated_genes = activated_genes,  activated_samples = activated_samples,
                         median_outliers_per_gene = median_outliers_per_gene, 
                         median_outliers_per_sample = median_outliers_per_sample)
  
  # Print summary
  message("Summary:")
  message("Expected values calculated using: ", expected_calc)
  message("Total outliers (padj or pval < 0.05): ", number_outliers)
  message("Observed > Expected: ", observed_higher_than_expected)
  message("Genes with some outlier: ",  activated_genes)
  message("Samples with some outlier: ",  activated_samples)
  message("Median of outliers per gene: ",  median_outliers_per_gene)
  message("Median of outliers per sample: ",  median_outliers_per_sample)
  
  
  list <- list("pval" = nb_probs, "padj" = adjusted_matrix, "mu" = mu_ij, "overview" = overview)
  return(list)
}