# GOOD FUNCTION SCRIPT

# Given a vector of theta and threshold, produces the long format results 

exploreParams <- function(out_obj, rarely_exp_rpkm, theta_vals, threshold_vals, 
                          interest_list, male_ID, too_large = FALSE){
  metric <- list()
  results_long_all <- data.frame()
  for (t in threshold_vals){
    for (i in theta_vals){
      print(paste0("threshold_", t, "theta_", i))
      results <- nb_act(out_obj, rarely_exp_rpkm, adj = "NO", threshold = t, theta = i)
      
      results_long <- summaryLong(results = results, interest_list = interest_list)
      results_long <- results_long$p_long
      results_long <- results_long[order(results_long$p_val_adj),]
      results_long$sex_male <- ifelse(results_long$SampleID %in% male_ID, TRUE, FALSE)
      results_long$positives <- ifelse(results_long$interest == TRUE & results_long$sex_male == TRUE, TRUE, FALSE)
      results_long$threshold<- rep(t, nrow(results_long))
      results_long$theta<- rep(i, nrow(results_long))
      results_long$cum <- cumsum(results_long$positives)
      results_long$index <- seq(1:nrow(results_long))
      
      # Determine N50 index:
      total_interest <- sum(results_long$positives)
      
      x_index <- tail(which(results_long$cum < 0.5*total_interest),1)
      method <- paste0("threshold_", t, "_theta_", i)
      
      # Save the important metrics
      metric[[method]] <- list("threshold" = t, "theta" = i, "N50" = x_index )
      
      # Leave only TRUE and join with results_long_all
      if (too_large == TRUE){
        results_long <- results_long[which(results_long$interest == TRUE),]
        results_long_all <- rbind(results_long_all, results_long)
      } else {
        results_long_all <- rbind(results_long_all, results_long)
      }
      
      results_long <- data.frame()
    }
  }
  list <- list("results_long_all" = results_long_all, "metric" = metric)
  return(list)
}


# Given a vector of theta and threshold, produces the long format results 

exploreParamsOnco <- function(out_obj, rarely_exp_rpkm, theta_vals, threshold_vals, 
                          interest_list, too_large = FALSE){
  metric <- list()
  results_long_all <- data.frame()
  for (t in threshold_vals){
    for (i in theta_vals){
      print(paste0("threshold_", t, "theta_", i))
      results <- nb_act(out_obj, rarely_exp_rpkm, adj = "NO", threshold = t, theta = i)
      
      results_long <- summaryLong(results = results, interest_list = interest_list)
      results_long <- results_long$p_long
      results_long <- results_long[order(results_long$p_val_adj),]
      results_long$positives <- ifelse(results_long$interest == TRUE, TRUE, FALSE)
      results_long$threshold<- rep(t, nrow(results_long))
      results_long$theta<- rep(i, nrow(results_long))
      results_long$cum <- cumsum(results_long$positives)
      results_long$index <- seq(1:nrow(results_long))
      
      # Determine N50 index:
      total_interest <- sum(results_long$positives)
      
      x_index <- tail(which(results_long$cum < 0.5*total_interest),1)
      method <- paste0("threshold_", t, "_theta_", i)
      
      # Save the important metrics
      metric[[method]] <- list("threshold" = t, "theta" = i, "N50" = x_index )
      
      # Leave only TRUE and join with results_long_all
      if (too_large == TRUE){
        results_long <- results_long[which(results_long$interest == TRUE),]
        results_long_all <- rbind(results_long_all, results_long)
      } else {
        results_long_all <- rbind(results_long_all, results_long)
      }
      
      results_long <- data.frame()
    }
  }
  list <- list("results_long_all" = results_long_all, "metric" = metric)
  return(list)
}

## Select subset of females (95%, males 5%)

# input OUTRIDER OBJECT
# params number of males, number females, column_name 
# output subset_counts

subsetRandom <- function(out_obj, num_males, num_females, groups = FALSE, seed = 123){
  counts <- as.data.frame(assays(out_obj)$counts)
  
  if (identical(groups, FALSE)){
    male_counts <- counts[,which(colData(out_obj)$SEX == "Male")] 
    female_counts <- counts[,which(colData(out_obj)$SEX == "Female")] 
  } else {
    male_counts <- counts[,which(colData(out_obj)[[groups[1]]] == groups[2])] 
    female_counts <- counts[,which(colData(out_obj)[[groups[1]]] == groups[3])] 
  }
  
  # Set a seed for reproducibility
  set.seed(seed)  # You can choose any integer
  
  random_males <- male_counts[,sample(ncol(male_counts), num_males) ]
  random_females <- female_counts[,sample(ncol(female_counts), num_females) ]
  
  subset_counts <- cbind(random_females, random_males)
  message("Resulting count matrix has ", nrow(subset_counts), " genes and ", ncol(subset_counts), " samples")
  
  return(subset_counts)
}

# Individual stats
summaryLong <- function(results, summary_of = "pval", sig_threshold = 0.05, interest_list){
  if (summary_of == "pval"){
    p <- as.data.frame(results$pval)
  } else if (summary_of == "padj") {
    p <- as.data.frame(results$padj)
  }
  
  # Convert to long format
  p <- rownames_to_column(as.data.frame(p), var = "Gene")
  
  
  p_long <-  p %>%
    pivot_longer(
      cols = -Gene,              # All columns except "Gene"
      names_to = "SampleID",     # New column for sample IDs
      values_to = "p_val_adj"         # New column for p-values
    )
  
  # Add Significance
  p_long$interest <- ifelse(sub("\\..*", "", p_long$Gene) %in% sub("\\..*", "", interest_list) , TRUE, FALSE)
  p_long$sig <- ifelse(p_long$p_val_adj < sig_threshold, TRUE, FALSE)
  
  # Individual stats
  total_interest <- nrow(p_long[which(p_long$interest == TRUE),])
  true_positives <- nrow(p_long[which(p_long$sig == TRUE & p_long$interest == TRUE), ])
  false_positives <- nrow(p_long[which(p_long$sig == TRUE & p_long$interest == FALSE),])
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / total_interest
  
  stats <- list(total_interest = total_interest, true_positives = true_positives,
                false_positives = false_positives, precision = precision, 
                recall = recall)
  
  # print message
  message("Outlier-wise stats:")
  message("There are ", total_interest, " gene x sample in the set of interest (potential outliers)")
  message(true_positives, " of are activated")
  message("Precision (TP / TP + FP): ", precision)
  message("Recall (TP / TP + FN): ", recall)
  
  list <- list("stats" = stats, "p_long" = p_long)
  return(list)
}


#' Filter count data frame to contain only rarely expressed genes:
#'
#' Selects those genes that are inactive (<1 FPKM) in more than 95% of the samples
#' They are the ones to be analysed by nb-act 
#'
#' @param subset_counts Optionally, you can enter a subset (of samples), for example selected by gender
#' @param out_obj Outrider object containing counts, sample and gene information
#' @return list containing the rpkm (all genes) and a filtered df with rpkm and 
#' counts of rarely expressed subset
filterRarely <- function(out_obj, subset_counts = FALSE){
  # extract counts
  if (identical(subset_counts, FALSE)) {
    counts <- as.data.frame(assays(out_obj)$counts)
  } else {
    counts <- subset_counts
  }
  
  message("count matrix has ", nrow(counts), " genes and ", ncol(counts), " samples")
  
  # Convert to RPKM
  gene_length_bp <- elementMetadata(out_obj)$basepairs # Extract gene length (in bp)
  gene_length_kb <- gene_length_bp/1000 # Convert to kb
  
  total_reads_sample <- colSums(as.matrix(counts)) # Obtain total reads per sample 
  
  rpkm <- t(t(counts/gene_length_kb)/total_reads_sample*1e6) # apply RPKM formula
  
  message("rpkm matrix has ", nrow(rpkm), " genes and ", ncol(rpkm), " samples")
  
  # Define rarely expressed using outrider criteria: 
  rarely_exp_rpkm <- rpkm[which(apply(rpkm, 1, common_exp) == FALSE),]
  rarely_exp_genes <- rownames(as.data.frame(rarely_exp_rpkm))
  message("There are ", length(rarely_exp_genes), " RE genes ")
  
  # Subset counts for rarely expressed genes
  subset_counts <- counts[rarely_exp_genes,, drop=FALSE]
  print(subset_counts[1:5,1:5])
  list <- list("rpkm" = rpkm, "rarely_exp_rpkm" = rarely_exp_rpkm , "subset_counts" = subset_counts, "rarely_exp_genes" = rarely_exp_genes)
  return(list)
}

#' Select rarely expressed genes following outrider criteria 
#'
#' Given a vector of rpkm values (for one gene across samples), determines if the
#' gene is commonly or rarely expressed following OUTRIDER criteria
#'
#' @param rpkm_row vector of rpkm values (for one gene across samples)
#' @return FALSE if the gene is rarely expressed
#' counts of rarely expressed subset
common_exp <- function(rpkm_row){
  samples <- length(rpkm_row)
  n <- 0.05*samples
  exp_threshold <- 1
  sum(rpkm_row > exp_threshold) > n 
}


# Function to calculate rarely expressed genes in a subset of samples
rarelyExpressedNumber <- function(rpkm_subset, threshold = 1, rare_threshold = 0.05) {
  samples <- ncol(rpkm_subset)
  max_expressed <- rare_threshold * samples
  # Count samples with RPKM > threshold for each gene
  rarely_exp <- rowSums(rpkm_subset > threshold) <= max_expressed
  return(rarely_exp) # Returns a logical vector indicating rarely expressed genes
}


# Function to optimize sample selection
findOptimalSubset <- function(rpkm, oncogenes, subset_size) {
  oncogene_rpkm <- rpkm[sub("\\..*", "", rownames(rpkm)) %in% sub("\\..*", "", rownames(oncogenes_df)), ] # Filter to oncogenes

  best_subset <- NULL
  max_rare_genes <- 0
  
  # Iterate over random subsets of samples to find the best subset
  for (i in 1:2000) { # Number of iterations
    sample_subset <- sample(colnames(oncogene_rpkm), size = subset_size) # Random sample
    subset_rpkm <- oncogene_rpkm[, sample_subset] # Subset RPKM matrix
    
    rare_genes <- sum(rarelyExpressedNumber(subset_rpkm)) # Count rarely expressed genes
    if (rare_genes > max_rare_genes) {
      max_rare_genes <- rare_genes
      best_subset <- sample_subset
    }
  }
  
  return(list("best_subset" = best_subset, "max_rare_genes" = max_rare_genes))
}




