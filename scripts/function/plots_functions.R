## Curve plot: 

curvePlot <- function(results_long_all, title, only_male = FALSE, male_ID){
  #results_long_all <- results_long_all[which(results_long_all$interest == TRUE),]
  if (only_male == TRUE){
    results_long_all <- results_long_all[which(results_long_all$SampleID %in% male_ID),]
  }
  
  plot <- ggplot(results_long_all, aes(x = index, y = cum, color = as.factor(threshold), linetype = as.factor(theta))) + 
    geom_line() + theme_minimal() + ggtitle(title) + xlab("Gene x sample sorted by pval") + ylab("proportion of positive")
  return(plot)
}

n50Plot <- function(results_long_all, plot_type = "boxplot"){
  # Extract the metrics list
  metrics_list <- results_long_all$metric
  
  # Convert the list to a data frame
  metric_df <- do.call(rbind, lapply(metrics_list, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  
  if (plot_type == "boxplot"){
    # Plot with ggplot2
    plot1 <- ggplot(metric_df, aes(x = threshold, y = N50)) +
      geom_boxplot(aes(group = threshold)) +  # Group by threshold to make boxplots for each threshold
      scale_x_log10() +  # Log scale for the x-axis 
      labs(
        x = "Expression Threshold (log scale)",
        y = "N50",
        title = "N50 vs. Threshold (Log Scale), Boxplots for Theta"
      ) +
      theme_minimal()
    
    
    # Plot with ggplot2
    plot2 <- ggplot(metric_df, aes(x = theta, y = N50)) +
      geom_boxplot(aes(group = theta)) +  # Group by threshold to make boxplots for each threshold
      scale_x_log10() +  # Log scale for the x-axis 
      labs(
        x = "Theta (log scale)",
        y = "N50",
        title = "N50 vs. Theta (Log Scale), Boxplots for Threshold"
      ) +
      theme_minimal()
    
  } else {
    plot1 <- ggplot(metric_df, aes(x = threshold, y = N50, color = as.factor(theta), group = theta)) +
      geom_line() +  # Line for each theta
      geom_point(size = 2) + # Optional: Add points to each line for clarity
      scale_x_log10() +  # Log scale for the x-axis
      labs(
        x = "Expression Threshold (log scale)",
        y = "N50",
        title = "N50 vs. Threshold, Lines for Each Theta",
        color = "Theta"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
    
    
    plot2 <- ggplot(metric_df, aes(x = theta, y = N50, color = as.factor(threshold), group = threshold)) +
      geom_line() +  # Line for each theta
      geom_point(size = 2) + # Optional: Add points to each line for clarity
      scale_x_log10() +  # Log scale for the x-axis
      labs(
        x = "Dispersion Theta (log scale)",
        y = "N50",
        title = "N50 vs. Theta, Lines for Each Theta",
        color = "Threshold"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  }
  
  combined_plot <- wrap_plots(list(plot1, plot2), ncol = 2) # Change ncol to adjust layout
  
  # Display combined plot
  return(combined_plot)
}
