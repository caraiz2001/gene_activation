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