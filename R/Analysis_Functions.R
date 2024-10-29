relative_error <- function(OM, EM){
  relative_error_all <- matrix(data = NA, nrow = length(OM), ncol = 70)
  for (i in 1:length(OM)) {
    true_values <- OM[[i]]$OM$SSB[26:95]
    estimate <- EM[[i]]$SD$value[1:70]
    relative_error_i <- (true_values - estimate) / true_values
    relative_error_all[i,] <- relative_error_i
  }
  return(relative_error_all)
}

plot_re <- function(re_df){
  library(ggplot2)
  library(reshape2)
  
  mean_re <- colMeans(re_df, na.rm = TRUE)
  std_error <- apply(re_df, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  
  mean_perf_re_df <- data.frame(
    Index = 1:length(mean_re),
    MeanRelativeError = mean_re,
    LowerCI = mean_re - (1.96 * std_error),
    UpperCI = mean_re + (1.96 * std_error)
  )
  
  figure <- ggplot(mean_perf_re_df, aes(x = Index, y = MeanRelativeError)) +
    geom_line() +
    geom_point() +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) +
    labs(title = "Mean Relative Error with 95% Confidence Intervals",
         x = "Index",
         y = "Mean Relative Error") +
    theme_minimal()
  
  return(figure)
}