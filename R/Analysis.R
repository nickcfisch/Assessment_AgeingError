library(matrixcalc)

#Write where you would like your output
#and .cpp file has to be in working directory
wd<-"C:/Users/fischn/Dropbox/"
#wd<-"C:/Users/Derek.Chamberlin/Work/Research/Age_Err_Simulation/Assessment_AgeingError/"

setwd(wd)

load("workspace.RData")
source(paste0(wd,"/R/Analysis_Functions.R"))

#Check Scenarios for Hessian Positive/Definite
hessian_check <- list()

for (i in 1:length(res_list_final)) {
  hessian_check[[i]] <- list()
  
  for (j in 1:length(res_list_final[[i]])) {
    hessian <- res_list_final[[i]][[j]]$hessian
    
    # Check if the Hessian is positive definite
    if (!is.null(hessian)) {
      #is_positive_definite <- is.positive.definite(hessian)
      is_positive_definite <- is.positive.definite(hessian) & res_list_final[[i]][[j]]$max_gradient <= 0.1
    } else {
      is_positive_definite <- FALSE
    }
    
    hessian_check[[i]][[j]] <- is_positive_definite
  }
}

summary_percent_true <- matrix(data = NA, nrow = length(hessian_check), ncol = 5)

for (i in 1:length(hessian_check)) {
  results <- hessian_check[[i]]
  
  # Calculate the percentage of TRUE values
  num_results <- length(results)
  num_true <- sum(unlist(results) == TRUE, na.rm = TRUE)

  percent_true <- (num_true / num_results) * 100

  summary_percent_true[i, 4] <- percent_true
}

# Initialize a matrix to store the mean convcounter values
mean_convcounter <- numeric(length(res_list_final))

for (i in 1:length(res_list_final)) {
  convcounter_values <- numeric(length(res_list_final[[i]]))
  
  for (j in 1:length(res_list_final[[i]])) {
    # Check if convcounter exists and is numeric
    if (!is.null(res_list_final[[i]][[j]]$convcounter) && is.numeric(res_list_final[[i]][[j]]$convcounter)) {
      convcounter_values[j] <- res_list_final[[i]][[j]]$convcounter
    }
  }
  
  # Calculate the mean, excluding NA values
  mean_convcounter[i] <- mean(convcounter_values, na.rm = TRUE)
}

# Add mean_convcounter to the summary_percent_true matrix
summary_percent_true[, 5] <- mean_convcounter

summary_percent_true[, 1] <- scenarios[, 1]
summary_percent_true[, 2] <- scenarios[, 2]
summary_percent_true[, 3] <- scenarios[, 3]

print(summary_percent_true)





#Time series plot of SSB relative error
load("./Output/GT_OM_perf_wdat.RData")
OM_perf <- OM_wdat
rm(OM_wdat)
perf_re <- relative_error(OM_perf,res_list_final[[1]])
perf_re_plot <- plot_re(perf_re)
print(perf_re_plot)