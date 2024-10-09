library(matrixcalc)

#Write where you would like your output
#and .cpp file has to be in working directory
wd<-"C:/Users/fischn/Dropbox/"
wd<-"C:/Users/Derek.Chamberlin/Work/Research/Age_Err_Simulation/Assessment_AgeingError"

load("workspace.RData")

#Check Scenarios for Hessian Positive/Definite
hessian_check <- list()

for (i in 1:length(res_list_final)) {
  hessian_check[[i]] <- list()
  
  for (j in 1:length(res_list_final[[i]])) {
    hessian <- res_list_final[[i]][[j]]$hessian
    
    # Check if the Hessian is positive definite
    if (!is.null(hessian)) {
      is_positive_definite <- is.positive.definite(hessian)
    } else {
      is_positive_definite <- FALSE
    }
    
    hessian_check[[i]][[j]] <- is_positive_definite
  }
}

summary_percent_true <- matrix(data = NA, nrow = length(hessian_check), ncol = 4)

for (i in 1:length(hessian_check)) {
  results <- hessian_check[[i]]
  
  # Calculate the percentage of TRUE values
  num_results <- length(results)
  num_true <- sum(unlist(results) == TRUE, na.rm = TRUE)

  percent_true <- (num_true / num_results) * 100

  summary_percent_true[i, 4] <- percent_true
}

summary_percent_true[, 1] <- scenarios[, 1]
summary_percent_true[, 2] <- scenarios[, 2]
summary_percent_true[, 3] <- scenarios[, 3]

print(summary_percent_true)
