library(matrixcalc)

#To Do:
# add in code to calculate F-ratio and B-ratio
# change plot_re function to things other than SSB to be plotted
# change grid_plot to allow y axes to be changed


#Write where you would like your output
#and .cpp file has to be in working directory
wd<-"C:/Users/Derek.Chamberlin/Work/Research/Age_Err_Simulation/Assessment_AgeingError/"
wd<-"C:/Users/fischn/Documents/GitHub/Assessment_AgeingError"

setwd(wd)

load("./Output/workspace.RData")
load("./Output/GT_OM_perf_wdat.RData")
GT_OM_perf_wdat <- OM_wdat
load("./Output/GT_OM_constant_wdat.RData")
GT_OM_constant_wdat <- OM_wdat
load("./Output/GT_OM_linear_wdat.RData")
GT_OM_linear_wdat <- OM_wdat
load("./Output/GT_OM_curvilinear_wdat.RData")
GT_OM_curvilinear_wdat <- OM_wdat
rm(OM_wdat)
wd<-"C:/Users/fischn/Documents/GitHub/Assessment_AgeingError"
source(paste0(wd,"/R/Analysis_Functions.R"))
source(paste0(wd,"/R/Functions.R"))

#Check Scenarios for Hessian Positive/Definite
hessian_check <- list()

for (i in 1:length(res_list_final)) {
  hessian_check[[i]] <- list()
  
  for (j in 1:length(res_list_final[[i]])) {
    hessian <- res_list_final[[i]][[j]]$hessian
    
    # Check if the Hessian is positive definite
    if (!is.null(hessian)) {
      is_positive_definite <- is.positive.definite(hessian)
      #is_positive_definite <- is.positive.definite(hessian) & res_list_final[[i]][[j]]$max_gradient <= 0.1
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



#Calculate RE in SSB
SSB_re <- list()
for (i in 1:nrow(scenarios)) {
  SSB_re[[i]] <- relative_error(get(scenarios[i , 2]),res_list_final[[i]]) #function on runs on hessian PD iterations
}


#Time series boxplots of RE in SSB
OM_Title <- c(rep("No AE", 4), rep("AE Constant Bias", 4), rep("AE Linear Bias", 4), rep("AE Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "AE Constant Bias", "AE Linear Bias", "AE Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(4.1, 4.6, 2.1, 1.1))
for (i in 1:nrow(scenarios)) {
  boxplot(SSB_re[[i]], ylim = c(-1.0, 1.0), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
  axis(1, at = 1, labels = 1, cex.axis = 1.1)
  axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1)
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = -2, cex = 1.3, outer = TRUE)
mtext("Relative Error in SSB", side = 2, line = -2, cex = 1.3, outer = TRUE)



#Time series plots of RE in SSB
grid_plot(SSB_re)




#F-ratio
#Finding fmsy for the true operating model
MSY <- find_msy()  #Looks to be 0.349

MSY_re <- vector("list", length(res_list_final))
for (k in 1:length(res_list_final)) {
  #add in if/else statement to skip over nonconverged iteration where res_list_final[[k]][[j]]$SD$par doesn't exist
  MSY_re[[k]] <- vector("list", length(res_list_final[[1]]))
  for (j in 1:length(res_list_final[[1]])) {
    hessian <- res_list_final[[k]][[j]]$hessian
    
    # Check if the Hessian is positive definite
    if (!is.null(hessian)) {
      MSY_re[[k]][[j]] <- find_msy(Mref=exp(res_list_final[[k]][[j]]$SD$par.fixed["log_M"]),                 
                      B1=res_list_final[[k]][[j]]$SD$par.fixed["B1"],                   
                      B2=res_list_final[[k]][[j]]$SD$par.fixed["B2"],                         
                      B3=res_list_final[[k]][[j]]$SD$par.fixed["B3"],
                      B4=res_list_final[[k]][[j]]$SD$par.fixed["B4"],
                      R0=exp(res_list_final[[k]][[j]]$SD$par.fixed["log_R0"]),
                      h=0.4593,
                      sd_rec=res_list_final[[k]][[j]]$SD$par.fixed["log_sigma_rec"])  
      MSY_re[[k]][[j]]$fratio_EM <- exp(res_list_final[[k]][[j]]$SD$par.fixed[9:77])/MSY_re[[k]][[j]]$fmsy 
      MSY_re[[k]][[j]]$bratio_EM <- res_list_final[[k]][[j]]$SD$value/MSY_re[[k]][[j]]$ssbmsy
      if (k >= 1 & k <= 4) {
        MSY_re[[k]][[j]]$fratio_OM <- GT_OM_perf_wdat[[j]]$OM$F_int[26:94]/MSY$fmsy
        MSY_re[[k]][[j]]$bratio_OM <- GT_OM_perf_wdat[[j]]$OM$SSB[26:95]/MSY$ssbmsy
      } else if (k >= 5 & k <= 8) {
        MSY_re[[k]][[j]]$fratio_OM <- GT_OM_constant_wdat[[j]]$OM$F_int[26:94]/MSY$fmsy
        MSY_re[[k]][[j]]$bratio_OM <- GT_OM_constant_wdat[[j]]$OM$SSB[26:95]/MSY$ssbmsy
      } else if (k >= 9 & k <= 12) {
        MSY_re[[k]][[j]]$fratio_OM <- GT_OM_linear_wdat[[j]]$OM$F_int[26:94]/MSY$fmsy
        MSY_re[[k]][[j]]$bratio_OM <- GT_OM_linear_wdat[[j]]$OM$SSB[26:95]/MSY$ssbmsy
      } else if(k >= 13 & k <= 16) {
        MSY_re[[k]][[j]]$fratio_OM <- GT_OM_curvilinear_wdat[[j]]$OM$F_int[26:94]/MSY$fmsy
        MSY_re[[k]][[j]]$bratio_OM <- GT_OM_curvilinear_wdat[[j]]$OM$SSB[26:95]/MSY$ssbmsy #would ssb_msy or fmsy be different for different scenarios?
      }
      MSY_re[[k]][[j]]$fratio_re <- (MSY_re[[k]][[j]]$fratio_EM - MSY_re[[k]][[j]]$fratio_OM) / MSY_re[[k]][[j]]$fratio_OM
      MSY_re[[k]][[j]]$bratio_re <- (MSY_re[[k]][[j]]$bratio_EM - MSY_re[[k]][[j]]$bratio_OM) / MSY_re[[k]][[j]]$bratio_OM
      MSY_re[[k]][[j]]$fratio_e <- (MSY_re[[k]][[j]]$fratio_EM - MSY_re[[k]][[j]]$fratio_OM)
      MSY_re[[k]][[j]]$bratio_e <- (MSY_re[[k]][[j]]$bratio_EM - MSY_re[[k]][[j]]$bratio_OM)
      
      MSY_re[[k]][[j]]$fmsy_re <- (MSY_re[[k]][[j]]$fmsy-MSY$fmsy)/MSY$fmsy
      MSY_re[[k]][[j]]$bmsy_re <- (MSY_re[[k]][[j]]$ssbmsy-MSY$ssbmsy)/MSY$ssbmsy
      MSY_re[[k]][[j]]$f_re <- (exp(res_list_final[[k]][[j]]$SD$par.fixed[9:77])-GT_OM_perf_wdat[[j]]$OM$F_int[26:94])/GT_OM_perf_wdat[[j]]$OM$F_int[26:94]
      MSY_re[[k]][[j]]$fmsy_e <- (MSY_re[[k]][[j]]$fmsy-MSY$fmsy)
      MSY_re[[k]][[j]]$bmsy_e <- (MSY_re[[k]][[j]]$ssbmsy-MSY$ssbmsy)
      MSY_re[[k]][[j]]$f_e <- (exp(res_list_final[[k]][[j]]$SD$par.fixed[9:77])-GT_OM_perf_wdat[[j]]$OM$F_int[26:94])
      
    } else {
      MSY_re[[k]][[j]] <- NA
    }
  }
}

save(MSY_re,file= "C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/Output/MSY_re.rds")

#load("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/Output/MSY_re.rds")

#Getting it into usable format
fratio_re<-bratio_re<-bratio_e<-fratio_e<-fmsy_re<-bmsy_re<-f_re<-f_e<-fmsy_e<-bmsy_e<-list()
for(i in 1:16){
  fratio_re[[i]]<-fratio_e[[i]]<-matrix(NA,nrow=100,ncol=69)
  bratio_re[[i]]<-bratio_e[[i]]<-matrix(NA,nrow=100,ncol=70)
  f_re[[i]]<-f_e[[i]]<-matrix(NA,nrow=100,ncol=69)
  fmsy_re[[i]]<-bmsy_re[[i]]<-fmsy_e[[i]]<-bmsy_e[[i]]<-NA
  for(j in 1:100){
    if (!is.null(res_list_final[[i]][[j]]$hessian)) {
    fratio_re[[i]][j,] <- MSY_re[[i]][[j]]$fratio_re
    bratio_re[[i]][j,] <- MSY_re[[i]][[j]]$bratio_re
    fratio_e[[i]][j,] <- MSY_re[[i]][[j]]$fratio_e
    bratio_e[[i]][j,] <- MSY_re[[i]][[j]]$bratio_e
    fmsy_re[[i]][j] <- MSY_re[[i]][[j]]$fmsy_re
    bmsy_re[[i]][j] <- MSY_re[[i]][[j]]$bmsy_re
    f_re[[i]][j,] <- MSY_re[[i]][[j]]$f_re
    f_e[[i]][j,] <- MSY_re[[i]][[j]]$f_e
    fmsy_e[[i]][j] <- MSY_re[[i]][[j]]$fmsy_e
    bmsy_e[[i]][j] <- MSY_re[[i]][[j]]$bmsy_e
    }
  }
}

#Time series boxplots of RE in Fratio and Bratio
OM_Title <- c(rep("No AE", 4), rep("AE Constant Bias", 4), rep("AE Linear Bias", 4), rep("AE Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "AE Constant Bias", "AE Linear Bias", "AE Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(4.1, 4.6, 2.1, 1.1))
for (i in 1:nrow(scenarios)) {
  boxplot(fratio_re[[i]], ylim = c(-1.0, 1.0), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
  axis(1, at = 1, labels = 1, cex.axis = 1.1)
  axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1)
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = -2, cex = 1.3, outer = TRUE)
mtext("Relative Error in Fratio", side = 2, line = -2, cex = 1.3, outer = TRUE)


OM_Title <- c(rep("No AE", 4), rep("AE Constant Bias", 4), rep("AE Linear Bias", 4), rep("AE Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "AE Constant Bias", "AE Linear Bias", "AE Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(4.1, 4.6, 2.1, 1.1))
for (i in 1:nrow(scenarios)) {
  boxplot(bratio_e[[i]], ylim = c(-1.0, 1.0), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
  axis(1, at = 1, labels = 1, cex.axis = 1.1)
  axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1)
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = -2, cex = 1.3, outer = TRUE)
mtext("Relative Error in Bratio", side = 2, line = -2, cex = 1.3, outer = TRUE)

######################################################
#Getting error in M and R0 to see what is going on
######################################################
M_R0_re <- vector("list", length(res_list_final))
for (i in 1:length(res_list_final)) {
  #add in if/else statement to skip over nonconverged iteration where res_list_final[[k]][[j]]$SD$par doesn't exist
  M_R0_re[[i]]<-matrix(NA, nrow=100,ncol=3)
  colnames( M_R0_re[[i]])<-c("M_re", "R0_re", "rSD_re")
  for (j in 1:length(res_list_final[[1]])) {
    hessian <- res_list_final[[i]][[j]]$hessian
    # Check if the Hessian is positive definite
    if (!is.null(hessian)) {

      M_R0_re[[i]][j,1]<-(exp(res_list_final[[i]][[j]]$SD$par.fixed["log_M"]) - 0.3015598)/0.3015598
      M_R0_re[[i]][j,2]<-(exp(res_list_final[[i]][[j]]$SD$par.fixed["log_R0"]) - exp(9.7608))/exp(9.7608)
      M_R0_re[[i]][j,3]<-(exp(res_list_final[[i]][[j]]$SD$par.fixed["log_sigma_rec"]) - 0.3582)/0.3582
      
      }}}



OM_Title <- c(rep("No AE", 4), rep("AE Constant Bias", 4), rep("AE Linear Bias", 4), rep("AE Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "AE Constant Bias", "AE Linear Bias", "AE Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(4.1, 4.6, 2.1, 1.1))
for (i in 1:nrow(scenarios)) {
  boxplot(M_R0_re[[i]], ylim = c(-1.0, 1.0), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  axis(1, at = 1:2, labels = c("M", "R0"), cex.axis = 1.1)
  axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1)
  abline(0, 0, lwd = 2)
}
mtext("Relative Error in M and R0", side = 2, line = -2, cex = 1.3, outer = TRUE)


#Fmsy, bmsy,
OM_Title <- c(rep("No AE", 4), rep("AE Constant Bias", 4), rep("AE Linear Bias", 4), rep("AE Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "AE Constant Bias", "AE Linear Bias", "AE Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(4.1, 4.6, 2.1, 1.1))
for (i in 1:nrow(scenarios)) {
  boxplot(cbind(fmsy_re[[i]],bmsy_re[[i]]), ylim = c(-1.0, 1.0), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  axis(1, at = 1:2, labels = c("Fmsy", "SSBmsy"), cex.axis = 1.1)
  axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1)
  abline(0, 0, lwd = 2)
}
mtext("Relative Error", side = 2, line = -2, cex = 1.3, outer = TRUE)

#Fs
OM_Title <- c(rep("No AE", 4), rep("AE Constant Bias", 4), rep("AE Linear Bias", 4), rep("AE Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "AE Constant Bias", "AE Linear Bias", "AE Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(4.1, 4.6, 2.1, 1.1))
for (i in 1:nrow(scenarios)) {
  boxplot(f_re[[i]], ylim = c(-0.3, 0.3), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
  axis(1, at = 1, labels = 1, cex.axis = 1.1)
  axis(2, at = seq(-0.3, 0.3, by=0.1), labels = seq(-0.3, 0.3, by=0.1), cex.axis = 1.1)
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = -2, cex = 1.3, outer = TRUE)
mtext("Relative Error in F", side = 2, line = -2, cex = 1.3, outer = TRUE)



