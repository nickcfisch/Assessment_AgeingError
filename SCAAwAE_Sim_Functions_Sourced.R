###############################################################################################
#Age structured population model simulator, and sampling model with Ageing error included
###############################################################################################

#Write where you would like your output
#and .cpp file has to be in working directory
wd<-"C:/Users/fischn/Dropbox/"
wd<-"C:/Users/Derek.Chamberlin/Work/Research/Age_Err_Simulation/Assessment_AgeingError"

source(paste0(wd,"/Functions.R"))

F_val_no_shrimp <- c(0.0021, 0.0008, 0.0044, 0.0073, 0.0104, 0.0144, 0.0179, 
                     0.0214, 0.0251, 0.029, 0.0331, 0.0359, 0.0388, 0.042, 
                     0.0453, 0.0493, 0.0514, 0.0533, 0.0546, 0.0566, 0.0589, 
                     0.061, 0.0636, 0.0641, 0.0659, 0.0675, 0.0713, 0.0791, 
                     0.0849, 0.0896, 0.0993, 0.1075, 0.1183, 0.1271, 0.1429, 
                     0.1534, 0.3418, 0.3873, 0.2956, 0.1406, 0.1144, 0.2267, 
                     0.1387, 0.3773, 0.44, 0.6917, 0.5739, 0.524, 0.537, 0.529, 
                     0.6349, 0.3841, 0.348, 0.3341, 0.3455, 0.2727, 0.2915, 
                     0.3712, 0.4446, 0.5625, 0.4175, 0.2842, 0.3183, 0.3144, 
                     0.2173, 0.1708, 0.2802, 0.1516, 0.2865)

#Example for Gray Triggerfish-like life-history
#For Triggerfish, based on stochastic runs, fmsy is 0.268, MSY is 2969, SSBmsy is 10046
N_sim<-100

Triggerfish_runs<-list()
for (s in 1:N_sim){
  Triggerfish_runs[[s]]<-SimPop(seed=s,
                             fage=0,
                             lage=10, #(Sedar 43, 2015)  p. 10
                             fyear=1,
                             #lyear=100,
                             lyear=length(F_val_no_shrimp)+25,
                             Linf=58.97, #(Sedar 43, 2015) p. 64
                             a3=0.5, 
                             L1=28.3, #(Sedar 43, 2015) p. 64
                             BK=0.14, #(Sedar 43, 2015) p. 64
                             Weight_scaling=2.16e-5, #(Sedar 43, 2015) p. 64
                             Weight_allometry=3.007, #(Sedar 43, 2015) p. 64
                             Mref=0.3015598,                 #Reference M for constant or lorenzen. 
                             M_pow=1.775641,                 #power for lorenzen M
                             Mat_50=31.0, #(Sedar 43, 2015) p. 64
                             Mat_slope=-0.065, #(Sedar 43, 2015) p. 64
                             Sel_50=28.9, # logistic selectivity, not used
                             Sel_slope=7, # logistic selectivity, not used
                             B1=4.374691,                   #Double normal selectivity parameters
                             B2=-3,                         #Nicks best approximation of trigger selectivity
                             B3=1.214063,
                             B4=1.582468,
                             R0=exp(9.7608), #(Sedar 43, 2015) p. 64
                             h=0.4593, #(Sedar 43, 2015) p. 64
                             sd_rec=0.3582, #(Sedar 43, 2015) p. 64
                             const_F=FALSE,
                             fint=0.0021, #based on cumulative F withouth shrimp fleet (Sedar 43)
                             fhigh=0.6917, #based on cumulative F withouth shrimp fleet (Sedar 43)
                             flow=0.0008, #based on cumulative F withouth shrimp fleet (Sedar 43)
                             F_man=TRUE,
                             F_val=F_val_no_shrimp,
                             stochastic=TRUE)
}

#save(Triggerfish_runs, file=paste0(wd,"/Triggerfish_Base.RData"))

#Looking at some plots of population
plot(1:95,Triggerfish_runs[[1]]$SSB/Triggerfish_runs[[1]]$SSB0, ylim=c(0,2.75), las=1, xlab="Year", ylab="SSB/SSB0", main="Triggerfish", type="l", col="grey50")
Triggerfish_Depl<-matrix(NA, nrow=N_sim,ncol=95)
Triggerfish_Depl[1,]<-Triggerfish_runs[[1]]$SSB/Triggerfish_runs[[1]]$SSB0
for(s in 2:N_sim){
  Triggerfish_Depl[s,]<-Triggerfish_runs[[s]]$SSB/Triggerfish_runs[[s]]$SSB0
  lines(1:95,Triggerfish_runs[[s]]$SSB/Triggerfish_runs[[s]]$SSB0, col="grey50")
}

#Ok looking at Inner 75% and inner 95% of simulations
plot(1:95,apply(Triggerfish_Depl,2,quantile,probs=0.975), type="l", lty=2, ylim=c(0,2), las=1, xlab="Year", ylab="SSB/SSB0", main="Gray Triggerfish")
lines(1:95,apply(Triggerfish_Depl,2,median),lty=1)
lines(1:95,apply(Triggerfish_Depl,2,quantile,probs=0.025),lty=2)
lines(1:95,apply(Triggerfish_Depl,2,quantile,probs=0.875),lty=3)
lines(1:95,apply(Triggerfish_Depl,2,quantile,probs=0.125),lty=3)

#############################################################
#TMB SCAAs fit to Fishery data
#############################################################

#TMB Section
library(TMB)

setwd(wd)
#Compile and load model 
compile("SCAA_forDerek_wAE.cpp")

#Ageing Error Definitions
#Need to refine and add in no bias but imprecision scenarios
{
  AE_mat<-diag(length(Triggerfish_runs[[s]]$fage:Triggerfish_runs[[s]]$lage))
  AE_mat_constant <- AE_mat
  AE_mat_linear <- AE_mat
  AE_mat_curvilinear <- AE_mat
  ages<-(1:nrow(AE_mat))-1
  
  #sd = 1
  sd = 0.2
  bias = 1
  plot(AE_mat[,3],col="white")
  for (i in 1:nrow(AE_mat)) {
    for(j in 1:nrow(AE_mat)){
      if(j==1){                      #if age=0 then integrate from 0.5 to 0
        AE_mat_constant[i,j]<-pnorm(ages[j]+0.5, mean = ages[i]-bias, sd = sd)
      }else if (j %in% 2:(nrow(AE_mat)-1)){        #integrate from age+0.5 to age-0.5
        AE_mat_constant[i,j]<-pnorm(ages[j]+0.5, mean = ages[i]-bias, sd = sd)-pnorm(ages[j]-0.5, mean = ages[i]-bias, sd = sd)
      }else if (j==nrow(AE_mat)){    # if you are in plus group integrate from age-0.5 to infinity
        AE_mat_constant[i,j]<-1-pnorm(ages[j]-0.5, mean = ages[i]-bias, sd = sd)
      }
    }
    lines(AE_mat_constant[i,])
  }
  
  
  bias = 0.25
  sd_def <- c(0.170711, 0.170711, 0.341422, 0.512133, 0.682844, 0.853555, 
              1.02427, 1.19498, 1.36569, 1.5364, 1.70711) #Based on GT otolith and spine relationship
  plot(AE_mat[,3],col="white")
  for (i in 1:nrow(AE_mat)) {
    for(j in 1:nrow(AE_mat)){
      if(j==1){                      #if age=0 then integrate from 0.5 to 0
        AE_mat_linear[i,j]<-pnorm(ages[j]+0.5, mean = ages[i]-(bias*ages[i]+0), sd = sd_def[i])
      }else if (j %in% 2:(nrow(AE_mat)-1)){        #integrate from age+0.5 to age-0.5
        AE_mat_linear[i,j]<-pnorm(ages[j]+0.5, mean = ages[i]-(bias*ages[i]+0), sd = sd_def[i])-pnorm(ages[j]-0.5, mean = ages[i]-(bias*ages[i]+0), sd = sd_def[i])
      }else if (j==nrow(AE_mat)){    # if you are in plus group integrate from age-0.5 to infinity
        AE_mat_linear[i,j]<-1-pnorm(ages[j]-0.5, mean = ages[i]-(bias*ages[i]+0), sd = sd_def[i])
      }
    }
    lines(AE_mat_linear[i,])
  }
  
  
  
  #From GT Oto, Old-spine comparison
  #Correct way to do it with the cumulative distribution function
  bias1 = -0.0329
  bias2 = 1.1207
  bias3 = 0.3772
  sd_def <- c(0.170711, 0.170711, 0.341422, 0.512133, 0.682844, 0.853555, 
              1.02427, 1.19498, 1.36569, 1.5364, 1.70711) #Based on GT otolith and spine relationship
  plot(AE_mat[,3],col="white")
  for (i in 1:nrow(AE_mat)) {
    for(j in 1:nrow(AE_mat)){
     if(j==1){                      #if age=0 then integrate from 0.5 to 0
      AE_mat_curvilinear[i,j]<-pnorm(ages[j]+0.5, mean = ((bias1*ages[i]^2)+(bias2*ages[i])-bias3), sd = sd_def[i])
     }else if (j %in% 2:(nrow(AE_mat)-1)){        #integrate from age+0.5 to age-0.5
       AE_mat_curvilinear[i,j]<-pnorm(ages[j]+0.5, mean = ((bias1*ages[i]^2)+(bias2*ages[i])-bias3), sd = sd_def[i])-pnorm(ages[j]-0.5, mean = ((bias1*ages[i]^2)+(bias2*ages[i])-bias3), sd = sd_def[i])
     }else if (j==nrow(AE_mat)){    # if you are in plus group integrate from age-0.5 to infinity
       AE_mat_curvilinear[i,j]<-1-pnorm(ages[j]-0.5, mean = ((bias1*ages[i]^2)+(bias2*ages[i])-bias3), sd = sd_def[i])
     }
    }
   lines(AE_mat_curvilinear[i,])
  }
}

#OMs run with different ageing error scenarios
{
  OM_Err(OM_text = "GT_OM_perf_wdat", AE_mat = AE_mat, N_sim = N_sim)
  OM_Err(OM_text = "GT_OM_constant_wdat", AE_mat = AE_mat_constant, N_sim = N_sim)
  OM_Err(OM_text = "GT_OM_linear_wdat", AE_mat = AE_mat_linear, N_sim = N_sim)
  OM_Err(OM_text = "GT_OM_curvilinear_wdat", AE_mat = AE_mat_curvilinear, N_sim = N_sim)
}


scenarios <- read.csv("Simulation Scenarios for model.csv") #data frame with columns Scenario #, OM_test, AE_mat
#scenarios <- scenarios[c(1,2),] #only OM, EM age error match scenarios

library(foreach)
library(doParallel)
library(parallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

parallel::clusterExport(cl = cl, varlist = c('AE_mat', 'AE_mat_constant', 'AE_mat_linear', 
                                             'AE_mat_curvilinear','scenarios', 'N_sim'), envir = .GlobalEnv)

res_list_final <- list()

res_list_final <- foreach(i=1:nrow(scenarios),.packages='TMB') %dopar% {
  sim_Fn(OM_text = as.character(scenarios[i,2]), 
                                N_sim = N_sim, AE_mat = get(scenarios[i,3]))
}

saveRDS(res_list_final, file=paste0(wd,"/SCAAfit_GT_All.RData"))
save.image("workspace.RData")
stopCluster(cl)





#res_list_perf <- sim_Fn(OM_text = "GT_OM_perf_wdat", N_sim = N_sim, AE_mat = AE_mat)
#saveRDS(res_list_perf, file=paste0(wd,"/SCAAfit_GT_perf.RData"))

#res_list_constant <- sim_Fn(OM_text = "GT_OM_constant_wdat", N_sim = N_sim, AE_mat = AE_mat_constant)
#saveRDS(res_list_constant, file=paste0(wd,"/SCAAfit_GT_constant.RData"))

#res_list_linear <- sim_Fn(OM_text = "GT_OM_linear_wdat", N_sim = N_sim, AE_mat = AE_mat_linear)
#saveRDS(res_list_linear, file=paste0(wd,"/SCAAfit_GT_linear.RData"))

#res_list_curvilinear <- sim_Fn(OM_text = "GT_OM_curvilinear_wdat", N_sim = N_sim, AE_mat = AE_mat_curvilinear)
#saveRDS(res_list_curvilinear, file=paste0(wd,"/SCAAfit_GT_curvilinear.RData"))


#plot(26:95,res_list_perf[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
#res_list_perf[[1]]$diagnostics

#plot(26:95,res_list_constant[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
#res_list_constant[[1]]$diagnostics

#plot(26:95,res_list_linear[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
#res_list_linear[[1]]$diagnostics

#plot(26:95,res_list_curvilinear[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
#res_list_curvilinear[[1]]$diagnostics
