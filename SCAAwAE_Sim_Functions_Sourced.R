###############################################################################################
#Age structured population model simulator, and sampling model with Ageing error included
###############################################################################################

#Write where you would like your output
#and .cpp file has to be in working directory
wd<-"C:/Users/fischn/Dropbox/"
wd<-"C:/Users/Derek.Chamberlin/Work/Research/GT_Simulation/Assessment_AgeingError"

source(paste0(wd,"/Functions.R"))

#Example for Gray Triggerfish-like life-history
#For Triggerfish, based on stochastic runs, fmsy is 0.268, MSY is 2969, SSBmsy is 10046
Nsim<-100
Trigger_runs<-list()
for (s in 1:Nsim){
  Trigger_runs[[s]]<-SimPop(seed=s,
                             fage=0,
                             lage=10, #(Sedar 43, 2015)  p. 10
                             fyear=1,
                             lyear=100,
                             Linf=58.97, #(Sedar 43, 2015) p. 64
                             a3=0.5, 
                             L1=28.3, #(Sedar 43, 2015) p. 64
                             BK=0.14, #(Sedar 43, 2015) p. 64
                             Weight_scaling=2.16e-5, #(Sedar 43, 2015) p. 64
                             Weight_allometry=3.007, #(Sedar 43, 2015) p. 64
                             Mref=0.28, #(Sedar 43, 2015) p. 32
                             Mat_50=31.0, #(Sedar 43, 2015) p. 64
                             Mat_slope=-0.065, #(Sedar 43, 2015) p. 64
                             Sel_50=28.9, # logistic selectivity, not used
                             Sel_slope=7, # logistic selectivity, not used
                             B1=4.2623060,                   #Double normal selectivity parameters
                             B2=-1.9183504,                  #Nicks best approximation of trigger selectivity
                             B3=0.9908788,
                             B4=0.4789121,
                             B5=-15.7304389,
                             B6=-13.3039320,
                             R0=exp(9.7608), #(Sedar 43, 2015) p. 64
                             h=0.4593, #(Sedar 43, 2015) p. 64
                             sd_rec=0.3582, #(Sedar 43, 2015) p. 64
                             const_F=FALSE,
                             fint=0.2, #need to update!!!!
                             fhigh=0.5425,#need to update!!!!
                             flow=0.1259,  #need to update!!!!
                             stochastic=TRUE)
}

#save(Trigger_runs, file=paste0(wd,"/Trigger_Base.RData"))

#Looking at some plots of population
plot(1:101,Trigger_runs[[1]]$SSB/Trigger_runs[[1]]$SSB0, ylim=c(0,1.25), las=1, xlab="Year", ylab="SSB/SSB0", main="Flatfish", type="l", col="grey50")
Trigger_Depl<-matrix(NA, nrow=Nsim,ncol=101)
Trigger_Depl[1,]<-Trigger_runs[[1]]$SSB/Trigger_runs[[1]]$SSB0
for(s in 2:Nsim){
  Trigger_Depl[s,]<-Trigger_runs[[s]]$SSB/Trigger_runs[[s]]$SSB0
  lines(1:101,Trigger_runs[[s]]$SSB/Trigger_runs[[s]]$SSB0, col="grey50")
}

#Ok looking at Inner 75% and inner 95% of simulations
plot(1:101,apply(Trigger_Depl,2,quantile,probs=0.975), type="l", lty=2, ylim=c(0,2), las=1, xlab="Year", ylab="SSB/SSB0", main="Flatfish")
lines(1:101,apply(Trigger_Depl,2,median),lty=1)
lines(1:101,apply(Trigger_Depl,2,quantile,probs=0.025),lty=2)
lines(1:101,apply(Trigger_Depl,2,quantile,probs=0.875),lty=3)
lines(1:101,apply(Trigger_Depl,2,quantile,probs=0.125),lty=3)
legend("top", c("Median","Inner 75%","Inner 95%"), lwd=1, lty=c(1,3,2))


#############################################################
#TMB SCAAs fit to Fishery data
#############################################################

#TMB Section
library(TMB)

setwd(wd)
#Compile and load model 
compile("SCAA_forDerek_wAE.cpp")

#Ageing Error Definitions
{
  AE_mat<-diag(length(Trigger_runs[[s]]$fage:Trigger_runs[[s]]$lage))
  AE_mat_constant <- AE_mat
  AE_mat_linear <- AE_mat
  AE_mat_curvilinear <- AE_mat
  
  
  sd = 1
  bias = 1
  plot(AE_mat[,3],col="white")
  for (i in 1:nrow(AE_mat)) {
    AE_mat_constant[,i]<-dnorm(1:nrow(AE_mat), mean = i+bias, sd = sd)/sum(dnorm(1:nrow(AE_mat), mean = i+bias, sd = sd))
    lines(AE_mat_constant[,i])
  }
  
  
  sd = 0.1
  bias = 0.25
  plot(AE_mat[,3],col="white")
  for (i in 1:nrow(AE_mat)) {
    AE_mat_linear[,i]<-dnorm(1:nrow(AE_mat), mean = i+(bias*i+0), sd = sd*i+0)/sum(dnorm(1:nrow(AE_mat), mean = i+(bias*i+0), sd = sd*i+0))
    lines(AE_mat_linear[,i])
  }
  
  
  sd = 0.1
  #From GT Oto, Old-spine comparison
  bias1 = -0.0329
  bias2 = 1.1207
  bias3 = 0.3772
  plot(AE_mat[,3],col="white")
  for (i in 1:nrow(AE_mat)) {
    AE_mat_curvilinear[,i]<-dnorm(1:nrow(AE_mat), mean = i+((bias1*i^2)+(bias2*i)-bias3), sd = sd*i+0)/sum(dnorm(1:nrow(AE_mat), mean = i+((bias1*i^2)+(bias2*i)-bias3), sd = sd*i+0))
    lines(AE_mat_curvilinear[,i])
  }
}

#OMs run with different ageing error scenarios
{
  OM_Err(OM_text = "GT_OM_perf_wdat", AE_mat = AE_mat)
  OM_Err(OM_text = "GT_OM_constant_wdat", AE_mat = AE_mat_constant)
  OM_Err(OM_text = "GT_OM_linear_wdat", AE_mat = AE_mat_linear)
  OM_Err(OM_text = "GT_OM_curvilinear_wdat", AE_mat = AE_mat_curvilinear)
}



N_sim <-1

#build a for loop to run based on a data frame of scenarios


res_list_perf <- sim_Fn(OM_text = "GT_OM_perf_wdat", N_sim = N_sim, AE_mat = AE_mat)
#saveRDS(res_list_perf, file=paste0(wd,"/SCAAfit_GT_perf.RData"))

res_list_constant <- sim_Fn(OM_text = "GT_OM_constant_wdat", N_sim = N_sim, AE_mat = AE_mat_constant)
#saveRDS(res_list_constant, file=paste0(wd,"/SCAAfit_GT_constant.RData"))

res_list_linear <- sim_Fn(OM_text = "GT_OM_linear_wdat", N_sim = N_sim, AE_mat = AE_mat_linear)
#saveRDS(res_list_linear, file=paste0(wd,"/SCAAfit_GT_linear.RData"))

res_list_curvilinear <- sim_Fn(OM_text = "GT_OM_curvilinear_wdat", N_sim = N_sim, AE_mat = AE_mat_curvilinear)
#saveRDS(res_list_curvilinear, file=paste0(wd,"/SCAAfit_GT_curvilinear.RData"))


plot(26:101,res_list_perf[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
res_list_perf[[1]]$diagnostics

plot(26:101,res_list_constant[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
res_list_constant[[1]]$diagnostics

plot(26:101,res_list_linear[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
res_list_linear[[1]]$diagnostics

plot(26:101,res_list_curvilinear[[1]]$SD$value, ylab="SSB", las=1, xlab="Year", type="b", pch=16, ylim=c(0,1e5))
res_list_curvilinear[[1]]$diagnostics
