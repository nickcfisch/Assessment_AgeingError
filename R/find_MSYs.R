
#Code to find the msy stuff with assessment estimates

#source("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/R/Functions.R")
source("C:/Users/Derek.Chamberlin/Work/Research/Age_Err_Simulation/Assessment_AgeingError/R/Functions.R")

find_msy<-function(fage=0,
                    lage=10, #(Sedar 43, 2015)  p. 10
                    fyear=1,
                    lyear=200,
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
                    B1=4.374691,                   #Double normal selectivity parameters
                    B2=-3,                         #Nicks best approximation of trigger selectivity
                    B3=1.214063,
                    B4=1.582468,
                    R0=exp(9.7608), #(Sedar 43, 2015) p. 64
                    h=0.4593, #(Sedar 43, 2015) p. 64
                    sd_rec=0.3582){ #(Sedar 43, 2015) p. 64

  
  Fs<-seq(0.001,1,0.001)
  res<-list()
  res_harv<-NA
  for(i in 1:1000){
   res[[i]]<-SimPop(seed=1,                         #seed to start random number generation for reproducibility
                    fage=fage,                         #first age of population
                    lage=lage,                        #last age/plus group of population
                    fyear=fyear,                        #first year of population
                    lyear=lyear,                      #last year of population
                    Linf=Linf,                        #asymptotic size 
                    a3=a3,                         #SS-like parameterization of growth, a3 parameter
                    L1=L1,                          #ditto^ L1 param
                    BK=BK,                         #brody growth coefficient
                    Weight_scaling=Weight_scaling,          #Weight-length a
                    Weight_allometry=Weight_allometry,           #Weight-length b
                    Mref=Mref,                 #Reference M for constant or lorenzen. 
                    M_pow=M_pow,                 #power for lorenzen M
                    Mat_50=Mat_50,                    #Mat param1, Age at which 50% maturity occurs
                    Mat_slope=Mat_slope,                 #Mat param2, Slope for logistic function of maturity
                    B1=B1,                    #Double normal selectivity parameters
                    B2=B2,                          #Nicks best approximation of trigger selectivity
                    B3=B3,
                    B4=B4,
                    R0=R0,                     #Unfished recruitment 
                    h=h,                         #Steepness
                    sd_rec=sd_rec,                    #Recruitment SD
                    const_F=FALSE,                  #Set constant fishing mortality??
                    F_man=TRUE,                    #Set F manually
                    F_val=Fs[i],                    #Vector of F value to manually set F
                    stochastic=FALSE)
   
   res_harv[i]<-rowSums(t(t(res[[i]]$Caa)*res[[i]]$Waa))[200]
  }
  
  #Getting ssbmsy
  ssb_msy<-SimPop(seed=1,                         #seed to start random number generation for reproducibility
         fage=fage,                         #first age of population
         lage=lage,                        #last age/plus group of population
         fyear=fyear,                        #first year of population
         lyear=lyear,                      #last year of population
         Linf=Linf,                        #asymptotic size 
         a3=a3,                         #SS-like parameterization of growth, a3 parameter
         L1=L1,                          #ditto^ L1 param
         BK=BK,                         #brody growth coefficient
         Weight_scaling=Weight_scaling,          #Weight-length a
         Weight_allometry=Weight_allometry,           #Weight-length b
         Mref=Mref,                 #Reference M for constant or lorenzen. 
         M_pow=M_pow,                 #power for lorenzen M
         Mat_50=Mat_50,                    #Mat param1, Age at which 50% maturity occurs
         Mat_slope=Mat_slope,                 #Mat param2, Slope for logistic function of maturity
         B1=B1,                    #Double normal selectivity parameters
         B2=B2,                          #Nicks best approximation of trigger selectivity
         B3=B3,
         B4=B4,
         R0=R0,                     #Unfished recruitment 
         h=h,                         #Steepness
         sd_rec=sd_rec,                    #Recruitment SD
         const_F=FALSE,                  #Set constant fishing mortality??
         F_man=TRUE,                    #Set F manually
         F_val=seq(0.001,1,0.001)[which.max(res_harv)],                    #Vector of F value to manually set F
         stochastic=FALSE)$SSB[200]
  
  
  
  return(list(fmsy=seq(0.001,1,0.001)[which.max(res_harv)], ssbmsy=ssb_msy))
}

#Finding fmsy for the true operating model
find_msy()  #Looks to be 0.349


#Ok lets look at some assessment results
#load("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/workspace.RData")
load("C:/Users/Derek.Chamberlin/Work/Research/Age_Err_Simulation/Assessment_AgeingError/workspace.RData")

k<-1
j<-1

find_msy(Mref=exp(res_list_final[[k]][[j]]$SD$par.fixed["log_M"]),                 
          B1=res_list_final[[k]][[j]]$SD$par.fixed["B1"],                   
          B2=res_list_final[[k]][[j]]$SD$par.fixed["B2"],                         
          B3=res_list_final[[k]][[j]]$SD$par.fixed["B3"],
          B4=res_list_final[[k]][[j]]$SD$par.fixed["B4"],
          R0=exp(res_list_final[[k]][[j]]$SD$par.fixed["log_R0"]),
          h=0.4593,
          sd_rec=res_list_final[[k]][[j]]$SD$par.fixed["log_sigma_rec"])  #sigma not really needed for msy


