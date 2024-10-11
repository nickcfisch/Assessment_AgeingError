#############################
#Population model Simulator
#############################
SimPop<-function(seed=1,                         #seed to start random number generation for reproducibility
                 fage=0,                         #first age of population
                 lage=10,                        #last age/plus group of population
                 fyear=1,                        #first year of population
                 lyear=94,                      #last year of population
                 Linf=25,                        #asymptotic size 
                 a3=0.5,                         #SS-like parameterization of growth, a3 parameter
                 L1=10,                          #ditto^ L1 param
                 BK=0.4,                         #brody growth coefficient
                 Weight_scaling=1.7e-5,          #Weight-length a
                 Weight_allometry=2.9,           #Weight-length b
                 Mref=0.3015598,                 #Reference M for constant or lorenzen. 
                 M_pow=1.775641,                 #power for lorenzen M
                 Mat_50=15.9,                    #Mat param1, Age at which 50% maturity occurs
                 Mat_slope=-0.9,                 #Mat param2, Slope for logistic function of maturity
                 Sel_50=15.9,                    #Sel param1, Age at which 50% selectivity occurs
                 Sel_slope=3.3,                  #Sel param2, Slope for logistic function of selectivity
                 B1=4.374691,                    #Double normal selectivity parameters
                 B2=-3,                          #Nicks best approximation of trigger selectivity
                 B3=1.214063,
                 B4=1.582468,
                 R0=exp(16),                     #Unfished recruitment 
                 h=0.59,                         #Steepness
                 sd_rec=0.73,                    #Recruitment SD
                 const_F=FALSE,                  #Set constant fishing mortality??
                 fint=0.25,                      #fully-selected fishing mortality to set each year if fishing mortality is constant
                 fhigh=0.25,                     #F high for fishing mortality ramp if const_F is false
                 flow=0.25,                      #F low for fishing mortality ramp if const_F is false
                 F_man=FALSE,                    #Set F manually
                 F_val=F_val,                    #Vector of F value to manually set F
                 stochastic=TRUE){               #Stochastic recruitment? If false then model is deterministic
  
  set.seed(seed)
  
  #Length at age
  Lmin<-0
  b<-(L1-Lmin)/a3
  Laa<-Lmin+b*0
  Laa<-NA
  Laa[1]<-Lmin+b*0
  Laa[2:(lage+1)]<-Linf+(L1-Linf)*exp(-BK*((1:lage)-a3))
#  Laa[2:lage]<-Linf+(L1-Linf)*exp(-BK*((1:(lage-1))-a3))
#  Laa[lage+1]<-Linf
  Laa[1]<-uniroot(f=function(x) x+(x-Linf)*(exp(-BK)-1)-Laa[2], interval=c(0.01,100))$root
  
  #Weight at age
  Waa<-Weight_scaling*Laa^Weight_allometry
  Fec<-51.357*Laa^2.8538
    
  #Natural mortality at age
  Maa<-Mref*(Laa/(Linf*0.75))^-M_pow
  #Maa<-rep(Mref,length(fage:lage))    #Constant M
  
  #Maturity
  Mat<-0.5*c(0, 0, 0.79, 0.91, 0.98, 0.99, 1, 1, 1, 1, 1) #Mat read in from assessment?
 # Mat<-0.5/(1+exp(Mat_slope*(Laa-Mat_50)))   #Changed mat to max at 0.5 to account for 50% females
  
  #Fishery Selectivity 
#  Sel<-1/(1+exp(-log(19)*(Laa-Sel_50)/Sel_slope))         #Logistic based on mean length at age
#  Sel<-1/(1+exp(-log(19)*((fage:lage)-Sel_50)/Sel_slope))  #Logistic based on age
  
  #Double Normal
  age<-fage:lage
  peak2<-B1+1+((0.99*lage-B1-1)/(1+exp(-B2)))
  j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
  j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
  asc<-exp(-(age-B1)^2/exp(B3))
  dsc<-exp(-(age-peak2)^2/exp(B4))
  Sel<-asc*(1-j1)+j1*((1-j2)+j2*dsc)
  
  #Fishing intensity, starts in year 25
  k_int<-0.15
  mid_int<-20
  F_int<-NA
  if (F_man==TRUE){
    F_int[fyear:25]<-0
    F_int[26:lyear]<-F_val
  } else if (F_man==FALSE){
  if (const_F==TRUE){
    F_int[fyear:lyear]<-fint
  } else if (const_F==FALSE){
    F_int[fyear:25]<-0
    F_int[26:85]<-fhigh/length(26:85)*(1:length(26:85))
    F_int[86:lyear]<-F_int[85]+(flow-fhigh)/length(86:lyear)*(1:length(86:lyear))
  }
  }
  
  #Now Pop Stuff
  lxo<-c(1,cumprod(exp(-Maa[1:lage])))   #survivorship
  lxo[lage+1]<-lxo[lage+1]/(1-exp(-Maa[lage+1]))   #plus group survivorship
  N0aa<-R0*lxo
  SSB0<-sum(N0aa*Mat*Fec)  #Unfished SSB calc
  
  #Population Simulation
  Faa<-Zaa<-matrix(NA, nrow=lyear, ncol=lage+1)
  Naa<-matrix(NA, nrow=lyear+1, ncol=lage+1)
  Naa[1,]<-N0aa
  SSB<-sum(Naa[1,]*Mat*Fec)
  lrecdevs<-0
  for(i in fyear:lyear){
    Faa[i,]<-F_int[i]*Sel
    Zaa[i,]<-Maa+Faa[i,]
    for(j in 1:lage){
      Naa[i+1,j+1]<-Naa[i,j]*exp(-Zaa[i,j])
    }
    Naa[i+1,lage+1]<-Naa[i+1,lage+1]+Naa[i,lage+1]*exp(-Zaa[i,lage+1]) #plus group
    
    SSB[i+1] <- sum(Naa[i+1,]*Mat*Fec, na.rm=TRUE)
    if (stochastic==TRUE){
      lrecdevs[i+1]<-rnorm(n=1,mean=0, sd=sd_rec)
      Naa[i+1,1] <- (4*h*R0*SSB[i+1])/(SSB0*(1-h)+SSB[i+1]*(5*h-1))*exp(lrecdevs[i+1]-0.5*sd_rec^2)
    } else {
      Naa[i+1,1] <- (4*h*R0*SSB[i+1])/(SSB0*(1-h)+SSB[i+1]*(5*h-1))
    }
  }
  
  Caa<-Faa/Zaa*Naa[1:lyear,]*(1-exp(-Zaa))
  
  return(list(fage=fage,lage=lage,seed=seed,fyear=fyear,lyear=lyear,Linf=Linf,a3=a3,L1=L1,BK=BK,Weight_scaling=Weight_scaling,
              Weight_allometry=Weight_allometry,Mref=Mref,Mat_50=Mat_50,Mat_slope=Mat_slope,Sel_50=Sel_50,Sel_slope=Sel_slope,
              B1=B1,B2=B2,B3=B3,B4=B4,R0=R0,h=h,sd_rec=sd_rec,fint=fint,fhigh=fhigh,flow=flow,stochastic=stochastic,
              lrecdevs=lrecdevs,
              Laa=Laa,
              Waa=Waa,
              Fec=Fec,
              Mat=Mat,
              F_int=F_int,
              lxo=lxo,
              Naa=Naa,
              Caa=Caa,
              SSB=SSB,
              SSB0=SSB0,
              N0aa=N0aa,
              Maa=Maa,
              Sel=Sel,
              Faa=Faa,
              Zaa=Zaa))
}


#########################################################################
#Data simulator, function which simulates data from a population model
#########################################################################
Get_Data<-function(OM=NA,              #Operating model from which to model
                   dat_seed=1,         #seed to start random number generation for reproducibility
                   fyear_dat=26,       #first year that has data
                   lyear_dat=94,      #last year with data
                   sd_catch=0.05,      #SD of catch data
                   #Composition sample size in each year (vector of length(fyear:lyear)). If zero then no comp data for that year
                   N_Comp=c(rep(0,25),30,rep(0,9),40,rep(0,9),50,rep(0,4),60,rep(0,4),70,rep(0,4),80,rep(0,4),90,rep(0,4),rep(100,24)),
                   q_index=0.0001,     #Catchability of fishery index 
                   sd_index=0.25,      #standard deviation of fishery index
                   #TRUE Ageing error matrix, if identity (diag), then no ageing error. This needs to have dimension length(fage:lage)*length(fage:lage)
                   AE_mat=diag(11))    #Given your true age i (row), this matrix defines the probability that you will be classified age j (column), #REMEMBER THIS NEEDS TO BE IN THE NUMBER OF AGES IN YOUR MODEL
{
  
  set.seed(dat_seed)
  #Getting Data
  Obs_Catch<-Obs_Index<-NA
  Obs_Catch_Comp<-Obs_Catch_Comp_noAE<-matrix(0,nrow=length(fyear_dat:lyear_dat),ncol=length(OM$fage:OM$lage))
  Obs_Catch_Comp_wiY<-array(0,dim=c(length(fyear_dat:lyear_dat),length(OM$fage:OM$lage),length(OM$fage:OM$lage)))
  for (d in fyear_dat:lyear_dat){
    Obs_Catch[d-(fyear_dat-1)]<-rlnorm(1, meanlog=log(sum(OM$Caa[d,]*OM$Waa)), sdlog=sd_catch)
    Obs_Index[d-(fyear_dat-1)]<-rlnorm(1, meanlog=log(sum(OM$Naa[d,]*((1-exp(-OM$Zaa[d,]))/OM$Zaa[d,])*OM$Sel*OM$Waa)*q_index), sdlog=sd_index)
    Obs_Catch_Comp_noAE[d-(fyear_dat-1),]<-rmultinom(n=1,size=N_Comp[d], prob=OM$Caa[d,])
    # *I think* this line code would be equivalent to the more kludgy and inefficient code below, but we won't run into efficiency issues so commented out in favor of double sampler
    #    Obs_Catch_Comp[d-(fyear_dat-1),]<-rmultinom(n=1,size=N_Comp[d], prob=(OM$Caa%*%AE_mat)[d,])
    
    #Getting observed data with Ageing error
    #Another sampler is needed to get data in integers, as opposed to commented out line below
    for (a in 1:length(OM$fage:OM$lage)){
      if(Obs_Catch_Comp_noAE[d-(fyear_dat-1),a]>0){
        Obs_Catch_Comp_wiY[d-(fyear_dat-1),a,]<-rmultinom(1,size=Obs_Catch_Comp_noAE[d-(fyear_dat-1),a],prob=AE_mat[a,])
      }
    }
    Obs_Catch_Comp[d-(fyear_dat-1),]<-colSums(Obs_Catch_Comp_wiY[d-(fyear_dat-1),,])
  }
  
  #This would be *expected* data with ageing error but in decimals. 
  #  Obs_Catch_Comp<-Obs_Catch_Comp_noAE%*%AE_mat  #Getting observed data with Ageing error
  
  return(list(OM=OM,dat_seed=dat_seed,sd_catch=sd_catch,N_Comp=N_Comp,q_index=q_index,sd_index=sd_index,fyear_dat=fyear_dat,lyear_dat=lyear_dat,
              Obs_Catch=Obs_Catch,
              Obs_Catch_CompnoAE=Obs_Catch_Comp_noAE,
              Obs_Catch_Comp=Obs_Catch_Comp,
              Obs_Index=Obs_Index))
}


####################################################################
#Getting Data from OM and applying Age Error the composition data
####################################################################
OM_Err <- function(OM_text, AE_mat, N_sim){
  OM_wdat<-list()
  N_comp<-c(rep(0,25),30,rep(0,9),40,rep(0,9),50,rep(0,4),60,rep(0,4),70,rep(0,4),80,rep(0,4),90,rep(0,4),rep(100,24))
  sd_catch<-0.05
  sd_index<-0.25
  fyear_dat<-26
  lyear_dat<-94
  for (s in N_sim){
    OM_wdat[[s]]<-Get_Data(OM=Triggerfish_runs[[s]],AE_mat=AE_mat,dat_seed=s,fyear_dat=fyear_dat,lyear_dat=lyear_dat,sd_catch=sd_catch,N_Comp=N_comp,q_index=0.0001,sd_index=sd_index)
  }
  save(OM_wdat,file=paste0(wd,"/",OM_text,".RData"))
}


#################################################
#EM including ageing error matrix 
#################################################
sim_Fn <- function(OM_text, N_sim, AE_mat, max_jitter){
  load(paste0(wd,"/",OM_text,".RData"))

  Triggerfish_OM<-OM_wdat

  #Doing N Simulations
  res_list<-list()
  for (s in N_sim){
    
    OM<-Triggerfish_OM[[s]]
    
    dat<-list(fyear=OM$OM$fyear, lyear=69, fage=OM$OM$fage, lage=OM$OM$lage, 
              years=OM$OM$fyear:69, ages=OM$OM$fage:OM$OM$lage,
              obs_harv=OM$Obs_Catch,
              obs_index=OM$Obs_Index,
              obs_fishery_comp=OM$Obs_Catch_Comp/rowSums(OM$Obs_Catch_Comp),
              SS_fishery=rowSums(OM$Obs_Catch_Comp),
              Mat=OM$OM$Mat,
              Laa=OM$OM$Laa,
              Waa=OM$OM$Waa,
              Fec=OM$OM$Fec,
              Lamda_Harvest=1,                           #Switch for whether to use a data source or not, 0=no, 1=yes
              Lamda_Comp=1,
              Lamda_Index=1,
              AE_mat=AE_mat)                             
    
    #Parameters
    set.seed(s)
    #Starting parameters drawing from uniform 35% below and above true value
    par <- list(log_M=log(runif(1,min=OM$OM$Mref-OM$OM$Mref*0.35,max=OM$OM$Mref+OM$OM$Mref*0.35)),
                log_q=log(runif(1,min=OM$q_index-OM$q_index*0.35,max=OM$q_index+OM$q_index*0.35)),
                log_recruit_devs_init=rep(0,dat$lage),
                log_recruit_devs=rep(0,dat$lyear),
                steepness=OM$OM$h,
                log_R0=log(runif(1,min=OM$OM$R0-OM$OM$R0*0.35,max=OM$OM$R0+OM$OM$R0*0.35)),
                log_sigma_rec=log(OM$OM$sd_rec),
                log_sd_catch=log(OM$sd_catch),
                log_sd_index=log(OM$sd_index),
#                Sel_logis_k=log(runif(1,min=OM$OM$Sel_slope-OM$OM$Sel_slope*0.35,max=OM$OM$Sel_slope+OM$OM$Sel_slope*0.35)),
#                Sel_logis_midpt=log(runif(1,min=OM$OM$Sel_50-OM$OM$Sel_50*0.35,max=OM$OM$Sel_50+OM$OM$Sel_50*0.35)),
                B1=runif(1,min=OM$OM$B1-abs(OM$OM$B1)*0.35,max=OM$OM$B1+abs(OM$OM$B1)*0.35),                       #Double normal selectivity parameters
                B2=runif(1,min=OM$OM$B2-abs(OM$OM$B2)*0.35,max=OM$OM$B2+abs(OM$OM$B2)*0.35),
                B3=runif(1,min=OM$OM$B3-abs(OM$OM$B3)*0.35,max=OM$OM$B3+abs(OM$OM$B3)*0.35),
                B4=runif(1,min=OM$OM$B4-abs(OM$OM$B4)*0.35,max=OM$OM$B4+abs(OM$OM$B4)*0.35),
                log_fint=log(runif(length(OM$OM$F_int[26:94]),min=OM$OM$F_int[26:94]-OM$OM$F_int[26:94]*0.35,max=OM$OM$F_int[26:94]+OM$OM$F_int[26:94]*0.35)))  

    ################
    #TMB stuff
    ################
    dyn.load(dynlib("SCAA_forDerek_wAE"))
    
    parm_names<-names(MakeADFun(dat, par, DLL="SCAA_forDerek_wAE")$par)
    
    fixed<-list(steepness=factor(NA),
                log_sd_catch=factor(NA),
                log_sd_index=factor(NA))
    
    #Bounds, need to be updated if you go back to logistic selectivity
    lower_bounds<-c(-5,-20,rep(-10,dat$lage),rep(-10,dat$lyear), 0, 5, -5,-5,-5, 1, -5, 0,-4,rep(-10,dat$lyear))
    upper_bounds<-c( 2,  1,rep( 10,dat$lage),rep( 10,dat$lyear), 1, 25, 2, 2, 2, 5,  5, 2, 2,rep(  0,dat$lyear))
    
    reffects=c("log_recruit_devs","log_recruit_devs_init")
    l<-lower_bounds[-which(parm_names %in% c(names(fixed),reffects))]
    u<-upper_bounds[-which(parm_names %in% c(names(fixed),reffects))]
    
    SCAA <- MakeADFun(dat, par, DLL="SCAA_forDerek_wAE", map=fixed, random=reffects)
    
counter<-1  
tryCatch({
  SCAA_fit <- TMBhelper::fit_tmb(obj=SCAA, startpar=SCAA$par, lower=l, upper=u, newtonsteps = 1,getsd=TRUE,bias.correct=TRUE,getHessian=TRUE)
}, error=function(e){
  counter<<-0
  SCAA_fit<<-list(NA)
})

convcounter<-1
if (
  if(length(SCAA_fit) == 2) {
    is.na(SCAA_fit$h[1,1]) || !matrixcalc::is.positive.definite(SCAA_fit$h) || SCAA_fit$opt$max_gradient > 0.1
  } else {
    is.na(SCAA_fit$hessian[1,1]) || !matrixcalc::is.positive.definite(SCAA_fit$hessian) || SCAA_fit$max_gradient > 0.1
  }
) {
  while (
    if(length(SCAA_fit) == 2) {
      (is.na(SCAA_fit$h[1,1]) || !matrixcalc::is.positive.definite(SCAA_fit$h) || SCAA_fit$opt$max_gradient > 0.1) && 
        convcounter < max_jitter && counter == 1
    } else {
      (is.na(SCAA_fit$hessian[1,1]) || !matrixcalc::is.positive.definite(SCAA_fit$hessian) || SCAA_fit$max_gradient > 0.1) && 
        convcounter < max_jitter && counter == 1
    }
  ){
    par <- list(log_M=log(runif(1,min=OM$OM$Mref-OM$OM$Mref*0.35,max=OM$OM$Mref+OM$OM$Mref*0.35)),
                log_q=log(runif(1,min=OM$q_index-OM$q_index*0.35,max=OM$q_index+OM$q_index*0.35)),
                log_recruit_devs_init=rep(0,dat$lage),
                log_recruit_devs=rep(0,dat$lyear),
                steepness=OM$OM$h,
                log_R0=log(runif(1,min=OM$OM$R0-OM$OM$R0*0.35,max=OM$OM$R0+OM$OM$R0*0.35)),
                log_sigma_rec=log(OM$OM$sd_rec),
                log_sd_catch=log(OM$sd_catch),
                log_sd_index=log(OM$sd_index),
                #                Sel_logis_k=log(runif(1,min=OM$OM$Sel_slope-OM$OM$Sel_slope*0.35,max=OM$OM$Sel_slope+OM$OM$Sel_slope*0.35)),
                #                Sel_logis_midpt=log(runif(1,min=OM$OM$Sel_50-OM$OM$Sel_50*0.35,max=OM$OM$Sel_50+OM$OM$Sel_50*0.35)),
                B1=runif(1,min=OM$OM$B1-abs(OM$OM$B1)*0.35,max=OM$OM$B1+abs(OM$OM$B1)*0.35),                       #Double normal selectivity parameters
                B2=runif(1,min=OM$OM$B2-abs(OM$OM$B2)*0.35,max=OM$OM$B2+abs(OM$OM$B2)*0.35),
                B3=runif(1,min=OM$OM$B3-abs(OM$OM$B3)*0.35,max=OM$OM$B3+abs(OM$OM$B3)*0.35),
                B4=runif(1,min=OM$OM$B4-abs(OM$OM$B4)*0.35,max=OM$OM$B4+abs(OM$OM$B4)*0.35),
                log_fint=log(runif(length(OM$OM$F_int[26:94]),min=OM$OM$F_int[26:94]-OM$OM$F_int[26:94]*0.35,max=OM$OM$F_int[26:94]+OM$OM$F_int[26:94]*0.35)))  
    
    
    SCAA <- MakeADFun(dat, par, DLL="SCAA_forDerek_wAE", map=fixed, random=reffects)
    tryCatch({
      SCAA_fit <- TMBhelper::fit_tmb(obj=SCAA, startpar=SCAA$par, lower=l, upper=u, newtonsteps = 1,getsd=TRUE,bias.correct=TRUE,getHessian=TRUE)
    }, error=function(e){
      counter<<-0
      SCAA_fit<<-list(NA)
    })
    convcounter<-sum(convcounter, 1)
  }
}
    SCAA_fit$convcounter <- convcounter
    #res_list saves all sorts of output related to the assessment, a couple of examples below. 
    res_list[[s]]<-SCAA_fit
  }
  return(res_list)
}

#x<-matrix(NA, nrow=100, ncol=length(res_list[[1]]$SD$value))
#for(i in 1:100){
# if(!is.null(res_list[[i]]$SD$value)){
#  x[i,]<-(res_list[[i]]$SD$value-Triggerfish_OM[[i]]$OM$SSB[26:95])/Triggerfish_OM[[i]]$OM$SSB[26:95]
# }
#}

#boxplot(x, ylim=c(-1,1))

#Example of extracting Data
#est_par<-as.list(SCAA_fit$SD, what="Est") #EXACT SAME STRUCTURE AS PARAMETER VECTOR
#sd_par<-as.list(SCAA_fit$SD, what="Std")

#est_der<-as.list(SCAA_fit$SD,report=TRUE, what="Est. (bias.correct)")
#sd_der<-as.list(SCAA_fit$SD,report=TRUE,  what="Std. Error")

#SCAA$report(SCAA$env$last.par.best)
