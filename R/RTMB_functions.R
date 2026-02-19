
##################################
#RTMB functions for Derek models
##################################

#Standard Von-Bertalanffy
get_VB<-function(dat){
  library(RTMB)
  
  par <- list(lLinf=log(58),lK=log(0.14),tnot=-1,lCV=log(0.1))
  
  VB <- function(par){
    getAll(dat,par)
    Linf<-exp(lLinf)
    K<-exp(lK)
    tnot<-tnot
    CV<-exp(lCV)
    
    pred<-Linf*(1-exp(-K*(age-tnot)))
    
    jnll <- -sum(dnorm(length,pred,CV*pred,log=TRUE),na.rm=TRUE)         #Fit assuming Normal dist of length at age
    
    pred_derived <- Linf*(1-exp(-K*((0:10)-tnot)))
    
    ADREPORT(pred_derived)
    ADREPORT(Linf)
    ADREPORT(K)
    ADREPORT(tnot)
    ADREPORT(CV)
    return(jnll)
  }
  
  obj <- MakeADFun(VB, par)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  newtonsteps <- 1
  for(i in seq_len(newtonsteps)) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite results
    opt <- nlminb(new_par, obj$fn, obj$gr,control = list(eval.max = 1e4, iter.max = 1e4))
  }
  Sdrep<-sdreport(obj)
  est_par<-as.list(Sdrep, what="Est") #EXACT SAME STRUCTURE AS PARAMETER VECTOR
  sd_par<-as.list(Sdrep, what="Std")
  
  est_der<-as.list(Sdrep,report=TRUE, what="Est")
  sd_der<-as.list(Sdrep,report=TRUE,  what="Std")
  
  return(list(Est=est_der, SEs=sd_der))
}

get_Mat<-function(dat){
  library(RTMB)
  
  val<-dat$age
  
  par <- list(lk=log(0.2), lx0=log(5))
  
  Mat <- function(par){
    getAll(dat,par)
    k<-exp(lk)
    x0<-exp(lx0)
    
    pred<-1/(1+exp(-k*(val-x0)))
    
    jnll <- -sum(dbinom(maturity,size=1,prob=pred,log=TRUE),na.rm=TRUE)
    #Force through the zero at age 2
    #    jnll <- -sum(dbinom(maturity,size=1,prob=pred,log=TRUE),na.rm=TRUE) + -dnorm(0,mean=1/(1+exp(-k*(2-x0))), sd=0.001,log=TRUE)
    
    pred_derived <- 1/(1+exp(-k*(seq(0,10)-x0)))
    ADREPORT(pred_derived)
    ADREPORT(k)
    ADREPORT(x0)
    return(jnll)
  }
  
  obj <- MakeADFun(Mat, par)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  newtonsteps <- 1
  for(i in seq_len(newtonsteps)) {
    g <- as.numeric(obj$gr(opt$par))
    h <- stats::optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
    new_par <- opt$par - solve(h, g)
    # rewrite results
    opt <- nlminb(new_par, obj$fn, obj$gr,control = list(eval.max = 1e4, iter.max = 1e4))
  }
  
  Sdrep<-sdreport(obj)
  est_par<-as.list(Sdrep, what="Est") #EXACT SAME STRUCTURE AS PARAMETER VECTOR
  sd_par<-as.list(Sdrep, what="Std")
  
  est_der<-as.list(Sdrep,report=TRUE, what="Est")
  sd_der<-as.list(Sdrep,report=TRUE,  what="Std")
  
  return(list(Est=est_der, SEs=sd_der))
  
}


