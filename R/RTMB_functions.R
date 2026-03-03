
##################################
#RTMB functions for Derek models
##################################

#Standard Von-Bertalanffy
get_VB<-function(dat){
  library(RTMB)
  
  par <- list(lLinf=log(58),lK=log(0.14),tnot=-4.16,lCV=log(0.2))
  
  VB <- function(par){
    getAll(dat,par)
    Linf<-exp(lLinf)
    K<-exp(lK)
    tnot<-tnot
    CV<-exp(lCV)
    
    pred<-Linf*(1-exp(-K*(0:10-tnot)))
    
    jnll<-0
    for(d in 1:length(age)){
     sumi<-0
      for(i in 1:11){
       sumi <- sumi + AE_mat[i,age[d]+1] * dnorm(length_cm[d],pred[i],CV*pred[i])         #Fit assuming Normal dist of length at age
      }
     jnll <- jnll - log(sumi)
    }
    
    pred_derived <- Linf*(1-exp(-K*(0:10-tnot)))
    
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
    h <- optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
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
  
  par <- list(lk=log(0.2), lx0=log(5))
  
  Mat <- function(par){
    getAll(dat,par)
    k<-exp(lk)
    x0<-exp(lx0)
    
    pred<-1/(1+exp(-k*(0:10-x0)))
    
    jnll<-0
    for(d in 1:length(age)){
      if(age[d]>1){
        sumi<-0
        for(i in 1:11){
          sumi <- sumi + AE_mat[i,age[d]+1] * dbinom(maturity[d],size=1,prob=pred[i])         #Fit assuming Normal dist of length at age
        }
        jnll <- jnll - log(sumi)
      }
    }
    
    pred_derived <- 1/(1+exp(-k*(0:10-x0)))
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
    h <- optimHess(opt$par, obj$fn, obj$gr) # Hessian matrix
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

get_Mat_VB<-function(dat){
  library(RTMB)
  age<-dat$age
  
  par <- list(lK=log(0.3),lb=log(3))
  
  Mat <- function(par){
    getAll(dat,par)
    K<-exp(lK) 
    b<-exp(lb) 

    pred <- (1-exp(-K*(0:10-1)))^b
    
    jnll<-0
    for(d in 1:length(age)){
      if(age[d]>1){
       sumi<-0
       for(i in 1:11){
        sumi <- sumi + AE_mat[i,age[d]+1] * sum(dbinom(maturity[d],size=1,prob=pred[i]),na.rm=TRUE)         #Fit assuming Normal dist of length at age
       }
       jnll <- jnll - log(sumi)
      }
    } 
    
    pred_derived <- (1-exp(-K*(seq(0,10)-1)))^b
    
    ADREPORT(pred_derived)
    ADREPORT(K)
    ADREPORT(b)
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

