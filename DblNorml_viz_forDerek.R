
#Double-normal selectivity

Amin<-0
Amax<-10
age<-Amin:Amax
B1 <- 4.374701
B2 <- -4
B3 <- 1.214091
B4 <- 1.582473

peak2<-B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
#Short form
asc<-exp(-(age-B1)^2/exp(B3))
dsc<-exp(-(age-peak2)^2/exp(B4))
Selvec<-asc*(1-j1)+j1*((1-j2)+j2*dsc)

plot(age, Selvec, las=1, pch=16, type="b",xlab="Age", ylab="Selectivity")


###########################################
#Now looking to see plausible bounds/prior
###########################################
for(i in 1:1000){
Amin<-0
Amax<-10
age<-Amin:Amax
B1 <- runif(1,1,5)
B2 <- runif(1,-5,5)
B3 <- runif(1,0,2)
B4 <- runif(1,-4,2)

peak2<-B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
#Short form
asc<-exp(-(age-B1)^2/exp(B3))
dsc<-exp(-(age-peak2)^2/exp(B4))
Selvec<-asc*(1-j1)+j1*((1-j2)+j2*dsc)

lines(age, Selvec, las=1, pch=16, xlab="Age", ylab="Selectivity")
}


#Finding params that approximate trigger sel

dat<-c(0,0.05,0.2,0.55,0.98,1,0.95,0.55,0.25,0.1,0) #approximate trigger sel
plot(age, dat, las=1, pch=16, type="b",xlab="Age", ylab="Selectivity")

opt<-function(theta){
Amin<-0
Amax<-10
age<-Amin:Amax
B1 <- theta[1]
B2 <- theta[2]
B3 <- theta[3]
B4 <- theta[4]

peak2<-B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
#Short form
asc<-exp(-(age-B1)^2/exp(B3))
dsc<-exp(-(age-peak2)^2/exp(B4))
Selvec<-asc*(1-j1)+j1*((1-j2)+j2*dsc)
ss<-sum((Selvec-dat)^2)
return(ss)
}

B1 <- 4.2623060
B2 <- -1.9183504
B3 <- 0.9908788
B4 <- 0.4789121
theta<-c(B1,B2,B3,B4)
optim(par=theta,fn=opt, lower=c(0,-15,0,0), method= "L-BFGS-B")



#########################
#Trying to approximate M
#########################
v<-function(theta){
pred<-exp(theta[1])*(Triggerfish_runs[[1]]$Laa/(Triggerfish_runs[[1]]$Linf*0.75))^theta[2]
sum((pred-M)^2)
}
optim(p=c(log(0.3),-1),fn=v)

points(0:10, exp(-1.198787)*(Triggerfish_runs[[1]]$Laa/(Triggerfish_runs[[1]]$Linf*0.75))^-1.775641, col=3, pch=16)
Mref<-0.3015598
M_pow<--1.775641
