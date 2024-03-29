
#Double-normal selectivity

Amin<-0
Amax<-10
age<-Amin:Amax
B1 <- 4.2623060
B2 <- -1.9183504
B3 <- 0.9908788
B4 <- 0.4789121
B5 <- -15.7304389
B6 <- -13.3039320

peak2<-B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)))
t1<-exp(-(Amin-B1)^2/exp(B3))
t2<-exp(-(Amax-peak2)^2/exp(B4))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
asc<-(1+exp(-B5))^-1+(1-(1+exp(-B5))^-1)*((exp(-(age-B1)^2/exp(B3))-t1)/(1-t1))
dsc<-1+(((1+exp(-B6))^-1)-1)*((exp(-(age-peak2)/exp(B4))-1)/(t2-1))
Selvec<-asc*(1-j1)+j1*((1-j2)+j2*dsc)

plot(age, Selvec, las=1, pch=16, type="b",xlab="Age", ylab="Selectivity")


###########################################
#Now looking to see plausible bounds/prior
###########################################
for(i in 1:1000){
Amin<-0
Amax<-10
age<-Amin:Amax
B1 <- runif(1,0,40)
B2 <- runif(1,-15,3)
B3 <- runif(1,0,10)
B4 <- runif(1,0,10)
B5 <- runif(1,-20,7)
B6 <- runif(1,-15,17)

peak2<-B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)))
t1<-exp(-(Amin-B1)^2/exp(B3))
t2<-exp(-(Amax-peak2)^2/exp(B4))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
asc<-(1+exp(-B5))^-1+(1-(1+exp(-B5))^-1)*((exp(-(age-B1)^2/exp(B3))-t1)/(1-t1))
dsc<-1+(((1+exp(-B6))^-1)-1)*((exp(-(age-peak2)/exp(B4))-1)/(t2-1))
Selvec<-asc*(1-j1)+j1*((1-j2)+j2*dsc)

lines(age, Selvec, las=1, pch=16, xlab="Age", ylab="Selectivity")
}

