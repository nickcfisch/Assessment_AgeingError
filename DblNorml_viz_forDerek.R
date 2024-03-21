
#Double-normal selectivity

Amin<-0
Amax<-15
age<-Amin:Amax
B1 <- 2.667
B2 <- -15.885
B3 <- 0.4
B4 <- 1.372
B5 <- -4.010
B6 <- 0.375

peak2<-B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)))
t1<-exp(-(Amin-B1)^2/exp(B3))
t2<-exp(-(Amax-peak2)^2/exp(B4))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
asc<-(1+exp(-B5))^-1+(1-(1+exp(-B5))^-1)*((exp(-(age-B1)^2/exp(B3))-t1)/(1-t1))
dsc<-1+(((1+exp(-B6))^-1)-1)*((exp(-(age-peak2)/exp(B4))-1)/(t2-1))
Selvec<-asc*(1-j1)+j1*((1-j2)+j2*dsc)

plot(age, Selvec, las=1, pch=16, type="b",xlab="Age", ylab="Selectivity")

