
#Some paper figures
fage=0
lage=10 #(Sedar 43, 2015)  p. 10
fyear=1
lyear=200
Linf=58.97 #(Sedar 43, 2015) p. 64
a3=0.5 
L1=28.3
B1=4.374691                   #Double normal selectivity parameters
B2=-3                        #Nicks best approximation of trigger selectivity
B3=1.214063
B4=1.582468

age<-fage:lage
peak2<-B1+1+((0.99*lage-B1-1)/(1+exp(-B2)))
j1<-(1+exp(-20*((age-B1)/(1+abs(age-B1)))))^-1
j2<-(1+exp(-20*((age-peak2)/(1+abs(age-peak2)))))^-1
asc<-exp(-(age-B1)^2/exp(B3))
dsc<-exp(-(age-peak2)^2/exp(B4))
Sel<-asc*(1-j1)+j1*((1-j2)+j2*dsc)

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


par(mfrow=c(1,2), mar=c(2,3,1,1), oma=c(3,2.5,1,1))
plot(fage:lage, Sel, ylim=c(0,1), las=1, type="b", xlab="Age", ylab="Selectivity", pch=16)
mtext(side=1, text="Age", line=3.5, font=2, cex=1.4)
mtext(side=2, text="Selectivity", line=3.5, font=2, cex=1.4)
plot(1:94, c(rep(0,25),F_val_no_shrimp), lty=1, xlab="Year", ylab="Fishing Mortality", las=1, type="l", lwd=2)
mtext(side=1, text="Year", line=3.5, font=2, cex=1.4)
mtext(side=2, text="Fishing Mortality", line=2.5, font=2, cex=1.4)

     
#load("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/Output/MSY_re.rds")

#Getting it into usable format
fratio_re<-bratio_re<-bratio_e<-fratio_e<-list()
for(i in 1:16){
  fratio_re[[i]]<-fratio_e[[i]]<-matrix(NA,nrow=100,ncol=69)
  bratio_re[[i]]<-bratio_e[[i]]<-matrix(NA,nrow=100,ncol=70)
  for(j in 1:100){
    if (!is.null(res_list_final[[i]][[j]]$hessian)) {
      fratio_re[[i]][j,] <- MSY_re[[i]][[j]]$fratio_re
      bratio_re[[i]][j,] <- MSY_re[[i]][[j]]$bratio_re
      fratio_e[[i]][j,] <- MSY_re[[i]][[j]]$fratio_e
      bratio_e[[i]][j,] <- MSY_re[[i]][[j]]$bratio_e
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

