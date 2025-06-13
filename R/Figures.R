
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

#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/Sel_and_F.tiff"), height=17, width=40, units='cm', compression="lzw", res=500)
par(mfrow=c(1,2), mar=c(2,3,1,1), oma=c(3,2.5,1,1))
plot(fage:lage, Sel, ylim=c(0,1), las=1, type="b", xlab="Age", ylab="Selectivity", pch=16)
mtext(side=1, text="Age", line=3.5, font=2, cex=1.4)
mtext(side=2, text="Selectivity", line=3.5, font=2, cex=1.4)
plot(1:94, c(rep(0,25),F_val_no_shrimp), lty=1, xlab="Year", ylab="Fishing Mortality", las=1, type="l", lwd=2)
mtext(side=1, text="Year", line=3.5, font=2, cex=1.4)
mtext(side=2, text="Fishing Mortality", line=2.5, font=2, cex=1.4)
#dev.off()

#Calculate RE in SSB
SSB_re <- list()
for (i in 1:nrow(scenarios)) {
  SSB_re[[i]] <- relative_error(get(scenarios[i , 2]),res_list_final[[i]]) #function on runs on hessian PD iterations
}


#Time series boxplots of RE in SSB
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/SSB_re.tiff"), height=35, width=50, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(4,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(SSB_re[[i]], ylim = c(-0.75, 0.75), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.7)
  if(i %in% 13:16){
    axis(1, at = seq(10, 70, by=10), labels = seq(10, 70, by=10), cex.axis = 1.75)
    axis(1, at = 1, labels = 1, cex.axis = 1.75)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-0.75, 0.75, by=0.25), labels = seq(-0.75, 0.75, by=0.25), cex.axis = 1.5, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("   Year", side = 1, line = 2, cex = 1.85, outer = TRUE)
mtext("Relative Error in SSB", side = 2, line = 2, cex = 1.75, outer = TRUE)
#dev.off()
     
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
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/F_ratio.tiff"), height=35, width=50, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(4,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(fratio_re[[i]], ylim = c(-1.0, 1.0), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.1)
  if(i %in% 13:16){
  axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
  axis(1, at = 1, labels = 1, cex.axis = 1.1)
  }
  if(i %in% c(1,5,9,13)){
  axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = 1, cex = 1.3, outer = TRUE)
mtext("Relative Error in F-ratio", side = 2, line = 1, cex = 1.3, outer = TRUE)
#dev.off()


#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/B_ratio.tiff"), height=35, width=50, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(4,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(bratio_re[[i]], ylim = c(-0.4, 0.4), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.1)
  if(i %in% 13:16){
    axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
    axis(1, at = 1, labels = 1, cex.axis = 1.1)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-0.4, 0.4, by=0.1), labels = seq(-0.4, 0.4, by=0.1), cex.axis = 1.1, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = 2, cex = 1.3, outer = TRUE)
mtext("Relative Error in B-ratio", side = 2, line = 2, cex = 1.3, outer = TRUE)
#dev.off()

#Now Raw
#Time series boxplots of E in Fratio and Bratio
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/F_ratio_raw.tiff"), height=35, width=50, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(4,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(fratio_e[[i]], ylim = c(-1.0, 1.0), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.1)
  if(i %in% 13:16){
    axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
    axis(1, at = 1, labels = 1, cex.axis = 1.1)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.1, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = 1, cex = 1.3, outer = TRUE)
mtext("Relative Error in F-ratio", side = 2, line = 1, cex = 1.3, outer = TRUE)
#dev.off()


#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/B_ratio_raw.tiff"), height=35, width=50, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(4,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(bratio_e[[i]], ylim = c(-0.75, 0.75), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.1)
  if(i %in% 13:16){
    axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
    axis(1, at = 1, labels = 1, cex.axis = 1.1)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-0.75, 0.75, by=0.25), labels = seq(-0.75, 0.75, by=0.25), cex.axis = 1.1, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = 2, cex = 1.3, outer = TRUE)
mtext("Relative Error in B-ratio", side = 2, line = 2, cex = 1.3, outer = TRUE)
#dev.off()

#Now only final year, raw ratio errors
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/fratio_bratio_raw_fy.tiff"), height=25, width=35, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(2.5,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(cbind(fratio_e[[i]][,69],bratio_e[[i]][,70]), ylim = c(-1.0, 1.0), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.2)
  if(i %in% 13:16){
    axis(1, at = 1:2, labels = c("F-ratio", "B_ratio"), cex.axis = 1.75)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.3, las=1)
  }
  abline(0, 0, lwd = 1, lty=2)
}
mtext("Error in Terminal Year", side = 2, line = 2, cex = 1.5, outer = TRUE)
#dev.off()


#M and R0
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/M_R0_re.tiff"), height=25, width=35, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(2.5,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(M_R0_re[[i]], ylim = c(-1.0, 1.0), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.2)
  if(i %in% 13:16){
    axis(1, at = 1:3, labels = c("M", "R0", "SD(R)"), cex.axis = 1.75)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.3, las=1)
  }
  abline(0, 0, lwd = 1, lty=2)
}
mtext("Relative Error", side = 2, line = 2, cex = 1.3, outer = TRUE)
#dev.off()

#msys
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/Fmsy_Bmsy_re.tiff"), height=25, width=35, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(2.5,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(cbind(fmsy_re[[i]],bmsy_re[[i]]), ylim = c(-1.0, 1.0), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.0)
  if(i %in% 13:16){
    axis(1, at = 1:2, labels = c("Fmsy", "SSBmsy"), cex.axis = 1.75)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-1.0, 1.0, by=0.5), labels = seq(-1.0, 1.0, by=0.5), cex.axis = 1.3, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("Relative Error", side = 2, line = 2, cex = 1.3, outer = TRUE)
#dev.off()

#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/F_re.tiff"), height=35, width=50, units='cm', compression="lzw", res=500)
OM_Title <- c(rep("No AE", 4), rep("Constant Bias", 4), rep("Linear Bias", 4), rep("Curvilinear Bias", 4))
EM_Title <- c(rep(c("No AE", "Constant Bias", "Linear Bias", "Curvilinear Bias"), 4))
par(mfrow = c(4,4), mar = c(2,2, 1, 1), oma=c(4,4,2,2))
for (i in 1:nrow(scenarios)) {
  boxplot(f_re[[i]], ylim = c(-0.3, 0.3), xlim = c(1, 70), axes = FALSE, 
          frame = TRUE, 
          main = paste0("TRUE = ", OM_Title[i], ", Model = ", EM_Title[i]),
          cex.main = 1.1)
  if(i %in% 13:16){
    axis(1, at = seq(5, 70, by=5), labels = seq(5, 70, by=5), cex.axis = 1.1)
    axis(1, at = 1, labels = 1, cex.axis = 1.1)
  }
  if(i %in% c(1,5,9,13)){
    axis(2, at = seq(-0.3, 0.3, by=0.1), labels = round(seq(-0.3, 0.3, by=0.1),1), cex.axis = 1.3, las=1)
  }
  abline(0, 0, lwd = 2)
}
mtext("     Year", side = 1, line = 1, cex = 1.3, outer = TRUE)
mtext("Relative Error in F", side = 2, line = 1, cex = 1.3, outer = TRUE)
#dev.off()


#Bubble plots for ageing error
AE_no<-diag(11)
library(PBSmodelling)
#tiff(filename=("C:/Users/fischn/Documents/GitHub/Assessment_AgeingError/figures/AE_Scenarios.tiff"), height=12, width=40, units='cm', compression="lzw", res=500)
par(mfrow=c(1,4), mar=c(2,2,1,1), oma=c(3,3,1,1))
plotBubbles(AE_no, xlab="True Age", ylab="Coded Age", las=1,xval=0:10, yval=0:10, main="No Ageing Error",size=0.175, cex.main=1.9, cex.axis=1.5)
abline(b=1,a=0, lty=2)
mtext(side=2, line=3, text="Coded Age", cex=1.5)

plotBubbles(t(AE_mat_constant), xlab="True Age", ylab="Coded Age", las=1,xval=0:10, yval=0:10, main="Constant Bias",size=0.175, cex.main=1.9, cex.axis=1.5)
abline(b=1,a=0, lty=2)

plotBubbles(t(AE_mat_linear), xlab="True Age", ylab="Coded Age", las=1,xval=0:10, yval=0:10, main="Linear Bias",size=0.175, cex.main=1.9, cex.axis=1.5)
abline(b=1,a=0, lty=2)
mtext(side=1, line=3, text="True Age", at=-1.5, cex=1.5)
plotBubbles(t(AE_mat_curvilinear), xlab="True Age", ylab="Coded Age", las=1,xval=0:10, yval=0:10, main="Curvilinear Bias",size=0.175, cex.main=1.9, cex.axis=1.5)
abline(b=1,a=0, lty=2)
#dev.off()
