relative_error <- function(OM, EM){
  relative_error_all <- matrix(data = NA, nrow = length(OM), ncol = 70)
  for (i in 1:length(OM)) {
    if (!is.null(EM[[i]]$hessian)) {
      if(is.positive.definite(EM[[i]]$hessian)){
        true_values <- OM[[i]]$OM$SSB[26:95]
        estimate <- EM[[i]]$SD$value[1:70]
        relative_error_i <- (estimate - true_values) / true_values
        relative_error_all[i,] <- relative_error_i
      }}}
  return(relative_error_all)
}

plot_re <- function(re_df){
  library(ggplot2)
  library(reshape2)
  
  median_re <- apply(re_df, 2, median, na.rm = TRUE)
  quantile_re <- apply(re_df, 2, function(x) quantile(x, probs = 0.95, na.rm = TRUE))
  
  mean_perf_re_df <- data.frame(
    Index = 1:length(median_re),
    MeanRelativeError = median_re,
    LowerCI = median_re - quantile_re, 
    UpperCI = median_re + quantile_re
  )
  
  figure <- ggplot(mean_perf_re_df, aes(x = Index, y = MeanRelativeError)) +
    geom_line(linetype = "solid", size = 1.25) +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.25) +
    labs(x = "Year",
         y = "Median Relative Error") +
    theme_void() + #theme void to remove axes
    scale_x_continuous(breaks = c(1, seq(5, 70, by = 5)), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
  
  return(figure)
}

grid_plot <- function(re_list){
  library(ggplot2)
  
  No_AE_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "No AE"), size = 8, color = "black") +
    theme_void() +  # Remove background and axes
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
  
  Constant_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Constant Bias"), size = 8, color = "black") +
    theme_void() +  # Remove background and axes
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
  
  Linear_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Linear Bias"), size = 8, color = "black") +
    theme_void() +  # Remove background and axes
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
  
  Curvilinear_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Curvilinear Bias"), size = 8, color = "black") +
    theme_void() +  # Remove background and axes
    theme(plot.margin = margin(0, 0, 0, 0),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 2))
  
  OM_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "OM"), size = 8, color = "black") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  EM_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "EM"), size = 8, color = "black") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  y_axis_label_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Year"), size = 8, color = "black") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  y_axis_plot <- ggplot() +
    geom_blank() +
    scale_x_continuous(limits = c(1, 70), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
    geom_text(aes(x = 10, y = 0.5, label = "10"), size = 8, color = "black")+
    geom_text(aes(x = 35, y = 0.5, label = "35"), size = 8, color = "black")+
    geom_text(aes(x = 60, y = 0.5, label = "60"), size = 8, color = "black")+
    geom_text(aes(x = 10, y = 1, label = "|"), size = 8, color = "black")+
    geom_text(aes(x = 35, y = 1, label = "|"), size = 8, color = "black")+
    geom_text(aes(x = 60, y = 1, label = "|"), size = 8, color = "black")+
    theme_void()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  x_axis_label_plot <- ggplot() +
    geom_text(aes(x = 1, y = 1, label = "Relative Error %"), size = 8, color = "black", angle = 270) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  x_axis_plot <- ggplot() +
    geom_blank() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + 
    coord_cartesian(ylim = c(-0.5, 0.5))+
    geom_text(aes(x = 0.5, y = 0.25, label = " 0.25"), size = 5, color = "black")+
    geom_text(aes(x = 0.5, y = -0.25, label = "-0.25"), size = 5, color = "black")+
    geom_text(aes(x = 0, y = 0.25, label = "-"), size = 8, color = "black", hjust = 0, vjust = 0.34)+
    geom_text(aes(x = 0, y = 0, label = "-"), size = 8, color = "black", hjust = 0, vjust = 0.34)+
    geom_text(aes(x = 0, y = -0.25, label = "-"), size = 8, color = "black", hjust = 0, vjust = 0.34)+
    theme_void()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  x_axis_plot2 <- ggplot() +
    geom_blank() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + 
    coord_cartesian(ylim = c(-0.5, 0.5))+
    geom_text(aes(x = 0, y = 0.25, label = "-"), size = 8, color = "black", hjust = 0, vjust = 0.34)+
    geom_text(aes(x = 0, y = 0, label = "-"), size = 8, color = "black", hjust = 0, vjust = 0.34)+
    geom_text(aes(x = 0, y = -0.25, label = "-"), size = 8, color = "black", hjust = 0, vjust = 0.34)+
    theme_void()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  scenario_1_plot <- plot_re(re_list[[1]])
  scenario_2_plot <- plot_re(re_list[[2]])
  scenario_3_plot <- plot_re(re_list[[3]])
  scenario_4_plot <- plot_re(re_list[[4]])
  scenario_5_plot <- plot_re(re_list[[5]])
  scenario_6_plot <- plot_re(re_list[[6]])
  scenario_7_plot <- plot_re(re_list[[7]])
  scenario_8_plot <- plot_re(re_list[[8]])
  scenario_9_plot <- plot_re(re_list[[9]])
  scenario_10_plot <- plot_re(re_list[[10]])
  scenario_11_plot <- plot_re(re_list[[11]])
  scenario_12_plot <- plot_re(re_list[[12]])
  scenario_13_plot <- plot_re(re_list[[13]])
  scenario_14_plot <- plot_re(re_list[[14]])
  scenario_15_plot <- plot_re(re_list[[15]])
  scenario_16_plot <- plot_re(re_list[[16]])
  
  blank_plot <- ggplot() + 
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) 
  
  library(gridExtra)
  
  final_figure <- grid.arrange(arrangeGrob(OM_plot, No_AE_plot, Constant_plot, Linear_plot, 
                                           Curvilinear_plot, blank_plot, blank_plot,
                                           ncol = 1, heights = c(1/19, 4/19, 4/19, 4/19, 4/19, 1/19, 1/19)), 
                               arrangeGrob(EM_plot, No_AE_plot, Constant_plot, Linear_plot, 
                                           Curvilinear_plot, No_AE_plot, Constant_plot, Linear_plot, 
                                           Curvilinear_plot, No_AE_plot, Constant_plot, Linear_plot, 
                                           Curvilinear_plot, No_AE_plot, Constant_plot, Linear_plot, 
                                           Curvilinear_plot, blank_plot, blank_plot,
                                           ncol = 1), 
                               arrangeGrob(blank_plot, scenario_1_plot, scenario_2_plot, 
                                           scenario_3_plot, scenario_4_plot, scenario_5_plot, 
                                           scenario_6_plot, scenario_7_plot, scenario_8_plot, 
                                           scenario_9_plot, scenario_10_plot, scenario_11_plot, 
                                           scenario_12_plot, scenario_13_plot, scenario_14_plot, 
                                           scenario_15_plot, scenario_16_plot, y_axis_plot, 
                                           y_axis_label_plot,
                                           ncol = 1),
                               arrangeGrob(blank_plot, x_axis_plot, x_axis_plot2, 
                                           x_axis_plot2, x_axis_plot2, x_axis_plot2, 
                                           x_axis_plot2, x_axis_plot2, x_axis_plot2, 
                                           x_axis_plot2, x_axis_plot2, x_axis_plot2, 
                                           x_axis_plot2, x_axis_plot2, x_axis_plot2, 
                                           x_axis_plot2, x_axis_plot2, blank_plot,
                                           blank_plot,
                                           ncol = 1),
                               arrangeGrob(x_axis_label_plot,
                                           ncol = 1),
                               ncol = 5, widths = c(1, 1, 1.4, 0.25, 0.15))
  
  
  return(final_figure)
}


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