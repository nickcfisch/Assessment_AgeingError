relative_error <- function(OM, EM) {
  relative_error_all <- matrix(data = NA, nrow = length(OM), ncol = 70)
  for (i in 1:length(EM)) {
    if (!is.null(EM[[i]]$hessian)) {
      if (is.positive.definite(EM[[i]]$hessian)) {
        true_values <- OM[[i]]$OM$SSB[26:95]
        estimate <- EM[[i]]$SD$value[1:70]
        relative_error_i <- (estimate - true_values) / true_values
        relative_error_all[i, ] <- relative_error_i
      }
    }
  }
  return(relative_error_all)
}

plot_re <- function(re_df) {
  library(ggplot2)
  library(reshape2)

  median_re <- apply(re_df, 2, median, na.rm = TRUE)
  quantile_re <- apply(re_df, 2, function(x) {
    quantile(x, probs = 0.95, na.rm = TRUE)
  })

  mean_perf_re_df <- data.frame(
    Index = 1:length(median_re),
    MeanRelativeError = median_re,
    LowerCI = median_re - quantile_re,
    UpperCI = median_re + quantile_re
  )

  figure <- ggplot(mean_perf_re_df, aes(x = Index, y = MeanRelativeError)) +
    geom_line(linetype = "solid", size = 1.25) +
    geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red",
      size = 1.25
    ) +
    labs(x = "Year", y = "Median Relative Error") +
    theme_void() + #theme void to remove axes
    scale_x_continuous(breaks = c(1, seq(5, 70, by = 5)), expand = c(0, 0)) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2)
    )

  return(figure)
}

SimPop_msy <- function(seed = 1, #seed to start random number generation for reproducibility
                       fage = 0, #first age of population
                       lage = 10, #last age/plus group of population
                       fyear = 1, #first year of population
                       lyear = 94, #last year of population
                       Linf = 25, #asymptotic size
                       a3 = 0.5, #SS-like parameterization of growth, a3 parameter
                       L1 = 10, #ditto^ L1 param
                       BK = 0.4, #brody growth coefficient
                       CV_LAA = 0.22, #CV of mean length at age
                       Weight_scaling = 1.7e-5, #Weight-length a
                       Weight_allometry = 2.9, #Weight-length b
                       Mref = 0.3015598, #Reference M for constant or lorenzen.
                       M_pow = 1.775641, #power for lorenzen M
                       Mat_50 = 15.9, #Mat param1, Age at which 50% maturity occurs
                       Mat_slope = -0.9, #Mat param2, Slope for logistic function of maturity
                       Sel_50 = 15.9, #Sel param1, Age at which 50% selectivity occurs
                       Sel_slope = 3.3, #Sel param2, Slope for logistic function of selectivity
                       B1 = 4.374691, #Double normal selectivity parameters
                       B2 = -3, #Nicks best approximation of trigger selectivity
                       B3 = 1.214063,
                       B4 = 1.582468,
                       R0 = exp(16), #Unfished recruitment
                       h = 0.59, #Steepness
                       sd_rec = 0.73, #Recruitment SD
                       const_F = FALSE, #Set constant fishing mortality??
                       fint = 0.25, #fully-selected fishing mortality to set each year if fishing mortality is constant
                       fhigh = 0.25, #F high for fishing mortality ramp if const_F is false
                       flow = 0.25, #F low for fishing mortality ramp if const_F is false
                       F_man = FALSE, #Set F manually
                       F_val = F_val, #Vector of F value to manually set F
                       stochastic = TRUE, #Stochastic recruitment? If false then model is deterministic
                       LAA = NA,
                       MAT = NA) {
  set.seed(seed)

  #Length at age
  Laa <- LAA

  #Weight at age
  Waa <- Weight_scaling * Laa^Weight_allometry
  Fec <- 51.357 * Laa^2.8538

  #Natural mortality at age
  Maa <- Mref * (Laa / (Linf * 0.75))^-M_pow
  #Maa<-rep(Mref,length(fage:lage))    #Constant M

  #Maturity
  Mat <- MAT
  # Mat<-0.5/(1+exp(Mat_slope*(Laa-Mat_50)))   #Changed mat to max at 0.5 to account for 50% females

  #Fishery Selectivity
  #  Sel<-1/(1+exp(-log(19)*(Laa-Sel_50)/Sel_slope))         #Logistic based on mean length at age
  #  Sel<-1/(1+exp(-log(19)*((fage:lage)-Sel_50)/Sel_slope))  #Logistic based on age

  #Double Normal
  age <- fage:lage
  peak2 <- B1 + 1 + ((0.99 * lage - B1 - 1) / (1 + exp(-B2)))
  j1 <- (1 + exp(-20 * ((age - B1) / (1 + abs(age - B1)))))^-1
  j2 <- (1 + exp(-20 * ((age - peak2) / (1 + abs(age - peak2)))))^-1
  asc <- exp(-(age - B1)^2 / exp(B3))
  dsc <- exp(-(age - peak2)^2 / exp(B4))
  Sel <- asc * (1 - j1) + j1 * ((1 - j2) + j2 * dsc)

  #Fishing intensity, starts in year 25
  k_int <- 0.15
  mid_int <- 20
  F_int <- NA
  if (F_man == TRUE) {
    F_int[fyear:25] <- 0
    F_int[26:lyear] <- F_val
  } else if (F_man == FALSE) {
    if (const_F == TRUE) {
      F_int[fyear:lyear] <- fint
    } else if (const_F == FALSE) {
      F_int[fyear:25] <- 0
      F_int[26:85] <- fhigh / length(26:85) * (1:length(26:85))
      F_int[86:lyear] <- F_int[85] +
        (flow - fhigh) / length(86:lyear) * (1:length(86:lyear))
    }
  }

  #Now Pop Stuff
  lxo <- c(1, cumprod(exp(-Maa[1:lage]))) #survivorship
  lxo[lage + 1] <- lxo[lage + 1] / (1 - exp(-Maa[lage + 1])) #plus group survivorship
  N0aa <- R0 * lxo
  SSB0 <- sum(N0aa * Mat * Fec) #Unfished SSB calc

  #Population Simulation
  Faa <- Zaa <- matrix(NA, nrow = lyear, ncol = lage + 1)
  Naa <- matrix(NA, nrow = lyear + 1, ncol = lage + 1)
  Naa[1, ] <- N0aa
  SSB <- sum(Naa[1, ] * Mat * Fec)
  lrecdevs <- 0
  for (i in fyear:lyear) {
    Faa[i, ] <- F_int[i] * Sel
    Zaa[i, ] <- Maa + Faa[i, ]
    for (j in 1:lage) {
      Naa[i + 1, j + 1] <- Naa[i, j] * exp(-Zaa[i, j])
    }
    Naa[i + 1, lage + 1] <- Naa[i + 1, lage + 1] + Naa[i, lage + 1] * exp(-Zaa[i, lage + 1]) #plus group

    SSB[i + 1] <- sum(Naa[i + 1, ] * Mat * Fec, na.rm = TRUE)
    if (stochastic == TRUE) {
      lrecdevs[i + 1] <- rnorm(n = 1, mean = 0, sd = sd_rec)
      Naa[i + 1, 1] <- (4 * h * R0 * SSB[i + 1]) / (SSB0 * (1 - h) + SSB[i + 1] * (5 * h - 1)) * exp(lrecdevs[i + 1] - 0.5 * sd_rec^2)
    } else {
      Naa[i + 1, 1] <- (4 * h * R0 * SSB[i + 1]) / (SSB0 * (1 - h) + SSB[i + 1] * (5 * h - 1))
    }
  }

  Caa <- Faa / Zaa * Naa[1:lyear, ] * (1 - exp(-Zaa))

  return(list(fage = fage,
              lage = lage,
              seed = seed,
              fyear = fyear,
              lyear = lyear,
              Linf = Linf,
              a3 = a3,
              L1 = L1,
              BK = BK,
              CV_LAA = CV_LAA,
              Weight_scaling = Weight_scaling,
              Weight_allometry = Weight_allometry,
              Mref = Mref,
              Mat_50 = Mat_50,
              Mat_slope = Mat_slope,
              Sel_50 = Sel_50,
              Sel_slope = Sel_slope,
              B1 = B1,
              B2 = B2,
              B3 = B3,
              B4 = B4,
              R0 = R0,
              h = h,
              sd_rec = sd_rec,
              fint = fint,
              fhigh = fhigh,
              flow = flow,
              stochastic = stochastic,
              lrecdevs = lrecdevs,
              Laa = Laa,
              Waa = Waa,
              Fec = Fec,
              Mat = Mat,
              F_int = F_int,
              lxo = lxo,
              Naa = Naa,
              Caa = Caa,
              SSB = SSB,
              SSB0 = SSB0,
              N0aa = N0aa,
              Maa = Maa,
              Sel = Sel,
              Faa = Faa,
              Zaa = Zaa))
}


find_msy <- function(fage = 0,
                     lage = 10, #(Sedar 43, 2015)  p. 10
                     fyear = 1,
                     lyear = 200,
                     Linf = 58.97, #(Sedar 43, 2015) p. 64
                     a3 = 0.5,
                     L1 = 28.3, #(Sedar 43, 2015) p. 64
                     BK = 0.14, #(Sedar 43, 2015) p. 64
                     Weight_scaling = 2.16e-5, #(Sedar 43, 2015) p. 64
                     Weight_allometry = 3.007, #(Sedar 43, 2015) p. 64
                     Mref = 0.3015598, #Reference M for constant or lorenzen.
                     M_pow = 1.775641, #power for lorenzen M
                     Mat_50 = 31.0, #(Sedar 43, 2015) p. 64
                     Mat_slope = -0.065, #(Sedar 43, 2015) p. 64
                     B1 = 4.374691, #Double normal selectivity parameters
                     B2 = -3, #Nicks best approximation of trigger selectivity
                     B3 = 1.214063,
                     B4 = 1.582468,
                     R0 = exp(9.7608), #(Sedar 43, 2015) p. 64
                     h = 0.4593, #(Sedar 43, 2015) p. 64
                     sd_rec = 0.3582, #(Sedar 43, 2015) p. 64
                     LAA = NA,
                     MAT = NA) {
  
  Fs <- seq(0.001, 1, 0.001)
  res <- list()
  res_harv <- NA
   for (i in 1:1000) {
     res[[i]] <- SimPop_msy(seed = 1, #seed to start random number generation for reproducibility
                            fage = fage, #first age of population
                            lage = lage, #last age/plus group of population
                            fyear = fyear, #first year of population
                            lyear = lyear, #last year of population
                            Linf = Linf, #asymptotic size
                            a3 = a3, #SS-like parameterization of growth, a3 parameter
                            L1 = L1, #ditto^ L1 param
                            BK = BK, #brody growth coefficient
                            Weight_scaling = Weight_scaling, #Weight-length a
                            Weight_allometry = Weight_allometry, #Weight-length b
                            Mref = Mref, #Reference M for constant or lorenzen.
                            M_pow = M_pow, #power for lorenzen M
                            Mat_50 = Mat_50, #Mat param1, Age at which 50% maturity occurs
                            Mat_slope = Mat_slope, #Mat param2, Slope for logistic function of maturity
                            B1 = B1, #Double normal selectivity parameters
                            B2 = B2, #Nicks best approximation of trigger selectivity
                            B3 = B3,
                            B4 = B4,
                            R0 = R0, #Unfished recruitment
                            h = h, #Steepness
                            sd_rec = sd_rec, #Recruitment SD
                            const_F = FALSE, #Set constant fishing mortality??
                            F_man = TRUE, #Set F manually
                            F_val = Fs[i], #Vector of F value to manually set F
                            stochastic = FALSE,
                            LAA = LAA,
                            MAT = MAT)

    res_harv[i] <- rowSums(t(t(res[[i]]$Caa) * res[[i]]$Waa))[200]
  }

  #Getting ssbmsy
  ssb_msy <- SimPop_msy(seed = 1, #seed to start random number generation for reproducibility
                        fage = fage, #first age of population
                        lage = lage, #last age/plus group of population
                        fyear = fyear, #first year of population
                        lyear = lyear, #last year of population
                        Linf = Linf, #asymptotic size
                        a3 = a3, #SS-like parameterization of growth, a3 parameter
                        L1 = L1, #ditto^ L1 param
                        BK = BK, #brody growth coefficient
                        Weight_scaling = Weight_scaling, #Weight-length a
                        Weight_allometry = Weight_allometry, #Weight-length b
                        Mref = Mref, #Reference M for constant or lorenzen.
                        M_pow = M_pow, #power for lorenzen M
                        Mat_50 = Mat_50, #Mat param1, Age at which 50% maturity occurs
                        Mat_slope = Mat_slope, #Mat param2, Slope for logistic function of maturity
                        B1 = B1, #Double normal selectivity parameters
                        B2 = B2, #Nicks best approximation of trigger selectivity
                        B3 = B3,
                        B4 = B4,
                        R0 = R0, #Unfished recruitment
                        h = h, #Steepness
                        sd_rec = sd_rec, #Recruitment SD
                        const_F = FALSE, #Set constant fishing mortality??
                        F_man = TRUE, #Set F manually
                        F_val = seq(0.001, 1, 0.001)[which.max(res_harv)], #Vector of F value to manually set F
                        stochastic = FALSE,
                        LAA = LAA,
                        MAT = MAT)$SSB[200]

  return(list(fmsy = seq(0.001, 1, 0.001)[which.max(res_harv)],ssbmsy = ssb_msy))
}
