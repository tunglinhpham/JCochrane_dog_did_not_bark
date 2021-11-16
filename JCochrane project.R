#PREPARING PACKAGES & LIBRARIES

  #Check packages and install if needed
    check_packages = c("readxl", "xtable", "data.table", "dyn", "kableExtra", "mvtnorm", "ggplot2", "rmarkdown")
    package.check = lapply(
      check_packages, FUN = function(x)
      {
        if (!require(x, character.only = TRUE)) {
          install.packages(x, dependencies = TRUE)
          library(x, character.only = TRUE)
        }
      }
    )
  
  #Load library
    library(readxl)
    library(xtable)
    library(data.table)
    library(dyn)
    library(kableExtra)
    library(mvtnorm)
    library(ggplot2)
    library(rmarkdown)
    
set.seed(234)


#GET RELEVANT DATA

  #Edit data source here
  data_link = "D:/University of Kent/01. Lectures/SPR - EC843 - Financial Econometrics_done/R/GWAnnual.csv"
  
  GWAnnual = read.csv(data_link)
  data = as.data.table(GWAnnual)
  
  GWAnnual = read.csv(data_link)
  data = as.data.table(GWAnnual)
  
  data = data[, DP := ((1+CRSP_SPvw)/(1+CRSP_SPvwx)-1)] # This means dividends are reinvested
  data = data[, Ret := (1+CRSP_SPvw)]
  data = data[, RF_Ret := (1+Rfree)]
  data = data[,list(Date, Ret, DP, RF_Ret, D12)]
  full_data = ts(data, start = 1926, end = 2017) #data is the whole data set
  

#REGRESSIONS
  
  #Estimate Return on Dividend Price Ratio
  OLS1 = dyn$lm(Ret ~ lag(DP, -1), data = full_data)
  OLS1_beta = summary(OLS1)$coef[2, 1]
  OLS1_sd = sd(OLS1[["fitted.values"]])*100
  OLS1_tstat = summary(OLS1)$coef[2, 3]
  OLS1_r2 = summary(OLS1)$r.squared*100
  
  #Estimate Equity Premium on Dividend Price Ratio
  OLS2 = dyn$lm((Ret-RF_Ret) ~ lag(DP, -1), data = full_data)
  OLS2_beta = summary(OLS2)$coef[2, 1]
  OLS2_sd = sd(OLS2[["fitted.values"]])*100
  OLS2_tstat = summary(OLS2)$coef[2, 3]
  OLS2_r2 = summary(OLS2)$r.squared*100
  
  #Estimate Dividend Growth on Dividend Price Ratio
  OLS3 = dyn$lm(D12/lag(D12, -1) ~ lag(DP, -1), data = full_data)
  OLS3_beta = summary(OLS3)$coef[2, 1]
  OLS3_sd = sd(OLS3[["fitted.values"]])*100
  OLS3_tstat = summary(OLS3)$coef[2, 3]
  OLS3_r2 = summary(OLS3)$r.squared*100
  
  #Estimate log return on log DP ratio
  OLS4 = dyn$lm(log(Ret) ~ log(lag(DP, -1)), data = full_data)
  OLS4_beta = summary(OLS4)$coef[2, 1]
  OLS4_sd = sd(OLS4[["fitted.values"]])*100
  OLS4_tstat = summary(OLS4)$coef[2, 3]
  OLS4_r2 = summary(OLS4)$r.squared*100
  OLS4_se = summary(OLS4)$coef[2, 2]
  
  #Estimate log dividend growth on log DP ratio
  OLS5 = dyn$lm(log(D12/lag(D12, -1)) ~ log(lag(DP, -1)), data = full_data)
  OLS5_beta = summary(OLS5)$coef[2, 1]
  OLS5_sd = sd(OLS5[["fitted.values"]])*100
  OLS5_tstat = summary(OLS5)$coef[2, 3]
  OLS5_r2 = summary(OLS5)$r.squared*100
  OLS5_se = summary(OLS5)$coef[2, 2]
  
  #Estimate log DP on its lag
  OLS6 = dyn$lm(log(DP) ~ log(lag(DP, -1)), data = full_data)
  OLS6_beta = summary(OLS6)$coef[2, 1]
  OLS6_sd = sd(OLS6$resid)
  

#TABLE 1

  #Table contents
  Table1 = matrix(
    c(
      OLS1_beta, OLS2_beta, OLS3_beta, OLS4_beta, OLS5_beta,
      OLS1_tstat, OLS2_tstat, OLS3_tstat, OLS4_tstat, OLS5_tstat,
      OLS1_r2, OLS2_r2, OLS3_r2, OLS4_r2, OLS5_r2,
      OLS1_sd, OLS2_sd, OLS3_sd, OLS4_sd, OLS5_sd
    ),
    nrow = 5,
    ncol = 4
  )
  colnames(Table1) = c(
    "$b$",
    "$t$",
    "$R^2(\\%)$",
    "$\\sigma(bx)(\\%)$"
  )
  rownames(Table1) = c(
    "$R_{t+1} = a + b(D_t / P_t ) + \\epsilon_{t+1}$",
    "$R_{t+1} - R^f_t = a + b(D_t / P_t ) + \\epsilon_{t+1}$",
    "$D_{t+1} / D_t = a + b(D_t / P_t ) + \\epsilon_{t+1}$",
    "$r_{t+1} = a + b(d_t - p_t ) + \\epsilon^r_{t+1}$",
    "$\\Delta d_{t+1} = a + b(d_t - p_t ) + \\epsilon^{dp}_{t+1}$"
  )
  
  Table1 %>%    
    kbl(row.names = TRUE,
        caption = "Table 1 - Forecasting regressions",
        escape = FALSE,
        digits = c(2, 2, 1, 1),
    ) %>%
    pack_rows("Regression", 1, 5) %>%
    pack_rows("", 1, 3) %>%
    pack_rows("", 4, 5) %>%
    kable_styling("striped", full_width = FALSE)
  

#TABLE 2
  
  #rho
  rho = exp(mean(log(1/full_data[,"DP"])))/(1+exp(mean(log(1/full_data[,"DP"]))))
  
  #Table contents
  Table2 = matrix(
    c(
      OLS4_beta, OLS5_beta, OLS6_beta,
      summary(OLS4)$coef[2, 2], summary(OLS5)$coef[2, 2], summary(OLS6)$coef[2, 2],
      1-rho*OLS6_beta + OLS5_beta, OLS4_beta-(1-rho*OLS6_beta), (1+OLS5_beta-OLS4_beta)/rho,
      sd(residuals(OLS4))*100, cor(residuals(OLS4),residuals(OLS5))*100, cor(residuals(OLS4),residuals(OLS6))*100,
      cor(residuals(OLS4),residuals(OLS5))*100, sd(residuals(OLS5))*100, cor(residuals(OLS5),residuals(OLS6))*100,
      cor(residuals(OLS4),residuals(OLS6))*100, cor(residuals(OLS5),residuals(OLS6))*100, sd(residuals(OLS6))*100,
      0, rho*OLS6_beta-1, OLS6_beta,
      0, rho*0.99-1, 0.99
    ),
    nrow = 3,
    ncol = 8
  )
  colnames(Table2) = c(
    "$\\hat{b}, \\hat{\\phi}$", "$\\sigma(\\hat{b})$", "$implied$",
    "$r$", "$\\Delta d$", "$dp$",
    "$b, \\phi$","$b, \\phi$"
  )
  rownames(Table2) = c("$r$", "$\\Delta d$", "$dp$")
  
  Table2 %>%    
    kbl(row.names = TRUE,
        caption = "Table 2 - Forecasting regressions and null hypothesis",
        escape = FALSE,
        digits = c(3, 3, 3, 1, 1, 1, 4, 4)) %>%
    add_header_above(c(" " = 1, "Estimates"=3, "$\\epsilon$ s.d. (diagonal) and correlation"=3, "Null 1"=1, "Null 2"=1)) %>%
    kable_styling("striped", full_width = FALSE)
  

#MONTE-CARLO SIMULATION FOR RETURN AND DIVIDEND GROWTH
  
  #Simulation function
  MCSim = function(iteration, phi, rho, data_length, cov_matrix)
  {
    sim_result = matrix(0, nrow = iteration, ncol = 5)
    
    for (i in 1:iteration)
    {
      #generate errors from bi-variate normal
      sim_err = rmvnorm(data_length, c(0, 0), cov_matrix)
      
      #generate dividend price
      sim_dp = rep(0, data_length)
      
      #generate initial observation
      if(phi >= 1) {sim_dp[1] = 0} else {sim_dp[1] = rnorm(1, 0, sqrt(cov_matrix[1, 1]/(1-phi^2)))}
      
      #generate the rest of the observation
      for (j in 2:data_length)
      {
        sim_dp[j] = phi * sim_dp[j-1] + sim_err[j,1]
      }
      sim_dp = ts(sim_dp)
      
      #generate dividend growth
      sim_dg = rep(0, data_length)
      
      #generate initial observation
      sim_dg[1] = rnorm(1, 0, sqrt(cov_matrix[2, 2]/(1-(phi*rho-1)^2)))
      
      #generate the rest of the observation
      for (j in 2:data_length)
      {
        sim_dg[j] = (phi*rho-1) * sim_dp[j-1] + sim_err[j,2]
      }
      sim_dg = ts(sim_dg)
      
      #generate return
      sim_ret = sim_err[, 2]-rho*sim_err[, 1]
      sim_ret = ts(sim_ret)
      
      #estimate the regressions using simulated data
      sim_reg_ret = dyn$lm(sim_ret ~ lag(sim_dp, -1))
      sim_reg_dg = dyn$lm(sim_dg ~ lag(sim_dp, -1))
      sim_reg_dp = dyn$lm(sim_dp ~ lag(sim_dp, -1))
      
      #get results from the regressions:
        #return
        sim_result[i, 1] = summary(sim_reg_ret)$coef[2, 1]
        sim_result[i, 2] = summary(sim_reg_ret)$coef[2, 3]
        
        #dividend growth
        sim_result[i, 3] = summary(sim_reg_dg)$coef[2, 1]
        sim_result[i, 4] = summary(sim_reg_dg)$coef[2, 3]
        
        #dividend price ratio
        sim_result[i, 5] = summary(sim_reg_dp)$coef[2, 1]
    }
    
    sim_result = as.data.frame(sim_result)
    colnames(sim_result) = c("br", "t_br", "bd", "t_bd", "phi")
    return(sim_result)
  }
  
  #Simulation parameters (general)
  iteration = 5000
  data_length = nrow(data)
  
    #Real
    cov_matrix = cov(cbind(OLS6$resid, OLS5$resid))
    
      #Run Simulation 1
      phi = OLS6_beta
      sim_real1 = MCSim(iteration, phi, rho, data_length, cov_matrix)
      
      #Run Simulation 2
      phi = 0.99
      sim_real2 = MCSim(iteration, phi, rho, data_length, cov_matrix)
  
    #Excess
      
      #Log return minus log risk-free return
      OLS7 = dyn$lm(log(Ret)-log(lag(RF_Ret, -1)) ~ log(lag(DP, -1)), data = full_data)
      OLS7_beta = summary(OLS7)$coef[2, 1]
      OLS7_tstat = summary(OLS7)$coef[2, 3]
      OLS7_se = summary(OLS7)$coef[2, 2]
      
      #log dividend growth minus log risk-free return
      OLS8 = dyn$lm(log(D12/lag(D12, -1))-log(lag(RF_Ret, -1)) ~ log(lag(DP, -1)), data = full_data)
      OLS8_beta = summary(OLS8)$coef[2, 1]
      OLS8_tstat = summary(OLS8)$coef[2, 3]
      cov_matrix = cov(cbind(OLS6$resid, OLS8$resid))
      
      phi = OLS6_beta
      
      #Run Simulation
      sim_excess1 = MCSim(iteration, phi, rho, data_length, cov_matrix)


#TABLE 3
      
  Table3 = matrix(
    c(
      sum(sim_real1$br > OLS4_beta)/iteration*100, sum(sim_excess1$br > OLS7_beta)/iteration*100,
      sum(sim_real1$t_br > OLS4_tstat)/iteration*100, sum(sim_excess1$t_br > OLS7_tstat)/iteration*100,
      sum(sim_real1$bd > OLS5_beta)/iteration*100, sum(sim_excess1$bd > OLS8_beta)/iteration*100,
      sum(sim_real1$t_bd > OLS5_tstat)/iteration*100, sum(sim_excess1$t_bd > OLS8_tstat)/iteration*100
    ),
    nrow = 2,
    ncol = 4
  )
  colnames(Table3) = c("$b_r$", "$t_r$", "$b_d$","$t_d$")
  rownames(Table3) = c("Real", "Excess")
  
  Table3 %>%
    kbl(row.names = TRUE,
        caption = "Table 3 - Percent probability values under the $\\phi = 0.939$ null",
        escape = FALSE,
        digits = 1) %>%
    kable_styling("striped", full_width = FALSE)
  

#FIGURE 1A
  
  ggplot(sim_real1, aes(x = br, y = bd)) +
    ggtitle(expression(paste("Figure 1a - Coefficients, ", phi, " = 0.94"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(b[r])) +
    ylab(expression(b[d])) +
    geom_jitter(size = 0.25) +
    geom_vline(xintercept = OLS4_beta, color = "blue") +
    geom_hline(yintercept = OLS5_beta, color = "blue") +
    geom_point(aes(x = 0, y = rho*OLS6_beta-1),shape = 18, size = 5, color = "red") +
    geom_point(aes(x = OLS4_beta, y = OLS5_beta),shape = 16, size = 5, color = "blue") +
    annotate("text", x = 0.3, y = 0.05, label = paste0(round(sum(sim_real1$br > OLS4_beta & sim_real1$bd > OLS5_beta)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = 0.3, y = -0.25, label = paste0(round(sum(sim_real1$br > OLS4_beta & sim_real1$bd <= OLS5_beta)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -0.1, y = 0.05, label = paste0(round(sum(sim_real1$br <= OLS4_beta & sim_real1$bd > OLS5_beta)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -0.1, y = -0.25, label = paste0(round(sum(sim_real1$br <= OLS4_beta & sim_real1$bd <= OLS5_beta)/iteration*100,1),"%"), color = "red")
  

#FIGURE 1B
  
  ggplot(sim_real1, aes(x = t_br, y = t_bd)) +
    ggtitle(expression(paste("Figure 1b - t-stats, ", phi, " = 0.94"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(paste("t, ", b[r]))) +
    ylab(expression(paste("t, ", b[d]))) +
    geom_jitter(size = 0.25) +
    geom_vline(xintercept = OLS4_tstat, color = "blue") +
    geom_hline(yintercept = OLS5_tstat, color = "blue") +
    geom_point(aes(x = OLS4_tstat, y = OLS5_tstat),shape = 16, size = 5, color = "blue") +
    annotate("text", x = 4, y = 1, label = paste0(round(sum(sim_real1$t_br > OLS4_tstat & sim_real1$t_bd > OLS5_tstat)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = 4, y = -7, label = paste0(round(sum(sim_real1$t_br > OLS4_tstat & sim_real1$t_bd <= OLS5_tstat)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -3, y = 1, label = paste0(round(sum(sim_real1$t_br <= OLS4_tstat & sim_real1$t_bd > OLS5_tstat)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -3, y = -7, label = paste0(round(sum(sim_real1$t_br <= OLS4_tstat & sim_real1$t_bd <= OLS5_tstat)/iteration*100,1),"%"), color = "red")
  
  
#FIGURE 1C
  
  ggplot(sim_real2, aes(x = br, y = bd)) +
    ggtitle(expression(paste("Figure 1c - Coefficients, ", phi, " = 0.99"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(b[r])) +
    ylab(expression(b[d])) +
    geom_jitter(size = 0.25) +
    geom_vline(xintercept = OLS4_beta, color = "blue") +
    geom_hline(yintercept = OLS5_beta, color = "blue") +
    geom_point(aes(x = 0, y = rho*0.99-1),shape = 18, size = 5, color = "red") +
    geom_point(aes(x = OLS4_beta, y = OLS5_beta),shape = 16, size = 5, color = "blue") +
    annotate("text", x = 0.3, y = 0.025, label = paste0(round(sum(sim_real2$br > OLS4_beta & sim_real2$bd > OLS5_beta)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = 0.3, y = -0.15, label = paste0(round(sum(sim_real2$br > OLS4_beta & sim_real2$bd <= OLS5_beta)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -0.05, y = 0.025, label = paste0(round(sum(sim_real2$br <= OLS4_beta & sim_real2$bd > OLS5_beta)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -0.05, y = -0.15, label = paste0(round(sum(sim_real2$br <= OLS4_beta & sim_real2$bd <= OLS5_beta)/iteration*100,1),"%"), color = "red")
  

#FIGURE 1D
  
  ggplot(sim_real2, aes(x = t_br, y = t_bd)) +
    ggtitle(expression(paste("Figure 1d - t-stats, ", phi, " = 0.99"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(paste("t, ", b[r]))) +
    ylab(expression(paste("t, ", b[d]))) +
    geom_jitter(size = 0.25) +
    geom_vline(xintercept = OLS4_tstat, color = "blue") +
    geom_hline(yintercept = OLS5_tstat, color = "blue") +
    geom_point(aes(x = OLS4_tstat, y = OLS5_tstat),shape = 16, size = 5, color = "blue") +
    annotate("text", x = 4, y = 0, label = paste0(round(sum(sim_real2$t_br > OLS4_tstat & sim_real2$t_bd > OLS5_tstat)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = 4, y = -6, label = paste0(round(sum(sim_real2$t_br > OLS4_tstat & sim_real2$t_bd <= OLS5_tstat)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -2, y = 0, label = paste0(round(sum(sim_real2$t_br <= OLS4_tstat & sim_real2$t_bd > OLS5_tstat)/iteration*100,1),"%"), color = "red") +
    annotate("text", x = -2, y = -6, label = paste0(round(sum(sim_real2$t_br <= OLS4_tstat & sim_real2$t_bd <= OLS5_tstat)/iteration*100,1),"%"), color = "red")
  
  
#TABLE 4
  
  Table4 = matrix(
    c(
      OLS4_beta/(1-rho*OLS6_beta), OLS5_beta/(1-rho*OLS6_beta), OLS7_beta/(1-rho*OLS6_beta),
      OLS4_se/(1-OLS6_beta^2), OLS5_se/(1-OLS6_beta^2), OLS7_se/(1-OLS6_beta^2),
      OLS4_beta/(1-rho*OLS6_beta)/(OLS4_se/(1-OLS6_beta^2)),
      (OLS5_beta/(1-rho*OLS6_beta)+1)/(OLS5_se/(1-OLS6_beta^2)),
      OLS7_beta/(1-rho*OLS6_beta)/(OLS7_se/(1-OLS6_beta^2)),
      2*pt(-abs(OLS4_beta/(1-rho*OLS6_beta)/(OLS4_se/(1-OLS6_beta^2))), df = data_length-1)*100,
      2*pt(-abs((OLS5_beta/(1-rho*OLS6_beta)+1)/(OLS5_se/(1-OLS6_beta^2))), df = data_length-1)*100,
      2*pt(-abs(OLS7_beta/(1-rho*OLS6_beta)/(OLS7_se/(1-OLS6_beta^2))), df = data_length-1)*100
    ),
    nrow = 3,
    ncol = 4
  )
  colnames(Table4) = c("$\\hat{b^{lr}}$", "$s. e.$", "$t$", paste("$\\%$", "$p$", "$value$"))
  rownames(Table4) = c("$r$", "$\\Delta d$", "Excess $r$")
  
  Table4 %>%
    kbl(row.names = TRUE,
        caption = "Table 4 - Long-run regression coefficients",
        escape = FALSE,
        digits = 2) %>%
    kable_styling("striped", full_width = FALSE)
  

#FIGURE 2A
  
  #phi = 0.94
  ggplot(sim_real1, aes(x = br/(1-rho*OLS6_beta))) +
    geom_histogram(binwidth = 0.1, color = "black", fill = "white") +
    ggtitle(expression(paste(phi, " = 0.94"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(paste(b[r], "/(1-", rho, phi, ")"))) +
    ylab("") +
    geom_vline(xintercept = OLS4_beta/(1-rho*OLS6_beta), color = "red") +
    annotate("text", x = 1.25, y = 500, label = "Data", color = "red")
  

#FIGURE 2B
  
  #phi = 0.99
  ggplot(sim_real2, aes(x = br/(1-rho*0.99))) +
    geom_histogram(binwidth = 0.1, color = "black", fill = "white") +
    ggtitle(expression(paste(phi, " = 0.99"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(paste(b[r], "/(1-", rho, phi, ")"))) +
    ylab("") +
    geom_vline(xintercept = OLS4_beta/(1-rho*0.99), color = "red") +
    annotate("text", x = 2.5, y = 300, label = "Data", color = "red")
  

#FIGURE 3A
  
  brlr = OLS4_beta/(1-rho*OLS6_beta)
  
  ggplot(sim_real1, aes(x = br, y = phi)) +
    ggtitle(expression(paste(b[r], " and ", phi, ", ", phi, " = 0.94"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(b[r])) +
    ylab(expression(phi)) +
    geom_jitter(size = 0.25) +
    geom_vline(xintercept = OLS4_beta, color = "blue") +
    geom_hline(yintercept = OLS6_beta, color = "blue") +
    geom_line(aes(x=1-rho*sim_real1$phi+OLS5_beta, y=sim_real1$phi), color = "red", linetype = "dashed") +
    geom_line(aes(x=brlr*(1-rho*sim_real1$phi), y=sim_real1$phi), color = "blue") +
    geom_point(aes(x = 0, y = OLS6_beta),shape = 18, size = 5, color = "red") +
    geom_point(aes(x = OLS4_beta, y = OLS6_beta),shape = 16, size = 5, color = "blue") +
    annotate("text", x = 0.4, y = 1, label = paste0(round(sum(sim_real1$br > OLS4_beta & sim_real1$phi > OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = 0.4, y = 0.5, label = paste0(round(sum(sim_real1$br > OLS4_beta & sim_real1$phi <= OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = -0.1, y = 1, label = paste0(round(sum(sim_real1$br <= OLS4_beta & sim_real1$phi > OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = -0.1, y = 0.5, label = paste0(round(sum(sim_real1$br <= OLS4_beta & sim_real1$phi <= OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = 0.35, y = 0.7, label = expression(b[r]^{lr}), color = "blue") +
    annotate("text", x = 0.3, y = 0.65, label = expression(b[d]), color = "red")
  

#FIGURE 3B
  
  brlr = OLS4_beta/(1-rho*OLS6_beta)
  
  ggplot(sim_real2, aes(x = br, y = phi)) +
    ggtitle(expression(paste(b[r], " and ", phi, ", ", phi, " = 0.99"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(b[r])) +
    ylab(expression(phi)) +
    geom_jitter(size = 0.25) +
    geom_vline(xintercept = OLS4_beta, color = "blue") +
    geom_hline(yintercept = OLS6_beta, color = "blue") +
    geom_point(aes(x = 0, y = 0.99),shape = 18, size = 5, color = "red") +
    geom_point(aes(x = OLS4_beta, y = OLS6_beta),shape = 16, size = 5, color = "blue") +
    geom_line(aes(x=1-rho*sim_real2$phi+OLS5_beta, y=sim_real2$phi), color = "red", linetype = "dashed") +
    geom_line(aes(x=brlr*(1-rho*sim_real2$phi), y=sim_real2$phi), color = "blue") +
    annotate("text", x = 0.4, y = 1, label = paste0(round(sum(sim_real2$br > OLS4_beta & sim_real2$phi > OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = 0.4, y = 0.5, label = paste0(round(sum(sim_real2$br > OLS4_beta & sim_real2$phi <= OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = -0.15, y = 1, label = paste0(round(sum(sim_real2$br <= OLS4_beta & sim_real2$phi > OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = -0.15, y = 0.5, label = paste0(round(sum(sim_real2$br <= OLS4_beta & sim_real2$phi <= OLS6_beta)/iteration*100,1),"%")) +
    annotate("text", x = 0.4, y = 0.65, label = expression(b[r]^{lr}), color = "blue") +
    annotate("text", x = 0.3, y = 0.65, label = expression(b[d]), color = "red")
  

#TABLE 5

  iteration = 5000
  data_length = nrow(data)
  
  #Real
  cov_matrix = cov(cbind(OLS6$resid, OLS5$resid))
  
    #Run Simulation 3
    phi = 0.9
    sim_real3 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 4
    phi = 0.96
    sim_real4 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 5
    phi = 0.98
    sim_real5 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 6
    phi = 1
    sim_real6 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 7
    phi = 1.01
    sim_real7 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    
  #Excess
  cov_matrix = cov(cbind(OLS6$resid, OLS8$resid))
  phi = OLS6_beta
  
    #Run Simulation 2
    phi = 0.9
    sim_excess2 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 3
    phi = 0.96
    sim_excess3 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 4
    phi = 0.98
    sim_excess4 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 5
    phi = 0.99
    sim_excess5 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 6
    phi = 1
    sim_excess6 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 7
    phi = 1.01
    sim_excess7 = MCSim(iteration, phi, rho, data_length, cov_matrix)
    
  Table5 = matrix(
    c(
      round(sum(sim_real3$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real1$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real4$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real5$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real2$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real6$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real7$br > OLS4_beta)/iteration*100,1),
      round(sum(sim_real3$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_real1$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_real4$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_real5$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_real2$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_real6$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_real7$bd > OLS5_beta)/iteration*100,1),
      round(sum(sim_excess2$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess1$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess3$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess4$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess5$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess6$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess7$br > OLS7_beta)/iteration*100,1),
      round(sum(sim_excess2$bd > OLS8_beta)/iteration*100,1),
      round(sum(sim_excess1$bd > OLS8_beta)/iteration*100,1),
      round(sum(sim_excess3$bd > OLS8_beta)/iteration*100,1),
      round(sum(sim_excess4$bd > OLS8_beta)/iteration*100,1),
      round(sum(sim_excess5$bd > OLS8_beta)/iteration*100,1),
      round(sum(sim_excess6$bd > OLS8_beta)/iteration*100,1),
      round(sum(sim_excess7$bd > OLS8_beta)/iteration*100,1),
      round(OLS6_sd/sqrt(1-0.9^2),2),
      round(OLS6_sd/sqrt(1-OLS6_beta^2),2),
      round(OLS6_sd/sqrt(1-0.96^2),2),
      round(OLS6_sd/sqrt(1-0.98^2),2),
      round(OLS6_sd/sqrt(1-0.99^2),2),
      "$\\infty$",
      "$\\infty$",
      round(log(0.5, base = 0.9),1),
      round(log(0.5, base = OLS6_beta),1),
      round(log(0.5, base = 0.96),1),
      round(log(0.5, base = 0.98),1),
      round(log(0.5, base = 0.99),1),
      "$\\infty$",
      "$\\infty$"
    ),
    nrow = 7,
    ncol = 6
  )
  
  colnames(Table5) = c("$b_r$", "$b_d$", "$b_r$", "$b_d$", "$\\sigma (dp)$", "1/2 life")
  rownames(Table5) = c("0.90", round(OLS6_beta, 3), "0.96", "0.98", "0.99", "1", "1.01")
  
  Table5 %>%
    kbl(row.names = TRUE,
        caption = paste("Table 5 - The effects of dividend-yield autocorrelation", "$\\phi$"),
        escape = FALSE) %>%
    add_header_above(c("Null $\\phi$"=1, "Real returns"=2, "Excess returns"=2, "Statistics"=2)) %>%
    add_header_above(c(" "=1, "Percentage probability values"=4, "Other"=2)) %>%
    kable_styling("striped", full_width = FALSE)

    
#FIGURE 5
    
  br5 = (1 + (rho*OLS6_beta)^1 + (rho*OLS6_beta)^2 + (rho*OLS6_beta)^3 + (rho*OLS6_beta)^4)*OLS4_beta
  br10 = (1 + (rho*OLS6_beta)^1 + (rho*OLS6_beta)^2 + (rho*OLS6_beta)^3 + (rho*OLS6_beta)^4 + (rho*OLS6_beta)^5 + (rho*OLS6_beta)^6 + (rho*OLS6_beta)^7 + (rho*OLS6_beta)^8 + (rho*OLS6_beta)^9)*OLS4_beta
  br20 = (1 + (rho*OLS6_beta)^1 + (rho*OLS6_beta)^2 + (rho*OLS6_beta)^3 + (rho*OLS6_beta)^4 + (rho*OLS6_beta)^5 + (rho*OLS6_beta)^6 + (rho*OLS6_beta)^7 + (rho*OLS6_beta)^8 + (rho*OLS6_beta)^9 + (rho*OLS6_beta)^10 + (rho*OLS6_beta)^11 + (rho*OLS6_beta)^12 + (rho*OLS6_beta)^13 + (rho*OLS6_beta)^14 + (rho*OLS6_beta)^15 + (rho*OLS6_beta)^16 + (rho*OLS6_beta)^17 + (rho*OLS6_beta)^18 + (rho*OLS6_beta)^19)*OLS4_beta
  brlr = OLS4_beta/(1-rho*OLS6_beta)
  brlr_unweighted = OLS4_beta/(1-OLS6_beta)
  
  ggplot(sim_real1, aes(x = br, y = phi)) +
    ggtitle(expression(paste(b[r], " and ", phi, ", with long run regressions"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(b[r])) +
    ylab(expression(phi)) +
    geom_jitter(size = 0.25) +
    geom_line(aes(x=OLS4_beta, y=sim_real1$phi), color = "blue") +
    geom_line(aes(x=br5/(1 + (rho*sim_real1$phi)^1 + (rho*sim_real1$phi)^2 + (rho*sim_real1$phi)^3 + (rho*sim_real1$phi)^4), y=sim_real1$phi), color = "blue") +
    geom_line(aes(x=br10/(1 + (rho*sim_real1$phi)^1 + (rho*sim_real1$phi)^2 + (rho*sim_real1$phi)^3 + (rho*sim_real1$phi)^4 + (rho*sim_real1$phi)^5 + (rho*sim_real1$phi)^6 + (rho*sim_real1$phi)^7 + (rho*sim_real1$phi)^8 + (rho*sim_real1$phi)^9), y=sim_real1$phi), color = "blue") +
    geom_line(aes(x=br20/(1 + (rho*sim_real1$phi)^1 + (rho*sim_real1$phi)^2 + (rho*sim_real1$phi)^3 + (rho*sim_real1$phi)^4 + (rho*sim_real1$phi)^5 + (rho*sim_real1$phi)^6 + (rho*sim_real1$phi)^7 + (rho*sim_real1$phi)^8 + (rho*sim_real1$phi)^9 + (rho*sim_real1$phi)^10 + (rho*sim_real1$phi)^11 + (rho*sim_real1$phi)^12 + (rho*sim_real1$phi)^13 + (rho*sim_real1$phi)^14 + (rho*sim_real1$phi)^15 + (rho*sim_real1$phi)^16 + (rho*sim_real1$phi)^17 + (rho*sim_real1$phi)^18 + (rho*sim_real1$phi)^19), y=sim_real1$phi), color = "blue") +
    geom_line(aes(x=brlr*(1-rho*sim_real1$phi), y=sim_real1$phi), color = "blue") +
    geom_line(aes(x=brlr_unweighted*(1-sim_real1$phi), y=sim_real1$phi), color = "blue", linetype = "dashed") +  
    geom_point(aes(x = 0, y = OLS6_beta),shape = 18, size = 5, color = "red") +
    geom_point(aes(x = OLS4_beta, y = OLS6_beta),shape = 16, size = 5, color = "blue") +
    annotate(geom="label", x = OLS4_beta, y = 0.6, label = "1", color = "blue", fill = "white") +
    annotate(geom="label", x = 0.165, y = 0.6, label = "5", color = "blue", fill = "white") +
    annotate(geom="label", x = 0.25, y = 0.6, label = "10", color = "blue", fill = "white") +
    annotate(geom="label", x = 0.345, y = 0.6, label = "20", color = "blue", fill = "white") +
    annotate(geom="label", x = 0.405, y = 0.6, label = expression(infinity), color = "blue", fill = "white") +
    annotate(geom="label", x = 0.45, y = 0.7, label = expression(paste(infinity, ", unweighted")), color = "blue", fill = "white")
  
  
#Goyal-Welch Statistics
  
  #Data's delta_RMSE
  FR = matrix(NA, data_length, 2)
  
  for (i in 20:data_length)
  {
    ret = ts(log(data[["Ret"]]), start = 1926, end = 1926 + i - 1)
    dp = ts(log(data[["DP"]]), start = 1926, end = 1926 + i - 1)
    
    #Forecast residual in the next period using mean:
    FR[i-19, 1] = log(data[["Ret"]][i+1]) - mean(ret)
    
    #Forecast residual in the next period using Prediction Model
    reg = dyn$lm(ret ~ lag(dp, -1))
    FR[i-19, 2] = log(data[["Ret"]][i+1]) - (reg$coef[1]+reg$coef[2]*log(data[["DP"]][i]))
  }
  
  delta_RMSE_sample = sqrt(mean(FR[,1]^2, na.rm = TRUE)) - sqrt(mean(FR[,2]^2, na.rm = TRUE))
  
  #Monte-Carlo Simulation
  
    #Simulation function
    MCSim_GW = function(iteration, phi, rho, data_length, cov_matrix)
    {
      delta_RMSE_sim = matrix(NA, nrow = iteration, ncol = 1)
      
      for (i in 1:iteration)
      {
        #Generate errors from bi-variate normal
        sim_err = rmvnorm(data_length, c(0, 0), cov_matrix)
        
        #Generate dividend price
        sim_dp = rep(0, data_length)
        
        #Generate initial observation
        sim_dp[1] = rnorm(1, 0, sqrt(cov_matrix[1, 1]/(1-phi^2)))
        
        #Generate the rest of the observation
        for (j in 2:data_length)
        {
          sim_dp[j] = phi * sim_dp[j-1] + sim_err[j,1]
        }
        
        #Generate dividend growth
        sim_dg = sim_err[, 2]
        
        #Generate return
        sim_ret = rep(0, data_length)
        
        #Generate initial observation
        sim_ret[1] = rnorm(1, 0, sd(OLS4$resid))
        
        #Generate the rest of the observation
        for (j in 2:data_length)
        {
          sim_ret[j] = (phi*rho-1) * sim_dp[j-1] + (sim_err[j, 1] - rho*sim_err[j, 2])
        }
        
        #Out of sample RMSE
        FR_sim = matrix(NA, data_length, 2)
        
        for (k in 20:data_length)
        {
          ret = ts(sim_ret, start = 1, end = k)
          dp = ts(sim_dp, start = 1, end = k)
          
          #Forecast residual in the next period using mean:
          FR_sim[k-19, 1] = sim_ret[k+1] - mean(ret)
          
          #Forecast residual in the next period using Prediction Model
          reg = dyn$lm(ret ~ lag(dp, -1))
          FR_sim[k-19, 2] = sim_ret[k+1] - (reg$coef[1]+reg$coef[2]*sim_dp[k])
        }
        
        delta_RMSE_sim[i] = sqrt(mean(FR_sim[,1]^2, na.rm = TRUE)) - sqrt(mean(FR_sim[,2]^2, na.rm = TRUE))
      }
      
      delta_RMSE_sim = as.data.frame(delta_RMSE_sim)
      colnames(delta_RMSE_sim) = "delta_RMSE_sim"
      return(delta_RMSE_sim)
    }
    
    #Simulation parameters (general)
    iteration = 5000
    data_length = nrow(data)
    cov_matrix = cov(cbind(OLS6$resid, OLS5$resid))
    
    #Run Simulation 1
    phi = OLS6_beta
    simGW_real1 = MCSim_GW(iteration, phi, rho, data_length, cov_matrix)
    #Run Simulation 2
    phi = 0.99
    simGW_real2 = MCSim_GW(iteration, phi, rho, data_length, cov_matrix)
    

#FIGURE 6A
    
  #phi = 0.94
  ggplot(simGW_real1, aes(x = delta_RMSE_sim)) +
    geom_histogram(color = "black", fill = "white") +
    ggtitle(expression(paste(phi, " = 0.94"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(paste(Delta, " rmse"))) +
    ylab("") +
    geom_vline(xintercept = delta_RMSE_sample, color = "red") +
    annotate("text", x = -0.0075, y = 700, label = "Data", color = "red") +
    annotate("text", x = -0.015, y = 200, label = paste(round(sum(simGW_real1 < delta_RMSE_sample)/iteration*100, 0), "%"), color = "red")
  
  
#FIGURE 6B  
  
  #phi = 0.99
  ggplot(simGW_real2, aes(x = delta_RMSE_sim)) +
    geom_histogram(color = "black", fill = "white") +
    ggtitle(expression(paste(phi, " = 0.99"))) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(expression(paste(Delta, " rmse"))) +
    ylab("") +
    geom_vline(xintercept = delta_RMSE_sample, color = "red") +
    annotate("text", x = -0.0075, y = 700, label = "Data", color = "red") +
    annotate("text", x = -0.015, y = 200, label = paste(round(sum(simGW_real2 < delta_RMSE_sample)/iteration*100, 0), "%"), color = "red")