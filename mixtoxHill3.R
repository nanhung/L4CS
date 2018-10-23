library(dplyr)
library(tidyr) #seperate
library(minpack.lm)

df1 <- read.csv("mixtoxtest.csv")
sheets <- readxl::excel_sheets("Mixture_Neuron.xlsx")
df2 <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
df1[,1] <- df2[,2]
colnames(df1)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", 
                 "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
DF <- df1 %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$dosen <- as.numeric(DF$dose)
DF1 <- DF %>% mutate(dose = 10^(dosen-3))
colnames(DF1)[4]<-"response"

df2 <- data.frame(df2[,c(2,6:13)])
names(df2) <- c("chemical","AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

chem <- df1[,1]



Fit <- function(i, init_n = 1){
  DF <- DF1 %>% filter(chemical == chem[i]) %>% group_by(round) 
  
  x <- DF$dose
  y <- DF$response/100
  
  n <- length(DF$dose)
  m <- 3
  
  Alpha <- 20
  Beta <- 1
  param <- c(Alpha, Beta)
  
  fun <- 'y ~ Beta / (1 + x/Alpha)'
  
  dframe <- data.frame(x, y)
  fit1 <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]),
                            control = nls.lm.control(maxiter = 1000))
  
  fun <- paste('y ~', coef(fit1)[2], '/ (1 + (x/Alpha)^Beta)')
  fit2 <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = coef(fit1)[1], Beta = init_n),
                            control = nls.lm.control(maxiter = 1000))
  
  paramHat <- c(coef(fit2)[1], coef(fit2)[2], coef(fit1)[2])
  names(paramHat) <- c("EC50", "n", "Emax")
  
  fitInfo <- summary(fit2) # fitting information	
  yhat <- paramHat[3] / (1 + (x/paramHat[1])^paramHat[2])
  
  sst <- sum((y - mean(y))^2) # total sum of squares
  sse <- sum((y - yhat)^2) # sum of squared errors
  r2 <- 1 - sse / sst # coefficient of determination
  adjr2 <- 1 - sse * (n - 1) / (sst * (n - m)) # adjusted coefficient of determination
  rmse <- sqrt(sse / (n - m)) # root-mean-square error
  mae <- sum(abs(y - yhat)) / n # mean absolute error
  lnL <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sse)))
  aic <- 2 * m - 2 * lnL # Akaike information criterion 
  aicc <- aic + 2 * m * (m + 1) / (n - m - 1)
  bic <- m * log(n) - 2 * lnL # Bayesian information criterion
  sta <- t(c(r2, adjr2, mae, rmse, aic, aicc, bic))
  colnames(sta) <- c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC')
  
  list(p = paramHat, sta = sta)  
}


plotFit <- function(i, init_n = 1){

  DF <- DF1 %>% filter(chemical == chem[i]) %>% group_by(round) 
  
  x <- DF$dose
  y <- DF$response/100
  
  n <- length(DF$dose)
  m <- 3
  
  Alpha <- 20
  Beta <- 1
  param <- c(Alpha, Beta)
  
  fun <- 'y ~ Beta / (1 + x/Alpha)'
  
  dframe <- data.frame(x, y)
  fit1 <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]),
                            control = nls.lm.control(maxiter = 1000))
  
  fun <- paste('y ~', coef(fit1)[2], '/ (1 + (x/Alpha)^Beta)')
  fit2 <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = coef(fit1)[1], Beta = init_n),
                            control = nls.lm.control(maxiter = 1000))
  
  paramHat <- c(coef(fit2)[1], coef(fit2)[2], coef(fit1)[2])
  names(paramHat) <- c("EC50", "n", "Emax")
  
  fitInfo <- summary(fit2) # fitting information	
  yhat <- paramHat[3] / (1 + (x/paramHat[1])^paramHat[2])
  
  sst <- sum((y - mean(y))^2) # total sum of squares
  sse <- sum((y - yhat)^2) # sum of squared errors
  r2 <- 1 - sse / sst # coefficient of determination
  adjr2 <- 1 - sse * (n - 1) / (sst * (n - m)) # adjusted coefficient of determination
  rmse <- sqrt(sse / (n - m)) # root-mean-square error
  mae <- sum(abs(y - yhat)) / n # mean absolute error
  lnL <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(sse)))
  aic <- 2 * m - 2 * lnL # Akaike information criterion 
  aicc <- aic + 2 * m * (m + 1) / (n - m - 1)
  bic <- m * log(n) - 2 * lnL # Bayesian information criterion
  sta <- t(c(r2, adjr2, mae, rmse, aic, aicc, bic))
  colnames(sta) <- c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC', 'AICc', 'BIC')
  
  list(p = paramHat, sta = sta)  
  
  px <- seq(-2,2,0.1)
  py <- paramHat[3] / (1 + (10^px/paramHat[1])^paramHat[2])
  
  sigLev = 0.05
  probT <- qt(1 - sigLev / 2, n - m) # the student t distribution
  mse <- rmse^2  # squared residual standard error
  jac <- mixtox::jacobian(eq = 'Hill_three', 10^px, paramHat)
  covPara <- mse * solve(t(jac) %*% jac)  # covariance matrix of the parameter estimates
  
  gap.PI <- sqrt(mse + diag(jac %*% covPara %*% t(jac))) # prediction intervals
  gap.CI <- sqrt(diag(jac %*% covPara %*% t(jac))) # confidence intervals
  
  PI.up <- py + probT * gap.PI # PI upper bound
  PI.low <- py - probT * gap.PI # PI lower bound
  CI.up <- py + probT * gap.CI # CI upper bound
  CI.low <- py - probT * gap.CI # CI lower bound
  crcInfo <- cbind(px, py, PI.low, PI.up, CI.low, CI.up)
  
  plot(log(DF$dose, 10), DF$response/100, pch = 19, col="darkgrey", main = paste(chem[i], "; r2 =", round(sta[1],2)),
       xlab=expression(paste("Log", "Conc. (", mu,"M)")), ylab="Response (%)")
  lines(px,py, col = 1, lwd = 2)
  lines(px,CI.up, col = 2, lwd = 1, lty = 2)
  lines(px,CI.low, col = 2, lwd = 1, lty = 2)
  lines(px,PI.up, col = 4, lwd = 1, lty = 2)
  lines(px,PI.low, col = 4, lwd = 1, lty = 2)
  #
  min <- log(df2[which(df2[,1] == chem[i]), "POD min"], 10)
  max <- log(df2[which(df2[,1] == chem[i]), "POD max"], 10)   
  #
  polygon(c(min, max, max, min), c(2 , 2, -1, -1), col=rgb(1, 0, 0,0.1), border=NA)

  print(list(p = paramHat, sta = sta))
}

png(file="ind-DR.png",width=4800,height=2800,res=300)
par(mfrow = c(4,6))
plotFit(4)
plotFit(5)
plotFit(6)
plotFit(7)
plotFit(12)
plotFit(15)
plotFit(16)
plotFit(18)
plotFit(19)
plotFit(20)
plotFit(22)
plotFit(23)
plotFit(25)
plotFit(27, init_n = 6)
plotFit(28)
plotFit(31)
plotFit(33)
plotFit(36)
plotFit(37)
plotFit(39)
plotFit(40)
plotFit(41)
plotFit(42)
dev.off()


###########################

source("mixECx.R")

effv_est <- function(i, init_n = 1){
  X <- Fit(i, init_n)
  C <- df2[i,'POD max']
  y <- X$p[3] / (1 + (C/X$p[1])^X$p[2])
  return(y)  
}

sheets <- readxl::excel_sheets("Mixture_Neuron.xlsx")
df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
df$ugl <- df$`Molecular weight`* df$POD.Highest
total_weight <- sum(df$ugl)
effv <- df$ugl/total_weight

effPoints <- rev((c(0.025, 0.03, 0.05, 0.1, 0.15, 0.2, 
                    0.25, 0.3, 0.35, 0.4, 0.45, 0.47, 0.5, 0.52, 
                    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.99)))
pctEcx <- t(t(effv/sum(effv)))

param0 <- c(1000, 1, 1.1)

model <- rep('Hill_three_rev', 42)
param <- matrix(c(rep(param0, 3),
                  Fit(4)$p,
                  Fit(5)$p,
                  Fit(6)$p,
                  Fit(7)$p,
                  rep(param0, 4),
                  Fit(12)$p,
                  rep(param0, 2),
                  Fit(15)$p,
                  Fit(16)$p,
                  rep(param0, 1),
                  Fit(18)$p,
                  Fit(19)$p,
                  Fit(20)$p,
                  rep(param0, 1),
                  Fit(22)$p,
                  Fit(23)$p,
                  rep(param0, 1),
                  Fit(25)$p,
                  rep(param0, 1),
                  Fit(27, init_n = 6)$p,
                  Fit(28)$p,
                  rep(param0, 2),
                  Fit(31)$p,
                  rep(param0, 1),
                  Fit(33)$p,
                  rep(param0, 2),
                  Fit(36)$p,
                  Fit(37)$p,
                  rep(param0, 1),
                  Fit(39)$p,
                  Fit(40)$p,
                  Fit(41)$p,
                  Fit(42)$p), byrow = T, ncol =3)

ECx(model, param, effPoints)

concAdd <- function(pctEcx, effPoints) {
  ecPoints <- ECx(model, param, effPoints)
  ca <- 1/(t(pctEcx) %*% (1/ecPoints))
  return(ca)
}

ca <- concAdd(pctEcx, rev(effPoints))
plot(ca, effPoints, log = "x", type = "l")
              