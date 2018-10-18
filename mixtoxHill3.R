library(dplyr)
library(tidyr) #seperate
library(minpack.lm)

df <- read.csv("mixtoxtest.csv")
colnames(df)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", 
                "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$dosen <- as.numeric(DF$dose)
DF1 <- DF %>% mutate(dose = 10^(dosen-3))
colnames(DF1)[4]<-"response"

chem <- c("DDT, O,P'-", "ALDRIN", "DIELDRIN", "CADMIUM(Chloride)", "HEPTACHLOR",
          "DDD, P,P'-", "MERCURIC CHLORIDE", "ENDOSULFAN", "DICOFOL", "CHLORPYRIFOS")

df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
df <- data.frame(df[,c(2,6:13)])
names(df) <- c("chemical","AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

df[37,1] <- "CADMIUM(Chloride)"

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
  
  plot(log(DF$dose, 10), DF$response/100, pch = 19, col="darkgrey", main = chem[i],
       xlab=expression(paste("Log", "Conc. (", mu,"M)")), ylab="Response (%)")
  lines(px,py, col = 1, lwd = 2)
  lines(px,CI.up, col = 2, lwd = 1, lty = 2)
  lines(px,CI.low, col = 2, lwd = 1, lty = 2)
  lines(px,PI.up, col = 4, lwd = 1, lty = 2)
  lines(px,PI.low, col = 4, lwd = 1, lty = 2)
  min <- log(df[which(df[,1] == chem[i]), "AC50 min"], 10)
  max <- log(df[which(df[,1] == chem[i]), "AC50 max"], 10)   
  polygon(c(min, max, max, min), c(2 , 2, -1, -1), col=rgb(1, 0, 0,0.1), border=NA)

  print(list(p = paramHat, sta = sta))
}

png(file="ind-DR.png",width=4800,height=2000,res=300)
par(mfrow = c(2,5))
plotFit(1)
plotFit(2)
plotFit(3)
plotFit(4)
plotFit(5)
plotFit(6)
plotFit(7)
plotFit(8, init_n = 6)
plotFit(9)
plotFit(10)
dev.off()




