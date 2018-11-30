ia_estimate <- function(data){
  eff <- data
  effv <- eff / sum(eff)
  pctEcx <- t(t(effv/sum(effv)))
  
  ia <- indAct(model, param, pctEcx, effPoints)
  ia <- c(1e-4, as.numeric(ia))
  return(ia)
}

indAct <- function(model, param, pctEcx, effPoints){
  # independent action
  ecPoints <- ECx(model, param, effPoints)
  fac <- nrow(pctEcx)
  lev <- ncol(pctEcx)
  iaFun <- as.character(rep(1, lev))
  
  for (i in seq(lev)){
    # IA equation construction
    # xx means x elsewhere
    for (j in seq(fac)){
      iaFun[i] <- paste(iaFun[i], '*', '(1 - (', param[j, 3] ,'/ (1 + (', param[j, 1], '/', pctEcx[j, i], '* xx)^', param[j, 2], ')))', sep = '')
    }
  }
  
  a <- 1e-9
  b <- 6
  eps <- 1e-10		
  root <- matrix(0, lev, ncol(ecPoints))
  
  for (i in seq(lev)){
    fia <-  iaFun[i]
    for (k in seq(ncol(ecPoints))){
      value <- 1 - effPoints[k]
      fun <- paste(value, '-',  fia, sep = '')
      f = function(xx) eval(parse(text = fun))
      root[i, k] <- uniroot(f, c(a, b), tol = eps)$root
    }
  }
  
  colName <- paste('EC', effPoints * 100, sep = '')
  colnames(root) <- colName
  return(root)
}

ca_estimate <- function(data){
  eff <- data
  effv <- eff / sum(eff)
  pctEcx <- t(t(effv/sum(effv)))
  
  concAdd <- function(pctEcx, effPoints) {
    ecPoints <- ECx(model, param, effPoints)
    ca <- 1/(t(pctEcx) %*% (1/ecPoints))
    return(ca)
  }
  
  ca <- concAdd(pctEcx, rev(effPoints))
  ca <- c(1e-4, as.numeric(ca))
  
  return(ca)
}


Fit <- function(i, init_n = 1, effect = sheets[3]){
  
  df <- readxl::read_xlsx("42_Chem_Neuron.xlsx", sheet = effect) 
  colnames(df)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", 
                   "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
  DF <- df %>% as.data.frame() %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
  DF$dosen <- as.numeric(DF$dose)
  DF2 <- DF %>% mutate(dose = 10^(dosen-3))
  colnames(DF2)[4]<-"response"
  
  DF <- DF2 %>% filter(chemical == chem[i,]) %>% group_by(round) 
  
  x <- DF$dose
  y <- DF$response/100
  
  n <- length(DF$dose)
  m <- 3
  
  Alpha <- 20
  Beta <- 1
  param <- c(Alpha, Beta)
  
  fun <- 'y ~ Beta / (1 + x/Alpha)'
  
  dframe <- data.frame(x, y)
  fit1 <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = param[1], Beta = param[2]))
  
  fun <- paste('y ~', coef(fit1)[2], '/ (1 + (x/Alpha)^Beta)')
  fit2 <- minpack.lm::nlsLM(fun, data = dframe, start = list(Alpha = coef(fit1)[1], Beta = init_n))
  
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





ECx <- function(model, param, effv, rtype = 'continuous', Scaled = TRUE){
  #calculate effect concentrations using associated inverse function
  if (missing(model) || missing (param)) stop('argument missing')
  if (missing(effv)) stop('error! input effv in ECx')
  if (is.vector(param)) param <- t(param)
  #effv <- sort(effv)
  
  if(max(effv) >= 1.0){
    rtype <- 'continuous'
    Scaled <- FALSE
  }
  
  ecx <- matrix(0, length(model), length(effv))
  
  if((rtype == 'continuous' || rtype == 'ctn') && Scaled == TRUE){
    rspnRange <- matrix(0, length(model), 2)
    rspnRange <- CEx(model, param, c(0, 1e20))
    effvAbs <- matrix(0, length(model), length(effv))
    
    for(j in seq(model)) effvAbs[j, ] <- rspnRange[j , 1] + (rspnRange[j , 2] - rspnRange[j , 1]) * effv
  }
  
  effv0 <- effv
  
  for (i in seq(model)){
    fun <- model[i]
    p <- param[i, ]
    
    if((rtype == 'continuous' || rtype == 'ctn') && Scaled == TRUE) effv0 <- effvAbs[i, ]
    
    ec <- switch(fun,
                 'Hill_three_rev' =  p[1] * ((p[3] / effv0 - 1)^(1 / p[2])),
                 'Hill' = p[1] / ((1 / effv0 - 1)^(1 / p[2])),
                 'Hill_two' = p[1] * effv0 / (p[2] - effv0),
                 'Hill_three' = p[1] / ((p[3] / effv0 - 1)^(1 / p[2])),			
                 'Hill_four' = p[1] / (((p[3] - p[4]) / (effv0 - p[4]) - 1)^(1 / p[2])),
                 'Weibull' = exp(-(-log(log(-1 / (-1 + effv0))) + p[1]) * log(10) / p[2]),
                 'Weibull_three' = exp(-(-log(log(p[3] / (p[3] - effv0))) + p[1]) * log(10) / p[2]),
                 'Weibull_four' = exp((log(log((-p[4] + p[3]) / (p[3] - effv0))) - p[1]) * log(10) / p[2]),
                 "Logit" = exp(-log(10) * (p[1] + log(-(-1 + effv0) / (effv0))) / p[2]),
                 'Logit_three' = exp(-log(10) * (p[1] + log((p[3] - effv0) / effv0)) / p[2]),
                 'Logit_four' = exp(-log(10) * (p[1] + log(-(p[3] - effv0) / (p[4] - effv0))) / p[2]),
                 "BCW" = exp(log(-(p[1] * p[3] - p[2] - log(-log(1 - effv0)) * p[3]) / p[2]) / p[3]),
                 "BCL" = exp(log(-(p[1] * p[3] - p[2] + log(-(-1 + effv0) / effv0) * p[3]) / p[2]) / p[3]),
                 "GL" = exp(-log(10) * (p[1] + log(exp(-log(effv0) / p[3]) - 1)) / p[2])
    )
    
    ecx[i, ] <- ec
  }
  
  if(rtype == 'quantal'){
    colName <- paste0('EC', effv * 100)
    
  }else if(rtype == 'continuous' || rtype == 'ctn'){
    
    if(Scaled == FALSE) colName <- paste0('EC', effv)
    
    if(Scaled == TRUE){
      colName <- paste0('EC', effv * 100)
      colNameAbs <- paste0('Abs_rspn@E', effv * 100)
      colnames(effvAbs) <- colNameAbs
      
      if(is.null(rownames(param))) rownames(effvAbs) <- model else rownames(effvAbs) <- rownames(param)
    }
  }
  
  colnames(ecx) <- colName
  
  if(is.null(rownames(param))) rownames(ecx) <- model else rownames(ecx) <- rownames(param)
  
  #if((rtype == 'continuous' || rtype == 'ctn') && Scaled == TRUE){
  #  list(ecx = ecx, effvAbs = effvAbs)
  #}else{
    return(ecx)
  #}
}

CEx <- function(model, param, conc){
  # calculate response based on concentration
  if (missing(model) || missing (param) || missing(conc)) stop('argument missing')
  #if (missing(conc)) conc = 0.00005
  if (is.vector(param)) param <- t(param)
  
  effv <- matrix(0, length(model), length(conc))
  
  for (i in seq(model)){
    fun <- model[i]
    p <- param[i, ]
    
    for (j in seq(conc)){
      if (fun == 'Hill')
        ev <- 1 / (1 + (p[1] / conc[j])^p[2])
      else if (fun == 'Hill_two')
        ev <- p[2] * conc[j] / (p[1] + conc[j])
      else if (fun == 'Hill_three')
        ev <- p[3] /(1 + (p[1] / conc[j])^p[2])
      else if (fun == 'Hill_three_rev')
        ev <- p[3] /(1 + (conc[j] / p[1])^p[2])
      else if (fun == 'Hill_four')
        ev <- p[4] + (p[3] - p[4]) / (1 + (p[1] / conc[j])^p[2])
      else if(fun == 'Weibull')
        ev <- 1 - exp(-exp(p[1] + p[2] * log10(conc[j])))
      else if(fun == 'Weibull_three')
        ev <- p[3] * (1 - exp(-exp(p[1] + p[2] * log10(conc[j]))))
      else if(fun == 'Weibull_four')
        ev <- p[3] + (p[4] - p[3]) * exp(-exp(p[1] + p[2] * log10(conc[j])))
      else if (fun == "Logit")
        ev <- 1 / (1 + exp(-p[1] - p[2] * log10(conc[j])))
      else if(fun == 'Logit_three')
        ev <- p[3] / (1 + exp((-p[1]) - p[2] * log10(conc[j])))
      else if(fun == 'Logit_four')
        ev <- p[4] + (p[3] - p[4]) / (1 + exp((-p[1]) - p[2] * log10(conc[j])))
      else if (fun == "BCW")
        ev <- 1 - exp(-exp(p[1] + p[2] * ((conc[j]^p[3] - 1) / p[3])))
      else if (fun == "BCL")
        ev <- 1 / (1 + exp(-p[1] - p[2]((conc[j]^p[3] - 1) / p[3])))
      else if (fun == "GL")
        ev <- 1 / (1 + exp(-p[1] - p[2] * log10(conc[j])))^p[3]
      else if (fun == "Brain_Consens") 
        ev <- 1 - (1 + p[1] * conc[j]) / (1 + exp(p[2] * p[3]) * conc[j]^p[2])
      else if(fun == "BCV") 
        ev <- 1 - p[1] * (1 + p[2] * conc[j]) / (1 + (1 + 2 * p[2] * p[3]) * (conc[j] / p[3])^p[4])
      else if(fun == "Cedergreen") 
        ev <- 1 - (1 + p[1] * exp(-1 / (conc[j]^p[2]))) / (1 + exp(p[3] * (log(conc[j]) - log(p[4]))))
      else if(fun == "Beckon") 
        ev <- (p[1] + (1 - (p[1]) / (1 + (p[2] / conc[j])^p[3]))) / (1 + (conc[j] / p[4])^p[5])
      else if(fun == "Biphasic") 
        ev <- p[1] - p[1] / (1 + 10^((conc[j] - p[2]) * p[3])) + (1 - p[1]) / (1 + 10^((p[4] - conc[j]) * p[5]))
      else if(fun == 'Hill_five')
        ev <- 1 - (1 + (p[3] - 1) / (1 + (p[1] / conc[j])^p[2])) * (1 - 1 / (1 + (p[4] / conc[j])^p[5]))
      effv[i, j] <- ev
    }
  }
  
  colName <- paste('Rspn_@_', conc, sep = '')
  colnames(effv) <- colName
  if(is.null(rownames(param))) rownames(effv) <- model else rownames(effv) <- rownames(param)
  return(effv)
}

mixtoxPlot <- function(x, y, ca, ia, effPoints, ...){
  plot(x, as.matrix(y), xlim = range(x, ca, ia),
       xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)",
       log = "x", pch = 19, ...)
  lines(ca, effPoints, lwd = 1.5, col = 2)
  lines(ia, effPoints, lwd = 1.5, col = 3)
}