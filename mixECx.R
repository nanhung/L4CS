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