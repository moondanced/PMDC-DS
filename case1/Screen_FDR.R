library("mvtnorm")
library('MASS')
library("knockoff")
start_time <- Sys.time()


PMDC_DS <- function(y, x, alpha, qreg, tau){
    source("mirrorstat_pmdc.R")
    source("DSscreen.R")
    n = nrow(x)
    p = ncol(x)

    # DS Method
    dsvalue = DSscreen(y, x, qreg, tau)
    t = sort(dsvalue, decreasing = FALSE, index.return = TRUE)$x
    t0 = t[t > 0]
    fdp = lapply(1:length(t0), function(k)
      (sum(dsvalue <= -t0[k]) + 1) / max(sum(dsvalue >= t0[k]), 1))
    
    cutoff_set <- max(abs(dsvalue))
    for(t in abs(dsvalue)){
      ps <- length(dsvalue [dsvalue > t])
      ng <- length(na.omit(dsvalue[dsvalue < -t]))
      rto <- (ng + 1)/max(ps, 1)
      if(rto <= alpha){
        cutoff_set <- c(cutoff_set, t)
      }
    }
    T.alpha<- min(cutoff_set) 
    
    ms = length(which(dsvalue >= T.alpha))
    pos= which(dsvalue >= rep(T.alpha, p))
    fdr = (sum((dsvalue <= -T.alpha))) / max(sum(dsvalue >= T.alpha), 1)
    
  return(list(
    "pos" = pos,
    "fdr" = fdr,
    "ms" = ms
  ))
}

PMDC_MDS <- function(y, x, alpha, qreg, tau){
  source("mirrorstat_pmdc.R")
  source("DSscreen.R")
  source("MDSscreen.R")
  n = nrow(x)
  p = ncol(x)
  msplit = 5
  
  # DS Method
  mdsvalue = MDSscreen(y,x,msplit,alpha, qreg, tau)
  IRvalue=sort(mdsvalue,decreasing=F, index.return=TRUE)$x 
  index.L=max(which(cumsum(IRvalue)<=alpha)) 
  T.alpha=IRvalue[index.L]
  
  ms = length(which(mdsvalue >= T.alpha))
  pos= which(mdsvalue >= rep(T.alpha, p))
  fdr = T.alpha
  
  return(list(
    "pos" = pos,
    "fdr" = fdr,
    "ms" = ms
  ))
}


REDS <- function(y, x, alpha, qreg, tau){
    source("REDS_pmdc.R")
    source("REDScreen.R")
    n = nrow(x)
    p = ncol(x)

    # DS Method
    dsvalue = DSscreen(y, x, qreg, tau)
    t = sort(dsvalue, decreasing = FALSE, index.return = TRUE)$x
    t0 = t[t > 0]
    fdp = lapply(1:length(t0), function(k)
      (sum(dsvalue <= -t0[k]) + 1) / max(sum(dsvalue >= t0[k]), 1))
    
    cutoff_set <- max(abs(dsvalue))
    for(t in abs(dsvalue)){
      ps <- length(dsvalue [dsvalue > t])
      ng <- length(na.omit(dsvalue[dsvalue < -t]))
      rto <- (ng + 1)/max(ps, 1)
      if(rto <= alpha){
        cutoff_set <- c(cutoff_set, t)
      }
    }
    T.alpha <- min(cutoff_set)
    
    ms = length(which(dsvalue >= T.alpha))
    pos= which(dsvalue >= rep(T.alpha, p))
    fdr = (sum((dsvalue <= -T.alpha))) / max(sum(dsvalue >= T.alpha), 1)
    

  return(list(
    "pos" = pos,
    "fdr" = fdr,
    "ms" = ms
  ))
}

Knockoff_PMDC <- function(y, x, alpha, qreg, tau){
    source("knockoff_pmdc.R")

    n = nrow(x)
    p = ncol(x)

    # knockoff Method
    pcvalue = lapply(1:p,function(j) pmdc(y,matrix(x[,j],n,1) ) ) 
    pcvalue = unlist(pcvalue)
    xtilde=knockoff::create.second_order(x,method = 'equi')
    #xtilde = knockoff::create.fixed(x,method = 'sdp',randomize = T)$Xk

    pctilde = lapply(1:p,function(j) pmdc(y,matrix(xtilde[,j],n,1) )  ) 
    pctilde = unlist(pctilde)
    dsvalue = pcvalue - pctilde

    t = sort(dsvalue, decreasing = FALSE, index.return = TRUE)$x
    t0 = t[t > 0]
    fdp = lapply(1:length(t0), function(k)
      (sum(dsvalue <= -t0[k]) + 1) / max(sum(dsvalue >= t0[k]), 1))
    
    cutoff_set <- max(abs(dsvalue))
    for(t in abs(dsvalue)){
      ps <- length(dsvalue [dsvalue > t])
      ng <- length(na.omit(dsvalue[dsvalue < -t]))
      rto <- (ng + 1)/max(ps, 1)
      if(rto <= alpha){
        cutoff_set <- c(cutoff_set, t)
      }
    }
    T.alpha<- min(cutoff_set)

    ms = length(which(dsvalue >= T.alpha))
    pos = which(dsvalue >= rep(T.alpha, p))
    fdr = (sum((dsvalue <= -T.alpha))) / max(sum(dsvalue >= T.alpha), 1)

  return(list(
    "pos" = pos,
    "fdr" = fdr,
    "ms" = ms
  ))
}


PMDC_BHq <- function(y, x, alpha, qreg, tau){
  source("pmdd_H.R")
  source("DSBHq.R")
  n = nrow(x)
  p = ncol(x)
  
  dsvalue <- DSscreen(y, x, qreg, tau)
  M <- dsvalue[[1]]  # Extract M values
  S <- dsvalue[[2]]  # Extract S values
  
  T1 = rep(0,p)
  T1 <- sqrt(n/2) * M / sqrt(S)
  p_values <- 2 * (1 - pnorm(abs(T1)))
  
  reject_H0 <- bh_fdr_control(p_values, alpha)
  TSM <- read.csv("data/phs000276/TSM.csv")
  X_true <- na.omit(unique(TSM[drug]))
  ms = sum(reject_H0)
  fdr = max(p_values[pos])*p/ms #
  pos= which(reject_H0==TRUE)
  
  return(list(
    "pos" = pos,
    "fdr" = fdr,
    "ms" = ms
  ))
}


Screen_process <- function(y, x, alpha = 0.2, method = 'PMDC-DS', qreg=FALSE, tau=0.5){
    if(method == "REDS"){
        result = REDS(y, x, alpha, qreg, tau)
    }else if(method == "PMDC-DS"){
        result = PMDC_DS(y, x, alpha, qreg, tau)
    }else if(method == "Knockoff"){
        result = Knockoff_PMDC(y, x, alpha, qreg, tau)
    }else if(method == "PMDC-MDS"){
        result = PMDC_MDS(y, x, alpha, qreg, tau)
    }else if(method == "PMDD-BHq"){
        result = PMDC_BHq(y, x, alpha, qreg, tau)
    }
    return(result)
}