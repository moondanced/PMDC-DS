
source("knockoff_pmdc.R")
source("generate_data.R")
start_time <- Sys.time()

Screen_process <- function(models, Nsim, para_list){
  n = para_list$n
  p = para_list$p
  s0 = para_list$s0
  d = para_list$d
  n1 = para_list$n1
  n2 = para_list$n2
  N = para_list$N
  Sig = para_list$Sig
  msplit = para_list$msplit
  C = para_list$C
  alpha = para_list$alpha
  qreg = para_list$qreg
  tau = para_list$tau
  
  pis = psm = matrix(0, Nsim, p)
  fdr = tpr = pa = ms = rep(0, Nsim)
  
  
  for (i in 1:Nsim) {
    set.seed(i)
    if(qreg){set.seed(i+tau*10000)}
    data = data_generate(models, para_list)
    x = data$x
    y = data$y
    index.true = data$index.true
    beta = rep(0, p)
    beta[index.true] = 1
    
    if(qreg){
      y=tau-(y<=quantile(y,tau))
    }
    
    # knockoff Method
    #source("PC_screen.r")
    Aset=1:d
    pcvalue = lapply(1:d,function(j) pmdc(y,matrix(x[,Aset[j]],n,1) ) ) 
    pcvalue = unlist(pcvalue)
    xtilde=knockoff::create.second_order(x,method = 'equi')[,Aset]
    pctilde = lapply(1:d,function(j) pmdc(y,matrix(xtilde[,Aset[j]],n,1) )  ) 
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
    
    ms[i] = length(which(dsvalue >= T.alpha))
    pis[i, ] = (dsvalue >= rep(T.alpha, p))
    psm[i, ] = apply(matrix(pis[1:i, ], i), 2, mean)
    pa[i] = (min(1 * (dsvalue[index.true] >= T.alpha)) >= 1)
    fdr[i] = (sum((dsvalue <= -T.alpha)) / max(sum(dsvalue >= T.alpha), 1))
    tpr[i] = (sum((dsvalue >= T.alpha) * (beta != 0))) / max(sum(beta != 0), 1)
    
    cat("Simulation", i, "for alpha =", alpha, "\n")
  }
  
  return(list(
    "pis" = pis,
    "pa" = pa,
    "psm" = psm,
    "fdr" = fdr,
    "ms" = ms,
    "tpr" = tpr
  ))
}
