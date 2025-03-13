
source("mirrorstat_pmdc.R")
source("DSscreen.R")
source("MDSscreen.R")
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
  
  pis = psm = matrix(0, Nsim, p)
  fdr = tpr = pa = ms = rep(0, Nsim)
  
  
  for (i in 1:Nsim) {
    set.seed(i)
    data = data_generate(models, para_list)
    x = data$x
    y = data$y
    index.true = data$index.true
    beta = rep(0, p)
    beta[index.true] = 1
    
    mdsvalue=MDSscreen(y,x,msplit,alpha, qreg, tau)
    IRvalue=sort(mdsvalue,decreasing=F, index.return=TRUE)$x 
    index.L=max(which(cumsum(IRvalue)<=alpha)) 
    T.alpha=IRvalue[index.L]
    
    plot(mdsvalue, main = 'mdsvalue')
    abline(h = T.alpha, col = "red")
    plot(mdsvalue[1:150], main = 'mdsvalue')
    abline(h = T.alpha, col = "red")
    
    print(which(mdsvalue >= T.alpha))
    ms[i] = length(which(mdsvalue >= T.alpha))
    pis[i, ] = (mdsvalue >= rep(T.alpha, p))
    psm[i, ] = apply(matrix(pis[1:i, ], i), 2, mean)
    pa[i] = (min(1 * (mdsvalue[index.true] >= T.alpha)) >= 1)
    fdr[i] = (sum((mdsvalue >= T.alpha) * (beta == 0))) / max(sum(mdsvalue >= T.alpha), 1)
    tpr[i] = (sum((mdsvalue >= T.alpha) * (beta != 0))) / max(sum(beta != 0), 1)
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