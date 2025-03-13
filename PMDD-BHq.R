source("pmdd_H.R")
source("DSBHq.R")
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

  pis = psm = T1 = matrix(0, Nsim, p)
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

    dsvalue <- DSscreen(y, x, qreg, tau)
    M <- dsvalue[[1]]  # Extract M values
    S <- dsvalue[[2]]  # Extract S values

    # Compute test statistics
    T1[i, ] <- sqrt(n1) * M / sqrt(S)
    p_values <- 2 * (1 - pnorm(abs(T1[i, ])))
    # Apply Benjamini-Hochberg FDR correction
    reject_H0 <- bh_fdr_control(p_values, alpha)
    ms[i] = sum(reject_H0)
    pis[i, ] = (reject_H0==TRUE)
    psm[i, ] = apply(matrix(pis[1:i, ], i), 2, mean)
    pa[i] = (min(1 * (reject_H0[index.true] == TRUE)) >= 1)
    fdr[i] = (sum(reject_H0 * (beta == 0))) / max(sum(reject_H0), 1)
    # fdr[i] = (sum((dsvalue <= -T.alpha))) / max(sum(dsvalue >= T.alpha), 1)
    tpr[i] = (sum(reject_H0 * (beta != 0))) / max(sum(beta != 0), 1)
    
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