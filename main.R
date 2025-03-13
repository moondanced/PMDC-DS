remove(list = ls())
library("mvtnorm")
library('MASS')
library("knockoff")
library("Matrix")

# 记录运行时间的示例代码块
start_time <- Sys.time()

#参数定义
for(n in c(100, 200)){
p = 5000
for (s0 in c(50, 100)){
d = 200
n2 = floor(n/2)
n1 = n-n2
N = n * p
alpha = 0.2
msplit = 5
qreg = FALSE
for(tau in c(0.5)){


    
Sig0 = toeplitz(seq(1,0,-1/(p/100-1)))*0.6
diag(Sig0) = rep(1,p/100)
S_list = list()
for (i in 1:100) {
  S_list = c(S_list, list(Sig0))
}
Sig = (bdiag(S_list))


beta = 1 * c(rep(1, s0), rep(0, p - s0))
index.true = (1:p)[beta != 0]
C = chol(Sig)
C = as.matrix(C,ncol=p)

# 模型设置
para_list = list(n=n,p=p,s0=s0,d=d,n1=n1,n2=n2,N=N,Sig=Sig,msplit=msplit,C=C,alpha=alpha,qreg=qreg,tau=tau)
for(models in c("linear")){ # "linear", "mixed-norm", "lin-t2tail", "hd-exp", "exp-additive", "poly-additive", "rsignal-lin"
Nsim = 500
screen_methods = c("PMDC-DS") # "PMDC-DS", "PMDC-MDS", "PMDD-BHq", "REDS", "Knockoff-PMDC"
results_list = list()

for (method in screen_methods) {
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
  
  source(paste(method, ".R", sep=''))
  
  pis = psm = matrix(0, Nsim, p)
  fdr = tpr = pa = ms = rep(0, Nsim)
  
  Screen_result <- Screen_process(models, Nsim, para_list)
  pis = Screen_result$pis
  psm = Screen_result$psm
  pa = Screen_result$pa
  fdr = Screen_result$fdr
  ms = Screen_result$ms
  tpr = Screen_result$tpr
  
  results_list[[method]] <- Screen_result
  if(qreg){
    # save(results_list, file = paste("tp_q", tau, method, models, n , s0, "results.RData", sep = "_"))
  }
  else{
    # save(results_list, file = paste("tp", method, models, n, s0, "results.RData", sep = "_"))
  }
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time

print(paste("run_time:", elapsed_time))

# Save all results
colnum <- length(screen_methods)
final_result <- data.frame(
  ps = rep(NaN, colnum),
  psm = rep(NaN, colnum),  
  pa = rep(NaN, colnum),  
  fdr = rep(NaN, colnum),  
  tpr = rep(NaN, colnum), 
  ms = rep(NaN, colnum) 
)
rownames(final_result) <- screen_methods

# Print summary results
for (i in seq_along(screen_methods)) {
  method <- screen_methods[i]
  cat("Summary Results for ", method, "\n")
  final_result[method, "fdr"] = mean(results_list[[method]]$fdr)
  final_result[method, "tpr"] = mean(results_list[[method]]$tpr)
  final_result[method, "ms"] = mean(results_list[[method]]$ms)
}

print(final_result[screen_methods, c("fdr", "tpr", "ms")])

}
}
}
}