data_generate <- function(models, para_list){
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
  if(models == "linear"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    x = matrix(rnorm(n * p), n, p) %*% C
    y = x %*% beta + 1 * rnorm(n, 0, 1) 
    index.true = 1:s0
  }
  if(models == "exp-additive"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    x = matrix(rnorm(n * p), n, p) %*% C
    y = matrix(x[, 1] + x[, 2] + x[, 3]^2 + exp((x[, 4] + x[, 5])/2) * rnorm(n,0,1),n,1)
    index.true = 1:s0
  }
  if(models == "poly-additive"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    x = matrix(rnorm(n * p), n, p) %*% C
    y = matrix(5*(x[, 1]+x[, 2]+x[, 3]+x[, 4]) + 5*((x[, 5]+x[, 6])^2+(x[, 7]+x[, 8])^3) + 5*((x[, 9]>0)+(x[, 10]>0)) + rnorm(n, 0, 1),n,1)
    index.true = 1:s0
  }
  if(models == "mixed-norm"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    x1=matrix(rnorm(n*p),n,p)
    x2=rt(n*p,df=2)
    x=(0.9*x1+0.1*x2)%*%C
    y=x%*%beta + rnorm(n)
    index.true = 1:s0
  }
  if(models == "mixed-t2tail"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    x1=matrix(rnorm(n*p),n,p)
    x2=rt(n*p,df=2)
    x=(0.9*x1+0.1*x2)%*%C
    y=x%*%beta + rnorm(n, 0, 1)
    index.true = 1:s0
  }
  if(models == "lin-t2tail"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    x = matrix(rnorm(n * p), n, p) %*% C
    y = x %*% beta + rt(n,df=2)
    index.true = 1:s0
  }
  if(models == "rsignal-lin"){
    beta = 1 * c(-seq(-floor(s0/2),-1), seq(1,ceiling(s0/2)), rep(0, p - s0))/10
    x = matrix(rnorm(n * p), n, p) %*% C
    y = x %*% beta + rnorm(n, 0, 1)
    index.true = 1:s0
  }
  if(models == "hd-additive"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    # beta1 = 1* c(rep(1, 3*s0/5), rep(0, p - 3*s0/5))
    beta1 = 1* c(rep(1, 1*s0/5), rep(0, p - 1*s0/5))
    beta2 = 1* c(rep(0, 1*s0/5), rep(1, 2*s0/5), rep(0, p - 3*s0/5))
    beta3 = 1* c(rep(0, 3*s0/5), rep(1, 2*s0/5), rep(0, p - 5*s0/5))
    # beta4 = 1* c(rep(0, 4*s0/5), rep(1, s0/5), rep(0, p - s0))
    x = matrix(rnorm(n * p), n, p) %*% C
    # print(dim(t(x)*beta4))
    # print(dim(rowSums(t(x) * beta4 > 0)))
    # y = matrix(5*(x %*% beta1) + 5*(x %*% beta2)^2 + 5*(x %*% beta3)^3 + rnorm(n, 0, 1),n,1) #+  5*matrix(rowSums((t(x) * beta4) > 0),n,1)
    y = matrix(5*(x %*% beta1) + 5*(x %*% beta2)^3 + 5*(x %*% beta3)^5 + rnorm(n, 0, 1),n,1) 
    index.true = 1:s0
  }
  if(models == "hd-exp"){
    beta = 1 * c(rep(1, s0), rep(0, p - s0))
    # beta1 = 1* c(rep(1, 10), rep(0, p - 10))
    beta1 = 1* c(rep(1, 40), rep(0, p - 40))
    # beta2 = 1* c(rep(0, 100), rep(1, 10), rep(0, p - 110))
    beta2 = 1* c(rep(0, 40), rep(1, 10), rep(0, p - 50))
    x = matrix(rnorm(n * p), n, p) %*% C
    # y = matrix(x %*% beta1/10 + exp(x %*% beta2/3) *rnorm(n,0,1),n,1)
    # y = matrix(x %*% beta1 + (abs(x %*% beta2) / 16)^(1/5)*rt(n,df=2),n,1)
    y = matrix(x %*% beta1 + (1 + (x %*% beta2) * 2)^(1)*rnorm(n,0,1),n,1)
    # y = matrix((x %*% beta1) + 1.5 * abs(x %*% beta2)*rt(n,df=2),n,1)
    # print((abs(x %*% beta2) / 16)^(1/5))
    index.true = 1:s0
    # index.true = c(1:50, 51:100)
  }
  
  return(list(x = x, y = y, index.true = index.true, beta = beta))
}