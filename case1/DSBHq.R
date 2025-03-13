DSscreen=function(y,x,qreg=FALSE,tau=0.5){
 

n=nrow(x); p=ncol(x);

train_num<-sample(1:n,n*0.5)  
y.spit1=y[train_num,] 
y.spit2=y[-train_num,] 

x.spit1=x[train_num,] 
n1=nrow(x.spit1)
x.spit2=x[-train_num,] 
n2=nrow(x.spit2)
 
if(qreg){
  y.spit1=tau-(y.spit1<=quantile(y.spit1,tau))
  y.spit2=tau-(y.spit2<=quantile(y.spit2,tau))
}

Hvalue=lapply(1:p,function(j) mirrorstat(y.spit1,matrix(x.spit1[,j],n1,1),
                  y.spit2,matrix(x.spit2[,j],n2,1) )  ) 

Mvalue = lapply(1:p,function(j) sum(Hvalue[[j]][[1]])/((n1*n2))/(2*pi)  ) 
Mvalue=unlist(Mvalue)

Svalue = lapply(1:p,function(j) var(colMeans(Hvalue[[j]][[2]]/(2*pi)))  )
Svalue=unlist(Svalue)

return(list(Mvalue,Svalue))
}                                          


bh_fdr_control <- function(p_values, alpha) {
  p <- length(p_values)
  sorted_indices <- order(p_values)  
  sorted_p <- p_values[sorted_indices]
  
  # Compute BH threshold
  threshold <- (1:p) / p * alpha #/ sum(1/(1:p))
  significant <- sorted_p <= threshold
  
  # Determine the largest k where condition holds
  if (any(significant)) {
    k <- max(which(significant))
    reject <- rep(FALSE, p)
    reject[sorted_indices[1:k]] <- TRUE  # Mark rejected hypotheses
  } else {
    reject <- rep(FALSE, p)
  }
  
  return(reject)
}