DSscreen<-function(y,x,qreg = FALSE,tau=0.5){
  
  K=3
  n=nrow(x); p=ncol(x);
  
  #train_num<-sample(1:n,n*0.5)  
  train_num<-sample(1:n,(K-1)*n/K)  
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

  Mvalue=lapply(1:p,function(j) mirrorstat(y.spit1,matrix(x.spit1[,j],n1,1),
                                           y.spit2,matrix(x.spit2[,j],n2,1) )  ) 
  
  Mvalue=unlist(Mvalue)
  
  return(Mvalue)
}