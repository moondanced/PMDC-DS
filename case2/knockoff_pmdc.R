library("knockoff")
pmdc=function(y1,Z1){

n1=nrow(Z1); p=ncol(Z1);
n=2*n1


a1= diag( c(1/sqrt(1+ Z1^2)) )
a=a1%*%(1+Z1%*%t(Z1))%*%a1

a[a>1]=1
a[a<(-1)]=-1
A=asin(a)

M=matrix(y1-mean(y1),1,n1)%*%A%*%matrix(y1-mean(y1),n1,1)/((n1*n1))/(2*pi)


B=var(y1)
Z=c(Z1)
a1= diag( c(1/sqrt(1+ Z^2)) )
a=a1%*%(1+Z%*%t(Z))%*%a1
a[a>1]=1
a[a<(-1)]=-1
A=asin(a)
C=1/4-sum(A)/(2*pi)/n^2

M_tilde=M/B/C
#print(M)

return(M_tilde)  
    
}                                          


mirrorstat <- function(y,x1){
  
  gamma = 1;
  n1=nrow(x1); p=ncol(x1);
    
  z1 = pmdc(y,x1)
  x2 =knockoff::create.fixed(x1,method = 'sdp',randomize = T)
  z2 = pmdc(y,x2)

  Mirror = z1-z2
  #Mirror = sign(w1-w2)*max(w1,w2)
  return(Mirror)
}