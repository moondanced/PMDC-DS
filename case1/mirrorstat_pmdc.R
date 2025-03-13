mirrorstat=function(y1,Z1, y2, Z2){

n1=nrow(Z1); p=ncol(Z1);
n2=nrow(Z2); n=n1+n2


#a1= diag( 1/sqrt(1+ diag(Z1%*%t(Z1))) )
a1= diag( c(1/sqrt(1+ Z1^2)) )
a2= diag( c(1/sqrt(1+ Z2^2)) )
a=a1%*%(1+Z1%*%t(Z2))%*%a2

a[a>1]=1
a[a<(-1)]=-1
A=asin(a)

M=matrix(y1-mean(y1),1,n1)%*%A%*%matrix(y2-mean(y2),n2,1)/((n1*n2))/(2*pi)

y=c(y1,y2)
B=var(y)
Z=c(Z1,Z2)
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

