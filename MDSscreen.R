MDSscreen=function(y,x,msplit,alpha,qreg=FALSE,tau=0.5){

n=nrow(x); p=ncol(x);
IR.mat=matrix(0, msplit, p)  #inclusion rate

for(i in 1:msplit){

Mirror=DSscreen(y,x,qreg,tau)

t=sort(Mirror,decreasing=F, index.return=TRUE)$x
t0=t[t>0]
fdp=lapply(1:length(t0),function(k)
           (sum(Mirror<= -t0[k])+1)/max( sum( Mirror>= t0[k]),1) ) 
       
ind.alpha=which(unlist(fdp)<=alpha)[1]

if( is.na(ind.alpha) ){ T.alpha=10^8;     ##if empty,T.alpha set +infinite 
}else{ T.alpha=t0[ind.alpha] } 
                                         
IR.mat[i,]=(Mirror>=T.alpha)/max(sum(Mirror>=T.alpha),1) 

}

IR=apply(IR.mat,2,mean)

return(IR)
}
 