remove(list = ls())
load("data/GSE5680.Rdata")
load("data/18976_x.Rdata")
library("mvtnorm")
library('MASS')
library("knockoff")
library("EDMeasure")
library("Matrix")

source("mirrorstat_pmdc.R")
source("DSscreen.R")
source("MDSscreen.R")
set.seed(123)
x<-x_18976
x=(scale(x))
y<-GSE5680['1389163_at',]
y=sapply(y, as.double)
dim(y)<-c(length(y),1)
p = ncol(x);
n = nrow(x);
alpha = 0.1;
remove(x_18976,GSE5680)


############################## MDS  ##############################
############################## MDS  ##############################
qreg = TRUE
tau=0.25
msplit=5
#y=tau-(y<=quantile(y,tau))

source("mirrorstat_pmdc.R")
source("DSscreen.R")

# DS Method
mdsvalue = MDSscreen(y, x, msplit, qreg, tau)
IRvalue=sort(mdsvalue,decreasing=F, index.return=TRUE)$x 
index.L=max(which(cumsum(IRvalue)<=alpha)) 
T.alpha=IRvalue[index.L]

load("GSE5680_025.RData")

mds.pos= which(mdsvalue >= rep(T.alpha, p))
mds.fdr = T.alpha
mds.tpr= length(intersect(mds.pos, idx_true))/length(idx_true)

results <- data.frame(
  Method = c("PMDC-MDS"),
  FDR = c(mds.fdr),
  TPR = c(mds.tpr)
)

# Print the dataframe
print(results)
