remove(list = ls())
load("data/GSE5680.Rdata")
load("data/18976_x.Rdata")
library("mvtnorm")
library('MASS')
library("knockoff")
library("EDMeasure")
library("Matrix")

source("REDS_pmdc.R")
source("mirrorstat_pmdc.R")

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

############################## DS  ##############################
############################## DS  ##############################
qreg = TRUE
tau=0.25
#y=tau-(y<=quantile(y,tau))

source("mirrorstat_pmdc.R")
source("DSscreen.R")

# DS Method
dsvalue = DSscreen(y, x, qreg, tau)
t = sort(dsvalue, decreasing = FALSE, index.return = TRUE)$x
t0 = t[t > 0]
fdp = lapply(1:length(t0), function(k)
  (sum(dsvalue <= -t0[k]) + 1) / max(sum(dsvalue >= t0[k]), 1))

cutoff_set <- max(abs(dsvalue))
for(t in abs(dsvalue)){
  ps <- length(dsvalue [dsvalue > t])
  ng <- length(na.omit(dsvalue[dsvalue < -t]))
  rto <- (ng + 1)/max(ps, 1)
  if(rto <= alpha){
    cutoff_set <- c(cutoff_set, t)
  }
}
T.alpha<- min(cutoff_set)
pmdc.pos= which(dsvalue >= rep(T.alpha, p))
pmdc.fdr = (sum((dsvalue <= -T.alpha))) / max(sum(dsvalue >= T.alpha), 1)


# REDS Method
source("REDS_pmdc.R")
source("REDScreen.R")

dsvalue = DSscreen(y, x, qreg, tau)
t = sort(dsvalue, decreasing = FALSE, index.return = TRUE)$x
t0 = t[t > 0]
fdp = lapply(1:length(t0), function(k)
  (sum(dsvalue <= -t0[k]) + 1) / max(sum(dsvalue >= t0[k]), 1))

cutoff_set <- max(abs(dsvalue))
for(t in abs(dsvalue)){
  ps <- length(dsvalue [dsvalue > t])
  ng <- length(na.omit(dsvalue[dsvalue < -t]))
  rto <- (ng + 1)/max(ps, 1)
  if(rto <= alpha){
    cutoff_set <- c(cutoff_set, t)
  }
}
T.alpha <- min(cutoff_set)


reds.pos= which(dsvalue >= rep(T.alpha, p))
reds.fdr = (sum((dsvalue <= -T.alpha))) / max(sum(dsvalue >= T.alpha), 1)

load("GSE5680_070.RData")
pmdc.tpr= length(intersect(pmdc.pos, idx_true))/length(idx_true)
reds.tpr= length(intersect(reds.pos, idx_true))/length(idx_true)

results <- data.frame(
  Method = c("PMDC", "REDS"),
  FDR = c(pmdc.fdr, reds.fdr),
  TPR = c(pmdc.tpr, reds.tpr)
)

# Print the dataframe
print(results)