#################MDS##################
######################################
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
#x<-GSE5680[-which(rownames(GSE5680)=='1389163_at'),]
x<-x_18976
#x=t(x)
x=(scale(x))
y<-GSE5680['1389163_at',]
#y=(scale(y))
y=sapply(y, as.double)
dim(y)<-c(length(y),1)
p = ncol(x);
n = nrow(x);
remove(x_18976,GSE5680)

tau=0.70
y=tau-(y<=quantile(y,tau))

mdcvalue = lapply(1:p, function(j) mdc(matrix(x[, j], n, 1), y, center = "D"))
mdcvalue = unlist(mdcvalue)
t1 = sort(mdcvalue, decreasing = TRUE)[d]
idx_mdc = sort(mdcvalue, decreasing = TRUE, index.return = TRUE)$ix[1:d]


dsvalue=lapply(1:p,function(j) pmdc(y, matrix(x[, j], n, 1))  ) 
dsvalue=unlist(dsvalue)
t2 = sort(dsvalue, decreasing = TRUE)[d]
idx_pmdc = sort(dsvalue, decreasing = TRUE, index.return = TRUE)$ix[1:d]

idx_true = intersect(idx_mdc, idx_pmdc)
save(idx_true, file="GSE5680_070.RData")
