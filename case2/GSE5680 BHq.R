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



############################## DS  ##############################
############################## DS  ##############################
qreg = TRUE
tau=0.75
#y=tau-(y<=quantile(y,tau))
source("pmdd_H.R")
source("DSBHq.R")
# BHq Method
dsvalue = DSscreen(y, x, qreg, tau)
M <- dsvalue[[1]]  # Extract M values
S <- dsvalue[[2]]  # Extract S values

# Compute test statistics
T1 = rep(0,p)
T1 <- sqrt(n/2) * M / sqrt(S)

# Compute two-sided p-values from normal distribution
p_values <- 2 * (1 - pnorm(abs(T1)))

# Apply Benjamini-Hochberg FDR correction
reject_H0 <- bh_fdr_control(p_values, alpha)

pmdc.pos= which(reject_H0)
load("GSE5680_075.RData")
pmdc.tpr= length(intersect(pmdc.pos, idx_true))/length(idx_true)
pmdc.fdr = max(p_values[pmdc.pos])*p/max(sum(reject_H0), 1)

results <- data.frame(
  Method = c("PMDD-BHq"),
  FDR = c(pmdc.fdr),
  TPR = c(pmdc.tpr)
)

# Print the dataframe
print(results)
