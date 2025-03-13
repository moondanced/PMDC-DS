# Remove any existing objects in the workspace
remove(list = ls())

# Load necessary libraries
library("mvtnorm")
library('MASS')
library("EDMeasure")
library("Matrix")

# Source the data generation script
source("generate_data.R")

# Record the start time of the simulation
start_time <- Sys.time()

# Function to perform the simulation for given n and s0
run_simulation <- function(n, s0) {
  # Parameter definitions
  p = 5000
  d = 30
  n1 = 40
  n2 = n - n1
  N = n * p
  alpha = 0.2
  msplit = 10
  
  Sig0 = toeplitz(seq(1,0,-1/(p/100-1)))*0.6
  diag(Sig0) = rep(1,p/100)
  S_list = list()
  for (i in 1:100) {
    S_list = c(S_list, list(Sig0))
  }
  Sig = (bdiag(S_list))
  
  beta = 1 * c(rep(1, s0), rep(0, p - s0))
  index.true = (1:p)[beta != 0]
  C = chol(Sig)
  
  # Model settings
  para_list = list(n = n, p = p, s0 = s0, d = d, n1 = n1, n2 = n2, N = N, Sig = Sig, msplit = msplit, C = C, alpha = alpha)
  models = "hd-additive"
  Nsim = 200
  d1 = floor(n / log(n))
  d2 = floor(2 * n / log(n))
  
  d1.pis = d1.psm = matrix(0, Nsim, p)
  d1.fdr = d1.tpr = d1.pa = d1.ms = rep(0, Nsim)
  d2.pis = d2.psm = matrix(0, Nsim, p)
  d2.fdr = d2.tpr = d2.pa = d2.ms = rep(0, Nsim)
  
  for (i in 1:Nsim) {
    set.seed(i)
    data = data_generate(models, para_list)
    x = data$x
    y = data$y
    index.true = data$index.true
    beta = rep(0, p)
    beta[index.true] = 1
    
    # DS Method
    dsvalue = lapply(1:p, function(j) mdc(matrix(x[, j], n, 1), y, center = "D"))
    dsvalue = unlist(dsvalue)
    t1 = sort(dsvalue, decreasing = TRUE)[d1]
    
    d1.ms[i] = length(which(dsvalue >= t1))
    d1.pis[i, ] = (dsvalue >= rep(t1, p))
    d1.psm[i, ] = apply(matrix(d1.pis[1:i, ], i), 2, mean)
    d1.pa[i] = (min(1 * (dsvalue[index.true] >= t1)) >= 1)
    d1.fdr[i] = (sum((dsvalue >= t1) * (beta == 0))) / max(sum(dsvalue >= t1), 1)
    d1.tpr[i] = (sum((dsvalue >= t1) * (beta != 0))) / max(sum(beta != 0), 1)
    
    t2 = sort(dsvalue, decreasing = TRUE)[d2]
    
    d2.ms[i] = length(which(dsvalue >= t2))
    d2.pis[i, ] = (dsvalue >= rep(t2, p))
    d2.psm[i, ] = apply(matrix(d2.pis[1:i, ], i), 2, mean)
    d2.pa[i] = (min(1 * (dsvalue[index.true] >= t2)) >= 1)
    d2.fdr[i] = (sum((dsvalue >= t2) * (beta == 0))) / max(sum(dsvalue >= t2), 1)
    d2.tpr[i] = (sum((dsvalue >= t2) * (beta != 0))) / max(sum(beta != 0), 1)
    
    cat("Simulation", i, "for alpha =", alpha, "\n")
  }
  
  results_list = list(
    d1.pis = d1.pis, d1.psm = d1.psm, d1.pa = d1.pa, d1.fdr = d1.fdr, d1.tpr = d1.tpr, d1.ms = d1.ms,
    d2.pis = d2.pis, d2.psm = d2.psm, d2.pa = d2.pa, d2.fdr = d2.fdr, d2.tpr = d2.tpr, d2.ms = d2.ms
  )
  
  return(results_list)
}

# Specify values for n and s0
n_values <- c(100, 200)
s0_values <- c(50, 100)

# Initialize an empty list to store results
all_results <- list()

# Run simulation for each combination of n and s0
for (n_val in n_values) {
  for (s0_val in s0_values) {
    results <- run_simulation(n_val, s0_val)
    all_results[[paste("n", n_val, "s0", s0_val)]] <- results
  }
}

save(all_results, file = "tp_hd-additive_MDC_results.RData")

# Combine results into a data frame
results_df1 <- data.frame(
  #n = rep(n_values, each = length(s0_values) * Nsim),
  #s0 = rep(rep(s0_values, each = Nsim), length(n_values)),
  fdr = unlist(lapply(all_results, function(res) mean(res$d1.fdr))),
  tpr = unlist(lapply(all_results, function(res) mean(res$d1.tpr))),
  pa = unlist(lapply(all_results, function(res) mean(res$d1.pa))),
  ms = unlist(lapply(all_results, function(res) mean(res$d1.ms)))
)

# Print the results data frame
print(results_df1)

# Combine results into a data frame
results_df2 <- data.frame(
  #n = rep(n_values, each = length(s0_values) * Nsim),
  #s0 = rep(rep(s0_values, each = Nsim), length(n_values)),
  fdr = unlist(lapply(all_results, function(res) mean(res$d2.fdr))),
  tpr = unlist(lapply(all_results, function(res) mean(res$d2.tpr))),
  pa = unlist(lapply(all_results, function(res) mean(res$d2.pa))),
  ms = unlist(lapply(all_results, function(res) mean(res$d2.ms)))
)

# Print the results data frame
print(results_df2)
