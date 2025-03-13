remove(list = ls())
screen_methods <- c('PMDD-BHq' )#'PMDC-DS', 'REDS', 'Knockoff-PMDC', 'PMDD-BHq','PMDC-MDS'
n <- c(100, 200)
s0 <- c(50, 100)
for (i in seq_along(n)) {
  for (j in seq_along(s0)) {
    n_val <- n[i]
    s0_val <- s0[j]
    final_result <- data.frame(
      ps = rep(NaN, length(screen_methods)),
      psm = rep(NaN, length(screen_methods)),
      pa = rep(NaN, length(screen_methods)),
      fdr = rep(NaN, length(screen_methods)),
      tpr = rep(NaN, length(screen_methods)),
      ms = rep(NaN, length(screen_methods))
    )
    rownames(final_result) <- screen_methods
    
    for (k in seq_along(screen_methods)) {
      method <- screen_methods[k]
      file_path <- paste("tp", method, "linear", n_val, s0_val, "results.RData", sep = "_")
      load(file = file_path)
      # cat("Summary Results for ", method, " (n =", n_val, ", s0 =", s0_val, ")\n")
      final_result[method, "fdr"] <- median(results_list[[k]]$fdr)
      final_result[method, "tpr"] <- mean(results_list[[k]]$tpr)
      final_result[method, "pa"] <- mean(results_list[[k]]$pa)
      final_result[method, "ms"] <- mean(results_list[[k]]$ms)
    }
    
    final_result_table <- final_result[, c("fdr", "tpr", "pa", "ms")]
    cat("\n")
    cat("Results for n =", n_val, ", s0 =", s0_val, ":\n")
    print(final_result_table)
    cat("\n\n")
  }
}
