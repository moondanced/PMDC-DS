remove(list = ls())
TSM <- read.csv("data/phs000276/TSM.csv")

sink("drug_output.log")

for (method in c("PMDD-BHq")) {
  for (drug in c('PI','NRTI')) {
    
    cat(paste("The screen method is:", method, "\n"))
    load(paste0(drug, "_", method, "_results.RData"))
    alpha <- 0.1
    qreg <- FALSE
    tau <- 0.5
    
    X_true <- na.omit(unique(TSM[drug]))
    if (drug == 'PI') {
      data <- read.csv("data/phs000276/PI_DATA.csv")
      response_variables <- c("IDV", "LPV", "NFV")
      response_data <- data[, response_variables]
    } else if (drug == 'NRTI') {
      data <- read.csv("data/phs000276/NRTI_DATA.csv")
      response_variables <- c("ABC", "AZT", "D4T")
      response_data <- data[, response_variables]
    }
    
    for (response in response_variables) {
      pos <- result_table$Result.pos[result_table$Response_Variable == response]
      df <- read.csv(paste0("data/phs000276_screen/", response, "_X_f.csv"))
      feature <- colnames(df)
      Aset <- feature[pos]
      G_pos <- as.numeric(gsub("P(\\d+).*", "\\1", Aset))
      G_pos <- unique(G_pos)
      fdr <- result_table$Result.fdr[result_table$Response_Variable == response][1]
      tpr <- length(intersect(G_pos, X_true[,1])) / length(X_true[,1])
      ms <- length(G_pos)
      
      cat("Response Variable:", response, "\n")
      cat("FDR:", fdr, "\n")
      cat("MS:", ms, "\n")
      cat("TPR:", tpr, "\n\n")
    }
    
    
  }
}

sink()

