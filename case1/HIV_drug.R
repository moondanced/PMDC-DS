remove(list = ls())

for (drug in c("PI", "NRTI")) {
  for (method in c("PMDD-BHq")) {
    TSM <- read.csv("data/phs000276/TSM.csv")
    
    alpha <- 0.1 #0.05
    qreg <- FALSE
    
    X_true <- na.omit(unique(TSM[[drug]]))
    X_f <- c()
    result_table <- data.frame(Response_Variable = character(),
                               Feature_Length = numeric(),
                               Result = character(),
                               stringsAsFactors = FALSE)
    
    if (drug == "PI") {
      data <- read.csv("data/phs000276/PI_DATA.csv")
      response_variables <- c("IDV", "LPV", "NFV")
    } else if (drug == "NRTI") {
      data <- read.csv("data/phs000276/NRTI_DATA.csv")
      response_variables <- c("ABC", "AZT", "D4T")
    }
    
    for (response in response_variables) {
      data_f <- data[complete.cases(data[[response]]), ]
       Y <- log(data_f[[response]])
      Y <- matrix(Y, length(Y), 1)

      
      if (drug == "PI") {
        X <- data_f[, 11:109]
      } else if (drug == "NRTI") {
        X <- data_f[, 10:249]
      }
      
      for (n in colnames(X)) {
        X[[n]][X[[n]] == "-"] <- 0
        X[[n]][X[[n]] == "."] <- 0  
        keep <- as.numeric(ave(X[[n]], X[[n]], FUN = length)) >= 3
        
        X[[n]] <- ifelse(keep, X[[n]], 0)
        g_add <- unique(X[[n]][X[[n]] != 0])
        for (g in g_add) {
          X_f <- c(X_f, as.numeric(X[[n]] == g))
        }
      }
      X_f <- c(X_f, 1*rbinom(800*nrow(data_f),1,0.5))
      X_f <- matrix(X_f, nrow = nrow(data_f))
      source("Screen_FDR.R")
      Result <- Screen_process(Y, X_f, alpha, method, qreg, tau)
      
      result_table <- rbind(result_table, data.frame(Response_Variable = response,
                                                     Feature_Length = length(Result),
                                                     Result = Result))
      
      X_f <- c()
    }
    
    if (drug == "PI") {
      save(result_table, file = paste0(drug, "_", method, "_results.RData"))
    } else if (drug == "NRTI") {
      save(result_table, file = paste0(drug, "_", method, "_results.RData"))
    }
    
  }
}
