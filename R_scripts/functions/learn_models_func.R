

# Multivariate Cox models -------------------------------------------------

# @description: learn models for elastic net, lasso, ridge and adaptive elastic net
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - y_cox: Surv object containing follow-up times and status
#   - method: "EN" for elastic net, "lasso" for lasso, "ridge" for ridge and "AEN"
#             for adaptive elastic net
#   - n_rep: number of models to learn
#
# @return:
#   - list of length 'n_rep' containing 'cv.glmnet' object
learn_models <- function(gene_data, clin, method){
  
  y_cox <- Surv(clin$time, clin$status)
  
  # Elastic Net
  if(method == "EN"){
    
    fit <- cv.glmnet(gene_data, y_cox, family = "cox",
                                alpha = 0.3, nfolds = 5, grouped = T, standardize = F)
    
  # Random Forest
  }else if(method == "RF"){
    
    data_ranger <- cbind(time = clin$time, status = clin$status, gene_data)
    
    colnames(data_ranger) <- str_replace_all(colnames(data_ranger), "-", "_")
    colnames(data_ranger) <- str_replace_all(colnames(data_ranger), "\\|", "_")
    colnames(data_ranger) <- str_replace_all(colnames(data_ranger), "\\?", "a")
    
    fit <- tuneRanger::tuneMtryFast(
      formula = NULL,
      data = data_ranger,
      dependent.variable.name = "time",
      status.variable.name = "status",
      num.treesTry = 50,
      doBest = TRUE,
      plot = F,
      trace = F) 
      
  }else if(method == "coxph"){
    
    clin_gene_data <- cbind(y_cox, gene_data)
    fit <- coxph(y_cox~., data = clin_gene_data)
  }
  
  print("Model learned")
  return(fit)
}

# @description: return a vector containing the names of the genes selected
#
# @parameters:
#   - fit: 'cv.glmnet' object
#   - lambda: weight of the penalty
#
# @return:
#   - vector containing the names of the genes selected
genes_selected_func <- function(fit, lambda){
  
  if(class(fit) == "cv.glmnet"){
    # coefficients and active coefficients
    coefs <- coef(fit, s = lambda)
    names_active_coefs <- row.names(coefs)[coefs@i+1]
    return(names_active_coefs)
  }else{
    return(NULL)
  }
  
  
  
}


# Univariate Cox selection ------------------------------------------------

# @description: univariate Cox for one feature
#
# @parameters:
#   - x: feature (one gene)
#   - y_cox: Surv object containing follow-up times and status
#
# @return:
#   - the p-value of the univariate Cox model
p_val_cox_func <- function(x, y_cox){
  
  res <- try(fit <- coxph(y_cox ~ x))
  
  if(!inherits(res, "try-error")){
    
    test <- summary(fit)
    return(test$logtest[3])
    
  }else{
    return(NA)
  }
  
}


# @description: univariate Cox for an entire matrix
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - y_cox: Surv object containing follow-up times and status
#
# @return:
#   - the p-value of the univariate Cox model
p_val_univCox_func <- function(gene_data, y_cox){
  
  # compute the p-values
  p_val_univCox_vect <- apply(gene_data, 2, function(x) p_val_cox_func(x, y_cox)) 
  
  print("Model learned")
  
  return(p_val_univCox_vect)
}

# @description: compute the p-value of the log-rank test with a 'survdiff' 
#               object as input
#
# @parameters:
#   - x: 'survdiff' object
#
# @return:
#   - the p-value of logrank test
pvalue.survdiff <- function(x, log_p=FALSE, ...){
  pchisq(x$chisq, length(x$n) - 1, lower.tail=FALSE, log.p=log_p)
}






# Elastic net -------------------------------------------------------------

# @description: learn models for elastic net, lasso, ridge and adaptive elastic net
#
# @parameters:
#   - gene_data: matrix of RNA-seq data (column: gene, row: patients)
#   - clin: data frame containing follow-up times and status
#   - alpha: weight in the elastic net penalizatinn between l1 and l2 norms
#   
# @return:
#   - 'cv.glmnet' object
learn_model_EN <- function(gene_data, clin, alpha){
  
  y_cox <- Surv(clin$time, clin$status)
  
  fit <- cv.glmnet(gene_data, y_cox, family = "cox",
                   alpha = alpha, nfolds = 5, grouped = T, standardize = F)
  
  print("Model learned")
  return(fit)
}