

# https://stackoverflow.com/questions/45104987/calculating-brier-score-and-integrated-brier-score-using-ranger-r-package
# prediction function for computing the IBS with pec
predictSurvProb.ranger <- function (object, newdata, times, ...) {
  ptemp <- ranger:::predict.ranger(object, data = newdata, importance = "none")$survival
  pos <- prodlim::sindex(jump.times = object$unique.death.times, 
                         eval.times = times)
  p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  p
}

# function to compute the C-index and the IBS
C_IBS_func <- function(fit, clin_train, clin_test, data_train, 
                            data_test, method){
  
  # test if there is at least 1 event
  if(sum(clin_test$status == 1) != 0){
    
    y_cox_test <- Surv(clin_test$time, clin_test$status)
    y_cox_train <- Surv(clin_train$time, clin_train$status)
    
  
    if(method == "EN"){
      
      if(length(genes_selected_func(fit, "lambda.min")) >= 1){
        
        # prognostic index
        beta <- as.numeric(coef(fit, "lambda.min"))
        names(beta) <- colnames(data_train) 
        PI_test <-  data_test %*% beta
        PI_train <- data_train %*% beta
        
        # C-index
        C <- concordance.index(PI_test, surv.time = clin_test$time, 
                               surv.event = clin_test$status)$c.index 
        
        # re-calibrate the model 
        df_coxph <- data.frame(time = clin_train$time, status = clin_train$status, PI = PI_train)
        fit_coxph <- coxph(Surv(time, status) ~ PI, data = df_coxph,  x = TRUE)
        PI_train <- PI_train * fit_coxph$coefficients
        PI_test <- PI_test * fit_coxph$coefficients
        
        # IBS
        data_pec <- data.frame(time = clin_test$time, status = clin_test$status, PI = PI_test)
        res <- try(PredError <- pec(object=fit_coxph,
                         formula = as.formula("Surv(time, status)~PI"), 
                         cens.model="cox",
                         data=data_pec, verbose=F))
        
        if(!inherits(res, "try-error")){
          IBS <- crps(PredError, start = 0, times = max(clin_test$time[clin_test$status==1]))[2]
        }else{
          IBS <- NA
        }
        
        return(c(C=C, IBS=IBS))
        
      }else{
        return(c(C = NA, IBS = NA))
      }
      
    }else if(method == "RF"){
      
      data_ranger_train <- cbind(time = clin_train$time, status = clin_train$status, data_train)
      data_ranger_test <- cbind(time = clin_test$time, status = clin_test$status, data_test)
      
      colnames(data_ranger_train) <- str_replace_all(colnames(data_ranger_train), "-", "_")
      colnames(data_ranger_train) <- str_replace_all(colnames(data_ranger_train), "\\|", "_")
      colnames(data_ranger_train) <- str_replace_all(colnames(data_ranger_train), "\\?", "a")
      
      colnames(data_ranger_test) <- str_replace_all(colnames(data_ranger_test), "-", "_")
      colnames(data_ranger_test) <- str_replace_all(colnames(data_ranger_test), "\\|", "_")
      colnames(data_ranger_test) <- str_replace_all(colnames(data_ranger_test), "\\?", "a")
      
      # risck score (mean of the cumulative hazard function) for the C-index ---
      PI_test <- ranger:::predict.ranger(fit, data = data_ranger_test)
      PI_test <- rowMeans(PI_test$chf)
      
      C <- concordance.index(PI_test, surv.time = clin_test$time, 
                             surv.event = clin_test$status)$c.index
      
      # IBS with pec ---
      
      # A formula to be inputted into the pec command 
      # https://stackoverflow.com/questions/45104987/calculating-brier-score-and-integrated-brier-score-using-ranger-r-package
      frm <- as.formula(paste("Surv(time, status)~",
                              paste(fit$forest$independent.variable.names, collapse="+")))
      
      # Using pec for IBS estimation
      PredError <- pec(object=fit,
                       formula = frm, cens.model="marginal",
                       data=as.data.frame(data_ranger_test), verbose=F)
      
      IBS <- crps(PredError, start = 0, times = max(clin_test$time[clin_test$status==1]))[2]
      
      return(c(C = C, IBS = IBS))
      
    }else if(method == "coxph"){
      
      PI_test <- as.vector(predict(fit, newdata = data_test, type = "lp"))
      PI_train <- as.vector(predict(fit, newdata = data_train, type = "lp"))
      
      df_coxph <- data.frame(time = clin_train$time, status = clin_train$status, PI = PI_train)
      
      # re-calibrate the model 
      fit_coxph <- coxph(Surv(time, status) ~ PI, data = df_coxph,  x = TRUE)
      PI_train <- PI_train * fit_coxph$coefficients
      PI_test <- PI_test * fit_coxph$coefficients
      
      C <- concordance.index(PI_test, surv.time = clin_test$time, 
                             surv.event = clin_test$status)$c.index 
      
      # IBS
      data_pec <- data.frame(time = clin_test$time, status = clin_test$status, PI = PI_test)
      PredError <- pec(object=fit_coxph,
                       formula = as.formula("Surv(time, status)~PI"), 
                       cens.model="marginal",
                       data=data_pec, verbose=F)
      
      IBS <- crps(PredError, times = max(clin_test$time[clin_test$status==1]))[2]
      
      return(c(C=C, IBS=IBS))
      
      return(c(C=C, IBS=IBS))
    }
    
  }else{
    return(c(C=NA, IBS=NA))
  }
}
