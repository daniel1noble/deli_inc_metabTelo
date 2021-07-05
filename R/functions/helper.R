#Model prediction function from imputation models
imp_mr_pred <- function(id, temp, samp_session, treatment){
  #Lizard deviation
  id_intercept <- paste0("r_id__lnmr[", id, ",", "Intercept]")
  id_slope <- paste0("r_id__lnmr[", id, ",", "temp]")
  
  #Samp_session deviation
  samp_session_intercept <- paste0("r_samp_session__lnmr[", samp_session, ",", "Intercept]")
  
  #Mean for each treatment
  mr_pred <- post["b_lnmr_Intercept"] +  #Intercept represents mean MR for cold
    post["b_lnmr_temp"] * temp +   #Temperature effect 
    post["bsp_lnmr_milnmass"] * -1.421746 + #Mass effect * mean age
    post["b_lnmr_z_age"] * 0 + #Age effect + mean age
    post["sd_id__lnmr_Intercept"] + 
    post["sd_id__lnmr_temp"] * temp + 
    post["sd_samp_session__lnmr_Intercept"] + 
    post[id_intercept] + #Lizard intercept
    post[samp_session_intercept] #Samp session intercept
  
  if(treatment == 29){
    mr_pred + 
      post[,"b_lnmr_treatment29"] * 1 + 
      post[,"b_lnmr_treatment29:temp"] * temp
  }
  #Get the posterior summary
  output <- posterior_summary(mr_pred)[,-2] %>% round(2)
  
  df <- data.frame(id = id,
                   treatment = treatment,
                   samp_session = samp_session,
                   temp = temp,
                   mean_lnmr = output[1],
                   lower = output[2],
                   upper = output[3])
}


#Model prediction function for complete data sets
mr_pred <- function(id, temp, samp_session, treatment){
  #Lizard deviation
  id_intercept <- paste0("r_id[", id, ",", "Intercept]")
  id_slope <- paste0("r_id[", id, ",", "temp]")
  
  #Samp_session deviation
  samp_session_intercept <- paste0("r_samp_session[", samp_session, ",", "Intercept]")
  
  #Mean for each treatment
  mr_pred <- post["b_Intercept"] +  #Intercept represents mean MR for cold
    post["b_temp"] * temp +   #Temperature effect 
    post["b_lnmass"] * -1.421746 + #Mass effect * mean age
    post["b_z_age"] * 0 + #Age effect + mean age
    post["sd_id__Intercept"] + #Average ID deviation
    post["sd_id__temp"] * temp  +  #Average ID deviation for each temp
    post["sd_samp_session__Intercept"] +  #Average session deviation
    post[id_intercept] + #Lizard intercept
    post[samp_session_intercept] #Samp session intercept
  
  if(treatment == 29){
    mr_pred + 
      post[,"b_treatment29"] * 1 + 
      post[,"b_treatment29:temp"] * 1
  }
  #Get the posterior summary
  output <- posterior_summary(mr_pred)[,-2] %>% round(2)
  
  df <- data.frame(id = id,
                   treatment = treatment,
                   samp_session = samp_session,
                   temp = temp,
                   mean_lnmr = output[1],
                   lower = output[2],
                   upper = output[3])
}

#Function to calculate temperature specific repeatability for complete_case models
brms_rpt <- function(x, model, type = "mean"){
  #Strings to search for relevant (co)variance components
  id_vars <- paste0("sd_","id")
  id_cors <- paste0("cor_","id")
  
  #Extract the relevant sd/(co)variance components for ID 
  SD <- posterior_samples(model, id_vars)
  COR <- posterior_samples(model, id_cors) 
  V <- SD^2
  names(V) <- str_replace(names(V), "sd", "var")
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV <- cbind(COR[1] * (SD[1] * SD[2]))
  names(COV) <- str_replace(names(COV), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  V_id <- V[1] + (x^2)*V[2] +  #The V of the intercept and linear slope
    2*x*COV[1]    # Covariance of intercept and linear slope
  
  #Samp session variance
  #Strings to search for relevant (co)variance components
  ss_vars <- paste0("sd_samp")
  
  #Extract the relevant sd/(co)variance components 
  SD_SS <- posterior_samples(model, ss_vars)
  
  #Squaring SD to get the variance
  V_SS <- (SD_SS)^2 
  
  #Residuals
  SD_e <- posterior_samples(model, "sigma") # Extract the variance of intercept, linear slope
  
  #Squaring SD to get the variance
  V_resid <- (SD_e)^2 
  
  #Calculate total phenotypic variance
  VtotalP <- V_id + V_SS + V_resid 
  
  #    # Calculate repeatability (proportion of variance explained by consistent ID variance)
  R <- V_id  / VtotalP
  
  if(type == "mean"){
    df <- data.frame(temp = x,
                     rpt = posterior_summary(R)[1] %>% round(2),
                     Lower =  posterior_summary(R)[3] %>% round(2),
                     Upper =  posterior_summary(R)[4] %>% round(2))
  }
  
  if(type == "mode"){
    df <- data.frame(temp = x,
                     rpt = MCMCglmm::posterior.mode(as.mcmc(R)) %>% round(2),
                     Lower =  coda::HPDinterval(as.mcmc(R))[1] %>% round(2),
                     Upper =  coda::HPDinterval(as.mcmc(R))[2] %>% round(2))
  }
  
  
  return(df)
}

#Function to calculate temperature specific repeatability for imputation models
brms_rpt_imp <- function(x, model, type = "mean"){
  #Strings to search for relevant (co)variance components
  id_vars <- paste0("sd_","id")
  id_cors <- paste0("cor_","id")
  
  #Extract the relevant sd/(co)variance components for ID 
  SD <- posterior_samples(model, id_vars)
  COR <- posterior_samples(model, id_cors) 
  V <- SD^2
  names(V) <- str_replace(names(V), "sd", "var")
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV <- cbind(COR[1] * (SD[1] * SD[2]))
  names(COV) <- str_replace(names(COV), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  V_id <- V[1] + (x^2)*V[2] +  #The V of the intercept and linear slope
    2*x*COV[1]    # Covariance of intercept and linear slope
  
  #Samp session variance
  #Strings to search for relevant (co)variance components
  ss_vars <- paste0("sd_samp")
  
  #Extract the relevant sd/(co)variance components 
  SD_SS <- posterior_samples(model, ss_vars)
  
  #Squaring SD to get the variance
  V_SS <- (SD_SS)^2 
  
  #Residuals
  SD_e <- posterior_samples(model, "sigma_lnmr") # Extract the variance of intercept, linear slope
  
  #Squaring SD to get the variance
  V_resid <- (SD_e)^2 
  
  #Calculate total phenotypic variance
  VtotalP <- V_id + V_SS + V_resid 
  
  #    # Calculate repeatability (proportion of variance explained by consistent ID variance)
  R <- V_id  / VtotalP
  
  if(type == "mean"){
  df <- data.frame(temp = x,
                   rpt = posterior_summary(R)[1] %>% round(2),
                   Lower =  posterior_summary(R)[3] %>% round(2),
                   Upper =  posterior_summary(R)[4] %>% round(2))
  }
  
  if(type == "mode"){
    df <- data.frame(temp = x,
                     rpt = MCMCglmm::posterior.mode(as.mcmc(R)) %>% round(2),
                     Lower =  coda::HPDinterval(as.mcmc(R))[1] %>% round(2),
                     Upper =  coda::HPDinterval(as.mcmc(R))[2] %>% round(2))
  }
  
  return(df)
}

#Repeatability of slope
brms_rptSlope_cc <- function(model, type = "mean"){
  idslope <- posterior_samples(model, "sd_id__temp") 
  seriesslope <- posterior_samples(model, "sd_series__temp") 
  
  r.slope <- ( idslope / ( idslope + seriesslope) )
  
  if(type == "mean"){
    df <- data.frame(rptSlope = posterior_summary(r.slope)[1] %>% round(2),
                     Lower =  posterior_summary(r.slope)[3] %>% round(2),
                     Upper =  posterior_summary(r.slope)[4] %>% round(2))
  }
  
  if(type == "mode"){
    df <- data.frame(rptSlope = MCMCglmm::posterior.mode(as.mcmc(r.slope)) %>% round(2),
                     Lower =  coda::HPDinterval(as.mcmc(r.slope))[1] %>% round(2),
                     Upper =  coda::HPDinterval(as.mcmc(r.slope))[2] %>% round(2))
  }
  
  return(df)
}


brms_rptSlope_imp <- function(model, transform = "z", type = "mean"){
  if(transform == "z"){
    idslope <- posterior_samples(model, "sd_id__lnmr_z_temp") 
    seriesslope <- posterior_samples(model, "sd_series__lnmr_z_temp")  
  }
  
  if(transform == "raw"){
  idslope <- posterior_samples(model, "sd_id__lnmr_temp") 
  seriesslope <- posterior_samples(model, "sd_series__lnmr_temp")
  }
  
  r.slope <- ( idslope / ( idslope + seriesslope) )
  
  if(type == "mean"){
    df <- data.frame(rptSlope = posterior_summary(r.slope)[1] %>% round(2),
                     Lower =  posterior_summary(r.slope)[3] %>% round(2),
                     Upper =  posterior_summary(r.slope)[4] %>% round(2))
    return(df)
  }
  
  if(type == "mode"){
    df <- data.frame(rptSlope = MCMCglmm::posterior.mode(as.mcmc(r.slope)) %>% round(2),
                     Lower =  coda::HPDinterval(as.mcmc(r.slope))[1] %>% round(2),
                     Upper =  coda::HPDinterval(as.mcmc(r.slope))[2] %>% round(2))
    return(df)
  }
  if(type == "posterior"){
    r.slope
    return(r.slope)
    }

  
}

imp_ph_pred <- function(id, temp, treatment){
  #Lizard deviation
  id_intercept <- paste0("r_id__lnmr[", id, ",", "Intercept]")
  
  #Mean for each treatment
  mr_pred <- post["b_lnmr_Intercept"] +  #Intercept represents mean MR for cold
    post["b_lnmr_temp"] * temp +   #Temperature effect 
    post["bsp_lnmr_milnmass"] * -1.421746 + #Mass effect * mean age
    post["b_lnmr_z_age"] * 0 + #Age effect + mean age
    post["sd_id__lnmr_Intercept"] + 
    post["sd_id__lnmr_temp"] * temp + 
    post[id_intercept] #Lizard intercept

  if(treatment == 29){
    mr_pred + 
      post[,"b_lnmr_treatment29"] * 1 + 
      post[,"b_lnmr_treatment29:temp"] * 1
  }
  #Get the posterior summary
  output <- posterior_summary(mr_pred)[,-2] %>% round(2)
  
  df <- data.frame(id = id,
                   treatment = treatment,
                   temp = temp,
                   mean_lnmr = output[1],
                   lower = output[2],
                   upper = output[3])
}


rpt_posterior <- function(x, model){
    #Strings to search for relevant (co)variance components
    id_vars <- paste0("sd_","id")
    id_cors <- paste0("cor_","id")
    
    #Extract the relevant sd/(co)variance components for ID 
    SD <- posterior_samples(model, id_vars)
    COR <- posterior_samples(model, id_cors) 
    V <- SD^2
    names(V) <- str_replace(names(V), "sd", "var")
    #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
    COV <- cbind(COR[1] * (SD[1] * SD[2]))
    names(COV) <- str_replace(names(COV), "cor", "cov")
    
    # Now, add everything together while accounting for covariances and their respective powers
    V_id <- V[1] + (x^2)*V[2] +  #The V of the intercept and linear slope
      2*x*COV[1]    # Covariance of intercept and linear slope
    
    #Samp session variance
    #Strings to search for relevant (co)variance components
    ss_vars <- paste0("sd_samp")
    
    #Extract the relevant sd/(co)variance components 
    SD_SS <- posterior_samples(model, ss_vars)
    
    #Squaring SD to get the variance
    V_SS <- (SD_SS)^2 
    
    #Residuals
    SD_e <- posterior_samples(model, "sigma_lnmr") # Extract the variance of intercept, linear slope
    
    #Squaring SD to get the variance
    V_resid <- (SD_e)^2 
    
    #Calculate total phenotypic variance
    VtotalP <- V_id + V_SS + V_resid 
    
    #    # Calculate repeatability (proportion of variance explained by consistent ID variance)
    R <- V_id  / VtotalP
    
    return(R)
}


var_comp_imp <- function(x, model, var_comp = "id"){
  #Strings to search for relevant (co)variance components
  id_vars <- paste0("sd_","id")
  id_cors <- paste0("cor_","id")
  
  #Extract the relevant sd/(co)variance components for ID 
  SD <- posterior_samples(model, id_vars)
  COR <- posterior_samples(model, id_cors) 
  V <- SD^2
  names(V) <- str_replace(names(V), "sd", "var")
  #Convert correlation to covariance #Cor(1,2) * (SD1 X SD2) = Cov(1,2)
  COV <- cbind(COR[1] * (SD[1] * SD[2]))
  names(COV) <- str_replace(names(COV), "cor", "cov")
  
  # Now, add everything together while accounting for covariances and their respective powers
  V_id <- V[1] + (x^2)*V[2] +  #The V of the intercept and linear slope
    2*x*COV[1]    # Covariance of intercept and linear slope
  
  #Samp session variance
  #Strings to search for relevant (co)variance components
  ss_vars <- paste0("sd_samp")
  
  #Extract the relevant sd/(co)variance components 
  SD_SS <- posterior_samples(model, ss_vars)
  
  #Squaring SD to get the variance
  V_SS <- (SD_SS)^2 
  
  #Residuals
  SD_e <- posterior_samples(model, "sigma_lnmr") # Extract the variance of intercept, linear slope
  
  #Squaring SD to get the variance
  V_resid <- (SD_e)^2 
  
  #Calculate total phenotypic variance
  VtotalP <- V_id + V_SS + V_resid 
  
  if(var_comp == "id"){
    df <- data.frame(temp = x,
                     var_comp = var_comp,
                     variance = posterior_summary(V_id)[1] %>% round(2),
                     Lower =  posterior_summary(V_id)[3] %>% round(2),
                     Upper =  posterior_summary(V_id)[4] %>% round(2))

  }
  
  if(var_comp == "sigma"){
    df <- data.frame(temp = x,
                     var_comp = var_comp,
                     variance = posterior_summary(V_resid)[1] %>% round(2),
                     Lower =  posterior_summary(V_resid)[3] %>% round(2),
                     Upper =  posterior_summary(V_resid)[4] %>% round(2))
    
  }
  
  return(df)
}
