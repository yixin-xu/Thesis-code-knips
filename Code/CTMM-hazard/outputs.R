generate_model_outputs <- function(hesim_dat, model_inputs) {
  

  # Fixed parameters
  mr <- model_inputs$mr
  log_rate_1st_revision_3_times <- model_inputs$log_rate_1st_revision_3_times
  log_rate_1st_revision_310_times <- model_inputs$log_rate_1st_revision_310_times
  log_rate_1st_revision_10_times <- model_inputs$log_rate_1st_revision_10_times
  
  
  revision_cost <- model_inputs$revision_cost
  
  # Random parameters
  input_parameters <- model_inputs$inputs_parameters
  log_rate_early_second <- input_parameters[, grep("log_rate_early_second", colnames(input_parameters))]
  log_rate_middle_second <- input_parameters[, grep("log_rate_middle_second", colnames(input_parameters))]
  log_rate_late_second <- input_parameters[, grep("log_rate_late_second", colnames(input_parameters))]
  log_rate_higher <- input_parameters[, grep("log_rate_higher", colnames(input_parameters))]
  
  log_rate_1st_revision_3 <- list()
  for(implant_name in implant_names) {
    temp <- input_parameters[, grep(paste0("log_rate_1st_revision_3_", implant_name), colnames(input_parameters))]
    
    colnames(temp) <- gsub(paste0("log_rate_1st_revision_3_", implant_name, "_"), "", colnames(temp))
    log_rate_1st_revision_3[[implant_name]] <- temp
  }
  
  log_rate_1st_revision_310 <- list()
  for(implant_name in implant_names) {
    temp <- input_parameters[, grep(paste0("log_rate_1st_revision_310_", implant_name), colnames(input_parameters))]
    
    colnames(temp) <- gsub(paste0("log_rate_1st_revision_310_", implant_name, "_"), "", colnames(temp))
    log_rate_1st_revision_310[[implant_name]] <- temp
  }
  
  log_rate_1st_revision_10 <- list()
  for(implant_name in implant_names) {
    temp <- input_parameters[, grep(paste0("log_rate_1st_revision_10_", implant_name), colnames(input_parameters))]
    
    colnames(temp) <- gsub(paste0("log_rate_1st_revision_10_", implant_name, "_"), "", colnames(temp))
    log_rate_1st_revision_10[[implant_name]] <- temp
  }
  
  if(is.infinite(age_range[2])) {
    log_rate_1st_revision_10[[implant_name]][, 2] <- -100}
  
  tkr_mortality_90d_probability <- input_parameters[, "tkr_mortality_90d_probability"]
  revision_mortality_90d_probability <- input_parameters[, "revision_mortality_90d_probability"]
  rerevision_mortality_90d_probability <- input_parameters[, "rerevision_mortality_90d_probability"]
  
  revision_disutility <- input_parameters[, "revision_disutility"]
  utility_post_tkr <- input_parameters[, "utility_post_tkr"]
  utility_post_1st_rev <- input_parameters[, "utility_post_1st_rev"]
  utility_post_2nd_rev <- input_parameters[, "utility_post_2nd_rev"]
  
  implant_costs <- input_parameters[, grep("implant_cost_", colnames(input_parameters))]
  colnames(implant_costs) <- gsub("implant_cost_", "", colnames(implant_costs))
  implant_costs <- t(implant_costs)
  
  ############################################################
  ## Simulate disease progression ############################
  ############################################################
  
  # Nested loop over implants
  # Allows independent spline models for each implant
  # Create sample parameters with each row repeated n_patients times
  # Simulate to get disprog
  # Simulate again using yrs_to_1st_rev_simulated
  # Relabel the resulting disprog so that samples go from 1:n_samples and patients from 1:n_patients
  
  
  
  # List of disease progression and state occupancies, one for each implant
  disprog <- list()
  stateprobs <- list()
  
  for(implant_name in implant_names) {
    print(paste0("Implant ", which(implant_names == implant_name), "/", length(implant_names)))
    
    # Treating patients as separate samples 
    n_samples_temp <- n_samples * n_patients
    
    # Simple function to convert vector to matrix with columns named constant ('cons')
    matrixv <- function(v, n = NULL){
      if (length(v) == 1) v <- rep(v, n_samples_temp) 
      m <- matrix(v)
      colnames(m) <- "cons"
      return(m)
    }
    
    surgery_d_second_early = log(-log(1-(1-exp(-exp(log_rate_early_second)))*rerevision_mortality_90d_probability))
    surgery_d_second_middle = log(-log(1-(1-exp(-exp(log_rate_middle_second)))*rerevision_mortality_90d_probability))
    surgery_d_second_late = log(-log(1-(1-exp(-exp(log_rate_late_second)))*rerevision_mortality_90d_probability))
    surgery_d_higher = log(-log(1-(1-exp(-exp(log_rate_higher)))*rerevision_mortality_90d_probability))
    

   ## combined for this implant

    # Define the transition rate models
    # Function to generate random samples of underlying model parameters
    # Create a list of randomly sampled survival models for transitions
    transition_model_params <- params_surv_list(
      # 1. Post TKR to Post early 1st revision
      params_surv(coefs = list(
        rep(log_rate_1st_revision_3[[implant_name]][, 1], each = n_patients),
        rep(log_rate_1st_revision_3[[implant_name]][, 2], each = n_patients)),
        aux = list(time = log_rate_1st_revision_3_times,
                   scale = "log_hazard"),
        dist = "pwexp"),
      
      # 2. Post TKR to Post middle 1st revision
      params_surv(coefs = list(
        rep(log_rate_1st_revision_310[[implant_name]][, 1], each = n_patients),
        rep(log_rate_1st_revision_310[[implant_name]][, 2], each = n_patients),
        rep(log_rate_1st_revision_310[[implant_name]][, 3], each = n_patients)),
        aux = list(time = log_rate_1st_revision_310_times,
                   scale = "log_hazard"),
        dist = "pwexp"),
      
      # 3. Post TKR to Post late 1st revision
      params_surv(coefs = list(
        rep(log_rate_1st_revision_10[[implant_name]][, 1], each = n_patients),
        rep(log_rate_1st_revision_10[[implant_name]][, 2], each = n_patients)),
        aux = list(time = log_rate_1st_revision_10_times,
                   scale = "log_hazard"),
        dist = "pwexp"),
      
      # 4. Post early revision to Post second revision
      params_surv(coefs = lapply(as.list(log_rate_early_second),matrixv),
        aux = list(scale = "log_hazard"),
        dist = "exp"),
      
      # 5. Post middle revision to Post second revision
      params_surv(coefs = lapply(as.list(log_rate_middle_second),matrixv),
                  aux = list(scale = "log_hazard"),
                  dist = "exp"),
      
      # 6. Post late revision to Post second revision
      params_surv(coefs = lapply(as.list(log_rate_late_second),matrixv),
                  aux = list(scale = "log_hazard"),
                  dist = "exp"),
      
      # 7-11. Post TKR to death
      params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                  aux = list(time = c(0:time_horizon),
                             scale = "log_hazard"),
                  dist = "pwexp"),
      params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                  aux = list(time = c(0:time_horizon),
                             scale = "log_hazard"),
                  dist = "pwexp"),
      params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                  aux = list(time = c(0:time_horizon),
                             scale = "log_hazard"),
                  dist = "pwexp"),
      params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                  aux = list(time = c(0:time_horizon),
                             scale = "log_hazard"),
                  dist = "pwexp"),
      params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                  aux = list(time = c(0:time_horizon),
                             scale = "log_hazard"),
                  dist = "pwexp"),
      
      # 12. Post TKR to surgery death
      params_surv(coefs = list(
        rep(log(-log(1-(1-exp(-exp(log_rate_1st_revision_3[[implant_name]][, 1])))*revision_mortality_90d_probability)), each = n_patients),
        rep(log(-log(1-(1-exp(-exp(log_rate_1st_revision_310[[implant_name]][, 2])))*revision_mortality_90d_probability)), each = n_patients),
        rep(log(-log(1-(1-exp(-exp(log_rate_1st_revision_10[[implant_name]][, 2])))*revision_mortality_90d_probability)), each = n_patients)),
        aux = list(time = c(0,3,10),
                   scale = "log_hazard"),
        dist = "pwexp"),
      
      
      # 13-15. Post 1st revision to surgery death
      params_surv(coefs = lapply(as.list(surgery_d_second_early),matrixv),
                   aux = list(scale = "log_hazard"),
                   dist = "exp"),
      params_surv(coefs = lapply(as.list(surgery_d_second_middle),matrixv),
                  aux = list(scale = "log_hazard"),
                  dist = "exp"),
      params_surv(coefs = lapply(as.list(surgery_d_second_late),matrixv),
                  aux = list(scale = "log_hazard"),
                  dist = "exp"),
      
      # 16. Post 2nd revision to surgery death
      # Depends on initial 90d period and ongoing risk of higher revision
      params_surv(coefs = lapply(as.list(surgery_d_higher),matrixv),
                  aux = list(scale = "log_hazard"),
                  dist = "exp")
      
    )
    

    # Now simulate the outcomes for all patients for this sample and this implant
    # Need temporary hesim data object with only one strategy and one patient
    hesim_dat_temp <- hesim_data(strategies = strategies[strategies$strategy_name == implant_name],
                                 patients = patients[1, ],  # All patients assumed identical at baseline
                                 states = states)
    
    transition_model_data <- expand(hesim_dat_temp, by = c("strategies", "patients"))
    
    
    transition_model_data[, cons := 1]
    transition_model_data[, x1 := 1]
    # Transition model
    transition_model <- create_IndivCtstmTrans(transition_model_params, 
                                               input_data = transition_model_data,
                                               trans_mat = tmat,
                                               clock = "forward",
                                               start_age = initial_age)

    
    set.seed(2243534)
    
    disprog[[implant_name]] <- transition_model$sim_disease(max_t = time_horizon, max_age = (starting_age + time_horizon))
    
    # Need to change sample and patient_id to match expected format for hesim
    disprog[[implant_name]]$patient_id <- disprog[[implant_name]]$sample %% n_patients
   
    # Also simulate average state probabilities
    stateprobs[[implant_name]] <- transition_model$sim_stateprobs(t = 0:(time_horizon),
                                                                  disprog = disprog[[implant_name]])
    
    
     # Replace patient 0 with patient n_patient 
    disprog[[implant_name]]$patient_id[which(disprog[[implant_name]]$patient_id == 0)] <- n_patients
    disprog[[implant_name]]$sample <- ceiling(disprog[[implant_name]]$sample / n_patients)
    
    # Add correct implant strategy number
    disprog[[implant_name]]$strategy_id <- which(implant_names == implant_name)
    
    
  } # End loop over implants
  
  disprog_combined <- rbindlist(disprog)
  stateprobs_combined <- rbindlist(stateprobs)
  
  # Correctly set the size attributes
  # Check if using attributes(disprog_combined)
  attributes_size <- c(n_samples, n_strategies, n_patients, n_states)
  names(attributes_size) <- c("n_samples", "n_strategies", "n_patients","n_states")
  setattr(disprog_combined, "size", attributes_size)

  
  # Create disease progression combined for this implant
  
  ############################################################
  ## Combined economic model #################################
  ############################################################
  
  
  cost_utility_models <- generate_cost_utility_models(n_samples = n_samples,
                                                      hesim_dat = hesim_dat,
                                                      model_inputs = model_inputs)
                                                      
  
  # Economic model combining everything
  economic_model <- IndivCtstm$new(trans_model = transition_model,
                            utility_model = cost_utility_models$utility_model,
                            cost_models = cost_utility_models$cost_models)
  
  
  # Use the combined disease progression results
  economic_model$disprog_ <- disprog_combined
  economic_model$stateprobs_ <- stateprobs_combined
  
  # Simulate QALYs and costs with discount rates
  # 3.5% discount rate
  economic_model$sim_qalys(dr = 0.035)
  economic_model$sim_costs(dr = 0.035)
  # Divide by number of patients
  economic_model$costs_[, c("costs")] <- economic_model$costs_[, c("costs")] / n_patients
  economic_model$qalys_[, c("qalys", "lys")] <- economic_model$qalys_[, c("qalys", "lys")] / n_patients
  
  # Return for analysis
  model_outputs <- list("economic_model" = economic_model)
  
}

