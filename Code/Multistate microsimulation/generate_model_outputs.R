generate_model_outputs <- function(hesim_dat, model_inputs) {
  
  # Extract requisite random and fixed parameters from model_inputs
  # This is to ensure input_parameters aligns with format required for BCEA-EVPPI
  
  # Fixed parameters
  ln_bhknots_first_revision <- model_inputs$ln_bhknots_first_revision
  ln_bhknots_second_revision <- model_inputs$ln_bhknots_second_revision
  mr <- model_inputs$mr
  tkr_mortality_times_90d <- model_inputs$tkr_mortality_times_90d
  revision_mortality_times_90d <- model_inputs$revision_mortality_times_90d
  rerevision_mortality_times_90d <- model_inputs$rerevision_mortality_times_90d
  revision_cost <- model_inputs$revision_cost
  
  # Random parameters
  input_parameters <- model_inputs$inputs_parameters
  rcs_second_revision <- input_parameters[, grep("rcs_second_revision", colnames(input_parameters))]
  colnames(rcs_second_revision) <- gsub("rcs_second_revision_", "", colnames(rcs_second_revision))
  
  probability_higher_revision <- input_parameters[, "probability_higher_revision"]
  
  rcs_first_revision <- list()
  for(implant_name in implant_names) {
    temp <- input_parameters[, grep(paste0("rcs_first_revision_", implant_name), colnames(input_parameters))]
    
    colnames(temp) <- gsub(paste0("rcs_first_revision_", implant_name, "_"), "", colnames(temp))
    rcs_first_revision[[implant_name]] <- temp
  }
  
  tkr_mortality_90d <- input_parameters[, grep("tkr_mortality_90d_", colnames(input_parameters))]
  revision_mortality_90d <- input_parameters[, grep("revision_mortality_90d_", colnames(input_parameters))]
  rerevision_mortality_90d <- input_parameters[, grep("rerevision_mortality_90d_", colnames(input_parameters))]
  
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
  yrs_to_1st_early_rev_implant <- list()
  yrs_to_1st_middle_rev_implant <- list()
  yrs_to_1st_late_rev_implant <- list()
  yrs_to_2nd_rev_implant <- list()
  
  
  
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
    
    # Times to 1st revision needed for 2nd revision rates and 1st revision disutility/cost
    yrs_to_1st_early_rev_simulated <- rep(0, n_samples_temp)
    yrs_to_1st_middle_rev_simulated <- rep(0, n_samples_temp)
    yrs_to_1st_late_rev_simulated <- rep(0, n_samples_temp)
    
    
    # Times to 2nd revision needed for 2nd revision disutility/cost
    yrs_to_2nd_rev_simulated <- rep(0, n_samples_temp)
    yrs_to_2nd_rev_simulated_value <- rep(0, n_samples_temp)
    
    # Function to add years to first revision before time in a clock reset state
    # so that correct rate from piecewise exponential is used
    offset_pwexp <- function(rates, start_times) {
      # Append empty values that aren't used as they're outside time horizon
      rates <- rbind(rates, matrix(rep(rates[length(rates)], max(start_times)), ncol = 1))
      start_times <- floor(start_times)
      output_list <- list()
      for(i_time in 1:length(rates)) {
        # Time is start_times later than i_time
        output_list[[i_time]] <- matrix(rates[start_times + i_time], ncol = 1)
        colnames(output_list[[i_time]]) <- "cons"
      }
      return(output_list)
    }
    
    
    # Define the transition rate models
    # Function to generate random samples of underlying model parameters
    # Create a list of randomly sampled survival models for transitions
    transition_model_params <- params_surv_list(
      # 1. Post TKR to Temp revision
      params_surv(coefs = list(
        gamma0 = data.frame(int = rep(rcs_first_revision[[implant_name]][, "cons"], each = n_patients)),
        gamma1 = data.frame(int = rep(rcs_first_revision[[implant_name]][, "rcs1"], each = n_patients)),
        gamma2 = data.frame(int = rep(rcs_first_revision[[implant_name]][, "rcs2"], each = n_patients)),
        gamma3 = data.frame(int = rep(rcs_first_revision[[implant_name]][, "rcs3"], each = n_patients))),  
        aux = list(
          knots = ln_bhknots_first_revision[[implant_name]],
          scale = "log_cumhazard",
          timescale = "log"),
        dist = "survspline"), 
      
      
      # 2. Temp revision to early revision
      params_surv(coefs = list(rep(100, each = n_patients*n_samples),
                               rep(0, each = n_patients*n_samples)),
                  aux = list(time = c(0,3)),
                  dist = "pwexp"),
      
      # 3. Temp revision to middle revision
      params_surv(coefs = list(rep(0, each = n_patients*n_samples),
                               rep(100, each = n_patients*n_samples),
                               rep(0, each = n_patients*n_samples)),
                  aux = list(time = c(0,3,10)),
                  dist = "pwexp"),
      
      # 4. Temp revision to late revision
      params_surv(coefs = list(rep(0, each = n_patients*n_samples),
                               rep(100, each = n_patients*n_samples)),
                  aux = list(time = c(0,10)),
                  dist = "pwexp"),
      
      # 5. Post early revision to post 2nd revision
      # Includes dependence on time from TKR to 1st revision
      params_surv(coefs = list(
        gamma0 = data.frame(int = rep(rcs_second_revision[, "cons"], each = n_patients) +
                              yrs_to_1st_early_rev_simulated*rep(rcs_second_revision[, "yrs_to_1st_rev"], each = n_patients)
        ),
        gamma1 = data.frame(int = rep(rcs_second_revision[, "rcs1"], each = n_patients)),
        gamma2 = data.frame(int = rep(rcs_second_revision[, "rcs2"], each = n_patients)),
        gamma3 = data.frame(int = rep(rcs_second_revision[, "rcs3"], each = n_patients))),  
        aux = list(
          knots = ln_bhknots_second_revision,
          scale = "log_cumhazard",
          timescale = "log"),
        dist = "survspline"),
      
      
      
      # 6. Post middle revision to post 2nd revision
      # Includes dependence on time from TKR to 1st revision
      params_surv(coefs = list(
        gamma0 = data.frame(int = rep(rcs_second_revision[, "cons"], each = n_patients) +
                              yrs_to_1st_middle_rev_simulated*rep(rcs_second_revision[, "yrs_to_1st_rev"], each = n_patients)
        ),
        gamma1 = data.frame(int = rep(rcs_second_revision[, "rcs1"], each = n_patients)),
        gamma2 = data.frame(int = rep(rcs_second_revision[, "rcs2"], each = n_patients)),
        gamma3 = data.frame(int = rep(rcs_second_revision[, "rcs3"], each = n_patients))),  
        aux = list(
          knots = ln_bhknots_second_revision,
          scale = "log_cumhazard",
          timescale = "log"),
        dist = "survspline"),
      
      # 7. Post late revision to post 2nd revision
      # Includes dependence on time from TKR to 1st revision
      params_surv(coefs = list(
        gamma0 = data.frame(int = rep(rcs_second_revision[, "cons"], each = n_patients) +
                              yrs_to_1st_late_rev_simulated*rep(rcs_second_revision[, "yrs_to_1st_rev"], each = n_patients)
        ),
        gamma1 = data.frame(int = rep(rcs_second_revision[, "rcs1"], each = n_patients)),
        gamma2 = data.frame(int = rep(rcs_second_revision[, "rcs2"], each = n_patients)),
        gamma3 = data.frame(int = rep(rcs_second_revision[, "rcs3"], each = n_patients))),  
        aux = list(
          knots = ln_bhknots_second_revision,
          scale = "log_cumhazard",
          timescale = "log"),
        dist = "survspline"),
      
      # 8. Post TKR to death
      params_surv(coefs = lapply(as.list(log(mr)), matrixv),
                  aux = list(time = c(0:time_horizon)),
                  dist = "pwexp"),
      
      # 9. Post early 1st revision to death
      params_surv(coefs = offset_pwexp(rates = log(mr), start_times = yrs_to_1st_early_rev_simulated),
                  aux = list(time = c(0:time_horizon)),
                  dist = "pwexp"),
      
      
      # 10. Post middle 1st revision to death
      params_surv(coefs = offset_pwexp(rates = log(mr), start_times = yrs_to_1st_middle_rev_simulated),
                  aux = list(time = c(0:time_horizon)),
                  dist = "pwexp"),
      
      # 11. Post late 1st revision to death
      params_surv(coefs = offset_pwexp(rates = log(mr), start_times = yrs_to_1st_late_rev_simulated),
                  aux = list(time = c(0:time_horizon)),
                  dist = "pwexp"),
      
      # 12. Post 2nd revision to death
      params_surv(coefs = offset_pwexp(rates = log(mr), start_times = yrs_to_2nd_rev_simulated),
                  aux = list(time = c(0:time_horizon)),
                  dist = "pwexp"),
      
      
      # 13. Post TKR to surgery death
      params_surv(coefs = list(rate = rep(log(tkr_mortality_90d[, 1]), each = n_patients),
                               as.matrix(rep(-100, n_samples_temp))),
                  aux = list(time = tkr_mortality_times_90d),
                  dist = "pwexp"),
      
      # 14, 15, 16. Post 1st revision to surgery death
      params_surv(coefs = list(rate = rep(log(revision_mortality_90d[, 1]), each = n_patients),
                               as.matrix(rep(-100, n_samples_temp))),
                  aux = list(time = revision_mortality_times_90d),
                  dist = "pwexp"),
      params_surv(coefs = list(rate = rep(log(revision_mortality_90d[, 1]), each = n_patients),
                               as.matrix(rep(-100, n_samples_temp))),
                  aux = list(time = revision_mortality_times_90d),
                  dist = "pwexp"),
      params_surv(coefs = list(rate = rep(log(revision_mortality_90d[, 1]), each = n_patients),
                               as.matrix(rep(-100, n_samples_temp))),
                  aux = list(time = revision_mortality_times_90d),
                  dist = "pwexp"),
      
      # 17. Post 2nd revision to surgery death
      # Depends on initial 90d period and ongoing risk of higher revision
      params_surv(coefs = list(rep(log(rerevision_mortality_90d[, 1]), each = n_patients),
                               rep(log(rerevision_mortality_90d[, 2]), each = n_patients)),
                  aux = list(time = rerevision_mortality_times_90d),
                  dist = "pwexp")
      
    )
    
    # Now simulate the outcomes for all patients for this sample and this implant
    # Need temporary hesim data object with only one strategy and one patient
    hesim_dat_temp <- hesim_data(strategies = strategies[strategies$strategy_name == implant_name],
                                 patients = patients[1, ],  # All patients assumed identical at baseline
                                 states = states)
    
    transition_model_data <- expand(hesim_dat_temp, by = c("strategies", "patients"))
    
    # Add numerical covariates
    # Need to add constants to match the spline and piecewise constant death parameters
    transition_model_data[, cons := 1]
    transition_model_data[, int := 1]
    transition_model_data[, x1 := 1]
    #transition_model_data[, yrs_to_1st_rev := 0]
    
    # Transition model
    transition_model <- create_IndivCtstmTrans(transition_model_params, 
                                               input_data = transition_model_data,
                                               trans_mat = tmat,
                                               clock = "mix",
                                               start_age = starting_age,
                                               reset_states = c(1))
    
    disprog[[implant_name]] <- transition_model$sim_disease(max_t = time_horizon, max_age = (starting_age + time_horizon))
    
    # Update the data using results of this initial simulation
    # Use the time to transition from state 1 to state 2
    # Not all patients on all transitions actually move from 1 to 2
    disprog_temp = disprog[[implant_name]][from == 2 & to == 3, ]
    yrs_to_1st_early_rev_simulated[disprog_temp$sample] <- disprog_temp$time_stop
    
    disprog_temp = disprog[[implant_name]][from == 2 & to == 4, ]
    yrs_to_1st_middle_rev_simulated[disprog_temp$sample] <- disprog_temp$time_stop
    
    disprog_temp = disprog[[implant_name]][from == 2 & to == 5, ]
    yrs_to_1st_late_rev_simulated[disprog_temp$sample] <- disprog_temp$time_stop
    
    disprog_temp = disprog[[implant_name]][from %in% c(3, 4, 5) & to == 6, ]
    yrs_to_2nd_rev_simulated[disprog_temp$sample] = disprog_temp$time_stop
    
    
    # Re-run simulation with the (correct) time to first revision
    # Transition model
    transition_model <- create_IndivCtstmTrans(transition_model_params, 
                                               input_data = transition_model_data,
                                               trans_mat = tmat,
                                               clock = "mix",
                                               start_age = starting_age,
                                               start_state = 1,
                                               reset_states = c(1,3,4,5,6))
    
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
    
    
    
    
    # Check that state probs make sense
    # with(stateprobs[[implant_name]], prob[strategy_id == which(implant_names == implant_name) & sample == 1 & t==20]) / n_patients
    
  } # End loop over implants
  
  # Create disease progression combined for this implant
  disprog_combined <- rbindlist(disprog)
  stateprobs_combined <- rbindlist(stateprobs)
  
  
  # Correctly set the size attributes
  # Check if using attributes(disprog_combined)
  attributes_size <- c(n_samples, n_strategies, n_patients, n_states + 2)
  names(attributes_size) <- c("n_samples", "n_strategies", "n_patients", "n_states")
  setattr(disprog_combined, "size", attributes_size)
  # Correctly set the absorbing attribute
  attributes_absorbing <- c(7, 8)
  names(attributes_absorbing) <- c("Death", "Surgery Death")
  setattr(disprog_combined, "absorbing", attributes_absorbing)
  
  # Check probabilities add to one
  # with(stateprobs[[]], prob[sample == 1 & strategy_id == 12 & t == 15])
  # with(stateprobs_combined_implant[[12]], prob[sample == 1 & strategy_id == 12 & t == 15])
  # with(stateprobs_combined, prob[sample == 1 & strategy_id == 12 & t == 15]) / n_patients
  
  
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
  
  # Look at simulations to ensure the combination was correct
  #head(economic_model$disprog_)
  
  # Simulate QALYs and costs with discount rates
  # 3.5% discount rate
  economic_model$sim_qalys(dr = discount_rate)
  economic_model$sim_costs(dr = discount_rate)
  # Divide by number of patients
  economic_model$costs_[, c("costs")] <- economic_model$costs_[, c("costs")] / n_patients
  economic_model$qalys_[, c("qalys", "lys")] <- economic_model$qalys_[, c("qalys", "lys")] / n_patients
  
  
  # Return for analysis
  model_outputs <- list("economic_model" = economic_model)
  
}

