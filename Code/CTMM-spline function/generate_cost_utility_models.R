 
generate_cost_utility_models <- function(n_samples, 
                                         hesim_dat,
                                         model_inputs) {
  
  # Extract requisite random and fixed parameters from model_inputs
  # This is to ensure input_parameters aligns with format required for BCEA-EVPPI
  
  revision_cost <- model_inputs$revision_cost
  
  # Random parameters
  input_parameters <- model_inputs$inputs_parameters
  probability_higher_revision <- input_parameters[, "probability_higher_revision"]
  ln_bhknots_first_revision <- model_inputs$ln_bhknots_first_revision

  rcs_first_revision <- list()
  for(implant_name in implant_names) {
    temp <- input_parameters[, grep(paste0("rcs_first_revision_", implant_name), colnames(input_parameters))]
    
    colnames(temp) <- gsub(paste0("rcs_first_revision_", implant_name, "_"), "", colnames(temp))
    rcs_first_revision[[implant_name]] <- temp
  }
  
  revision_disutility <- input_parameters[, "revision_disutility"]
  utility_post_tkr <- input_parameters[, "utility_post_tkr"]
  utility_post_1st_rev <- input_parameters[, "utility_post_1st_rev"]
  utility_post_2nd_rev <- input_parameters[, "utility_post_2nd_rev"]
  
  implant_costs <- input_parameters[, grep("implant_cost_", colnames(input_parameters))]
  colnames(implant_costs) <- gsub("implant_cost_", "", colnames(implant_costs))
  implant_costs <- t(implant_costs)
  
  lograte_revision <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "revision_log_rate")
  lograte_revision[,2:4]= log(lograte_revision[,2:4])
  lograte_revision <- as.data.frame(lograte_revision)
  rownames(lograte_revision) <- lograte_revision[, "parameter"]
  log_rate_2nd_revision_early_mean <- as.numeric(lograte_revision["early_revision", "mean"]) 
  log_rate_2nd_revision_early_se = ((as.numeric(lograte_revision["early_revision", "UL"])-as.numeric(lograte_revision["early_revision", "LL"]))/2*1.96)
  log_rate_2nd_revision_middle_mean <- as.numeric(lograte_revision["middle_revision", "mean"]) 
  log_rate_2nd_revision_middle_se = ((as.numeric(lograte_revision["middle_revision", "UL"])-as.numeric(lograte_revision["middle_revision", "LL"]))/2*1.96) 
  log_rate_2nd_revision_late_mean <- as.numeric(lograte_revision["late_revision", "mean"]) 
  log_rate_2nd_revision_late_se = ((as.numeric(lograte_revision["late_revision", "UL"])-as.numeric(lograte_revision["late_revision", "LL"]))/2*1.96) 
  other_rates_raw <- as.matrix(read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "other_rates"))
  rownames(other_rates_raw) <- other_rates_raw[, "parameter"]
  lograte_higher_revision_mean <-  as.numeric(other_rates_raw["lograte_higher_revision_mean", "value"])
  lograte_higher_revision_se <-  as.numeric(other_rates_raw["lograte_higher_revision_se", "value"])
  
  
  n_times = time_horizon
  
  
  # Need probabilities of 1st revision in each 1-year time interval
  # Probability 1st revision depends on strategy, sample, and time (not patient)
  probability_1st_revision <- data.table(strategy_id = rep(strategies$strategy_id, each = (n_samples * n_times)),
                                         sample = rep(rep(1:n_samples, times = n_strategies), each = n_times),
                                         time_start = rep(0:(n_times - 1), times = (n_strategies * n_samples)),
                                         value = rep(0, times = n_strategies * n_samples * n_times))
  
  
  # For efficiency, create a single matrix with spline parameters
  rcs_first_revision_temp <- rbindlist(rcs_first_revision)
  # And a matrix with one set of knots for each sample
  ln_bhknots_first_revision_temp <- rbindlist(lapply(ln_bhknots_first_revision, as.data.frame))
  ln_bhknots_first_revision_temp <- ln_bhknots_first_revision_temp[rep(1:12, each = n_samples), ]
  
  # Intervals are 1 year if time_horizon is equal to n_time_intervals
  for(time_start_ in 0:(n_times - 1)) {
    # Varies by strategy and sample
    probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"] <- 
      (exp(-Hsurvspline(
        x = time_start_,
        gamma = as.matrix(rcs_first_revision_temp),
        knots = as.matrix(ln_bhknots_first_revision_temp))) - 
      exp(-Hsurvspline(
        x = time_start_ + time_horizon/n_times,
        gamma = as.matrix(rcs_first_revision_temp),
        knots = as.matrix(ln_bhknots_first_revision_temp)
      )))/exp(-Hsurvspline(
        x = time_start_,
        gamma = as.matrix(rcs_first_revision_temp),
        knots = as.matrix(ln_bhknots_first_revision_temp)))
  }
  
  
  utility_tbl <- stateval_tbl(
    data.table(strategy_id = rep(c(1:n_strategies), each = n_samples * n_states * n_times),
               sample = rep(rep(1:n_samples, each = n_states * n_times), times = n_strategies),
               state_id = rep(rep(c(1:n_states), each = n_times), times = n_strategies * n_samples),
               time_start = rep(c(0: (n_times-1)), times = n_strategies * n_samples * n_states),
               value = 0),
    dist = "custom"
  )
  
  # Utilities in Post TKR are dependent on time (clock-forward as all start in Post TKR)
  for(time_start_ in 0:(n_times - 1)) {
    utility_tbl[utility_tbl$state_id == 1 & utility_tbl$time_start == time_start_, "value"] <- rep(utility_post_tkr, each = n_strategies) +# general utility
      # Disutility times probability of revision during this time interval
      rep(revision_disutility, each = n_strategies) *
      unlist(probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"])
  }
  
  
  
  utility_tbl[state_id == 2, "value"] <- 0
  utility_tbl[state_id == 3, "value"] <- rep(utility_post_tkr, each = n_strategies) + revision_disutility*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                     mean = log_rate_2nd_revision_early_mean,
                                                                                                     sd = log_rate_2nd_revision_early_se))))
  utility_tbl[state_id == 4, "value"] <- rep(utility_post_tkr, each = n_strategies) + revision_disutility*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                     mean = log_rate_2nd_revision_middle_mean,
                                                                                                     sd = log_rate_2nd_revision_middle_se))))
  utility_tbl[state_id == 5, "value"] <- rep(utility_post_tkr, each = n_strategies) + revision_disutility* (1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                      mean = log_rate_2nd_revision_late_mean,
                                                                                                      sd = log_rate_2nd_revision_late_se))))
  # Utility post 2nd revision 
  # Need to subtract disutilities of higher revision multiplied by (annual) probability of having revision
  # Don't depend on time, time to second revision, or strategy
  utility_tbl[state_id == 6, "value"] <- rep(utility_post_tkr, each = n_strategies) + revision_disutility*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                     mean = lograte_higher_revision_mean,
                                                                                                     sd = lograte_higher_revision_se)))) 
            
  
  # Correct the utility table so that no utilities applied after the time_horizon
  utility_tbl$time_stop[is.infinite(utility_tbl$time_stop)] <- time_horizon
  utility_tbl$time_stop[utility_tbl$time_stop > time_horizon] <- time_horizon
  
  # Implant costs only in "Post TKR" state and depends on strategy
  implantcost_tbl <- stateval_tbl(
    data.table(strategy_id = rep(strategies$strategy_id, each = n_states * n_samples),
               sample = rep(rep(1:n_samples, each = n_states), times = n_strategies),
               state_id = rep(states$state_id, times = n_strategies * n_samples),
               value = 0),
    dist = "custom"
  )
  
  # Only non-zero for first state
  for(i_implant in 1:n_strategies) {
    implantcost_tbl[implantcost_tbl$strategy_id == i_implant &
                      implantcost_tbl$state_id == 1, "value"] = implant_costs[i_implant, ] 
  }
  
  
  # Medical costs
  medcost_tbl <- stateval_tbl(
    data.table(strategy_id = rep(c(1:n_strategies), each = n_samples * n_states * n_times),
               sample = rep(rep(1:n_samples, each = n_states * n_times), times = n_strategies),
               state_id = rep(rep(c(1:n_states), each = n_times), times = n_strategies * n_samples),
               time_start = rep(c(0: (n_times-1)), times = n_strategies * n_samples * n_states),
               value = 0),
    dist = "custom"
  )
  
  # Costs in Post TKR are dependent on time (clock-forward as all start in Post TKR)
  for(time_start_ in 0:(n_times - 1)) {
    medcost_tbl[medcost_tbl$state_id == 1 & medcost_tbl$time_start == time_start_, "value"] <-  rep(revision_cost, each = n_strategies ) *
      unlist(rep(probability_1st_revision[probability_1st_revision$time_start == time_start_, "value"]))
  }
  
  medcost_tbl[state_id == 2, "value"] <- 0
  # No cost in post 1st revision
  medcost_tbl[state_id == 3, "value"] <- revision_cost*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                            mean = log_rate_2nd_revision_early_mean,
                                                                            sd = log_rate_2nd_revision_early_se))))
  medcost_tbl[state_id == 4, "value"] <- revision_cost*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                            mean = log_rate_2nd_revision_middle_mean,
                                                                            sd = log_rate_2nd_revision_middle_se))))
  medcost_tbl[state_id == 5, "value"] <- revision_cost* (1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                             mean = log_rate_2nd_revision_late_mean,
                                                                             sd = log_rate_2nd_revision_late_se))))
  # Cost post 2nd revision 
  # Cost of a higher revision multiplied by (annual) probability of having revision
  # Don't depend on time, time to second revision, or strategy
  medcost_tbl[state_id == 6, "value"] <-  revision_cost*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                             mean = lograte_higher_revision_mean,
                                                                             sd = lograte_higher_revision_se)))) 
  
  # Correct the utility table so that no utilities applied after the time_horizon
  medcost_tbl$time_stop[is.infinite(medcost_tbl$time_stop)] <- time_horizon
  medcost_tbl$time_stop[medcost_tbl$time_stop > time_horizon] <- time_horizon
  ############################################################
  ## Create cost and utility models
  ############################################################
  
  # Utility and cost models
  # Utility
  utility_model <- create_StateVals(utility_tbl, n = n_samples, 
                                 time_reset = TRUE,
                                 hesim_data = hesim_dat) # Not totally sure this argument is needed
  
  # Costs
  # The 'starting' option means costs are only applied when patients enter the state
  implant_cost_model <- create_StateVals(implantcost_tbl, n = n_samples,
                                     method = "starting",
                                     hesim_data = hesim_dat) # Expand by patients
  
  medical_cost_model <- create_StateVals(medcost_tbl, n = n_samples,
                                 time_reset = TRUE,
                                 hesim_data = hesim_dat) 
  
  cost_models <- list(Drug = implant_cost_model,
                   Medical = medical_cost_model)
  
  return(list("utility_model" = utility_model, 
              "cost_models" = cost_models))
  
}

