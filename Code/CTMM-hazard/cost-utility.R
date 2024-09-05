generate_cost_utility_models <- function(n_samples, 
                                         hesim_dat,
                                         model_inputs) {
  
  # Extract requisite random and fixed parameters from model_inputs
  # This is to ensure input_parameters aligns with format required for BCEA-EVPPI
  hazard_rate_1st_revision <- read_excel((paste0(data_directory, "/", paste0(sample_gender, "-", age_range[1],"--",  "rate.xlsx"))),sheet = "Sheet1")
  lograte_revision <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "revision_log_rate")
  log_rate_1st_revision = hazard_rate_1st_revision
  log_rate_1st_revision[,2:10]= log(hazard_rate_1st_revision[,2:10])
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
  
  revision_cost <- model_inputs$revision_cost
  
  # Random parameters
  input_parameters <- model_inputs$inputs_parameters
  probability_higher_revision <- input_parameters[, "probability_higher_revision"]
  
  revision_disutility <- input_parameters[, "revision_disutility"]
  utility_post_tkr <- input_parameters[, "utility_post_tkr"] 
  utility_post_1st_rev <- input_parameters[, "utility_post_1st_rev"]
  utility_post_2nd_rev <- input_parameters[, "utility_post_2nd_rev"]
  
  
  implant_costs <- input_parameters[, grep("implant_cost_", colnames(input_parameters))]
  colnames(implant_costs) <- gsub("implant_cost_", "", colnames(implant_costs))
  implant_costs <- t(implant_costs)
  
  n_times <- 3
  
  
  prob1= data.frame(matrix(ncol = n_samples, nrow = length(implant_names)))
  prob2 = data.frame(implant_names)
  prob_3=cbind(prob2, prob1)
  prob_310=cbind(prob2, prob1)
  prob_10=cbind(prob2, prob1)
  
  for(i_strategy in 1:n_strategies){
    prob_3[i_strategy,2:(n_samples+1)]= (1 - exp(- exp(with(log_rate_1st_revision,
                                                            rnorm(n_samples,
                                                                  mean = mean1[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]],
                                                                  sd = ((UL1[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]]-LL1[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]])/2*1.96))))))
  }
  
  for(i_strategy in 1:n_strategies){
    prob_310[i_strategy,2:(n_samples+1)]= (1 - exp(- exp(with(log_rate_1st_revision,
                                                              rnorm(n_samples,
                                                                    mean = mean2[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]],
                                                                    sd = ((UL2[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]]-LL2[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]])/2*1.96))))))
  }
  
  for(i_strategy in 1:n_strategies){
    prob_10[i_strategy,2:(n_samples+1)]= (1 - exp(- exp(with(log_rate_1st_revision,
                                                             rnorm(n_samples,
                                                                   mean = mean3[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]],
                                                                   sd = ((UL3[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]]-LL3[treatment == strategies$strategy_name[which(strategies$strategy_id == i_strategy)]])/2*1.96))))))
  }
  
  if(is.infinite(age_range[2])) {
   for(i_strategy in 1:n_strategies){
    prob_10[i_strategy,2:(n_samples+1)]= rep(5.136108e-04, times = n_samples)
  }
  }
  
  utility_tbl <- stateval_tbl(
    data.table(strategy_id = rep(c(1:n_strategies), each = n_samples * n_states * n_times),
               sample = rep(rep(1:n_samples, each = n_states * n_times), times = n_strategies),
               state_id = rep(rep(c(1:n_states), each = n_times), times = n_strategies * n_samples),
               time_start = rep(c(0, 3, 10), times = n_strategies * n_samples * n_states),
               value = 0),
    dist = "custom"
  )
  
  
  revision_disutility <- input_parameters[, "revision_disutility"]
  multi_disutilities = rlnorm(n_samples, log(4.5), log(1.5))
  
  for(i_strategy in 1:n_strategies) { 
    utility_tbl[strategy_id==i_strategy & state_id == 1 & time_id == 1, "value"] <- utility_post_tkr + revision_disutility*prob_3[i_strategy,2:(n_samples+1)]
    utility_tbl[strategy_id==i_strategy & state_id == 1 & time_id == 2, "value"] <- utility_post_tkr + revision_disutility*prob_310[i_strategy,2:(n_samples+1)]
    utility_tbl[strategy_id==i_strategy & state_id == 1 & time_id == 3, "value"] <- utility_post_tkr + revision_disutility*prob_10[i_strategy,2:(n_samples+1)]}
  utility_tbl[state_id == 2, "value"] <- utility_post_tkr + revision_disutility*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                     mean = log_rate_2nd_revision_early_mean,
                                                                                                     sd = log_rate_2nd_revision_early_se))))
  utility_tbl[state_id == 3, "value"] <- utility_post_tkr + revision_disutility*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                     mean = log_rate_2nd_revision_middle_mean,
                                                                                                     sd = log_rate_2nd_revision_middle_se))))
  utility_tbl[state_id == 4, "value"] <- utility_post_tkr + revision_disutility* (1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                      mean = log_rate_2nd_revision_late_mean,
                                                                                                      sd = log_rate_2nd_revision_late_se))))
  
  # Utility post 2nd revision 
  # Need to subtract disutilities of higher revision multiplied by (annual) probability of having revision
  # Don't depend on time, time to second revision, or strategy
  utility_tbl[state_id == 5, "value"] <- utility_post_tkr + revision_disutility*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                                                     mean = lograte_higher_revision_mean,
                                                                                                     sd = lograte_higher_revision_se)))) 
 
  # Correct the utility table so that no utilities applied after the time_horizon
  utility_tbl$time_stop[is.infinite(utility_tbl$time_stop)] <- time_horizon
  utility_tbl$time_stop[utility_tbl$time_stop > time_horizon] <- time_horizon
  
  utility_tbl[state_id == 6, "value"] <- 0
  utility_tbl[state_id == 7, "value"] <- 0
  
  
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
               time_start = rep(c(0, 3, 10), times = n_strategies * n_samples * n_states),
               value = 0),
    dist = "custom"
  )
  
  
  
  for(i_strategy in 1:n_strategies) {
    medcost_tbl[strategy_id==i_strategy & state_id == 1 & time_id == 1, "value"] <- revision_cost*prob_3[i_strategy,2:(n_samples+1)]
    medcost_tbl[strategy_id==i_strategy & state_id == 1 & time_id == 2, "value"] <- revision_cost*prob_310[i_strategy,2:(n_samples+1)]
    medcost_tbl[strategy_id==i_strategy & state_id == 1 & time_id == 3, "value"] <- revision_cost*prob_10[i_strategy,2:(n_samples+1)]}
  
  
  medcost_tbl[state_id == 2, "value"] <- revision_cost*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                            mean = log_rate_2nd_revision_early_mean,
                                                                            sd = log_rate_2nd_revision_early_se))))
  medcost_tbl[state_id == 3, "value"] <- revision_cost*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                            mean = log_rate_2nd_revision_middle_mean,
                                                                            sd = log_rate_2nd_revision_middle_se))))
  medcost_tbl[state_id == 4, "value"] <- revision_cost* (1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                             mean = log_rate_2nd_revision_late_mean,
                                                                             sd = log_rate_2nd_revision_late_se))))
  medcost_tbl[state_id == 5, "value"] <- revision_cost*(1 - exp(- exp(rnorm(n = n_times *n_samples*n_strategies,
                                                                            mean = lograte_higher_revision_mean,
                                                                            sd = lograte_higher_revision_se))))
  
  
  
  medcost_tbl$time_stop[is.infinite(medcost_tbl$time_stop)] <- time_horizon
  medcost_tbl$time_stop[medcost_tbl$time_stop > time_horizon] <- time_horizon
 
  ############################################################
  ## Create cost and utility models
  ############################################################
  
  # Utility and cost models
  # Utility
  
  utility_model <- create_StateVals(utility_tbl, n = n_samples, 
                                    hesim_data = hesim_dat) # Not totally sure this argument is needed
  
  
  # Costs
  # The 'starting' option means costs are only applied when patients enter the state
  implant_cost_model <- create_StateVals(implantcost_tbl, n = n_samples,
                                         method = "starting",
                                         hesim_data = hesim_dat) # Expand by patients
  
  medical_cost_model <- create_StateVals(medcost_tbl, n = n_samples,
                                         hesim_data = hesim_dat) 
  
  cost_models <- list(Drug = implant_cost_model,
                      Medical = medical_cost_model)
  
  return(list("utility_model" = utility_model, 
              "cost_models" = cost_models))
  
}

