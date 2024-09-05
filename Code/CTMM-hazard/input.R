# Generate an inputs list for the economic model
# This includes a matrix of random sampled parameters for use in BCEA-EVPPI calculations
# Fixed values are not included in the input_parameters matrix but as separate elements in model_inputs

generate_model_inputs <- function(n_samples, age_range, sample_gender) {
  
  lograte_revision <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "revision_log_rate")
  hazard_rate_1st_revision <- read_excel((paste0(data_directory, "/", paste0(sample_gender, "-", age_range[1],"--",  "rate.xlsx"))),sheet = "Sheet1")
  log_rate_1st_revision = hazard_rate_1st_revision
  log_rate_1st_revision[,2:10]= log(hazard_rate_1st_revision[,2:10])
  lifetables <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "uk_lifetables")
  
  ############################################################
  ## Transition rate parameters ##############################
  ############################################################
  
  # Transition 1
  # From Post TKR to Post early revision
  log_rate_1st_revision_3 = list()
  for(implant_name in implant_names) {
    log_rate_1st_revision_3[[implant_name]] <- cbind(as.matrix(with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == implant_name],
                                                                                                 sd = ((UL1[treatment == implant_name]-LL1[treatment == implant_name])/2*1.96)))),-100)
    colnames(log_rate_1st_revision_3[[implant_name]]) <- paste0("log_rate_1st_revision_3_", c(1, 2)) }
  log_rate_1st_revision_3_times <- c(0, 3)
  
  # Transition 2 
  # From Post TKR to Post middle revision
  log_rate_1st_revision_310 = list()
  for(implant_name in implant_names) {
    log_rate_1st_revision_310[[implant_name]] <- cbind(-100, as.matrix(with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == implant_name],
                                                                                                      sd = ((UL2[treatment == implant_name]-LL2[treatment == implant_name])/2*1.96)))),-100)
    colnames(log_rate_1st_revision_310[[implant_name]]) <- paste0("log_rate_1st_revision_310_", c(1, 2, 3))}
  log_rate_1st_revision_310_times <- c(0, 3,10)
  
  
  # Transition 3
  # From Post TKR to Post late revision
  log_rate_1st_revision_10 = list()
  for(implant_name in implant_names) {
    log_rate_1st_revision_10[[implant_name]] <- cbind(-100, as.matrix(with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == implant_name],
                                                                                                     sd = ((UL3[treatment == implant_name]-LL3[treatment == implant_name])/2*1.96)))) )
    colnames(log_rate_1st_revision_10[[implant_name]]) <- paste0("log_rate_1st_revision_10_", c(1, 2))}
  log_rate_1st_revision_10_times <- c(0, 10)
  
  
  
  # Transition 4
  # From Post early revision to Post second revision
  lograte_revision[,2:4]= log(lograte_revision[,2:4])
  lograte_revision <- as.data.frame(lograte_revision)
  rownames(lograte_revision) <- lograte_revision[, "parameter"]
  log_rate_2nd_revision_early_mean <- as.numeric(lograte_revision["early_revision", "mean"]) 
  log_rate_2nd_revision_early_se = ((as.numeric(lograte_revision["early_revision", "UL"])-as.numeric(lograte_revision["early_revision", "LL"]))/2*1.96) 
  log_rate_early_second <- rnorm(n = n_samples, mean = log_rate_2nd_revision_early_mean, sd = log_rate_2nd_revision_early_se)
  
  # Transition 5
  # From Post middle revision to Post second revision
  log_rate_2nd_revision_middle_mean <- as.numeric(lograte_revision["middle_revision", "mean"]) 
  log_rate_2nd_revision_middle_se = ((as.numeric(lograte_revision["middle_revision", "UL"])-as.numeric(lograte_revision["middle_revision", "LL"]))/2*1.96) 
  log_rate_middle_second <- rnorm(n = n_samples, mean = log_rate_2nd_revision_middle_mean, sd = log_rate_2nd_revision_middle_se)
  
  # Transition 6
  # From Post late revision to Post second revision
  log_rate_2nd_revision_late_mean <- as.numeric(lograte_revision["late_revision", "mean"]) 
  log_rate_2nd_revision_late_se = ((as.numeric(lograte_revision["late_revision", "UL"])-as.numeric(lograte_revision["late_revision", "LL"]))/2*1.96) 
  log_rate_late_second <- rnorm(n = n_samples, mean = log_rate_2nd_revision_late_mean, sd = log_rate_2nd_revision_late_se)
  
  # Transition 7-11
  # Post TKR to death 
  # Currenly background only but need to include KNIPS data
  # Annualised mortality probability
  colnames(lifetables) = c( "age",  "males", "females")
  mr <- as.matrix(lifetables[starting_age:(starting_age + time_horizon), paste0(sample_gender, "s")])
  mr_times <- c(0:time_horizon)
  
  # Transition 12
  # Hazard rate of Post TKR to surgery death 
  # 90-day mortality from NJR
  mortality_90d_raw <- as.matrix(read_excel(paste0(data_directory,"/KNIPS_Main_input_data.xlsx"), sheet = "90d_mortality"))
  mortality_row_index <- which(mortality_90d_raw[, "age_upper"] == age_range[2] &
                                 mortality_90d_raw[, "gender"] == sample_gender)
  
  # Mortality probability is sampled from a beta distribution
  # Each sample is converted to a rate
  tkr_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "tkr_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "tkr_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "tkr_prop_dead"]))
  tkr_mortality_90d_probability <- rbeta(n_samples, shape1 = tkr_mortality_90d_beta_params[1], shape2 = tkr_mortality_90d_beta_params[2])
  tkr_mortality_90d <- as.matrix(-log(1 - tkr_mortality_90d_probability))
  
  # Rate is zero after 90 days
  #tkr_mortality_90d <- cbind(tkr_mortality_90d, 0)
  #colnames(tkr_mortality_90d) <- paste0("tkr_mortality_90d_", c(1, 2))
  #tkr_mortality_times_90d <- c(0, 90/365.25)
  
  # Transition 13-15
  # Hazard rate of Post 1st revision to death
  # 90-day mortality from NJR
  
  # Mortality probability is sampled from a beta distribution
  # Each sample is converted to a rate
  revision_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "revision_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "revision_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "revision_prop_dead"]))
  revision_mortality_90d_probability <- rbeta(n_samples, shape1 = revision_mortality_90d_beta_params[1], shape2 = revision_mortality_90d_beta_params[2])
  revision_mortality_90d <- as.matrix(-log(1 - revision_mortality_90d_probability))
  
  # Rate is zero after 90 days
  #revision_mortality_90d <- cbind(revision_mortality_90d, 0)
  #colnames(revision_mortality_90d) <- paste0("revision_mortality_90d_", c(1, 2))
  #revision_mortality_times_90d <- c(0, 90/365.25)
  
  # Rates of 3rd or higher revision
  # Used for surgery mortality after 2nd revision (transition 22) and in costs and utilities
  # Inverse variance meta-analysis of logrates of 3rd to 8th revision from NJR
  other_rates_raw <- as.matrix(read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "other_rates"))
  rownames(other_rates_raw) <- other_rates_raw[, "parameter"]
  
  lograte_higher_revision_mean <-  as.numeric(other_rates_raw["lograte_higher_revision_mean", "value"])
  lograte_higher_revision_se <-  as.numeric(other_rates_raw["lograte_higher_revision_se", "value"])
  log_rate_higher <- rnorm(n = n_samples, mean = lograte_higher_revision_mean, sd = lograte_higher_revision_se)
  # Assumed not to depend on time, time to 2nd revision, or implant
  # Need annualised probabilities for multiplying by 90d mortality and in costs and utilities in "Post 2nd revision"
  probability_higher_revision <- 1 - exp(- exp(rnorm(n = n_samples,
                                                     mean = lograte_higher_revision_mean,
                                                     sd = lograte_higher_revision_se)))
  
  # Transition 16
  # Hazard rate of Post 2nd revision to death
  rerevision_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "rerevision_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "rerevision_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "rerevision_prop_dead"]))
  rerevision_mortality_90d_probability <- rbeta(n_samples, shape1 = rerevision_mortality_90d_beta_params[1], shape2 = rerevision_mortality_90d_beta_params[2])
  rerevision_mortality_90d <- as.matrix(-log(1 - rerevision_mortality_90d_probability))
  
  # Beyond initial 90d period there is an ongoing risk of higher revision and surgical mortality
  # Probability of dying (within 90 days) times annual probability higher revision
  higher_revision_annual_mortality_probability <- rerevision_mortality_90d_probability * probability_higher_revision
  higher_revision_mortality <- as.matrix(-log(1 - higher_revision_annual_mortality_probability))
  
  # Initial 90d risk from 2nd revision, ongoing risk from further revision
  #rerevision_mortality_90d <- cbind(rerevision_mortality_90d, higher_revision_mortality)
  #colnames(rerevision_mortality_90d) <- paste0("rerevision_mortality_90d_", c(1, 2))
  #rerevision_mortality_times_90d <- c(0, 90/365.25)
  
  
  ############################################################
  ## Costs and utilities
  ############################################################
  utilities_raw <- as.matrix(read_excel(paste0(data_directory,"/KNIPS_Main_input_data.xlsx"), sheet = "utilities"))
  utility_row_index <- which(utilities_raw[, "age_upper"] == age_range[2] &
                               utilities_raw[, "gender"] == sample_gender)
  
  revision_disutility_mean <- as.numeric(utilities_raw[utility_row_index, "disutilities"])
  revision_disutility_se <-  as.numeric(utilities_raw[utility_row_index, "disutilities SE"])
  utility_post_tkr_mean <-  as.numeric(utilities_raw[utility_row_index, "6 months after primary"])
  utility_post_1st_rev_mean <-  as.numeric(utilities_raw[utility_row_index, "6 months after revision"])
  utility_post_2nd_rev_mean <-   as.numeric(utilities_raw[utility_row_index, "6 months after revision"])
  utility_post_tkr_se <-  as.numeric(utilities_raw[utility_row_index, "6 months after primary SE"])
  utility_post_1st_rev_se <-  as.numeric(utilities_raw[utility_row_index, "6 months after revision SE"])
  utility_post_2nd_rev_se <-  as.numeric(utilities_raw[utility_row_index, "6 months after revision SE"])
  
  # Sample the utilities
  revision_disutility <- rnorm(n_samples, mean = revision_disutility_mean, sd = revision_disutility_se)
  utility_post_tkr <- rnorm(n_samples, mean = utility_post_tkr_mean, sd = utility_post_tkr_se)
  utility_post_1st_rev <- rnorm(n_samples, mean = utility_post_1st_rev_mean, sd = utility_post_1st_rev_se)
  utility_post_2nd_rev <- rnorm(n_samples, mean = utility_post_2nd_rev_mean, sd = utility_post_2nd_rev_se)
  
  # Implant costs 
  implant_costs_raw <- as.data.frame(read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "implant_costs"))
  
  implant_costs <- array(dim = c(n_strategies, n_samples), dimnames = list(implant_names, NULL))
  for(implant_name in implant_names) {
    implant_costs[implant_name, ] <- rep(implant_costs_raw$mean[implant_costs_raw$implant_name == implant_name],times = n_samples)
  }
  
  # Load other costs
  other_costs_raw <- as.data.frame(read_excel(paste0(data_directory,"/KNIPS_Main_input_data.xlsx"), sheet = "other_costs"))
  # Cost of revision surgery
  revision_cost <- other_costs_raw$value[other_costs_raw$parameter == "revision_cost"]
  
  ############################################################
  ## Return the inputs #######################################
  ############################################################
  
  # List of inputs for the economic model
  # This is a matrix of sample random parameters and some fixed parameters (e.g., knot locations, background mortality)
  model_inputs <- list()
  input_parameters <- matrix()
  
  # Fixed parameters
  model_inputs$mr <- mr
  #model_inputs$tkr_mortality_times_90d <- tkr_mortality_times_90d
  #model_inputs$revision_mortality_times_90d <- revision_mortality_times_90d
  #model_inputs$rerevision_mortality_times_90d <- rerevision_mortality_times_90d
  model_inputs$revision_cost <- revision_cost
  model_inputs$log_rate_1st_revision_3_times <- log_rate_1st_revision_3_times
  model_inputs$log_rate_1st_revision_310_times <- log_rate_1st_revision_310_times
  model_inputs$log_rate_1st_revision_10_times <- log_rate_1st_revision_10_times
  
  # Random parameters
  input_parameters <- log_rate_early_second
  input_parameters <- cbind(input_parameters, log_rate_middle_second, log_rate_late_second, log_rate_higher)
  input_parameters <- cbind(input_parameters, probability_higher_revision)
  colnames(input_parameters) <- c("log_rate_early_second", "log_rate_middle_second", "log_rate_late_second", "log_rate_higher","probability_higher_revision")
  for(implant_name in implant_names) {
    temp3 <- log_rate_1st_revision_3[[implant_name]]
    colnames(temp3) <- paste0("log_rate_1st_revision_3_", implant_name, "_", colnames(log_rate_1st_revision_3[[implant_name]]))
    input_parameters <- cbind(input_parameters, temp3)
    temp310 <- log_rate_1st_revision_310[[implant_name]]
    colnames(temp310) <- paste0("log_rate_1st_revision_310_", implant_name, "_", colnames(log_rate_1st_revision_310[[implant_name]]))
    input_parameters <- cbind(input_parameters, temp310)
    temp10 <- log_rate_1st_revision_10[[implant_name]]
    colnames(temp10) <- paste0("log_rate_1st_revision_10_", implant_name, "_", colnames(log_rate_1st_revision_10[[implant_name]]))
    input_parameters <- cbind(input_parameters, temp10)
  }
  
  
  input_parameters <- cbind(input_parameters, tkr_mortality_90d)
  colnames(input_parameters)[ncol(input_parameters)] <- "tkr_mortality_90d"
  
  input_parameters <- cbind(input_parameters, tkr_mortality_90d_probability)
  colnames(input_parameters)[ncol(input_parameters)] <- "tkr_mortality_90d_probability"
  
  input_parameters <- cbind(input_parameters, revision_mortality_90d)
  colnames(input_parameters)[ncol(input_parameters)] <- "revision_mortality_90d"
  
  input_parameters <- cbind(input_parameters, revision_mortality_90d_probability)
  colnames(input_parameters)[ncol(input_parameters)] <- "revision_mortality_90d_probability"
  
  input_parameters <- cbind(input_parameters, rerevision_mortality_90d)
  colnames(input_parameters)[ncol(input_parameters)] <- "rerevision_mortality_90d"
  
  input_parameters <- cbind(input_parameters, rerevision_mortality_90d_probability)
  colnames(input_parameters)[ncol(input_parameters)] <- "rerevision_mortality_90d_probability"
  
  input_parameters <- cbind(input_parameters, revision_disutility)
  input_parameters <- cbind(input_parameters, utility_post_tkr)
  input_parameters <- cbind(input_parameters, utility_post_1st_rev)
  input_parameters <- cbind(input_parameters, utility_post_2nd_rev)
  
  temp <- t(implant_costs)
  colnames(temp) <- paste0("implant_cost_", colnames(t(implant_costs)))
  input_parameters <- cbind(input_parameters, temp)
  
  # Add random parameters to the list
  model_inputs$inputs_parameters <- input_parameters
  
  return(model_inputs)
}

