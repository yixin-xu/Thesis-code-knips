# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to age-dependent transition matrices
# Output is a 3-dimensional array of n_samples x n_treatments x n_states x n_states


generate_transition_matrices <- function(input_parameters, n_cycles=50, sensitivity = NULL) {
  # create transition matrix 
  transition_matrices <- array(0, dim = c(n_cycles, n_treatments, n_samples, n_states, n_states),
                               dimnames = list(NULL, treatment_names, NULL, state_names, state_names))
  
  lifetime <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "uk_lifetables")
  mortality <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "90d_mortality")
  
  
  # Mortality probability 
  mortality_row_index <- which(mortality[, "age_upper"] == finial_age &
                                 mortality[, "gender"] == gender)
  
  # Mortality probability is sampled from a beta distribution
  # Each sample is converted to a rate
  # First mortality
  tkr_mortality_90d_beta_params <- as.numeric(mortality[mortality_row_index, "tkr_number"]) * c(
    as.numeric(mortality[mortality_row_index, "tkr_prop_dead"]),
    1 - as.numeric(mortality[mortality_row_index, "tkr_prop_dead"]))
  first_mortality <- rbeta(n_samples, shape1 = tkr_mortality_90d_beta_params[1], shape2 = tkr_mortality_90d_beta_params[2])
  
  # Revision mortality
  revision_mortality_90d_beta_params <- as.numeric(mortality[mortality_row_index, "revision_number"]) * c(
    as.numeric(mortality[mortality_row_index, "revision_prop_dead"]),
    1 - as.numeric(mortality[mortality_row_index, "revision_prop_dead"]))
  second_mortality <- rbeta(n_samples, shape1 = revision_mortality_90d_beta_params[1], shape2 = revision_mortality_90d_beta_params[2])
  
  # Second the higher revision mortality
  rerevision_mortality_90d_beta_params <- as.numeric(mortality[mortality_row_index, "rerevision_number"]) * c(
    as.numeric(mortality[mortality_row_index, "rerevision_prop_dead"]),
    1 - as.numeric(mortality[mortality_row_index, "rerevision_prop_dead"]))
  higher_mortality <- rbeta(n_samples, shape1 = rerevision_mortality_90d_beta_params[1], shape2 = rerevision_mortality_90d_beta_params[2])
  
  # Transition probabilities between states
  for(i_cycle in 1:n_cycles) {
    for(treatment_name in treatment_names) {
      transition_matrices[i_cycle, treatment_name, , "State Post TKR <3 years", "State Early revision"] <- 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))*(1-first_mortality)
      transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State middle revision"] <- 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))*(1-first_mortality)
      transition_matrices[i_cycle, treatment_name, , "State Post TKR >=10 years", "State late revision"] <- 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))*(1-first_mortality)
      
      
      transition_matrices[i_cycle, treatment_name, , "State Early revision", "State second revision"] <- 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))*(1-second_mortality)
      transition_matrices[i_cycle, treatment_name, , "State middle revision", "State second revision"] <- 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))*(1-second_mortality)
      transition_matrices[i_cycle, treatment_name, , "State late revision", "State second revision"] <- 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))*(1-second_mortality)
      
      # Annual probability of death the same for all states
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))*first_mortality
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))*first_mortality
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=10 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))*first_mortality
      
      transition_matrices[ i_cycle, treatment_name, , "State Early revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))*second_mortality
      transition_matrices[ i_cycle, treatment_name, , "State middle revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))*second_mortality
      transition_matrices[ i_cycle, treatment_name, , "State late revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))*second_mortality
      
      transition_matrices[ i_cycle, treatment_name, , "State second revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, 2])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))*higher_mortality
      
      # Tunnel states
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Post TKR >=3 years < 10 years"] = (1-transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Early revision"]-
                                                                                                                             transition_matrices[i_cycle, treatment_name, , "State Post TKR <3 years", "State Death"])/4
      
      transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Post TKR >=10 years"] = (1-transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State middle revision"]-
                                                                                                                              transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Death"])/8
      
      
      # Ensure remaining patients stay in state
      # Sum probabilities of transitions to other states
      # This ensures probabilities sum to 1.
      for(i_state in 1:length(state_names)) {
        transition_matrices[i_cycle, treatment_name, , i_state, i_state] <- 1 - 
          apply(transition_matrices[i_cycle, treatment_name, , i_state, -i_state], c(1), sum, na.rm=TRUE)
      } # End loop over states
      
      transition_matrices[is.na(transition_matrices)] <- 0
      
    } # End loop over implant_names
  } # End loop over cycles
  
  # At this point do a test and throw and error if it fails
  if(prod(rowSums(transition_matrices[sample(1:n_cycles, 1),
                                      sample(treatment_names, 1),
                                      sample(1:n_samples, 1), , ]) == 1) != 1) {
    stop("Rows must sum to 1!")
  }
  
  return(transition_matrices)
}
