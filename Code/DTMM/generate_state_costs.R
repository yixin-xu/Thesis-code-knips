# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to state costs


generate_state_costs<-function(input_parameters,sensitivity=NULL) {
  
  
  # Construct array of state costs with default value of zero
  state_costs <- array(0, dim = c(n_samples, n_treatments, n_states),
                       dimnames = list( NULL, treatment_names, state_names))
  
  # State costs are a mixture of treatment costs, acute event costs, event managment costs, 
  # and costs of transient events
  
  for(treatment_name in treatment_names) {
    state_costs[ ,treatment_name , "State Post TKR <3 years"] <- input_parameters[ , "cost_revision"] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))
    state_costs[ ,treatment_name , "State Post TKR >=3 years < 10 years"] <- input_parameters[ , "cost_revision"] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))
    state_costs[ ,treatment_name , "State Post TKR >=10 years"] <- input_parameters[ , "cost_revision"] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))
    
  }
  
  state_costs[ , , "State Early revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))
  state_costs[ , , "State middle revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))
  state_costs[ , , "State late revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))
  state_costs[ , , "State second revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))
  
  
  # End loop over treatments
  
  return(state_costs)
} # End function
