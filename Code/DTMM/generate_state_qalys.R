# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters state utilities

generate_state_qalys <- function(input_parameters, sensitivity = NULL) {
  

  # Construct array of state utilities with default value 0
  state_utilities <- array(dim = c(n_samples, n_treatments, n_states), dimnames = list(NULL, treatment_names,state_names))
  event_disutilities <- array(dim = c(n_samples, n_treatments, n_states), dimnames = list(NULL, treatment_names,state_names))
 
   # State utilities are a mixture of state event  utilities, and disutilities of transient events
  
  for(treatment_name in treatment_names) {
    event_disutilities[ ,treatment_name , "State Post TKR <3 years"] <- input_parameters[ ,"revision_disutility"]  * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))
    event_disutilities[ ,treatment_name , "State Post TKR >=3 years < 10 years"] <- input_parameters[ ,"revision_disutility"]  * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))
    event_disutilities[ ,treatment_name , "State Post TKR >=10 years"] <- input_parameters[ ,"revision_disutility"]  * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))
    
  }
  
  event_disutilities[ , , "State Early revision"] <- input_parameters[ ,"revision_disutility"]  * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))
  event_disutilities[ , , "State middle revision"] <- input_parameters[ ,"revision_disutility"]  * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))
  event_disutilities[ , , "State late revision"] <- input_parameters[ ,"revision_disutility"]  * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))
  event_disutilities[ , , "State second revision"] <- input_parameters[ ,"revision_disutility"]  * (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))
  
  
  state_utilities[, , "State Post TKR <3 years"] = rep(input_parameters[ ,"qalys_State Post TKR <3 years"], n= n_treatments) + event_disutilities[, , "State Post TKR <3 years"]
  state_utilities[, , "State Post TKR >=3 years < 10 years"] = rep(input_parameters[ ,"qalys_State Post TKR >=3 years < 10 years"], n= n_treatments) + event_disutilities[, , "State Post TKR >=3 years < 10 years"]
  state_utilities[, , "State Post TKR >=10 years"] =rep(input_parameters[ ,"qalys_State Post TKR >=10 years"], n= n_treatments) + event_disutilities[, , "State Post TKR >=10 years"]
  state_utilities[, , "State Early revision"] = rep(input_parameters[ ,"qalys_State Early revision"], n= n_treatments) + event_disutilities[, , "State Early revision"]
  state_utilities[, , "State middle revision"] = rep(input_parameters[ ,"qalys_State middle revision"], n= n_treatments) + event_disutilities[, , "State middle revision"]
  state_utilities[, , "State late revision"] = rep(input_parameters[ ,"qalys_State late revision"], n= n_treatments) + event_disutilities[, , "State late revision"]
  state_utilities[, , "State second revision"] = rep(input_parameters[ ,"qalys_State second revision"], n= n_treatments) + event_disutilities[, , "State second revision"]
  state_utilities[, , "State Death"] = rep(input_parameters[ ,"qalys_State Death"], n= n_treatments)
  
  
  return(state_utilities)
} # End function
