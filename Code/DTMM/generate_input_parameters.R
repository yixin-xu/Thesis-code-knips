generate_input_parameters <- function(n_samples, sensitivity = NULL) {
  
  
  parameter_names <- c(paste0("log_rate_1st_revision_<3", treatment_names),
                       paste0("log_rate_1st_revision_3-10", treatment_names),
                       paste0("log_rate_1st_revision_>10", treatment_names),
                       "log_rate_2nd_revision_early", "log_rate_2nd_revision_middle", 
                       "log_rate_2nd_revision_late", "log_rate_higher_revision",
                       "cost_revision",
                       paste0("implant_cost_", treatment_names),
                       paste0("qalys_", state_names),
                       "revision_disutility")
  n_parameters <- length(parameter_names)
  input_parameters <- array(dim = c(n_samples, n_parameters), dimnames = list(NULL, parameter_names))
  
  # data: female 0-55 years old group
  
  lograte_revision <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "revision_log_rate")
  implant_costs <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "implant_costs")
  other_costs <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "other_costs")
  utilities <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "utilities")
  
  if(is.infinite(finial_age)) {
    # Above age_range[1]
    log_rate_1st_revision = read_excel(paste0(data_directory,"/",gender,"-", initial_age, "-",  "rate.xlsx"))  
  }else{
    log_rate_1st_revision = read_excel(paste0(data_directory,"/",gender, "-", initial_age,"-",  "rate.xlsx"))
  }
  
  log_rate_1st_revision[,2:10]= log(log_rate_1st_revision[,2:10])
  
  # Impute log hazard rate for each implant for first revision less than 3 years
  if ("log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cem CR_Fix Mod"],
                                                                                                                    sd = ((UL1[treatment == "Cem CR_Fix Mod"]-LL1[treatment == "Cem CR_Fix Mod"])/2*1.96))) }     
  if ("log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mono" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mono"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cem CR_Fix Mono"],
                                                                                                                   sd = ((UL1[treatment == "Cem CR_Fix Mono"]-LL1[treatment == "Cem CR_Fix Mono"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem CR_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cem CR_Mob Mod"],
                                                                                                                   sd = ((UL1[treatment == "Cem CR_Mob Mod"]-LL1[treatment == "Cem CR_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cem PS_Fix Mod"],
                                                                                                                   sd = ((UL1[treatment == "Cem PS_Fix Mod"]-LL1[treatment == "Cem PS_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem PS_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cem PS_Mob Mod"],
                                                                                                                   sd = ((UL1[treatment == "Cem PS_Mob Mod"]-LL1[treatment == "Cem PS_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Cem Con_Con Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem Con_Con Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Cem Con_Con Mod"],
                                                                                                                    sd = ((UL1[treatment == "Cem Con_Con Mod"]-LL1[treatment == "Cem Con_Con Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Unc CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Unc CR_Fix Mod"],
                                                                                                                   sd = ((UL1[treatment == "Unc CR_Fix Mod"]-LL1[treatment == "Unc CR_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Unc CR_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Unc CR_Mob Mod"],
                                                                                                                   sd = ((UL1[treatment == "Unc CR_Mob Mod"]-LL1[treatment == "Unc CR_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Unc PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Unc PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Unc PS_Fix Mod"],
                                                                                                                   sd = ((UL1[treatment == "Unc PS_Fix Mod"]-LL1[treatment == "Unc PS_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Hyb CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "Hyb CR_Fix Mod"],
                                                                                                                   sd = ((UL1[treatment == "Hyb CR_Fix Mod"]-LL1[treatment == "Hyb CR_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_<3Implant OX Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant OX Cem CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "OX Cem CR_Fix Mod"],
                                                                                                                  sd = ((UL1[treatment == "OX Cem CR_Fix Mod"]-LL1[treatment == "OX Cem CR_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_<3Implant OX Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant OX Cem PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean1[treatment == "OX Cem PS_Fix Mod"],
                                                                                                                  sd = ((UL1[treatment == "OX Cem PS_Fix Mod"]-LL1[treatment == "OX Cem PS_Fix Mod"])/2*1.96)))}                                                                                       

  # Log hazard rate for first revision more than 3 years and less than 10 years
  if ("log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mono" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mono"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cem CR_Fix Mono"],
                                                                                                                      sd = ((UL2[treatment == "Cem CR_Fix Mono"]-LL2[treatment == "Cem CR_Fix Mono"])/2*1.96)))}
  if ("log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cem CR_Fix Mod"],
                                                                                                                     sd = ((UL2[treatment == "Cem CR_Fix Mod"]-LL2[treatment == "Cem CR_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_3-10Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem CR_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cem CR_Mob Mod"],
                                                                                                                     sd = ((UL2[treatment == "Cem CR_Mob Mod"]-LL2[treatment == "Cem CR_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cem PS_Fix Mod"],
                                                                                                                     sd = ((UL2[treatment == "Cem PS_Fix Mod"]-LL2[treatment == "Cem PS_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem PS_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cem PS_Mob Mod"],
                                                                                                                     sd = ((UL2[treatment == "Cem PS_Mob Mod"]-LL2[treatment == "Cem PS_Mob Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_3-10Implant MoP Cem Con_Con Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem Con_Con Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Cem Con_Con Mod"],
                                                                                                                      sd = ((UL2[treatment == "Cem Con_Con Mod"]-LL2[treatment == "Cem Con_Con Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_3-10Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Unc CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Unc CR_Fix Mod"],
                                                                                                                     sd = ((UL2[treatment == "Unc CR_Fix Mod"]-LL2[treatment == "Unc CR_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Unc CR_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Unc CR_Mob Mod"],
                                                                                                                     sd = ((UL2[treatment == "Unc CR_Mob Mod"]-LL2[treatment == "Unc CR_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Unc PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Unc PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Unc PS_Fix Mod"],
                                                                                                                     sd = ((UL2[treatment == "Unc PS_Fix Mod"]-LL2[treatment == "Unc PS_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Hyb CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "Hyb CR_Fix Mod"],
                                                                                                                     sd = ((UL2[treatment == "Hyb CR_Fix Mod"]-LL2[treatment == "Hyb CR_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_3-10Implant OX Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant OX Cem CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "OX Cem CR_Fix Mod"],
                                                                                                                    sd = ((UL2[treatment == "OX Cem CR_Fix Mod"]-LL2[treatment == "OX Cem CR_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_3-10Implant OX Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant OX Cem PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean2[treatment == "OX Cem PS_Fix Mod"],
                                                                                                                    sd = ((UL2[treatment == "OX Cem PS_Fix Mod"]-LL2[treatment == "OX Cem PS_Fix Mod"])/2*1.96)))}   

  # Log hazard rate for first revision more than 10 years 
  if ("log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mono" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mono"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cem CR_Fix Mono"],
                                                                                                                     sd = ((UL3[treatment == "Cem CR_Fix Mono"]-LL3[treatment == "Cem CR_Fix Mono"])/2*1.96)))}
  if ("log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cem CR_Fix Mod"],
                                                                                                                    sd = ((UL3[treatment == "Cem CR_Fix Mod"]-LL3[treatment == "Cem CR_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem CR_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cem CR_Mob Mod"],
                                                                                                                    sd = ((UL3[treatment == "Cem CR_Mob Mod"]-LL3[treatment == "Cem CR_Mob Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cem PS_Fix Mod"],
                                                                                                                    sd = ((UL3[treatment == "Cem PS_Fix Mod"]-LL3[treatment == "Cem PS_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem PS_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cem PS_Mob Mod"],
                                                                                                                    sd = ((UL3[treatment == "Cem PS_Mob Mod"]-LL3[treatment == "Cem PS_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_>10Implant MoP Cem Con_Con Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem Con_Con Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Cem Con_Con Mod"],
                                                                                                                     sd = ((UL3[treatment == "Cem Con_Con Mod"]-LL3[treatment == "Cem Con_Con Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Unc CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Unc CR_Fix Mod"],
                                                                                                                    sd = ((UL3[treatment == "Unc CR_Fix Mod"]-LL3[treatment == "Unc CR_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Unc CR_Mob Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Unc CR_Mob Mod"],
                                                                                                                    sd = ((UL3[treatment == "Unc CR_Mob Mod"]-LL3[treatment == "Unc CR_Mob Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_>10Implant MoP Unc PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Unc PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Unc PS_Fix Mod"],
                                                                                                                    sd = ((UL3[treatment == "Unc PS_Fix Mod"]-LL3[treatment == "Unc PS_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_>10Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Hyb CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "Hyb CR_Fix Mod"],
                                                                                                                    sd = ((UL3[treatment == "Hyb CR_Fix Mod"]-LL3[treatment == "Hyb CR_Fix Mod"])/2*1.96)))}                                                                                       
  if ("log_rate_1st_revision_>10Implant OX Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant OX Cem CR_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "OX Cem CR_Fix Mod"],
                                                                                                                   sd = ((UL3[treatment == "OX Cem CR_Fix Mod"]-LL3[treatment == "OX Cem CR_Fix Mod"])/2*1.96))) }                                                                                      
  if ("log_rate_1st_revision_>10Implant OX Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant OX Cem PS_Fix Mod"] <- with(log_rate_1st_revision, rnorm(n_samples, mean = mean3[treatment == "OX Cem PS_Fix Mod"],
                                                                                                                   sd = ((UL3[treatment == "OX Cem PS_Fix Mod"]-LL3[treatment == "OX Cem PS_Fix Mod"])/2*1.96))) }          
 # Re-revision rate 
  lograte_revision[,2:4]= log(lograte_revision[,2:4])
  
  input_parameters[ , "log_rate_2nd_revision_early"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "early_revision"],
                                                                                     sd = ((UL[parameter == "early_revision"]-LL[parameter == "early_revision"])/2*1.96)))
  input_parameters[ , "log_rate_2nd_revision_middle"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "middle_revision"],
                                                                                      sd = ((UL[parameter == "middle_revision"]-LL[parameter == "middle_revision"])/2*1.96)))
  input_parameters[ , "log_rate_2nd_revision_late"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "late_revision"],
                                                                                    sd = ((UL[parameter == "late_revision"]-LL[parameter == "late_revision"])/2*1.96)))
  # third revision
  sd_third = (as.numeric(lograte_revision[5,4]) - as.numeric(lograte_revision[5,3]))/2*1.96
  se_third = sd_third/sqrt(as.numeric(lograte_revision[5,5]))
  tau_third = 1/se_third^2
  
  # fourth revision
  sd_fourth = (as.numeric(lograte_revision[6,4]) - as.numeric(lograte_revision[6,3]))/2*1.96
  se_fourth = sd_fourth/sqrt(as.numeric(lograte_revision[6,5]))
  tau_fourth = 1/se_fourth^2
  
  # fifth revision
  sd_fifth = (as.numeric(lograte_revision[7,4]) - as.numeric(lograte_revision[7,3]))/2*1.96
  se_fifth = sd_fifth/sqrt(as.numeric(lograte_revision[7,5]))
  tau_fifth = 1/se_fifth^2
  
  # sixth revision
  sd_sixth = (as.numeric(lograte_revision[8,4]) - as.numeric(lograte_revision[8,3]))/2*1.96
  se_sixth = sd_sixth/sqrt(as.numeric(lograte_revision[8,5]))
  tau_sixth = 1/se_sixth^2
  
  # seventh revision
  sd_seventh = (as.numeric(lograte_revision[9,4]) - as.numeric(lograte_revision[9,3]))/2*1.96
  se_seventh = sd_seventh/sqrt(as.numeric(lograte_revision[9,5]))
  tau_seventh = 1/se_seventh^2
  
  
  # weight of ecah revision, using in the calculation of weighted average of log rates
  w_third = tau_third/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh )
  w_fourth = tau_fourth/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  w_fifth = tau_fifth/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  w_sixth = tau_sixth/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  w_seventh = tau_seventh/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  
  average_lograte = w_third*as.numeric(lograte_revision[5,3]) + w_fourth*as.numeric(lograte_revision[6,3]) + w_fifth*as.numeric(lograte_revision[7,3]) + 
    w_sixth*as.numeric(lograte_revision[8,3]) + w_seventh*as.numeric(lograte_revision[9,3])
  average_se = sqrt(w_third^2*se_third^2 + w_fourth^2*se_fourth^2 + w_fifth^2*se_fifth^2 + w_sixth^2*se_sixth^2 + w_seventh^2*se_seventh^2)
  
  input_parameters[ , "log_rate_higher_revision"] = rnorm(n_samples, average_lograte, average_se)
  
  ############################################################
  ## Costs and utilities
  ############################################################
  
  # Input utilities
  utility_row_index <- which(utilities[, "age_upper"] == finial_age &
                               utilities[, "gender"] == gender)
  
  revision_disutility_mean <- as.numeric(utilities[utility_row_index, "disutilities"])
  revision_disutility_se <-  as.numeric(utilities[utility_row_index, "disutilities SE"])
  utility_post_tkr_mean <-  as.numeric(utilities[utility_row_index, "6 months after primary"])
  utility_post_1st_rev_mean <-  as.numeric(utilities[utility_row_index, "6 months after revision"])
  utility_post_2nd_rev_mean <-   as.numeric(utilities[utility_row_index, "6 months after revision"])
  utility_post_tkr_se <-  as.numeric(utilities[utility_row_index, "6 months after primary SE"])
  utility_post_1st_rev_se <-  as.numeric(utilities[utility_row_index, "6 months after revision SE"])
  utility_post_2nd_rev_se <-  as.numeric(utilities[utility_row_index, "6 months after revision SE"])
  
  # Sample the utilities
  input_parameters[ ,"revision_disutility"] <- rnorm(n_samples, mean = revision_disutility_mean, sd = revision_disutility_se)
  
  input_parameters[ ,"qalys_State Post TKR <3 years"] <- input_parameters[ ,"qalys_State Post TKR >=3 years < 10 years"] <- input_parameters[ ,"qalys_State Post TKR >=10 years"] <- 
    rnorm(n_samples, mean = utility_post_tkr_mean, sd = utility_post_tkr_se)
  
  input_parameters[ ,"qalys_State Early revision"] <- input_parameters[ ,"qalys_State middle revision"] <-  input_parameters[ ,"qalys_State late revision"] <- 
    rnorm(n_samples, mean = utility_post_1st_rev_mean, sd = utility_post_1st_rev_se)
  
  input_parameters[ ,"qalys_State second revision"] <- rnorm(n_samples, mean = utility_post_2nd_rev_mean, sd = utility_post_2nd_rev_se)
  
  input_parameters[ ,"qalys_State Death"]<- rep(rep(0, each = n_samples), n= n_treatments)
  
  # Input costs for each implant
  if ("implant_cost_Implant MoP Cem CR_Fix Mono" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem CR_Fix Mono"] <- rep(as.numeric(implant_costs[1,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem CR_Fix Mod"] <- rep(as.numeric(implant_costs[2,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem CR_Mob Mod"] <-  rep(as.numeric(implant_costs[3,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem PS_Fix Mod"] <- rep(as.numeric(implant_costs[4,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem PS_Mob Mod"] <- rep(as.numeric(implant_costs[5,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Cem Con_Con Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem Con_Con Mod"] <- rep(as.numeric(implant_costs[6,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Unc CR_Fix Mod"] <- rep(as.numeric(implant_costs[7,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Unc CR_Mob Mod"] <- rep(as.numeric(implant_costs[8,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Unc PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Unc PS_Fix Mod"] <- rep(as.numeric(implant_costs[9,2]),times = n_samples) }
  if ("implant_cost_Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Hyb CR_Fix Mod"] <- rep(as.numeric(implant_costs[10,2]),times = n_samples) }
  if ("implant_cost_Implant OX Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant OX Cem CR_Fix Mod"] <- rep(as.numeric(implant_costs[11,2]),times = n_samples) }
  if ("implant_cost_Implant OX Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant OX Cem PS_Fix Mod"] <- rep(as.numeric(implant_costs[12,2]),times = n_samples) }
  
  
  input_parameters[ , "cost_revision"] <- rep(as.numeric(other_costs[2,2]), n = n_samples) 
  
  
  
  return(input_parameters)
} # End function

