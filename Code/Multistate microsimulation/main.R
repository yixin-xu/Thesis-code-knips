library(hesim)
library(data.table)
library(kableExtra)
library(flexsurv)
library(ggplot2)
library(MASS)
library(readxl)
library(BCEA)
library(dplyr)

data_directory <-  "Thesis code/Dummy data"
code_directory <- "Thesis code/Code/Multistate microsimulation"

# Load requisite functions
source(paste0(code_directory, "/knips_functions.R"))
source(paste0(code_directory, "/generate_model_inputs.R"))
source(paste0(code_directory, "/generate_cost_utility_models.R"))
source(paste0(code_directory, "/generate_model_outputs.R"))


lifetables <- read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "uk_lifetables")


############################################################
## Model specification #####################################
############################################################

n_samples <- 10
n_patients <- 10

discount_rate <- 0.035

# To specify <55 use c(0, 55) and to specify 85+ use c(85, Inf)
age_range <- c(0,55)
starting_age <- 53
final_age <- 100

# Specify the gender
sample_gender <- "female"

#sys.time <-system.time({


implant_costs_raw <- as.data.frame(read_excel(paste0(data_directory, "/KNIPS_Main_input_data.xlsx"), sheet = "implant_costs"))  

n_implants <- dim(implant_costs_raw)[1]
implant_names <- implant_costs_raw$implant_name


# Time horizon is based on final age and starting age. The latter is age/gender dependent
if(!is.infinite(age_range[2])){
  time_horizon <- final_age - starting_age
}

if(is.infinite(age_range[2])) {
  time_horizon = 15}

# Number of time intervals only used for costs and utilities
# as these are discretized
n_time_intervals <- time_horizon

# All state names
state_names <- c("Post TKR", "temp revision",
                 "Early revision",  "middle revision", "late revision", "second revision", "Death", "Surgery Death")

############################################################
## hesim specification #####################################
############################################################

# Treatment strategies - one for each implant
strategies <- data.table(
  strategy_id = 1:n_implants,
  strategy_name = implant_names 
)
n_strategies <- nrow(strategies)

# Patients  
# Randomly sample their ages between the lower and upper age range limits

if(!is.infinite(age_range[2])){
  age = rep(age_range[2] - 5, n_patients)}
if(is.infinite(age_range[2])){
  age = rep(90, n_patients)}
patients <- data.table(
  patient_id = 1:n_patients,
  gender = sample_gender,
  age = age)



# Non-death States - total hip replacement (THR) and revisions
states <- data.table( # Non-other cause death health states
  state_id = 1:6,
  state_name = state_names[1:6]#, "Surgery Death")
) 
n_states <- nrow(states)

# "hesim data"
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states)

# Nice labels - added the two death states manually
labs <- get_labels(hesim_dat) 
labs$state_id <- c(1:length(state_names))
names(labs$state_id) <- state_names


# Possible transitions
tmat <- rbind(c(NA, 1, NA, NA, NA, NA, 8, 13), # From "Post TKR"
              c(NA, NA, 2, 3, 4, NA, NA, NA), 
              c(NA, NA, NA, NA, NA, 5, 9, 14), 
              c(NA, NA, NA, NA, NA, 6, 10, 15),
              c(NA, NA, NA, NA, NA, 7, 11, 16),
              c(NA, NA, NA, NA, NA, NA, 12, 17),
              c(NA, NA, NA, NA, NA, NA, NA, NA),
              c(NA, NA, NA, NA, NA, NA, NA, NA)) 
colnames(tmat) <- rownames(tmat) <- names(labs$state_id)
############################################################
## Run model ###############################################
############################################################

model_inputs <- generate_model_inputs(n_samples, 
                                      age_range = age_range, 
                                      sample_gender = sample_gender)

system.time({
  model_outputs <- generate_model_outputs(hesim_dat = hesim_dat, 
                                          model_inputs = model_inputs)
})

############################################################
## Analyse results #########################################
############################################################
economic_model <- model_outputs$economic_model

reference_implant <- "Cem CR_Fix Mod"

# First use the hesim summary functions as a check for results below
ce_sim <- economic_model$summarize()
summary(ce_sim, labels = labs) %>%
  format()

# Calculate total costs, QALYs, and net benefit manually in format for BCEA
total_costs  <- total_qalys <- matrix(nrow = n_samples, ncol = n_strategies)
colnames(total_costs) <- colnames(total_qalys) <- implant_names
for(implant_name in implant_names) {
  total_costs[, implant_name] <-  as.matrix(ce_sim$costs[ce_sim$costs$category == "total" & 
                                                           ce_sim$costs$strategy_id == which(implant_names == implant_name), "costs"])
  total_qalys[, implant_name] <-  as.matrix(ce_sim$qalys[ce_sim$qalys$strategy_id == which(implant_names == implant_name), "qalys"]) 
}

knips_bcea <- bcea(e = total_qalys, 
                   c = total_costs, ref = which(implant_names == reference_implant), 
                   interventions = implant_names) 

summary(knips_bcea, wtp = 20000)

#plot the cost-effectiveness plane
ceplane.plot(knips_bcea, wtp = 20000, xlim = c(-10, 10), ylim = c(-5000, 5000))


results_table <- summarise_results(total_costs = total_costs,
                                   total_qalys = total_qalys,
                                   intervention_names = implant_names,
                                   reference_intervention = reference_implant)


write.csv(results_table, file = paste0("Thesis code/Results/Multistate microsimulation/"
                                       , sample_gender, "_", age_range[1], "-",age_range[2],"_s",n_samples, "_p", n_patients,".csv"))

# Plots a CEAC but only include interventions with a probability > 5% of being most cost-effective
# And which have >5% probability
lambdas <- c(1:50) * 1000
net_benefit <- array(NA, dim = c(n_samples, n_implants, length(lambdas)))
ceac <- matrix(NA, nrow = n_implants, ncol = length(lambdas))
# Use loops as not computationally intensive
for(i_lambda in 1:50) {
  net_benefit[, , i_lambda] <- total_qalys * lambdas[i_lambda] - total_costs
  which_max_net_benefit <- apply(net_benefit[, , i_lambda], c(1), which.max)
  for(i_implant in 1:n_implants) {
    ceac[i_implant, i_lambda] <- mean(which_max_net_benefit == i_implant)
  }
}

optimal_implants <- which(rowSums(ceac > 0.05) > 0)
knips_bcea_optimal <- bcea(e = total_qalys[, optimal_implants], 
                           c = total_costs[, optimal_implants], 
                           ref = which(implant_names[optimal_implants] == reference_implant), 
                           interventions = implant_names[optimal_implants]) 

# Plot the CEAC
knips_multi_ce <- multi.ce(knips_bcea_optimal)


jpeg(filename = paste0("Thesis code/Results/Multistate microsimulation/",
                       sample_gender, "_", age_range[1], "-",age_range[2],"_s",n_samples, "_p", n_patients,".jpeg"))
ceac.plot(knips_multi_ce, graph = "ggplot",
          line = list(colors = c(1:length(optimal_implants))),
          pos = c(0, 1.00))
dev.off()
