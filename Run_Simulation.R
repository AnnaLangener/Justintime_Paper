####################################
######### Load Packages ############
####################################

library(randomForest)
library(foreach)
library(doParallel)
library(dplyr)
library(lme4)
library(caret)
library(pROC)

source("Simulation_Functions_notshiny.R")

# Basic Parameters from Study
n_features <- 10  # Features
n_samples <- 30  # Timepoints
n_subjects <- 10 # Subjects

# Generating Outcome
overall_prob_outcome <- 0.2  # Overall probability of 1 (for the entire dataset) (e.g., tracking depression 0.146)
sd_outcome <- 0.2 # Controls variability BETWEEN different subject (nneds to be smaller than sqrt(overall_prob_outcome*(1-overall_prob_outcome)))

# Generating Features
A <- 0.05 # Relationship between features and outcome
feature_std <- 0.1 # population level feature generating process
B <- 0.8  # Cross-subject variability ("random effect") (added per participants for all timepoints)
C <- 0.2 # Within-subject variability (added within participant for each timepoint)

time_effect = 0
timevariability = 0


features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B, C, overall_prob_outcome, sd_outcome, time_effect,timevariability)
features_sample_Astd <- features_sample[[2]]
features_sample <- features_sample[[1]]

test = features_sample  %>%
  group_by(subject) %>%
  summarize(mean =mean(prob))

mean(test$mean)
sd(test$mean)
####################################################
################ Run Simulation ####################
####################################################
source("Simulation_Functions_notshiny.R")

# n_features <- 10 
# n_samples <- 90  
# n_subjects <- 150 
# time_effect = FALSE
# 
# A <- 0 # Relationship between features and outcome, we also ran this for A = 0.05 and A = 0.2, and then merged the results
# feature_std <- seq(0.1,0.5, by = 0.15) # population level feature generating process
# 
# # Initialize parameters
# overall_prob_outcome <- seq(0.1, 0.9, by = 0.2)  
# sd_outcome <- seq(0.05, 0.29, by = 0.1)            
# B <- seq(0.1, 0.9, by = 0.1)       
# C <-  seq(0.05,0.45, by = 0.1)   
# 
# time_effect <- c(0.1,0.2,0.4,0.6,0.8)
# timevariability <- c(2,4,8,10)

### Old

n_features <- 10 
n_samples <- 200  
n_subjects <- 150 
time_effect = FALSE

A <- 0.2 # Relationship between features and outcome, we also ran this for A = 0.05 and A = 0.2, and then merged the results
feature_std <- 0.1 # population level feature generating process

# Initialize parameters
overall_prob_outcome <- seq(0.1, 0.9, by = 0.1)  
sd_outcome <- seq(0.05, 0.29, by = 0.05)            
B <- seq(0.1, 0.9, by = 0.05)       
C <-  seq(0.05,0.45, by = 0.05)   

time_effect <- 0
timevariability <- 0


# Create a grid of all parameter combinations
param_grid <- expand.grid(overall_prob_outcome = overall_prob_outcome,
                          sd_outcome = sd_outcome,
                          B = B,
                          C = C,
                          time_effect = time_effect,
                          timevariability = timevariability,
                          feature_std = feature_std)

# Set up parallel backend
n_cores <- parallel::detectCores() - 1  
cl <- makeCluster(n_cores)
registerDoParallel(cl)
old <- Sys.time() # get start time

# Parallel loop with foreach
result <- foreach(i = 0:nrow(param_grid), .combine = rbind, .packages = c("dplyr", "lme4", "pROC", "randomForest"), .errorhandling = "pass") %dopar% {
  # Extract parameter combination for this iteration
  params <- param_grid[i, ]
  overall_prob <- params$overall_prob_outcome
  sd <- params$sd_outcome
  B_value <- params$B
  C_value <- params$C
  time_effect <- params$time_effect
  timevariability <- params$timevariability
  feature_std <- params$feature_std
  
  # Generate data for this iteration
  features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B_value, C_value, overall_prob, sd, time_effect,timevariability)
  features_sample_Astd <- features_sample[[2]]
  features_sample <- features_sample[[1]]
  
  ### Descriptives ####
  test <- features_sample %>% group_by(subject) %>% summarise(prob = mean(y))
  
  # ICC outcome
  model <- lmer(y ~ 1 + (1 | subject), data = features_sample)
  var_random <- as.data.frame(VarCorr(model))$vcov[1]  # Random intercept variance
  var_residual <- attr(VarCorr(model), "sc")^2         # Residual variance
  icc = var_random / (var_random + var_residual)
  
  # ICC predictor
  model_pred <- lmer(V1 ~ 1 + ( 1  | subject), data = features_sample) # In our example we simulated the same relationships for all features
  var_random <- as.data.frame(VarCorr(model_pred))$vcov[1]  
  var_residual <- attr(VarCorr(model_pred), "sc")^2        
  icc_pred <- var_random / (var_random + var_residual)
  
  sim <- run_simulation(features_sample, "row-wise", 1, testsize = 0.3)
  sim3 <- run_simulation_centering(features_sample,"row-wise",1,testsize = 0.3)
  
  # Create results
  result_row <- data.frame(
    n_features = n_features,
    n_samples = n_samples,
    n_subjects = n_subjects,
    feature_std = feature_std,
    A = A,
    B = B_value,
    C = C_value,
    feature_std = feature_std,
    time_effect = time_effect,
    timevariability = timevariability,
    overall_prob_outcome = overall_prob,
    sd_outcome = sd,
    mean = mean(features_sample$y),
    sd = sd(test$prob),
    icc = icc,
    icc_pred = icc_pred,
    auc_value_base = sim$auc_value_base,
    auc = mean(sim$auc_value),
    auc_individual = sim$overall_summary$mean,
    auc_individual_sd = sim$overall_summary$sd,
    total_n = sim$overall_summary$total_n,
    auc_c = mean(sim3$auc_value),
    auc_c_individual = sim3$overall_summary$mean,
    auc_c_individual_sd = sim3$overall_summary$sd,
    sd_intercept = as.data.frame(VarCorr(model_pred))$sdcor[1], # B recovered
    sd_residual = as.data.frame(VarCorr(model_pred))$sdcor[2] # C recovered
  )
  
  # Save the result to a file in real time
  write.table(result_row, file = "simulation_results_02.csv", append = TRUE, sep = ",", col.names = !file.exists("simulation_results_02.csv"), row.names = FALSE)
  
  # Return the result 
  result_row
}


# Stop the cluster
stopCluster(cl)

# print elapsed time
Sys.time() - old # calculate difference



####################################################
############### Baseline model #####################
####################################################

source("Simulation_Functions_notshiny.R")

n_features <- 10 
n_samples <- 90 
n_subjects <- 150
time_effect = FALSE

A <- 0 
feature_std <- 0.1 

# Initialize parameters
overall_prob_outcome <- seq(0.1, 0.9, by = 0.01)  
sd_outcome <- seq(0.05, 0.29, by = 0.01) 
B <- 0.1
C <- 0.1

# Create a grid of all parameter combinations
param_grid <- expand.grid(overall_prob_outcome = overall_prob_outcome,
                          sd_outcome = sd_outcome,
                          n_samples = n_samples,
                          n_subjects = n_subjects)

# Set up parallel backend
n_cores <- parallel::detectCores() - 1  # Use all but two core
cl <- makeCluster(n_cores)
registerDoParallel(cl)
old <- Sys.time() # get start time

# Parallel loop with foreach
result <- foreach(i = 1:nrow(param_grid), .combine = rbind, .packages = c("dplyr", "lme4", "pROC", "randomForest"), .errorhandling = "remove") %dopar% {
  # Extract parameter combination for this iteration
  params <- param_grid[i, ]
  overall_prob <- params$overall_prob_outcome
  sd <- params$sd_outcome
  n_samples <- params$n_samples
  n_subjects <- params$n_subjects
  
  # Generate data for this iteration
  features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B, C, overall_prob, sd, time_effect)
  features_sample_Astd <- features_sample[[2]]
  features_sample <- features_sample[[1]]

  
  ### Descriptives ####
  test <- features_sample %>% group_by(subject) %>% summarise(prob = mean(y))
  
  # ICC outcome
  model <- lmer(y ~ 1 + (1 | subject), data = features_sample)
  var_random <- as.data.frame(VarCorr(model))$vcov[1]  # Random intercept variance
  var_residual <- attr(VarCorr(model), "sc")^2         # Residual variance
  icc = var_random / (var_random + var_residual)
  
  # ICC predictor
  model_pred <- lmer(V1 ~ 1 + ( 1  | subject), data = features_sample) # In our example we simulated the same relationships for all features
  var_random <- as.data.frame(VarCorr(model_pred))$vcov[1]  
  var_residual <- attr(VarCorr(model_pred), "sc")^2        
  icc_pred <- var_random / (var_random + var_residual)
  
  sim <- run_simulation(features_sample, "row-wise", 1, testsize = 0.3)
  
  # Create results
  result_row <- data.frame(
    n_features = n_features,
    n_samples = n_samples,
    n_subjects = n_subjects,
    feature_std = feature_std,
    A = A,
    B = B,
    C = C,
    overall_prob_outcome = overall_prob,
    sd_outcome = sd,
    time_effect = FALSE,
    mean = mean(features_sample$y),
    sd = sd(test$prob),
    icc = icc,
    icc_pred = icc_pred,
    auc_value_base = sim$auc_value_base,
    auc = mean(sim$auc_value),
    auc_individual = sim$overall_summary$mean,
    auc_individual_sd = sim$overall_summary$sd,
    total_n = sim$overall_summary$total_n,
    sd_intercept = as.data.frame(VarCorr(model_pred))$sdcor[1], # B recovered
    sd_residual = as.data.frame(VarCorr(model_pred))$sdcor[2] # C recovered
  )
  
  # Save the result to a file in real time
  write.table(result_row, file = "simulation_results_baseline.csv", append = TRUE, sep = ",", col.names = !file.exists("simulation_results_baseline.csv"), row.names = FALSE)
  
  # Return the result 
  result_row
}


# Stop the cluster
stopCluster(cl)
