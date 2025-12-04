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
library(CorBin)
library(tidyr)



source("Simulation_continous.R")
source("Simulation_Functions_notshiny.R")


# Basic Parameters from Study
n_features <- 10  # Features
n_samples <- 90  # Timepoints
n_subjects <- 150 # Subjects

# Generating Outcome
overall_prob_outcome <- 0.5  # 0.5 | Overall probability of 1 (for the entire dataset) (e.g., tracking depression 0.146)
sd_outcome <- 0.25 # 0| Controls variability BETWEEN different subject (needs to be smaller than sqrt(overall_prob_outcome*(1-overall_prob_outcome)))

# Generating Features
A <- 0.1 # Relationship between features and outcome
feature_std <- 0.1 # population level feature generating process
B <- 0.8 #0.9  # Cross-subject variability ("random effect") (added per participants for all timepoints)
C <- 0.45 #0.45 # Within-subject variability (added within participant for each timepoint)

time_type = "DCP" #DCP, 1-dependent
time_effect = 0.2
time_effect =  seq(0.9, 0.01, length.out = (n_samples-1)) #0.2

time_effect = c(rep(0.5,n_samples-1))

#time_effect = c(rep(0.9,45),rep(0,44))


time_effect = NULL
# 
# lags <- c(1:(n_samples-1))
# values <- c(0.1, 0.05, 0.05, 0.05, rep(0.01,n_samples - 5))
# 
# 
# time_effect <- lapply(seq_along(lags), function(i) {
#   rep(c(values[i]), n_samples - lags[i])
# })



features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B, C, overall_prob_outcome, sd_outcome, time_effect, time_type)
features_sample_Astd <- features_sample[[2]]
features_sample <- features_sample[[1]]

test = run_simulation(features_sample, "subject-wise", 1, testsize = 0.3)

sub_id = unique(features_sample$subject)[1]

auc_value_ind = c()
auc_value_base = c()
for(k in 1:length(sub_id))
if(length(unique(features_sample$y[features_sample$subject == sub_id])) > 1){
  res = run_simulation_slidingwindow(features_sample[features_sample$subject == sub_id,],1,windowsize = 14)
  auc_value_ind[i] = res$auc_value
  auc_value_base[i] = res$auc_value_base
}


test = run_simulation(features_sample, "subject-wise", 1, testsize = 0.3)

mean(test$auc_value)

auc_value <- list()
for(i in 1:100000){
features_sample_shuffle <- shuffle_data(features_sample,subject_var = "subject", outcome_var = "y", n_features_upload = c("V1","V2","V3","V4", "V5","V6","V7","V8","V9", "V10"), seed = sample(1:1000000000,1))
test_shuffle = run_simulation(features_sample_shuffle, "subject-wise", 1, testsize = 0.3,seed = sample(1:1000000000,1))
auc_value[i] = test_shuffle$auc_value
}

mean(unlist(auc_value))

sum(auc_value > 0.5)/length(auc_value)

#####

test_row = run_simulation(features_sample, "row-wise", 10, testsize = 0.3)
mean(test$auc_value)


features_sample_shuffle <- shuffle_data(features_sample,subject_var = "subject", outcome_var = "y", n_features_upload = c("V1","V2","V3","V4", "V5","V6","V7","V8","V9", "V10"))
test_shuffle_row = run_simulation(features_sample_shuffle, "row-wise", 10, testsize = 0.3)


mean(test$auc_value)
mean(test_shuffle$auc_value)
mean(test_row$auc_value)
mean(test_shuffle_row$auc_value)

####



test = features_sample  %>%
  group_by(subject) %>%
  summarize(mean =mean(prob))

mean(test$mean)
sd(test$mean)

library(nlme)

model <-lme(
  y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10,                 
  random = ~1 | subject, # random intercept per participant
  correlation = corAR1(form = ~ time | subject),  # AR(1) within participant
  data = features_sample
)


model <-glm(
  y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10,                
  data = features_sample,
)


model <-glmer(
  y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + (1 | subject),                
  data = features_sample,
  family = "binomial"
)

performance::r2(model)

# r2/1-r2

# library(performance)
# performance::r2(model)
# 
# 0.43/(1-0.43)


##

model <-lme(
  y ~ time,                  # fixed effects
  random = ~1 | subject, # random intercept per participant
  correlation = corAR1(form = ~ time | subject),  # AR(1) within participant
  data = features_sample
)


phi <- coef(model$modelStruct$corStruct, unconstrained = FALSE)
phi



library(geepack)
model <- geeglm(
  y ~ time,
  id = subject,
  corstr = "ar1",
  family = binomial,
  data = features_sample
)

summary(mapVarsummary(model))

# Compute autocorrelation function for each participant
acf_df <- features_sample %>%
  group_by(subject) %>%
  summarise(acf_vals = list(acf(y, plot=FALSE)$acf),
            lag_vals = list(acf(y, plot=FALSE)$lag)) %>%
  unnest(c(lag_vals, acf_vals))

# Plot autocorrelation per participant
ggplot(acf_df[acf_df$subject %in% sample(unique(acf_df$subject),5),], aes(x=lag_vals, y=acf_vals, color=factor(subject))) +
  geom_line() +
  geom_point() +
  labs(x="Lag", y="Autocorrelation", color="Participant") +
  theme_minimal()


###
ggplot(features_sample[features_sample$subject %in% sample(unique(features_sample$subject),5),], aes(x = time, y = y)) +
  geom_line() +
  geom_point() +  
  labs(
    x = "Time",
    y = "Y Variable",
    color = "Participant"
  ) +
  facet_grid(~subject) +
  theme_minimal()

acf(features_sample$y)
pacf(features_sample$y)  


####################################################
################ Signal to Noise ###################
####################################################
source("Simulation_Functions_notshiny.R")


n_features <- 10
n_samples <- 200
n_subjects <- 150
time_effect = FALSE

A <- 0.05 # Relationship between features and outcome, we also ran this for A = 0.05 and A = 0.2, and then merged the results
feature_std <- 0.1 # population level feature generating process

# Initialize parameters
overall_prob_outcome <- seq(0.1, 0.9, by = 0.1)
sd_outcome <- seq(0.05, 0.29, by = 0.05)
B <- seq(0.1, 0.9, by = 0.05)
C <-  seq(0.05,0.45, by = 0.05)

time_effect <- 0



# Create a grid of all parameter combinations
param_grid <- expand.grid(overall_prob_outcome = overall_prob_outcome,
                          sd_outcome = sd_outcome,
                          B = B,
                          C = C,
                          time_effect_short = time_effect,
                          feature_std = feature_std)

# Set up parallel backend
n_cores <- parallel::detectCores() - 1  
cl <- makeCluster(n_cores)
registerDoParallel(cl)
old <- Sys.time() # get start time

# Parallel loop with foreach
result <- foreach(i = 1:nrow(param_grid), .combine = rbind, .packages = c("dplyr", "lme4", "pROC", "randomForest", "CorBin", "nlme", "performance"), .errorhandling = "pass") %dopar% {
  
  # Extract parameter combination for this iteration
  params <- param_grid[i, ]
  overall_prob <- params$overall_prob_outcome
  sd <- params$sd_outcome
  B_value <- params$B
  C_value <- params$C
  time_effect <- c(rep(params$time_effect_short,n_samples-1))
  #time_effect <- seq(0.5, 0.01, length.out = (n_samples-1))
  #time_effect = c(rep(0.9,45),rep(0,44))
  feature_std <- params$feature_std
  time_type = "DCP"
  
  # Generate data for this iteration
  features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B_value, C_value, overall_prob, sd, time_effect,time_type)
  features_sample_Astd <- features_sample[[2]]
  features_sample <- features_sample[[1]]
  
  ### Descriptives ####
  test <- features_sample %>% group_by(subject) %>% summarise(prob = mean(y))
  
  
  # Signal to noise ratio
  
  model <-glmer(
    y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + (1 | subject),                
    data = features_sample,
    family = "binomial"
  )
  
  s_n <- performance::r2(model)
  s_n$R2_marginal
  
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
    overall_prob_outcome = overall_prob,
    sd_outcome = sd,
    mean = mean(features_sample$y),
    sd = sd(test$prob),
    R2_conditional =  s_n$R2_conditional,
    R2_marginal = s_n$R2_marginal
  )
  
  # Save the result to a file in real time
  write.table(result_row, file = "signal_to_noise.csv", append = TRUE, sep = ",", col.names = !file.exists("signal_to_noise.csv"), row.names = FALSE)
  
  # Return the result 
  result_row
}


# Stop the cluster
stopCluster(cl)

# print elapsed time
Sys.time() - old # calculate difference


s_n <- read.csv("signal_to_noise.csv")

s_n$R2_marginal

mean(s_n$R2_marginal)
median(s_n$R2_marginal)
sd(s_n$R2_marginal)
min(s_n$R2_marginal)
max(s_n$R2_marginal)

####################################################
################ Run Simulation ####################
####################################################
source("Simulation_Functions_notshiny.R")

# Reported in Manuscript
n_features <- 10
n_samples <- 90
n_subjects <- 150
A <- 0 

# # Initialize parameters
overall_prob_outcome <- seq(0.1, 0.9, by = 0.1)
sd_outcome <- seq(0.05, 0.29, by = 0.05)
B <- seq(0.1, 0.9, by = 0.05)
C <-  seq(0.05,0.45, by = 0.05)

# time_effect <- 0
# timevariability <- 0
# feature_std <- 0.1 # population level feature generating process



# Figure A.1.1
time_effect <- 0.25
time_type <- "1-dependent"

# Create a grid of all parameter combinations
param_grid <- expand.grid(overall_prob_outcome = overall_prob_outcome,
                          sd_outcome = sd_outcome,
                          B = B,
                          C = C,
                          feature_std = feature_std)

#Figure A.1.2
time_type = "DCP"
time_effect =  seq(0.9, 0.01, length.out = (n_samples-1)) #0.2


# Figure A.1.3
time_type = 0
time_effect =  0

feature_std = c(0.5,1,2)

# Create a grid of all parameter combinations
param_grid <- expand.grid(overall_prob_outcome = overall_prob_outcome,
                          sd_outcome = sd_outcome,
                          B = B,
                          C = C,
                          feature_std = feature_std)


# Set up parallel backend
n_cores <- parallel::detectCores() - 1  
cl <- makeCluster(n_cores)
registerDoParallel(cl)
old <- Sys.time() # get start time

# Parallel loop with foreach
result <- foreach(i = 1:nrow(param_grid), .combine = rbind, .packages = c("dplyr", "lme4", "pROC", "randomForest", "CorBin", "nlme"), .errorhandling = "pass") %dopar% {
 
   # Extract parameter combination for this iteration
  params <- param_grid[i, ]
  overall_prob <- params$overall_prob_outcome
  sd <- params$sd_outcome
  B_value <- params$B
  C_value <- params$C
  time_effect <- time_effect

  feature_std <- params$feature_std
  time_type = time_type
  
  # Generate data for this iteration
  features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B_value, C_value, overall_prob, sd, time_effect,time_type)
  features_sample_Astd <- features_sample[[2]]
  features_sample <- features_sample[[1]]
  
  ### Descriptives ####
  test <- features_sample %>% group_by(subject) %>% summarise(prob = mean(y))
  
  # phi
  model <-lme(
    y ~ time,                  # fixed effects
    random = ~1 | subject, # random intercept per participant
    correlation = corAR1(form = ~ time | subject),  # AR(1) within participant
    data = features_sample
  )

  phi <- coef(model$modelStruct$corStruct, unconstrained = FALSE)

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
    time_effect = phi, #phi, # recovered
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
  #write.table(result_row, file = "simulation_results_02.csv", append = TRUE, sep = ",", col.names = !file.exists("simulation_results_02.csv"), row.names = FALSE)
  write.table(result_row, file = "simulation_results_time_A13.csv", append = TRUE, sep = ",", col.names = !file.exists("simulation_results_time_A13.csv"), row.names = FALSE)
  
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

time_effect <- 0
time_type <- 0




# Figure A.1.1
time_effect <- 0.25
time_type <- "1-dependent"



# #Figure A.2.2
time_type = "DCP"
time_effect =  seq(0.9, 0.01, length.out = (n_samples-1)) #0.2




# Set up parallel backend
n_cores <- parallel::detectCores() - 1  # Use all but two core
cl <- makeCluster(n_cores)
registerDoParallel(cl)
old <- Sys.time() # get start time

# Parallel loop with foreach
result <- foreach(i = 1:nrow(param_grid), .combine = rbind, .packages = c("dplyr", "lme4", "pROC", "randomForest", "CorBin", "nlme"), .errorhandling = "remove") %dopar% {
  # Extract parameter combination for this iteration
  params <- param_grid[i, ]
  overall_prob <- params$overall_prob_outcome
  sd <- params$sd_outcome
  n_samples <- params$n_samples
  n_subjects <- params$n_subjects
  time_effect <- time_effect
  time_type <- time_type
  
  # Generate data for this iteration
  features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B_value, C_value, overall_prob, sd, time_effect,time_type)
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
  write.table(result_row, file = "simulation_results_baseline_A22.csv", append = TRUE, sep = ",", col.names = !file.exists("simulation_results_baseline_A22.csv"), row.names = FALSE)
  
  # Return the result 
  result_row
}


# Stop the cluster
stopCluster(cl)
