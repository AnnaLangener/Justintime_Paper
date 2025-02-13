######### Load Packages ############
####################################

library(randomForest)
library(R.matlab)
library(foreach)
library(doParallel)
library(dplyr)
library(lme4)
library(tseries)
library(caret)
library(pROC)
library(boot)
library(plotly)


#set.seed(1234)  # Setting the seed for random number generation
##################################################
########### Set Simualtion Parameters ############
##################################################
# Basic Parameters from Study
n_features <- 10 # FIXED IN SAEB
n_samples <- 90  # Timepoints/ RECORDS (fixed in saeb) Number of samples per subject
n_subjects <- 150 # SUBJECTS (VARIED 4 TO 32) # Needs to be an even number

# Generating Outcome
overall_prob_outcome <- 0.3  # Overall probability of 1 (for the entire dataset) (e.g., tracking depression 0.146)
sqrt(overall_prob_outcome*(1-overall_prob_outcome))
sd_outcome <- 0.05 # Controls variability BETWEEN different subject (e.g., tracking depression 0.21)
# within_variability_outcome has to be: sd_outcome < sqrt(overall_prob_outcome*(1-overall_prob_outcome))
# (variance < mean(1-mean)

time_effect = FALSE

# Generating Features
A <- 0.05 # Relationship between features and outcome
feature_std <- 0.1 # population level feature generating process
B <- 0.8  # Cross-subject variability ("random effect") (added per participants for all timepoints)
# Add varying B (B[i], C[i] should work)
C <- 0.1 # Within-subject variability (added within participant for each timepoint)

##################################################
################# Data Simulation ################
##################################################
source("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/Simulation_Functions_notshiny.R")

features_sample <- create_data(n_features,n_samples,n_subjects,A,feature_std,B,C,overall_prob_outcome,sd_outcome,time_effect)
features_sample_Astd <- features_sample[[2]]
#features_sample_centered <- features_sample[[3]]
features_sample <- features_sample[[1]]


sim <- run_simulation_slidingwindow(features_sample,1,windowsize = 14)

auc(roc(as.numeric(as.character(sim$ind[[1]]$true)), as.numeric(as.character(sim$ind[[1]]$pred)), quiet = TRUE))

#features_sample_rep <- create_data(n_features,n_samples,n_subjects = n_subjects_rep,n_test,n_train,A,feature_std,B,C,overall_prob_outcome,sd_outcome,time_effect)
#features_sample <- read.csv("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/feature_sample.csv")

sim <- run_simulation(features_sample,"row-wise",1, testsize = 0.3)
write.csv(features_sample,"/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/feature_sample_005_low.csv" )

sim <- run_simulation_centering(features_sample,"row-wise",1,testsize = 0.3)

sim$ind[[1]]$true

sim <- run_simulation(features_sample,"subject-wise",1, testsize = 0.3)


features_sample_centered %>% 
  ggplot(aes(x=V1, fill=as.factor(y))) +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#FC9D9A","#83AF9B")) +
  theme_minimal() +
  labs(fill="")
mean(features_sample_centered$V1[features_sample_centered$y == 1])
mean(features_sample_centered$V1[features_sample_centered$y == 0])

features_sample %>% 
  ggplot(aes(x=V1, fill=as.factor(y))) +
  geom_histogram(alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#FC9D9A","#83AF9B")) +
  theme_minimal() +
  labs(fill="")

mean(features_sample$V1[features_sample$y == 1])
mean(features_sample$V1[features_sample$y == 0])


####################################################
############ Test different parameters #############
####################################################
source("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/Simulation_Functions_notshiny.R")

n_features <- 10 # FIXED IN SAEB
n_samples <- 90  # Timepoints/ RECORDS (fixed in saeb) Number of samples per subject
n_subjects <- 150 # SUBJECTS (VARIED 4 TO 32) # Needs to be an even number
time_effect = FALSE

A <- 0 # Relationship between features and outcome try 0.15 again
feature_std <- 0.1 # population level feature generating process

# Initialize parameters
overall_prob_outcome <- seq(0.1, 0.9, by = 0.1)  
overall_prob_outcome <- 0.5
sd_outcome <- seq(0.05, 0.29, by = 0.05)            
B <- seq(0.1, 1, by = 0.05)       
C <-  seq(0.05,0.5, by = 0.1)   


# Create a grid of all parameter combinations
param_grid <- expand.grid(overall_prob_outcome = overall_prob_outcome,
                          sd_outcome = sd_outcome,
                          B = B,
                          C = C)

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
  B_value <- params$B
  C_value <- params$C
  
  # Generate data for this iteration
  features_sample <- create_data(n_features, n_samples, n_subjects, A, feature_std, B_value, C_value, overall_prob, sd, time_effect)
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
    n_features <- n_features,
    n_samples <- n_samples,
    n_subjects <- n_subjects,
    feature_std = feature_std,
    A = A,
    B = B_value,
    C = C_value,
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
    auc_c = mean(sim3$auc_value),
    auc_c_individual = sim3$overall_summary$mean,
    auc_c_individual_sd = sim3$overall_summary$sd,
    sd_intercept = as.data.frame(VarCorr(model_pred))$sdcor[1], # B recovered
    sd_residual = as.data.frame(VarCorr(model_pred))$sdcor[2] # C recovered
  )
  
  # Save the result to a file in real time
  write.table(result_row, file = "/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_0_other.csv", append = TRUE, sep = ",", col.names = !file.exists("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_0_other.csv"), row.names = FALSE)
  
  # Return the result 
  result_row
}


# Stop the cluster
stopCluster(cl)

# print elapsed time
Sys.time() - old # calculate difference


result <- read.csv("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_test_new.csv")
ggplot(result, aes(x=icc, y= auc, color = sd)) +
  geom_point(alpha=0.5, size = 3) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
 # ggtitle(paste("A = ", A)) +
  # ggtitle(paste("A = ", A, "AUC_c min =",round(min(result$auc_c),2), "AUC_c max = ", round(max(result$auc_c),2),"AUC min =",round(min(result$auc),2), "AUC max = ", round(max(result$auc),2))) +
  geom_hline(yintercept = 0) 

####### Overall Performance ##########
library(patchwork)

p1 <- ggplot(result_0, aes(x=icc, y=auc_c_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 3) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
  ggtitle(paste("A = ", result_0$A)) +
  # ggtitle(paste("A = ", A, "AUC_c min =",round(min(result$auc_c),2), "AUC_c max = ", round(max(result$auc_c),2),"AUC min =",round(min(result$auc),2), "AUC max = ", round(max(result$auc),2))) +
  geom_hline(yintercept = 0.5) +
  ylim(range= c(0.5,0.9))

p3 <- ggplot(result_0.1, aes(x=icc, y=auc_c_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 3) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
  ggtitle(paste("A = ", result_0.1$A)) +
  # ggtitle(paste("A = ", A, "AUC_c min =",round(min(result$auc_c),2), "AUC_c max = ", round(max(result$auc_c),2),"AUC min =",round(min(result$auc),2), "AUC max = ", round(max(result$auc),2))) +
  geom_hline(yintercept = 0.5) +
  ylim(range= c(0.5,0.9))

combined_plot <- (p1 / p3) +
  plot_annotation(
    title = "Overall Performance (not centered)",
  ) 

########## Within-Person Perfromance ##########
p1 <-  ggplot(result, aes(x=icc, y=auc_c_individual-auc_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 2) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC centered - AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
  ggtitle(paste(
    "Within person results:",
    "A = ", A,
    "\nAUC centered mean = ", round(mean(result$auc_c_individual), 2),",",
    "AUC not centered mean = ", round(mean(result$auc_individual), 2)
  )) +
  geom_hline(yintercept = 0) 


p2 <- ggplot(result, aes(x=icc, y=auc_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 2) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
  geom_hline(yintercept = 0.5) +
  ylim(range= c(0.5,0.9))

p3 <- ggplot(result, aes(x=icc, y=auc_c_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 2) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
  geom_hline(yintercept = 0.5) +
  ylim(range= c(0.5,0.9))
  

combined_plot <- (p1 / p2/ p3) +
  plot_annotation(
    title = "Within-Person Performance",
  ) 



############ Relationship AUC & AUC individual ###########
p1 <- ggplot(result, aes(x=auc_c, y=auc_c_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 3) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC Individual centered") +
  xlab("AUC centered") +
  guides(col = guide_colourbar(title = "ICC predictor")) +
  ggtitle("A = 0.05") +
 # ggtitle(paste("A = ", A, "AUC_c min =",round(min(result$auc_c),2), "AUC_c max = ", round(max(result$auc_c),2),"AUC min =",round(min(result$auc),2), "AUC max = ", round(max(result$auc),2))) +
  geom_hline(yintercept = 0.5)
  
p2 <- ggplot(result, aes(x=auc, y=auc_individual, color = icc_pred)) +
  geom_point(alpha=0.5, size = 3) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC Individual") +
  xlab("AUC") +
  guides(col = guide_colourbar(title = "ICC predictor")) +
  ggtitle("A = 0.05") +
# ggtitle(paste("A = ", A, "AUC_c min =",round(min(result$auc_c),2), "AUC_c max = ", round(max(result$auc_c),2),"AUC min =",round(min(result$auc),2), "AUC max = ", round(max(result$auc),2))) +
  geom_hline(yintercept = 0.5) +
  ylim(range= c(0.5,0.9))



combined_plot <- (p2 / p1) +
  plot_annotation(
    title = "Relationship Overall AUC and within person AUC",
  ) 



plot_ly(result, 
        x = ~icc_pred, 
        y = ~icc, 
        z = ~auc, 
        type = 'scatter3d', 
        mode = 'markers',
        marker = list(size = 3),
        text = ~paste('ICC Pred:', icc_pred, '<br>ICC:', icc, '<br>AUC:', auc, '<br>SD:', sd),  # Custom hover text
        hoverinfo = 'text') 


#####################################
### Baseline model
source("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/Simulation_Functions_notshiny.R")

n_features <- 10 # FIXED IN SAEB
n_samples <- c(30,60,90,200,500)  # Timepoints/ RECORDS (fixed in saeb) Number of samples per subject
n_subjects <- c(30,60,90,200,500)# SUBJECTS (VARIED 4 TO 32) # Needs to be an even number
time_effect = FALSE

A <- 0 # Relationship between features and outcome try 0.15 again
feature_std <- 0.1 # population level feature generating process

# Initialize parameters
overall_prob_outcome <- seq(0.1, 0.9, by = 0.05)  
sd_outcome <- seq(0.05, 0.29, by = 0.05) 
B <- 0.1
C <- 0.1
#B <- seq(0.1, 1, by = 0.05)       
#C <-  seq(0.05,0.5, by = 0.1)   


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
    n_features <- n_features,
    n_samples <- n_samples,
    n_subjects <- n_subjects,
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
  write.table(result_row, file = "/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_baseline_new.csv", append = TRUE, sep = ",", col.names = !file.exists("/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_baseline_new.csv"), row.names = FALSE)
  
  # Return the result 
  result_row
}


# Stop the cluster
stopCluster(cl)
result$n_samples....n_samples
# print elapsed time
Sys.time() - old # calculate difference
ggplot(result, aes(x=icc, y=auc_value_base, color = n_samples....n_samples)) +
  geom_point(alpha=0.5, size = 1.2) +
  #scale_size(range = c(.001, 6), name="") +
  theme_minimal() +
  scale_colour_gradientn(colours = c("#2A363B","#83AF9B","#C8C8A9","#F9CDAD","#FC9D9A","#FE4365")) +
  ylab("AUC") +
  xlab("ICC outcome") +
  guides(col = guide_colourbar()) +
  ggtitle(paste("A = ", result$A)) +
  # ggtitle(paste("A = ", A, "AUC_c min =",round(min(result$auc_c),2), "AUC_c max = ", round(max(result$auc_c),2),"AUC min =",round(min(result$auc),2), "AUC max = ", round(max(result$auc),2))) +
  geom_hline(yintercept = 0.5) +
  ylim(range= c(0.5,1))



######## Run Diagnostics ############
#####################################

### To be PROVIDED by the USER ###

#features_sample <- read data in
feature_names <- colnames(features_sample)[1:n_features] # name of columns that contain features
outcome_variable <- "y" # name of outcome variable
id_variable <- "subject" # name of id/ subject variable
time_variable <- "time" # name of time variable
# add feature_std and check again if this is really needed
# testsize can be added
# number of ML repitions can be added

### Basic Descriptives (extract prevalence and SD) ###
overall_prob_outcome = mean(features_sample$y)
grouped <- features_sample %>% group_by(subject) %>% summarise(prob = mean(y))
sd_outcome = sd(grouped$prob)

n_subjects  <- nrow(unique(features_sample[id_variable]))  
n_samples <- nrow(unique(features_sample[time_variable]))  

# - To do: add plot from below

### ICC Predictor (extract B, C) ###
feature_names <- colnames(features_sample)[1:n_features]
id_variable <- "subject"
time_variable <- "time"

predictors_des <- data.frame(variable = rep(NA,n_features), icc = rep(NA,n_features), B = rep(NA,n_features), C = rep(NA,n_features))

for(i in 1:n_features){
  predictors_des$variable[i] <- feature_names[i]
  
  model_pred <- lmer(as.formula(paste0(feature_names[i], "~ 1 + ( 1  |", id_variable,")")), data = features_sample) 
  var_random <- as.data.frame(VarCorr(model_pred))$vcov[1]  
  var_residual <- attr(VarCorr(model_pred), "sc")^2     
  
  predictors_des$icc[i] <- var_random / (var_random + var_residual)
  predictors_des$B[i] <- as.data.frame(VarCorr(model_pred))$sdcor[1]
  predictors_des$C[i] <- as.data.frame(VarCorr(model_pred))$sdcor[2]
}


### Run Simulation (to get "baseline comparisons) ###


A <- 0 # Relationship between features and outcome
feature_std <- 0.1 # CHECK AGAIN population level feature generating process
B <- predictors_des$B  # Cross-subject variability ("random effect") (added per participants for all timepoints)
# Add varying B (B[i], C[i] should work)
C <- predictors_des$C # Within-subject variability (added within participant for each timepoint)



features_sample <- create_data(n_features,n_samples,n_subjects,A,feature_std,B,C,overall_prob_outcome,sd_outcome,time_effect)
sim <- run_simulation(features_sample,"row-wise",1, testsize = 0.3)
sim <- run_simulation(features_sample,"subject-wise",1, testsize = 0.3)

#### ICC Outcome (maybe delete?)
model <- lmer(y ~ 1 + (1 | subject), data = features_sample)
var_random <- as.data.frame(VarCorr(model))$vcov[1]  # Random intercept variance
var_residual <- attr(VarCorr(model), "sc")^2         # Residual variance
icc <- var_random / (var_random + var_residual)


# person center predictor variables to get model to generalize
# To enter condition in which model is unable to do it (person centering) > information gets thrown away
# Based on within person process that generalzie 

# 1. Person centered
# 2. Custom loss (keras, R)


# run same plot with some relationship

# difference in AUC plot (AUC RF, AUC Baseline)

# Low ICC, check true signal to noise, check how to quantify, start with low ICC
# within-person: signal, between-person: noise

# add some individual loss function > and check how it performs
####################################################


#### Data Visualization Outome #####
stat <- features_sample %>%
  group_by(subject) %>%
  summarise(stationary_p = adf.test(y)$p.value)

# Merge the p-values with the original dataset for plotting
data <- features_sample %>%
  left_join(stat, by = "subject")

# Plot the data and add p-values as annotations (TODO remove NA pvalue)
p2 <- ggplot(data, aes(x = time, y = y)) +
  geom_point(color = "#83AF9B") +                       
  facet_wrap(~subject) +               
  labs(x = "Time", y = "Outcome") +
  theme_minimal() +
  guides(fill = "none") 

p2

#geom_text(data = stat,               # Add p-values
#                  aes(x = Inf, y = Inf, label = paste0("p = ", round(stationary_p, 3))),
#          inherit.aes = FALSE,       # Do not inherit x and y aesthetics
#          hjust = 1.1, vjust = 1.8) +


#########################################################
##########################################################
# Other code

# Code to extract individual AUC from sim study
test <- sim$ind[[1]] %>%
  group_by(subject) %>%
  filter(length(unique(true)) > 1) %>%
  filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
  summarise(
    auc_val = auc(roc(as.numeric(as.character(true)), as.numeric(as.character(pred)), quiet = TRUE))[1],
    .groups = "drop" # Ungroup after summarizing
  )




#model <- glmer(y ~ 0 + V1 + V2 + (V1 + V2| subject), data = features_sample,family =  binomial(link = "logit"))
#summary(model) 
#anova(model)
#coef(summary(model))

#model <- lm(y ~ 0 + v1, data = features_sample)
#coef(summary(model))



# Stationarity Binary
stat <- data %>%
  group_by(subject) %>%
  summarise(stationary_p = summary(lm(y ~ time))$coefficients[2, 4])

# test time effect (randomness needs to be deleted)
test = data %>% group_by(time) %>% summarise(n = sum(y))

p1 <- ggplot(test, aes(x = time, y = n)) +
  geom_line() +
  labs(x = "Time", y = "Count") +  # Add axis labels (optional, based on context)
  theme_minimal()

#hist(data$y)

################ Empirical Example ###############
exp <- read.csv('/Users/f007qrc/Documents/Unishare_Saved/Tracking Depression Study/ESM_data_rollingrmssd.csv')

binarize_rmssd <- function(df) {
  # Calculate mean and standard deviation for each participant within the normalized day range (-4 to 24)
  stats <- df %>%
    filter(normalized_day >= -4 & normalized_day <= 90) %>%
    group_by(uid) %>%
    summarise(
      std = sd(rolling_rmssd, na.rm = TRUE),
      mean = mean(rolling_rmssd, na.rm = TRUE),
    )
  
  df <- df %>%
    left_join(stats, by = "uid")
  
  df["treshhold"] = df["mean"] + df["std"]
  df["rmmssd_binary"] = ifelse(df["rolling_rmssd"] > df["mean"] + df["std"], 1,0)
  
  return(df)
}


#exp <- binarize_rmssd(exp)
rmssd = as.numeric(unlist(na.omit(exp["rolling_rmssd"])))

exp$rmssd_binary = ifelse(exp$rolling_rmssd > mean(rmssd)+  sd(rmssd), 1,0)

mean(exp$rmssd_binary, na.rm = TRUE)
test <- exp %>% group_by(uid) %>% summarise(prob = mean(rmssd_binary, na.rm = TRUE))
sd(test$prob, na.rm = TRUE)



mean(exp$rolling_rmssd, na.rm = TRUE)
test <- exp %>% group_by(uid) %>% summarise(prob = mean(rolling_rmssd, na.rm = TRUE))
sd(test$prob, na.rm = TRUE)

#######################
#######################
exp <- read.csv('/Users/f007qrc/projects/Justintime/example_justintime.csv')


origin = 14
horizon = 7
exp = exp[,-1]
colnames(exp)[colnames(exp) == "uid"]  <- "subject"
colnames(exp)[colnames(exp) == "normalized_day"]  <- "time"
colnames(exp)[colnames(exp) == "rolling_rmssd"]  <- "y"

exp <- exp %>%
  mutate(subject = gsub("[^0-9]", "", subject))

exp$subject <- as.numeric(exp$subject)

features_sample = exp

n_bootstrap = 1

windowsize = 14

n_samples = length(unique(features_sample$time))
n_features = 10
acc <- numeric(n_bootstrap)
auc_value <- numeric(n_bootstrap)
acc_rep <- numeric(n_bootstrap)
true_list <- list()  # To store true_sw for each iteration
pred_list <- list() 
subject <- list() 
ind <- list()
auc_value_meansw <- numeric(n_bootstrap)

acc_base <- numeric(n_bootstrap)
auc_value_base <- numeric(n_bootstrap)
acc_rep_base <- numeric(n_bootstrap)
true_list_base <- list()  # To store true_sw for each iteration
pred_list_base <- list() 
ind_base <- list()

print("sliding-window")

i = 1

#  Prepare training and testing sets
timeSlices <- SlidingWindow_CV(n_samples,windowsize,1) #window size, prediction horizon
timeSlices <- SlidingWindow_CV(n_samples,windowsize,7) #window size, prediction horizon
timeSlices <- createTimeSlices(seq(0,n_samples),14,1,FALSE)
trainSlices <- timeSlices[[1]]
testSlices <- timeSlices[[2]]

# # Example usage
# result <- sliding_window_cv(maxlookback = 10 + 6)
# trainSlices <- result$train_index
# testSlices <- result$val_index
# #testSlices <- result$test_index


lengthtest <- length(unlist(testSlices))

acc_sw <- numeric(lengthtest)
auc_value_sw <- numeric(lengthtest)
acc_sw_base <- numeric(lengthtest)
auc_value_base_sw <- numeric(lengthtest)
ind[[i]] = data.frame(subject = as.character(), true = as.numeric(), pred = as.numeric())
ind_t <- data.frame()
ind_t_base <- data.frame()


features_sample <- features_sample[!is.na(features_sample$y),]
    
for(k in 1:length(trainSlices)){
  train_X = features_sample[features_sample$time %in% trainSlices[[k]], 1:n_features]
  test_X = features_sample[features_sample$time %in% testSlices[[k]], 1:n_features]
  train_Y = as.factor(features_sample$y[features_sample$time %in% trainSlices[[k]]])
  test_Y = as.factor(features_sample$y[features_sample$time %in% testSlices[[k]]])
  
  # Train and evaluate the Random Forest model 
  #### Baseline
  Baseline <- glmer(
    "y ~ 1 + (1|subject)", # Random intercept for 'subject'
    data = features_sample[features_sample$time %in% trainSlices[[k]],], 
    family = "binomial" # Logistic regression
  )
  
  class_pred_base <- predict(Baseline,re.form = ~(1 | subject), newdata = features_sample[features_sample$time %in% testSlices[[k]],], type = "response")
 # class_pred_base <- ifelse(class_pred_base > 0.5, 1, 0)
 # acc_sw_base[k] <- mean(as.numeric(as.character(class_pred_base)) == test_Y)
  
  if (length(unique(test_Y)) == 2) {
    roc_curve_base <- roc(features_sample$y[features_sample$time %in% testSlices[[k]]],  as.numeric(as.character(class_pred_base)),quiet = TRUE)
    auc_value_base_sw[k] <- auc(roc_curve_base)
  } else {
    auc_value_base_sw[k] <- NA
  }
  
  ind_t_base <- rbind(
    ind_t_base,
    data.frame(
      subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
      true = features_sample$y[features_sample$time %in% testSlices[[k]]],                                   
      pred = class_pred_base
    )
  )
  
  
  print(paste("Window ",k,": ", " Baseline: ",round(auc_value_base_sw[k],2), sep = ""))
}




print(ind_t$true)
auc_value[i] <- auc(roc(as.numeric(as.character(ind_t$true)), as.numeric(as.character(ind_t$pred)), quiet = TRUE))
#auc_value_meansw[i] <- mean(auc_value_sw,na.rm = TRUE)
ind[[i]] <- ind_t
auc_value_base[i] <- auc(roc(ind_t_base$true, as.numeric(as.character(ind_t_base$pred)), quiet = TRUE))


results_summary <- lapply(ind, function(ind_result) {
  processed <- ind_result %>%
    group_by(subject) %>%
    filter(length(unique(true)) > 1) %>%
    filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
    summarise(
      auc_val = auc(roc(as.numeric(as.character(true)), as.numeric(as.character(pred)), quiet = TRUE))[1],
      .groups = "drop" # Ungroup after summarizing
    )
  
  # Bootstrap-level statistics
  mean_auc <- mean(processed$auc_val, na.rm = TRUE)
  sd_auc <- sd(processed$auc_val, na.rm = TRUE)
  percent_above_0_5 <- (sum(processed$auc_val > 0.5, na.rm = TRUE) / nrow(processed)) * 100
  n <- nrow(processed) # Number of observations in this bootstrap
  
  # Return a summary as a list
  list(
    mean_auc = mean_auc,
    sd_auc = sd_auc,
    percent_above_0_5 = percent_above_0_5,
    total_n = n,
    bootstrap_results = processed # Include detailed bootstrap-level data
  )
})
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean = x$mean_auc, sd = x$sd_auc, percent_above_0_5 = x$percent_above_0_5, total_n = x$total_n)
  }))
  
  results <- list(
    auc_value = auc_value,
    ind = ind,
    overall_summary = overall_summary
  )
  
  print(paste("Baseline Mean AUC:",mean(auc_value_base)))
  print("Model Results:")
  print(paste("Mean AUC:",mean(auc_value)))
  print("Individual Results Summary:")
  print(overall_summary)




#######################
######################



# Intraclass correlation
model <- lmer(rmssd_binary ~ 1 + (1 | uid), data = exp)
var_random <- as.data.frame(VarCorr(model))$vcov[1]  # Random intercept variance
var_residual <- attr(VarCorr(model), "sc")^2         # Residual variance
icc <- var_random / (var_random + var_residual)

icc



r1 <- read.csv('/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results05.csv')
r2 <- read.csv('/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_0.csv')
r3 <- read.csv('/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_08.csv')
r4 <- read.csv('/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results_02.csv')



r <- rbind(r1,r2,r3,r4)

write.csv(r,'/Users/f007qrc/Library/CloudStorage/GoogleDrive-anna.m.langener@dartmouth.edu/My Drive/Darmouth Drive/ML CrossValidation Project/simulation_results.csv')
####


# Define train, val, & test set split with ROLLING WINDOW FOR TRAIN/VAL/TEST
sliding_window_cv <- function(maxlookback) {
  train_index <- list()
  val_index <- list()
  test_index <- list()
  
  # Initial indices (adapt starting values)
  train_index[[1]] <- seq(-4 + maxlookback, -4 + maxlookback + 14 - 1)
  val_index[[1]] <- seq(-4 + maxlookback + 20, -4 + maxlookback + 27 - 1)
  test_index[[1]] <- c(-4 + maxlookback + 33)
  
  counter <- 1
  index <- test_index[[1]]
  value <- -4 + maxlookback + 14
  
  while (index < 83 ) {
    # Expand training window by one day
    value <- value + 1
    train_index[[counter + 1]] <- seq(-4 + maxlookback, value)
    
    # Shift validation and test set by one day
    val_index[[counter + 1]] <- val_index[[counter]] + 1
    test_index[[counter + 1]] <- test_index[[counter]] + 1
    
    # Calculate new "index"
    index <- test_index[[counter + 1]]
    counter <- counter + 1
  }
  
  return(list(train_index = train_index, val_index = val_index, test_index = test_index))
}

# Example usage
result <- sliding_window_cv(maxlookback = 10 + 6)
train_index <- result$train_index 
val_index <- result$val_index
test_index <- result$test_index

run_simulation_slidingwindow(features_sample,1,windowsize = 14)
