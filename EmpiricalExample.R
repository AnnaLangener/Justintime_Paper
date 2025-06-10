
library(dplyr)
library(lme4)
library(tibble)
library(dplyr)
library(pROC)
library(tidyr)
library(ggplot2)
library(purrr)
library(randomForest)

run_simulation_example_icc <- function(data = data, days_inlcuded = 50,participants_included = 30, icc_treshhold = 0, direction = ">"){
  
  source("Simulation_UploadData.R")
  set.seed(1236)
  
  data = data  %>% group_by %>%
    group_by(user_id) %>%
    mutate(day = row_number())
  
  data <- data %>%
    group_by(user_id) %>%
    filter(n_distinct(day) >= days_inlcuded) %>%
    ungroup()
  data = data[data$day <= days_inlcuded,]
  print(length(unique(data$user_id)))

  sampled_ids <- sample(unique(data$user_id), size = participants_included, replace = FALSE)

  data <- data %>%
    filter(user_id %in% sampled_ids)

  
  results_list <- list()     # For summary stats & ICC
  model_list   <- list()     # For model results on shuffled data
  
  
  id_variable = "user_id"   
  time_variable = "day"
  n_features_upload = c("wantingAnything", "wantingIntensity", 
                        "withOthers", "feelingBadToGood", "energy", 
                        "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed")
  
  binary_vars = c("wantedToFeel_pleasant", "wantedToFeel_control", "wantedToFeel_goodAboutMe", 
                  "wantedToFeel_connected", "wantedToFeel_comfort", "wantedToFeel_calm", 
                  "wantedToFeel_excited", "wantedToFeel_competent", "wantedToFeel_acknowledged", 
                  "wantedToFeel_alert", "physically_hungry", 
                  "physically_tired", "physically_uncomfortable", "physically_energized", 
                  "physically_good", "feeling_angry", "feeling_anxious", 
                  "feeling_awe", "feeling_compassion", "feeling_content", "feeling_enthusiastic", 
                  "feeling_grateful", "feeling_guilty", "feeling_happy", "feeling_resentful", 
                  "feeling_restless", "feeling_sad")
  
  binary_vars = c("wantedToFeel_pleasant")
  
  non_binar <- names(data)[sapply(data, function(x) length(unique(na.omit(x))) > 2)]
  
  
  
  model_list <- data.frame(outcome = rep(NA,length(binary_vars)),icc = rep(NA,length(binary_vars)))
  
  for (i in seq_along(binary_vars)) {
    outcome_variable <- binary_vars[i]
    print(outcome_variable)
    
    ## --- Summary statistics ---
    overall_prob_outcome <- mean(data[[outcome_variable]], na.rm = TRUE)
    
    grouped <- data %>%
      group_by_at(id_variable) %>%
      summarise(prob = mean(.data[[outcome_variable]], na.rm = TRUE), .groups = "drop")
    
    sd_outcome <- sd(grouped$prob, na.rm = TRUE)
    
    model <- lmer(as.formula(paste0(outcome_variable, " ~ 1 + (1|", id_variable, ")")), data = data)
    var_random   <- as.data.frame(VarCorr(model))$vcov[1]
    var_residual <- attr(VarCorr(model), "sc")^2
    icc_data <- var_random / (var_random + var_residual)
    model_list$outcome[i] = outcome_variable
    model_list$icc[i] = icc_data
  }
  print(max(model_list$icc))
  return(model_list)
}


df_1 = run_simulation_example_icc(read.csv("esm_clean.csv"),200,10,0, ">") # most extreme example (data was already named esm_clean, no cleaning steps were taken)
df_2 = run_simulation_example_icc(read.csv("esm_clean.csv"),150,20,0, ">") # 
df_3 = run_simulation_example_icc(read.csv("esm_clean.csv"), 100,30,0, ">")
df_4 = run_simulation_example_icc(read.csv("esm_clean.csv"),60,100,0, ">") # ...
df_5 = run_simulation_example_icc(read.csv("esm_clean.csv"),50,241,0, ">") # all participants

df_all <- df_1 %>%
  left_join(df_2, by = "outcome") %>%
  left_join(df_3, by = "outcome") %>%
  left_join(df_4, by = "outcome") %>%
  left_join(df_5, by = "outcome")



#####################################
########## Full function ############
#####################################
library(furrr)
library(future)
plan(multisession, workers = parallel::detectCores()-1)  


data <-read.csv("esm_clean.csv")
run_simulation_example_bin <- function(data = data, days_inlcuded = 200,participants_included = 10, icc_treshhold = 0, direction = ">"){
    
  source("Simulation_UploadData.R")
  set.seed(1234)
  
  data = data  %>% group_by %>%
    group_by(user_id) %>%
    mutate(day = row_number())
  
  data <- data %>%
    group_by(user_id) %>%
    filter(n_distinct(day) >= days_inlcuded) %>%
    ungroup()
  data = data[data$day <= days_inlcuded,]
  
  sampled_ids <- sample(unique(data$user_id), size = participants_included, replace = FALSE)

  data <- data %>% 
    filter(user_id %in% sampled_ids)
  
  data$user_id <- sapply(data$user_id, match, unique(unlist(data$user_id)))
  
  
  results_list <- list()     # For summary stats & ICC
  model_list   <- list()     # For model results on shuffled data
  
  
  id_variable = "user_id"   
  time_variable = "day"

  
  n_features_upload = c("feelingBadToGood", "energy",
                        "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed")

  # n_features_upload = c( "physically_uncomfortable", "physically_tired", 
  #                        "physically_good")
  
  binary_vars = c("wantedToFeel_pleasant", "wantedToFeel_control", "wantedToFeel_goodAboutMe", 
    "wantedToFeel_connected", "wantedToFeel_comfort", "wantedToFeel_calm", 
    "wantedToFeel_excited", "wantedToFeel_competent", "wantedToFeel_acknowledged", 
    "wantedToFeel_alert", "physically_hungry", 
    "physically_tired", "physically_uncomfortable", "physically_energized", 
    "physically_good", "feeling_angry", "feeling_anxious", 
    "feeling_awe", "feeling_compassion", "feeling_content", "feeling_enthusiastic", 
    "feeling_grateful", "feeling_guilty", "feeling_happy", "feeling_resentful", 
    "feeling_restless", "feeling_sad")

  ## --- Features ---
  icc_data_var <- data.frame(variable = n_features_upload, icc = NA)
  for (i in seq_along(n_features_upload)) {
    model <- lmer(as.formula(paste0(n_features_upload[i], " ~ 1 + (1|", id_variable, ")")), data = data)
    var_random   <- as.data.frame(VarCorr(model))$vcov[1]
    var_residual <- attr(VarCorr(model), "sc")^2
    icc_data_var$icc[i] <- var_random / (var_random + var_residual)
  }

    n_features_upload = n_features_upload
    print(paste("ICC predictors mean:",mean(icc_data_var$icc)))
 
    ##### Set up parallel
    for (i in seq_along(binary_vars)) {
      
  # model_list <- future_map(binary_vars, function(outcome_variable) {
    
      
    suppressMessages({
    outcome_variable <- binary_vars[i]
    print(outcome_variable)
    

    data <- data %>%
      group_by(user_id) %>%
      filter(n_distinct(day) >= 3) %>%
      ungroup()
    
    ## --- Summary statistics ---
    overall_prob_outcome <- mean(data[[outcome_variable]], na.rm = TRUE)
    
    grouped <- data %>%
      group_by_at(id_variable) %>%
      summarise(prob = mean(.data[[outcome_variable]], na.rm = TRUE), .groups = "drop")
    
    sd_outcome <- sd(grouped$prob, na.rm = TRUE)
    
    model <- lmer(as.formula(paste0(outcome_variable, " ~ 1 + (1|", id_variable, ")")), data = data)
    var_random   <- as.data.frame(VarCorr(model))$vcov[1]
    var_residual <- attr(VarCorr(model), "sc")^2
    icc_data <- var_random / (var_random + var_residual)
    
    
    ## --- Shuffle data and run model ---
    data_shuffle <- shuffle_data(data, subject_var = id_variable, outcome_var = outcome_variable)
    colnames(data_shuffle)[colnames(data_shuffle) == id_variable]       <- "subject"
    colnames(data_shuffle)[colnames(data_shuffle) == time_variable]     <- "time"
   # colnames(data_shuffle)[colnames(data_shuffle) == outcome_variable]  <- "y"
    
    data_shuffle <- data_shuffle %>%
      group_by(subject) %>%
      mutate(y = lag(.data[[outcome_variable]], n = 1))
    
    data_shuffle <- data_shuffle[!is.na(data_shuffle$y),]
    data_shuffle$subject <- sapply(data_shuffle$subject, match, unique(unlist(data_shuffle$subject)))
    
    result_shuffle <- t(run_simulation_own(
      features_sample = data_shuffle,
      cv              = "record-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
      
    ))
    result_shuffle_df <- as.data.frame(result_shuffle)
    colnames(result_shuffle_df) <- paste0("shuffle_", colnames(result_shuffle_df))
    
    ## - Run Sliding Window 
    result_sliding_window_shuffle <- run_simulation_slidingwindow_own(
      features_sample = data_shuffle,
      n_bootstrap     = 1,
      windowsize      = 14,
      n_features      = n_features_upload,
      seed = "1234"
    )
    result_sliding_window_df_shuffle <- as.data.frame(t(result_sliding_window_shuffle))
    colnames(result_sliding_window_df_shuffle) <- paste0("movingwindowshuffle_", colnames(result_sliding_window_df_shuffle))
    
    ## --- Prepare real data for modeling ---
    data_model <- data
    colnames(data_model)[colnames(data_model) == id_variable]       <- "subject"
    colnames(data_model)[colnames(data_model) == time_variable]     <- "time"
    #colnames(data_model)[colnames(data_model) == outcome_variable]  <- "y"
    data_model <- data_model %>%
      group_by(subject) %>%
      mutate(y = lag(.data[[outcome_variable]], n = 1))
    
    data_model <- data_model[!is.na(data_model$y),]
    data_model$subject <- sapply(data_model$subject, match, unique(unlist(data_model$subject)))
    
    ## --- Run model on real data ---
    result_real <- t(run_simulation_own(
      features_sample = data_model,
      cv              = "record-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
      
      ))
    result_real_df <- as.data.frame(result_real)
    colnames(result_real_df) <- paste0("real_", colnames(result_real_df))
    
    ## --- Run centered model ---
    result_centered <- t(run_simulation_centering_own(
      features_sample = data_model,
      cv              = "record-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
    ))
    result_centered_df <- as.data.frame(result_centered)
    colnames(result_centered_df) <- paste0("centered_", colnames(result_centered_df))
    
    
    ## --- Split-Wise ---
    result_subject_wise <- t(run_simulation_own(
      features_sample = data_model,
      cv              = "subject-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed ="1234"
      
    ))
    result_subject_wise_df <- as.data.frame(result_subject_wise)
    colnames(result_subject_wise_df) <- paste0("subjectwise_", colnames(result_subject_wise_df))
    
    
    ## --- Run Sliding Window ---
    result_sliding_window <- run_simulation_slidingwindow_own(
      features_sample = data_model,
      n_bootstrap     = 1,
      windowsize      = 14,
      n_features      = n_features_upload,
      seed = "1234"
    )
    result_sliding_window_df <- as.data.frame(t(result_sliding_window))
    colnames(result_sliding_window_df) <- paste0("movingwindow_", colnames(result_sliding_window_df))
    
    
    ## --- Combine everything ---
    combined_df <- tibble(
      outcome = outcome_variable,
      overall_prob_outcome = overall_prob_outcome,
      sd_outcome = sd_outcome,
      num_subjects = n_distinct(data[[id_variable]]),
      num_samples = length(unique(data[[time_variable]])),
      icc_data = icc_data,
      icc_pred = mean(icc_data_var$icc)
    ) %>%
      bind_cols(result_shuffle_df, result_real_df, result_centered_df,result_subject_wise_df,result_sliding_window_df,result_sliding_window_df_shuffle)
    model_list[[i]] <- combined_df
        })
  }
  
  # Combine all results into one final dataframe

  final_results_df <- bind_rows(model_list)
  return(final_results_df)
}

###########################
library(ggplot2)
library(cowplot)
library(tidyr)
library(directlabels)

final_results_df_2_b = run_simulation_example_bin(read.csv("esm_clean.csv"),150,20,0, ">") # 
final_results_df_3_b = run_simulation_example_bin(read.csv("esm_clean.csv"), 100,30,0, ">")
final_results_df_4_b = run_simulation_example_bin(read.csv("esm_clean.csv"),60,100,0, ">") # ...
final_results_df_5_b = run_simulation_example_bin(read.csv("esm_clean.csv"),50,241,0, ">") # all participants

write.csv(final_results_df_2_b,"final_results_df_2_b.csv")
write.csv(final_results_df_3_b,"final_results_df_3_b.csv")
write.csv(final_results_df_4_b,"final_results_df_4_b.csv")
write.csv(final_results_df_5_b,"final_results_df_5_b.csv")


###### Between Plot #######

plot_data = function(data_viz,spanpar){
  # Reshape to long format
  df_long <- data_viz %>%
    pivot_longer(
      cols = c(`shuffle_AUC:`, `real_AUC random intercept only:`, `real_AUC:`,`movingwindowshuffle_AUC:`,`movingwindow_AUC:`,`subjectwise_AUC:`,`movingwindow_AUC random intercept only:`),
      names_to = "AUC_type",
      values_to = "AUC"
    )
  p = ggplot(df_long, aes(x = icc_data, y = AUC, color = AUC_type)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_smooth(se = FALSE, method = "loess",span = spanpar) +
    theme_minimal() +
    scale_color_manual(values = c(
      "shuffle_AUC:" = "#40916C",
      "real_AUC random intercept only:" ="#B5E48C" ,
      "real_AUC:" = "#74C69D",
      "movingwindowshuffle_AUC:" ="#4D8FAC",
      "movingwindow_AUC:" ="#9EC5E7",
      "movingwindow_AUC random intercept only:" = "#A2D2FF",
      "subjectwise_AUC:" ="#FFB5A7"
    ), 
    labels = c("shuffle_AUC:" = "No true Relationship (Record-Wise)", 
               "real_AUC random intercept only:"="Random Intercept only (Record-Wise)",
               "real_AUC:" ="Prediction (Record-Wise) ",
               "movingwindowshuffle_AUC:" ="No true Relationship (Moving-Window)",
               "movingwindow_AUC:" ="Prediction (Moving-Window)",
               "subjectwise_AUC:" ="Prediction (Subject-Wise)",
               "movingwindow_AUC random intercept only:" = "Random Intercept only (Moving-Window)")) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0.5, color = "grey", size = 0.8) +
    labs(color = "AUC Type") +
    ylim(0.4,0.95) +
    ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ","\n",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))
  return(p)
}


plot_data = function(data_viz,spanpar){
  # Reshape to long format
  df_long <- data_viz %>%
    pivot_longer(
      cols = c(`shuffle_AUC:`,`movingwindowshuffle_AUC:`),
      names_to = "AUC_type",
      values_to = "AUC"
    )
  p = ggplot(df_long, aes(x = icc_data, y = AUC, color = AUC_type)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_smooth(se = FALSE, method = "loess",span = spanpar) +
    theme_minimal() +
    scale_color_manual(values = c(
      "shuffle_AUC:" = "#40916C",
      "movingwindowshuffle_AUC:" ="#4D8FAC"
    ), 
    labels = c("shuffle_AUC:" = "No true Relationship (Record-Wise)",
               "movingwindowshuffle_AUC:" ="No true Relationship (Moving-Window)")) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0.5, color = "grey", size = 0.8) +
    labs(color = "AUC Type") +
    ylim(0.4,0.95) +
    ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ","\n",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))
  return(p)
}

spanpar = 0.6
#p1 = plot_data(final_results_df_1_b, spanpar)
p2 = plot_data(final_results_df_2_b, spanpar)
p3 = plot_data(final_results_df_3_b, spanpar)
p4 = plot_data(final_results_df_4_b, spanpar)
p5 = plot_data(final_results_df_5_b, spanpar)

legend_plot <- get_legend(
  p2 + theme(legend.position = "right")
)

# Convert legend to a plot object
legend_as_plot <- ggpubr::as_ggplot(legend_plot)
layout <- "
ABC
DE
"

final_plot <-  p2 + p3 + p4 + p5 + legend_as_plot +
  plot_layout(design = layout)

print(final_plot)


#### Wihtin Plot
plot_data = function(data_viz,spanpar){
  # Reshape to long format
  df_long <- data_viz %>%
    pivot_longer(
      cols = c(`shuffle_Mean AUC within-person:`, `real_Mean AUC within-person:`, `movingwindowshuffle_Mean AUC within-person:`,`subjectwise_Mean AUC within-person:`,`movingwindow_Mean AUC within-person:`,`movingwindow_Mean AUC random intercept only:`),
      names_to = "AUC_type",
      values_to = "AUC"
    )
  p = ggplot(df_long, aes(x = icc_data, y = AUC, color = AUC_type)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_smooth(se = FALSE, method = "loess",span = spanpar) +
    theme_minimal() +
    scale_color_manual(values = c(
      "shuffle_Mean AUC within-person:" = "#40916C",
      "real_Mean AUC within-person:" ="#74C69D" ,
      "movingwindowshuffle_Mean AUC within-person:" = "#4D8FAC",
      "movingwindow_Mean AUC within-person:" ="#9EC5E7",
      "subjectwise_Mean AUC within-person:" ="#A2D2FF",
      "movingwindow_Mean AUC random intercept only:" = "#FFB5A7"
    ),
    labels = c( "shuffle_Mean AUC within-person:" = "No true Relationship (Record Wise)", 
                "real_Mean AUC within-person:" = "Record-Wise (Random Forest)",
                "movingwindowshuffle_Mean AUC within-person:" ="No true Relationship (Moving Window)",
                "movingwindow_Mean AUC within-person:"="Moving-Window (Random Forest)",
                "subjectwise_Mean AUC within-person:"="Subject-Wise (Random Forest)",
                "movingwindow_Mean AUC random intercept only:" = "Random Intercept Only (Moving Window)"
              )) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 0.5, color = "grey", size = 0.8) +
    labs(color = "AUC Type") +
    ylim(0.4,0.95) +
    ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ","\n",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))
  return(p)
}



spanpar = 0.3
p1 = plot_data(final_results_df_1_b, spanpar)
p2 = plot_data(final_results_df_2_b, spanpar)
p3 = plot_data(final_results_df_3_b, spanpar)
p4 = plot_data(final_results_df_4_b, spanpar)
p5 = plot_data(final_results_df_5_b, spanpar)

legend_plot <- get_legend(
  p1 + theme(legend.position = "right")
)

# Convert legend to a plot object
legend_as_plot <- ggpubr::as_ggplot(legend_plot)
layout <- "
ABC
DEL
"

final_plot <- p1 + p2 + p3 + p4 + p5 + legend_as_plot +
  plot_layout(design = layout)

print(final_plot)

############## Individual vs Non individual Performance

data = read.csv("esm_clean.csv")

source("Simulation_UploadData.R")
set.seed(1234)

data = data  %>% group_by %>%
  group_by(user_id) %>%
  mutate(day = row_number())

data <- data %>%
  group_by(user_id) %>%
  filter(n_distinct(day) >= 150) %>%
  ungroup()
data = data[data$day <= 150,]

sampled_ids <- sample(unique(data$user_id), size = 20, replace = FALSE)

data <- data %>% 
  filter(user_id %in% sampled_ids)

data$user_id <- sapply(data$user_id, match, unique(unlist(data$user_id)))


id_variable = "user_id"   
time_variable = "day"



colnames(data)[colnames(data) == "user_id"] <- "subject"
colnames(data)[colnames(data) == "day"]     <- "time"
#colnames(data)[colnames(data) == "wantedToFeel_goodAboutMe"]  <- "y"

data <- data %>%
  group_by(subject) %>%
  mutate(y = lag(wantedToFeel_goodAboutMe, n = 1))

data <- data[!is.na(data$y),]
data$subject <- sapply(data$subject, match, unique(unlist(data$subject)))
## - Run Sliding Window
result <- run_simulation_slidingwindow_own(
  features_sample = data,
  n_bootstrap     = 1,
  windowsize      = 14,
  n_features      = c("feelingBadToGood", "energy",
                      "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed"),
  seed = 1234
)

pred_true = as.data.frame(result$ind)


roc_obj <- pROC::roc(as.numeric(as.character(pred_true$true)), 
                     as.numeric(as.character(pred_true$pred)), 
                     direction = "<")

# Create a data frame for plotting
roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

# Plot the ROC curve
a <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1D3557", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  coord_equal() +
  theme_minimal() +
  labs(title = paste("Overall ROC Curve - AUC:", round(roc_obj$auc, 2)),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme(text = element_text(family = 'Futura', size = 11))

# Compute ROC curve for each subject
roc_list <- pred_true %>%
  group_by(subject) %>%
  filter(length(unique(true)) > 1) %>%
  filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
  group_split() %>%
  map(function(df) {
    true_vals <- as.numeric(as.character(df$true))
    pred_vals <- as.numeric(as.character(df$pred))
    
    r <- pROC::roc(response = true_vals, predictor = pred_vals, direction = "<")
    
    data.frame(
      subject = unique(df$subject),
      specificity = r$specificities,
      sensitivity = r$sensitivities,
      auc = as.numeric(r$auc),
      n_pos = sum(true_vals == 1))
    
  }) %>%
  bind_rows()


# Step 1: Create labels
roc_list <- roc_list %>%
  mutate(subject_label = paste0("\nAUC = ", round(auc, 2), ",\nPositives = ", n_pos))

# Step 2: Set subject_label as factor and order by AUC
roc_list$subject_label <- factor(roc_list$subject_label, 
                                 levels = roc_list %>%
                                   distinct(subject_label, auc) %>%
                                   arrange(desc(auc)) %>%
                                   pull(subject_label))

# Step 3: Plot
b <- ggplot(roc_list, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1D3557") +
  geom_abline(linetype = "dashed", color = "gray") +
  facet_wrap(~ subject_label) +
  theme_minimal() +
  # coord_equal() +  # Optional, uncomment if you want square axes
  labs(
    title = paste("Individual ROC Curves, mean AUC = ", round(mean(roc_list$auc),2)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme(
    text = element_text(family = 'Futura', size = 11),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank()
  )



a + b


ind_auc <- pred_true %>%
  group_by(subject) %>%
  filter(length(unique(true)) > 1) %>%
  filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
  summarise(
    # threshold = {
    #   roc_obj <- pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)), direction = "<")
    #   coords_val <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = "threshold")
    #   as.numeric(coords_val[[1]])  # Safely extract numeric threshold
    # },
    auc_val = as.numeric(pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)), direction = "<"))),
    .groups = "drop"
  )

pred_true <- pred_true %>%
  left_join(ind_auc, by = "subject")

roc_obj <- pROC::roc(as.numeric(as.character(pred_true$true)), as.numeric(as.character(pred_true$pred)), direction = "<")
coords_val <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = "threshold")
pred_true$threshold = 0.5#coords_val$threshold

pred_true <- pred_true %>%
  mutate(pred = ifelse(pred > threshold, 1, 0))

# Add prediction correctness column
pred_true <- pred_true %>%
  mutate(prediction_correct = ifelse(true == pred, "Correct", "Incorrect"))

pred_true <- pred_true %>%
  arrange(subject,prediction_correct, true) %>%
  group_by(subject) %>%
  mutate(day = row_number())


library(pROC)
library(dplyr)


# Merge AUCs back into pred_true
pred_true$auc_val = round(pred_true$auc_val,2)
pred_true$auc_val[is.na(pred_true$auc_val)] = "-"


plot = ggplot(pred_true, aes(day, as.factor(subject), color = factor(true), shape = prediction_correct)) + 
  geom_point(size = 2, alpha = 0.8) +
  geom_text(aes(label = ifelse(!duplicated(subject), paste("AUC: ", auc_val), "")),
            x = -Inf, hjust = 0, size = 3, color = "gray40") +  # AUC on the left side
  scale_color_manual(values = c("0" = "#FFB5A7", "1" = "#40916C")) +
  scale_shape_manual(values = c("Correct" = 15, "Incorrect" = 0)) +
  labs(color = NULL, shape = NULL) +  # Removes legend titles
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.text.x = element_blank()     # Removes x-axis tick labels
  ) +
  xlim(-6,max(pred_true$day)) +
  labs(y = paste0("Within-Person"), x = "") +
  theme(text = element_text('Futura'),axis.title=element_text(size=14)) +
  ggtitle(paste("AUC evaluated over the whole dataset:", round(result$auc_value,2), "(Between-Person)"))
  
plot
######### ROC CURVES ######

library(pROC)
library(ggplot2)


result = run_simulation_own(
  features_sample = data,
  cv              = "record-wise",
  n_bootstrap     = 1,
  testsize        = 0.3,
  n_features      =  c("feelingBadToGood", "energy", 
                       "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed"),
  seed = 1234
)


result <- run_simulation_slidingwindow_own(
  features_sample = data,
  n_bootstrap     = 1,
  windowsize      = 14,
  n_features      = c("feelingBadToGood", "energy",
                      "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed"),
  seed = 1234
)

pred_true = as.data.frame(result$ind)
# Compute ROC curve across all data
roc_obj <- pROC::roc(as.numeric(as.character(pred_true$true)), 
               as.numeric(as.character(pred_true$pred)), 
               direction = "<")

# Create a data frame for plotting
roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

# Plot the ROC curve
a <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1D3557", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  coord_equal() +
  theme_minimal() +
  labs(title = paste("Overall ROC Curve - AUC:", round(roc_obj$auc, 2)),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme(text = element_text(family = 'Futura', size = 11))

# Compute ROC curve for each subject
roc_list <- pred_true %>%
  group_by(subject) %>%
  filter(length(unique(true)) > 1) %>%
  filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
  group_split() %>%
  map(function(df) {
    true_vals <- as.numeric(as.character(df$true))
    pred_vals <- as.numeric(as.character(df$pred))
    
    r <- pROC::roc(response = true_vals, predictor = pred_vals, direction = "<")
    
    data.frame(
      subject = unique(df$subject),
      specificity = r$specificities,
      sensitivity = r$sensitivities,
      auc = as.numeric(r$auc),
      n_pos = sum(true_vals == 1))
    
  }) %>%
  bind_rows()


# Step 1: Create labels
roc_list <- roc_list %>%
  mutate(subject_label = paste0("\nAUC = ", round(auc, 2), ",\nPositives = ", n_pos))

# Step 2: Set subject_label as factor and order by AUC
roc_list$subject_label <- factor(roc_list$subject_label, 
                                 levels = roc_list %>%
                                   distinct(subject_label, auc) %>%
                                   arrange(desc(auc)) %>%
                                   pull(subject_label))

# Step 3: Plot
b <- ggplot(roc_list, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1D3557") +
  geom_abline(linetype = "dashed", color = "gray") +
  facet_wrap(~ subject_label) +
  theme_minimal() +
  # coord_equal() +  # Optional, uncomment if you want square axes
  labs(
    title = paste("Individual ROC Curves, mean AUC = ", round(mean(roc_list$auc),2)),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme(
    text = element_text(family = 'Futura', size = 11),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank()
  )



a + b



### Get Table for example 1 #####

sim_1 <-  run_simulation_own(
  features_sample = shuffle_data(data, subject_var = "subject", outcome_var = "y"),
  cv              = "record-wise",
  n_bootstrap     = 1,
  testsize        = 0.3,
  n_features      =  c("feelingBadToGood", "energy", 
                       "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed")
)

sim_2 <- run_simulation_centering(features_sample,"row-wise",1,testsize = 0.3)
sim_3 <- run_simulation(features_sample,"subject-wise",1, testsize = 0.3)
sim_4 <- run_simulation_slidingwindow_own(
  features_sample = data,
  n_bootstrap     = 1,
  windowsize      = 14,
  n_features      = c("feelingBadToGood", "energy",
                     "thinkingOverAndOver", "lonely", "selfWorth", "appreciating", "stressed"))

table <- rbind(sim_1$overall_summary,sim_3$overall_summary, sim_4$overall_summary)
table$'AUC Overall' <- c(sim_1$auc_value,sim_3$auc_value,sim_4$auc_value)
colnames(table) <-  c("Mean AUC within-person:", "SD AUC within-person:", "% of AUC > 0.5 within-person:", "N included within-person:", "AUC Overall")
rownames(table) <- c("row-wise","subject-wise","rolling-window")
table <- round(table[,c(5,1,2,3,4)], 2)
table


########################
data = read.csv("esm_clean.csv") 
data = data  %>% group_by %>%
  group_by(user_id) %>%
  mutate(day = row_number())
ggplot(data, aes(x = day, y = wantedToFeel_pleasant)) +
  geom_point(color = "#FC9D9A", size = 0.5) +                       
  facet_wrap(~user_id) +               
  labs(
    x = "Time", 
    y = "Outcome", 
    title = paste("Within-person variability over time: \n88 participants (58.6%) do not experience within-person variability. ")
  ) +
  scale_y_continuous(breaks = c(0, 1)) +  # Ensure only 0 and 1 appear on the y-axis
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  guides(fill = "none") 

#########################

### Between person
data_viz = final_results_df
# Reshape to long format
df_long <- data_viz %>%
  pivot_longer(
    cols = c(`shuffle_AUC:`, `real_AUC random intercept only:`, `real_AUC:`, `centered_AUC:`),
    names_to = "AUC_type",
    values_to = "AUC_value"
  )

ggplot(df_long, aes(x = icc_data, y = AUC_value, color = AUC_type)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_smooth(se = FALSE, method = "loess") +
  theme_minimal() +
  scale_color_manual(values = c(
    "shuffle_AUC:" = "#B5E48C",
    "real_AUC random intercept only:" = "#40916C",
    "real_AUC:" = "#1B4332",
    "centered_AUC:" = "#FFB5A7"
  )) +
  labs(color = "AUC Type") +
  ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))

# within
# Reshape to long format
df_long <- data_viz %>%
  pivot_longer(
    cols = c(`shuffle_Mean AUC within-person:`, `real_Mean AUC within-person:`, `centered_Mean AUC within-person:`),
    names_to = "AUC_type",
    values_to = "AUC_value"
  )

ggplot(df_long, aes(x = icc_data, y = AUC_value, color = AUC_type)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +  # optional smoothing line
  theme_minimal() +
  scale_color_manual(values = c(
    "shuffle_Mean AUC within-person:" = "#B5E48C",
    "real_Mean AUC within-person:" = "#40916C",
    "centered_Mean AUC within-person:" = "#1B4332"
  )) +
  labs(color = "AUC Type") +
  ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))



## Centered vs Uncentered

final_results_df$diff_between = final_results_df$`centered_AUC:` - final_results_df$`real_AUC:`
final_results_df$diff_within = final_results_df$`centered_Mean AUC within-person:` - final_results_df$`real_Mean AUC within-person:` 



ggplot(final_results_df,aes(x = icc_data, y = diff_between )) +
  geom_point()



p3 <- ggplot(final_results_df, aes(x = icc_data, y = diff_between)) +
  geom_point(aes(color = ifelse(diff_between > 0, "Above Zero", "Below or Equal to Zero")), 
             alpha = 1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Above Zero" = "#FC9D9A", "Below or Equal to Zero" = "#83AF9B")) +
  theme_minimal() +
  ylab("Performance Difference \n(centered - uncentered)") +
  xlab("% of Total Variation Explained by Between Person Differences in the Outcome (ICC Outcome)") +
  ggtitle("Between-Person Performance") +
  geom_hline(yintercept = 0) 

p4 <- ggplot(final_results_df, aes(x = icc_data, y = diff_within)) +
  geom_point(aes(color = ifelse(diff_within > 0, "Above Zero", "Below or Equal to Zero")), 
             alpha =1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Above Zero" = "#FC9D9A", "Below or Equal to Zero" = "#83AF9B")) +
  theme_minimal() +
  ylab("Performance Difference \n(centered - uncentered)") +
  xlab("% of Total Variation Explained by Between Person Differences in the Outcome (ICC Outcome)") +
  ggtitle("Within-Person Performance") +
  geom_hline(yintercept = 0) 



combined_plot <- (p3 / p4) +
  plot_annotation(
    title = "Performance difference: centered vs. not-centered",
  ) 
combined_plot


###########################
library(dplyr)



ema <- read.csv("/Users/f007qrc/projects/Justintime_Paper/general_ema.csv")


# ema <- ema[!is.na(ema$stress),]
# ema <- ema[ema$stress != 3,]
# ema$stress[ema$stress == 1 | ema$stress == 2 ] <- 0
# ema$stress[ema$stress == 4 | ema$stress == 5] <- 1


ema <- ema[!is.na(ema$phq4_score),]
ema$phq4_score[ema$phq4_score < 8] <- 0
ema$phq4_score[ema$phq4_score>= 8 ] <- 1
ema <- ema[!is.na(ema$phq4_score),]

length(unique(ema$uid))

write.csv(ema,"/Users/f007qrc/projects/Justintime_Paper/ema.csv")


sensing <- read.csv("/Users/f007qrc/projects/Justintime_Paper/sensing.csv")

test = left_join(ema,sensing)


test = test[colnames(test) %in% c("uid", "day", "pam", "phq4.1", "phq4.2", "phq4.3", "phq4.4", "phq4_score", "social_level", "sse3.1", "sse3.2", "sse3.3", "sse3.4", "stress", "audio_amp_mean_ep_0", "audio_convo_duration_ep_0", "audio_voice_ep_0", "call_in_duration_ep_0", "loc_dist_ep_0", "loc_max_dis_from_campus_ep_0", "loc_home_still", 
                                  "loc_home_dur", "loc_leisure_dur", "loc_visit_num_ep_0", "sleep_duration","loc_leisure_dur","unlock_duration_ep_0")]

length(unique(test$uid))

test = test  %>% group_by %>%
  group_by(uid) %>%
  mutate(day = row_number())

test <- test %>%
  group_by(uid) %>%
  filter(n_distinct(day) >= 35) %>%
  ungroup()

length(unique(test$uid))

# 
# test <- test %>%
#   group_by(uid) %>%
#   filter(n_distinct(day) < 150) %>%
#   ungroup()


test <- test[test$day < 150,]

length(unique(test$uid))

write.csv(test,"/Users/f007qrc/projects/Justintime_Paper/all_studentlife_phq6.csv")
# 
# test = na.omit(test)
# length(test$uid)


###############
library(dplyr)
test = read.csv("test_new.csv")
test <- test %>%
  group_by(uid) %>%
  mutate(count = row_number() - 1)

test$stress[test$stress == 2] <- 1

write.csv(test,"test_new2.csv")
