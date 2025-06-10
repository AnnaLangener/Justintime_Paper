
library(dplyr)
library(lme4)
library(tibble)
library(dplyr)
library(pROC)
library(tidyr)
library(ggplot2)
library(purrr)
library(randomForest)

########## Binary Outcome ###########
#####################################

run_simulation_example_bin <- function(data = data, days_inlcuded = 200,participants_included = 10, icc_treshhold = 0, direction = ">"){ # direction for AUC
    
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






########### Other code to look at the results that was not directly used in the paper ###########
########### In the paper we presented separate plots, see Paper_Results.Rmd #####################


###### Quick overview of results #######

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

############## ROC CURVES PLOTS ######

library("purrr")

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

