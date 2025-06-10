
download.OSF.file <- function(GUID,Access_Token=NULL,file_name)
{
  require(httr)
  require(rjson)
  #search for file private/public status
  GETurl <- paste0("https://api.osf.io/v2/files/",GUID)
  req <- GET(GETurl, write_disk("test",overwrite=T))
  json_data <- fromJSON(file = "test")
  if (length(json_data$data) > 0){
    req1 <- GET(json_data$data$links$download,
                write_disk(file_name, overwrite = TRUE))
    print(paste0("The file has been downloaded to your working directory as: ",
                 file_name))
  }
  else if (length(Access_Token) == 1){
    if (grepl("https://osf.io",Access_Token)==TRUE){
      req1 <- GET(paste0("https://api.osf.io/v2/files/",GUID,"/",gsub(".*/","",Access_Token)),
                  write_disk("test", overwrite = TRUE))
      json_data <- fromJSON(file = "test")
      if (length(json_data$data) > 0){
        req1 <- GET(json_data$data$links$download,
                    write_disk(file_name, overwrite = TRUE))
        print(paste0("The file has been downloaded to your working directory as: ",
                     file_name))
      }
      else{
        print(json_data$errors[[1]]$detail[1])
      }
    }
    else if (grepl("https://osf.io",Access_Token)==FALSE){
      req1 <- GET(paste0("https://api.osf.io/v2/files/",GUID),
                  write_disk("test", overwrite = TRUE),
                  add_headers("Authorization" = paste0("Bearer ",Access_Token)))
      json_data <- fromJSON(file = "test")
      if (length(json_data$data) > 0){
        req1 <- GET(json_data$data$links$download,
                    write_disk(file_name, overwrite = TRUE),
                    add_headers("Authorization" = paste0("Bearer ",Access_Token)))
        print(paste0("The file has been downloaded to your working directory as: ",
                     file_name))
      }
      else{
        print(json_data$errors[[1]]$detail[1])
      }
    }
    else{
      print(json_data$errors[[1]]$detail[1])
    }
  }
  else{
    print(json_data$errors[[1]]$detail[1])
  }
}
download.OSF.file(GUID="grm4t",file_name="Clinical EMA data.xlsx")



data = data  %>% group_by %>%
  group_by(user_id) %>%
  mutate(day = row_number())

data <- data %>%
  group_by(user_id) %>%
  filter(n_distinct(day) >= days_inlcuded) %>%
  ungroup()
data = data[data$day <= days_inlcuded,]

sampled_ids <- sample(unique(data$user_id), size = participants_included, replace = FALSE)




data[data == 999] <- NA
############### Run Empirical Example #########
###############################################

data_all = read.csv("esm_clean.csv")
run_simulation_example <- function(data_all = data, days_inlcuded = 50,participants_included = 30, icc_treshhold = 0, direction = ">"){

  data_all = data_all  %>% group_by %>%
    group_by(user_id) %>%
    mutate(day = row_number())
  
  data_all <- data_all %>%
    group_by(user_id) %>%
    filter(n_distinct(day) >= days_inlcuded) %>%
    ungroup()
  data_all = data_all[data_all$day <= days_inlcuded,]
  
  source("Simulation_continous.R")
  

  
  sampled_ids <- sample(unique(data_all$user_id), size = participants_included, replace = FALSE)

  data_all <- data_all %>%
    filter(user_id %in% sampled_ids)

  
  results_list <- list()     # For summary stats & ICC
  model_list   <- list()     # For model results on shuffled data
  
  
  id_variable = "user_id"   
  time_variable = "day"

  n_features_all =  c("feelingBadToGood", "energy",
                      "thinkingOverAndOver", "lonely", "selfWorth", "appreciating","stressed")

  binary_vars =  c("feelingBadToGood", "energy",
                   "thinkingOverAndOver", "lonely", "selfWorth", "appreciating","stressed")
  
  
  ## --- Features ---
 
  for (i in seq_along(binary_vars)) {   
    set.seed(12341688)
    suppressMessages({
    ## --- Features ---
    outcome_variable <- binary_vars[i]
    print(outcome_variable)
    
    n_features_upload = n_features_all[n_features_all != outcome_variable]
    icc_data_var <- data.frame(variable = n_features_upload, icc = NA)
    
    data = data_all[!is.na(data_all[outcome_variable]),]
    
    
    for (k in seq_along(n_features_upload)) {
      model <- lmer(as.formula(paste0(n_features_upload[k], " ~ 1 + (1|", id_variable, ")")), data = data)
      var_random   <- as.data.frame(VarCorr(model))$vcov[1]
      var_residual <- attr(VarCorr(model), "sc")^2
      icc_data_var$icc[k] <- var_random / (var_random + var_residual)
    }
    
    print(paste("ICC predictors mean:",mean(icc_data_var$icc)))
    
    
   
    # 
    # data <- data %>%
    #   group_by(ID) %>%
    #   filter(n_distinct(Day) >= 3) %>%
    #   ungroup()
    
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
    #colnames(data_shuffle)[colnames(data_shuffle) == outcome_variable]  <- "y"
    
    data_shuffle <- data_shuffle %>%
      group_by(subject) %>%
      mutate(y = lag(.data[[outcome_variable]], n = 1))
    
    data_shuffle <- data_shuffle[!is.na(data_shuffle$y),]
    data_shuffle$subject <- sapply(data_shuffle$subject, match, unique(unlist(data_shuffle$subject)))
    
    
    result_shuffle <- run_simulation_own(
      features_sample = data_shuffle,
      cv              = "record-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
    )
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
    result_sliding_window_df_shuffle <- as.data.frame(result_sliding_window_shuffle)
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
    result_real <- run_simulation_own(
      features_sample = data_model,
      cv              = "record-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
    )
    result_real_df <- as.data.frame(result_real)
    colnames(result_real_df) <- paste0("real_", colnames(result_real_df))
    
    ## --- Run centered model ---
    result_centered <- run_simulation_centering_own(
      features_sample = data_model,
      cv              = "record-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
    )
    result_centered_df <- as.data.frame(result_centered)
    colnames(result_centered_df) <- paste0("centered_", colnames(result_centered_df))
    
    
    ## --- Split-Wise ---
    result_subject_wise <- run_simulation_own(
      features_sample = data_model,
      cv              = "subject-wise",
      n_bootstrap     = 1,
      testsize        = 0.3,
      n_features      = n_features_upload,
      seed = "1234"
    )
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
    result_sliding_window_df <-  as.data.frame(result_sliding_window)
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
      bind_cols(result_shuffle_df, result_real_df, result_centered_df,result_subject_wise_df,result_sliding_window_df, result_sliding_window_df_shuffle)
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
library(patchwork)



#final_results_df_1 = run_simulation_example(read.csv("esm_clean.csv"),200,10,0, ">") # most extreme example
final_results_df_2 = run_simulation_example(read.csv("esm_clean.csv"),150,20,0, ">") # 
final_results_df_3 = run_simulation_example(read.csv("esm_clean.csv"), 100,30,0, ">")
final_results_df_4 = run_simulation_example(read.csv("esm_clean.csv"),60,100,0, ">") # ...
final_results_df_5 = run_simulation_example(read.csv("esm_clean.csv"),50,241,0, ">") # all participants


write.csv(final_results_df_2,"final_results_df_2.csv")
write.csv(final_results_df_3,"final_results_df_3.csv")
write.csv(final_results_df_4,"final_results_df_4.csv")
write.csv(final_results_df_5,"final_results_df_5.csv")

plot_data = function(data_viz,spanpar){
  # Reshape to long format
  df_long <- data_viz %>%
    pivot_longer(
    #  cols = c(`shuffle_mae`, `real_mae_base`, `real_mae`,`centered_mae`,`subjectwise_mae`),
      cols = c(`shuffle_cor`, `real_cor_base`, `real_cor`,`movingwindow_cor`,`subjectwise_cor`,movingwindowshuffle_cor,movingwindow_cor_base),
      names_to = "type",
      values_to = "Correlation"
    )
  p = ggplot(df_long, aes(x = icc_data, y = Correlation, color = type)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_smooth(se = FALSE, method = "loess",span = spanpar) +
    theme_minimal() +
    scale_color_manual(values = c(
      "shuffle_cor" = "#40916C",
      "real_cor_base" ="#B5E48C" ,
      "real_cor" = "#74C69D",
      "movingwindowshuffle_cor" = "#4D8FAC",
      "movingwindow_cor" = "#9EC5E7",
      "movingwindow_cor_base" = "#A2D2FF",
      "subjectwise_cor" ="#FFB5A7"
    ), labels = c("shuffle_cor" = "No true Relationship (Record-Wise)", 
                  "real_cor_base"="Random Intercept only (Record-Wise)",
                  "real_cor" ="Record-Wise (Random Forest)",
                  "movingwindowshuffle_cor" ="No true Relationship (Moving-Window)",
                  "movingwindow_cor" ="Moving-Window (Random Forest)",
                  "subjectwise_cor" = "Subject-Wise (Random Forest)",
                  "movingwindow_cor_base" = "Random Intercept only (Moving-Window)"
    )) +
    theme(legend.position = "none") +
    labs(color = "RSQ Value") +
    ylim(0,1) +
   # xlim(0.1,1) +
    ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ","\n",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))
  return(p)
}


plot_data = function(data_viz,spanpar){
  # Reshape to long format
  df_long <- data_viz %>%
    pivot_longer(
      #  cols = c(`shuffle_mae`, `real_mae_base`, `real_mae`,`centered_mae`,`subjectwise_mae`),
      cols = c(`shuffle_cor`,movingwindowshuffle_cor),
      names_to = "type",
      values_to = "Correlation"
    )
  p = ggplot(df_long, aes(x = icc_data, y = Correlation, color = type)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_smooth(se = FALSE, method = "loess",span = spanpar) +
    theme_minimal() +
    scale_color_manual(values = c(
      "shuffle_cor" = "#40916C",
      "movingwindowshuffle_cor" = "#4D8FAC"
    ), labels = c("shuffle_cor" = "No true Relationship (Record-Wise)", 
                  "movingwindowshuffle_cor" ="No true Relationship (Moving-Window)"
    )) +
    theme(legend.position = "none") +
    labs(color = "RSQ Value") +
    ylim(0,1) +
    # xlim(0.1,1) +
    ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ","\n",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))
  return(p)
}


spanpar = 0.3
#p1 = plot_data(final_results_df_1, spanpar)
p2 = plot_data(final_results_df_2, spanpar)
p3 = plot_data(final_results_df_3, spanpar)
p4 = plot_data(final_results_df_4, spanpar)
p5 = plot_data(final_results_df_5, spanpar)

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
final_plot


plot_data = function(data_viz,spanpar){
  # Reshape to long format
  df_long <- data_viz %>%
    pivot_longer(
      cols = c(`shuffle_overall_summary.mean_cor`,`real_overall_summary.mean_cor`,movingwindowshuffle_overall_summary.mean_cor,`centered_overall_summary.mean_cor`,`subjectwise_overall_summary.mean_cor`,`movingwindow_overall_summary.mean_cor`,movingwindow_overall_summary_base.mean_cor_base),
      names_to = "Type",
      values_to = "Correlation"
    )
  p = ggplot(df_long, aes(x = icc_data, y = Correlation, color = Type)) +
    geom_point(alpha = 0.8, size = 1) +
    geom_smooth(se = FALSE, method = "loess",span = spanpar) +
    theme_minimal() +
    scale_color_manual(values = c(
      "shuffle_overall_summary.mean_cor" = "#40916C",
      "real_overall_summary.mean_cor" = "#74C69D",
      "movingwindowshuffle_overall_summary.mean_cor" ="#4D8FAC",
      "subjectwise_overall_summary.mean_cor" ="#9EC5E7",
      "movingwindow_overall_summary.mean_cor:" ="#A2D2FF",
      "movingwindow_overall_summary_base.mean_cor_base" = "#FFB5A7"
    ), labels = c("shuffle_overall_summary.mean_cor" = "No true Relationship (Record-Wise)", 
                  "real_overall_summary.mean_cor" ="Record-Wise (Random Forest)",
                  "movingwindowshuffle_overall_summary.mean_cor" ="No true Relationship (Moving-Window)",
                  "movingwindow_overall_summary.mean_cor" ="Moving-Window (Random Forest)",
                  "subjectwise_overall_summary.mean_cor" = "Subject-Wise (Random Forest)",
                  "movingwindow_overall_summary_base.mean_cor_base" = "Random Intercept only (Moving-Window)"
    )) +
    theme(legend.position = "none") +
    labs(color = "AUC Type") +
    #ylim(-10,10) +
    ggtitle(paste0("N: ", max(data_viz$num_subjects),", ", "Timepoints: ", max(data_viz$num_samples),", ","\n",  "ICC Predictors (mean): ", round(data_viz$icc_pred,2)))
  return(p)
}
#p1 = plot_data(final_results_df_1, spanpar)
p2 = plot_data(final_results_df_2, spanpar)

p3 = plot_data(final_results_df_3, spanpar)

p4 = plot_data(final_results_df_4, spanpar)
p5 = plot_data(final_results_df_5, spanpar)

legend_plot <- get_legend(
  p1 + theme(legend.position = "right")
)

# Convert legend to a plot object
legend_as_plot <- ggpubr::as_ggplot(legend_plot)
layout <- "
ABC
DEL
"

final_plot <- p2 + p3 + p4 + p5 + legend_as_plot

final_plot
