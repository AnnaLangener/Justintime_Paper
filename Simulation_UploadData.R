##############################
### Simulation Functions #####
##############################

shuffle_data <- function(data, subject_var, outcome_var) {
  # Shuffle participant IDs:
  unique_subjects <- unique(data[[subject_var]])
  shuffled_subjects <- sample(unique_subjects)
  mapping <- data.frame(old_subject = unique_subjects, new_subject = shuffled_subjects,
                        stringsAsFactors = FALSE)

  # Merge mapping with data and update the subject column:
  data <- data %>%
    left_join(mapping, by = setNames("old_subject", subject_var)) %>%
    mutate(!!sym(subject_var) := new_subject) %>%
    select(-new_subject)
  
  # Shuffle the outcome variable within each subject:
  data <- data %>%
    group_by(!!sym(subject_var)) %>%
    mutate(!!sym(outcome_var) := sample(!!sym(outcome_var))) %>%
    ungroup()
  
  return(data)
}


########################################################
############## Cross Validation Strategies #############
########################################################

############# Subject Wise and Row Wise ################

run_simulation_own <- function(features_sample,cv,n_bootstrap,testsize, seed = "12361488",n_features){
  set.seed(seed)
  acc <- numeric(n_bootstrap)
  auc_value <- numeric(n_bootstrap)
  acc_rep <- numeric(n_bootstrap)
  true_list <- list()  # To store true_sw for each iteration
  pred_list <- list() 
  subject <- list() 
  ind <- list()
  
  acc_base <- numeric(n_bootstrap)
  auc_value_base <- numeric(n_bootstrap)
  acc_rep_base <- numeric(n_bootstrap)
  true_list_base <- list()  # To store true_sw for each iteration
  pred_list_base <- list() 
  ind_base <- list()
  
  for (i in 1:n_bootstrap) {
    n_subjects  = length(unique(features_sample$subject))
    n_samples  = length(unique(features_sample$time))

    
    #### Row wise and subject wise ####
    #  Prepare training and testing sets
    if(cv == "record-wise"){
      samples_test <- sample(x =nrow(features_sample), size = floor(testsize * nrow(features_sample))) 
      samples_train <- setdiff(1:(nrow(features_sample)), samples_test) 
      
      train_X = features_sample[samples_train, which(colnames(features_sample) %in% n_features)]
      test_X = features_sample[samples_test, which(colnames(features_sample) %in% n_features)]
      
      train_Y =  as.factor(features_sample$y[samples_train])
      test_Y = as.factor(features_sample$y[samples_test])
      
      subject[[i]] = features_sample$subject[samples_test]
      
      
    }else if(cv == "subject-wise"){
      subjects_test <- sample(n_subjects -round(testsize*n_subjects) , round(testsize*n_subjects)) 
      subjects_train <- setdiff(1:n_subjects, subjects_test)  
      
      train_X = features_sample[features_sample$subject %in% subjects_train,  which(colnames(features_sample) %in% n_features)]
      test_X = features_sample[features_sample$subject %in% subjects_test, which(colnames(features_sample) %in% n_features)]
      
      train_Y = as.factor(features_sample$y[features_sample$subject %in% subjects_train])
      test_Y = as.factor(features_sample$y[features_sample$subject %in% subjects_test])
      
      subject[[i]] = features_sample$subject[features_sample$subject %in% subjects_test]
    }  
    
    # Train and evaluate the Random Forest model  
    if(cv == "subject-wise" | cv == "record-wise"){
      
      
      RF <- randomForest(na.roughfix(train_X),train_Y)
      class_pred <- predict(RF, na.roughfix(test_X),type="response")
      class_pred_auc <- predict(RF, na.roughfix(test_X),type="prob")[, 2]
      acc[i] <- mean(as.numeric(as.character(class_pred)) == test_Y)
      
      true_list[[i]]<- test_Y
      pred_list[[i]] <- class_pred_auc
      roc_curve <- pROC::roc(test_Y,  class_pred_auc,quiet = TRUE,direction = "<",positive = "1", levels = c(0,1))
      auc_value[i] <- pROC::auc(roc_curve)
      
      ind[[i]] = data.frame(subject = subject[[i]], true = true_list[[i]], pred = pred_list[[i]])
      
      # Baseline model for record-wise CV
      if(cv == "record-wise"){ 
        Baseline <- glmer(
          train_Y ~ 1 + (1 | subject),
          data = features_sample[samples_train,], 
          family = binomial(link = "logit") 
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject), type = "response")
        #class_pred_base_auc <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject), type = "prob")[, 2]
        
        true_list_base[[i]]<- test_Y
        pred_list_base[[i]] <- as.numeric(as.character(class_pred_base))
        roc_curve_base <- pROC::roc(test_Y,  class_pred_base, quiet = TRUE,direction = "<")
        auc_value_base[i] <- pROC::auc(roc_curve_base)
        
        class_pred_base2 <- ifelse(class_pred_base > 0.5, 1, 0)
        acc_base[i] <- mean(as.numeric(as.character(class_pred_base2)) == test_Y)
        
        ind_base[[i]] = data.frame(subject = subject[[i]], true = true_list_base[[i]], pred = pred_list_base[[i]])
      }
    }
  } # end bootstrap
  
  if(cv == "record-wise"){ 

  }
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      filter(length(unique(true)) > 1) %>%
      filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
      summarise(
        auc_val = pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)),direction = "<", quiet = TRUE))[1],
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
    acc = acc,
    auc_value = auc_value,
    ind = ind,
    auc_value_base = auc_value_base,
    overall_summary = overall_summary
  )
  results_shiny <- list(
    auc_value_base = auc_value_base,
    acc_base = acc_base,
    Mean_AUC = mean(auc_value),
    Mean_Accuracy = mean(acc, na.rm = TRUE),
    overall_summary = overall_summary
  #  ind = ind
  )
  
  results_shiny <- t(as.data.frame(results_shiny))
  colnames(results_shiny) <- NULL
  rownames(results_shiny) <- c("AUC random intercept only:", "Accuracy random intercept only:", "AUC:", "Accuracy:", "Mean AUC within-person:", "SD AUC within-person:", "% of AUC > 0.5 within-person:", "N included within-person:")
  
  return(results)
} 


#features_sample <- read.csv("/Users/f007qrc/projects/Justintime_Paper/all_studentlife.csv")
############################################
run_simulation_centering_own <- function(features_sample,cv,n_bootstrap,testsize, seed = "12361488",n_features){
  set.seed(seed)
  acc <- numeric(n_bootstrap)
  auc_value <- numeric(n_bootstrap)
  acc_rep <- numeric(n_bootstrap)
  true_list <- list()  # To store true_sw for each iteration
  pred_list <- list() 
  subject <- list() 
  ind <- list()
  
  acc_base <- numeric(n_bootstrap)
  auc_value_base <- numeric(n_bootstrap)
  acc_rep_base <- numeric(n_bootstrap)
  true_list_base <- list()  # To store true_sw for each iteration
  pred_list_base <- list() 
  ind_base <- list()
  

  
  for (i in 1:n_bootstrap) {
    #### Row wise and subject wise ####
    #  Prepare training and testing sets
    
    samples_test <- sample(x = nrow(features_sample), size = floor(testsize * nrow(features_sample))) 
    samples_train <- setdiff(1:(nrow(features_sample)), samples_test) 
    
    train_X = features_sample[samples_train, which(colnames(features_sample) %in% n_features)]
    test_X = features_sample[samples_test, which(colnames(features_sample) %in% n_features)]
    
    train_Y =  as.factor(features_sample$y[samples_train])
    test_Y = as.factor(features_sample$y[samples_test])
    
    # train_Y =  as.factor(features_sample$stress[samples_train])
    # test_Y = as.factor(features_sample$stress[samples_test])
    
    subject[[i]] = features_sample$subject[samples_test]
    
    #subject[[i]] = features_sample$uid[samples_test]
    
    # Center the training set
    subject_means <- features_sample[samples_train, ] %>%
      group_by(subject) %>%
      summarise(across(all_of(names(train_X)), ~mean(.x, na.rm = TRUE)), .groups = "drop")
    
    # Center the training set
    train_X <- features_sample[samples_train, ] %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(train_X)[-1]), ~ . - get(paste0(cur_column(), "_mean")))) %>%
      select(all_of(names(train_X)))
    
    # Center the test set using training set means
    test_X <- features_sample[samples_test, ] %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(train_X)[-1]), ~ . - get(paste0(cur_column(), "_mean")))) %>%
      select(all_of(names(train_X)))
    
    # Train and evaluate the Random Forest model  
    if(cv == "record-wise"){
      RF <- randomForest(na.roughfix(train_X),train_Y)
      class_pred <- predict(RF, na.roughfix(test_X))
      class_pred_auc <- predict(RF, na.roughfix(test_X), type = "prob")[,2]
      
      acc[i] <- mean(as.numeric(as.character(class_pred)) == test_Y)
      
      true_list[[i]]<- test_Y
      pred_list[[i]] <- as.numeric(as.character(class_pred_auc))
      roc_curve <- pROC::roc(test_Y,  class_pred_auc,quiet = TRUE,direction = "<")
      auc_value[i] <- pROC::auc(roc_curve)
      
      ind[[i]] = data.frame(subject = subject[[i]], true = true_list[[i]], pred = pred_list[[i]])
      
      # Baseline model for record-wise CV
      if(cv == "record-wise"){ 
        Baseline <- glmer(
          train_Y ~ 1 + (1 | subject),
          data = features_sample[samples_train,], 
          family = binomial(link = "logit") 
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject), type = "response")
        #class_pred_base_auc <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject), type = "prob")[,2]
        
        true_list_base[[i]]<- test_Y
        pred_list_base[[i]] <- as.numeric(as.character(class_pred_base))
        roc_curve_base <- pROC::roc(test_Y,  class_pred_base,quiet = TRUE,direction = "<")
        auc_value_base[i] <- pROC::auc(roc_curve_base)
        
        class_pred_base2 <- ifelse(class_pred_base > 0.5, 1, 0)
        acc_base[i] <- mean(as.numeric(as.character(class_pred_base2)) == test_Y)
        
        ind_base[[i]] = data.frame(subject = subject[[i]], true = true_list_base[[i]], pred = pred_list_base[[i]])
      }
    }
  } # end bootstrap

  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      filter(length(unique(true)) > 1) %>%
      filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
      summarise(
        auc_val = pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)),direction = "<", quiet = TRUE))[1],
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
  
  results_shiny <- list(
    auc_value_base = auc_value_base,
    acc_base = acc_base,
    Mean_AUC = mean(auc_value),
    Mean_Accuracy = mean(acc, na.rm = TRUE),
    overall_summary = overall_summary
  )
  
  results_shiny <- t(as.data.frame(results_shiny))
  colnames(results_shiny) <- NULL
  rownames(results_shiny) <- c("AUC random intercept only:", "Accuracy random intercept only:", "AUC:", "Accuracy:", "Mean AUC within-person:", "SD AUC within-person:", "% of AUC > 0.5 within-person:", "N included within-person:")
  

  return(results)
}      




############# Sliding Window ################

SlidingWindow_CV <- function(data, origin, horizon){
  samplesize <- data
  trainindex <- list()
  testindex <- list()
  
  trainindex[[1]] <- 1:origin
  testindex[[1]] <- (origin + 1):(origin+horizon)
  
  counter <- 1
  index <- testindex[[1]][horizon]+1
  
  while(testindex[[counter]][horizon] < (samplesize-horizon)){
    index <- testindex[[counter]][horizon]+1
    trainindex[[counter+1]] <- (index-origin+1):index-1
    testindex[[counter+1]] <- (index):(index+horizon-1)
    counter <- counter + 1
  }
  return(list(trainindex,testindex))
}


#############################

run_simulation_slidingwindow_own <- function(features_sample,n_bootstrap,windowsize,n_features = n_features,seed = "1234"){
  set.seed(seed) 
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
  
  n_samples  = length(unique(features_sample$time))

  for (i in 1:n_bootstrap) {
    #  Prepare training and testing sets
    timeSlices <- SlidingWindow_CV(n_samples,windowsize,1) #window size, prediction horizon
    #timeSlices <- createTimeSlices(seq(0,n_samples),14,1,FALSE)
    trainSlices <- timeSlices[[1]]
    testSlices <- timeSlices[[2]]
    
    lengthtest <- length(unlist(testSlices))
    
    acc_sw <- numeric(lengthtest)
    auc_value_sw <- numeric(lengthtest)
    acc_sw_base <- numeric(lengthtest)
    auc_value_base_sw <- numeric(lengthtest)
    ind[[i]] = data.frame(subject = as.character(), true = as.numeric(), pred = as.numeric())
    ind_t <- data.frame()
    ind_t_base <- data.frame()
    
    for(k in 1:length(trainSlices)){
      tryCatch({
      train_X = features_sample[features_sample$time %in% trainSlices[[k]], which(colnames(features_sample) %in% n_features)]
      test_X = features_sample[features_sample$time %in% testSlices[[k]], which(colnames(features_sample) %in% n_features)]
      train_Y = as.factor(features_sample$y[features_sample$time %in% trainSlices[[k]]])
      test_Y = as.factor(features_sample$y[features_sample$time %in% testSlices[[k]]])
      
      # Train and evaluate the Random Forest model 
      
      RF <- randomForest(na.roughfix(train_X),train_Y)
      class_pred <- predict(RF, na.roughfix(test_X))
      class_pred_auc <- predict(RF, na.roughfix(test_X), type = "prob")[,2]
      
      
      acc_sw[k] <- mean(as.numeric(as.character(class_pred)) == test_Y)
      
      if (length(unique(test_Y)) == 2) {
        roc_curve <- pROC::roc(test_Y, class_pred_auc, quiet = TRUE,direction = "<")
        auc_value_sw[k] <- pROC::auc(roc_curve)
      } else {
        auc_value_sw[k] <- NA
      }
      
      ind_t <- rbind(
        ind_t,
        data.frame(
          subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
          true = test_Y,                                   
          pred = class_pred_auc
        )
      )
      
      #### Baseline
      Baseline <- glmer(
        "y ~ 1 + (1|subject)", # Random intercept for 'subject'
        data = features_sample[features_sample$time %in% trainSlices[[k]],], 
        family = "binomial" # Logistic regression
      )
      
      class_pred_base <- predict(Baseline, newdata = features_sample[features_sample$time %in% testSlices[[k]],], re.form = ~(1 | subject), type = "response")
      #class_pred_base_auc <- predict(Baseline, newdata = features_sample[features_sample$time %in% testSlices[[k]],], re.form = ~(1 | subject), type = "prob")[,2]
      
      if (length(unique(test_Y)) == 2) {
        roc_curve_base <- pROC::roc(test_Y,  class_pred_base,quiet = TRUE,direction = "<")
        auc_value_base_sw[k] <- pROC::auc(roc_curve_base)
      } else {
        auc_value_base_sw[k] <- NA
      }
      
      class_pred_base2 <- ifelse(class_pred_base > 0.5, 1, 0)
      acc_sw_base[k] <- mean(as.numeric(as.character(class_pred_base2)) == test_Y)
      
      
      ind_t_base <- rbind(
        ind_t_base,
        data.frame(
          subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
          true = test_Y,                                   
          pred = class_pred_base
        )
      )
      }, error = function(e) {
        message(paste("Skipping iteration", k, "due to error:", e$message))
      })
    }

    auc_value[i] <- pROC::auc(pROC::roc(as.numeric(as.character(ind_t$true)), as.numeric(as.character(ind_t$pred)),direction = "<", quiet = TRUE))
    auc_value_meansw[i] <- mean(auc_value_sw,na.rm = TRUE)
    ind[[i]] <- ind_t
    ind_base[[i]] <- ind_t_base
    auc_value_base[i] <- pROC::auc(pROC::roc(as.numeric(as.character(ind_t_base$true)), as.numeric(as.character(ind_t_base$pred)),direction = "<", quiet = TRUE))
    acc_sw <- mean(acc_sw, na.rm = TRUE)
    acc_base <- mean(acc_sw_base, na.rm = TRUE)
  } # end bootstrap
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      filter(length(unique(true)) > 1) %>%
      filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
      summarise(
        auc_val = pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)),direction = "<", quiet = TRUE))[1],
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
  
  results_summary_base <- lapply(ind_base, function(ind_result_base) {
    processed_base <- ind_result_base %>%
      group_by(subject) %>%
      filter(length(unique(true)) > 1) %>%
      filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
      summarise(
        auc_val = pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)),direction = "<", quiet = TRUE))[1],
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_auc_base <- mean(processed_base$auc_val, na.rm = TRUE)
    sd_auc_base <- sd(processed_base$auc_val, na.rm = TRUE)
    percent_above_0_5_base <- (sum(processed_base$auc_val > 0.5, na.rm = TRUE) / nrow(processed_base)) * 100
    n_base <- nrow(processed_base) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_auc_base = mean_auc_base,
      sd_auc_base = sd_auc_base,
      percent_above_0_5_base = percent_above_0_5_base,
      total_n_base = n_base
    )
  })
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean = x$mean_auc, sd = x$sd_auc, percent_above_0_5 = x$percent_above_0_5, total_n = x$total_n)
  }))
  
  overall_summary_base <- do.call(rbind, lapply(results_summary_base, function(x) {
    data.frame(mean_auc = x$mean_auc_base, sd = x$sd_auc_base, percent_above_0_5 = x$percent_above_0_5_base, total_n = x$total_n_base)
  }))
  
  results <- list(
    auc_value = auc_value,
    ind = ind,
    ind_base = ind_base,
    overall_summary = overall_summary,
    overall_summary_base = overall_summary_base
  )
  
  results_shiny <- list(
    auc_value_base = mean(auc_value_base, na.rm = TRUE),
    acc_base = mean(acc_base, na.rm = T),
    Mean_AUC = mean(auc_value, na.rm = T),
    Mean_Accuracy = mean(acc_sw, na.rm = TRUE),
    overall_summary = overall_summary,
    overall_summary_base = overall_summary_base
   # ind = ind
  )
  
  results_shiny <- t(as.data.frame(results_shiny))
  colnames(results_shiny) <- NULL
  rownames(results_shiny) <- c("AUC random intercept only:", "Accuracy random intercept only:", "AUC:", "Accuracy:", "Mean AUC within-person:", "SD AUC within-person:", "% of AUC > 0.5 within-person:", "N included within-person:", "Mean AUC random intercept only:","SD AUC within-person base:", "% of AUC > 0.5 within-person base:", "N included within-person base:")
  
  
  return(results)
}      


run_simulation_slidingwindow_own_centering <- function(features_sample,n_bootstrap,windowsize,n_features = n_features, seed = "1234"){
 set.seed(seed)
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
  
  n_samples  = length(unique(features_sample$time))
  

  for (i in 1:n_bootstrap) {
    #  Prepare training and testing sets
    timeSlices <- SlidingWindow_CV(n_samples,windowsize,1) #window size, prediction horizon
    #timeSlices <- createTimeSlices(seq(0,n_samples),14,1,FALSE)
    trainSlices <- timeSlices[[1]]
    testSlices <- timeSlices[[2]]
    
    lengthtest <- length(unlist(testSlices))
    
    acc_sw <- numeric(lengthtest)
    auc_value_sw <- numeric(lengthtest)
    acc_sw_base <- numeric(lengthtest)
    auc_value_base_sw <- numeric(lengthtest)
    ind[[i]] = data.frame(subject = as.character(), true = as.numeric(), pred = as.numeric())
    ind_t <- data.frame()
    ind_t_base <- data.frame()
    
    for (k in 1:length(trainSlices)) {
      tryCatch({
        train_X <- features_sample[features_sample$time %in% trainSlices[[k]], which(colnames(features_sample) %in% n_features)]
        test_X <- features_sample[features_sample$time %in% testSlices[[k]], which(colnames(features_sample) %in% n_features)]
        train_Y <- as.factor(features_sample$y[features_sample$time %in% trainSlices[[k]]])
        test_Y <- as.factor(features_sample$y[features_sample$time %in% testSlices[[k]]])
        
        subject_means <- features_sample[features_sample$time %in% trainSlices[[k]],] %>%
          group_by(subject) %>%
          summarise(across(all_of(names(train_X)), mean), .groups = "drop")
        
        train_X <- features_sample[features_sample$time %in% trainSlices[[k]],] %>%
          left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
          mutate(across(all_of(names(train_X)), ~ . - get(paste0(cur_column(), "_mean")))) %>%
          select(all_of(names(train_X)))
        
        test_X <- features_sample[features_sample$time %in% testSlices[[k]],] %>%
          left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
          mutate(across(all_of(names(train_X)), ~ . - get(paste0(cur_column(), "_mean")))) %>%
          select(all_of(names(train_X)))
        
        RF <- randomForest(na.roughfix(train_X), train_Y)
        class_pred <- predict(RF, na.roughfix(test_X))
        acc_sw[k] <- mean(as.numeric(as.character(class_pred)) == test_Y)
        
        if (length(unique(test_Y)) == 2) {
          roc_curve <- pROC::roc(test_Y, as.numeric(as.character(class_pred)), quiet = TRUE,direction = "<")
          auc_value_sw[k] <- pROC::auc(roc_curve)
        } else {
          auc_value_sw[k] <- NA
        }
        
        ind_t <- rbind(
          ind_t,
          data.frame(
            subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
            true = test_Y,
            pred = class_pred
          )
        )
        
        Baseline <- glmer(
          "y ~ 1 + (1|subject)",
          data = features_sample[features_sample$time %in% trainSlices[[k]],],
          family = "binomial"
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[features_sample$time %in% testSlices[[k]],], re.form = ~(1 | subject), type = "response")
        
        if (length(unique(test_Y)) == 2) {
          roc_curve_base <- pROC::roc(test_Y, as.numeric(as.character(class_pred_base)), quiet = TRUE,direction = "<")
          auc_value_base_sw[k] <- pROC::auc(roc_curve_base)
        } else {
          auc_value_base_sw[k] <- NA
        }
        
        class_pred_base2 <- ifelse(class_pred_base > 0.5, 1, 0)
        acc_sw_base[k] <- mean(as.numeric(as.character(class_pred_base2)) == test_Y)
        
        ind_t_base <- rbind(
          ind_t_base,
          data.frame(
            subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
            true = test_Y,
            pred = class_pred_base
          )
        )
      }, error = function(e) {
        message(paste("Skipping iteration", k, "due to error:", e$message))
      })
    }
    
    
    auc_value[i] <- pROC::auc(pROC::roc(as.numeric(as.character(ind_t$true)), as.numeric(as.character(ind_t$pred)),direction = "<", quiet = TRUE))
    auc_value_meansw[i] <- mean(auc_value_sw,na.rm = TRUE)
    ind[[i]] <- ind_t
    auc_value_base[i] <- pROC::auc(pROC::roc(as.numeric(as.character(ind_t_base$true)), as.numeric(as.character(ind_t_base$pred)),direction = "<", quiet = TRUE))
    acc_sw <- mean(acc_sw, na.rm = TRUE)
    acc_base <- mean(acc_sw_base, na.rm = TRUE)
  } # end bootstrap
  
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      filter(length(unique(true)) > 1) %>%
      filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
      summarise(
        auc_val = pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)),direction = "<", quiet = TRUE))[1],
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
  
  results_summary_base <- lapply(ind_base, function(ind_result_base) {
    processed_base <- ind_result_base %>%
      group_by(subject) %>%
      filter(length(unique(true)) > 1) %>%
      filter(sum(true == 1) > 0, sum(true == 0) > 0) %>%
      summarise(
        auc_val = pROC::auc(pROC::roc(as.numeric(as.character(true)), as.numeric(as.character(pred)),direction = "<", quiet = TRUE))[1],
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_auc_base <- mean(processed$auc_val, na.rm = TRUE)
    sd_auc_base <- sd(processed$auc_val, na.rm = TRUE)
    percent_above_0_5_base <- (sum(processed$auc_val > 0.5, na.rm = TRUE) / nrow(processed)) * 100
    n_base <- nrow(processed) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_auc_base = mean_auc_base,
      sd_auc_base = sd_auc_base,
      percent_above_0_5_base = percent_above_0_5_base,
      total_n_base = n_base,
    )
  })
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean = x$mean_auc, sd = x$sd_auc, percent_above_0_5 = x$percent_above_0_5, total_n = x$total_n)
  }))
  
  overall_summary_base <- do.call(rbind, lapply(results_summary_base, function(x) {
    data.frame(mean_auc = x$mean_auc_base, sd = x$sd_auc_base, percent_above_0_5 = x$percent_above_0_5_base, total_n = x$total_n_base)
  }))
  
  results <- list(
    auc_value = auc_value,
    ind = ind,
    overall_summary = overall_summary,
    overall_summary_base = overall_summary_base
  )
  
  results_shiny <- list(
    auc_value_base = mean(auc_value_base, na.rm = TRUE),
    acc_base = mean(acc_base, na.rm = T),
    Mean_AUC = mean(auc_value, na.rm = T),
    Mean_Accuracy = mean(acc_sw, na.rm = TRUE),
    overall_summary = overall_summary,
    overall_summary_base = overall_summary_base
  )
  
  results_shiny <- t(as.data.frame(results_shiny))
  colnames(results_shiny) <- NULL
  rownames(results_shiny) <- c("AUC random intercept only:", "Accuracy random intercept only:", "AUC:", "Accuracy:", "Mean AUC within-person:", "SD AUC within-person:", "% of AUC > 0.5 within-person:", "N included within-person:", "Mean AUC random intercept only:","SD AUC within-person base:", "% of AUC > 0.5 within-person base:", "N included within-person base:")
  
  
  return(results)
}      

