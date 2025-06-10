library("Metrics")
library("TSrepr") #mape

r2 <- function(true,pred){
  rss <- sum((pred- true) ^ 2,na.rm = TRUE)  # residual sum of squares
  tss <- sum((true- mean(true, na.rm = TRUE))^ 2,na.rm = TRUE)  # total sum of squares
  rsq <- 1 - rss/tss
  return(rsq)
}

##############################
### Simulation Functions #####
##############################

shuffle_data <- function(data, subject_var, outcome_var) {
  # Shuffle participant IDs:
  unique_subjects <- unique(data[[subject_var]])
  shuffled_subjects <- sample(unique_subjects) # shuffles order of subjects without replacement
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
  
  cor <- numeric(n_bootstrap)
  rsq <- numeric(n_bootstrap)
  mape <- numeric(n_bootstrap)
  mae <- numeric(n_bootstrap)
  rae <- numeric(n_bootstrap)
 
  true_list <- list()  # To store true_sw for each iteration
  pred_list <- list() 
  subject <- list() 
  ind <- list()
  
  cor_base <- numeric(n_bootstrap)
  rsq_base <- numeric(n_bootstrap)
  mape_base <- numeric(n_bootstrap)
  mae_base <- numeric(n_bootstrap)
  rae_base <- numeric(n_bootstrap)
  
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
      
      train_Y = features_sample$y[samples_train]
      test_Y = features_sample$y[samples_test]
      
      subject[[i]] = features_sample$subject[samples_test]
      
      
    }else if(cv == "subject-wise"){
      subjects_test <- sample(n_subjects -round(testsize*n_subjects) , round(testsize*n_subjects)) 
      subjects_train <- setdiff(1:n_subjects, subjects_test)  
      
      train_X = features_sample[features_sample$subject %in% subjects_train,  which(colnames(features_sample) %in% n_features)]
      test_X = features_sample[features_sample$subject %in% subjects_test, which(colnames(features_sample) %in% n_features)]
      
      train_Y = features_sample$y[features_sample$subject %in% subjects_train]
      test_Y = features_sample$y[features_sample$subject %in% subjects_test]
      
      subject[[i]] = features_sample$subject[features_sample$subject %in% subjects_test]
    }  
    
    # Train and evaluate the Random Forest model  
    if(cv == "subject-wise" | cv == "record-wise"){
      
      
      RF <- randomForest(na.roughfix(train_X),train_Y)
      class_pred <- predict(RF, na.roughfix(test_X))

      rsq[i] <- r2(test_Y,class_pred)
      cor[i] <- cor(test_Y,class_pred, use = "pairwise.complete.obs")
      mape[i] <- mape(test_Y,class_pred)
      mae[i] <- mae(test_Y,class_pred)
      rae[i] <- rae(test_Y,class_pred)
      
      true_list[[i]]<- test_Y
      pred_list[[i]] <- class_pred

      ind[[i]] = data.frame(subject = subject[[i]], true = true_list[[i]], pred = pred_list[[i]])
      
      # Baseline model for record-wise CV
      if(cv == "record-wise"){ 
        Baseline <- lmer(
          train_Y ~ 1 + (1 | subject),
          data = features_sample[samples_train,]
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject))
        
        rsq_base[i] <- r2(test_Y,class_pred_base)
        cor_base[i] <- cor(test_Y,class_pred_base, use = "pairwise.complete.obs")
        mape_base[i] <- mape(test_Y,class_pred_base)
        mae_base[i] <- mae(test_Y,class_pred_base)
        rae_base[i] <- rae(test_Y,class_pred_base)
        
        true_list_base[[i]]<- test_Y
        pred_list_base[[i]] <- class_pred_base
        
        ind_base[[i]] = data.frame(subject = subject[[i]], true = true_list_base[[i]], pred = pred_list_base[[i]])
      }
    }
  } # end bootstrap
  
  if(cv == "record-wise"){ 
    
  }
  
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      summarise(
        r2_value = r2(true, pred),
        correlation =  cor(true, pred),
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_r2 <- mean(processed$r2_value, na.rm = TRUE)
    mean_cor <- mean(processed$correlation, na.rm = TRUE)
    sd_r2 <- sd(processed$r2_value, na.rm = TRUE)
    n <- nrow(processed) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_r2 = mean_r2,
      mean_cor = mean_cor,
      sd_r2 = sd_r2,
      total_n = n,
      bootstrap_results = processed # Include detailed bootstrap-level data
    )
  })
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean_r2 = x$mean_r2,mean_cor = x$mean_cor,  sd = x$sd_r2, total_n = x$total_n)
  }))
  
  results <- list(
    rsq = rsq,
    cor = cor,
    mape = mape,
    mae = mae,
    rae = rae,
    #   ind = ind,
    rsq_base = rsq_base,
    cor_base = cor_base,
    mape_base = mape_base,
    mae_base = mae_base,
    rae_base = rae_base,
    overall_summary = overall_summary
  )
  
  results <- as.data.frame(results)
  return(results)
} 


############################################
run_simulation_centering_own <- function(features_sample,cv,n_bootstrap,testsize, seed = "12361488",n_features){
  set.seed(seed)

  cor <- numeric(n_bootstrap)
  rsq <- numeric(n_bootstrap)
  mape <- numeric(n_bootstrap)
  mae <- numeric(n_bootstrap)
  rae <- numeric(n_bootstrap)
  
  true_list <- list()  # To store true_sw for each iteration
  pred_list <- list() 
  subject <- list() 
  ind <- list()
  
  cor_base <- numeric(n_bootstrap)
  rsq_base <- numeric(n_bootstrap)
  mape_base <- numeric(n_bootstrap)
  mae_base <- numeric(n_bootstrap)
  rae_base <- numeric(n_bootstrap)
  
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
    
    train_Y =  features_sample$y[samples_train]
    test_Y = features_sample$y[samples_test]
    
    subject[[i]] = features_sample$subject[samples_test]
    
    #subject[[i]] = features_sample$uid[samples_test]
    
    # Center the training set
    subject_means <- features_sample[samples_train, ] %>%
      group_by(subject) %>%
      summarise(across(all_of(names(train_X)), ~mean(.x, na.rm = TRUE)), .groups = "drop")
    
    # Center the training set
    train_X <- features_sample[samples_train, ] %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(train_X[-1])), ~ . - get(paste0(cur_column(), "_mean")))) %>%
      select(all_of(names(train_X)))
    
    # Center the test set using training set means
    test_X<- features_sample[samples_test, ] %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(train_X[-1])), ~ . - get(paste0(cur_column(), "_mean")))) %>%
      select(all_of(names(train_X)))
    
    # Train and evaluate the Random Forest model  
    if(cv == "record-wise"){
      RF <- randomForest(na.roughfix(train_X),train_Y)
      class_pred <- predict(RF, na.roughfix(test_X))
      
      rsq[i] <- r2(test_Y,class_pred)
      cor[i] <- cor(test_Y,class_pred, use = "pairwise.complete.obs")
      mape[i] <- mape(test_Y,class_pred)
      mae[i] <- mae(test_Y,class_pred)
      rae[i] <- rae(test_Y,class_pred)
      
      true_list[[i]]<- test_Y
      pred_list[[i]] <- class_pred
      
      ind[[i]] = data.frame(subject = subject[[i]], true = true_list[[i]], pred = pred_list[[i]])
      
      
      # Baseline model for record-wise CV
      if(cv == "record-wise"){ 
        Baseline <- lmer(
          train_Y ~ 1 + (1 | subject),
          data = features_sample[samples_train,]
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject))
        
        rsq_base[i] <- r2(test_Y,class_pred_base)
        cor_base[i] <- cor(test_Y,class_pred_base, use = "pairwise.complete.obs")
        mape_base[i] <- mape(test_Y,class_pred_base)
        mae_base[i] <- mae(test_Y,class_pred_base)
        rae_base[i] <- rae(test_Y,class_pred_base)
        
        true_list_base[[i]]<- test_Y
        pred_list_base[[i]] <- class_pred_base
        
        ind_base[[i]] = data.frame(subject = subject[[i]], true = true_list_base[[i]], pred = pred_list_base[[i]])
        
      }
    }
  } # end bootstrap
  
  
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      summarise(
        r2_value = r2(true, pred),
        correlation =  cor(true, pred),
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_r2 <- mean(processed$r2_value, na.rm = TRUE)
    mean_cor <- mean(processed$correlation, na.rm = TRUE)
    sd_r2 <- sd(processed$r2_value, na.rm = TRUE)
    n <- nrow(processed) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_r2 = mean_r2,
      mean_cor = mean_cor,
      sd_r2 = sd_r2,
      total_n = n,
      bootstrap_results = processed # Include detailed bootstrap-level data
    )
  })
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean_r2 = x$mean_r2,mean_cor = x$mean_cor,  sd = x$sd_r2, total_n = x$total_n)
  }))
  
  results <- list(
    rsq = rsq,
    cor = cor,
    mape = mape,
    mae = mae,
    rae = rae,
    #   ind = ind,
    rsq_base = rsq_base,
    cor_base = cor_base,
    mape_base = mape_base,
    mae_base = mae_base,
    rae_base = rae_base,
    overall_summary = overall_summary
  )
  
  results <- as.data.frame(results)
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

run_simulation_slidingwindow_own <- function(features_sample,n_bootstrap,windowsize,n_features = n_features, seed = "1234"){
  set.seed(seed)
  cor <- numeric(n_bootstrap)
  rsq <- numeric(n_bootstrap)
  mape <- numeric(n_bootstrap)
  mae <- numeric(n_bootstrap)
  rae <- numeric(n_bootstrap)
  
  true_list <- list()  # To store true_sw for each iteration
  pred_list <- list() 
  subject <- list() 
  ind <- list()
  
  cor_base <- numeric(n_bootstrap)
  rsq_base <- numeric(n_bootstrap)
  mape_base <- numeric(n_bootstrap)
  mae_base <- numeric(n_bootstrap)
  rae_base <- numeric(n_bootstrap)
  
  true_list_base <- list()  # To store true_sw for each iteration
  pred_list_base <- list() 
  ind_base <- list()
  
  
  r2_value_meansw <- numeric(n_bootstrap)
  
  
  n_samples  = length(unique(features_sample$time))
  
  for (i in 1:n_bootstrap) {
    #  Prepare training and testing sets
    timeSlices <- SlidingWindow_CV(n_samples,windowsize,1) #window size, prediction horizon
    #timeSlices <- createTimeSlices(seq(0,n_samples),14,1,FALSE)
    trainSlices <- timeSlices[[1]]
    testSlices <- timeSlices[[2]]
    
    lengthtest <- length(unlist(testSlices))
    
    cor_sw <- numeric(lengthtest)
    rsq_sw <- numeric(lengthtest)
    mape_sw <- numeric(lengthtest)
    mae_sw <- numeric(lengthtest)
    rae_sw <- numeric(lengthtest)
    cor_base_sw <- numeric(lengthtest)
    rsq_base_sw <- numeric(lengthtest)
    mape_base_sw <- numeric(lengthtest)
    mae_base_sw <- numeric(lengthtest)
    rae_base_sw <- numeric(lengthtest)
    
    ind[[i]] = data.frame(subject = as.character(), true = as.numeric(), pred = as.numeric())
    ind_t <- data.frame()
    ind_t_base <- data.frame()
    
    for(k in 1:length(trainSlices)){
      tryCatch({
        train_X = features_sample[features_sample$time %in% trainSlices[[k]], which(colnames(features_sample) %in% n_features)]
        test_X = features_sample[features_sample$time %in% testSlices[[k]], which(colnames(features_sample) %in% n_features)]
        train_Y = features_sample$y[features_sample$time %in% trainSlices[[k]]]
        test_Y = features_sample$y[features_sample$time %in% testSlices[[k]]]
        
        # Train and evaluate the Random Forest model 
        
        RF <- randomForest(na.roughfix(train_X),train_Y)
        class_pred <- predict(RF, na.roughfix(test_X))
        
        rsq_sw[k] <- r2(test_Y,class_pred)
        cor_sw[k] <- cor(test_Y,class_pred, use = "pairwise.complete.obs")
        mape_sw[k] <- mape(test_Y,class_pred)
        mae_sw[k] <- mae(test_Y,class_pred)
        rae_sw[k] <- rae(test_Y,class_pred)
        
        ind_t <- rbind(
          ind_t,
          data.frame(
            subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
            true = test_Y,                                   
            pred = class_pred
          )
        )
        
        #### Baseline
        Baseline <- lmer(
          "y ~ 1 + (1|subject)", # Random intercept for 'subject'
          data = features_sample[features_sample$time %in% trainSlices[[k]],]
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[features_sample$time %in% testSlices[[k]],], re.form = ~(1 | subject))
        
        
        rsq_base_sw[k] <- r2(test_Y,class_pred_base)
        cor_base_sw[k] <- cor(test_Y,class_pred_base, use = "pairwise.complete.obs")
        mape_base_sw[k] <- mape(test_Y,class_pred_base)
        mae_base_sw[k] <- mae(test_Y,class_pred_base)
        rae_base_sw[k] <- rae(test_Y,class_pred_base)
        
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
    ind[[i]] <- ind_t
    ind_base[[i]] <- ind_t_base
    
    rsq[i] <- r2(ind_t$true,ind_t$pred)
    mean_rsq_sw <- mean(rsq_sw,na.rm = TRUE)
    cor[i] <- cor(ind_t$true,ind_t$pred, use = "pairwise.complete.obs")
    mape[i] <- mape(ind_t$true,ind_t$pred)
    mae[i] <- mae(ind_t$true,ind_t$pred)
    rae[i] <- rae(ind_t$true,ind_t$pred)
    
    
    rsq_base[i] <- r2(ind_t_base$true,ind_t_base$pred)
    mean_rsq_base_sw <- mean(rsq_base_sw,na.rm = TRUE)
    cor_base[i] <- cor(ind_t_base$true,ind_t_base$pred, use = "pairwise.complete.obs")
    mape_base[i] <- mape(ind_t_base$true,ind_t_base$pred)
    mae_base[i] <- mae(ind_t_base$true,ind_t_base$pred)
    rae_base[i] <- rae(ind_t_base$true,ind_t_base$pred)
    
  } # end bootstrap
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      summarise(
        r2_value = r2(true, pred),
        correlation =  cor(true, pred),
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_r2 <- mean(processed$r2_value, na.rm = TRUE)
    mean_cor <- mean(processed$correlation, na.rm = TRUE)
    sd_r2 <- sd(processed$r2_value, na.rm = TRUE)
    n <- nrow(processed) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_r2 = mean_r2,
      mean_cor = mean_cor,
      sd_r2 = sd_r2,
      total_n = n,
      bootstrap_results = processed # Include detailed bootstrap-level data
    )
  })
  
  results_summary_base <- lapply(ind_base, function(ind_result_base) {
    processed_base <- ind_result_base %>%
      group_by(subject) %>%
      summarise(
        r2_value = r2(true, pred),
        correlation =  cor(true, pred),
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_r2_base <- mean(processed_base$r2_value, na.rm = TRUE)
    mean_cor_base <- mean(processed_base$correlation, na.rm = TRUE)
    sd_r2_base <- sd(processed_base$r2_value, na.rm = TRUE)
    n_base <- nrow(processed_base) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_r2_base = mean_r2_base,
      mean_cor_base = mean_cor_base,
      sd_r2_base = sd_r2_base,
      total_n_base = n_base,
      bootstrap_results_base = processed_base # Include detailed bootstrap-level data
    )
  })
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean_r2 = x$mean_r2,mean_cor = x$mean_cor,  sd = x$sd_r2, total_n = x$total_n)
  }))
  
  overall_summary_base <- do.call(rbind, lapply(results_summary_base, function(x) {
    data.frame(mean_r2_base = x$mean_r2_base, mean_cor_base = x$mean_cor_base,  sd_base = x$sd_r2_base, total_n_base = x$total_n_base)
  }))
  
  results <- list(
    rsq = rsq,
    cor = cor,
    mape = mape,
    mae = mae,
    rae = rae,
    #   ind = ind,
    rsq_base = rsq_base,
    cor_base = cor_base,
    mape_base = mape_base,
    mae_base = mae_base,
    rae_base = rae_base,
    overall_summary = overall_summary,
    overall_summary_base = overall_summary_base
  )
  
  results <- as.data.frame(results)
}      


run_simulation_slidingwindow_own_centering <- function(features_sample,n_bootstrap,windowsize,n_features = n_features, seed = "1234"){
 set.seed(1234)
  cor <- numeric(n_bootstrap)
  rsq <- numeric(n_bootstrap)
  mape <- numeric(n_bootstrap)
  mae <- numeric(n_bootstrap)
  rae <- numeric(n_bootstrap)
  
  true_list <- list()  # To store true_sw for each iteration
  pred_list <- list() 
  subject <- list() 
  ind <- list()
  
  cor_base <- numeric(n_bootstrap)
  rsq_base <- numeric(n_bootstrap)
  mape_base <- numeric(n_bootstrap)
  mae_base <- numeric(n_bootstrap)
  rae_base <- numeric(n_bootstrap)
  
  true_list_base <- list()  # To store true_sw for each iteration
  pred_list_base <- list() 
  ind_base <- list()
  
  r2_value_meansw <- numeric(n_bootstrap)
  n_samples  = length(unique(features_sample$time))
  
  
  for (i in 1:n_bootstrap) {
    #  Prepare training and testing sets
    timeSlices <- SlidingWindow_CV(n_samples,windowsize,1) #window size, prediction horizon
    #timeSlices <- createTimeSlices(seq(0,n_samples),14,1,FALSE)
    trainSlices <- timeSlices[[1]]
    testSlices <- timeSlices[[2]]
    
    lengthtest <- length(unlist(testSlices))
    
    cor_sw <- numeric(lengthtest)
    rsq_sw <- numeric(lengthtest)
    mape_sw <- numeric(lengthtest)
    mae_sw <- numeric(lengthtest)
    rae_sw <- numeric(lengthtest)
    cor_base_sw <- numeric(lengthtest)
    rsq_base_sw <- numeric(lengthtest)
    mape_base_sw <- numeric(lengthtest)
    mae_base_sw <- numeric(lengthtest)
    rae_base_sw <- numeric(lengthtest)
    
    ind[[i]] = data.frame(subject = as.character(), true = as.numeric(), pred = as.numeric())
    ind_t <- data.frame()
    ind_t_base <- data.frame()
    
    for (k in 1:length(trainSlices)) {
      tryCatch({
        train_X <- features_sample[features_sample$time %in% trainSlices[[k]], which(colnames(features_sample) %in% n_features)]
        test_X <- features_sample[features_sample$time %in% testSlices[[k]], which(colnames(features_sample) %in% n_features)]
        train_Y <- features_sample$y[features_sample$time %in% trainSlices[[k]]]
        test_Y <- features_sample$y[features_sample$time %in% testSlices[[k]]]
      
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
        class_pred <- predict(RF, na.roughfix(test_X))
        
        rsq_sw[k] <- r2(test_Y,class_pred)
        cor_sw[k] <- cor(test_Y,class_pred, use = "pairwise.complete.obs")
        mape_sw[k] <- mape(test_Y,class_pred)
        mae_sw[k] <- mae(test_Y,class_pred)
        rae_sw[k] <- rae(test_Y,class_pred)
        
        ind_t <- rbind(
          ind_t,
          data.frame(
            subject = as.character(features_sample$subject[features_sample$time %in% testSlices[[k]]]),
            true = test_Y,                                   
            pred = class_pred
          )
        )
        
        # Baseline
        Baseline <- lmer(
          "y ~ 1 + (1|subject)",
          data = features_sample[features_sample$time %in% trainSlices[[k]],]
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[features_sample$time %in% testSlices[[k]],], re.form = ~(1 | subject))
        
        rsq_base_sw[k] <- r2(test_Y,class_pred_base)
        cor_base_sw[k] <- cor(test_Y,class_pred_base, use = "pairwise.complete.obs")
        mape_base_sw[k] <- mape(test_Y,class_pred_base)
        mae_base_sw[k] <- mae(test_Y,class_pred_base)
        rae_base_sw[k] <- rae(test_Y,class_pred_base)
        
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
    
    ind[[i]] <- ind_t
    ind_base[[i]] <- ind_t_base
    
    rsq[i] <- r2(ind_t$true,ind_t$pred)
    mean_rsq_sw <- mean(rsq_sw,na.rm = TRUE)
    cor[i] <- cor(ind_t$true,ind_t$pred, use = "pairwise.complete.obs")
    mape[i] <- mape(ind_t$true,ind_t$pred)
    mae[i] <- mae(ind_t$true,ind_t$pred)
    rae[i] <- rae(ind_t$true,ind_t$pred)
    
    
    rsq_base[i] <- r2(ind_t_base$true,ind_t_base$pred)
    mean_rsq_base_sw <- mean(rsq_base_sw,na.rm = TRUE)
    cor_base[i] <- cor(ind_t_base$true,ind_t_base$pred, use = "pairwise.complete.obs")
    mape_base[i] <- mape(ind_t_base$true,ind_t_base$pred)
    mae_base[i] <- mae(ind_t_base$true,ind_t_base$pred)
    rae_base[i] <- rae(ind_t_base$true,ind_t_base$pred)
  } # end bootstrap
  
  
  results_summary <- lapply(ind, function(ind_result) {
    processed <- ind_result %>%
      group_by(subject) %>%
      summarise(
        r2_value = r2(true, pred),
        correlation =  cor(true, pred),
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_r2 <- mean(processed$r2_value, na.rm = TRUE)
    mean_cor <- mean(processed$correlation, na.rm = TRUE)
    sd_r2 <- sd(processed$r2_value, na.rm = TRUE)
    n <- nrow(processed) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_r2 = mean_r2,
      mean_cor = mean_cor,
      sd_r2 = sd_r2,
      total_n = n,
      bootstrap_results = processed # Include detailed bootstrap-level data
    )
  })
  
  results_summary_base <- lapply(ind_base, function(ind_result_base) {
    processed_base <- ind_result_base %>%
      group_by(subject) %>%
      summarise(
        r2_value = r2(true, pred),
        correlation =  cor(true, pred),
        .groups = "drop" # Ungroup after summarizing
      )
    
    # Bootstrap-level statistics
    mean_r2_base <- mean(processed$r2_value, na.rm = TRUE)
    mean_cor_base <- mean(processed$correlation, na.rm = TRUE)
    sd_r2_base <- sd(processed$r2_value, na.rm = TRUE)
    n_base <- nrow(processed) # Number of observations in this bootstrap
    
    # Return a summary as a list
    list(
      mean_r2_base = mean_r2,
      mean_cor_base = mean_cor,
      sd_r2_base = sd_r2,
      total_n_base = n,
      bootstrap_results_base = processed # Include detailed bootstrap-level data
    )
  })
  
  # Combine results for overall summary
  overall_summary <- do.call(rbind, lapply(results_summary, function(x) {
    data.frame(mean_r2 = x$mean_r2,mean_cor = x$mean_cor,  sd = x$sd_r2, total_n = x$total_n)
  }))
  
  overall_summary_base <- do.call(rbind, lapply(results_summary_base, function(x) {
    data.frame(mean_r2_base = x$mean_r2_base,mean_cor_base = x$mean_cor_base,  sd_base = x$sd_r2_base, total_n_base = x$total_n_base)
  }))
  
  results <- list(
    rsq = rsq,
    cor = cor,
    mape = mape,
    mae = mae,
    rae = rae,
  #   ind = ind,
    rsq_base = rsq_base,
    cor_base = cor_base,
    mape_base = mape_base,
    mae_base = mae_base,
    rae_base = rae_base,
   overall_summary = overall_summary,
   overall_summary_base = overall_summary_base
  )
  
  results <- as.data.frame(results)
  
  
  return(results)
  }      

