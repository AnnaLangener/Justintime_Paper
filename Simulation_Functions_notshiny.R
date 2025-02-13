##############################
### Simulation Functions #####
##############################

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
} #taken from https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance/12239#12239

create_data <- function(n_features,n_samples,n_subjects,A,feature_std,B,C,overall_prob_outcome,sd_outcome,time_effect){
  
  y <- vector("numeric", length = n_subjects * n_samples) 
  subject_id <- rep(1:n_subjects, each = n_samples)  # Repeat subject IDs for each timepoint
  time_id <- rep(1:n_samples, times = n_subjects)  # Repeat timepoints for each subject
  data <- data.frame(subject = subject_id, time = time_id, y = y)
  
  # Simulate Outcome
  para <- estBetaParams(overall_prob_outcome,sd_outcome^2)
  prob_l <- rbeta(n=n_subjects ,shape1 = para$alpha, shape2 = para$beta) # sample probability for each subject
  
  print(paste("Mean prob:",mean(prob_l)))
  print(paste("SD prob:",sd(prob_l)))
  
  prob_l = rep(prob_l, each = n_samples) 
  
  if (time_effect) {
    #Can be added in the future
  } else {
    prob_matrix = prob_l
  }
  
  data$y <- rbinom(n_subjects * n_samples, size = 1, prob = prob_matrix) #sample 0 or 1 for each subject & timepoint, prob_matrix indicates the probability for each subject
  data$A = data$y

  # Simulate Features (adapted from Saeb et al. and https://aosmith.rbind.io/2018/04/23/simulate-simulate-part-2/)
  data$A[data$y == 0] = -A # relationship to outcome
  data$A[data$y == 1] = A
  
  features_sample <- list()
  
  features_sample_ABCstd <- as.data.frame(array(0, dim = c(n_subjects * n_samples, n_features)))
  features_sample_ABCstd <- cbind(features_sample_ABCstd,data)
  
  features_sample_Astd <- as.data.frame(array(0, dim = c(n_subjects * n_samples, n_features)))
  features_sample_Astd <- cbind(features_sample_Astd,data)
  
  features_sample[[1]] <- features_sample_ABCstd
  features_sample[[2]] <- features_sample_Astd
  
  for (i in 1:n_features) {
    
    if(length(B) == 1){
    feature_noise <- rep(feature_std * rnorm(1,mean = 0, sd = 1), each = n_subjects*n_samples) #d_e POPULATION LEVEL FEATURE GENERATING PROCESS
    subjecteff = B *rnorm(n_subjects, mean = 0, sd = 1) #B * rnorm(n_subjects) 
    subjecteff = rep(subjecteff, each = n_samples)
    timeeff = C * rnorm(n_subjects*n_samples, mean = 0, sd = 1) #C * rnorm(n_subjects*n_samples)
    maineff = features_sample[[1]]$A
    
    features_sample[[1]][,i] = maineff + feature_noise + subjecteff + timeeff 
    features_sample[[1]][, paste("V", i, "_mean", sep = "")] <- subjecteff
    features_sample[[2]][,i] = maineff + feature_noise + timeeff 
    
    }else{
      feature_noise <- rep(feature_std * rnorm(1,mean = 0, sd = 1), each = n_subjects*n_samples) #d_e POPULATION LEVEL FEATURE GENERATING PROCESS
      subjecteff = B[i] *rnorm(n_subjects, mean = 0, sd = 1) #B * rnorm(n_subjects) 
      subjecteff = rep(subjecteff, each = n_samples)
      timeeff = C[i] * rnorm(n_subjects*n_samples, mean = 0, sd = 1) #C * rnorm(n_subjects*n_samples)
      maineff = features_sample[[1]]$A
      
      features_sample[[1]][,i] = maineff + feature_noise + subjecteff + timeeff 
      features_sample[[1]][, paste("V", i, "_mean", sep = "")] <- subjecteff
      features_sample[[2]][,i] = maineff + feature_noise + timeeff 
    }
    
  }
  return(features_sample)
}


########################################################
############## Cross Validation Strategies #############
########################################################

############# Subject Wise and Row Wise ################

run_simulation <- function(features_sample,cv,n_bootstrap,testsize, seed = "12361488"){
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
  
  print(cv)
  
  for (i in 1:n_bootstrap) {
    
  #### Row wise and subject wise ####
  #  Prepare training and testing sets
    if(cv == "row-wise"){
      samples_test <- base::sample(x = n_subjects * n_samples,size = testsize * n_subjects * n_samples) 
      samples_train <- setdiff(1:(n_subjects * n_samples), samples_test) 
      
      train_X = features_sample[samples_train, 1:n_features]
      test_X = features_sample[samples_test, 1:n_features]
      
      train_Y =  as.factor(features_sample$y[samples_train])
      test_Y = as.factor(features_sample$y[samples_test])
      
      subject[[i]] = features_sample$subject[samples_test]
      
    }else if(cv == "subject-wise"){
      subjects_test <- sample(n_subjects -round(testsize*n_subjects) , round(testsize*n_subjects)) 
      subjects_train <- setdiff(1:n_subjects, subjects_test)  
      
      train_X = features_sample[features_sample$subject %in% subjects_train, 1:n_features]
      test_X = features_sample[features_sample$subject %in% subjects_test,1:n_features]
      
      train_Y = as.factor(features_sample$y[features_sample$subject %in% subjects_train])
      test_Y = as.factor(features_sample$y[features_sample$subject %in% subjects_test])
      
      subject[[i]] = features_sample$subject[features_sample$subject %in% subjects_test]
    }  
    
    # Train and evaluate the Random Forest model  
    if(cv == "subject-wise" | cv == "row-wise"){
      RF <- randomForest(train_X,train_Y)
      class_pred <- predict(RF, test_X)
      acc[i] <- mean(as.numeric(as.character(class_pred)) == test_Y)
      
      true_list[[i]]<- test_Y
      pred_list[[i]] <- as.numeric(as.character(class_pred))
      roc_curve <- roc(test_Y,  as.numeric(as.character(class_pred)),quiet = TRUE)
      auc_value[i] <- auc(roc_curve)
      
      ind[[i]] = data.frame(subject = subject[[i]], true = true_list[[i]], pred = pred_list[[i]])
      
      # Baseline model for row-wise CV
      if(cv == "row-wise"){ 
        Baseline <- glmer(
          train_Y ~ 1 + (1 | subject),
          data = features_sample[samples_train,], 
          family = binomial(link = "logit") 
        )
      
      class_pred_base <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject), type = "response")
      
      true_list_base[[i]]<- test_Y
      pred_list_base[[i]] <- as.numeric(as.character(class_pred_base))
      roc_curve_base <- roc(test_Y,  as.numeric(as.character(class_pred_base)),quiet = TRUE)
      auc_value_base[i] <- auc(roc_curve_base)
      
      class_pred_base2 <- ifelse(class_pred_base > 0.5, 1, 0)
      acc_base[i] <- mean(as.numeric(as.character(class_pred_base2)) == test_Y)
      
      ind_base[[i]] = data.frame(subject = subject[[i]], true = true_list_base[[i]], pred = pred_list_base[[i]])
      }
    }
  } # end bootstrap

  if(cv == "row-wise"){ 
    print(paste("Baseline Mean AUC:",mean(auc_value_base)))
    print(paste("Baseline Mean Accuracy:",mean(acc_base)))
    print("Model Results:")
    
  }
  
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
    acc = acc,
    auc_value = auc_value,
    ind = ind,
    auc_value_base = auc_value_base,
    overall_summary = overall_summary,
    ind_base = ind_base
  )
  
  print(paste("Mean AUC:",mean(auc_value)))
  print(paste("Mean Accuracy:",mean(acc)))
  print("Individual Results Summary:")
  print(overall_summary)
  return(results)
} 
############################################
run_simulation_centering <- function(features_sample,cv,n_bootstrap,testsize, seed = "12361488"){
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
  
  print(cv)
  
  for (i in 1:n_bootstrap) {
    
    #### Row wise and subject wise ####
    #  Prepare training and testing sets
    
    if(cv == "row-wise"){
    samples_test <- base::sample(x = n_subjects * n_samples,size = testsize * n_subjects * n_samples) 
    samples_train <- setdiff(1:(n_subjects * n_samples), samples_test) 
    
    id_variable = "subject"
    
    subject_means <- data.frame(subject = unique(features_sample[[id_variable]]))
    
    features_sample <- features_sample %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(features_sample[,1:10])), ~ . - get(paste0(cur_column(), "_mean")))) 
    
    train_X = features_sample[samples_train, 1:n_features]
    test_X = features_sample[samples_test, 1:n_features]
    
    train_Y =  as.factor(features_sample$y[samples_train])
    test_Y = as.factor(features_sample$y[samples_test])
    
    subject[[i]] = features_sample$subject[samples_test]
    
    feature_names = colnames(features_sample)[1:10]
   
  
    # Center the training set
    subject_means <- features_sample[samples_train, ] %>%
     group_by(subject) %>%
     summarise(across(all_of(names(train_X)), mean), .groups = "drop")

    # Center the training set
    train_X <- features_sample[samples_train, ] %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(train_X)), ~ . - get(paste0(cur_column(), "_mean")))) %>%
      select(all_of(names(train_X)))

    # Center the test set using training set means
    test_X <- features_sample[samples_test, ] %>%
      left_join(subject_means, by = "subject", suffix = c("", "_mean")) %>%
      mutate(across(all_of(names(train_X)), ~ . - get(paste0(cur_column(), "_mean")))) %>%
      select(all_of(names(train_X)))
    }
    
    # Train and evaluate the Random Forest model  
    if(cv == "row-wise"){
      RF <- randomForest(train_X,train_Y)
      class_pred <- predict(RF, test_X)
      acc[i] <- mean(as.numeric(as.character(class_pred)) == test_Y)
      
      true_list[[i]]<- test_Y
      pred_list[[i]] <- as.numeric(as.character(class_pred))
      roc_curve <- roc(test_Y,  as.numeric(as.character(class_pred)),quiet = TRUE)
      auc_value[i] <- auc(roc_curve)
      
      ind[[i]] = data.frame(subject = subject[[i]], true = true_list[[i]], pred = pred_list[[i]])
      
      # Baseline model for row-wise CV
      if(cv == "row-wise"){ 
        Baseline <- glmer(
          train_Y ~ 1 + (1 | subject),
          data = features_sample[samples_train,], 
          family = binomial(link = "logit") 
        )
        
        class_pred_base <- predict(Baseline, newdata = features_sample[samples_test,], re.form = ~(1 | subject), type = "response")
        
        true_list_base[[i]]<- test_Y
        pred_list_base[[i]] <- as.numeric(as.character(class_pred_base))
        roc_curve_base <- roc(test_Y,  as.numeric(as.character(class_pred_base)),quiet = TRUE)
        auc_value_base[i] <- auc(roc_curve_base)
        
        class_pred_base2 <- ifelse(class_pred_base > 0.5, 1, 0)
        acc_base[i] <- mean(as.numeric(as.character(class_pred_base2)) == test_Y)
        
        ind_base[[i]] = data.frame(subject = subject[[i]], true = true_list_base[[i]], pred = pred_list_base[[i]])
      }
    }
  } # end bootstrap
  
  if(cv == "row-wise"){ 
    print(paste("Baseline Mean AUC:",mean(auc_value_base)))
    print(paste("Baseline Mean Accuracy:",mean(acc_base)))
    print("Model Results:")
    
  }
  
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
    acc = acc,
    auc_value = auc_value,
    ind = ind,
    overall_summary = overall_summary,
    ind_base = ind_base
  )
  
  print(paste("Mean AUC:",mean(auc_value)))
  print(paste("Mean Accuracy:",mean(acc)))
  print("Individual Results Summary:")
  print(overall_summary)
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

run_simulation_slidingwindow <- function(features_sample,n_bootstrap,windowsize){
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
  
  for (i in 1:n_bootstrap) {
    #  Prepare training and testing sets
    timeSlices <- SlidingWindow_CV(n_samples,windowsize,1) #window size, prediction horizon
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
      train_X = features_sample[features_sample$time %in% trainSlices[[k]], 1:n_features]
      test_X = features_sample[features_sample$time %in% testSlices[[k]], 1:n_features]
      train_Y = as.factor(features_sample$y[features_sample$time %in% trainSlices[[k]]])
      test_Y = as.factor(features_sample$y[features_sample$time %in% testSlices[[k]]])
      
      # Train and evaluate the Random Forest model 
      
      RF <- randomForest(train_X,train_Y)
      class_pred <- predict(RF, test_X)
      acc_sw[k] <- mean(as.numeric(as.character(class_pred)) == test_Y)
      
      if (length(unique(test_Y)) == 2) {
        roc_curve <- roc(test_Y, as.numeric(as.character(class_pred)), quiet = TRUE)
        auc_value_sw[k] <- auc(roc_curve)
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
      
      #### Baseline
      Baseline <- glmer(
        "y ~ 1 + (1|subject)", # Random intercept for 'subject'
        data = features_sample[features_sample$time %in% trainSlices[[k]],], 
        family = "binomial" # Logistic regression
      )
      
      class_pred_base <- predict(Baseline, newdata = features_sample[features_sample$time %in% testSlices[[k]],], re.form = ~(1 | subject), type = "response")

      if (length(unique(test_Y)) == 2) {
        roc_curve_base <- roc(test_Y,  as.numeric(as.character(class_pred_base)),quiet = TRUE)
        auc_value_base_sw[k] <- auc(roc_curve_base)
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
      
      
      print(paste("Window ",k,": ","Model: ", round(auc_value_sw[k],2), " Baseline: ",round(auc_value_base_sw[k],2), sep = ""))
      
    }
    
    auc_value[i] <- auc(roc(as.numeric(as.character(ind_t$true)), as.numeric(as.character(ind_t$pred)), quiet = TRUE))
    auc_value_meansw[i] <- mean(auc_value_sw,na.rm = TRUE)
    ind[[i]] <- ind_t
    ind_base[[i]] <- ind_t_base
    auc_value_base[i] <- auc(roc(as.numeric(as.character(ind_t_base$true)), as.numeric(as.character(ind_t_base$pred)), quiet = TRUE))
  } # end bootstrap
 
  
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
    overall_summary = overall_summary,
    ind_base = ind_base,
    auc_value_base = auc_value_base
  )
  
  print(paste("Baseline Mean AUC:",mean(auc_value_base)))
  print("Model Results:")
  print(paste("Mean AUC:",mean(auc_value)))
  print("Individual Results Summary:")
  print(overall_summary)
  return(results)
}      

