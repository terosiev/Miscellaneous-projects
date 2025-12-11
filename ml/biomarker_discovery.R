# Script: Nested cross-validation for logistic and Cox regresion with feature selection
# Purpose: To test if blood-based features predict treatment response in cancer patients
# Date: 11.12.2025
# Note: Claude.ai was used to clean the original script

# Load libraries
library(glmnet)
library(survival)
library(caret)
library(pROC)


## FEATURE SELECTION ----

#' Select top features based on selection frequency across CV folds
#' 
#' @param feature_freq Named vector of feature selection frequencies
#' @param n Maximum number of features to return
#' @param exclude Features to exclude from selection
#' @param min_prop Minimum proportion of folds a feature must appear in
#' @return Character vector of top features
top_n_features <- function(feature_freq, n = NULL, exclude = NULL, min_prop = 0.5) {
  # Convert list to numeric vector (if needed)
  freq_table <- unlist(feature_freq)
  
  # Remove any features specified in exclude list
  freq_table <- freq_table[!(names(freq_table) %in% exclude)]
  
  # Return empty if no features remain
  if (length(freq_table) == 0) return(character(0))
  
  # Calculate stability threshold: feature must appear in at least min_prop of folds
  n_folds <- max(freq_table)
  keep <- names(freq_table[freq_table >= min_prop * n_folds])
  
  # Return top n features (or all if n is NULL)
  head(keep, n)
}


## NESTED CV FOR BINARY OUTCOMES ----

#' Nested cross-validation for logistic regression with feature selection
#' 
#' Performs nested CV using penalized logistic regression (lasso/elastic net)
#' for feature selection in the inner loop, then fits unpenalized models
#' for unbiased performance estimation in the outer loop.
#' 
#' @param X Matrix of predictors (samples x features)
#' @param y Binary outcome vector
#' @param covariate_df Data frame of covariates to force into models
#' @param n_features Maximum number of features to select
#' @param seed Random seed
#' @param alpha Elastic net mixing (1 = lasso, 0 = ridge)
#' @param repeats Number of inner CV repeats for stability
#' @return List with results (fold-level AUCs), feature_freq, all_selected
nested_cv_log <- function(X, 
                          y, 
                          covariate_df, 
                          n_features = 10, 
                          seed = 123, 
                          alpha = 1,
                          repeats = 10) {   
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Create outer folds for testing (3-fold CV)
  # Each fold will be held out once as a test set
  outer_folds <- createFolds(y, k = 3, returnTrain = FALSE)
  
  # Initialize storage for results and selected features
  results <- list()
  feature_list <- list()
  
  # OUTER LOOP: iterate over each test fold
  for (i in seq_along(outer_folds)) {
    # Initialize empty feature list for this fold
    feature_list[[i]] <- character(0)
    
    # Split data into train and test sets
    test_idx  <- outer_folds[[i]]
    train_idx <- setdiff(seq_along(y), test_idx)
    
    # Skip fold if training set has only one class (can't train model)
    if (length(unique(y[train_idx])) < 2) next
    
    # Split predictors into train/test
    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Split covariates into train/test
    cov_train <- covariate_df[train_idx, , drop = FALSE]
    cov_test  <- covariate_df[test_idx, , drop = FALSE]
    
    # Drop unused factor levels to prevent errors
    cov_train <- droplevels(cov_train)
    cov_test  <- droplevels(cov_test)
    
    # Build design matrices (convert categorical variables to dummy codes)
    # Combine predictors and covariates
    train_all <- cbind(X_train, cov_train)
    test_all  <- cbind(X_test, cov_test)
    
    # Create model matrix (handles factor encoding automatically)
    terms_obj <- terms(~ ., data = train_all)
    mm_train <- model.matrix(terms_obj, data = train_all)[, -1, drop = FALSE]  # Remove intercept
    mm_test  <- model.matrix(terms_obj, data = test_all)[, -1, drop = FALSE]
    
    # Identify which columns are covariates (should not be penalized)
    mm_cov   <- model.matrix(~ ., data = cov_train)[, -1, drop = FALSE]
    cov_cols <- intersect(colnames(mm_train), colnames(mm_cov))
    
    # Create penalty factor: 0 for covariates (forced in), 1 for predictors (penalized)
    penalty <- ifelse(colnames(mm_train) %in% cov_cols, 0, 1)
    
    # INNER LOOP: Repeat CV multiple times for stable lambda selection
    # Each repeat uses different random folds
    cvfit_list <- lapply(1:repeats, function(rep) {
      # Use different seed for each repeat to get different fold splits
      set.seed(seed + rep)
      
      # Fit penalized logistic regression with cross-validation
      tryCatch(
        cv.glmnet(
          x = mm_train,           # Training predictors
          y = y_train,            # Training outcome
          family = "binomial",    # Logistic regression
          alpha = alpha,          # 1 = lasso, 0 = ridge
          nfolds = 3,            # 3-fold CV for lambda tuning
          penalty.factor = penalty,  # Don't penalize covariates
          standardize = TRUE,     # Standardize predictors
          maxit = 1e6            # Max iterations
        ),
        error = function(e) NULL  # Return NULL if fit fails
      )
    })
    
    # Remove any failed fits (NULL values)
    cvfit_list <- Filter(Negate(is.null), cvfit_list)
    
    # If all fits failed, skip this fold
    if (length(cvfit_list) == 0) {
      results[[i]] <- data.frame(fold = i, auc = NA_real_)
      feature_list[[i]] <- character(0)
      next
    }
    
    # Select the best model: lowest cross-validation error across all repeats
    best_fit <- cvfit_list[[which.min(sapply(cvfit_list, function(f) min(f$cvm)))]]
    
    # Extract coefficients at optimal lambda (lambda.min)
    coefs <- as.matrix(coef(best_fit, s = "lambda.min"))
    
    # Find non-zero coefficients (selected features)
    nz <- which(coefs != 0)
    nz_names <- rownames(coefs)[nz]
    
    # Keep only predictor features (exclude intercept and covariates)
    selected <- setdiff(nz_names, c("(Intercept)", cov_cols))
    
    # Limit to top n_features
    selected <- head(selected, n_features)
    feature_list[[i]] <- selected

    # FINAL MODEL: Refit unpenalized logistic regression on selected features
    # This gives unbiased estimates of performance
    use_cols <- unique(c(selected, cov_cols))  # Include selected features + covariates
    
    # Fit standard logistic regression (no penalty)
    fit <- glm(y_train ~ ., 
               data = as.data.frame(mm_train[, use_cols, drop = FALSE]), 
               family = binomial())
    
    # Predict probabilities on held-out test set
    probs <- predict(fit, 
                    newdata = as.data.frame(mm_test[, use_cols, drop = FALSE]), 
                    type = "response")
    
    # Calculate AUC (area under ROC curve) for this fold
    auc <- tryCatch(auc(y_test, probs), error = function(e) NA_real_)
    
    # Store fold results
    results[[i]] <- data.frame(fold = i, auc = as.numeric(auc))
  }
  
  # Combine results from all folds
  results_df <- dplyr::bind_rows(results)
  
  # Count how often each feature was selected across folds (stability measure)
  all_features <- unlist(feature_list)
  feature_freq <- sort(table(all_features), decreasing = TRUE)
  
  # Return performance metrics and feature selection info
  list(
    results = results_df,           # AUC for each fold
    feature_freq = feature_freq,    # How often each feature was selected
    all_selected = unique(all_features)  # All unique features selected
  )
}


## NESTED CV FOR SURVIVAL OUTCOMES ----

#' Nested cross-validation for Cox regression with feature selection
#' 
#' Performs nested CV using penalized Cox regression (lasso/elastic net)
#' for feature selection, with optional stratification variable.
#' 
#' @param X Matrix of predictors (samples x features)
#' @param time Survival time vector
#' @param event Event indicator (1 = event, 0 = censored)
#' @param covariate_df Data frame of covariates (can include strata variable)
#' @param strata_var Name of stratification variable in covariate_df (optional)
#' @param n_features Maximum number of features to select
#' @param seed Random seed
#' @param alpha Elastic net mixing parameter
#' @param repeats Number of inner CV repeats
#' @return List with results (fold-level C-indices), feature_freq, all_selected
nested_cv_cox <- function(X, 
                          time, 
                          event, 
                          covariate_df,
                          strata_var = NULL,
                          n_features = 10, 
                          seed = 123, 
                          alpha = 1,
                          repeats = 10) {   
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Create outer folds for testing (3-fold CV)
  outer_folds <- createFolds(event, k = 3, returnTrain = FALSE)
  
  # Initialize storage
  results <- list()
  feature_list <- list()
  
  # OUTER LOOP: iterate over each test fold
  for (i in seq_along(outer_folds)) {
    # Initialize empty feature list for this fold
    feature_list[[i]] <- character(0)
    
    # Split data into train and test indices
    test_idx  <- outer_folds[[i]]
    train_idx <- setdiff(seq_along(event), test_idx)
    
    # Split predictors into train/test
    X_train <- as.data.frame(X[train_idx, , drop = FALSE])
    X_test  <- as.data.frame(X[test_idx, , drop = FALSE])
    
    # Split covariates into train/test
    cov_train <- as.data.frame(covariate_df[train_idx, , drop = FALSE])
    cov_test  <- as.data.frame(covariate_df[test_idx, , drop = FALSE])
    
    # Split survival outcomes into train/test
    time_train  <- time[train_idx]
    event_train <- event[train_idx]
    time_test   <- time[test_idx]
    event_test  <- event[test_idx]
    
    # Drop unused factor levels to prevent errors
    cov_train <- droplevels(cov_train)
    cov_test  <- droplevels(cov_test)

    # Handle stratification variable if provided
    # Stratification allows baseline hazards to differ between groups
    if (!is.null(strata_var) && strata_var %in% colnames(cov_train)) {
      # Ensure strata variable is a factor with consistent levels
      cov_train[[strata_var]] <- factor(cov_train[[strata_var]])
      cov_test[[strata_var]]  <- factor(cov_test[[strata_var]], 
                                        levels = levels(cov_train[[strata_var]]))
      
      # Remove strata variable from penalized regression (will add back for final model)
      cov_train_nostrata <- cov_train[, setdiff(colnames(cov_train), strata_var), drop = FALSE]
      cov_test_nostrata  <- cov_test[,  setdiff(colnames(cov_test),  strata_var), drop = FALSE]
    } else {
      # No stratification: use all covariates
      cov_train_nostrata <- cov_train
      cov_test_nostrata <- cov_test
    }
    
    # Build design matrices for penalized regression
    # Combine predictors and covariates (excluding strata)
    data_all <- cbind(X_train, cov_train_nostrata)
    
    # Create model matrix (handles factor encoding)
    mm_all   <- model.matrix(~ ., data = data_all)[, -1, drop = FALSE]
    
    # Identify covariate columns (will not be penalized)
    mm_cov   <- model.matrix(~ ., data = cov_train_nostrata)[, -1, drop = FALSE]
    cov_cols <- intersect(colnames(mm_all), colnames(mm_cov))
    
    # Create penalty factor: 0 for covariates, 1 for predictors
    penalty <- ifelse(colnames(mm_all) %in% cov_cols, 0, 1)
    
    # Create survival object for glmnet
    y_surv <- Surv(time_train, event_train)
    
    # INNER LOOP: Repeated CV for stable lambda selection
    cvfit_list <- lapply(1:repeats, function(rep) {
      # Different seed for each repeat
      set.seed(seed + rep)
      
      # Fit penalized Cox regression with cross-validation
      tryCatch(
        cv.glmnet(
          x = mm_all,              # Predictors
          y = y_surv,              # Survival outcome
          family = "cox",          # Cox proportional hazards
          alpha = alpha,           # 1 = lasso, 0 = ridge
          nfolds = 3,             # Inner CV folds
          standardize = TRUE,      # Standardize predictors
          penalty.factor = penalty # Don't penalize covariates
        ),
        error = function(e) NULL   # Return NULL if fit fails
      )
    })
    
    # Remove failed fits
    cvfit_list <- Filter(Negate(is.null), cvfit_list)
    
    # If all fits failed, skip this fold
    if (length(cvfit_list) == 0) {
      results[[i]] <- data.frame(fold = i, cindex = NA_real_)
      next
    }
    
    # Select best model: lowest CV error across repeats
    best_fit <- cvfit_list[[which.min(sapply(cvfit_list, function(f) min(f$cvm)))]]
    
    # Extract non-zero coefficients (selected features)
    coefs <- as.matrix(coef(best_fit, s = "lambda.min"))
    nz_idx <- which(coefs != 0)
    
    # If no features selected, skip this fold
    if (length(nz_idx) == 0) {
      results[[i]] <- data.frame(fold = i, cindex = NA_real_)
      next
    }
    
    # Get names of selected features
    nz_names <- rownames(coefs)[nz_idx]
    
    # Keep only predictor features (exclude covariates)
    feat_names <- setdiff(nz_names, cov_cols)
    
    # Select top n features by absolute coefficient size
    selected <- head(feat_names[order(abs(coefs[feat_names, 1]), decreasing = TRUE)], 
                    n_features)
    feature_list[[i]] <- c(feature_list[[i]], selected)

    # FINAL MODEL: Refit unpenalized Cox model with selected features
    # Build complete design matrices for final model
    train_all <- cbind(X_train, cov_train_nostrata)
    test_all  <- cbind(X_test,  cov_test_nostrata)
    terms_obj <- terms(~ ., data = train_all)
    mm_train  <- model.matrix(terms_obj, data = train_all)[, -1, drop = FALSE]
    mm_test   <- model.matrix(terms_obj, data = test_all)[,  -1, drop = FALSE]
    
    # Include selected features + covariates
    use_cols <- unique(c(selected, cov_cols))
    mm_train_use <- mm_train[, use_cols, drop = FALSE]
    mm_test_use  <- mm_test[,  use_cols, drop = FALSE]
    
    # If stratification is used, add strata variable back
    if (!is.null(strata_var)) {
      # Prepare training data with strata
      df_train <- data.frame(time = time_train, 
                            event = event_train, 
                            mm_train_use,
                            strata = cov_train[[strata_var]])
      
      # Prepare test data with strata
      df_test  <- data.frame(time = time_test,  
                            event = event_test,  
                            mm_test_use,
                            strata = cov_test[[strata_var]])
      
      # Remove test observations with unseen strata levels
      keep <- !is.na(df_test$strata)
      if (!any(keep)) {
        results[[i]] <- data.frame(fold = i, cindex = NA_real_)
        next
      }
      df_test <- df_test[keep, , drop = FALSE]
      
      # Fit stratified Cox model
      # strata() allows different baseline hazards for each stratum
      fit <- coxph(Surv(time, event) ~ . - strata + strata(strata),
                  data = df_train, 
                  model = FALSE,  # Don't store model frame (save memory)
                  x = FALSE,      # Don't store design matrix
                  y = FALSE)      # Don't store response
    } else {
      # No stratification: standard Cox model
      df_train <- data.frame(time = time_train, event = event_train, mm_train_use)
      df_test  <- data.frame(time = time_test,  event = event_test,  mm_test_use)
      
      fit <- coxph(Surv(time, event) ~ ., 
                  data = df_train,
                  model = FALSE, 
                  x = FALSE, 
                  y = FALSE)
    }
    
    # Predict risk scores (linear predictor) on test set
    risk <- predict(fit, newdata = df_test, type = "lp")
    
    # Calculate concordance index (C-index) for this fold
    # C-index measures discriminative ability (similar to AUC)
    cidx <- concordance(Surv(df_test$time, df_test$event) ~ risk)$concordance
    
    # Store fold results
    results[[i]] <- data.frame(fold = i, cindex = cidx)
  }
  
  # Combine results from all folds
  results_df <- dplyr::bind_rows(results)
  
  # Count feature selection frequency (stability measure)
  all_features <- unlist(feature_list)
  feature_freq <- sort(table(all_features), decreasing = TRUE)
  
  # Return performance metrics and feature selection info
  list(
    results = results_df,           # C-index for each fold
    feature_freq = feature_freq,    # Feature selection frequency
    all_selected = unique(all_features)  # All unique selected features
  )
}


## EXAMPLE USAGE ----

# Example 1: Binary outcome (treatment response)
if (FALSE) {
  # Simulate data
  set.seed(123)
  n <- 200  # Number of samples
  p <- 50   # Number of features
  
  # Create predictor matrix (e.g., gene expression data)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("gene_", 1:p)
  
  # Create binary outcome with signal from first 3 genes
  # Response depends on gene_1, gene_2, and gene_3
  y <- rbinom(n, 1, prob = plogis(X[,1] + 0.5*X[,2] - 0.3*X[,3]))
  
  # Create covariate data frame (clinical variables)
  covariates <- data.frame(
    age = rnorm(n, 60, 10),  # Patient age
    gender = factor(sample(c("M", "F"), n, replace = TRUE)),  # Gender
    batch = factor(sample(paste0("batch", 1:3), n, replace = TRUE))  # Technical batch
  )
  
  # Run nested cross-validation
  cv_results <- nested_cv_log(
    X = X,                    # Predictor matrix
    y = y,                    # Binary outcome
    covariate_df = covariates,  # Clinical covariates
    n_features = 10,          # Maximum features to select
    seed = 123,               # Random seed
    alpha = 1,                # Lasso (alpha=1), Ridge (alpha=0)
    repeats = 10              # Number of inner CV repeats
  )
  
  # View results
  print(cv_results$results)  # AUC for each fold
  print(cv_results$feature_freq)  # How many times each feature was selected
  
  # Select most stable features (appeared in >50% of folds)
  top_features <- top_n_features(cv_results$feature_freq, n = 5, min_prop = 0.5)
  print(top_features)
}


# Example 2: Survival outcome (time to event)
if (FALSE) {
  # Simulate survival data
  set.seed(456)
  n <- 200  # Number of patients
  p <- 50   # Number of features
  
  # Create predictor matrix (e.g., protein expression)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("protein_", 1:p)
  
  # Create survival outcome with signal from first 2 proteins
  # Higher protein_1 → shorter survival
  # Higher protein_2 → longer survival
  time <- rexp(n, rate = exp(0.5*X[,1] - 0.3*X[,2]))
  event <- rbinom(n, 1, 0.7)  # 70% event rate (30% censored)
  
  # Create covariates including a stratification variable
  covariates <- data.frame(
    age = rnorm(n, 65, 8),  # Patient age
    stage = factor(sample(c("I", "II", "III"), n, replace = TRUE)),  # Disease stage
    disease_type = factor(sample(c("TypeA", "TypeB"), n, replace = TRUE))  # Disease subtype
  )
  
  # Run nested CV with stratification by disease type
  # Stratification allows different baseline hazards for TypeA vs TypeB
  cv_results <- nested_cv_cox(
    X = X,                      # Predictor matrix
    time = time,                # Survival time
    event = event,              # Event indicator (1=event, 0=censored)
    covariate_df = covariates,  # Clinical covariates
    strata_var = "disease_type",  # Stratify by disease type
    n_features = 10,            # Maximum features to select
    seed = 456,                 # Random seed
    alpha = 1,                  # Lasso penalty
    repeats = 10                # Inner CV repeats
  )
  
  # View results
  print(cv_results$results)  # C-index for each fold (higher is better)
  print(cv_results$feature_freq)  # Feature selection stability
  
  # Select stable biomarkers (appeared in >60% of folds)
  biomarkers <- top_n_features(cv_results$feature_freq, n = 5, min_prop = 0.6)
  print(biomarkers)
}
