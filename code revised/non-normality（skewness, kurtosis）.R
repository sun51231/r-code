




library(sjPlot)
library(MASS)
library(mvtnorm)
library(lavaan) 
library(CompQuadForm) 
library(BAMMtools) 
library(Ckmeans.1d.dp) 
library(expm)
library(Matrix)
library(gtools)  
library(regsem)
library(semTests)    
library(future.apply)


# ======================
# 1. Constructing the CFA model covariance matrix function    
# ======================

create_cfa_cov <- function(n_factors, lambda = 0.7, factor_cor = 0.3) {
  n_indicators <- 5  # 5 indicators for each factor
  total_vars <- n_factors * n_indicators
  
  # Create factor covariance matrix
  Phi <- matrix(factor_cor, nrow = n_factors, ncol = n_factors)
  diag(Phi) <- 1
  
  # Create factor loading matrix
  Lambda <- matrix(0, nrow = total_vars, ncol = n_factors)
  
  for (i in 1:n_factors) {
    start_idx <- (i-1)*n_indicators + 1
    end_idx <- i*n_indicators
    Lambda[start_idx:end_idx, i] <- lambda
  }
  
  # Residual Matrix(1 - λ²)
  Theta <- diag(1 - lambda^2, total_vars)
  
  # Observed variable covariance matrix: Σ = ΛΦΛ' + Θ
  Sigma <- Lambda %*% Phi %*% t(Lambda) + Theta
  
  # Add variable name
  var_names <- paste0("x", 1:total_vars)
  colnames(Sigma) <- rownames(Sigma) <- var_names
  
  return(Sigma)
}

# ======================
# 2. Generates moderately nonnormal data
# ======================

###    normal (skewness=0, kurtosis=3)
###   Moderate non-normality (skewness=1, kurtosis=3)
###   Substantial non-normality (skewness=2, kurtosis=7)
####  Severe non-normality (skewness=3, kurtosis=21)
generate_cfa_data <- function(n_factors, n, skew = 3, kurt = 21) {
  # Create a covariance matrix.
  Sigma <- create_cfa_cov(n_factors)
  total_vars <- n_factors * 5
  
  # Generate data using the Vale-Maurelli method
  data_nonnorm <- semTools::mvrnonnorm(
    n = n,
    mu = rep(0, total_vars),
    Sigma = Sigma,
    skewness = rep(skew, total_vars),
    kurtosis = rep(kurt, total_vars)
  )
  
  colnames(data_nonnorm) <- paste0("x", 1:total_vars)
  return(as.data.frame(data_nonnorm))
}

# ======================
# 3. Defining CFA model syntax
# ======================

create_cfa_syntax <- function(n_factors) {
  syntax <- ""
  
  # Defining latent variables
  for (i in 1:n_factors) {
    start_idx <- (i-1)*5 + 1
    indicators <- paste0("x", start_idx:(start_idx+4), collapse = " + ")
    syntax <- paste0(syntax, sprintf("f%d =~ %s\n", i, indicators))
  }
  
  # Adding correlations between factors
  if (n_factors > 1) {
    factors <- paste0("f", 1:n_factors)
    corrs <- combn(factors, 2, function(pair) paste(pair, collapse = " ~~ "))
    syntax <- paste0(syntax, paste(corrs, collapse = "\n"))
  }
  
  return(syntax)
}


### Changing the number of factors
m<-n_factors<-12





############  #  Revised 


# Convergence count statistics
fit_non_converged_count <- 0  


# Check if the matrix is positive definite
is_positive_definite <- function(matrix) {
  if (is.null(matrix) || any(is.na(matrix))) {
    return(FALSE)
  }
  
  # Check if it is a numerical matrix
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    return(FALSE)
  }
  
  # Check if it is a square matrix
  if (nrow(matrix) != ncol(matrix)) {
    return(FALSE)
  }
  
  # Check symmetry (numerical tolerance)
  symmetric_tolerance <- 1e-8
  is_symmetric <- all(abs(matrix - t(matrix)) < symmetric_tolerance)
  
  if (!is_symmetric) {
    # Attempt symmetrization
    matrix <- (matrix + t(matrix)) / 2
  }
  
  # Compute eigenvalues (values only)
  eigenvalues <- tryCatch({
    eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  }, error = function(e) {
    return(NA)
  })
  
  if (any(is.na(eigenvalues))) {
    return(FALSE)
  }
  
  # Positive definite condition: All eigenvalues > 0
  tolerance <- 0             
  return(all(eigenvalues > tolerance))
}

# Function to check model convergence (using more practical methods)
check_model_convergence <- function(fit) {
  # Check 1: Does the optimizer report convergence?
  if (!fit@optim$converged) {
    return(FALSE)
  }
  
  # Check 2: Attempt to check the positive definiteness of the information matrix
  tryCatch({
    # First attempt to obtain the Hessian matrix
    hessian_matrix <- lavInspect(fit, "hessian")
    
    if (!is.null(hessian_matrix)) {
      if (is_positive_definite(hessian_matrix)) {
        return(TRUE)
      } 
    }
    
    # If Hessian is unavailable, attempt the information matrix
    info_matrix <- lavInspect(fit, "information")
    if (!is.null(info_matrix)) {
      if (is_positive_definite(info_matrix)) {
        return(TRUE)
      }
    }
    
    # If both fail, check the variance-covariance matrix
    vcov_matrix <- vcov(fit)
    if (!is.null(vcov_matrix)) {
      if (is_positive_definite(vcov_matrix)) {
        return(TRUE)
      }
    }
    
  }, error = function(e) {
    # Return FALSE if any check fails
    return(FALSE)
  })
  
  # Check 3: Ensure valid parameter estimates exist
  param_estimates <- tryCatch({
    coef(fit)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(param_estimates) || any(is.na(param_estimates)) || any(is.infinite(param_estimates))) {
    return(FALSE)
  }
  
  # Check 4: Ensure standard errors are reasonable
  standard_errors <- tryCatch({
    fit@Fit@se
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(standard_errors) || any(is.na(standard_errors))) {
    return(FALSE)
  }
  
  # For standard errors equal to zero, check for a reasonable explanation (e.g., fixed parameters)
  if (any(standard_errors <= 0, na.rm = TRUE)) {
    # If standard errors are zero or negative, check if they are fixed parameters
    free_params <- fit@ParTable$free > 0
    if (any(standard_errors[free_params] <= 0, na.rm = TRUE)) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# Modified get_pvalues function
get_pvalues_optimized <- function(n) {
  tryCatch({
    # 
    model <- create_cfa_syntax(n_factors)
    dat<- generate_cfa_data(n_factors = n_factors, n = n)
    
    # 
    yn <- ncol(dat)/nrow(dat)
    tau <- ncol(dat) - ncol(dat)*(1-1/yn)*log(1-yn) - 0.5*log(1-yn)
    v <- sqrt(-2*log(1-yn) - 2*yn)
    
    mu_y1 <- 2.015910*(m/ncol(dat)) + 1.291412*(ncol(dat)/nrow(dat)) - 0.278377*m + 0.036066*ncol(dat) - 2.393643
    mu_y2 <- -49.62*(ncol(dat)/(nrow(dat))^2) + 4.270*(ncol(dat)/nrow(dat)) + 0.003475*nrow(dat) - 0.1819*m - 2.626
    mu_y3 <- 1.531018*(ncol(dat)/nrow(dat)) - 0.237301*m + 0.029336*ncol(dat) - 2.035224
    mu_y4 <- -50.735178*(ncol(dat)/(nrow(dat))^2) + 4.788150*(ncol(dat)/nrow(dat)) + 0.003553*nrow(dat) - 0.183627*m - 2.639546
    
    # Initialize result vector (14 NA values)
    result <- rep(NA, 14)
    names(result) <- c("p_mlb", "p_rmlb", "p_rls", "p_CTrls", 
                       "pT_F", "pT_F_c1", "pT_F_c2", "pT_F_cr",
                       "pT_CsF", "pT_CsF_c1", "pT_CsF_c2", "pT_CsF_cr",
                       "pEBA4_rls", "pols_2_rls")
    
    # ==================== Type I: Standard statistics ====================
    fit <- NULL
    fit_converged <- FALSE
    
    tryCatch({
      # The default optimization algorithm, “NLMINB”, is robust when handling problems with boundary constraints
      fit <- cfa(model, data = dat, estimator = "MLM",test = c("standard", "browne.residual.nt.model"))
      
      # Check for convergence
      fit_converged <- check_model_convergence(fit)
      
      if (!fit_converged) {
        # Record if fit fails to converge
        assign("fit_non_converged_count", fit_non_converged_count + 1, envir = .GlobalEnv)
      }
      
    }, error = function(e) {
      # Error during fitting
      assign("fit_non_converged_count", fit_non_converged_count + 1, envir = .GlobalEnv)
    })
    
    # If fit converges, compute Type I p-value
    if (fit_converged && !is.null(fit)) {
      # Compute the original statistic
      cb <- (1 - (2*ncol(dat) + 4*m + 11)/(6*n)) * n / (n - 1)
      chisq_mlb <- cb * fitmeasures(fit, "chisq")
      chisq_rmlb <- cb * fitmeasures(fit, "chisq.scaled")
      df <- fitmeasures(fit, "df")
      
      result["p_mlb"] <- 1 - pchisq(chisq_mlb, df)
      result["p_rmlb"] <- 1 - pchisq(chisq_rmlb, df)
      
      # Check if browne.residual.nt.model exists
      if (!is.null(fit@test$browne.residual.nt.model)) {
        result["p_rls"] <- fit@test$browne.residual.nt.model$pvalue
        
        T_rls <- fit@test$browne.residual.nt.model$stat
        CT_rls <- T_rls / fitMeasures(fit, "chisq.scaling.factor")
        result["p_CTrls"] <- 1 - pchisq(CT_rls, df)
      }
      
      # Newly proposed statistic
      chisq_val <- fitmeasures(fit, "chisq")
      chisq_scaled_val <- fitmeasures(fit, "chisq.scaled")
      
      T_F <- (chisq_val/nrow(dat) - tau)/v
      result["pT_F"] <- 1 - pnorm(T_F, 0, 1)
      
      T_F_c1 <- (chisq_val/nrow(dat) - tau - mu_y1*v)/v
      result["pT_F_c1"] <- 1 - pnorm(T_F_c1, 0, 1)
      
      T_F_c2 <- (chisq_val/nrow(dat) - tau - mu_y2*v)/v
      result["pT_F_c2"] <- 1 - pnorm(T_F_c2, 0, 1)
      
      T_CsF <- (chisq_scaled_val/nrow(dat) - tau)/v
      result["pT_CsF"] <- 1 - pnorm(T_CsF, 0, 1)
      
      T_CsF_c1 <- (chisq_scaled_val/nrow(dat) - tau - mu_y3*v)/v
      result["pT_CsF_c1"] <- 1 - pnorm(T_CsF_c1, 0, 1)
      
      T_CsF_c2 <- (chisq_scaled_val/nrow(dat) - tau - mu_y4*v)/v
      result["pT_CsF_c2"] <- 1 - pnorm(T_CsF_c2, 0, 1)
    }
    
    # ==================== Type II: Ridge statistics. ====================
    # 
    fit1 <- NULL
    fit1_converged <- FALSE
    
    tryCatch({
      # 
      fit1 <- cfa(model, data = dat, ridge = TRUE, ridge.constant = 1e-3,estimator = "MLM")
      
      # 
      fit1_converged <- check_model_convergence(fit1)
      
    }, error = function(e) {
      # Error during fit1 fitting
    })
    
    # Conditions for calculating Type II p-values:
    # 1.fit fails to converge but fit1 converges
    # 2.fit converges and fit1 also converges
    if ((!fit_converged && fit1_converged) || (fit_converged && fit1_converged)) {
      if (!is.null(fit1) && fit1_converged) {
        T_F_cr <- (fitmeasures(fit1, "chisq")/nrow(dat) - tau - mu_y2*v)/v
        result["pT_F_cr"] <- 1 - pnorm(T_F_cr, 0, 1)
        
        T_CsF_cr <- (fitmeasures(fit1, "chisq.scaled")/nrow(dat) - tau - mu_y4*v)/v
        result["pT_CsF_cr"] <- 1 - pnorm(T_CsF_cr, 0, 1)
      }
    }
    
    # ==================== Type III: EBA statistic.====================
    # 
    if (!is.null(fit)) {
      tryCatch({
        pv_EBA <- semTests::pvalues(fit)
        # Ensure pv_EBA has sufficient length.
        if (length(pv_EBA) >= 6) {
          result["pEBA4_rls"] <- pv_EBA[5]
          result["pols_2_rls"] <- pv_EBA[6]
        }
      }, error = function(e) {
        # If calculation fails, retain NA
      })
    }
    
    return(result)
    
  }, error = function(e) {
    # If error occurs, return all NA
    result <- rep(NA, 14)
    names(result) <- c("p_mlb", "p_rmlb", "p_rls", "p_CTrls", 
                       "pT_F", "pT_F_c1", "pT_F_c2", "pT_F_cr",
                       "pT_CsF", "pT_CsF_c1", "pT_CsF_c2", "pT_CsF_cr",
                       "pEBA4_rls", "pols_2_rls")
    return(result)
  })
}




plan(multisession)
set.seed(123)
n_simulations <- 1000
# 
fit_non_converged_count <- 0
# 
res <- future_replicate(n_simulations, get_pvalues_optimized(200), future.seed = TRUE)

# Transform results into a matrix (14 rows x n_simulations columns).
result_matrix <- matrix(unlist(res), nrow = 14, byrow = FALSE)
rownames(result_matrix) <- c("p_mlb", "p_rmlb", "p_rls", "p_CTrls", 
                             "pT_F", "pT_F_c1", "pT_F_c2", "pT_F_cr",
                             "pT_CsF", "pT_CsF_c1", "pT_CsF_c2", "pT_CsF_cr",
                             "pEBA4_rls", "pols_2_rls")

# 
cat("==================================================\n")
cat(sprintf("Total simulation count:                   %d\n", n_simulations))
cat(sprintf("Number of fit failures:                %d (%.1f%%)\n", 
            fit_non_converged_count, 
            fit_non_converged_count/n_simulations*100))
cat("==================================================\n\n")

# 
cat("\n============== Table output ===============\n")
L <- data.frame(matrix(nrow = 1, ncol = 14))
colnames(L) <- rownames(result_matrix)
rejection_rates <- rowMeans(result_matrix < 0.05, na.rm = TRUE)
L[1, ] <- rejection_rates
print(tab_df(L, digits = 3))








