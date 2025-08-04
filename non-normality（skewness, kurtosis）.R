




library(sjPlot)
library(MASS)
library(mvtnorm)
library(lavaan) 
library(CompQuadForm) # Version1.4.2
library(BAMMtools) # Version 2.1.6
library(Ckmeans.1d.dp) # Version 4.0.1
library(expm)
library(Matrix)
library(gtools)  
library(regsem)
library(semTools)

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
generate_cfa_data <- function(n_factors, n, skew = 1, kurt = 3) {
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
m<-n_factors<-3



get_pvalues <- function(n) {
  model <- create_cfa_syntax(n_factors)
  dat<- generate_cfa_data(n_factors = n_factors, n = n)
  fit <- cfa(model, data = dat, estimator = "MLM", test = c("standard", "browne.residual.nt.model"))
  cb <- (1 - (2*ncol(dat) + 4*m + 11)/(6*n)) * n / (n - 1)
  chisq_mlb <- cb * fitmeasures(fit, "chisq")
  chisq_rmlb <- cb * fitmeasures(fit, "chisq.scaled")
  df <- fitmeasures(fit, "df")
  p_mlb <- 1 - pchisq(chisq_mlb, df)                  # 1  pvalue of  TML^B
  p_rmlb <- 1 - pchisq(chisq_rmlb, df)                # 2  pvalue of  TRML^B
  p_rls <- fit@test$browne.residual.nt.model$pvalue   # 3  pvalue of  RLS  (Hayakawa, 2019)
  p_ml <- fitmeasures(fit, "pvalue")                  # 4  pvalue of  TML
  
  
  T_rls <- fit@test$browne.residual.nt.model$stat          ##  RLS  (Hayakawa, 2019)
  CT_rls<-T_rls/fitMeasures(fit, "chisq.scaling.factor")  #    scaled RLS
  p_CTrls <- 1 - pchisq(CT_rls, df)                   #5     pvalue of  scaled RLS
  
  
  #extract test statistic T                            ####  Foldnes Statistics (2017)
  T = fitmeasures(fit, "chisq")
  #Extract U*Gamma
  UG <- lavInspect(fit, "UGamma")
  #The estimated eigenvalues
  #df = fitmeasures(fit,"DF")
  eig_hat = Re(eigen(UG)$values[1:df])
  ########
  ## p values for various tests, based on eigenvalues
  ########
  #NTML
  pNTML <- imhof(T,rep(1, df))$Qq                    # 6
  #SB
  pSB <-imhof(T,rep(mean(eig_hat), df))$Qq          # 7
  #EBAF
  pEBAF <-imhof(T,eig_hat)$Qq                      # 8
  #EBA2
  eigs <- c(rep(mean(eig_hat[1:ceiling(df/2)]), ceiling(df/2)),
            rep(mean(eig_hat[(ceiling(df/2) +1):df]), df - ceiling(df/2))) 
  pEBA2 <- imhof(T, eigs)$Qq                       # 9
  #Jenks EBA2
  breaks <- getJenksBreaks(eig_hat, k = 3)
  block1 <- eig_hat[eig_hat <= breaks [2]] ; block2 = eig_hat[eig_hat > breaks[2]]
  eigs <- c(rep (mean(block1), length(block1)), rep(mean(block2), length(block2)))
  pEBA2J <- imhof(T, eigs)$Qq                    # 10
  #EBAA
  t = Ckmeans.1d.dp(eig_hat)
  means <- t$centers
  clusters <- t$cluster
  eigs <- sapply(clusters, function (x) means [x])
  pEBAA <- imhof(T, eigs) $Qq                  #  11
  
  
  ####  our proposed statistics
  yn<-ncol(dat)/nrow(dat)
  tau<-ncol(dat)-ncol(dat)*(1-1/yn)*log(1-yn)-0.5*log(1-yn)
  v<-sqrt(-2*log(1-yn)-2*yn)
  T_F<-(fitmeasures(fit, "chisq")/nrow(dat)-tau)/v         ##  T_F
  pT_F<-1-pnorm(T_F,0,1)                                   #   pvalue of  T_F
  ##
  mu_y1<-2.015910*(m/ncol(dat))+1.291412*(ncol(dat)/nrow(dat))-0.278377*m+0.036066*ncol(dat)-2.393643  ## g
  T_F_c1<-(fitmeasures(fit, "chisq")/nrow(dat)-tau-mu_y1*v)/v
  pT_F_c1<-1-pnorm(T_F_c1,0,1) 
  
  mu_y2<- -49.62*(ncol(dat)/(nrow(dat))^2)+4.270*(ncol(dat)/nrow(dat))+0.003475*nrow(dat)-0.1819*m-2.626  ## g'
  T_F_c2<-(fitmeasures(fit, "chisq")/nrow(dat)-tau-mu_y2*v)/v
  pT_F_c2<-1-pnorm(T_F_c2,0,1)
  
  ##
  T_CsF<-(fitmeasures(fit, "chisq.scaled")/nrow(dat)-tau)/v   ##  T_CsF   
  pT_CsF<-1-pnorm(T_CsF,0,1)                                 #    pvalue of  T_CsF 
  ##
  mu_y3<-1.531018*(ncol(dat)/nrow(dat))-0.237301*m+0.029336*ncol(dat)-2.035224       ## g1
  T_CsF_c1<-(fitmeasures(fit, "chisq.scaled")/nrow(dat)-tau-mu_y3*v)/v
  pT_CsF_c1<-1-pnorm(T_CsF_c1,0,1) 
  
  mu_y4<--50.735178*(ncol(dat)/(nrow(dat))^2)+4.788150*(ncol(dat)/nrow(dat))+0.003553*nrow(dat)-0.183627*m-2.639546  ## g'1
  T_CsF_c2<-(fitmeasures(fit, "chisq.scaled")/nrow(dat)-tau-mu_y4*v)/v
  pT_CsF_c2<-1-pnorm(T_CsF_c2,0,1) 
  
  ###
  fit1<-cfa(model,data =dat,ridge = TRUE, ridge.constant = 1e-3,estimator="MLM")  ###   ridge
  T_F_cr<-(fitmeasures(fit1, "chisq")/nrow(dat)-tau-mu_y2*v)/v
  pT_F_cr<-1-pnorm(T_F_cr,0,1)  
  T_CsF_cr<-(fitmeasures(fit1, "chisq.scaled")/nrow(dat)-tau-mu_y4*v)/v
  pT_CsF_cr<-1-pnorm(T_CsF_cr,0,1) 
  
  
  
  c(chisq_mlb,chisq_rmlb,T_rls,CT_rls,T_F,T_F_c1,T_F_c2,T_F_cr,T_CsF,T_CsF_c1,T_CsF_c2,T_CsF_cr,p_mlb,p_rmlb,p_rls,p_CTrls,pT_F,pT_F_c1,pT_F_c2,pT_F_cr,pT_CsF,pT_CsF_c1,pT_CsF_c2,pT_CsF_cr,p_ml,pNTML,pSB,pEBAF,pEBA2,pEBA2J,pEBAA)
}

library(future.apply)
plan(multisession)
set.seed(123)
res <- future_replicate(1000, get_pvalues(100))
rowMeans(res[13:31,]< 0.05)
rowMeans(res[1:12,])
###
L<-data.frame(matrix(nrow = 1,ncol = 19))
colnames(L)<-c("p_mlb","p_rmlb","p_rls","p_CTrls","pT_F","pT_F_c1","pT_F_c2","pT_F_cr","pT_CsF","pT_CsF_c1","pT_CsF_c2","pT_CsF_cr","p_ml","pNTML","pSB","pEBAF","pEBA2","pEBA2J","pEBAA")
L[1,]<-rowMeans(res[13:31,]< 0.05)

tab_df(L,digits = 3)









library(future.apply)
plan(multisession)
set.seed(123)

#####
#
failure_count <- 0  

get_pvalues_safe <- function(n) {
  tryCatch({
    get_pvalues(n)  
  }, error = function(e) {
    assign("failure_count", failure_count + 1, envir = .GlobalEnv)
    return(NA)  
  })
}



res <- future_replicate(1000, get_pvalues_safe(100))
valid_results <- res[!sapply(res, function(x) all(is.na(x)))]
failure_count <- sum(sapply(res, function(x) all(is.na(x))))
cat("Number of samples that did not converge：", failure_count, "/1000\n")
result_matrix <- do.call(rbind, valid_results)
cat("Effective result dimension：", dim(result_matrix), "\n")

result_matrix_t<-t(result_matrix)

#####

rowMeans(result_matrix_t[13:31,]< 0.05)
rowMeans(result_matrix_t[1:12,])
###
L<-data.frame(matrix(nrow = 1,ncol = 19))
colnames(L)<-c("p_mlb","p_rmlb","p_rls","p_CTrls","pT_F","pT_F_c1","pT_F_c2","pT_F_cr","pT_CsF","pT_CsF_c1","pT_CsF_c2","pT_CsF_cr","p_ml","pNTML","pSB","pEBAF","pEBA2","pEBA2J","pEBAA")
L[1,]<-rowMeans(result_matrix_t[13:31,]< 0.05)

tab_df(L,digits = 3)











