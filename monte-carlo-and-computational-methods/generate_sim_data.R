
# ============================================================
# General Simulation Data Generator
# File: /Data Generation/generate_sim_data.R
# Author: Mark Fulginiti
#
# Purpose:
#   Generate flexible synthetic datasets for simulation studies in
#   causal inference, regression, and predictive modeling under
#   controlled departures from ideal assumptions.
#
# Main features:
#   - Multiple outcome types: continuous, binary, count, multinomial
#   - Correlated baseline covariates
#   - Optional nonlinearity and interaction structure
#   - Configurable treatment assignment with confounding and overlap control
#   - Optional heteroskedasticity and measurement error
#   - Optional MAR-style missingness
#   - Returns observed data, metadata, truth, and optionally potential outcomes
#
# Notes:
#   This is the general-purpose simulation engine for the repository.
#   Note-specific scripts should wrap this function and impose only the
#   settings needed for a given technical note.
#
# Main function:
#   generate_sim_data()
# ============================================================



generate_sim_data <- function(
    n = 1000,
    seed = NULL,
    outcome_type = c("continuous", "binary", "count", "multinomial"),
    correlation = c("moderate", "none", "strong"),
    nonlinearity = c("none", "moderate", "strong"),
    interactions = c("none", "moderate", "strong"),
    confounding = c("none", "moderate", "strong"),
    overlap = c("good", "moderate", "poor"),
    heteroskedasticity = c("none", "moderate", "strong"),
    missingness = c("none", "mar_moderate", "mar_strong"),
    measurement_error = c("none", "moderate", "strong"),
    return_potential_outcomes = TRUE
) {
  outcome_type <- match.arg(outcome_type)
  correlation <- match.arg(correlation)
  nonlinearity <- match.arg(nonlinearity)
  interactions <- match.arg(interactions)
  confounding <- match.arg(confounding)
  overlap <- match.arg(overlap)
  heteroskedasticity <- match.arg(heteroskedasticity)
  missingness <- match.arg(missingness)
  measurement_error <- match.arg(measurement_error)
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- 0. Helper functions ---------------------------------------------------
  
  inv_logit <- function(z) 1 / (1 + exp(-z))
  
  softmax_vec <- function(eta_mat) {
    eta_shift <- eta_mat - apply(eta_mat, 1, max)
    exp_eta <- exp(eta_shift)
    exp_eta / rowSums(exp_eta)
  }
  
  # --- 1. Base covariates ----------------------------------------------------
  rho <- switch(
    correlation,
    "none" = 0.0,
    "moderate" = 0.35,
    "strong" = 0.7
  )
  
  Sigma <- matrix(
    c(1, rho, rho,
      rho, 1, rho,
      rho, rho, 1),
    nrow = 3, byrow = TRUE
  )
  
  # Generate 3 correlated continuous covariates
  Z <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  X_cont <- Z %*% chol(Sigma) 
  
  x1 <- as.numeric(X_cont[, 1])
  x2 <- as.numeric(X_cont[, 2])
  x3 <- as.numeric(X_cont[, 3])
  
  x4 <- rbinom(n, size = 1, prob = 0.5)                # binary
  x5 <- sample(0:2, size = n, replace = TRUE,          # ordinal-ish
               prob = c(0.3, 0.45, 0.25))
  x6 <- rnorm(n)                                       # noise variable
  
  X_true <- data.frame(
    id = seq_len(n),
    x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6
  )
  
  # --- 2. Optional measurement error ----------------------------------------
  X_obs <- X_true
  if (measurement_error != "none") {
    me_sd <- switch(
      measurement_error,
      "moderate" = 0.4,
      "strong" = 0.8
    )
    X_obs$x1 <- X_true$x1 + rnorm(n, 0, me_sd)
    X_obs$x2 <- X_true$x2 + rnorm(n, 0, me_sd)
  }
  
  # --- 3. True signal for outcomes ------------------------------------------
  f_base <- 0.8 * x1 - 0.6 * x2 + 0.5 * x3 + 0.7 * x4 - 0.3 * x5
  
  f_nl <- switch(
    nonlinearity,
    "none" = rep(0, n),
    "moderate" = 0.5 * x1^2 - 0.4 * sin(x2),
    "strong" = 0.9 * x1^2 - 0.7 * sin(x2) + 0.5 * pmax(x3, 0)
  )
  
  f_int <- switch(
    interactions,
    "none" = rep(0, n),
    "moderate" = 0.6 * x1 * x4 - 0.4 * x2 * x3,
    "strong" = 1.0 * x1 * x4 - 0.7 * x2 * x3 + 0.5 * x1 * x5
  )
  
  # Untreated latent mean structure
  mu0 <- f_base + f_nl + f_int  
  
  # Heterogeneous treatment effect
  tau <- 1.0 + 0.5 * x1 - 0.4 * x4 + 0.3 * (x5 == 2)
  
  # Potential-outcome means on latent scale
  mu1 <- mu0 + tau
  
  # --- 4. Treatment assignment ----------------------------------------------
  conf_strength <- switch(
    confounding,
    "none" = 0,
    "moderate" = 0.8,
    "strong" = 1.5
  )
  
  overlap_scale <- switch(
    overlap,
    "good" = 0.6,
    "moderate" = 1.2,
    "poor" = 2.5
  )
  
  # If confounding = none, treatment is randomized around 0.5
  lp_treat <- if (confounding == "none") {
    rep(0, n)
  } else {
    overlap_scale * (
      conf_strength * (0.8 * x1 - 0.6 * x2 + 0.5 * x4 - 0.3 * x5)
    )
  }
  
  # Propensity score and treatment assignment
  ps <- inv_logit(lp_treat)
  D <- rbinom(n, size = 1, prob = ps)
  
  # --- 5. Realized potential outcomes ---------------------------------------
  if (outcome_type == "continuous") {
    sigma_eps <- switch(
      heteroskedasticity,
      "none" = rep(1, n),
      "moderate" = exp(0.35 * x1),
      "strong" = exp(0.7 * x1)
    )
    
    Y0 <- mu0 + rnorm(n, 0, sigma_eps)
    Y1 <- mu1 + rnorm(n, 0, sigma_eps)
    
  } else if (outcome_type == "binary") {
    p0 <- inv_logit(mu0)
    p1 <- inv_logit(mu1)
    Y0 <- rbinom(n, size = 1, prob = p0)
    Y1 <- rbinom(n, size = 1, prob = p1)
    
  } else if (outcome_type == "count") {
    lambda0 <- exp(mu0 / 3)
    lambda1 <- exp(mu1 / 3)
    Y0 <- rpois(n, lambda0)
    Y1 <- rpois(n, lambda1)
    
  } else if (outcome_type == "multinomial") {
    
    # 3-category multinomial outcome:
    # class 1 is reference, class 2 and 3 have linear predictors
    
    eta0_mat <- cbind(
      rep(0, n),
      0.6 * mu0,
      -0.4 * mu0 + 0.5 * x2 - 0.3 * x4
    )
    
    eta1_mat <- cbind(
      rep(0, n),
      0.6 * mu1,
      -0.4 * mu1 + 0.5 * x2 - 0.3 * x4
    )
    
    p0_mat <- softmax_vec(eta0_mat)
    p1_mat <- softmax_vec(eta1_mat)
    
    draw_multinom <- function(prob_mat) {
      apply(prob_mat, 1, function(p) sample(1:3, size = 1, prob = p))
    }
    
    Y0 <- draw_multinom(p0_mat)
    Y1 <- draw_multinom(p1_mat)
  }
  
  Y <- ifelse(D == 1, Y1, Y0)
  
  # --- 6. Optional missingness (MAR-ish) ------------------------------------
  miss_y <- rep(FALSE, n)
  miss_x2 <- rep(FALSE, n)
  
  if (missingness != "none") {
    miss_strength <- switch(
      missingness,
      "mar_moderate" = 0.6,
      "mar_strong" = 1.1
    )
    
    p_miss_y <- inv_logit(-2 + miss_strength * (0.8 * x1 + 0.6 * D))
    p_miss_x2 <- inv_logit(-2.2 + miss_strength * (0.7 * x4 + 0.5 * x5))
    
    miss_y <- rbinom(n, 1, p_miss_y) == 1
    miss_x2 <- rbinom(n, 1, p_miss_x2) == 1
  }
  
  dat <- X_obs
  dat$D <- D
  dat$ps_true <- ps
  dat$Y <- Y
  
  dat$Y[miss_y] <- NA
  dat$x2[miss_x2] <- NA
  
  # --- 7. True estimands -----------------------------------------------------
  if (outcome_type == "multinomial") {
    true_ate <- NA_real_
    true_att <- NA_real_
    true_atc <- NA_real_
  } else {
    true_ate <- mean(Y1 - Y0)
    true_att <- mean((Y1 - Y0)[D == 1])
    true_atc <- mean((Y1 - Y0)[D == 0])
  }
  
  out <- list(
    data = dat,
    metadata = list(
      n = n,
      outcome_type = outcome_type,
      correlation = correlation,
      nonlinearity = nonlinearity,
      interactions = interactions,
      confounding = confounding,
      overlap = overlap,
      heteroskedasticity = heteroskedasticity,
      missingness = missingness,
      measurement_error = measurement_error
    ),
    truth = list(
      ate = true_ate,
      att = true_att,
      atc = true_atc
    )
  )
  
  if (return_potential_outcomes) {
    out$full_data <- transform(
      X_true,
      D = D,
      ps_true = ps,
      Y0 = Y0,
      Y1 = Y1,
      Y = ifelse(D == 1, Y1, Y0),
      tau = if (outcome_type == "multinomial") NA_real_ else Y1 - Y0
    )
  }
  
  return(out)
}

