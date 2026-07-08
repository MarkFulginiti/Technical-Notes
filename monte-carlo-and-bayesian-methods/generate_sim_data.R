
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
  # Resolve each user-supplied option to one allowed value.
  # This guards against misspecified arguments and standardizes downstream logic.
  outcome_type <- match.arg(outcome_type)
  correlation <- match.arg(correlation)
  nonlinearity <- match.arg(nonlinearity)
  interactions <- match.arg(interactions)
  confounding <- match.arg(confounding)
  overlap <- match.arg(overlap)
  heteroskedasticity <- match.arg(heteroskedasticity)
  missingness <- match.arg(missingness)
  measurement_error <- match.arg(measurement_error)
  
  # Set a random seed if requested so simulation results are reproducible.  
  if (!is.null(seed)) set.seed(seed)
  
  # ---------------------------------------------------------------------------
  # 0. Helper functions
  # ---------------------------------------------------------------------------
  
  # Logistic inverse link.
  # Used for treatment probabilities, binary outcomes, and missingness models.
  inv_logit <- function(z) 1 / (1 + exp(-z))
  
  # Row-wise softmax transformation.
  # Converts a matrix of category-specific linear predictors into multinomial
  # probabilities that sum to 1 within each row.
  # Subtracting the row maximum improves numerical stability.   
  softmax_vec <- function(eta_mat) {
    eta_shift <- eta_mat - apply(eta_mat, 1, max)
    exp_eta <- exp(eta_shift)
    exp_eta / rowSums(exp_eta)
  }
  
  # ---------------------------------------------------------------------------
  # 1. Generate base covariates
  # ---------------------------------------------------------------------------
  
  # Map the requested correlation level to the pairwise correlation rho
  # used in the 3 x 3 covariance matrix for the continuous covariates.
  rho <- switch(
    correlation,
    "none" = 0.0,
    "moderate" = 0.35,
    "strong" = 0.7
  )
  
  # Covariance matrix for (x1, x2, x3).
  Sigma <- matrix(
    c(1, rho, rho,
      rho, 1, rho,
      rho, rho, 1),
    nrow = 3, byrow = TRUE
  )
  
  # Draw independent standard normals and then induce the target covariance
  # structure using a Cholesky transformation.
  Z <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  X_cont <- Z %*% chol(Sigma) 
  
  # Store each continuous covariate as a numeric vector.  
  x1 <- as.numeric(X_cont[, 1])
  x2 <- as.numeric(X_cont[, 2])
  x3 <- as.numeric(X_cont[, 3])
  
  
  # Additional covariates:
  # x4: binary covariate
  # x5: ordinal-style categorical covariate with three levels
  # x6: auxiliary covariate that is noise for the outcome model but enters
  # the treatment-assignment mechanism to vary overlap independently of confounding
  x4 <- rbinom(n, size = 1, prob = 0.5)                
  x5 <- sample(0:2, size = n, replace = TRUE,        
               prob = c(0.3, 0.45, 0.25))
  x6 <- rnorm(n)                                      
  
  # True covariate data before any measurement error is introduced.  
  X_true <- data.frame(
    id = seq_len(n),
    x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6
  )
  
  # ---------------------------------------------------------------------------
  # 2. Add optional measurement error to observed covariates
  # ---------------------------------------------------------------------------
  
  # Start the observed data as a copy of the truth, then perturb selected
  # covariates if measurement error is requested.
  X_obs <- X_true
  if (measurement_error != "none") {
    # Standard deviation of additive measurement error.
    me_sd <- switch(
      measurement_error,
      "moderate" = 0.4,
      "strong" = 0.8
    )
    
    # Add classical measurement error to x1 and x2 in the observed data only.
    # The true values remain available in X_true and, optionally, full_data.    
    X_obs$x1 <- X_true$x1 + rnorm(n, 0, me_sd)
    X_obs$x2 <- X_true$x2 + rnorm(n, 0, me_sd)
  }
  
  # ---------------------------------------------------------------------------
  # 3. Construct the latent outcome signal
  # ---------------------------------------------------------------------------
  
  # Baseline linear signal shared across outcome families.
  f_base <- 0.8 * x1 - 0.6 * x2 + 0.5 * x3 + 0.7 * x4 - 0.3 * x5
  
  # Optional nonlinear component.
  # Introduces departures from a purely linear signal through quadratic,
  # smooth oscillating, and threshold-type nonlinear terms.
  f_nl <- switch(
    nonlinearity,
    "none" = rep(0, n),
    "moderate" = 0.5 * x1^2 - 0.4 * sin(x2),
    "strong" = 0.9 * x1^2 - 0.7 * sin(x2) + 0.5 * pmax(x3, 0)
  )
  
  # Optional interaction component.
  # This is useful when studying omitted interactions or partial misspecification.  
  f_int <- switch(
    interactions,
    "none" = rep(0, n),
    "moderate" = 0.6 * x1 * x4 - 0.4 * x2 * x3,
    "strong" = 1.0 * x1 * x4 - 0.7 * x2 * x3 + 0.5 * x1 * x5
  )
  
  # Untreated latent mean.
  mu0 <- f_base + f_nl + f_int  
  
  # Heterogeneous treatment effect.
  # Treatment effects vary across units as a function of covariates.
  tau <- 1.0 + 0.5 * x1 - 0.4 * x4 + 0.3 * (x5 == 2)
  
  # Treated latent mean.
  mu1 <- mu0 + tau
  
  # ---------------------------------------------------------------------------
  # 4. Generate treatment assignment
  # ---------------------------------------------------------------------------
  
  # Confounding controls how strongly treatment depends on outcome-relevant
  # covariates. This creates imbalance on variables that also affect outcomes.
  conf_strength <- switch(
    confounding,
    "none" = 0,
    "moderate" = 0.8,
    "strong" = 1.5
  )
  
  # Overlap controls how extreme the treatment probabilities become.
  # Larger values push propensity scores closer to 0 or 1, even when the
  # confounding component is weak or absent.
  overlap_strength <- switch(
    overlap,
    "good" = 0,
    "moderate" = 1.1,
    "poor" = 1.8
  )
  
  # Confounding component:
  # depends on outcome-relevant covariates and therefore induces selection bias
  # when treatment is not randomized.
  lp_conf <- conf_strength * (0.8 * x1 - 0.6 * x2 + 0.5 * x4 - 0.3 * x5)
  
  # Overlap component:
  # depends primarily on a variable excluded from the outcome model, so it can
  # worsen common support without directly creating outcome confounding.
  # Because x6 is mean-zero in expectation, this component perturbs treatment
  # probabilities without systematically shifting them toward treatment or control.
  lp_overlap <- overlap_strength * x6
  
  # Total treatment-assignment linear predictor.
  lp_treat <- lp_conf + lp_overlap
  
  # True propensity score and realized treatment assignment.
  ps <- inv_logit(lp_treat)
  D <- rbinom(n, size = 1, prob = ps)
  
  # ---------------------------------------------------------------------------
  # 5. Generate potential outcomes by outcome family
  # ---------------------------------------------------------------------------
  if (outcome_type == "continuous") {
    # Outcome-specific error scale.
    # Heteroskedasticity is only used for continuous outcomes.    
    # An exponential form is used for the conditional error standard deviation 
    # so that heteroskedasticity varies smoothly with x1 while remaining 
    # strictly positive for all units.
    sigma_eps <- switch(
      heteroskedasticity,
      "none" = rep(1, n),
      "moderate" = exp(0.35 * x1),
      "strong" = exp(0.7 * x1)
    )
    
    # Gaussian potential outcomes around latent means.    
    Y0 <- mu0 + rnorm(n, 0, sigma_eps)
    Y1 <- mu1 + rnorm(n, 0, sigma_eps)
    
  } else if (outcome_type == "binary") {
    # Convert latent means to probabilities using the logistic link.    
    p0 <- inv_logit(mu0)
    p1 <- inv_logit(mu1)
    
    # Bernoulli potential outcomes.    
    Y0 <- rbinom(n, size = 1, prob = p0)
    Y1 <- rbinom(n, size = 1, prob = p1)
    
  } else if (outcome_type == "count") {
    # Convert latent means to Poisson rates.
    # Dividing by 3 tempers the rate scale so counts remain in a useful range
    # across a variety of signal settings.    
    lambda0 <- exp(mu0 / 3)
    lambda1 <- exp(mu1 / 3)
    Y0 <- rpois(n, lambda0)
    Y1 <- rpois(n, lambda1)
    
  } else if (outcome_type == "multinomial") {
    # Three-category multinomial outcome.
    # Category 1 is the reference category with linear predictor 0.    
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
    
    # Convert category-specific linear predictors into probabilities.
    p0_mat <- softmax_vec(eta0_mat)
    p1_mat <- softmax_vec(eta1_mat)
    
    # Draw one category per row according to its multinomial probability vector.    
    draw_multinom <- function(prob_mat) {
      apply(prob_mat, 1, function(p) sample(1:3, size = 1, prob = p))
    }
    
    Y0 <- draw_multinom(p0_mat)
    Y1 <- draw_multinom(p1_mat)
  }
  
  # Realized observed outcome under the switching equation.
  Y <- ifelse(D == 1, Y1, Y0)
  
  # ---------------------------------------------------------------------------
  # 6. Introduce optional missingness
  # ---------------------------------------------------------------------------
  
  # Initialize missingness indicators.
  miss_y <- rep(FALSE, n)
  miss_x2 <- rep(FALSE, n)
  
  if (missingness != "none") {
    # Missingness strength controls how strongly covariates/treatment predict
    # whether Y and x2 are observed.    
    miss_strength <- switch(
      missingness,
      "mar_moderate" = 0.6,
      "mar_strong" = 1.1
    )
    
    # Outcome missingness depends on x1 and treatment.
    # Covariate missingness for x2 depends on x4 and x5.
    p_miss_y <- inv_logit(-2 + miss_strength * (0.8 * x1 + 0.6 * D))
    p_miss_x2 <- inv_logit(-2.2 + miss_strength * (0.7 * x4 + 0.5 * x5))
    
    # Realize missingness indicators.    
    miss_y <- rbinom(n, 1, p_miss_y) == 1
    miss_x2 <- rbinom(n, 1, p_miss_x2) == 1
  }
  
  # Build the observed dataset seen by the analyst.  
  dat <- X_obs
  dat$D <- D
  dat$ps_true <- ps
  dat$Y <- Y
  
  # Apply missingness after the observed data object is assembled.  
  dat$Y[miss_y] <- NA
  dat$x2[miss_x2] <- NA
  
  # ---------------------------------------------------------------------------
  # 7. Compute true causal estimands
  # ---------------------------------------------------------------------------
  
  # For multinomial outcomes, a scalar treatment effect defined as Y1 - Y0
  # is not directly meaningful in the same way as for the other outcome types,
  # so the current implementation returns NA for these summary estimands.
  if (outcome_type == "multinomial") {
    true_ate <- NA_real_
    true_att <- NA_real_
    true_atc <- NA_real_
  } else {
    true_ate <- mean(Y1 - Y0)
    true_att <- mean((Y1 - Y0)[D == 1])
    true_atc <- mean((Y1 - Y0)[D == 0])
  }
  
  # Main return object: observed data, settings, and truth.  
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
  
  
  # Optionally return the underlying full data, including true covariates and
  # both potential outcomes. This is useful for validation and simulation studies.  
  if (return_potential_outcomes) {
    out$full_data <- transform(
      X_true,
      D = D,
      ps_true = ps,
      Y0 = Y0,
      Y1 = Y1,
      Y = ifelse(D == 1, Y1, Y0),
      ite = if (outcome_type == "multinomial") NA_real_ else Y1 - Y0
    )
  }
  
  return(out)
}

