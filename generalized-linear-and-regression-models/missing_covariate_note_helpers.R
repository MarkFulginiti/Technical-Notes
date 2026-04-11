# ============================================================
# Missing Clinical Covariates: helper functions
# Built on top of generate_sim_data()
# ============================================================

# Logistic inverse-link function used to convert a linear predictor
# into a probability for the binary outcome.
inv_logit <- function(z) 1 / (1 + exp(-z))

# ------------------------------------------------------------
# 1. Note-specific wrapper around generate_sim_data()
# ------------------------------------------------------------
# This wrapper uses generate_sim_data() only as a covariate engine.
# We keep the moderate correlation structure from the general DGF,
# but we do not use its generated binary outcome. Instead, we
# rescale the covariates into more application-like values and then
# generate the outcome from the note-specific logistic model.
generate_note_base_data <- function(n = 1000, seed = NULL) {
  sim <- generate_sim_data(
    n = n,
    seed = seed,
    outcome_type = "binary",   # placeholder only; Y from the general DGF is not used
    correlation = "moderate",
    nonlinearity = "none",
    interactions = "none",
    confounding = "none",
    overlap = "good",
    heteroskedasticity = "none",
    missingness = "none",
    measurement_error = "none",
    return_potential_outcomes = FALSE
  )
  
  # Rescale latent covariates into more application-like variables.
  # age is centered around 60 years with moderate spread.
  # lab is centered around 100 in plausible laboratory units.
  # bmi is centered around 25 with moderate spread.
  # smoking remains binary.
  # comorbidity remains an ordinal 0/1/2 burden measure.
  #
  # Smoking prevalence remains relatively high because x4 comes from
  # the general DGF Bernoulli(0.5) generator. For this note, that is
  # interpreted as a somewhat higher-risk clinical cohort rather than
  # a general-population sample.
  dat <- dplyr::transmute(
    sim$data,
    id = id,
    age = round(60 + 12 * x1),
    lab = 100 + 15 * x2,
    bmi = 25 + 4 * x3,
    smoking = x4,
    comorbidity = x5
  )
  
  # Standardize the rescaled continuous covariates for modeling.
  # The note uses the standardized versions for outcome generation,
  # missingness mechanisms, and regression estimation, while the raw
  # variables are retained for application-scale interpretation.
  dat <- dplyr::mutate(
    dat,
    age_std = (age - 60) / 12,
    lab_std = (lab - 100) / 15,
    bmi_std = (bmi - 25) / 4,
    risk_linear_predictor =
      -1.2 +
      0.6 * age_std +
      0.5 * comorbidity +
      0.8 * lab_std +
      0.35 * bmi_std +
      0.55 * smoking,
    y = rbinom(dplyr::n(), size = 1, prob = inv_logit(risk_linear_predictor))
  )
  
  # Truth object for the simulation study.
  # The parameter of primary interest is the standardized lab coefficient
  # in the true logistic regression.
  truth <- list(
    beta_lab = 0.8
  )
  
  list(data = dat, truth = truth)
}

# ------------------------------------------------------------
# 2. Missingness mechanisms for the note
# ------------------------------------------------------------
# This function imposes note-specific missingness on the partially
# observed covariates lab, bmi, and smoking.
#
# MCAR:
#   Missingness is assigned independently of patient characteristics.
#
# MAR:
#   Missingness depends on observed predictors already available in the
#   analytic dataset, primarily age_std and comorbidity.
#
# MNAR:
#   Missingness depends partly on the latent value of the covariate
#   itself through lab_std or bmi_std. This is used only as a stress case,
#   not as a full sensitivity-analysis framework.
#
# Missingness is applied to both the raw variables (lab, bmi) and the
# standardized versions (lab_std, bmi_std) so that the fitted model and
# the stored application-scale variables remain synchronized.
impose_note_missingness <- function(dat,
                                    mechanism = c("MCAR", "MAR", "MNAR"),
                                    severity = c("mild", "moderate", "strong"),
                                    seed = NULL) {
  mechanism <- match.arg(mechanism)
  severity <- match.arg(severity)
  
  # Set a seed if provided so the missingness pattern is reproducible.
  if (!is.null(seed)) set.seed(seed)
  
  # Severity increases the strength of the covariate-driven component
  # in the MAR and MNAR-style missingness models.
  sev <- switch(
    severity,
    mild = 0.5,
    moderate = 0.9,
    strong = 1.3
  )
  
  n <- nrow(dat)
  
  # MCAR benchmark:
  # fixed missingness probabilities, independent of the covariates and outcome.
  if (mechanism == "MCAR") {
    p_lab <- rep(switch(severity, mild = 0.10, moderate = 0.20, strong = 0.35), n)
    p_bmi <- rep(switch(severity, mild = 0.08, moderate = 0.15, strong = 0.25), n)
    p_smoking <- rep(switch(severity, mild = 0.05, moderate = 0.10, strong = 0.18), n)
  }
  
  # MAR setting:
  # Missingness depends only on predictors that are observed in the analysis.
  # The intercept terms (-1.6, -2.0, -2.3) set baseline missingness rates so
  # the probabilities are not unrealistically high for everyone before the
  # severity adjustment is applied.
  #
  # The severity multiplier then strengthens or weakens how much the observed
  # predictors matter. Larger severity means stronger dependence of missingness
  # on patient characteristics.
  #
  # Signs are chosen to create simple, interpretable patterns:
  # - lab missingness decreases as age/comorbidity increase, representing more
  #   selective ordering in lower-risk patients
  # - bmi missingness increases with age/comorbidity, representing less complete
  #   documentation in older or sicker patients
  # - smoking missingness also increases with age/comorbidity, but at a lower
  #   baseline rate
  if (mechanism == "MAR") {
    p_lab <- inv_logit(-1.6 - sev * (0.8 * dat$age_std + 0.9 * dat$comorbidity))
    p_bmi <- inv_logit(-2.0 + sev * (0.7 * dat$age_std + 0.6 * dat$comorbidity))
    p_smoking <- inv_logit(-2.3 + sev * (0.5 * dat$age_std + 0.5 * dat$comorbidity))
  }
  
  # MNAR-style stress case:
  # Missingness depends partly on the covariate's own latent value before that
  # value is set to missing. This is why the code can look recursive: in a
  # simulation we generate the full data first, compute the missingness
  # probabilities using those full values, and only then blank out selected
  # entries.
  #
  # The intercept terms again control the baseline missingness rate, while the
  # severity multiplier controls how strongly the observed predictors contribute.
  # The extra terms involving lab_std or bmi_std are what make this an MNAR-style
  # mechanism rather than MAR, because missingness is no longer explained only by
  # variables that remain observed in the final dataset.
  #
  # This is a stress-test design, not a full MNAR theory model. Its purpose is
  # to show what happens once the missingness process depends on information that
  # the imputation model does not fully observe.
  if (mechanism == "MNAR") {
    p_lab <- inv_logit(-1.6 - sev * (0.6 * dat$age_std + 0.7 * dat$comorbidity) - 1.2 * dat$lab_std)
    p_bmi <- inv_logit(-2.0 + sev * (0.5 * dat$age_std + 0.5 * dat$comorbidity) - 0.9 * dat$bmi_std)
    p_smoking <- inv_logit(-2.3 + sev * (0.4 * dat$age_std + 0.4 * dat$comorbidity) + 1.0 * dat$smoking)
  }
  
  # Truncate probabilities away from 0 and 1 to avoid degenerate edge cases
  # in very extreme covariate configurations.
  miss_lab <- rbinom(n, 1, pmin(pmax(p_lab, 0.001), 0.999)) == 1
  miss_bmi <- rbinom(n, 1, pmin(pmax(p_bmi, 0.001), 0.999)) == 1
  miss_smoking <- rbinom(n, 1, pmin(pmax(p_smoking, 0.001), 0.999)) == 1
  
  out <- dat
  
  # Apply missingness to the raw variables and to the standardized versions
  # used in the fitted regression model.
  out$lab[miss_lab] <- NA
  out$lab_std[miss_lab] <- NA
  
  out$bmi[miss_bmi] <- NA
  out$bmi_std[miss_bmi] <- NA
  
  out$smoking[miss_smoking] <- NA
  
  # Store scenario-level missingness diagnostics as an attribute so later
  # code can summarize average missingness rates and complete-case retention.
  attr(out, "missing_rates") <- c(
    lab = mean(miss_lab),
    bmi = mean(miss_bmi),
    smoking = mean(miss_smoking),
    complete_case_retained = mean(
      complete.cases(out[, c("age_std", "comorbidity", "lab_std", "bmi_std", "smoking", "y")])
    )
  )
  
  out
}

# ------------------------------------------------------------
# 3. Analysis methods
# ------------------------------------------------------------

# Complete-case analysis:
# fit the logistic regression only on observations with no missing
# covariates used in the model. This is the benchmark deletion strategy
# the note is evaluating.
fit_complete_case <- function(dat) {
  cc <- dat |>
    dplyr::select(y, age_std, comorbidity, lab_std, bmi_std, smoking) |>
    dplyr::mutate(
      # Treat smoking as a binary factor so the fitted model handles it
      # explicitly as categorical rather than as a numeric 0/1 covariate.
      smoking = factor(smoking, levels = c(0, 1))
    ) |>
    stats::na.omit()
  
  # Guard against pathological scenarios where too few complete cases
  # remain to fit the model sensibly.
  if (nrow(cc) < 50) {
    return(list(
      est = NA_real_,
      se = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      n_used = nrow(cc),
      converged = FALSE
    ))
  }
  
  fit <- stats::glm(
    y ~ age_std + comorbidity + lab_std + bmi_std + smoking,
    data = cc,
    family = stats::binomial()
  )
  
  # Extract the standardized lab coefficient, since that is the target
  # truth quantity tracked across simulation replicates.
  co <- summary(fit)$coefficients
  est <- unname(co["lab_std", "Estimate"])
  se <- unname(co["lab_std", "Std. Error"])
  
  list(
    est = est,
    se = se,
    lower = est - 1.96 * se,
    upper = est + 1.96 * se,
    n_used = nrow(cc),
    converged = isTRUE(fit$converged)
  )
}

# Multiple imputation analysis:
# Use MICE (multiple imputation by chained equations) to fill in the
# partially missing covariates, then fit the logistic regression in each
# completed dataset and combine the results using Rubin's rules.
#
# In this note:
# - age_std and comorbidity are fully observed, so they are used as predictors
#   in the imputation models but are not themselves imputed.
# - lab_std and bmi_std are continuous and are imputed with predictive mean
#   matching (PMM), which preserves plausible observed values and avoids
#   purely model-based normal imputations.
# - smoking is binary and is imputed with logistic regression.
#
# "Chained equations" means MICE cycles through the incomplete variables one
# at a time, updating each variable conditionally on the others. After several
# iterations, this produces one completed dataset. Repeating that process m
# times produces m completed datasets, which reflect uncertainty from the
# missing values rather than pretending there is only one best fill-in.
fit_multiple_imputation <- function(dat, m = 20, seed = NULL, maxit = 10) {
  work <- dat |>
    dplyr::select(y, age_std, comorbidity, lab_std, bmi_std, smoking) |>
    dplyr::mutate(
      # Treat smoking as a binary factor so MICE uses the intended
      # logistic imputation model rather than treating it as continuous.
      smoking = factor(smoking, levels = c(0, 1))
    )
  
  # Build MICE method and predictor specifications.
  # make.method() starts from defaults, then we override variable-specific
  # choices below. make.predictorMatrix() defines which variables can be used
  # to impute which other variables.
  meth <- mice::make.method(work)
  pred <- mice::make.predictorMatrix(work)
  
  # Tell MICE which variables should and should not be imputed.
  # An empty string means "do not impute this variable."
  #
  # y is left unimputed because this note focuses on missing covariates,
  # not missing outcomes.
  # age_std and comorbidity are fully observed by design.
  #
  # PMM for lab_std and bmi_std:
  # MICE first fits a regression-type model, then matches each missing case
  # to observed cases with similar predicted values, and finally draws one
  # of those observed values as the imputation. This tends to preserve the
  # observed distribution better than plugging in a raw regression prediction.
  #
  # logreg for smoking:
  # MICE fits a logistic model for the binary smoking indicator and draws
  # imputations from the resulting Bernoulli probabilities.
  meth["y"] <- ""
  meth["age_std"] <- ""
  meth["comorbidity"] <- ""
  meth["lab_std"] <- "pmm"
  meth["bmi_std"] <- "pmm"
  meth["smoking"] <- "logreg"
  
  # Run the chained-equations algorithm.
  #
  # m      = number of completed datasets to generate
  # maxit  = number of Gibbs-like cycling iterations within each imputation run
  #
  # Intuitively:
  # 1. Start with rough initial fills for missing values
  # 2. Update lab_std conditional on the other variables
  # 3. Update bmi_std conditional on the latest values
  # 4. Update smoking conditional on the latest values
  # 5. Repeat this cycle maxit times
  # 6. Save the resulting completed dataset
  # 7. Repeat until m completed datasets are produced
  imp <- mice::mice(
    work,
    m = m,
    method = meth,
    predictorMatrix = pred,
    maxit = maxit,
    printFlag = FALSE,
    seed = seed
  )
  
  # Fit the substantive logistic regression separately within each completed
  # dataset. This gives m slightly different coefficient estimates because
  # each completed dataset reflects different plausible missing values.
  fit <- with(
    imp,
    stats::glm(y ~ age_std + comorbidity + lab_std + bmi_std + smoking,
               family = stats::binomial())
  )
  
  # Pool the m fitted models using Rubin's rules.
  # This combines:
  # - within-imputation variability (ordinary model uncertainty inside each fit)
  # - between-imputation variability (how much estimates change across the
  #   m completed datasets)
  #
  # The pooled standard error is therefore larger than what you would get
  # from pretending one imputed dataset were the truth.
  pooled <- mice::pool(fit)
  summ <- summary(pooled)
  
  # Extract the pooled estimate and standard error for the lab coefficient,
  # which is the parameter of primary interest in the note.
  row_lab <- summ[summ$term == "lab_std", , drop = FALSE]
  
  est <- row_lab$estimate
  se <- row_lab$std.error
  
  list(
    est = est,
    se = se,
    lower = est - 1.96 * se,
    upper = est + 1.96 * se,
    n_used = nrow(work),
    converged = TRUE
  )
}

# ------------------------------------------------------------
# 4. One replicate
# ------------------------------------------------------------
# Run one complete simulation replicate for a given mechanism/severity pair.
# This function:
#   1. generates the full dataset
#   2. imposes missingness
#   3. fits complete-case and MI analyses
#   4. stores coefficient estimates and basic diagnostics

run_one_replicate <- function(n = 1000,
                              mechanism = c("MCAR", "MAR", "MNAR"),
                              severity = c("mild", "moderate", "strong"),
                              seed = NULL,
                              mi_m = 20) {
  mechanism <- match.arg(mechanism)
  severity <- match.arg(severity)
  
  # Generate the complete data before any missingness is imposed.
  # This gives the full target cohort for the current replicate.
  full <- generate_note_base_data(n = n, seed = seed)
  
  # Impose the chosen missingness mechanism on the complete data.
  # A seed offset is used so the missingness draw is reproducible but not
  # identical to the seed used for data generation itself.
  dat_miss <- impose_note_missingness(
    dat = full$data,
    mechanism = mechanism,
    severity = severity,
    seed = if (is.null(seed)) NULL else seed + 1000
  )
  
  # Fit the deletion-based benchmark analysis.
  cc <- fit_complete_case(dat_miss)
  
  # Fit the multiple-imputation analysis.
  # A different seed offset is used so the imputation step is also reproducible
  # but remains distinct from the data-generation and missingness steps.
  mi <- fit_multiple_imputation(
    dat_miss,
    m = mi_m,
    seed = if (is.null(seed)) NULL else seed + 2000
  )
  
  # Recover the missingness diagnostics stored by impose_note_missingness().
  # These are later averaged across replicates within each scenario.
  mr <- attr(dat_miss, "missing_rates")
  
  # Identify the rows retained by complete-case analysis so we can compare
  # the retained subset to the full cohort through the risk-score diagnostic.
  # This is one way to quantify the compositional distortion caused by deletion.
  cc_idx <- complete.cases(dat_miss[, c("age_std", "comorbidity", "lab_std", "bmi_std", "smoking", "y")])
  
  # Return one tidy row for this replicate.
  # Each row contains:
  #   - scenario labels (mechanism, severity)
  #   - truth value for the target coefficient
  #   - CC and MI estimates with standard errors and intervals
  #   - missingness rates and retained sample fraction
  #   - a simple diagnostic comparing the full cohort to the CC subset
  tibble::tibble(
    mechanism = mechanism,
    severity = severity,
    truth = full$truth$beta_lab,
    
    cc_est = cc$est,
    cc_se = cc$se,
    cc_lower = cc$lower,
    cc_upper = cc$upper,
    cc_n = cc$n_used,
    
    mi_est = mi$est,
    mi_se = mi$se,
    mi_lower = mi$lower,
    mi_upper = mi$upper,
    mi_n = mi$n_used,
    
    miss_lab = mr["lab"],
    miss_bmi = mr["bmi"],
    miss_smoking = mr["smoking"],
    cc_retained = mr["complete_case_retained"],
    
    # Mean risk score in the full dataset and in the retained CC subset.
    # Their difference is a compact proxy for how deletion changes who remains.
    full_mean_risk = mean(full$data$risk_linear_predictor),
    cc_mean_risk = mean(dat_miss$risk_linear_predictor[cc_idx])
  )
}

# ------------------------------------------------------------
# 5. Repeated simulation
# ------------------------------------------------------------
# Run all mechanism/severity combinations across many replicates.
#
# Seeds are offset by scenario and replicate so the full simulation is
# reproducible while each run remains distinct.
run_simulation <- function(n_sims = 500,
                           n = 1000,
                           mechanisms = c("MCAR", "MAR", "MNAR"),
                           severities = c("mild", "moderate", "strong"),
                           mi_m = 20,
                           base_seed = 20260327) {
  # Create the full grid of scenario settings.
  # For this note, each scenario is defined by one missingness mechanism
  # and one severity level.
  settings <- expand.grid(
    mechanism = mechanisms,
    severity = severities,
    stringsAsFactors = FALSE
  )
  
  # Pre-allocate a list large enough to hold every replicate from every scenario.
  # This is more efficient than growing an object repeatedly inside the loops.
  results <- vector("list", length = nrow(settings) * n_sims)
  k <- 1
  
  # Outer loop over scenarios.
  for (s in seq_len(nrow(settings))) {
    mech <- settings$mechanism[s]
    sev <- settings$severity[s]
    
    # Inner loop over replicates within a scenario.
    for (r in seq_len(n_sims)) {
      # Each scenario-replicate pair gets its own seed derived from base_seed.
      # This makes the entire study reproducible while preventing accidental
      # reuse of the same random-number stream across runs.
      results[[k]] <- run_one_replicate(
        n = n,
        mechanism = mech,
        severity = sev,
        seed = base_seed + 10000 * s + r,
        mi_m = mi_m
      )
      k <- k + 1
    }
  }
  
  # Stack all replicate-level rows into one long data frame.
  # This is the raw simulation output that will later be summarized.
  dplyr::bind_rows(results)
}

# ------------------------------------------------------------
# 6. Performance summaries
# ------------------------------------------------------------
# Aggregate replicate-level output into the main performance table.
# The note focuses on:
#   - mean estimate
#   - bias
#   - RMSE
#   - 95% interval coverage
#   - average retained sample size
#   - average missingness rates
#   - average complete-case risk shift

summarize_performance <- function(res) {
  res |>
    # Create replicate-level derived quantities first.
    # Bias is estimate minus truth.
    # Coverage records whether the nominal 95% interval contains the truth.
    dplyr::mutate(
      cc_bias_rep = cc_est - truth,
      mi_bias_rep = mi_est - truth,
      cc_cover = cc_lower <= truth & cc_upper >= truth,
      mi_cover = mi_lower <= truth & mi_upper >= truth
    ) |>
    # Summarize within each mechanism/severity scenario.
    dplyr::group_by(mechanism, severity) |>
    dplyr::summarise(
      # Truth should be constant within scenario; the mean is just a convenient extractor.
      truth = mean(truth),
      
      # Complete-case performance.
      # cc_mean_est: average estimated lab coefficient
      # cc_bias: average deviation from truth
      # cc_rmse: square-root of average squared bias across replicates
      # cc_coverage: proportion of replicates whose 95% interval covered the truth
      # cc_avg_n: average usable sample size after deletion
      # cc_retained: average proportion of rows retained after deletion
      cc_mean_est = mean(cc_est, na.rm = TRUE),
      cc_bias = mean(cc_bias_rep, na.rm = TRUE),
      cc_rmse = sqrt(mean(cc_bias_rep^2, na.rm = TRUE)),
      cc_coverage = mean(cc_cover, na.rm = TRUE),
      cc_avg_n = mean(cc_n, na.rm = TRUE),
      cc_retained = mean(cc_retained, na.rm = TRUE),
      
      # Multiple-imputation performance.
      # mi_avg_n is the nominal analysis size because MI uses all rows,
      # but it is still included for parallelism with the CC output.
      mi_mean_est = mean(mi_est, na.rm = TRUE),
      mi_bias = mean(mi_bias_rep, na.rm = TRUE),
      mi_rmse = sqrt(mean(mi_bias_rep^2, na.rm = TRUE)),
      mi_coverage = mean(mi_cover, na.rm = TRUE),
      mi_avg_n = mean(mi_n, na.rm = TRUE),
      
      # Average missingness rates by variable across replicates in the scenario.
      miss_lab = mean(miss_lab, na.rm = TRUE),
      miss_bmi = mean(miss_bmi, na.rm = TRUE),
      miss_smoking = mean(miss_smoking, na.rm = TRUE),
      
      # Average shift in mean risk score induced by complete-case deletion.
      # Positive or negative values indicate that the retained CC subset differs
      # systematically from the full cohort rather than being a harmless random subset.
      risk_shift_cc = mean(cc_mean_risk - full_mean_risk, na.rm = TRUE),
      .groups = "drop"
    )
}