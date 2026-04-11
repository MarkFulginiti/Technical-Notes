# Technical Notes

This repository contains compact technical notes, short papers, simulation studies, and reproducible code focused on how statistical methods behave when real-world data violate ideal assumptions. The emphasis is on clear explanation, method comparison, practical consequences, and reproducible analysis for everyday practitioners across applied statistics, causal inference, modeling, computation, probability, and related quantitative workflows.

The notes are organized primarily by field and method, with an additional industry map below highlighting common problem settings across domains.

## Field Map

- [Causal Inference](./causal-inference/)
- [Generalized Linear and Regression Models](./generalized-linear-and-regression-models/)
- [Machine Learning](./machine-learning/)
- [Statistical Learning and Validation](./statistical-learning-and-validation/)
- [Statistical Modeling](./statistical-modeling/)
- [Sampling and Survey Inference](./sampling-and-survey-inference/)
- [Probability and Stochastic Processes](./probability-and-stochastic-processes/)
- [Monte Carlo and Computational Methods](./monte-carlo-and-computational-methods/)

## Industry Map

### Healthcare
Applied problems may include readmission risk, treatment-effect estimation, count modeling for utilization, missing clinical data, survival outcomes, calibration, and related healthcare decision problems.

- [Missing Clinical Covariates](./generalized-linear-and-regression-models/Missing-Clinical-Covariates.pdf)
  - [Supporting R code: missing_covariate_note_helpers.R](./generalized-linear-and-regression-models/missing_covariate_note_helpers.R)
    
### Policy
Applied problems may include program evaluation, treatment heterogeneity, regression adjustment, weighting, difference-in-differences, overlap, selection bias, and related policy analysis workflows.

### Screening / Risk
Applied problems may include fraud screening, sanctions or denied-party screening, threshold selection, class imbalance, feature leakage, calibration, covariate shift, and related risk-modeling workflows.

### Customer Behavior
Applied problems may include churn prediction, promotion response, spend modeling, zero-inflated purchase behavior, interaction effects, and related customer analytics workflows.

### Operations
Applied problems may include queueing behavior, service congestion, workflow bottlenecks, threshold-crossing problems, demand uncertainty, and related operational decision workflows.

### Manufacturing
Applied problems may include defect counts, quality control, contamination, measurement error, mixture structure, maintenance timing, and related industrial modeling workflows.

### Education
Applied problems may include tutoring effects, school-level heterogeneity, performance bands, attendance missingness, effect modification, and related educational measurement workflows.

### Labor
Applied problems may include training effects, wage modeling, hiring-related selection issues, earnings heterogeneity, measurement error, and related labor-market workflows.

### Public-Sector Measurement
Applied problems may include survey design, weighting, stratification, clustering, nonresponse bias, frame mismatch, domain estimation, and related official-statistics workflows.


## Foundational Notes
This section highlights core reference notes, reusable simulation infrastructure, and theory-oriented technical material that supports multiple applied domains. These notes are not tied to a single industry setting; instead, they provide shared methodological foundations for later application papers.

Examples include identifiability foundations, data-generating functions, estimator reference notes, validation diagnostics, and other reusable technical building blocks.

- [Data Generating Function for Simulation Studies](./monte-carlo-and-computational-methods/Data_Generating_Function.pdf)
  - [Supporting R code: generate_sim_data.R](./monte-carlo-and-computational-methods/generate_sim_data.R)
- [Identifiability Foundations for Causal Inference](./causal-inference/Identifiability_Foundations.pdf)


## Notes

- The repository is designed to grow as new technical notes and application areas are added.
- Notes are intended to be modular, reusable, and easy to reference across related projects.
- When appropriate, notes will include both R and Python implementations.
