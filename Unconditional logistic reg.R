# Load analyses dataset 
index <- read.csv("filterd_disease_pairs.csv")  
data <- read.csv("main_dataset.csv")

# Define covariates
covs <- c("sex", "age", "ethnicity", "employment", "tdi", "education",
                "diet_score", "alcohol_status", "smoking_status",
                "sleep_duration", "physical_activity", "bmi")

# Initialize results dataframe -------------------------------------------------
results_df <- data.frame(
  exposure = character(),
  outcome = character(),
  or = numeric(),
  ci_low = numeric(),
  ci_high = numeric(),
  p_value = numeric(),
  fdr_p_value = numeric(),
  stringsAsFactors = FALSE
)

# Main analysis loop -----------------------------------------------------------
for (i in 1:nrow(index)) {
  # Get current exposure-outcome pair
  current_exposure <- index[i, "disease_1"]
  current_outcome <- index[i, "disease_2"]
  
  # Build unconditional logistic regression formula -----------------------------
  model_formula <- paste(
    current_exposure, "~", 
    current_outcome, "+", 
    paste(covariates, collapse = " + ")
  )
  
  # Fit regression model ----------------------------------------------
  model <- glm(
    formula = as.formula(model_formula),
    data = data,
    family = binomial
  )
  
  # Extract model results ------------------------------------------------------
  model_summary <- summary(model)
  coef_index <- 2  # Index for primary predictor
  odds_ratio <- exp(coef(model)[coef_index])
  conf_int <- exp(confint(model)[coef_index, ])
  p_val <- model_summary$coefficients[coef_index, 4]
  
  # Store results --------------------------------------------------------------
  results_df <- rbind(
    results_df,
    data.frame(
      exposure = current_exposure,
      outcome = current_outcome,
      or = round(odds_ratio, 2),
      ci_low = round(conf_int[1], 2),
      ci_high = round(conf_int[2], 2),
      p_value = p_val,
      fdr_q_value = NA
    )
  )
}

# Adjust for multiple comparisons ----------------------------------------------
results_df <- results_df %>% 
  mutate(fdr_p_value = p.adjust(p_value, method = "fdr"))
