# ============================================================
# THEORY-PRACTICE COMPARISON
# P+ LIVSCYKLUS VERSUS THE THEORETICAL OPTIMUM
# ============================================================
#
# This script must be run after the Chapter 5 calibration script.
#
# It uses the final five-annual-premium calibration sample
# gamma_5yr and the calibrated coefficients of relative risk
# aversion gamma_hat_i.
#
# Notation:
#
#   X_i             = current pension wealth
#   h_i             = present value of remaining premiums
#   gamma_hat_i     = calibrated coefficient of relative
#                     risk aversion
#   theta_i         = product-assigned P+ portfolio
#   w_i^M           = Merton risky-asset allocation
#   w_i^*           = theoretical optimal risky-asset allocation
#
# Before retirement:
#
#   w_i^M =
#     (1 / gamma_hat_i)
#     Sigma^{-1}(alpha - r 1)
#
# and
#
#   w_i^* =
#     (1 + h_i / X_i) w_i^M.
#
# After retirement:
#
#   h(t) = 0
#
# and therefore
#
#   w^*(t,x) = w^M.
#
# The post-retirement theoretical benchmark also uses
#
#   c^*(t,x) = x / abar(t).
#
# Retirement boundary in this thesis:
#
#   T = 67.
#
# ============================================================


# ============================================================
# 1. REQUIRED OBJECTS
# ============================================================

required_objects <- c(
  "gamma_5yr",
  "returns",
  "variance_matrix",
  "pf",
  "retirement_age",
  "calibration_r",
  "figure_dir",
  "results_dir"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "Run the Chapter 5 calibration script first. Missing objects: ",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
}


# ------------------------------------------------------------
# Confirm thesis retirement boundary
# ------------------------------------------------------------

if (!isTRUE(all.equal(retirement_age, 67))) {
  stop(
    "The theory-practice comparison requires the thesis retirement boundary T = 67."
  )
}


# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  results_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. COMMON MARKET QUANTITIES
# ============================================================

# The comparison uses the same market parameters as the
# Chapter 5 inverse calibration.

r <- calibration_r

alpha <- as.numeric(
  returns(1)
)

Sigma <- variance_matrix(1)

n_assets <- length(alpha)

one_vec <- rep(
  1,
  n_assets
)

excess_return <- as.numeric(
  alpha -
    r * one_vec
)

Sigma_inv <- solve(
  Sigma
)


# ------------------------------------------------------------
# Squared market price of risk
# ------------------------------------------------------------
#
# Theta^2 =
#   (alpha - r 1)' Sigma^{-1} (alpha - r 1)

Theta_sq <- as.numeric(
  t(excess_return) %*%
    Sigma_inv %*%
    excess_return
)


# ------------------------------------------------------------
# Matrix square root used in simulations
# ------------------------------------------------------------

L <- t(
  chol(
    Sigma
  )
)


# ============================================================
# 3. THEORETICAL PORTFOLIO FUNCTIONS
# ============================================================

# ------------------------------------------------------------
# Merton allocation
# ------------------------------------------------------------
#
# w^M =
#   (1 / gamma)
#   Sigma^{-1}(alpha - r 1)

merton_allocation <- function(gamma) {
  
  if (
    !is.finite(gamma) ||
    gamma <= 0
  ) {
    return(
      rep(
        NA_real_,
        n_assets
      )
    )
  }
  
  as.numeric(
    (1 / gamma) *
      (
        Sigma_inv %*%
          excess_return
      )
  )
}


# ------------------------------------------------------------
# Optimal pre-retirement risky-asset allocation
# ------------------------------------------------------------
#
# w^*(t,x) =
#   (1 + h(t)/x) w^M

optimal_pre_retirement_allocation <- function(
    X,
    h,
    gamma
) {
  
  if (
    !is.finite(X) ||
    X <= 0 ||
    !is.finite(h) ||
    h < 0 ||
    !is.finite(gamma) ||
    gamma <= 0
  ) {
    
    return(
      rep(
        NA_real_,
        n_assets
      )
    )
  }
  
  w_M <- merton_allocation(
    gamma
  )
  
  as.numeric(
    (
      1 +
        h / X
    ) *
      w_M
  )
}


# ------------------------------------------------------------
# Aggregate asset classes 5-10
# ------------------------------------------------------------

growth_asset_weight <- function(weights) {
  
  if (
    length(weights) != 10 ||
    any(!is.finite(weights))
  ) {
    return(
      NA_real_
    )
  }
  
  sum(
    weights[5:10]
  )
}


# ============================================================
# 3A. FEASIBLE PORTFOLIO PROJECTION FOR ROBUSTNESS
# ============================================================
#
# Chapters 2-4 derive an unconstrained theoretical optimum.
# That unconstrained solution remains the primary benchmark.
#
# For the product-design realism discussion, this function constructs
# a separate long-only, no-borrowing projection satisfying
#
#   w_j >= 0,     sum_j w_j <= 1.
#
# If the non-negative theoretical weights already sum to at most one,
# they are retained. Otherwise the vector is projected onto the unit
# simplex. This is a feasibility projection only. It is NOT presented
# as the solution of a newly constrained stochastic-control problem.

project_to_simplex <- function(weights) {
  w <- as.numeric(weights)
  if (any(!is.finite(w))) return(rep(NA_real_, length(w)))
  w_nonnegative <- pmax(w, 0)
  if (sum(w_nonnegative) <= 1) return(w_nonnegative)
  u <- sort(w_nonnegative, decreasing = TRUE)
  cssv <- cumsum(u)
  rho_candidates <- which(
    u > (cssv - 1) / seq_along(u)
  )
  rho_index <- max(rho_candidates)
  tau <- (cssv[rho_index] - 1) / rho_index
  pmax(w_nonnegative - tau, 0)
}


# ============================================================
# 4. FINAL CALIBRATION SAMPLE
# ============================================================

comparison_sample <- gamma_5yr %>%
  filter(
    is.finite(gamma),
    gamma > 0,
    is.finite(depot),
    depot > 0,
    is.finite(future_premiums)
  ) %>%
  mutate(
    X_i = depot,
    h_i = future_premiums,
    gamma_hat_i = gamma,
    ell_i = annual_premium
  )


cat(
  "\n==================================================\n"
)

cat(
  "THEORY-PRACTICE COMPARISON\n"
)

cat(
  "==================================================\n"
)

cat(
  "Retirement boundary T:",
  retirement_age,
  "\n"
)

cat(
  "Risk-free rate r:",
  r,
  "\n"
)

cat(
  "Final calibrated sample:",
  nrow(comparison_sample),
  "\n"
)


# ============================================================
# 5. PRE-RETIREMENT SAMPLE
# ============================================================

# The pre-retirement problem is defined for age < T.

pre_retirement_sample <- comparison_sample %>%
  filter(
    alder < retirement_age
  )


cat(
  "Pre-retirement comparison sample (age < T):",
  nrow(pre_retirement_sample),
  "\n"
)


# ============================================================
# 6. POLICY-LEVEL PORTFOLIO COMPARISON
# ============================================================

# For each member i:
#
#   theta_i = product-assigned P+ allocation
#
#   w_i^M =
#     (1 / gamma_hat_i)
#     Sigma^{-1}(alpha - r 1)
#
#   w_i^* =
#     (1 + h_i/X_i) w_i^M.
#
# The comparison is evaluated at the member's recorded age.

pre_retirement_sample$growth_weight_pplus <- NA_real_
pre_retirement_sample$growth_weight_optimal <- NA_real_

pre_retirement_sample$
  total_risky_weight_optimal <- NA_real_

pre_retirement_sample$
  riskfree_weight_optimal <- NA_real_

pre_retirement_sample$
  future_premium_multiplier <- NA_real_

pre_retirement_sample$
  h_over_X <- NA_real_

pre_retirement_sample$
  allocation_gap_5_10 <- NA_real_

pre_retirement_sample$
  growth_weight_feasible <- NA_real_

pre_retirement_sample$
  total_risky_weight_feasible <- NA_real_

pre_retirement_sample$
  riskfree_weight_feasible <- NA_real_


# Store complete theoretical portfolios for diagnostics.

optimal_weight_matrix <- matrix(
  NA_real_,
  nrow = nrow(pre_retirement_sample),
  ncol = n_assets
)

pplus_weight_matrix <- matrix(
  NA_real_,
  nrow = nrow(pre_retirement_sample),
  ncol = n_assets
)

feasible_weight_matrix <- matrix(
  NA_real_,
  nrow = nrow(pre_retirement_sample),
  ncol = n_assets
)


for (
  i in seq_len(
    nrow(
      pre_retirement_sample
    )
  )
) {
  
  age_i <-
    pre_retirement_sample$alder[i]
  
  X_i <-
    pre_retirement_sample$X_i[i]
  
  h_i <-
    pre_retirement_sample$h_i[i]
  
  gamma_hat_i <-
    pre_retirement_sample$gamma_hat_i[i]
  
  profile_i <-
    pre_retirement_sample$investeringsProfil[i]
  
  
  # ----------------------------------------------------------
  # Product-assigned P+ portfolio
  # ----------------------------------------------------------
  
  theta_i <- as.numeric(
    pf(
      age_i,
      profile_i
    )
  )
  
  
  # ----------------------------------------------------------
  # Theoretical Merton allocation
  # ----------------------------------------------------------
  
  w_M_i <- merton_allocation(
    gamma_hat_i
  )
  
  
  # ----------------------------------------------------------
  # Theoretical optimal pre-retirement allocation
  # ----------------------------------------------------------
  
  w_star_i <-
    optimal_pre_retirement_allocation(
      X = X_i,
      h = h_i,
      gamma = gamma_hat_i
    )
  
  # Feasible long-only/no-borrowing projection used only as a
  # robustness and product-design diagnostic.
  w_star_feasible_i <- project_to_simplex(w_star_i)
  
  
  # ----------------------------------------------------------
  # Store complete portfolios
  # ----------------------------------------------------------
  
  pplus_weight_matrix[i, ] <-
    theta_i
  
  optimal_weight_matrix[i, ] <-
    w_star_i
  
  feasible_weight_matrix[i, ] <-
    w_star_feasible_i
  
  
  # ----------------------------------------------------------
  # Stocks and illiquid asset weights, asset classes 5-10
  # ----------------------------------------------------------
  
  pre_retirement_sample$
    growth_weight_pplus[i] <-
    growth_asset_weight(
      theta_i
    )
  
  pre_retirement_sample$
    growth_weight_optimal[i] <-
    growth_asset_weight(
      w_star_i
    )
  
  pre_retirement_sample$
    growth_weight_feasible[i] <-
    growth_asset_weight(
      w_star_feasible_i
    )
  
  pre_retirement_sample$h_over_X[i] <-
    h_i / X_i
  
  pre_retirement_sample$allocation_gap_5_10[i] <-
    pre_retirement_sample$growth_weight_optimal[i] -
    pre_retirement_sample$growth_weight_pplus[i]
  
  
  # ----------------------------------------------------------
  # Total theoretical risky exposure
  # ----------------------------------------------------------
  
  pre_retirement_sample$
    total_risky_weight_optimal[i] <-
    sum(
      w_star_i
    )
  
  
  # ----------------------------------------------------------
  # Implied theoretical risk-free allocation
  # ----------------------------------------------------------
  
  pre_retirement_sample$
    riskfree_weight_optimal[i] <-
    1 -
    sum(
      w_star_i
    )
  
  pre_retirement_sample$
    total_risky_weight_feasible[i] <-
    sum(w_star_feasible_i)
  
  pre_retirement_sample$
    riskfree_weight_feasible[i] <-
    1 - sum(w_star_feasible_i)
  
  
  # ----------------------------------------------------------
  # Future-premium multiplier
  # ----------------------------------------------------------
  
  pre_retirement_sample$
    future_premium_multiplier[i] <-
    1 +
    h_i / X_i
}

# ============================================================
# 6A. CALIBRATION-CANCELLATION CHECK
# ============================================================
#
# The inverse calibration gives
#
#   gamma_hat_i
#     =
#   M_i * (1 + h_i / X_i),
#
# where
#
#   M_i
#     =
#   (alpha_theta_i - r) / sigma_theta_i^2.
#
# Substitution into the theoretical allocation gives
#
#   w_i^*
#     =
#   (1 / M_i) Sigma^{-1}(alpha - r1).
#
# Therefore, at the calibration date, h_i / X_i cancels from the
# recomputed unconstrained theoretical allocation. This check confirms
# that the numerical implementation satisfies the algebraic identity.

calibration_M <-
  (
    pre_retirement_sample$portfolio_return -
      r
  ) /
  pre_retirement_sample$portfolio_variance

market_direction <-
  as.numeric(
    Sigma_inv %*%
      excess_return
  )

cancelled_weight_matrix <-
  outer(
    1 / calibration_M,
    market_direction
  )

maximum_cancellation_error <-
  max(
    abs(
      optimal_weight_matrix -
        cancelled_weight_matrix
    ),
    na.rm = TRUE
  )

cat(
  "Maximum calibration-cancellation error:",
  format(
    maximum_cancellation_error,
    scientific = TRUE
  ),
  "\n"
)

if (
  !is.finite(
    maximum_cancellation_error
  ) ||
  maximum_cancellation_error > 1e-10
) {
  
  stop(
    "The recomputed theoretical allocation fails the calibration-cancellation check."
  )
}

# ============================================================
# 7. THEORETICAL PORTFOLIO DIAGNOSTICS
# ============================================================

cat(
  "\nTheoretical portfolios with total risky weight > 1:",
  sum(
    pre_retirement_sample$
      total_risky_weight_optimal > 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Theoretical portfolios with negative risk-free weight:",
  sum(
    pre_retirement_sample$
      riskfree_weight_optimal < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum theoretical risk-free weight:",
  round(
    min(
      pre_retirement_sample$
        riskfree_weight_optimal,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Maximum theoretical total risky weight:",
  round(
    max(
      pre_retirement_sample$
        total_risky_weight_optimal,
      na.rm = TRUE
    ),
    4
  ),
  "\n"
)

cat(
  "Feasible projections with total risky weight > 1:",
  sum(
    pre_retirement_sample$total_risky_weight_feasible > 1 + 1e-10,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Feasible projections with negative risk-free weight:",
  sum(
    pre_retirement_sample$riskfree_weight_feasible < -1e-10,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 8. AGE-BAND SUMMARY
# ============================================================

pre_retirement_sample$age_band <- cut(
  pre_retirement_sample$alder,
  breaks = c(
    50,
    55,
    60,
    retirement_age
  ),
  include.lowest = TRUE,
  right = FALSE
)


portfolio_comparison_summary <-
  pre_retirement_sample %>%
  
  group_by(
    investeringsProfil,
    age_band
  ) %>%
  
  summarise(
    
    N = n(),
    
    Mean_PPlus_Growth_Weight =
      mean(
        growth_weight_pplus,
        na.rm = TRUE
      ),
    
    Mean_Optimal_Growth_Weight =
      mean(
        growth_weight_optimal,
        na.rm = TRUE
      ),
    
    Mean_Feasible_Growth_Weight =
      mean(
        growth_weight_feasible,
        na.rm = TRUE
      ),
    
    Mean_Allocation_Gap_5_10 =
      mean(
        allocation_gap_5_10,
        na.rm = TRUE
      ),
    
    Mean_Optimal_Riskfree_Weight =
      mean(
        riskfree_weight_optimal,
        na.rm = TRUE
      ),
    
    Mean_Future_Premium_Multiplier =
      mean(
        future_premium_multiplier,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


print(
  portfolio_comparison_summary
)


# ============================================================
# 9. SAVE POLICY-LEVEL PORTFOLIO WEIGHTS
# ============================================================

colnames(
  pplus_weight_matrix
) <- paste0(
  "theta_",
  1:n_assets
)

colnames(
  optimal_weight_matrix
) <- paste0(
  "w_star_",
  1:n_assets
)

colnames(
  feasible_weight_matrix
) <- paste0(
  "w_feasible_",
  1:n_assets
)


portfolio_weight_output <- cbind(
  
  pre_retirement_sample %>%
    select(
      policeNr,
      alder,
      investeringsProfil,
      X_i,
      h_i,
      gamma_hat_i,
      future_premium_multiplier,
      h_over_X,
      growth_weight_pplus,
      growth_weight_optimal,
      growth_weight_feasible,
      allocation_gap_5_10,
      total_risky_weight_optimal,
      riskfree_weight_optimal,
      total_risky_weight_feasible,
      riskfree_weight_feasible
    ),
  
  as.data.frame(
    pplus_weight_matrix
  ),
  
  as.data.frame(
    optimal_weight_matrix
  ),
  
  as.data.frame(
    feasible_weight_matrix
  )
)


write.csv(
  portfolio_weight_output,
  file.path(
    results_dir,
    "theory_practice_pre_retirement_portfolios.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    portfolio_comparison_summary
  ),
  file.path(
    results_dir,
    "theory_practice_portfolio_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. PORTFOLIO COMPARISON FIGURE
# ============================================================

portfolio_plot_data <- rbind(
  
  data.frame(
    Age =
      pre_retirement_sample$alder,
    
    Profile =
      pre_retirement_sample$investeringsProfil,
    
    Growth_Weight =
      pre_retirement_sample$
      growth_weight_pplus,
    
    Strategy =
      "Product-assigned P+ strategy"
  ),
  
  data.frame(
    Age =
      pre_retirement_sample$alder,
    
    Profile =
      pre_retirement_sample$investeringsProfil,
    
    Growth_Weight =
      pre_retirement_sample$
      growth_weight_optimal,
    
    Strategy =
      "Theoretical optimal strategy"
  )
)


portfolio_plot_data$Profile_Label <- recode(
  portfolio_plot_data$Profile,
  hoj = "Low risk aversion",
  mellem = "Moderate risk aversion",
  lav = "High risk aversion"
)

portfolio_plot_data$Profile_Label <- factor(
  portfolio_plot_data$Profile_Label,
  levels = c(
    "Low risk aversion",
    "Moderate risk aversion",
    "High risk aversion"
  )
)

portfolio_plot_data$Strategy <- factor(
  portfolio_plot_data$Strategy,
  levels = c(
    "Product-assigned P+ strategy",
    "Theoretical optimal strategy"
  )
)


portfolio_comparison_plot <- ggplot(
  portfolio_plot_data,
  aes(
    x = Age,
    y = Growth_Weight,
    shape = Strategy
  )
) +
  
  geom_point(
    alpha = 0.65,
    size = 1.8,
    position = position_jitter(
      width = 0.05,
      height = 0
    )
  ) +
  
  facet_wrap(
    ~ Profile_Label
  ) +
  
  labs(
    x = "Age",
    y = TeX("$w_{5-10}$"),
    shape = NULL
  ) +
  
  scale_y_continuous(
    breaks = scales::pretty_breaks(
      n = 7
    ),
    labels = function(x) {
      sprintf(
        "%.2f",
        x
      )
    }
  ) +
  
  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 7
    )
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "bottom"
  )


portfolio_comparison_plot


ggsave(
  filename = file.path(
    figure_dir,
    "Theory vs PPlus Growth-Asset Weight.png"
  ),
  plot = portfolio_comparison_plot,
  width = 9,
  height = 5,
  dpi = 300
)


# ============================================================
# 10A. WEALTH-RATIO DIAGNOSTIC
# ============================================================
#
# In the general pre-retirement model, the optimal allocation is
# wealth dependent through h_i / X_i when gamma is held fixed.
#
# In the present empirical comparison, gamma_hat_i is calibrated from
# the same observed P+ portfolio and the same value of h_i / X_i that
# are subsequently used to calculate w_i^*. Consequently, h_i / X_i
# cancels algebraically from the recomputed allocation at the
# calibration date.
#
# This figure therefore reports a descriptive cross-sectional
# association between h_i / X_i and the theory-practice allocation
# gap. It does not identify a direct causal effect of h_i / X_i on
# the initial recomputed theoretical allocation. Any association can
# also reflect differences in age, P+ profile and assigned portfolio.

state_dependence_results <- pre_retirement_sample %>%
  filter(
    is.finite(h_over_X),
    is.finite(allocation_gap_5_10)
  ) %>%
  select(
    policeNr,
    alder,
    investeringsProfil,
    X_i,
    h_i,
    h_over_X,
    future_premium_multiplier,
    growth_weight_pplus,
    growth_weight_optimal,
    allocation_gap_5_10
  )

state_dependence_summary <- state_dependence_results %>%
  summarise(
    N = n(),
    Mean_h_over_X = mean(h_over_X, na.rm = TRUE),
    Median_h_over_X = median(h_over_X, na.rm = TRUE),
    Q25_h_over_X = as.numeric(quantile(h_over_X, 0.25, na.rm = TRUE)),
    Q75_h_over_X = as.numeric(quantile(h_over_X, 0.75, na.rm = TRUE)),
    Min_h_over_X = min(h_over_X, na.rm = TRUE),
    Max_h_over_X = max(h_over_X, na.rm = TRUE),
    Mean_Allocation_Gap_5_10 = mean(allocation_gap_5_10, na.rm = TRUE),
    Median_Allocation_Gap_5_10 = median(allocation_gap_5_10, na.rm = TRUE),
    Mean_Absolute_Allocation_Gap_5_10 = mean(abs(allocation_gap_5_10), na.rm = TRUE)
  )

print(state_dependence_summary)

write.csv(
  state_dependence_results,
  file.path(results_dir, "theory_practice_state_dependence.csv"),
  row.names = FALSE
)

write.csv(
  state_dependence_summary,
  file.path(results_dir, "theory_practice_state_dependence_summary.csv"),
  row.names = FALSE
)

state_dependence_plot <- ggplot(
  state_dependence_results,
  aes(
    x = h_over_X,
    y = allocation_gap_5_10
  )
) +
  geom_point(
    alpha = 0.70,
    size = 2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = TeX("$h_i/X_i$"),
    y = TeX("$w^{*}_{5-10} - w^{P+}_{5-10}$")
  ) +
  scale_x_continuous(
    breaks = seq(0, 2.5, by = 0.5),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_y_continuous(
    breaks = seq(-0.4, 0, by = 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_bw()

state_dependence_plot

ggsave(
  filename = file.path(
    figure_dir,
    "Theory-Practice Allocation Gap by Future Premiums Relative to Wealth.png"
  ),
  plot = state_dependence_plot,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ============================================================
# 10B. FEASIBLE-PROJECTION ROBUSTNESS COMPARISON
# ============================================================
#
# This compares the observed P+ allocation with both the primary
# unconstrained theoretical optimum and the feasible projection.
# The feasible projection is a product-design diagnostic only.

feasibility_plot_data <- rbind(
  data.frame(
    Age = pre_retirement_sample$alder,
    Profile = pre_retirement_sample$investeringsProfil,
    Weight_5_10 = pre_retirement_sample$growth_weight_pplus,
    Strategy = "Product-assigned P+ strategy"
  ),
  data.frame(
    Age = pre_retirement_sample$alder,
    Profile = pre_retirement_sample$investeringsProfil,
    Weight_5_10 = pre_retirement_sample$growth_weight_optimal,
    Strategy = "Theoretical unconstrained strategy"
  ),
  data.frame(
    Age = pre_retirement_sample$alder,
    Profile = pre_retirement_sample$investeringsProfil,
    Weight_5_10 = pre_retirement_sample$growth_weight_feasible,
    Strategy = "Feasible theoretical projection"
  )
)

feasibility_plot_data$Profile_Label <- recode(
  feasibility_plot_data$Profile,
  hoj = "Low risk aversion",
  mellem = "Moderate risk aversion",
  lav = "High risk aversion"
)

feasibility_plot_data$Profile_Label <- factor(
  feasibility_plot_data$Profile_Label,
  levels = c(
    "Low risk aversion",
    "Moderate risk aversion",
    "High risk aversion"
  )
)

feasibility_plot_data$Age_Band <- floor(feasibility_plot_data$Age)

feasibility_age_summary <- feasibility_plot_data %>%
  group_by(Profile_Label, Age_Band, Strategy) %>%
  summarise(
    Median_Weight_5_10 = median(Weight_5_10, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

feasibility_comparison_plot <- ggplot(
  feasibility_age_summary,
  aes(
    x = Age_Band,
    y = Median_Weight_5_10,
    linetype = Strategy,
    group = Strategy
  )
) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 1.6
  ) +
  facet_wrap(
    ~ Profile_Label
  ) +
  labs(
    x = "Observed age",
    y = TeX("$w_{5-10}$"),
    linetype = NULL
  ) +
  scale_x_continuous(
    breaks = c(50, 55, 60, 65),
    labels = c("50", "55", "60", "65")
  ) +
  scale_y_continuous(
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
    labels = c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7")
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

feasibility_comparison_plot

ggsave(
  filename = file.path(
    figure_dir,
    "Observed and Theoretical Stocks and Illiquid Asset Allocation.png"
  ),
  plot = feasibility_comparison_plot,
  width = 9,
  height = 5.5,
  dpi = 300
)

write.csv(
  feasibility_age_summary,
  file.path(results_dir, "theory_practice_feasible_projection_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 10C. TERMINAL P+ PORTFOLIO CALIBRATION AT AGE 70
# ============================================================
#
# The supplied P+ glide path becomes constant from age 70.
#
# This section does NOT use individual policy observations aged
# 70 or above. Instead, it asks a product-level question:
#
#   If the fixed age-70 P+ portfolio for each profile were treated
#   as one customised mutual fund, what coefficient of relative
#   risk aversion would make full investment in that fund optimal
#   under the post-retirement Merton model?
#
# Since future premiums have ended,
#
#   h = 0,
#
# so the customised-fund calibration becomes
#
#   gamma_hat_70
#     =
#   (alpha_theta - r) / sigma_theta^2.
#
# This removes the policy-specific h_i / X_i term that enters the
# pre-retirement inverse calibration.


terminal_pplus_age <- 70

terminal_profiles <- c(
  "hoj",
  "mellem",
  "lav"
)

terminal_profile_labels <- c(
  hoj = "Low risk aversion",
  mellem = "Moderate risk aversion",
  lav = "High risk aversion"
)


# ------------------------------------------------------------
# Construct fixed terminal P+ portfolios
# ------------------------------------------------------------

terminal_pplus_weight_matrix <- t(
  sapply(
    terminal_profiles,
    function(profile) {
      
      as.numeric(
        pf(
          terminal_pplus_age,
          profile
        )
      )
    }
  )
)


rownames(
  terminal_pplus_weight_matrix
) <- terminal_profiles

colnames(
  terminal_pplus_weight_matrix
) <- paste0(
  "Asset_",
  seq_len(
    n_assets
  )
)


# ------------------------------------------------------------
# Calibrate profile-specific terminal gamma
# ------------------------------------------------------------

terminal_gamma_results <- data.frame(
  
  Profile =
    terminal_profiles,
  
  Profile_Label =
    unname(
      terminal_profile_labels[
        terminal_profiles
      ]
    ),
  
  Age =
    terminal_pplus_age,
  
  PPlus_Expected_Return =
    NA_real_,
  
  PPlus_Excess_Return =
    NA_real_,
  
  PPlus_Variance =
    NA_real_,
  
  PPlus_Volatility =
    NA_real_,
  
  Gamma_Hat_70 =
    NA_real_,
  
  stringsAsFactors = FALSE
)


for (
  i in seq_along(
    terminal_profiles
  )
) {
  
  theta_i <-
    as.numeric(
      terminal_pplus_weight_matrix[
        i,
      ]
    )
  
  
  alpha_theta_i <-
    as.numeric(
      sum(
        theta_i *
          alpha
      )
    )
  
  
  excess_theta_i <-
    alpha_theta_i -
    r
  
  
  variance_theta_i <-
    as.numeric(
      t(theta_i) %*%
        Sigma %*%
        theta_i
    )
  
  
  volatility_theta_i <-
    sqrt(
      variance_theta_i
    )
  
  
  gamma_hat_70_i <-
    excess_theta_i /
    variance_theta_i
  
  
  terminal_gamma_results$
    PPlus_Expected_Return[i] <-
    alpha_theta_i
  
  terminal_gamma_results$
    PPlus_Excess_Return[i] <-
    excess_theta_i
  
  terminal_gamma_results$
    PPlus_Variance[i] <-
    variance_theta_i
  
  terminal_gamma_results$
    PPlus_Volatility[i] <-
    volatility_theta_i
  
  terminal_gamma_results$
    Gamma_Hat_70[i] <-
    gamma_hat_70_i
}


cat(
  "\n==================================================\n"
)

cat(
  "TERMINAL P+ PORTFOLIO CALIBRATION AT AGE 70\n"
)

cat(
  "==================================================\n"
)

print(
  terminal_gamma_results
)


write.csv(
  terminal_gamma_results,
  file.path(
    results_dir,
    "terminal_pplus_age70_implied_risk_aversion.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10D. REVERSE TEN-ASSET TEST USING TERMINAL GAMMA
# ============================================================
#
# The terminal calibration above treats each P+ portfolio as one
# customised mutual fund and finds gamma_hat_70 such that full
# investment in that fund is optimal.
#
# This section reverses the exercise.
#
# Holding gamma_hat_70 fixed, return to the unrestricted
# ten-asset Merton problem:
#
#   w^M
#     =
#   (1 / gamma_hat_70)
#   Sigma^{-1}(alpha - r 1).
#
# The resulting theoretical ten-asset allocation is then compared
# with the actual fixed P+ age-70 portfolio.
#
# Equality is NOT imposed. The purpose is to test how closely the
# scalar risk-aversion coefficient inferred from the customised
# fund reproduces the internal P+ asset allocation when the model
# is allowed to choose freely across all ten assets.


terminal_theoretical_weight_matrix <- matrix(
  NA_real_,
  nrow = length(
    terminal_profiles
  ),
  ncol = n_assets
)


for (
  i in seq_along(
    terminal_profiles
  )
) {
  
  gamma_i <-
    terminal_gamma_results$
    Gamma_Hat_70[i]
  
  
  terminal_theoretical_weight_matrix[
    i,
  ] <-
    merton_allocation(
      gamma_i
    )
}


rownames(
  terminal_theoretical_weight_matrix
) <- terminal_profiles

colnames(
  terminal_theoretical_weight_matrix
) <- paste0(
  "Asset_",
  seq_len(
    n_assets
  )
)


# ------------------------------------------------------------
# Long-form comparison table
# ------------------------------------------------------------

terminal_asset_comparison <- do.call(
  rbind,
  lapply(
    seq_along(
      terminal_profiles
    ),
    function(i) {
      
      data.frame(
        
        Profile =
          terminal_profiles[i],
        
        Profile_Label =
          unname(
            terminal_profile_labels[
              terminal_profiles[i]
            ]
          ),
        
        Asset =
          seq_len(
            n_assets
          ),
        
        PPlus_Weight =
          as.numeric(
            terminal_pplus_weight_matrix[
              i,
            ]
          ),
        
        Theoretical_Weight =
          as.numeric(
            terminal_theoretical_weight_matrix[
              i,
            ]
          ),
        
        Weight_Difference =
          as.numeric(
            terminal_theoretical_weight_matrix[
              i,
            ] -
              terminal_pplus_weight_matrix[
                i,
              ]
          ),
        
        Absolute_Weight_Difference =
          abs(
            as.numeric(
              terminal_theoretical_weight_matrix[
                i,
              ] -
                terminal_pplus_weight_matrix[
                  i,
                ]
            )
          )
      )
    }
  )
)


print(
  terminal_asset_comparison
)


write.csv(
  terminal_asset_comparison,
  file.path(
    results_dir,
    "terminal_pplus_vs_theoretical_ten_asset_weights.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10E. THESIS-READY TERMINAL TEN-ASSET TABLE
# ============================================================

terminal_asset_table <- data.frame(
  Asset = seq_len(n_assets)
)


for (
  i in seq_along(
    terminal_profiles
  )
) {
  
  profile_name <-
    terminal_profiles[i]
  
  label_short <-
    switch(
      profile_name,
      hoj = "Low",
      mellem = "Moderate",
      lav = "High"
    )
  
  
  terminal_asset_table[[
    paste0(
      "PPlus_",
      label_short
    )
  ]] <-
    as.numeric(
      terminal_pplus_weight_matrix[
        i,
      ]
    )
  
  
  terminal_asset_table[[
    paste0(
      "Theory_",
      label_short
    )
  ]] <-
    as.numeric(
      terminal_theoretical_weight_matrix[
        i,
      ]
    )
}


print(
  terminal_asset_table
)


write.csv(
  terminal_asset_table,
  file.path(
    results_dir,
    "terminal_age70_ten_asset_thesis_table.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10F. TERMINAL PORTFOLIO REPRODUCTION DIAGNOSTICS
# ============================================================

terminal_portfolio_diagnostics <- data.frame(
  
  Profile =
    terminal_profiles,
  
  Profile_Label =
    unname(
      terminal_profile_labels[
        terminal_profiles
      ]
    ),
  
  Gamma_Hat_70 =
    terminal_gamma_results$
    Gamma_Hat_70,
  
  PPlus_Total_Risky_Weight =
    NA_real_,
  
  Theory_Total_Risky_Weight =
    NA_real_,
  
  PPlus_Riskfree_Weight =
    NA_real_,
  
  Theory_Riskfree_Weight =
    NA_real_,
  
  PPlus_Weight_5_10 =
    NA_real_,
  
  Theory_Weight_5_10 =
    NA_real_,
  
  Mean_Absolute_Asset_Gap =
    NA_real_,
  
  Maximum_Absolute_Asset_Gap =
    NA_real_,
  
  Euclidean_Asset_Gap =
    NA_real_
)


for (
  i in seq_along(
    terminal_profiles
  )
) {
  
  theta_i <-
    as.numeric(
      terminal_pplus_weight_matrix[
        i,
      ]
    )
  
  theory_i <-
    as.numeric(
      terminal_theoretical_weight_matrix[
        i,
      ]
    )
  
  difference_i <-
    theory_i -
    theta_i
  
  
  terminal_portfolio_diagnostics$
    PPlus_Total_Risky_Weight[i] <-
    sum(
      theta_i
    )
  
  terminal_portfolio_diagnostics$
    Theory_Total_Risky_Weight[i] <-
    sum(
      theory_i
    )
  
  terminal_portfolio_diagnostics$
    PPlus_Riskfree_Weight[i] <-
    1 -
    sum(
      theta_i
    )
  
  terminal_portfolio_diagnostics$
    Theory_Riskfree_Weight[i] <-
    1 -
    sum(
      theory_i
    )
  
  terminal_portfolio_diagnostics$
    PPlus_Weight_5_10[i] <-
    sum(
      theta_i[
        5:10
      ]
    )
  
  terminal_portfolio_diagnostics$
    Theory_Weight_5_10[i] <-
    sum(
      theory_i[
        5:10
      ]
    )
  
  terminal_portfolio_diagnostics$
    Mean_Absolute_Asset_Gap[i] <-
    mean(
      abs(
        difference_i
      )
    )
  
  terminal_portfolio_diagnostics$
    Maximum_Absolute_Asset_Gap[i] <-
    max(
      abs(
        difference_i
      )
    )
  
  terminal_portfolio_diagnostics$
    Euclidean_Asset_Gap[i] <-
    sqrt(
      sum(
        difference_i^2
      )
    )
}


print(
  terminal_portfolio_diagnostics
)


write.csv(
  terminal_portfolio_diagnostics,
  file.path(
    results_dir,
    "terminal_age70_portfolio_reproduction_diagnostics.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 10G. TERMINAL TEN-ASSET COMPARISON FIGURE
# ============================================================

terminal_asset_plot_data <-
  terminal_asset_comparison %>%
  
  select(
    Profile_Label,
    Asset,
    PPlus_Weight,
    Theoretical_Weight
  ) %>%
  
  pivot_longer(
    cols = c(
      PPlus_Weight,
      Theoretical_Weight
    ),
    names_to = "Strategy",
    values_to = "Weight"
  )


terminal_asset_plot_data$Strategy <- recode(
  terminal_asset_plot_data$Strategy,
  PPlus_Weight = "P+ age-70 portfolio",
  Theoretical_Weight = "Unrestricted theoretical portfolio"
)


terminal_asset_plot_data$Strategy <- factor(
  terminal_asset_plot_data$Strategy,
  levels = c(
    "P+ age-70 portfolio",
    "Unrestricted theoretical portfolio"
  )
)


terminal_asset_comparison_plot <- ggplot(
  terminal_asset_plot_data,
  aes(
    x = factor(
      Asset
    ),
    y = Weight,
    linetype = Strategy,
    group = Strategy
  )
) +
  
  geom_line(
    linewidth = 0.85
  ) +
  
  geom_point(
    size = 1.8
  ) +
  
  facet_wrap(
    ~ Profile_Label
  ) +
  
  labs(
    x = "Asset class",
    y = "Portfolio weight",
    linetype = NULL
  ) +
  
  theme_bw() +
  
  theme(
    legend.position = "bottom"
  )


terminal_asset_comparison_plot


ggsave(
  filename = file.path(
    figure_dir,
    "Terminal PPlus and Theoretical Ten-Asset Allocation.png"
  ),
  plot = terminal_asset_comparison_plot,
  width = 10,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 10H. TERMINAL CALIBRATION INTERNAL CHECK
# ============================================================
#
# For each customised P+ fund, the calibrated gamma_hat_70
# should imply an optimal scalar exposure of exactly one:
#
#   pi^*
#     =
#   (alpha_theta-r) /
#   (gamma_hat_70 sigma_theta^2)
#     =
#   1.


terminal_gamma_results$
  Recovered_Optimal_Fund_Exposure <-
  
  terminal_gamma_results$
  PPlus_Excess_Return /
  
  (
    terminal_gamma_results$
      Gamma_Hat_70 *
      
      terminal_gamma_results$
      PPlus_Variance
  )


terminal_fund_exposure_error <-
  max(
    abs(
      terminal_gamma_results$
        Recovered_Optimal_Fund_Exposure -
        1
    ),
    na.rm = TRUE
  )


cat(
  "\nMaximum terminal customised-fund calibration error:",
  format(
    terminal_fund_exposure_error,
    scientific = TRUE
  ),
  "\n"
)


if (
  terminal_fund_exposure_error >
  1e-10
) {
  
  warning(
    paste0(
      "Terminal customised-fund calibration does not recover ",
      "unit optimal fund exposure."
    )
  )
}

# ============================================================
# 11. CRRA CERTAINTY-EQUIVALENT RETIREMENT WEALTH
# ============================================================

# For gamma != 1:
#
#   CE(X_T) =
#     [E(X_T^(1-gamma))]^(1/(1-gamma)).
#
# For gamma = 1:
#
#   CE(X_T) =
#     exp(E[log(X_T)]).
#
# The calculation is performed in log space for numerical
# stability.

certainty_equivalent_wealth <- function(
    gamma,
    wealth
) {
  
  if (
    !is.finite(gamma) ||
    gamma <= 0 ||
    any(!is.finite(wealth)) ||
    any(wealth <= 0)
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  if (
    abs(
      gamma - 1
    ) < 1e-8
  ) {
    
    return(
      exp(
        mean(
          log(
            wealth
          )
        )
      )
    )
  }
  
  
  log_wealth <-
    log(
      wealth
    )
  
  exponent <-
    (
      1 -
        gamma
    ) *
    log_wealth
  
  m <-
    max(
      exponent
    )
  
  
  log_mean_power <-
    m +
    log(
      mean(
        exp(
          exponent -
            m
        )
      )
    )
  
  
  exp(
    log_mean_power /
      (
        1 -
          gamma
      )
  )
}


# ============================================================
# 12. PRE-RETIREMENT WELFARE COMPARISON
# ============================================================

# Both strategies:
#
#   1. start from the same current pension wealth X_i,
#   2. receive the same deterministic premium stream,
#   3. use the same market parameters,
#   4. use the same calibrated gamma_hat_i,
#   5. are evaluated at the same retirement boundary T.
#
# The only difference is the investment strategy.
#
# P+:
#
#   theta(a,p)
#
# Theoretical optimum:
#
#   w^*(t,x)
#     =
#   (1 + h(t)/x) w^M.
#
# Common random numbers are used for the two strategies to
# reduce Monte Carlo noise in the policy-level comparison.


set.seed(
  2026
)

N_sim <- 10000

dt <- 1 / 12


welfare_results <- data.frame(
  
  policeNr =
    pre_retirement_sample$policeNr,
  
  alder =
    pre_retirement_sample$alder,
  
  investeringsProfil =
    pre_retirement_sample$
    investeringsProfil,
  
  gamma_hat =
    pre_retirement_sample$
    gamma_hat_i,
  
  X_initial =
    pre_retirement_sample$X_i,
  
  h_initial =
    pre_retirement_sample$h_i,
  
  annual_premium =
    pre_retirement_sample$ell_i,
  
  Median_XT_PPlus = NA_real_,
  Median_XT_Optimal = NA_real_,
  Median_Balance_Difference = NA_real_,
  Q10_XT_PPlus = NA_real_,
  Q10_XT_Optimal = NA_real_,
  Q10_Balance_Difference = NA_real_,
  Probability_PPlus_Below_Optimal = NA_real_,
  
  CE_PPlus =
    NA_real_,
  
  CE_Optimal =
    NA_real_,
  
  Additional_Initial_Wealth =
    NA_real_,
  
  Additional_Initial_Wealth_Percent =
    NA_real_
)


for (
  i in seq_len(
    nrow(
      pre_retirement_sample
    )
  )
) {
  
  # ----------------------------------------------------------
  # Member-specific quantities
  # ----------------------------------------------------------
  
  age_i <-
    pre_retirement_sample$alder[i]
  
  X_i <-
    pre_retirement_sample$X_i[i]
  
  h_i <-
    pre_retirement_sample$h_i[i]
  
  ell_i <-
    pre_retirement_sample$ell_i[i]
  
  gamma_hat_i <-
    pre_retirement_sample$gamma_hat_i[i]
  
  profile_i <-
    pre_retirement_sample$
    investeringsProfil[i]
  
  
  # ----------------------------------------------------------
  # Exact simulation grid ending at T
  # ----------------------------------------------------------
  
  age_grid <- seq(
    from = age_i,
    to = retirement_age,
    by = dt
  )
  
  if (
    tail(
      age_grid,
      1
    ) <
    retirement_age
  ) {
    
    age_grid <- c(
      age_grid,
      retirement_age
    )
  }
  
  age_grid[
    length(
      age_grid
    )
  ] <- retirement_age
  
  
  step_dt <- diff(
    age_grid
  )
  
  n_steps <- length(
    step_dt
  )
  
  
  # ----------------------------------------------------------
  # Common Brownian innovations
  # ----------------------------------------------------------
  
  Z_array <- array(
    rnorm(
      N_sim *
        n_assets *
        n_steps
    ),
    dim = c(
      N_sim,
      n_assets,
      n_steps
    )
  )
  
  
  # ==========================================================
  # P+ STRATEGY
  # ==========================================================
  
  # M_t is the growth factor on initial wealth.
  #
  # C_t is the accumulated value of deterministic premiums.
  
  M_pplus <- rep(
    1,
    N_sim
  )
  
  C_pplus <- rep(
    0,
    N_sim
  )
  
  
  for (
    s in seq_len(
      n_steps
    )
  ) {
    
    age_now <-
      age_grid[s]
    
    dt_s <-
      step_dt[s]
    
    
    theta_t <- as.numeric(
      pf(
        age_now,
        profile_i
      )
    )
    
    
    portfolio_drift <-
      r +
      sum(
        theta_t *
          excess_return
      )
    
    
    portfolio_diffusion <-
      as.numeric(
        t(theta_t) %*%
          L
      )
    
    
    portfolio_variance <-
      sum(
        portfolio_diffusion^2
      )
    
    
    Z_s <-
      Z_array[
        ,
        ,
        s,
        drop = FALSE
      ]
    
    dim(
      Z_s
    ) <- c(
      N_sim,
      n_assets
    )
    
    
    shock_pplus <-
      as.numeric(
        Z_s %*%
          portfolio_diffusion
      ) *
      sqrt(
        dt_s
      )
    
    
    log_return_pplus <-
      (
        portfolio_drift -
          0.5 *
          portfolio_variance
      ) *
      dt_s +
      shock_pplus
    
    
    growth_factor <-
      exp(
        log_return_pplus
      )
    
    
    premium_increment <-
      ell_i *
      dt_s
    
    
    M_pplus <-
      M_pplus *
      growth_factor
    
    
    C_pplus <-
      (
        C_pplus +
          0.5 *
          premium_increment
      ) *
      growth_factor +
      0.5 *
      premium_increment
  }
  
  
  X_PPlus_T <- function(
    initial_wealth
  ) {
    
    initial_wealth *
      M_pplus +
      C_pplus
  }
  
  
  # ==========================================================
  # THEORETICAL OPTIMAL STRATEGY
  # ==========================================================
  
  # Define total effective resources
  #
  #   Y(t) = X(t) + h(t).
  #
  # Under the optimal strategy:
  #
  #   w^M =
  #     (1/gamma_hat)
  #     Sigma^{-1}(alpha-r1).
  #
  # The diffusion exposure of Y is therefore generated by w^M.
  
  w_M_i <- merton_allocation(
    gamma_hat_i
  )
  
  
  merton_diffusion <-
    as.numeric(
      t(w_M_i) %*%
        L
    )
  
  
  merton_variance <-
    as.numeric(
      t(w_M_i) %*%
        Sigma %*%
        w_M_i
    )
  
  
  merton_expected_excess_return <-
    sum(
      w_M_i *
        excess_return
    )
  
  
  Y_optimal <-
    rep(
      X_i +
        h_i,
      N_sim
    )
  
  
  for (
    s in seq_len(
      n_steps
    )
  ) {
    
    dt_s <-
      step_dt[s]
    
    
    Z_s <-
      Z_array[
        ,
        ,
        s,
        drop = FALSE
      ]
    
    dim(
      Z_s
    ) <- c(
      N_sim,
      n_assets
    )
    
    
    shock_optimal <-
      as.numeric(
        Z_s %*%
          merton_diffusion
      ) *
      sqrt(
        dt_s
      )
    
    
    log_return_optimal <-
      (
        r +
          merton_expected_excess_return -
          0.5 *
          merton_variance
      ) *
      dt_s +
      shock_optimal
    
    
    Y_optimal <-
      Y_optimal *
      exp(
        log_return_optimal
      )
  }
  
  
  # At retirement h(T) = 0, hence Y(T) = X(T).
  
  X_Optimal_T <-
    Y_optimal
  
  X_PPlus_T_current <- X_PPlus_T(X_i)
  
  # ==========================================================
  # RETIREMENT-BALANCE DISTRIBUTION METRICS
  # ==========================================================
  # Median terminal balance, 10th percentile terminal balance,
  # and the matched-path probability that P+ finishes below the
  # theoretical benchmark.
  
  median_XT_PPlus <- median(X_PPlus_T_current, na.rm = TRUE)
  median_XT_Optimal <- median(X_Optimal_T, na.rm = TRUE)
  q10_XT_PPlus <- as.numeric(quantile(X_PPlus_T_current, 0.10, na.rm = TRUE))
  q10_XT_Optimal <- as.numeric(quantile(X_Optimal_T, 0.10, na.rm = TRUE))
  probability_PPlus_below_Optimal <- mean(
    X_PPlus_T_current < X_Optimal_T,
    na.rm = TRUE
  )
  
  # ==========================================================
  # CERTAINTY-EQUIVALENT RETIREMENT WEALTH
  # ==========================================================
  
  CE_PPlus_i <-
    certainty_equivalent_wealth(
      gamma_hat_i,
      X_PPlus_T_current
    )
  
  
  CE_Optimal_i <-
    certainty_equivalent_wealth(
      gamma_hat_i,
      X_Optimal_T
    )
  
  
  # ==========================================================
  # ADDITIONAL INITIAL WEALTH REQUIRED UNDER P+
  # ==========================================================
  #
  # Delta_i solves
  #
  #   CE_i^P+(X_i + Delta_i)
  #     =
  #   CE_i^*(X_i).
  
  welfare_difference <- function(
    Delta
  ) {
    
    certainty_equivalent_wealth(
      gamma_hat_i,
      X_PPlus_T(
        X_i +
          Delta
      )
    ) -
      CE_Optimal_i
  }
  
  
  lower_bound <- 0
  
  upper_bound <-
    50 *
    max(
      X_i,
      1
    )
  
  
  root <- tryCatch(
    
    uniroot(
      welfare_difference,
      lower = lower_bound,
      upper = upper_bound,
      extendInt = "yes",
      tol = 1
    ),
    
    error = function(e) {
      NULL
    }
  )
  
  
  Delta_i <-
    if (
      is.null(
        root
      )
    ) {
      
      NA_real_
      
    } else {
      
      root$root
    }
  
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  welfare_results$Median_XT_PPlus[i] <- median_XT_PPlus
  welfare_results$Median_XT_Optimal[i] <- median_XT_Optimal
  welfare_results$Median_Balance_Difference[i] <-
    median_XT_PPlus - median_XT_Optimal
  welfare_results$Q10_XT_PPlus[i] <- q10_XT_PPlus
  welfare_results$Q10_XT_Optimal[i] <- q10_XT_Optimal
  welfare_results$Q10_Balance_Difference[i] <-
    q10_XT_PPlus - q10_XT_Optimal
  welfare_results$Probability_PPlus_Below_Optimal[i] <-
    probability_PPlus_below_Optimal
  
  welfare_results$
    CE_PPlus[i] <-
    CE_PPlus_i
  
  welfare_results$
    CE_Optimal[i] <-
    CE_Optimal_i
  
  welfare_results$
    Additional_Initial_Wealth[i] <-
    Delta_i
  
  welfare_results$
    Additional_Initial_Wealth_Percent[i] <-
    100 *
    Delta_i /
    X_i
  
  
  cat(
    sprintf(
      paste0(
        "[%d/%d] policeNr=%s ",
        "age=%.2f gamma_hat=%.3f ",
        "CE_PPlus=%.0f ",
        "CE_Optimal=%.0f ",
        "additional_initial_wealth=%.2f%%\n"
      ),
      i,
      nrow(
        pre_retirement_sample
      ),
      welfare_results$policeNr[i],
      age_i,
      gamma_hat_i,
      CE_PPlus_i,
      CE_Optimal_i,
      welfare_results$
        Additional_Initial_Wealth_Percent[i]
    )
  )
}

# ============================================================
# PRE-RETIREMENT WELFARE VALIDATION CHECKS
# ============================================================

if (
  any(
    !is.finite(
      welfare_results$CE_PPlus
    )
  ) ||
  any(
    !is.finite(
      welfare_results$CE_Optimal
    )
  ) ||
  any(
    welfare_results$CE_PPlus <= 0
  ) ||
  any(
    welfare_results$CE_Optimal <= 0
  )
) {
  
  stop(
    "Non-finite or non-positive certainty-equivalent wealth detected."
  )
}


if (
  any(
    !is.finite(
      welfare_results$Additional_Initial_Wealth
    )
  )
) {
  
  warning(
    "Additional initial wealth could not be calculated for at least one member."
  )
}


monte_carlo_tolerance <- 0.01

welfare_dominance_failures <-
  which(
    welfare_results$CE_Optimal <
      (
        1 -
          monte_carlo_tolerance
      ) *
      welfare_results$CE_PPlus
  )

if (
  length(
    welfare_dominance_failures
  ) > 0
) {
  
  warning(
    paste(
      length(
        welfare_dominance_failures
      ),
      "members have an optimal certainty equivalent more than",
      100 * monte_carlo_tolerance,
      "percent below the P+ certainty equivalent."
    )
  )
}

# Khemka-style terminology: the same quantity is the extra starting
# balance required under P+ to match the theoretical benchmark.
welfare_results$Extra_Starting_Balance <-
  welfare_results$Additional_Initial_Wealth
welfare_results$Extra_Starting_Balance_Percent <-
  welfare_results$Additional_Initial_Wealth_Percent


# ============================================================
# 12A. RETIREMENT-BALANCE DISTRIBUTION SUMMARY
# ============================================================

retirement_balance_summary <- welfare_results %>%
  summarise(
    N = n(),
    Mean_Median_XT_PPlus = mean(Median_XT_PPlus, na.rm = TRUE),
    Mean_Median_XT_Optimal = mean(Median_XT_Optimal, na.rm = TRUE),
    Mean_Median_Balance_Difference = mean(Median_Balance_Difference, na.rm = TRUE),
    Mean_Q10_XT_PPlus = mean(Q10_XT_PPlus, na.rm = TRUE),
    Mean_Q10_XT_Optimal = mean(Q10_XT_Optimal, na.rm = TRUE),
    Mean_Q10_Balance_Difference = mean(Q10_Balance_Difference, na.rm = TRUE),
    Mean_Probability_PPlus_Below_Optimal = mean(
      Probability_PPlus_Below_Optimal,
      na.rm = TRUE
    ),
    Median_Probability_PPlus_Below_Optimal = median(
      Probability_PPlus_Below_Optimal,
      na.rm = TRUE
    )
  )

retirement_balance_summary_by_profile <- welfare_results %>%
  group_by(investeringsProfil) %>%
  summarise(
    N = n(),
    Mean_Median_XT_PPlus = mean(Median_XT_PPlus, na.rm = TRUE),
    Mean_Median_XT_Optimal = mean(Median_XT_Optimal, na.rm = TRUE),
    Mean_Q10_XT_PPlus = mean(Q10_XT_PPlus, na.rm = TRUE),
    Mean_Q10_XT_Optimal = mean(Q10_XT_Optimal, na.rm = TRUE),
    Mean_Probability_PPlus_Below_Optimal = mean(
      Probability_PPlus_Below_Optimal,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(retirement_balance_summary)
print(retirement_balance_summary_by_profile)

write.csv(
  retirement_balance_summary,
  file.path(results_dir, "theory_practice_retirement_balance_summary.csv"),
  row.names = FALSE
)

write.csv(
  retirement_balance_summary_by_profile,
  file.path(results_dir, "theory_practice_retirement_balance_summary_by_profile.csv"),
  row.names = FALSE
)


# ============================================================
# 13. PRE-RETIREMENT WELFARE SUMMARY
# ============================================================

welfare_summary_by_profile <-
  welfare_results %>%
  
  group_by(
    investeringsProfil
  ) %>%
  
  summarise(
    
    N = n(),
    
    Mean_Gamma_Hat =
      mean(
        gamma_hat,
        na.rm = TRUE
      ),
    
    Mean_CE_PPlus =
      mean(
        CE_PPlus,
        na.rm = TRUE
      ),
    
    Mean_CE_Optimal =
      mean(
        CE_Optimal,
        na.rm = TRUE
      ),
    
    Mean_CE_Ratio =
      mean(
        CE_Optimal /
          CE_PPlus,
        na.rm = TRUE
      ),
    
    Median_Additional_Initial_Wealth_Percent =
      median(
        Additional_Initial_Wealth_Percent,
        na.rm = TRUE
      ),
    
    Mean_Additional_Initial_Wealth_Percent =
      mean(
        Additional_Initial_Wealth_Percent,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


print(
  welfare_summary_by_profile
)


write.csv(
  welfare_results,
  file.path(
    results_dir,
    "theory_practice_pre_retirement_welfare.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    welfare_summary_by_profile
  ),
  file.path(
    results_dir,
    "theory_practice_pre_retirement_welfare_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. PRE-RETIREMENT WELFARE FIGURE
# ============================================================

welfare_plot_data <-
  welfare_results %>%
  
  mutate(
    
    Profile_Label = recode(
      investeringsProfil,
      hoj = "Low",
      mellem = "Moderate",
      lav = "High"
    )
  )


welfare_plot_data$Profile_Label <- factor(
  welfare_plot_data$Profile_Label,
  levels = c(
    "Low",
    "Moderate",
    "High"
  )
)


additional_wealth_plot <- ggplot(
  welfare_plot_data,
  aes(
    x = Profile_Label,
    y = Additional_Initial_Wealth_Percent
  )
) +
  
  geom_boxplot(
    alpha = 0.6
  ) +
  
  labs(
    x = "Risk aversion profile",
    y = "Additional initial wealth required under P+ (%)"
  ) +
  
  theme_bw()


additional_wealth_plot


ggsave(
  filename = file.path(
    figure_dir,
    "Theory vs PPlus Additional Initial Wealth Required.png"
  ),
  plot = additional_wealth_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 15. POST-RETIREMENT THEORETICAL BENCHMARK
# ============================================================
#
# After retirement:
#
#   h(t) = 0
#
# and the optimal risky-asset allocation is
#
#   w^M =
#     (1/gamma)
#     Sigma^{-1}(alpha-r1).
#
# The optimal consumption feedback rule is
#
#   c^*(t,x) = x / abar(t).
#
# The numerical implementation uses terminal age N and
# subjective discount rate rho.


terminal_age <- 100

rho <- 0


# ------------------------------------------------------------
# Effective rate
# ------------------------------------------------------------
#
# nu =
#
#   [rho -
#    (1-gamma)
#    (r + Theta^2/(2 gamma))]
#   / gamma

nu_function <- function(
    gamma,
    rho_value = rho
) {
  
  (
    rho_value -
      (
        1 -
          gamma
      ) *
      (
        r +
          Theta_sq /
          (
            2 *
              gamma
          )
      )
  ) /
    gamma
}


# ------------------------------------------------------------
# Effective annuity factor
# ------------------------------------------------------------
#
# abar(t) =
#   integral_t^N exp[-nu(s-t)] ds
#
# If nu != 0:
#
#   abar(t) =
#     [1-exp(-nu(N-t))]/nu.

annuity_factor <- function(
    t,
    gamma,
    N = terminal_age,
    rho_value = rho
) {
  
  remaining_horizon <-
    N -
    t
  
  
  if (
    remaining_horizon <= 0
  ) {
    return(
      0
    )
  }
  
  
  nu <-
    nu_function(
      gamma = gamma,
      rho_value = rho_value
    )
  
  
  if (
    abs(
      nu
    ) < 1e-8
  ) {
    
    return(
      remaining_horizon
    )
  }
  
  
  (
    1 -
      exp(
        -nu *
          remaining_horizon
      )
  ) /
    nu
}


# ------------------------------------------------------------
# Post-retirement value function
# ------------------------------------------------------------

post_retirement_value <- function(
    t,
    x,
    gamma,
    N = terminal_age,
    rho_value = rho
) {
  
  abar_t <-
    annuity_factor(
      t = t,
      gamma = gamma,
      N = N,
      rho_value = rho_value
    )
  
  
  exp(
    -rho_value *
      t
  ) *
    abar_t^gamma *
    x^(
      1 -
        gamma
    ) /
    (
      1 -
        gamma
    )
}


# ============================================================
# 16. POST-RETIREMENT SAMPLE
# ============================================================

post_retirement_sample <-
  comparison_sample %>%
  
  filter(
    alder >= retirement_age,
    alder < terminal_age
  )


cat(
  "\nPost-retirement comparison sample (age >= T):",
  nrow(
    post_retirement_sample
  ),
  "\n"
)


# ============================================================
# 17. P+ INVESTMENT WITH THEORETICAL CONSUMPTION RULE
# ============================================================
#
# This is not interpreted as the actual P+ consumption policy.
#
# The P+ investment allocation theta(t) is combined with the
# same theoretical consumption feedback rule
#
#   c^*(t,x) = x / abar(t)
#
# used in the theoretical benchmark.
#
# Holding the consumption rule fixed isolates the effect of the
# investment allocation.


expected_utility_pplus_investment <- function(
    t,
    gamma,
    profile,
    N = terminal_age,
    rho_value = rho,
    dt = 1 / 12
) {
  
  if (
    t >= N
  ) {
    return(
      NA_real_
    )
  }
  
  
  time_grid <- seq(
    from = t,
    to = N,
    by = dt
  )
  
  
  if (
    tail(
      time_grid,
      1
    ) <
    N
  ) {
    
    time_grid <- c(
      time_grid,
      N
    )
  }
  
  
  time_grid[
    length(
      time_grid
    )
  ] <- N
  
  
  step_dt <- diff(
    time_grid
  )
  
  
  n_steps <- length(
    step_dt
  )
  
  
  portfolio_mean <- numeric(
    n_steps
  )
  
  portfolio_variance <- numeric(
    n_steps
  )
  
  
  for (
    s in seq_len(
      n_steps
    )
  ) {
    
    time_mid <-
      0.5 *
      (
        time_grid[s] +
          time_grid[s + 1]
      )
    
    
    theta_t <- as.numeric(
      pf(
        time_mid,
        profile
      )
    )
    
    
    portfolio_mean[s] <-
      r +
      sum(
        theta_t *
          excess_return
      )
    
    
    portfolio_variance[s] <-
      as.numeric(
        t(theta_t) %*%
          Sigma %*%
          theta_t
      )
  }
  
  
  nu <-
    nu_function(
      gamma = gamma,
      rho_value = rho_value
    )
  
  
  abar_initial <-
    annuity_factor(
      t = t,
      gamma = gamma,
      N = N,
      rho_value = rho_value
    )
  
  
  # ----------------------------------------------------------
  # Dynamics of Y_u = X_u / abar(u)
  # ----------------------------------------------------------
  
  log_mean_increment <-
    (
      portfolio_mean -
        nu -
        0.5 *
        portfolio_variance
    ) *
    step_dt
  
  
  log_variance_increment <-
    portfolio_variance *
    step_dt
  
  
  cumulative_log_mean <-
    c(
      0,
      cumsum(
        log_mean_increment
      )
    )
  
  
  cumulative_log_variance <-
    c(
      0,
      cumsum(
        log_variance_increment
      )
    )
  
  
  log_Y_initial <-
    -log(
      abar_initial
    )
  
  
  mean_log_Y <-
    log_Y_initial +
    cumulative_log_mean
  
  
  variance_log_Y <-
    cumulative_log_variance
  
  
  expected_utility_integrand <-
    exp(
      (
        1 -
          gamma
      ) *
        mean_log_Y +
        0.5 *
        (
          1 -
            gamma
        )^2 *
        variance_log_Y
    ) /
    (
      1 -
        gamma
    )
  
  
  expected_utility_integrand <-
    exp(
      -rho_value *
        time_grid
    ) *
    expected_utility_integrand
  
  
  sum(
    0.5 *
      (
        expected_utility_integrand[-1] +
          expected_utility_integrand[
            -length(
              expected_utility_integrand
            )
          ]
      ) *
      step_dt
  )
}


# ============================================================
# 18. POST-RETIREMENT WELFARE COMPARISON
# ============================================================

post_retirement_results <- data.frame(
  
  policeNr =
    post_retirement_sample$policeNr,
  
  alder =
    post_retirement_sample$alder,
  
  investeringsProfil =
    post_retirement_sample$
    investeringsProfil,
  
  gamma_hat =
    post_retirement_sample$
    gamma_hat_i,
  
  X_initial =
    post_retirement_sample$X_i,
  
  V_Optimal_Unit_Wealth =
    NA_real_,
  
  J_PPlus_Investment_Unit_Wealth =
    NA_real_,
  
  Additional_Initial_Wealth_Percent =
    NA_real_
)


if (
  nrow(
    post_retirement_sample
  ) > 0
) {
  
  for (
    i in seq_len(
      nrow(
        post_retirement_sample
      )
    )
  ) {
    
    age_i <-
      post_retirement_sample$alder[i]
    
    gamma_hat_i <-
      post_retirement_sample$gamma_hat_i[i]
    
    profile_i <-
      post_retirement_sample$
      investeringsProfil[i]
    
    
    # --------------------------------------------------------
    # Theoretical value from unit wealth
    # --------------------------------------------------------
    
    V_optimal_i <-
      post_retirement_value(
        t = age_i,
        x = 1,
        gamma = gamma_hat_i,
        N = terminal_age,
        rho_value = rho
      )
    
    
    # --------------------------------------------------------
    # P+ investment combined with the same theoretical
    # consumption feedback rule
    # --------------------------------------------------------
    
    J_pplus_i <-
      expected_utility_pplus_investment(
        t = age_i,
        gamma = gamma_hat_i,
        profile = profile_i,
        N = terminal_age,
        rho_value = rho
      )
    
    
    # --------------------------------------------------------
    # Additional initial wealth required under P+
    # --------------------------------------------------------
    #
    # CRRA homogeneity gives
    #
    #   J(t,kx) = k^(1-gamma) J(t,x).
    
    additional_initial_wealth_pct <-
      100 *
      (
        (
          V_optimal_i /
            J_pplus_i
        )^(
          1 /
            (
              1 -
                gamma_hat_i
            )
        ) -
          1
      )
    
    
    post_retirement_results$
      V_Optimal_Unit_Wealth[i] <-
      V_optimal_i
    
    post_retirement_results$
      J_PPlus_Investment_Unit_Wealth[i] <-
      J_pplus_i
    
    post_retirement_results$
      Additional_Initial_Wealth_Percent[i] <-
      additional_initial_wealth_pct
    
    
    cat(
      sprintf(
        paste0(
          "[%d/%d] policeNr=%s ",
          "age=%.2f gamma_hat=%.3f ",
          "additional_initial_wealth=%.2f%%\n"
        ),
        i,
        nrow(
          post_retirement_sample
        ),
        post_retirement_results$policeNr[i],
        age_i,
        gamma_hat_i,
        additional_initial_wealth_pct
      )
    )
  }
}


print(
  post_retirement_results
)


write.csv(
  post_retirement_results,
  file.path(
    results_dir,
    "theory_practice_post_retirement_welfare.csv"
  ),
  row.names = FALSE
)


if (
  nrow(
    post_retirement_results
  ) > 0
) {
  
  cat(
    "\nMean additional initial wealth required post-retirement:",
    round(
      mean(
        post_retirement_results$
          Additional_Initial_Wealth_Percent,
        na.rm = TRUE
      ),
      2
    ),
    "%\n"
  )
  
  
  cat(
    "Median additional initial wealth required post-retirement:",
    round(
      median(
        post_retirement_results$
          Additional_Initial_Wealth_Percent,
        na.rm = TRUE
      ),
      2
    ),
    "%\n"
  )
}


# ============================================================
# 19. REPRESENTATIVE POST-RETIREMENT MEMBER
# ============================================================
#
# The following section illustrates the Chapter 4 optimal
# consumption rule
#
#   c^*(t,x) = x / abar(t).
#
# This is a theoretical illustration. It is separate from the
# policy-level welfare comparison above.
#
# The representative member is the post-retirement observation
# whose calibrated gamma_hat is closest to the median of the
# post-retirement subsample.


if (
  nrow(
    post_retirement_sample
  ) > 0
) {
  
  median_gamma_post <-
    median(
      post_retirement_sample$
        gamma_hat_i,
      na.rm = TRUE
    )
  
  
  representative_index <-
    which.min(
      abs(
        post_retirement_sample$
          gamma_hat_i -
          median_gamma_post
      )
    )
  
  
  representative_member <-
    post_retirement_sample[
      representative_index,
      ,
      drop = FALSE
    ]
  
  
  t0 <-
    as.numeric(
      representative_member$
        alder[1]
    )
  
  X0 <-
    as.numeric(
      representative_member$
        X_i[1]
    )
  
  gamma_hat_rep <-
    as.numeric(
      representative_member$
        gamma_hat_i[1]
    )
  
  profile_rep <-
    as.character(
      representative_member$
        investeringsProfil[1]
    )
  
  police_rep <-
    representative_member$
    policeNr[1]
  
  
  cat(
    "\n==================================================\n"
  )
  
  cat(
    "REPRESENTATIVE POST-RETIREMENT MEMBER\n"
  )
  
  cat(
    "==================================================\n"
  )
  
  cat(
    "Police number:",
    police_rep,
    "\n"
  )
  
  cat(
    "Current age:",
    round(
      t0,
      2
    ),
    "\n"
  )
  
  cat(
    "Current pension wealth:",
    format(
      round(
        X0,
        0
      ),
      big.mark = ",",
      scientific = FALSE
    ),
    "\n"
  )
  
  cat(
    "Calibrated gamma_hat:",
    round(
      gamma_hat_rep,
      4
    ),
    "\n"
  )
  
  cat(
    "Risk aversion profile:",
    profile_rep,
    "\n"
  )
  
  
  # ==========================================================
  # 20. REPRESENTATIVE MEMBER PARAMETERS
  # ==========================================================
  
  nu_rep <-
    nu_function(
      gamma = gamma_hat_rep,
      rho_value = rho
    )
  
  
  abar0 <-
    annuity_factor(
      t = t0,
      gamma = gamma_hat_rep,
      N = terminal_age,
      rho_value = rho
    )
  
  
  if (
    !is.finite(
      abar0
    ) ||
    abar0 <= 0
  ) {
    
    stop(
      "The initial effective annuity factor is non-positive or undefined."
    )
  }
  
  
  initial_optimal_consumption <-
    X0 /
    abar0
  
  
  cat(
    "Effective rate nu:",
    round(
      nu_rep,
      6
    ),
    "\n"
  )
  
  cat(
    "Initial effective annuity factor abar(t0):",
    round(
      abar0,
      4
    ),
    "\n"
  )
  
  cat(
    "Initial optimal annual consumption:",
    format(
      round(
        initial_optimal_consumption,
        0
      ),
      big.mark = ",",
      scientific = FALSE
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Merton allocation
  # ----------------------------------------------------------
  
  w_M_rep <-
    merton_allocation(
      gamma_hat_rep
    )
  
  
  merton_excess_return_rep <-
    sum(
      w_M_rep *
        excess_return
    )
  
  
  merton_variance_rep <-
    as.numeric(
      t(w_M_rep) %*%
        Sigma %*%
        w_M_rep
    )
  
  
  cat(
    "Merton expected excess return:",
    round(
      merton_excess_return_rep,
      6
    ),
    "\n"
  )
  
  cat(
    "Merton portfolio volatility:",
    round(
      sqrt(
        merton_variance_rep
      ),
      6
    ),
    "\n"
  )
  
  
  # ==========================================================
  # 21. OPTIMAL POST-RETIREMENT CONSUMPTION PATH
  # ==========================================================
  
  set.seed(
    2026
  )
  
  
  dt_consumption <-
    1 / 12
  
  
  age_grid <-
    seq(
      from = t0,
      to = terminal_age,
      by = dt_consumption
    )
  
  
  if (
    tail(
      age_grid,
      1
    ) <
    terminal_age
  ) {
    
    age_grid <- c(
      age_grid,
      terminal_age
    )
  }
  
  
  age_grid[
    length(
      age_grid
    )
  ] <- terminal_age
  
  
  n_grid <-
    length(
      age_grid
    )
  
  
  elapsed_time <-
    age_grid -
    t0
  
  
  abar_grid <-
    sapply(
      age_grid,
      function(age) {
        
        annuity_factor(
          t = age,
          gamma = gamma_hat_rep,
          N = terminal_age,
          rho_value = rho
        )
      }
    )
  
  
  # ----------------------------------------------------------
  # Closed-loop stochastic factor
  # ----------------------------------------------------------
  
  merton_variance_rep <-
    Theta_sq /
    gamma_hat_rep^2
  
  
  merton_volatility_rep <-
    sqrt(
      merton_variance_rep
    )
  
  
  stochastic_factor_drift <-
    r +
    Theta_sq /
    gamma_hat_rep -
    Theta_sq /
    (
      2 *
        gamma_hat_rep^2
    ) -
    nu_rep
  
  
  stochastic_factor <-
    rep(
      1,
      n_grid
    )
  
  
  if (
    n_grid > 1
  ) {
    
    for (
      j in 2:n_grid
    ) {
      
      dt_j <-
        age_grid[j] -
        age_grid[j - 1]
      
      Z_j <-
        rnorm(
          1
        )
      
      
      stochastic_factor[j] <-
        stochastic_factor[j - 1] *
        exp(
          stochastic_factor_drift *
            dt_j +
            merton_volatility_rep *
            sqrt(
              dt_j
            ) *
            Z_j
        )
    }
  }
  
  
  # ----------------------------------------------------------
  # Optimal wealth
  # ----------------------------------------------------------
  
  optimal_wealth <-
    X0 *
    (
      abar_grid /
        abar0
    ) *
    stochastic_factor
  
  
  optimal_wealth[
    length(
      optimal_wealth
    )
  ] <- 0
  
  
  # ----------------------------------------------------------
  # Realised optimal consumption
  # ----------------------------------------------------------
  
  realised_optimal_consumption <-
    (
      X0 /
        abar0
    ) *
    stochastic_factor
  
  
  # ----------------------------------------------------------
  # Expected optimal consumption
  # ----------------------------------------------------------
  
  expected_stochastic_factor <-
    exp(
      (
        r +
          Theta_sq /
          gamma_hat_rep -
          nu_rep
      ) *
        elapsed_time
    )
  
  
  expected_optimal_consumption <-
    (
      X0 /
        abar0
    ) *
    expected_stochastic_factor
  
  
  # ----------------------------------------------------------
  # Optimal consumption-to-wealth ratio
  # ----------------------------------------------------------
  
  consumption_wealth_ratio <-
    rep(
      NA_real_,
      n_grid
    )
  
  
  non_terminal <-
    age_grid <
    terminal_age
  
  
  consumption_wealth_ratio[
    non_terminal
  ] <-
    1 /
    abar_grid[
      non_terminal
    ]
  
  
  # ----------------------------------------------------------
  # Store representative-member paths
  # ----------------------------------------------------------
  
  consumption_path_data <- data.frame(
    
    Age =
      age_grid,
    
    Effective_Annuity_Factor =
      abar_grid,
    
    Optimal_Wealth =
      optimal_wealth,
    
    Realised_Optimal_Consumption =
      realised_optimal_consumption,
    
    Expected_Optimal_Consumption =
      expected_optimal_consumption,
    
    Consumption_Wealth_Ratio =
      consumption_wealth_ratio
  )
  
  
  write.csv(
    consumption_path_data,
    file.path(
      results_dir,
      "post_retirement_optimal_consumption_representative_member.csv"
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 22. OPTIMAL CONSUMPTION FIGURE
  # ==========================================================
  
  consumption_plot_base <-
    subset(
      consumption_path_data,
      Age <
        terminal_age
    )
  
  
  consumption_plot_data <- rbind(
    
    data.frame(
      Age =
        consumption_plot_base$Age,
      
      Consumption =
        consumption_plot_base$
        Realised_Optimal_Consumption,
      
      Series =
        "One simulated optimal path"
    ),
    
    data.frame(
      Age =
        consumption_plot_base$Age,
      
      Consumption =
        consumption_plot_base$
        Expected_Optimal_Consumption,
      
      Series =
        "Expected optimal consumption"
    )
  )
  
  
  consumption_plot_data$Series <- factor(
    consumption_plot_data$Series,
    levels = c(
      "One simulated optimal path",
      "Expected optimal consumption"
    )
  )
  
  
  optimal_consumption_plot <- ggplot(
    consumption_plot_data,
    aes(
      x = Age,
      y = Consumption,
      linetype = Series,
      group = Series
    )
  ) +
    
    geom_line(
      linewidth = 0.85
    ) +
    
    scale_linetype_manual(
      name = NULL,
      values = c(
        "One simulated optimal path" = "solid",
        "Expected optimal consumption" = "dashed"
      )
    ) +
    
    scale_y_continuous(
      labels = scales::comma,
      breaks = scales::pretty_breaks(
        n = 6
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::pretty_breaks(
        n = 8
      )
    ) +
    
    labs(
      x = "Age",
      y = "Optimal annual consumption (DKK)"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  
  optimal_consumption_plot
  
  
  ggsave(
    filename = file.path(
      figure_dir,
      "Optimal Post-Retirement Consumption - Representative Member.png"
    ),
    plot = optimal_consumption_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  
  # ==========================================================
  # 23. CONSUMPTION-TO-WEALTH RATIO
  # ==========================================================
  
  # Since
  #
  #   c^*(t,x)/x = 1/abar(t),
  #
  # the ratio is deterministic for a given gamma and rho.
  #
  # The final year is excluded because abar(t) tends to zero
  # as t approaches N.
  
  ratio_plot_end_age <-
    terminal_age -
    1
  
  
  consumption_ratio_plot_data <-
    subset(
      consumption_path_data,
      Age <=
        ratio_plot_end_age &
        is.finite(
          Consumption_Wealth_Ratio
        )
    )
  
  
  consumption_ratio_plot_data$
    Consumption_Wealth_Percent <-
    100 *
    consumption_ratio_plot_data$
    Consumption_Wealth_Ratio
  
  
  optimal_consumption_ratio_plot <- ggplot(
    consumption_ratio_plot_data,
    aes(
      x = Age,
      y = Consumption_Wealth_Percent
    )
  ) +
    
    geom_line(
      linewidth = 0.85
    ) +
    
    scale_x_continuous(
      breaks = scales::pretty_breaks(
        n = 8
      )
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks(
        n = 6
      ),
      labels = scales::label_number(
        accuracy = 1
      )
    ) +
    
    labs(
      x = "Age",
      y = "Optimal annual consumption-to-wealth ratio (%)"
    ) +
    
    theme_bw()
  
  
  optimal_consumption_ratio_plot
  
  
  ggsave(
    filename = file.path(
      figure_dir,
      "Optimal Post-Retirement Consumption-to-Wealth Ratio.png"
    ),
    plot = optimal_consumption_ratio_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  
  # ==========================================================
  # 24. SENSITIVITY TO SUBJECTIVE DISCOUNT RATE RHO
  # ==========================================================
  
  # gamma_hat is held fixed.
  #
  # Only rho is varied.
  #
  # For each rho:
  #
  #   c^*(t,x)/x =
  #     1/abar_rho(t).
  
  rho_values <- c(
    0.00,
    0.05,
    0.10,
    0.50
  )
  
  
  rho_ratio_age_grid <-
    age_grid[
      age_grid <=
        ratio_plot_end_age
    ]
  
  
  rho_ratio_data <- do.call(
    rbind,
    
    lapply(
      rho_values,
      
      function(
    rho_value
      ) {
        
        abar_values <- sapply(
          rho_ratio_age_grid,
          
          function(
    age
          ) {
            
            annuity_factor(
              t = age,
              gamma = gamma_hat_rep,
              N = terminal_age,
              rho_value = rho_value
            )
          }
        )
        
        
        data.frame(
          Age =
            rho_ratio_age_grid,
          
          Rho =
            rho_value,
          
          Consumption_Wealth_Percent =
            100 /
            abar_values
        )
      }
    )
  )
  
  
  rho_ratio_data$Rho_Label <- factor(
    rho_ratio_data$Rho,
    levels = rho_values,
    labels = c(
      "0",
      "0.05",
      "0.10",
      "0.50"
    )
  )
  
  
  optimal_consumption_ratio_rho_plot <- ggplot(
    rho_ratio_data,
    aes(
      x = Age,
      y = Consumption_Wealth_Percent,
      linetype = Rho_Label,
      group = Rho_Label
    )
  ) +
    
    geom_line(
      linewidth = 0.9
    ) +
    
    scale_linetype_manual(
      name = expression(rho),
      values = c(
        "solid",
        "dashed",
        "dotted",
        "dotdash"
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::pretty_breaks(
        n = 8
      )
    ) +
    
    scale_y_continuous(
      breaks = scales::pretty_breaks(
        n = 6
      ),
      labels = scales::label_number(
        accuracy = 1
      )
    ) +
    
    labs(
      x = "Age",
      y = "Optimal annual consumption-to-wealth ratio (%)"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "bottom"
    )
  
  
  optimal_consumption_ratio_rho_plot
  
  
  ggsave(
    filename = file.path(
      figure_dir,
      "Optimal Post-Retirement Consumption-to-Wealth Ratio - Rho Sensitivity.png"
    ),
    plot = optimal_consumption_ratio_rho_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  
  write.csv(
    rho_ratio_data,
    file.path(
      results_dir,
      "post_retirement_consumption_ratio_rho_sensitivity.csv"
    ),
    row.names = FALSE
  )
  
  # ==========================================================
  # 24A. RHO SENSITIVITY AT SELECTED RETIREMENT AGES
  # ==========================================================
  #
  # This table makes the economic interpretation of rho explicit.
  # It reports the optimal annual consumption-to-wealth ratio at
  # selected ages for each value of rho.
  
  
  selected_consumption_ages <- c(
    67,
    70,
    75,
    80,
    85,
    90,
    95,
    99
  )
  
  
  rho_selected_age_summary <- do.call(
    rbind,
    lapply(
      rho_values,
      function(rho_value) {
        
        abar_values <-
          sapply(
            selected_consumption_ages,
            function(age) {
              
              annuity_factor(
                t = age,
                gamma = gamma_hat_rep,
                N = terminal_age,
                rho_value = rho_value
              )
            }
          )
        
        
        data.frame(
          
          Rho =
            rho_value,
          
          Age =
            selected_consumption_ages,
          
          Effective_Annuity_Factor =
            abar_values,
          
          Consumption_Wealth_Ratio =
            1 /
            abar_values,
          
          Consumption_Wealth_Percent =
            100 /
            abar_values
        )
      }
    )
  )
  
  
  print(
    rho_selected_age_summary
  )
  
  
  write.csv(
    rho_selected_age_summary,
    file.path(
      results_dir,
      "post_retirement_rho_selected_age_summary.csv"
    ),
    row.names = FALSE
  )  
  
  # ==========================================================
  # 24B. EFFECT OF RHO ON INITIAL RETIREMENT CONSUMPTION
  # ==========================================================
  
  rho_initial_consumption_summary <- data.frame(
    
    Rho =
      rho_values,
    
    Initial_Age =
      t0,
    
    Initial_Wealth =
      X0,
    
    Effective_Annuity_Factor =
      NA_real_,
    
    Initial_Consumption_Wealth_Ratio =
      NA_real_,
    
    Initial_Consumption_Wealth_Percent =
      NA_real_,
    
    Initial_Annual_Consumption =
      NA_real_
  )
  
  
  for (
    i in seq_along(
      rho_values
    )
  ) {
    
    rho_i <-
      rho_values[i]
    
    
    abar_i <-
      annuity_factor(
        t = t0,
        gamma = gamma_hat_rep,
        N = terminal_age,
        rho_value = rho_i
      )
    
    
    rho_initial_consumption_summary$
      Effective_Annuity_Factor[i] <-
      abar_i
    
    rho_initial_consumption_summary$
      Initial_Consumption_Wealth_Ratio[i] <-
      1 /
      abar_i
    
    rho_initial_consumption_summary$
      Initial_Consumption_Wealth_Percent[i] <-
      100 /
      abar_i
    
    rho_initial_consumption_summary$
      Initial_Annual_Consumption[i] <-
      X0 /
      abar_i
  }
  
  
  print(
    rho_initial_consumption_summary
  )
  
  
  write.csv(
    rho_initial_consumption_summary,
    file.path(
      results_dir,
      "post_retirement_rho_initial_consumption_summary.csv"
    ),
    row.names = FALSE
  )  
  
  
  # ==========================================================
  # 25. OPTIMAL WEALTH FIGURE
  # ==========================================================
  
  optimal_wealth_plot <- ggplot(
    consumption_path_data,
    aes(
      x = Age,
      y = Optimal_Wealth
    )
  ) +
    
    geom_line(
      linewidth = 0.85
    ) +
    
    scale_y_continuous(
      labels = scales::comma,
      breaks = scales::pretty_breaks(
        n = 6
      )
    ) +
    
    scale_x_continuous(
      breaks = scales::pretty_breaks(
        n = 8
      )
    ) +
    
    labs(
      x = "Age",
      y = "Optimal pension wealth (DKK)"
    ) +
    
    theme_bw()
  
  
  optimal_wealth_plot
  
  
  ggsave(
    filename = file.path(
      figure_dir,
      "Optimal Post-Retirement Wealth - Representative Member.png"
    ),
    plot = optimal_wealth_plot,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  
  # ==========================================================
  # 26. REPRESENTATIVE-MEMBER SUMMARY
  # ==========================================================
  
  representative_consumption_summary <- data.frame(
    
    Police_Number =
      police_rep,
    
    Initial_Age =
      t0,
    
    Initial_Wealth =
      X0,
    
    Risk_Profile =
      profile_rep,
    
    Gamma_Hat =
      gamma_hat_rep,
    
    Rho =
      rho,
    
    Effective_Rate_Nu =
      nu_rep,
    
    Initial_Effective_Annuity_Factor =
      abar0,
    
    Initial_Optimal_Consumption =
      initial_optimal_consumption,
    
    Initial_Consumption_Wealth_Ratio =
      1 /
      abar0,
    
    Terminal_Age =
      terminal_age
  )
  
  
  print(
    representative_consumption_summary
  )
  
  
  write.csv(
    representative_consumption_summary,
    file.path(
      results_dir,
      "post_retirement_representative_member_summary.csv"
    ),
    row.names = FALSE
  )
  
} else {
  
  representative_consumption_summary <-
    data.frame()
  
  cat(
    "\nNo post-retirement observations are available ",
    "for the representative-member consumption illustration.\n",
    sep = ""
  )
}


# ============================================================
# 27. FINAL INTERNAL CHECKS
# ============================================================

# ------------------------------------------------------------
# All pre-retirement observations must have age < T
# ------------------------------------------------------------

if (
  any(
    pre_retirement_sample$alder >=
    retirement_age
  )
) {
  
  warning(
    "The pre-retirement comparison sample contains observations at or after T."
  )
}


# ------------------------------------------------------------
# All post-retirement observations must have age >= T
# ------------------------------------------------------------

if (
  nrow(
    post_retirement_sample
  ) > 0 &&
  any(
    post_retirement_sample$alder <
    retirement_age
  )
) {
  
  warning(
    "The post-retirement comparison sample contains observations before T."
  )
}


# ------------------------------------------------------------
# h_i must equal zero after retirement
# ------------------------------------------------------------

post_retirement_nonzero_h <-
  subset(
    post_retirement_sample,
    abs(
      h_i
    ) >
      1e-10
  )


if (
  nrow(
    post_retirement_nonzero_h
  ) > 0
) {
  
  warning(
    "At least one post-retirement observation has non-zero h_i."
  )
}


# ------------------------------------------------------------
# Theoretical allocation must converge to Merton allocation
# when h = 0
# ------------------------------------------------------------

test_gamma <-
  median(
    comparison_sample$
      gamma_hat_i,
    na.rm = TRUE
  )


test_w_M <-
  merton_allocation(
    test_gamma
  )


test_w_star_at_T <-
  optimal_pre_retirement_allocation(
    X = 1,
    h = 0,
    gamma = test_gamma
  )


merton_boundary_error <-
  max(
    abs(
      test_w_M -
        test_w_star_at_T
    )
  )


if (
  merton_boundary_error >
  1e-12
) {
  
  warning(
    "The pre-retirement allocation does not converge numerically to the Merton allocation when h = 0."
  )
}


# ============================================================
# 28. SAVE KEY RESULTS
# ============================================================

capture.output({
  
  cat(
    "THEORY-PRACTICE COMPARISON RESULTS\n"
  )
  
  cat(
    "==================================\n\n"
  )
  
  
  cat(
    "Retirement boundary T:",
    retirement_age,
    "\n"
  )
  
  cat(
    "Terminal age N:",
    terminal_age,
    "\n"
  )
  
  cat(
    "Risk-free rate r:",
    r,
    "\n"
  )
  
  cat(
    "Baseline subjective discount rate rho:",
    rho,
    "\n"
  )
  
  cat(
    "Theta squared:",
    Theta_sq,
    "\n\n"
  )
  
  
  cat(
    "Final calibrated sample:",
    nrow(
      comparison_sample
    ),
    "\n"
  )
  
  cat(
    "Pre-retirement comparison sample:",
    nrow(
      pre_retirement_sample
    ),
    "\n"
  )
  
  cat(
    "Post-retirement comparison sample:",
    nrow(
      post_retirement_sample
    ),
    "\n\n"
  )
  
  
  cat(
    "Portfolio comparison summary:\n"
  )
  
  print(
    as.data.frame(
      portfolio_comparison_summary
    ),
    row.names = FALSE
  )
  
  
  cat(
    "\nPre-retirement welfare summary:\n"
  )
  
  print(
    as.data.frame(
      welfare_summary_by_profile
    ),
    row.names = FALSE
  )
  
  
  cat(
    "\nPost-retirement welfare results:\n"
  )
  
  print(
    post_retirement_results,
    row.names = FALSE
  )
  
  
  cat(
    "\nPortfolio diagnostics:\n"
  )
  
  cat(
    "\nTerminal P+ age-70 implied risk aversion:\n"
  )
  
  print(
    terminal_gamma_results,
    row.names = FALSE
  )
  
  
  cat(
    "\nTerminal P+ versus unrestricted theoretical portfolio diagnostics:\n"
  )
  
  print(
    terminal_portfolio_diagnostics,
    row.names = FALSE
  )
  
  
  cat(
    "\nTerminal customised-fund calibration check:\n"
  )
  
  cat(
    "Maximum |recovered optimal fund exposure - 1|:",
    terminal_fund_exposure_error,
    "\n"
  )
  
  cat(
    "Theoretical portfolios with total risky weight > 1:",
    sum(
      pre_retirement_sample$
        total_risky_weight_optimal > 1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Theoretical portfolios with negative risk-free weight:",
    sum(
      pre_retirement_sample$
        riskfree_weight_optimal < 0,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "\nBoundary checks:\n"
  )
  
  cat(
    "Post-retirement observations with non-zero h_i:",
    nrow(
      post_retirement_nonzero_h
    ),
    "\n"
  )
  
  cat(
    "Maximum difference between w*(h=0) and w^M:",
    merton_boundary_error,
    "\n"
  )
  
  
  if (
    nrow(
      representative_consumption_summary
    ) > 0
  ) {
    
    cat(
      "\nRepresentative-member consumption summary:\n"
    )
    
    print(
      representative_consumption_summary,
      row.names = FALSE
    )
  }
  
  cat(
    "\nRho sensitivity at selected ages:\n"
  )
  
  print(
    rho_selected_age_summary,
    row.names = FALSE
  )
  
  
  cat(
    "\nEffect of rho on initial optimal consumption:\n")
  
  print(
    rho_initial_consumption_summary,
    row.names = FALSE
  )
  
  
  cat(
    "\nSession information:\n"
  )
  
  print(
    sessionInfo()
  )
  
},
file = file.path(
  results_dir,
  "theory_practice_key_results.txt"
))


# ============================================================
# 29. RESULTS MANIFEST
# ============================================================

write.csv(
  
  data.frame(
    file =
      sort(
        list.files(
          results_dir
        )
      ),
    stringsAsFactors = FALSE
  ),
  
  file.path(
    results_dir,
    "results_manifest.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 30. FINAL CONSOLE OUTPUT
# ============================================================

cat(
  "\n==================================================\n"
)

cat(
  "THEORY-PRACTICE COMPARISON COMPLETE\n"
)

cat(
  "==================================================\n"
)

cat(
  "Retirement boundary T:",
  retirement_age,
  "\n"
)

cat(
  "Final calibrated sample:",
  nrow(
    comparison_sample
  ),
  "\n"
)

cat(
  "Pre-retirement comparison sample:",
  nrow(
    pre_retirement_sample
  ),
  "\n"
)

cat(
  "Post-retirement comparison sample:",
  nrow(
    post_retirement_sample
  ),
  "\n"
)

cat(
  "Theoretical portfolios with negative risk-free weight:",
  sum(
    pre_retirement_sample$
      riskfree_weight_optimal < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Boundary error |w*(h=0)-w^M|:",
  format(
    merton_boundary_error,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Figures written to '",
  figure_dir,
  "'.\n",
  sep = ""
)

cat(
  "Numerical results written to '",
  results_dir,
  "'.\n",
  sep = ""
)

cat(
  "==================================================\n"
)
