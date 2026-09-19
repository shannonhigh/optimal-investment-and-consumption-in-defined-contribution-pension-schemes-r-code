# ============================================================
# CHAPTER 5
# EMPIRICAL CALIBRATION OF IMPLIED RISK AVERSION
# ============================================================
#
# Initial calibration implementation supplied by the project supervisor. The accompanying anonymised P+ observations and
# market inputs relate to the 2023 project data.
#
# Market assumptions:
# Rådet for Afkastforventninger, 2023Q1-2.
#
# The implementation below follows the theoretical framework
# developed in Chapters 2-4.
#
# Before retirement, the deterministic present value of remaining
# premiums is
#
#   h(t) = l/r * [1 - exp(-r(T-t))].
#
# The implied coefficient of relative risk aversion is obtained
# from
#
#   gamma_hat =
#     (alpha_theta - r) / sigma_theta^2
#     * (X + h(t)) / X,
#
# where P+ invests the full pension savings account in the
# constructed risky mutual fund, so pi_hat = 1.
#
# Retirement age in this thesis:
#
#   T = 67.
#
# The empirical calibration is evaluated at each investor's
# currently recorded age. The current pension balance, remaining
# premiums and P+ portfolio therefore refer to the same age.
#
# Mortality is not included in the valuation of future premiums.
# Future premiums are deterministic and stop at retirement.
# ============================================================


# ============================================================
# 1. SETTINGS
# ============================================================

library(scales)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(latex2exp)


# ------------------------------------------------------------
# General options
# ------------------------------------------------------------

options(scipen = 999)
options(digits = 16)

scale <- 100

linetypes <- c(
  "solid",
  "dashed",
  "dotted",
  "dotdash"
)


# ------------------------------------------------------------
# Model and calibration settings
# ------------------------------------------------------------

retirement_age <- 67
calibration_r <- 0.02


# ------------------------------------------------------------
# Files and output directories
# ------------------------------------------------------------

data_file <- "Samlet information til JP.csv"

figure_dir <- "figures"
results_dir <- "results"

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

if (!file.exists(data_file)) {
  stop(
    paste0(
      "Input file not found: '",
      data_file,
      "'. Set the working directory to the project root ",
      "and place the anonymised CSV there before running this script."
    )
  )
}


# ============================================================
# 2. MARKET PARAMETERS
# ============================================================

# Aggregate the ten asset classes into:
#
#   1. bonds     = asset classes 1-4
#   2. stocks    = asset classes 5-6
#   3. illiquid  = asset classes 7-10

pf_fordelingsfunktion <- function(A) {
  rbind(
    sum(A[1:4, 1]),
    sum(A[5:6, 1]),
    sum(A[7:10, 1])
  )
}


# ------------------------------------------------------------
# Correlation matrix
# Rådet for Afkastforventninger
# ------------------------------------------------------------

correlation_matrix <- matrix(
  c(
    1.0,  0.6,  0.1,  0.3, -0.1, -0.1, -0.2, -0.1, -0.1, -0.1,
    0.6,  1.0,  0.6,  0.6,  0.2,  0.2,  0.2,  0.1,  0.1,  0.3,
    0.1,  0.6,  1.0,  0.7,  0.7,  0.6,  0.6,  0.4,  0.3,  0.7,
    0.3,  0.6,  0.7,  1.0,  0.5,  0.6,  0.4,  0.2,  0.2,  0.5,
    -0.1,  0.2,  0.7,  0.5,  1.0,  0.7,  0.8,  0.4,  0.4,  0.8,
    -0.1,  0.2,  0.6,  0.6,  0.7,  1.0,  0.7,  0.4,  0.4,  0.7,
    -0.2,  0.2,  0.6,  0.4,  0.8,  0.7,  1.0,  0.4,  0.4,  0.7,
    -0.1,  0.1,  0.4,  0.2,  0.4,  0.4,  0.4,  1.0,  0.3,  0.4,
    -0.1,  0.1,  0.3,  0.2,  0.4,  0.4,  0.4,  0.3,  1.0,  0.4,
    -0.1,  0.3,  0.7,  0.5,  0.8,  0.7,  0.7,  0.4,  0.4,  1.0
  ),
  nrow = 10,
  ncol = 10,
  byrow = TRUE
)


# ------------------------------------------------------------
# Volatilities
# ------------------------------------------------------------

std_dev <- matrix(
  data = c(
    3.7,
    5.5,
    11.9,
    10.7,
    15.3,
    21.7,
    20.4,
    14.0,
    10.8,
    9.4
  ) / scale,
  nrow = 10,
  ncol = 1
)


# ------------------------------------------------------------
# Variance-covariance matrix
# ------------------------------------------------------------

variance_matrix <- function(t = 0) {
  
  if (t < 11) {
    
    covariance_matrix <- matrix(
      0,
      nrow = 10,
      ncol = 10
    )
    
    for (i in 1:10) {
      for (j in 1:10) {
        
        covariance_matrix[i, j] <-
          correlation_matrix[i, j] *
          std_dev[i, 1] *
          std_dev[j, 1]
      }
    }
    
    return(covariance_matrix)
  }
  
  return(
    matrix(
      data = c(
        0.08^2, 0,       0,
        0,       0.18^2, 0,
        0,       0,       0.12^2
      ),
      nrow = 3,
      ncol = 3,
      byrow = TRUE
    )
  )
}


# ------------------------------------------------------------
# Expected returns
# ------------------------------------------------------------

returns <- function(t = 0) {
  
  if (t < 6) {
    
    return(
      matrix(
        data = c(
          1.9,
          2.2,
          4.9,
          4.3,
          6.1,
          8.3,
          10.2,
          5.6,
          4.1,
          3.8
        ) / scale,
        nrow = 10,
        ncol = 1
      )
    )
  }
  
  if (t < 11) {
    
    return(
      matrix(
        data = c(
          2.1,
          3.1,
          5.0,
          4.7,
          6.3,
          8.8,
          10.5,
          5.9,
          5.1,
          4.3
        ) / scale,
        nrow = 10,
        ncol = 1
      )
    )
  }
  
  return(
    matrix(
      data = c(
        3.5,
        6.5,
        6.5
      ) / scale,
      nrow = 3,
      ncol = 1
    )
  )
}


# ============================================================
# 3. LOAD P+ DATA
# ============================================================

data <- read.csv(
  data_file,
  header = TRUE,
  sep = ";",
  dec = ",",
  colClasses = c(
    "character",
    "character",
    "character",
    "character",
    "character",
    "character",
    "numeric",
    "numeric",
    "numeric",
    "character"
  )
)


# ============================================================
# 4. RISK-AVERSION PROFILE CODING
# ============================================================

# Internal coding:
#
# "hoj"    = low risk aversion / higher-risk allocation
# "mellem" = moderate risk aversion
# "lav"    = high risk aversion / lower-risk allocation

data$investeringsProfil[
  data$investeringsProfil != "P+ Livscyklus mellem" &
    data$investeringsProfil != "P+ Livscyklus lav"
] <- "hoj"

data$investeringsProfil[
  data$investeringsProfil == "P+ Livscyklus mellem"
] <- "mellem"

data$investeringsProfil[
  data$investeringsProfil == "P+ Livscyklus lav"
] <- "lav"


# ============================================================
# 5. DESCRIPTIVE P+ DATA
# ============================================================

# Active P211R policyholders with positive pension savings.

test_data <- subset(
  data,
  gruppeNavn == "P211R" &
    policeStatus == "Betalende" &
    depot > 0
)


# ------------------------------------------------------------
# Age distribution
# ------------------------------------------------------------

d_age <- ggplot(
  test_data,
  aes(x = alder)
) +
  geom_histogram(
    bins = 30,
    colour = "grey40",
    fill = "grey75"
  ) +
  ylab("Frequency") +
  xlab("Age") +
  theme_bw()

d_age

ggsave(
  filename = file.path(
    figure_dir,
    "Age Distribution of P211R Policyholders.png"
  ),
  plot = d_age,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# Risk-aversion profile distribution
# ------------------------------------------------------------

d_profile <- ggplot(
  test_data,
  aes(
    x = fct_relevel(
      investeringsProfil,
      "hoj",
      "mellem",
      "lav"
    )
  )
) +
  geom_bar(
    colour = "grey40",
    fill = "grey75"
  ) +
  ylab("Frequency") +
  xlab("Risk aversion profile") +
  scale_x_discrete(
    labels = c(
      hoj = "Low",
      mellem = "Moderate",
      lav = "High"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 15,
      hjust = 1
    )
  )

d_profile

ggsave(
  filename = file.path(
    figure_dir,
    "Distribution of Risk Profiles, P+ Livscyklus.png"
  ),
  plot = d_profile,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# Pension savings distribution
# ------------------------------------------------------------

d_savings <- ggplot(
  test_data,
  aes(x = depot / 1000000)
) +
  geom_histogram(
    binwidth = 0.25,
    boundary = 0,
    colour = "grey40",
    fill = "grey75"
  ) +
  coord_cartesian(
    xlim = c(0, 5)
  ) +
  labs(
    x = "Pension savings (DKK millions)",
    y = "Frequency"
  ) +
  theme_bw()

d_savings

ggsave(
  filename = file.path(
    figure_dir,
    "Distribution of Pension Savings (Depot).png"
  ),
  plot = d_savings,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 6. P+ LIFE-CYCLE PORTFOLIOS
# ============================================================

# The supplied portfolio weights do not always sum exactly to one.
# Each allocation is therefore rescaled to sum to one.

rescale_function <- function(A) {
  
  total_weight <- sum(A)
  
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("Portfolio weights must have a finite positive sum.")
  }
  
  as.matrix(
    A / total_weight
  )
}


# ------------------------------------------------------------
# Low risk-aversion profile ("hoj")
# ------------------------------------------------------------

pension_om_30_hoj_original <- matrix(
  c(
    5.0, 14.0, 25.0, 6.0, 69.0,
    11.0, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_30_hoj <-
  rescale_function(pension_om_30_hoj_original)


pension_om_15_hoj_original <- matrix(
  c(
    5.0, 14.0, 25.0, 6.0, 69.0,
    11.0, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_15_hoj <-
  rescale_function(pension_om_15_hoj_original)


pension_om_5_hoj_original <- matrix(
  c(
    14.0, 16.8, 22.2, 6.0, 50.4,
    7.9, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_5_hoj <-
  rescale_function(pension_om_5_hoj_original)


pension_efter_5_hoj_original <- matrix(
  c(
    18.3, 18.1, 20.9, 6.0, 41.4,
    6.4, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_efter_5_hoj <-
  rescale_function(pension_efter_5_hoj_original)


# ------------------------------------------------------------
# Moderate risk-aversion profile ("mellem")
# ------------------------------------------------------------

pension_om_30_mellem_original <- matrix(
  c(
    10.8, 15.8, 23.2, 6.0, 57.0,
    9.0, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_30_mellem <-
  rescale_function(pension_om_30_mellem_original)


pension_om_15_mellem_original <- matrix(
  c(
    10.8, 15.8, 23.2, 6.0, 57.0,
    9.0, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_15_mellem <-
  rescale_function(pension_om_15_mellem_original)


pension_om_5_mellem_original <- matrix(
  c(
    21.8, 19.2, 19.8, 6.0, 34.3,
    5.2, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_5_mellem <-
  rescale_function(pension_om_5_mellem_original)


pension_efter_5_mellem_original <- matrix(
  c(
    25.2, 20.3, 18.7, 6.0, 27.1,
    4.0, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_efter_5_mellem <-
  rescale_function(pension_efter_5_mellem_original)


# ------------------------------------------------------------
# High risk-aversion profile ("lav")
# ------------------------------------------------------------

pension_om_30_lav_original <- matrix(
  c(
    23.3, 19.7, 19.3, 6.0, 31.2,
    4.7, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_30_lav <-
  rescale_function(pension_om_30_lav_original)


pension_om_15_lav_original <- matrix(
  c(
    23.3, 19.7, 19.3, 6.0, 31.2,
    4.7, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_15_lav <-
  rescale_function(pension_om_15_lav_original)


pension_om_5_lav_original <- matrix(
  c(
    30.4, 21.9, 17.1, 6.0, 16.5,
    2.3, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_om_5_lav <-
  rescale_function(pension_om_5_lav_original)


pension_efter_5_lav_original <- matrix(
  c(
    32.9, 22.7, 16.3, 6.0, 11.2,
    1.4, 4.0, 6.0, 10.0, 6.0
  ) / scale,
  nrow = 10,
  ncol = 1
)

pension_efter_5_lav <-
  rescale_function(pension_efter_5_lav_original)


# ============================================================
# 7. P+ PORTFOLIO AS A FUNCTION OF AGE
# ============================================================

# P+ allocations are linearly interpolated between the supplied
# age points. The portfolio continues to de-risk until age 70.
#
# The age-70 portfolio rule is separate from the retirement-age
# assumption T = 67 used to value future premiums.

pf <- function(alder, risiko) {
  
  if (
    !is.finite(alder)
  ) {
    stop("Age must be finite.")
  }
  
  if (
    !risiko %in%
    c(
      "hoj",
      "mellem",
      "lav"
    )
  ) {
    stop(
      paste(
        "Unknown risk-aversion profile:",
        risiko
      )
    )
  }
  
  
  profile_column <- switch(
    risiko,
    hoj = 1,
    mellem = 2,
    lav = 3
  )
  
  
  portfolios <- matrix(
    data = list(
      pension_om_30_hoj,
      pension_om_15_hoj,
      pension_om_5_hoj,
      pension_efter_5_hoj,
      
      pension_om_30_mellem,
      pension_om_15_mellem,
      pension_om_5_mellem,
      pension_efter_5_mellem,
      
      pension_om_30_lav,
      pension_om_15_lav,
      pension_om_5_lav,
      pension_efter_5_lav
    ),
    nrow = 4,
    ncol = 3
  )
  
  
  if (alder < 50) {
    
    lower_row <- 1
    upper_row <- 1
    lower_age <- 0
    upper_age <- 50
    
  } else if (alder < 60) {
    
    lower_row <- 2
    upper_row <- 3
    lower_age <- 50
    upper_age <- 60
    
  } else if (alder < 70) {
    
    lower_row <- 3
    upper_row <- 4
    lower_age <- 60
    upper_age <- 70
    
  } else {
    
    lower_row <- 4
    upper_row <- 4
    lower_age <- 70
    upper_age <- 100
  }
  
  
  lower_portfolio <-
    portfolios[
      lower_row,
      profile_column
    ][[1]]
  
  upper_portfolio <-
    portfolios[
      upper_row,
      profile_column
    ][[1]]
  
  
  if (lower_row == upper_row) {
    
    return(
      as.matrix(
        lower_portfolio
      )
    )
  }
  
  
  interpolation_weight <-
    (alder - lower_age) /
    (upper_age - lower_age)
  
  
  interpolated_portfolio <-
    (1 - interpolation_weight) *
    lower_portfolio +
    interpolation_weight *
    upper_portfolio
  
  
  as.matrix(
    interpolated_portfolio
  )
}


# ============================================================
# 8. CHECK SUPPLIED PORTFOLIOS
# ============================================================

# Check that interpolated P+ portfolios sum to one over the
# age range used for plotting.

portfolio_check_ages <- seq(
  25,
  80,
  by = 1
)

portfolio_check_profiles <- c(
  "hoj",
  "mellem",
  "lav"
)

portfolio_sum_check <- expand.grid(
  age = portfolio_check_ages,
  profile = portfolio_check_profiles,
  stringsAsFactors = FALSE
)

portfolio_sum_check$weight_sum <- mapply(
  function(age, profile) {
    sum(
      pf(
        age,
        profile
      )
    )
  },
  portfolio_sum_check$age,
  portfolio_sum_check$profile
)

if (
  any(
    abs(
      portfolio_sum_check$weight_sum - 1
    ) > 1e-10
  )
) {
  warning(
    "At least one interpolated P+ portfolio does not sum to one."
  )
}


# ============================================================
# 9. PLOT P+ LIFE-CYCLE STRATEGY
# ============================================================

stocks <- data.frame(
  age = seq(
    25,
    80,
    by = 1
  )
)


stocks$Low <- sapply(
  stocks$age,
  function(x) {
    sum(
      pf(
        x,
        "hoj"
      )[5:10]
    )
  }
)


stocks$Moderate <- sapply(
  stocks$age,
  function(x) {
    sum(
      pf(
        x,
        "mellem"
      )[5:10]
    )
  }
)


stocks$High <- sapply(
  stocks$age,
  function(x) {
    sum(
      pf(
        x,
        "lav"
      )[5:10]
    )
  }
)


scaleFUN <- function(x) {
  sprintf(
    "%.2f",
    x
  )
}


data_long <- melt(
  stocks,
  id = "age"
)


gfg_plot <- ggplot(
  data_long,
  aes(
    x = age,
    y = round(
      value,
      digits = 4
    ),
    group = variable
  )
) +
  geom_line(
    aes(
      linetype = variable
    )
  ) +
  ylab(
    TeX("$w_{5-10}$")
  ) +
  xlab("Age") +
  scale_linetype_manual(
    values = linetypes,
    name = "Risk aversion profile"
  ) +
  theme_bw() +
  scale_y_continuous(
    breaks = scales::pretty_breaks(
      n = 5
    ),
    labels = scaleFUN
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 5
    )
  )


gfg_plot


ggsave(
  filename = file.path(
    figure_dir,
    "P+ Livscyklus Growth-Asset Weight by Age and Risk Profile.png"
  ),
  plot = gfg_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# ============================================================
# 10. CALIBRATION SAMPLE
# ============================================================

# The broad calibration sample contains P211R observations:
#
#   age >= 50,
#   pension savings > 0.
#
# No upper-age restriction is imposed at this stage.
#
# This preserves the broad sample of 569 observations and permits
# diagnostics around the retirement boundary. The P+ portfolio
# continues to change until age 70, while deterministic premiums
# stop at the thesis retirement age T = 67.

thepensiondata <- subset(
  data,
  gruppeNavn == "P211R" &
    depot > 0 &
    alder >= 50
)


# ============================================================
# 11. PRESENT VALUE OF FUTURE PREMIUMS
# ============================================================

# For age a < T:
#
#                 l
# h(a) = ---------------- [1 - exp(-r(T-a))]
#                 r
#
# where:
#
#   l = annual premium rate
#   r = risk-free rate
#   T = retirement age
#
# For age a >= T:
#
#   h(a) = 0.
#
# The premium stream is deterministic. No mortality adjustment
# is included.

future_premiums <- function(
    age,
    annual_premium,
    retirement_age = 67,
    r = 0.02
) {
  
  if (!is.finite(age)) {
    return(NA_real_)
  }
  
  if (!is.finite(annual_premium)) {
    return(NA_real_)
  }
  
  if (!is.finite(retirement_age)) {
    stop("Retirement age must be finite.")
  }
  
  if (!is.finite(r)) {
    stop("Risk-free rate must be finite.")
  }
  
  
  years_to_retirement <-
    max(
      retirement_age - age,
      0
    )
  
  
  if (
    annual_premium <= 0 ||
    years_to_retirement <= 0
  ) {
    
    return(0)
  }
  
  
  if (abs(r) < 1e-12) {
    
    return(
      annual_premium *
        years_to_retirement
    )
  }
  
  
  (
    annual_premium / r
  ) *
    (
      1 -
        exp(
          -r *
            years_to_retirement
        )
    )
}


# ============================================================
# 12. INVERSE CALIBRATION OF RISK AVERSION
# ============================================================

# For investor i:
#
#   X_i       = current pension savings
#   P_i       = current annual premium rate
#   h_i       = PV of remaining premiums at current age
#   theta_i   = P+ portfolio at current age
#
# The customised risky mutual fund has
#
#   alpha_theta,i =
#     theta_i' alpha
#
# and
#
#   sigma_theta,i^2 =
#     theta_i' Sigma theta_i.
#
# The theoretical pre-retirement risky allocation is
#
#   pi_i =
#
#     (alpha_theta,i - r)
#     ---------------------
#     gamma_i sigma_theta,i^2
#
#     * (X_i + h_i) / X_i.
#
# In the inverse calibration, the observed P+ portfolio itself is
# treated as the customised risky mutual fund. The full savings
# account is invested in this fund, so pi_i = 1.
#
# Therefore
#
#   gamma_hat_i =
#
#     (alpha_theta,i - r)
#     ---------------------
#       sigma_theta,i^2
#
#     * (X_i + h_i) / X_i.
#
# Current pension savings, future premiums and the P+ portfolio
# are all evaluated at the investor's currently recorded age.

calibrate_gamma <- function(
    y,
    r = 0.02,
    retirement_age = 67
) {
  
  y$annual_premium <- NA_real_
  y$future_premiums <- NA_real_
  y$portfolio_return <- NA_real_
  y$portfolio_variance <- NA_real_
  y$gamma <- NA_real_
  
  
  for (row in seq_len(nrow(y))) {
    
    # --------------------------------------------------------
    # Individual data
    # --------------------------------------------------------
    
    age <-
      y$alder[row]
    
    monthly_premium <-
      y$maanedligPraemie[row]
    
    depot <-
      y$depot[row]
    
    risk_profile <-
      y$investeringsProfil[row]
    
    
    # --------------------------------------------------------
    # Validate individual data
    # --------------------------------------------------------
    
    if (
      !is.finite(age) ||
      !is.finite(depot) ||
      depot <= 0
    ) {
      
      warning(
        paste(
          "Invalid age or pension savings for row",
          row
        )
      )
      
      next
    }
    
    
    if (!is.finite(monthly_premium)) {
      
      warning(
        paste(
          "Missing or non-finite premium for row",
          row
        )
      )
      
      next
    }
    
    
    annual_premium <-
      12 *
      monthly_premium
    
    
    # --------------------------------------------------------
    # Present value of future premiums at current age
    # --------------------------------------------------------
    
    h <- future_premiums(
      age = age,
      annual_premium = annual_premium,
      retirement_age = retirement_age,
      r = r
    )
    
    
    if (!is.finite(h)) {
      
      warning(
        paste(
          "Invalid present value of future premiums for row",
          row
        )
      )
      
      next
    }
    
    
    # --------------------------------------------------------
    # P+ portfolio at current age
    # --------------------------------------------------------
    
    theta <- as.matrix(
      pf(
        age,
        risk_profile
      )
    )
    
    
    # --------------------------------------------------------
    # Portfolio checks
    # --------------------------------------------------------
    
    if (
      nrow(theta) != 10 ||
      ncol(theta) != 1
    ) {
      
      warning(
        paste(
          "Invalid portfolio dimension for row",
          row
        )
      )
      
      next
    }
    
    
    if (
      abs(
        sum(theta) - 1
      ) > 1e-10
    ) {
      
      warning(
        paste(
          "Portfolio weights do not sum to one for row",
          row
        )
      )
    }
    
    
    # --------------------------------------------------------
    # Expected return
    # --------------------------------------------------------
    
    alpha_theta <- as.numeric(
      t(theta) %*%
        returns(1)
    )
    
    
    # --------------------------------------------------------
    # Variance
    # --------------------------------------------------------
    
    sigma2_theta <- as.numeric(
      t(theta) %*%
        variance_matrix(1) %*%
        theta
    )
    
    
    # --------------------------------------------------------
    # Validate portfolio moments
    # --------------------------------------------------------
    
    if (!is.finite(alpha_theta)) {
      
      warning(
        paste(
          "Invalid portfolio expected return for row",
          row
        )
      )
      
      next
    }
    
    
    if (
      !is.finite(sigma2_theta) ||
      sigma2_theta <= 0
    ) {
      
      warning(
        paste(
          "Invalid portfolio variance for row",
          row
        )
      )
      
      next
    }
    
    
    # --------------------------------------------------------
    # Implied coefficient of relative risk aversion
    # --------------------------------------------------------
    
    gamma_hat <-
      (
        (alpha_theta - r) /
          sigma2_theta
      ) *
      (
        (depot + h) /
          depot
      )
    
    
    # --------------------------------------------------------
    # Store results
    # --------------------------------------------------------
    
    y$annual_premium[row] <-
      annual_premium
    
    y$future_premiums[row] <-
      h
    
    y$portfolio_return[row] <-
      alpha_theta
    
    y$portfolio_variance[row] <-
      sigma2_theta
    
    y$gamma[row] <-
      gamma_hat
  }
  
  
  y
}


# ------------------------------------------------------------
# Run calibration
# ------------------------------------------------------------

thepensiondatainclgamma <- calibrate_gamma(
  thepensiondata,
  r = calibration_r,
  retirement_age = retirement_age
)


# ============================================================
# 13. PREMIUM AND WEALTH DIAGNOSTICS
# ============================================================

# Ratio of current pension savings to current annual premium.

thepensiondatainclgamma$depot_to_annual_premium <-
  thepensiondatainclgamma$depot /
  thepensiondatainclgamma$annual_premium


# The ratio is undefined for observations with a zero or negative
# annual premium.

thepensiondatainclgamma$depot_to_annual_premium[
  !is.finite(
    thepensiondatainclgamma$annual_premium
  ) |
    thepensiondatainclgamma$annual_premium <= 0
] <- NA_real_


zero_premium_count <- sum(
  is.finite(
    thepensiondatainclgamma$annual_premium
  ) &
    thepensiondatainclgamma$annual_premium <= 0,
  na.rm = TRUE
)


missing_premium_count <- sum(
  !is.finite(
    thepensiondatainclgamma$annual_premium
  )
)


cat(
  "Zero-premium observations:",
  zero_premium_count,
  "\n"
)

cat(
  "Missing/non-finite premium observations:",
  missing_premium_count,
  "\n"
)


# ============================================================
# 14. INTERNAL CONSISTENCY CHECKS
# ============================================================

# ------------------------------------------------------------
# Future premiums at retirement
# ------------------------------------------------------------

# h(a) must equal zero for a >= T.

post_retirement_premium_check <- subset(
  thepensiondatainclgamma,
  alder >= retirement_age &
    is.finite(future_premiums) &
    abs(future_premiums) > 1e-10
)


if (
  nrow(
    post_retirement_premium_check
  ) > 0
) {
  
  warning(
    paste(
      nrow(
        post_retirement_premium_check
      ),
      "observations at or after retirement have non-zero future premiums."
    )
  )
}


# ------------------------------------------------------------
# Non-negative future premiums
# ------------------------------------------------------------

negative_future_premiums <- subset(
  thepensiondatainclgamma,
  is.finite(future_premiums) &
    future_premiums < -1e-10
)


if (
  nrow(
    negative_future_premiums
  ) > 0
) {
  
  warning(
    "Negative present values of future premiums detected."
  )
}


# ------------------------------------------------------------
# Portfolio weights
# ------------------------------------------------------------

portfolio_weight_check <- sapply(
  seq_len(
    nrow(
      thepensiondatainclgamma
    )
  ),
  function(row) {
    
    sum(
      pf(
        thepensiondatainclgamma$alder[row],
        thepensiondatainclgamma$investeringsProfil[row]
      )
    )
  }
)


if (
  any(
    abs(
      portfolio_weight_check - 1
    ) > 1e-10
  )
) {
  
  warning(
    "At least one P+ portfolio does not sum to one."
  )
}


# ------------------------------------------------------------
# Non-finite gamma
# ------------------------------------------------------------

non_finite_gamma <- subset(
  thepensiondatainclgamma,
  !is.finite(gamma)
)


if (
  nrow(
    non_finite_gamma
  ) > 0
) {
  
  warning(
    paste(
      nrow(
        non_finite_gamma
      ),
      "observations have non-finite calibrated gamma values."
    )
  )
}


# ------------------------------------------------------------
# Non-positive gamma
# ------------------------------------------------------------

non_positive_gamma <- subset(
  thepensiondatainclgamma,
  is.finite(gamma) &
    gamma <= 0
)


if (
  nrow(
    non_positive_gamma
  ) > 0
) {
  
  warning(
    paste(
      nrow(
        non_positive_gamma
      ),
      "observations have non-positive calibrated gamma values."
    )
  )
}


# ============================================================
# 15. AGE COMPOSITION
# ============================================================

age_summary <- thepensiondatainclgamma %>%
  summarise(
    
    N_total = n(),
    
    N_age_50_to_66 =
      sum(
        alder >= 50 &
          alder < 67
      ),
    
    N_age_67_to_69 =
      sum(
        alder >= 67 &
          alder < 70
      ),
    
    N_age_70_plus =
      sum(
        alder >= 70
      )
  )


age_summary


# ============================================================
# 16. UNRESTRICTED GAMMA DISTRIBUTION
# ============================================================

res_unrestricted <-
  as.data.frame(
    thepensiondatainclgamma
  )


res_unrestricted$investeringsProfil <- factor(
  res_unrestricted$investeringsProfil,
  levels = c(
    "hoj",
    "mellem",
    "lav"
  )
)


d_gamma_unrestricted <- ggplot(
  res_unrestricted,
  aes(
    x = gamma,
    fill = investeringsProfil
  )
) +
  geom_histogram(
    bins = 30,
    position = "identity",
    alpha = 0.6,
    colour = "white"
  ) +
  ylab("Frequency") +
  xlab(
    TeX("$\\gamma$")
  ) +
  scale_fill_manual(
    name = "Risk aversion profile",
    breaks = c(
      "hoj",
      "mellem",
      "lav"
    ),
    labels = c(
      "Low",
      "Moderate",
      "High"
    ),
    values = c(
      "#858484",
      "gray",
      "black"
    )
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 8
    )
  )


d_gamma_unrestricted


ggsave(
  filename = file.path(
    figure_dir,
    "Calibrated Risk Aversion - No Depot Restriction.png"
  ),
  plot = d_gamma_unrestricted,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 17. WEALTH-TO-PREMIUM DISTRIBUTION
# ============================================================

summary(
  thepensiondatainclgamma$
    depot_to_annual_premium
)


quantile(
  thepensiondatainclgamma$
    depot_to_annual_premium,
  
  probs = c(
    0,
    0.01,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  ),
  
  na.rm = TRUE
)
# ------------------------------------------------------------
# Wealth-to-premium ratio histogram
# ------------------------------------------------------------

wealth_premium_plot_data <- thepensiondatainclgamma %>%
  filter(
    is.finite(depot_to_annual_premium),
    depot_to_annual_premium >= 0
  )


d_wealth_premium <- ggplot(
  wealth_premium_plot_data,
  aes(
    x = depot_to_annual_premium
  )
) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    colour = "grey40",
    fill = "grey75"
  ) +
  geom_vline(
    xintercept = 5,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  coord_cartesian(
    xlim = c(0, 20)
  ) +
  labs(
    x = TeX("$R_i=X_i/P_i$"),
    y = "Frequency"
  ) +
  theme_bw()


d_wealth_premium


ggsave(
  filename = file.path(
    figure_dir,
    "Distribution of Wealth-to-Premium Ratio.png"
  ),
  plot = d_wealth_premium,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 18. ALTERNATIVE SAVINGS RESTRICTIONS
# ============================================================

# The restriction compares current pension savings X_i with
# multiples of the current annual premium P_i.
#
# It does not compare X_i with the present value h_i.

gamma_no_restriction <-
  thepensiondatainclgamma


gamma_1yr <- subset(
  thepensiondatainclgamma,
  is.finite(annual_premium) &
    depot >
    1 *
    annual_premium
)


gamma_2yr <- subset(
  thepensiondatainclgamma,
  is.finite(annual_premium) &
    depot >
    2 *
    annual_premium
)


gamma_5yr <- subset(
  thepensiondatainclgamma,
  is.finite(annual_premium) &
    depot >
    5 *
    annual_premium
)


N_unrestricted <-
  nrow(
    gamma_no_restriction
  )


restriction_summary <- data.frame(
  
  Restriction = c(
    "No restriction",
    "Depot > 1 annual premium",
    "Depot > 2 annual premiums",
    "Depot > 5 annual premiums"
  ),
  
  N = c(
    nrow(
      gamma_no_restriction
    ),
    nrow(
      gamma_1yr
    ),
    nrow(
      gamma_2yr
    ),
    nrow(
      gamma_5yr
    )
  )
)


restriction_summary$Removed <-
  N_unrestricted -
  restriction_summary$N


restriction_summary$Removed_Percent <-
  100 *
  restriction_summary$Removed /
  N_unrestricted


restriction_summary


# ============================================================
# 19. GAMMA STATISTICS BY RESTRICTION
# ============================================================

gamma_summary_function <- function(
    dataset,
    restriction_name
) {
  
  valid_gamma <-
    dataset$gamma[
      is.finite(
        dataset$gamma
      )
    ]
  
  
  if (
    length(
      valid_gamma
    ) == 0
  ) {
    
    return(
      data.frame(
        Restriction = restriction_name,
        N = nrow(dataset),
        N_Gamma = 0,
        Mean_Gamma = NA_real_,
        Median_Gamma = NA_real_,
        SD_Gamma = NA_real_,
        Min_Gamma = NA_real_,
        Q25_Gamma = NA_real_,
        Q75_Gamma = NA_real_,
        Max_Gamma = NA_real_
      )
    )
  }
  
  
  data.frame(
    
    Restriction =
      restriction_name,
    
    N =
      nrow(
        dataset
      ),
    
    N_Gamma =
      length(
        valid_gamma
      ),
    
    Mean_Gamma =
      mean(
        valid_gamma
      ),
    
    Median_Gamma =
      median(
        valid_gamma
      ),
    
    SD_Gamma =
      if (
        length(
          valid_gamma
        ) > 1
      ) {
        sd(
          valid_gamma
        )
      } else {
        NA_real_
      },
    
    Min_Gamma =
      min(
        valid_gamma
      ),
    
    Q25_Gamma =
      as.numeric(
        quantile(
          valid_gamma,
          0.25
        )
      ),
    
    Q75_Gamma =
      as.numeric(
        quantile(
          valid_gamma,
          0.75
        )
      ),
    
    Max_Gamma =
      max(
        valid_gamma
      )
  )
}


gamma_restriction_summary <- rbind(
  
  gamma_summary_function(
    gamma_no_restriction,
    "No restriction"
  ),
  
  gamma_summary_function(
    gamma_1yr,
    "Depot > 1 annual premium"
  ),
  
  gamma_summary_function(
    gamma_2yr,
    "Depot > 2 annual premiums"
  ),
  
  gamma_summary_function(
    gamma_5yr,
    "Depot > 5 annual premiums"
  )
)


gamma_restriction_summary


# ============================================================
# 20. COMPARE UNRESTRICTED AND FIVE-PREMIUM SAMPLES
# ============================================================

gamma_no_restriction_plot <-
  gamma_no_restriction

gamma_5yr_plot <-
  gamma_5yr


gamma_no_restriction_plot$Sample <-
  "No depot restriction"

gamma_5yr_plot$Sample <-
  "Depot > 5 annual premiums"


gamma_compare <- rbind(
  gamma_no_restriction_plot,
  gamma_5yr_plot
)


gamma_compare$Sample <- factor(
  gamma_compare$Sample,
  levels = c(
    "No depot restriction",
    "Depot > 5 annual premiums"
  )
)


gamma_compare$investeringsProfil <- factor(
  gamma_compare$investeringsProfil,
  levels = c(
    "hoj",
    "mellem",
    "lav"
  )
)


d_gamma_compare <- ggplot(
  gamma_compare,
  aes(
    x = gamma,
    fill = investeringsProfil
  )
) +
  geom_histogram(
    bins = 30,
    position = "identity",
    alpha = 0.6,
    colour = "white"
  ) +
  facet_wrap(
    ~ Sample,
    scales = "free",
    ncol = 1
  ) +
  ylab("Frequency") +
  xlab(
    TeX("$\\gamma$")
  ) +
  scale_fill_manual(
    name = "Risk aversion profile",
    breaks = c(
      "hoj",
      "mellem",
      "lav"
    ),
    labels = c(
      "Low",
      "Moderate",
      "High"
    ),
    values = c(
      "#858484",
      "gray",
      "black"
    )
  ) +
  theme_bw()


d_gamma_compare


ggsave(
  filename = file.path(
    figure_dir,
    "Calibrated Risk Aversion - No Restriction vs 5 Annual Premiums.png"
  ),
  plot = d_gamma_compare,
  width = 7,
  height = 7,
  dpi = 300
)


# ============================================================
# 21. FINAL FIVE-ANNUAL-PREMIUM SAMPLE
# ============================================================

thepensiondatainclgammayeardeposit <-
  gamma_5yr


res <- as.data.frame(
  thepensiondatainclgammayeardeposit
)


res$investeringsProfil <- factor(
  res$investeringsProfil,
  levels = c(
    "hoj",
    "mellem",
    "lav"
  )
)


# ============================================================
# 22. GAMMA STATISTICS BY RISK-AVERSION PROFILE
# ============================================================

gamma_profile_summary <- gamma_5yr %>%
  
  group_by(
    investeringsProfil
  ) %>%
  
  summarise(
    
    N = n(),
    
    N_Gamma =
      sum(
        is.finite(
          gamma
        )
      ),
    
    Mean_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 0,
        mean(
          gamma[
            is.finite(
              gamma
            )
          ]
        ),
        NA_real_
      ),
    
    Median_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 0,
        median(
          gamma[
            is.finite(
              gamma
            )
          ]
        ),
        NA_real_
      ),
    
    SD_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 1,
        sd(
          gamma[
            is.finite(
              gamma
            )
          ]
        ),
        NA_real_
      ),
    
    Q25_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 0,
        as.numeric(
          quantile(
            gamma[
              is.finite(
                gamma
              )
            ],
            0.25
          )
        ),
        NA_real_
      ),
    
    Q75_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 0,
        as.numeric(
          quantile(
            gamma[
              is.finite(
                gamma
              )
            ],
            0.75
          )
        ),
        NA_real_
      ),
    
    Min_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 0,
        min(
          gamma[
            is.finite(
              gamma
            )
          ]
        ),
        NA_real_
      ),
    
    Max_Gamma =
      ifelse(
        sum(
          is.finite(
            gamma
          )
        ) > 0,
        max(
          gamma[
            is.finite(
              gamma
            )
          ]
        ),
        NA_real_
      ),
    
    .groups = "drop"
  )


gamma_profile_summary


gamma_profile_medians <- gamma_5yr %>%
  
  group_by(
    investeringsProfil
  ) %>%
  
  summarise(
    
    N = n(),
    
    Median_Gamma =
      ifelse(
        any(
          is.finite(
            gamma
          )
        ),
        median(
          gamma[
            is.finite(
              gamma
            )
          ]
        ),
        NA_real_
      ),
    
    .groups = "drop"
  )


gamma_profile_medians


# ============================================================
# 23. FINAL CALIBRATION HISTOGRAM
# ============================================================

d_gamma_final <- ggplot(
  res,
  aes(
    x = gamma,
    fill = investeringsProfil
  )
) +
  geom_histogram(
    bins = 30,
    position = "identity",
    alpha = 0.6,
    colour = "white"
  ) +
  ylab("Frequency") +
  xlab(
    TeX("$\\gamma$")
  ) +
  scale_fill_manual(
    name = "Risk aversion profile",
    breaks = c(
      "hoj",
      "mellem",
      "lav"
    ),
    labels = c(
      "Low",
      "Moderate",
      "High"
    ),
    values = c(
      "#858484",
      "gray",
      "black"
    )
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 8
    )
  )


d_gamma_final


ggsave(
  filename = file.path(
    figure_dir,
    "Calibrated Risk Aversion by Risk Profile - 5 Annual Premium Restriction.png"
  ),
  plot = d_gamma_final,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 24. RETIREMENT-AGE DIAGNOSTICS
# ============================================================

# These diagnostics distinguish the theoretical retirement age
# T = 67 from the age-70 endpoint of the supplied P+ glide path.

retirement_summary <- thepensiondatainclgamma %>%
  summarise(
    
    N_total = n(),
    
    N_pre_retirement =
      sum(
        alder <
          retirement_age
      ),
    
    N_at_or_after_retirement =
      sum(
        alder >=
          retirement_age
      ),
    
    N_age_50_to_66 =
      sum(
        alder >= 50 &
          alder < 67
      ),
    
    N_age_67_to_69 =
      sum(
        alder >= 67 &
          alder < 70
      ),
    
    N_age_70_plus =
      sum(
        alder >= 70
      )
  )


retirement_summary


# Same diagnostic for the final five-premium sample.

final_retirement_summary <- gamma_5yr %>%
  summarise(
    
    N_total = n(),
    
    N_pre_retirement =
      sum(
        alder <
          retirement_age
      ),
    
    N_at_or_after_retirement =
      sum(
        alder >=
          retirement_age
      ),
    
    N_age_50_to_66 =
      sum(
        alder >= 50 &
          alder < 67
      ),
    
    N_age_67_to_69 =
      sum(
        alder >= 67 &
          alder < 70
      ),
    
    N_age_70_plus =
      sum(
        alder >= 70
      )
  )


final_retirement_summary


# ============================================================
# 25. FINAL SAMPLE CHECKS
# ============================================================

# Confirm that every observation in the final sample satisfies
# the five-annual-premium restriction.

final_restriction_check <- subset(
  gamma_5yr,
  !is.finite(annual_premium) |
    depot <=
    5 *
    annual_premium
)


if (
  nrow(
    final_restriction_check
  ) > 0
) {
  
  warning(
    "At least one observation in the final sample fails the five-annual-premium restriction."
  )
}


# Check gamma values in the final sample.

final_non_finite_gamma <- subset(
  gamma_5yr,
  !is.finite(
    gamma
  )
)


if (
  nrow(
    final_non_finite_gamma
  ) > 0
) {
  
  warning(
    paste(
      nrow(
        final_non_finite_gamma
      ),
      "observations in the final sample have non-finite gamma."
    )
  )
}

# ============================================================
# 26. CROSS-SECTIONAL IMPLIED RISK AVERSION BY OBSERVED AGE
# ============================================================
#
# Research Question 3 asks whether the estimated values of
# relative risk aversion are consistent over time.
#
# The available P+ data are cross-sectional rather than
# longitudinal. The same investor is not observed repeatedly
# through time. A true within-investor time-series test of
# stability is therefore not possible from these data.
#
# The analysis below instead examines whether the calibrated
# gamma estimates differ systematically across the observed
# ages in the final five-annual-premium sample.
#
# No new calibration equation is introduced. Each plotted
# gamma is the current-age estimate already calculated in
# Section 12:
#
#   gamma_hat_i =
#
#     (alpha_theta,i - r)
#     ---------------------
#       sigma_theta,i^2
#
#     * (X_i + h_i) / X_i.
#
# The analysis is descriptive and should not be interpreted as
# the evolution of an individual investor's risk aversion over
# time.
# ============================================================


gamma_age_sample <- gamma_5yr %>%
  filter(
    is.finite(gamma),
    is.finite(alder),
    alder >= 50,
    alder < retirement_age
  )


# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

gamma_age_summary <- gamma_age_sample %>%
  summarise(
    N = n(),
    Mean_Age = mean(alder),
    Min_Age = min(alder),
    Max_Age = max(alder),
    Mean_Gamma = mean(gamma),
    Median_Gamma = median(gamma),
    SD_Gamma = sd(gamma),
    Min_Gamma = min(gamma),
    Q25_Gamma = as.numeric(
      quantile(
        gamma,
        0.25
      )
    ),
    Q75_Gamma = as.numeric(
      quantile(
        gamma,
        0.75
      )
    ),
    Max_Gamma = max(gamma)
  )


gamma_age_summary


# ------------------------------------------------------------
# Observed gamma against observed age
# ------------------------------------------------------------

d_gamma_by_age <- ggplot(
  gamma_age_sample,
  aes(
    x = alder,
    y = gamma
  )
) +
  geom_point(
    alpha = 0.70,
    size = 2
  ) +
  labs(
    x = "Observed age",
    y = TeX("Implied $\\gamma$"),
    title = "Implied Risk Aversion by Observed Age"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(
      50,
      retirement_age - 1,
      by = 2
    )
  )


d_gamma_by_age


ggsave(
  filename = file.path(
    figure_dir,
    "Implied Risk Aversion by Observed Age.png"
  ),
  plot = d_gamma_by_age,
  width = 8,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# Save cross-sectional age output
# ------------------------------------------------------------

write.csv(
  gamma_age_sample,
  file.path(
    results_dir,
    "calibration_gamma_by_observed_age.csv"
  ),
  row.names = FALSE
)

write.csv(
  gamma_age_summary,
  file.path(
    results_dir,
    "calibration_gamma_observed_age_summary.csv"
  ),
  row.names = FALSE
)
# ============================================================
# 26A. IMPLIED RISK AVERSION AT AND AFTER RETIREMENT
# ============================================================
#
# The direct inverse investment calibration is also informative
# at and after the theoretical retirement age.
#
# For age >= T = 67:
#
#   h_i = 0,
#
# so
#
#   gamma_hat_i =
#     (alpha_theta,i - r) / sigma_theta,i^2.
#
# Between ages 67 and 70, the P+ portfolio continues along the
# supplied de-risking path.
#
# From age 70 onwards, the P+ portfolio is held constant at its
# age-70 allocation. Therefore, under fixed market parameters,
# observations in the same P+ profile aged 70+ should imply the
# same gamma apart from numerical precision.
# ============================================================


gamma_age_extended <- thepensiondatainclgamma %>%
  filter(
    is.finite(gamma),
    is.finite(alder),
    alder >= 50
  ) %>%
  mutate(
    Age_Group = case_when(
      alder < retirement_age ~ "50 to <67",
      alder < 70 ~ "67 to <70",
      TRUE ~ "70+"
    )
  )


gamma_age_extended$Age_Group <- factor(
  gamma_age_extended$Age_Group,
  levels = c(
    "50 to <67",
    "67 to <70",
    "70+"
  )
)


# ------------------------------------------------------------
# Summary by age group
# ------------------------------------------------------------

gamma_age_group_summary <- gamma_age_extended %>%
  group_by(
    Age_Group
  ) %>%
  summarise(
    N = n(),
    Mean_Age = mean(alder),
    Min_Age = min(alder),
    Max_Age = max(alder),
    Mean_Gamma = mean(gamma),
    Median_Gamma = median(gamma),
    SD_Gamma = ifelse(
      n() > 1,
      sd(gamma),
      NA_real_
    ),
    Min_Gamma = min(gamma),
    Q25_Gamma = as.numeric(
      quantile(
        gamma,
        0.25
      )
    ),
    Q75_Gamma = as.numeric(
      quantile(
        gamma,
        0.75
      )
    ),
    Max_Gamma = max(gamma),
    .groups = "drop"
  )


gamma_age_group_summary


# ------------------------------------------------------------
# Age 70+ observations
# ------------------------------------------------------------

gamma_70_plus <- gamma_age_extended %>%
  filter(
    alder >= 70
  )


gamma_70_plus


# ------------------------------------------------------------
# Age 70+ summary by P+ profile
# ------------------------------------------------------------

gamma_70_plus_profile_summary <- gamma_70_plus %>%
  group_by(
    investeringsProfil
  ) %>%
  summarise(
    N = n(),
    Mean_Age = mean(alder),
    Min_Age = min(alder),
    Max_Age = max(alder),
    Mean_Gamma = mean(gamma),
    Median_Gamma = median(gamma),
    SD_Gamma = ifelse(
      n() > 1,
      sd(gamma),
      NA_real_
    ),
    Min_Gamma = min(gamma),
    Max_Gamma = max(gamma),
    .groups = "drop"
  )


gamma_70_plus_profile_summary


# ------------------------------------------------------------
# Check h = 0 for all age 70+ observations
# ------------------------------------------------------------

gamma_70_plus_nonzero_h <- gamma_70_plus %>%
  filter(
    !is.finite(future_premiums) |
      abs(future_premiums) > 1e-10
  )


if (
  nrow(
    gamma_70_plus_nonzero_h
  ) > 0
) {
  
  warning(
    "At least one age-70+ observation has non-zero future premiums."
  )
}


# ------------------------------------------------------------
# Check gamma is constant within profile after age 70
# ------------------------------------------------------------

gamma_70_plus_constancy_check <- gamma_70_plus %>%
  group_by(
    investeringsProfil
  ) %>%
  summarise(
    N = n(),
    Min_Gamma = min(gamma),
    Max_Gamma = max(gamma),
    Gamma_Range =
      Max_Gamma -
      Min_Gamma,
    .groups = "drop"
  )


gamma_70_plus_constancy_check


# ------------------------------------------------------------
# Plot gamma across the full direct-calibration age range
# ------------------------------------------------------------

d_gamma_by_age_extended <- ggplot(
  gamma_age_extended,
  aes(
    x = alder,
    y = gamma
  )
) +
  geom_point(
    alpha = 0.70,
    size = 2
  ) +
  geom_vline(
    xintercept = retirement_age,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_vline(
    xintercept = 70,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  labs(
    x = "Observed age",
    y = TeX("Implied $\\gamma$"),
    title = "Implied Risk Aversion by Observed Age"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = scales::pretty_breaks(
      n = 8
    )
  )


d_gamma_by_age_extended


ggsave(
  filename = file.path(
    figure_dir,
    "Implied Risk Aversion by Observed Age - Extended.png"
  ),
  plot = d_gamma_by_age_extended,
  width = 8,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# Save age-group outputs
# ------------------------------------------------------------

write.csv(
  gamma_age_group_summary,
  file.path(
    results_dir,
    "calibration_gamma_summary_by_age_group.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_70_plus,
  file.path(
    results_dir,
    "calibration_gamma_age_70_plus.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_70_plus_profile_summary,
  file.path(
    results_dir,
    "calibration_gamma_age_70_plus_by_profile.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_70_plus_constancy_check,
  file.path(
    results_dir,
    "calibration_gamma_age_70_plus_constancy_check.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 26B. AGE GROUPS IN THE FIVE-PREMIUM SAMPLE
# ============================================================

gamma_5yr_age_extended <- gamma_5yr %>%
  filter(
    is.finite(gamma),
    is.finite(alder),
    alder >= 50
  ) %>%
  mutate(
    Age_Group = case_when(
      alder < retirement_age ~ "50 to <67",
      alder < 70 ~ "67 to <70",
      TRUE ~ "70+"
    )
  )


gamma_5yr_age_extended$Age_Group <- factor(
  gamma_5yr_age_extended$Age_Group,
  levels = c(
    "50 to <67",
    "67 to <70",
    "70+"
  )
)


gamma_5yr_age_group_summary <- gamma_5yr_age_extended %>%
  group_by(
    Age_Group
  ) %>%
  summarise(
    N = n(),
    Mean_Age = mean(alder),
    Min_Age = min(alder),
    Max_Age = max(alder),
    Mean_Gamma = mean(gamma),
    Median_Gamma = median(gamma),
    SD_Gamma = ifelse(
      n() > 1,
      sd(gamma),
      NA_real_
    ),
    Min_Gamma = min(gamma),
    Q25_Gamma = as.numeric(
      quantile(
        gamma,
        0.25
      )
    ),
    Q75_Gamma = as.numeric(
      quantile(
        gamma,
        0.75
      )
    ),
    Max_Gamma = max(gamma),
    .groups = "drop"
  )


gamma_5yr_age_group_summary


gamma_5yr_70_plus <- gamma_5yr_age_extended %>%
  filter(
    alder >= 70
  )


gamma_5yr_70_plus_profile_summary <- gamma_5yr_70_plus %>%
  group_by(
    investeringsProfil
  ) %>%
  summarise(
    N = n(),
    Mean_Age = mean(alder),
    Min_Age = min(alder),
    Max_Age = max(alder),
    Mean_Gamma = mean(gamma),
    Median_Gamma = median(gamma),
    SD_Gamma = ifelse(
      n() > 1,
      sd(gamma),
      NA_real_
    ),
    Min_Gamma = min(gamma),
    Max_Gamma = max(gamma),
    .groups = "drop"
  )


gamma_5yr_70_plus_profile_summary


write.csv(
  gamma_5yr_age_group_summary,
  file.path(
    results_dir,
    "calibration_final_gamma_summary_by_age_group.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_5yr_70_plus,
  file.path(
    results_dir,
    "calibration_final_gamma_age_70_plus.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_5yr_70_plus_profile_summary,
  file.path(
    results_dir,
    "calibration_final_gamma_age_70_plus_by_profile.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 27. SAVE NUMERICAL OUTPUTS
# ============================================================

write.csv(
  restriction_summary,
  file.path(
    results_dir,
    "calibration_restriction_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_restriction_summary,
  file.path(
    results_dir,
    "calibration_gamma_summary_by_restriction.csv"
  ),
  row.names = FALSE
)


write.csv(
  as.data.frame(
    gamma_profile_summary
  ),
  file.path(
    results_dir,
    "calibration_gamma_summary_by_profile.csv"
  ),
  row.names = FALSE
)


write.csv(
  gamma_5yr,
  file.path(
    results_dir,
    "calibration_final_sample.csv"
  ),
  row.names = FALSE
)


write.csv(
  retirement_summary,
  file.path(
    results_dir,
    "calibration_retirement_age_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  final_retirement_summary,
  file.path(
    results_dir,
    "calibration_final_retirement_age_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  thepensiondatainclgamma,
  file.path(
    results_dir,
    "calibration_full_sample_with_gamma.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. SAVE KEY RESULTS
# ============================================================

capture.output({
  
  cat(
    "CHAPTER 5 CALIBRATION RESULTS\n"
  )
  
  cat(
    "=============================\n\n"
  )
  
  
  cat(
    "Retirement age:",
    retirement_age,
    "\n"
  )
  
  
  cat(
    "Calibration risk-free rate:",
    calibration_r,
    "\n\n"
  )
  
  
  cat(
    "Calibration timing:",
    "current recorded age\n\n"
  )
  
  
  cat(
    "Zero-premium observations:",
    zero_premium_count,
    "\n"
  )
  
  
  cat(
    "Missing/non-finite premium observations:",
    missing_premium_count,
    "\n\n"
  )
  
  
  cat(
    "Broad calibration sample:\n"
  )
  
  print(
    retirement_summary,
    row.names = FALSE
  )
  
  
  cat(
    "\nFinal five-premium sample:\n"
  )
  
  print(
    final_retirement_summary,
    row.names = FALSE
  )
  
  
  cat(
    "\nSample restrictions:\n"
  )
  
  print(
    restriction_summary,
    row.names = FALSE
  )
  
  
  cat(
    "\nGamma statistics by restriction:\n"
  )
  
  print(
    gamma_restriction_summary,
    row.names = FALSE
  )
  
  
  cat(
    "\nFinal gamma statistics by profile:\n"
  )
  
  print(
    as.data.frame(
      gamma_profile_summary
    ),
    row.names = FALSE
  )
  
  
  cat(
    "\nCross-sectional variation by observed age:\n"
  )
  
  cat(
    "The available data are cross-sectional, so longitudinal ",
    "within-investor stability cannot be tested directly.\n"
  )
  
  print(
    as.data.frame(
      gamma_age_summary
    ),
    row.names = FALSE
  )
  
  
  
  cat(
    "\nInternal checks:\n"
  )
  
  cat(
    "Post-retirement observations with non-zero h:",
    nrow(
      post_retirement_premium_check
    ),
    "\n"
  )
  
  cat(
    "Observations with negative h:",
    nrow(
      negative_future_premiums
    ),
    "\n"
  )
  
  cat(
    "Non-finite gamma values:",
    nrow(
      non_finite_gamma
    ),
    "\n"
  )
  
  cat(
    "Non-positive gamma values:",
    nrow(
      non_positive_gamma
    ),
    "\n"
  )
  
},
file = file.path(
  results_dir,
  "calibration_key_results.txt"
))


# ============================================================
# 28. FINAL CONSOLE OUTPUT
# ============================================================

cat(
  "\nChapter 5 calibration complete.\n",
  "Retirement age: ",
  retirement_age,
  "\n",
  "Calibration timing: current recorded age\n",
  "Broad sample size: ",
  nrow(
    thepensiondatainclgamma
  ),
  "\n",
  "Final five-premium sample size: ",
  nrow(
    gamma_5yr
  ),
  "\n",
  "Updated figures are in '",
  figure_dir,
  "'.\n",
  "Numerical outputs are in '",
  results_dir,
  "'.\n",
  sep = ""
)