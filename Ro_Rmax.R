library(ggplot2)

rm(list = ls())
gc()

Tmax_dev <- 45

# 1. Model development Definitions

# A - development Rates Egg (Brière-1)
dev_rate_egg <- function(T, a = 0.0001004294, T0 = 5, Tmax = Tmax_dev){
  rate <- a * T * (T - T0) * sqrt(pmax(0, Tmax - T))
  return(ifelse(T >= T0 & T <= Tmax, rate, 0))
}

# B - development Rates larva (Brière-2)
dev_rate_larva <- function(T) {
  a    <- 0.0006138710
  T0   <- 5
  Tmax <- Tmax_dev
  m    <- 7.358106
  
  rate <- a * T * (T - T0) * (pmax(0, Tmax - T))^(1 / m)
  return(ifelse(T >= T0 & T <= Tmax, rate, 0))
}

# C - development Rates Nymphs (quadratic)
dev_rate_nymphs <- function(T) {
  beta0 <- 0.4985068
  beta1 <- -0.02980347
  beta2 <- 0.0007908584
  
  rate <- beta0 + (beta1 * T) + (beta2 * (T^2))
  return(pmin(pmax(rate, 0.01), 0.999))
}

# 2. Survival model

# A - Egg Hatchability Function (Quadratic Fit)
egg_hatch_T <- function(T) {
  pmin(
    pmax(
      0.8216286 + 0.01255143 * T - 0.0002314286 * T^2,
      0.01
    ),
    0.999
  )
}

# 3. Mortality model

# A - Egg mort
mort_egg_T <- function(T) {
  -log(egg_hatch_T(T)) * dev_rate_egg(T)  # DAILY egg mortality rate
}

# B - larva and Nymphs mort
compound_stage_partition <- function(T) {
  raw_surv <- 0.7334286 + 0.02027143 * T - 0.0004314286 * (T^2)
  tot_surv <- pmin(pmax(raw_surv, 0.01), 0.999)
  
  dur_larva <- 1.85   # Larva block (larva + q1)
  dur_nymph <- 4.225  # Nymph block ((protonymph + q2) + (deutonymph + q3))
  dur_total <- dur_larva + dur_nymph
  
  w_larva <- dur_larva / dur_total  
  w_nymph <- dur_nymph / dur_total  
  
  surv_larva <- tot_surv ^ w_larva
  surv_nymph <- tot_surv ^ w_nymph
  
  data.frame(
    Temperature     = T,
    Total_Survival  = tot_surv,
    Larva_Survival  = surv_larva,
    Nymph_Survival  = surv_nymph,
    Larva_Mortality = 1 - surv_larva,
    Nymph_Mortality = 1 - surv_nymph
  )
}

# C - Adult mort
mort_adult_T <- function(T) {
  pmin(
    pmax(
      -0.03134242 + 0.003521138 * T,
      0.001
    ),
    0.999
  )
}

# 4. Fecundity & Sex Ratio (Ratkowsky-2 model)
fecundity_adult_T <- function(T) {
  T0   <- 15.000000
  b    <- 1.7
  c    <- 0.1
  Tmax <- 35
  
  rate <- b * (T - T0) * (1 - exp(c * (T - Tmax)))
  return(ifelse(T >= T0 & T <= Tmax, pmax(rate, 0), 0))
}

sex_ratio <- c(0.831, 0.871, 0.872, 0.849, 0.865)
sr_mean <- mean(sex_ratio)

# =========================================================
# Plant Resource Capacity R_max(T) Function (Gent, 1986)
# =========================================================
get_plant_Rmax <- function(T) {
  I <- 370 # High irradiance level from Gent (1986)
  rgr <- -0.2949 + 0.000329 * I + 0.04001 * T - 0.000797 * (T^2)
  return(pmin(pmax(rgr, 0.01), 1.0))
}

# Attack / Consumption Rate alpha(T) - Brière-2 Model
alpha_T <- function(T) {
  a  <- 9.939369e-06
  T0 <- 10.39128
  m  <- 1.0
  Tmax <- 38.0
  if (is.na(T) || T < T0 || T > Tmax) return(0.001)
  rate <- a * T * (T - T0) * (Tmax - T)^(1/m)
  return(max(rate, 0.001))
}

# Conversion Efficiency epsilon(T) - Brière-2 Model
epsilon_T <- function(T) {
  a  <- 3.605811e-05
  T0 <- 5.216056
  m  <- 1.0
  Tmax <- 38.0
  if (is.na(T) || T < T0 || T > Tmax) return(0.01)
  rate <- a * T * (T - T0) * (Tmax - T)^(1/m)
  return(max(rate, 0.01))
}

# =========================================================
# Manuscript-Based Resource-Coupled R0 Function
# with Constant Test Values for Alpha and Epsilon
# =========================================================
compute_suitability_manuscript <- function(T) {
  if (is.na(T) || T < 11 || T > 45) return(NA_real_)
  
  d1 <- dev_rate_egg(T)
  d2 <- dev_rate_larva(T)
  d3 <- dev_rate_nymphs(T) 
  
  mu_e <- mort_egg_T(T)
  
  immature_mort <- compound_stage_partition(T)
  mu_l <- immature_mort$Larva_Mortality
  mu_n <- immature_mort$Nymph_Mortality
  
  mu_a <- mort_adult_T(T) 
  
  phi <- fecundity_adult_T(T) * sr_mean
  
  # Get Temperature-dependent R_max from Gent (1986)
  R_max <- get_plant_Rmax(T)
  
  # Constant test values for Alpha and Epsilon (as requested)
  alpha   <- alpha_T(T)  # Constant consumption/attack rate test value
  epsilon <- epsilon_T(T)  # Constant conversion efficiency test value
  
  # Manuscript Denominator: (mu_4 - epsilon * alpha * R_max)
  net_adult_loss <- mu_a - (epsilon * alpha * R_max)
  
  if (!is.finite(net_adult_loss) || net_adult_loss <= 0) {
    net_adult_loss <- 0.0001 
  }
  
  numerator   <- phi * d1 * d2 * d3
  denominator <- (d1 + mu_e) * (d2 + mu_l) * (d3 + mu_n) * net_adult_loss
  
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  
  return(numerator / denominator)
}

# =========================================================
# Compute and Plot Manuscript R_0 Curve
# =========================================================
temps <- seq(11, 45, by = 0.1)
Ro_values <- sapply(temps, compute_suitability_manuscript)

plot_data <- data.frame(
  temp = temps,
  Ro   = Ro_values
)

plot_data <- na.omit(plot_data)

p_thermal_manuscript <- ggplot(plot_data, aes(x = temp, y = Ro)) +
  geom_line(color = "#084594", linewidth = 1.3) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.9) +
  labs(x = "Temperature (°C)", y = expression(R[0] ~ "(Manuscript Model)")) +
  theme_classic(base_size = 14)

print(p_thermal_manuscript)