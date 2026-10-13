# ==============================================================================
# Probability of extinction using BN inference. Real data of Abakaliki smallpox 
# outbreak
# ==============================================================================

if (!require("outbreaks")) install.packages("outbreaks")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")

library(outbreaks)
library(ggplot2)
library(patchwork)

data("smallpox_abakaliki_1967")
fechas <- smallpox_abakaliki_1967$date_of_onset
fechas <- fechas[!is.na(fechas)]

# Interval grouping (~14 days)
dias_inicio <- as.numeric(fechas - min(fechas))
generacion_id <- floor(dias_inicio / 14) + 1

# Population (Z_0, Z_1, ..., Z_k)
Z <- as.numeric(table(factor(generacion_id, levels = 1:max(generacion_id))))
N_gen <- length(Z)

# ------------------------------------------------------------------------------
# FORWARD SAMPLING OVER THE BAYESIAN NETWORK
# ------------------------------------------------------------------------------
# FUNCTION FOR PROPAGATING THE TRAJECTORIES Z_{k+1}, ..., Z_{k+m}
propagar_extincion_dbn <- function(Z_k, alpha_k, beta_k, horizonte_m = 3, 
                                   n_lambda = 1000, n_trayectorias = 200) {
  
  # If the population has extinguished, the future probability of extinction is 1
  if (Z_k == 0) {
    return(c(media = 1, ic_inf = 1, ic_sup = 1, map = 1))
  }
  
  # Sample the parameter lambda ~ Gamma(alpha_k, beta_k)
  lambda_samples <- rgamma(n_lambda, shape = alpha_k, rate = beta_k)
  p_ext_por_lambda <- numeric(n_lambda)
  
  # Propagation in the BN
  for (s in 1:n_lambda) {
    l <- lambda_samples[s]
    # Initial parallel trajectories
    z_estado <- rep(Z_k, n_trayectorias)
    
    # Transitions Z_{j} ~ Poisson(lambda * Z_{j-1})
    for (step in 1:horizonte_m) {
      z_estado <- rpois(n_trayectorias, lambda = l * z_estado)
    }
    
    # Proportion of trajectories extinguishes at horizon m
    p_ext_por_lambda[s] <- mean(z_estado == 0)
  }
  
  # Plug-in MAP estimator
  lambda_map <- alpha_k / beta_k
  z_map <- rep(Z_k, 10000)
  for (step in 1:horizonte_m) {
    z_map <- rpois(10000, lambda = lambda_map * z_map)
  }
  p_ext_map <- mean(z_map == 0)
  
  return(c(
    media  = mean(p_ext_por_lambda),
    ic_inf = as.numeric(quantile(p_ext_por_lambda, 0.025)),
    ic_sup = as.numeric(quantile(p_ext_por_lambda, 0.975)),
    map    = p_ext_map
  ))
}

# ------------------------------------------------------------------------------
# MAIN LOOP. Inference and extinction ovaluation at m = 3 steps ahead
# ------------------------------------------------------------------------------
horizonte_n <- 3 
alpha_0 <- 1
beta_0  <- 1

sum_Z_obs  <- 0
sum_Z_prev <- 0

matriz_extincion <- matrix(0, nrow = N_gen - 1, ncol = 4)
colnames(matriz_extincion) <- c("Media_Posterior", "IC_95_Inf", "IC_95_Sup", "Estimador_MAP")

set.seed(123) 

for (k in 1:(N_gen - 1)) {
  sum_Z_obs  <- sum_Z_obs + Z[k + 1]
  sum_Z_prev <- sum_Z_prev + Z[k]
  
  alpha_k <- alpha_0 + sum_Z_obs
  beta_k  <- beta_0 + sum_Z_prev
  
  Z_actual <- Z[k + 1]
  
  # Propagation on the BN
  matriz_extincion[k, ] <- propagar_extincion_dbn(
    Z_k = Z_actual, 
    alpha_k = alpha_k, 
    beta_k = beta_k, 
    horizonte_m = horizonte_n
  )
}

df_res <- data.frame(
  Gen_Actual = 1:(N_gen - 1),
  Poblacion_Zk = Z[2:N_gen],
  Prob_Extincion_n3 = matriz_extincion[, "Media_Posterior"],
  IC_95_Inf = matriz_extincion[, "IC_95_Inf"],
  IC_95_Sup = matriz_extincion[, "IC_95_Sup"],
  Estimado_MAP = matriz_extincion[, "Estimador_MAP"]
)

# ------------------------------------------------------------------------------
#  ggplot2 visualization
# ------------------------------------------------------------------------------

# Panel A: Cases
p1 <- ggplot(df_res, aes(x = Gen_Actual, y = Poblacion_Zk)) +
  geom_line(color = "#08306b", linewidth = 0.9) +
  geom_point(color = "#08519c", size = 2.5) +
  scale_x_continuous(breaks = df_res$Gen_Actual) +
  labs(
    title = "",
    x = "Generation (k)",
    y = expression("Active Cases (" * Z[k] * ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Panel B: Extinction probability
p2 <- ggplot(df_res, aes(x = Gen_Actual)) +
  geom_ribbon(aes(ymin = IC_95_Inf, ymax = IC_95_Sup, fill = "95% CI"), alpha = 0.25) +
  geom_line(aes(y = Prob_Extincion_n3, color = "Predictive Distribution"), linewidth = 0.9) +
  geom_point(aes(y = Prob_Extincion_n3, color = "Predictive Distribution"), size = 2.5) +
  geom_line(aes(y = Estimado_MAP, color = "MAP Plug-in Estimator"), linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "firebrick", linewidth = 0.8) +
  scale_x_continuous(breaks = df_res$Gen_Actual) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Predictive Distribution" = "#006d2c",
      "MAP Plug-in Estimator" = "#e6550d"
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("95% CI" = "#74c476")
  ) +
  labs(
    title = "",
    x = "Generation (k)",
    y = expression("P(" * Z[k + 3] * " = 0 | " * Z[k] * ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  )

figura_extincion_dbn <- p1 | p2
print(figura_extincion_dbn)
