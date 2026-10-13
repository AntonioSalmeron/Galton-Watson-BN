# ==============================================================================
# Real data: Ebola outbreak
# ==============================================================================

if (!require("outbreaks")) install.packages("outbreaks")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")

library(outbreaks)
library(ggplot2)
library(patchwork)

data("ebola_sim")
linelist <- ebola_sim$linelist

# Dates where symptoms started
fechas <- linelist$date_of_onset
fechas <- fechas[!is.na(fechas)]

# Grouping by infection generation (Serial interval ~ 14 days)
dias_inicio <- as.numeric(fechas - min(fechas))
generacion_id <- floor(dias_inicio / 14) + 1

# Individuals per generation (Z_0, Z_1, Z_2, ...)
Z <- as.numeric(table(factor(generacion_id, levels = 1:max(generacion_id))))
N_gen <- length(Z)

cat("--- OBSERVED INDIVIDUAL BY GENERATION (Z_i) ---\n")
print(data.frame(Generacion = 0:(N_gen - 1), Casos = Z))

# BAYESIAN NETWORK INFERENCE
alpha_0 <- 1
beta_0  <- 1

alpha_post  <- numeric(N_gen - 1)
beta_post   <- numeric(N_gen - 1)
lambda_mean <- numeric(N_gen - 1)
ic_lower    <- numeric(N_gen - 1)
ic_upper    <- numeric(N_gen - 1)

sum_Z_obs  <- 0 
sum_Z_prev <- 0 

for (k in 1:(N_gen - 1)) {
  sum_Z_obs  <- sum_Z_obs + Z[k + 1]
  sum_Z_prev <- sum_Z_prev + Z[k]    
  
  alpha_post[k] <- alpha_0 + sum_Z_obs
  beta_post[k]  <- beta_0 + sum_Z_prev
  
  lambda_mean[k] <- alpha_post[k] / beta_post[k]
  ic_lower[k]    <- qgamma(0.025, shape = alpha_post[k], rate = beta_post[k])
  ic_upper[k]    <- qgamma(0.975, shape = alpha_post[k], rate = beta_post[k])
}

# Data frames for ggplot2
df_casos <- data.frame(
  Generacion = 0:(N_gen - 1),
  Casos = Z
)

df_lambda <- data.frame(
  Gen = 1:(N_gen - 1),
  Lambda_MQR = lambda_mean,
  IC_Inf = ic_lower,
  IC_Sup = ic_upper
)

# Plots generation for ggplot2

# Panel A: Curve by generation
p1 <- ggplot(df_casos, aes(x = Generacion, y = Casos)) +
  geom_line(color = "#a50f15", linewidth = 0.9) +
  geom_point(color = "#67000d", size = 2) +
  labs(
    title = "",
    x = "Generation (14-days intervals)",
    y = expression("No. of cases (" * Z[i] * ")")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Panel B: Bayesian network inference over Lambda
p2 <- ggplot(df_lambda, aes(x = Gen, y = Lambda_MQR)) +
  geom_ribbon(aes(ymin = IC_Inf, ymax = IC_Sup), fill = "#6baed6", alpha = 0.35) +
  geom_line(color = "#08519c", linewidth = 0.9) +
  geom_point(color = "#08306b", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  annotate("text", x = N_gen - 6, y = 1.08, label = "", 
           color = "firebrick", fontface = "italic", size = 3.5) +
  labs(
    title = "",
    x = "Generation",
    y = expression("Estimation " * hat(lambda)[MQR] ~ "(E[" * lambda * " | Z])")
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Figure composition
figura_final <- p1 | p2

# Output
print(figura_final)
