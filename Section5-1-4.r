# ==============================================================================
# GW with immigration and evolving parameter
# ==============================================================================

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")

library(ggplot2)
library(patchwork)

# GW process with immigration simulation
set.seed(42)

n_gen <- 10
lambda_real <- 0.85 # Sub-critical regime
z <- numeric(n_gen)
z_0 <- 5            # Initial population

# Generaction 1
z[1] <- rpois(1, lambda = z_0 * lambda_real)

# Generations 2 to 4: Normal evolution
for (i in 2:4) {
  z[i] <- rpois(1, lambda = z[i-1] * lambda_real)
}

# Massive immigration at generation 5
z[5] <- rpois(1, lambda = z[4] * lambda_real) + 30 

# Generations 6 to 10: The process goes on
for (i in 6:n_gen) {
  z[i] <- rpois(1, lambda = z[i-1] * lambda_real)
}

# Bayesian estimation (Theorem 7)
lambda_est <- numeric(n_gen)
alpha_prior <- 1
beta_prior  <- 10

# Initial inference
lambda_est[1] <- (alpha_prior + z[1]) / (beta_prior + z_0)

# Iterative update (Gen 2 to 10)
for (j in 2:n_gen) {
  alpha_post <- z[j] + 2 * lambda_est[j-1]
  beta_post  <- z[j-1] + 2
  lambda_est[j] <- alpha_post / beta_post # MQR estimation
}

# ggplot2 panel construction
df_sim <- data.frame(
  Generacion = 1:n_gen,
  Casos_Z = z,
  Lambda_Est = lambda_est,
  Lambda_Real = lambda_real
)

# ==============================================================================
# 5. VISUALIZATION 
# ==============================================================================

# Panel A: Time series
p1 <- ggplot(df_sim, aes(x = Generacion, y = Casos_Z)) +
  geom_line(color = "#2b5c8f", linewidth = 0.9) +
  geom_point(color = "#1a365d", size = 2.5) +
  geom_point(data = subset(df_sim, Generacion == 5), 
             aes(x = Generacion, y = Casos_Z), 
             color = "#c0392b", size = 4) +
  annotate("text", x = 5.2, y = max(z) - 2, 
           label = "Immigration (+30 individals)", 
           color = "#c0392b", fontface = "bold", size = 3.3, hjust = 0) +
  scale_x_continuous(breaks = 1:n_gen) +
  labs(
    title = "",
    x = "Generation (j)",
    y = expression(Z[j])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Panel B: Adaptation of the Lambda parameter
p2 <- ggplot(df_sim, aes(x = Generacion)) +
  # Real value of Lambda
  geom_hline(aes(yintercept = lambda_real, color = "True value"), 
             linetype = "dashed", linewidth = 0.8) +
  # MQR estimates
  geom_line(aes(y = Lambda_Est, color = "MQR estimation"), 
            linewidth = 0.9) +
  geom_point(aes(y = Lambda_Est, color = "MQR estimation"), 
             size = 2.5) +
  # Vertical line where the immigration occurs
  geom_vline(xintercept = 5, linetype = "dotted", color = "#c0392b", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:n_gen) +
  scale_y_continuous(limits = c(0, max(lambda_est) + 0.3)) +
  scale_color_manual(
    name = NULL,
    values = c(
      "MQR estimation" = "#27ae60",
      "True value" = "#2c3e50"
    )
  ) +
  labs(
    title = "",
    x = "Generation (j)",
    y = expression(hat(lambda)[j])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  )

# multi-panel figure
figura_inmigracion <- p1 | p2

# Output
print(figura_inmigracion)


