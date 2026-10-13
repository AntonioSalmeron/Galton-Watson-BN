# ==============================================================================
# Missing data habdling
# ==============================================================================

if (!require("rjags")) install.packages("rjags")
if (!require("coda")) install.packages("coda")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")

library(rjags)
library(coda)
library(ggplot2)
library(patchwork)

# 1. Data simulation (Ground Truth)
set.seed(42)
n_gen <- 10
lambda_real <- 1.15
z_clean <- numeric(n_gen)
z_0 <- 5 
z_clean[1] <- rpois(1, lambda = z_0 * lambda_real)

for(i in 2:n_gen) {
  z_clean[i] <- rpois(1, lambda = z_clean[i-1] * lambda_real)
}

# 2. Including missing data (generations 4 and 5)
z_obs <- z_clean
z_obs[4] <- NA 
z_obs[5] <- NA 

# 3. Model specification in BUGS
gw_model_code <- "
model {
  lambda ~ dgamma(alpha, beta)
  Z[1] ~ dpois(Z0 * lambda)
  for (i in 2:N) {
    Z[i] ~ dpois(Z[i-1] * lambda)
  }
}
"
writeLines(gw_model_code, "gw_model.jags")

# 4. Model compilation and MCMC inference
data_jags <- list(Z = z_obs, N = n_gen, Z0 = z_0, alpha = 1, beta = 10)
z_inits <- rep(NA, n_gen)
z_inits[4] <- 10
z_inits[5] <- 10
inits_list <- list(Z = z_inits, lambda = 1)

jags_model <- jags.model("gw_model.jags", data = data_jags, inits = inits_list, n.chains = 3, quiet = TRUE)
update(jags_model, 1000)
samples <- coda.samples(jags_model, variable.names = c("Z[4]", "Z[5]", "lambda"), n.iter = 5000)

# 5. MCMC sample extraction
mcmc_mat <- as.matrix(samples)
z4_samples <- mcmc_mat[, "Z[4]"]
z5_samples <- mcmc_mat[, "Z[5]"]

z4_mean <- mean(z4_samples)
z4_ci   <- quantile(z4_samples, probs = c(0.025, 0.975))

z5_mean <- mean(z5_samples)
z5_ci   <- quantile(z5_samples, probs = c(0.025, 0.975))

# ==============================================================================
# 6. Data frame construction for ggplot2
# ==============================================================================

# Data frame with t he complete trajectory (Ground Truth)
df_truth <- data.frame(
  Generacion = 1:n_gen,
  Z_real = z_clean
)

# Data frame for observed blocks 
df_obs1 <- data.frame(Generacion = 1:3, Z_obs = z_clean[1:3])
df_obs2 <- data.frame(Generacion = 6:10, Z_obs = z_clean[6:10])

# Data frame for the imputation of the latent variables (Z[4] and Z[5])
df_imp <- data.frame(
  Generacion = c(4, 5),
  Media = c(z4_mean, z5_mean),
  IC_Inf = c(z4_ci[1], z5_ci[1]),
  IC_Sup = c(z4_ci[2], z5_ci[2])
)

# Data frame para for the posterior MCMC densities
df_mcmc <- rbind(
  data.frame(Variable = "Z[4]", Muestra = z4_samples),
  data.frame(Variable = "Z[5]", Muestra = z5_samples)
)

df_vlines <- data.frame(
  Variable = c("Z[4]", "Z[5]"),
  Real = c(z_clean[4], z_clean[5])
)

# ==============================================================================
# 7. VISUALIZATION 
# ==============================================================================

# Panel A: Time series and imputation
p1 <- ggplot() +
  # Unknown trajectory
  geom_line(data = df_truth, aes(x = Generacion, y = Z_real, linetype = "Unknown true values"), 
            color = "gray50", linewidth = 0.8) +
  # Observed data
  geom_line(data = df_obs1, aes(x = Generacion, y = Z_obs), color = "#08519c", linewidth = 1) +
  geom_point(data = df_obs1, aes(x = Generacion, y = Z_obs, shape = "Observed values"), 
             color = "#08519c", size = 2.5) +
  geom_line(data = df_obs2, aes(x = Generacion, y = Z_obs), color = "#08519c", linewidth = 1) +
  geom_point(data = df_obs2, aes(x = Generacion, y = Z_obs, shape = "Observed values"), 
             color = "#08519c", size = 2.5) +
  # 95% credible intervals and imputation
  geom_errorbar(data = df_imp, aes(x = Generacion, ymin = IC_Inf, ymax = IC_Sup), 
                color = "#cb181d", width = 0.2, linewidth = 0.9) +
  geom_point(data = df_imp, aes(x = Generacion, y = Media, shape = "Imputation (Mean + 95% CI)"), 
             color = "#cb181d", size = 3.2) +
  scale_x_continuous(breaks = 1:n_gen) +
  scale_shape_manual(
    name = NULL,
    values = c("Observed values" = 19, "Imputation (Mean + 95% CI)" = 18)
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Unknown true values" = "dashed")
  ) +
  labs(
    title = "",
    x = "Generation (i)",
    y = expression(Z[i])
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  )

# Panel B: Posterior densities
p2 <- ggplot(df_mcmc, aes(x = Muestra, fill = Variable, color = Variable)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.35, bins = 25) +
  geom_density(alpha = 0.1, linewidth = 0.8) +
  geom_vline(data = df_vlines, aes(xintercept = Real, color = Variable), 
             linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(
    name = "Latent variable:",
    labels = c(expression(Z[4]), expression(Z[5])),
    values = c("Z[4]" = "#cb181d", "Z[5]" = "#08519c")
  ) +
  scale_color_manual(
    name = "Latent variable:",
    labels = c(expression(Z[4]), expression(Z[5])),
    values = c("Z[4]" = "#a50f15", "Z[5]" = "#08306b")
  ) +
  labs(
    title = "",
    x = expression(Z[i]),
    y = "MCMC density"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.box = "horizontal"
  )

# Panel combination
figura_missing_data <- p1 | p2

# Plot
print(figura_missing_data)

