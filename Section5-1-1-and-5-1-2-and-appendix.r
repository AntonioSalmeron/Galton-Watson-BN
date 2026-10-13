# ==============================================================================
# "Bayesian Networks with Branching Processes: Modeling, Inference and Learning"
# ==============================================================================

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("patchwork")) install.packages("patchwork")

library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# Auxiliary functions
# ------------------------------------------------------------------------------

# ------------------
# GW tree generation
# ------------------
generate_gw_tree <- function(offspring_distribution = rpois, param = 1, max_generations = 10) {
  tree <- list()
  tree[[1]] <- list(id = 1, parent = 0)
  generations <- list()
  id_current <- 2
  current_generation <- list(1)
  
  for (gen in 1:max_generations) {
    new_generation <- list()
    for (id_parent in current_generation) {
      n_children <- offspring_distribution(1, param)
      if (n_children > 0) {
        for (i in 1:n_children) {
          tree[[id_current]] <- list(id = id_current, parent = id_parent)
          new_generation <- c(new_generation, id_current)
          id_current <- id_current + 1
        }
      }
    }
    if (length(new_generation) == 0) break
    current_generation <- new_generation
    generations[[gen]] <- new_generation
  }
  return(generations)
}

# -----------------
# Sample generation
# -----------------
generate_sample <- function(n, off_dist, param_value) {
  z <- 0
  while (length(z) < 2) {
    tree <- generate_gw_tree(offspring_distribution = off_dist, param = param_value, max_generations = n)
    for (i in 1:length(tree)) {
      if (length(tree) > 0) {
        z[i] <- sum(1:length(tree[[i]]))
      }
    }
  }
  l <- length(z)
  if (l < n) {
    for (i in (l + 1):n) {
      z[[i]] <- 0
    }
  }
  return(z)
}

# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------

# Panel A: Parameter estimation (MAP, MQR vs Real)
plot_estimates_gg <- function(df_est, exact_val = NULL, is_dynamic = FALSE, title = "A) Estimación de Parámetros") {
  p <- ggplot(df_est, aes(x = Generation)) +
    geom_line(aes(y = MAP, color = "MAP"), linewidth = 0.8) +
    geom_point(aes(y = MAP, color = "MAP"), size = 2) +
    geom_line(aes(y = MQR, color = "MQR"), linewidth = 0.8) +
    geom_point(aes(y = MQR, color = "MQR"), size = 2) +
    scale_x_continuous(breaks = min(df_est$Generation):max(df_est$Generation)) +
    scale_color_manual(
      name = NULL,
      values = c("MAP" = "#08519c", "MQR" = "#cb181d", "Exact" = "#238b45")
    ) +
    labs(title = title, x = "Generation", y = "Parameter value") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position = "top")
  
  if (is_dynamic) {
    p <- p + geom_line(aes(y = Exact, color = "Exact"), linewidth = 0.8, linetype = "dashed") +
      geom_point(aes(y = Exact, color = "Exact"), size = 2)
  } else if (!is.null(exact_val)) {
    p <- p + geom_hline(aes(yintercept = exact_val, color = "Exact"), linetype = "dashed", linewidth = 0.8)
  }
  return(p)
}

# Panel B: Evolution of the posterior densities
plot_posterior_gg <- function(df_dens, title = "B) Evolución de la Distribución A Posteriori") {
  ggplot(df_dens, aes(x = x, y = Density, color = Etapa, group = Etapa)) +
    geom_line(aes(linewidth = Etapa, linetype = Etapa)) +
    scale_color_manual(
      name = "",
      values = c("Prior" = "#08519c", "Intermediate" = "gray60", "Posterior" = "#cb181d")
    ) +
    scale_linewidth_manual(
      name = "",
      values = c("Prior" = 1, "Intermediate" = 0.4, "Posterior" = 1)
    ) +
    scale_linetype_manual(
      name = "",
      values = c("Prior" = "dashed", "Intermediate" = "solid", "Posterior" = "solid")
    ) +
    labs(title = title, x = "Parameter value", y = "Density value") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position = "top")
}

# Panel C and D: Probability of extingtion and population size
plot_series_gg <- function(df_series, y_var, y_label, title, color_hex = "#08519c") {
  ggplot(df_series, aes(x = Generation, y = .data[[y_var]])) +
    geom_line(color = color_hex, linewidth = 0.8) +
    geom_point(color = color_hex, size = 2) +
    scale_x_continuous(breaks = min(df_series$Generation):max(df_series$Generation)) +
    labs(title = title, x = "Generation", y = y_label) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
}

# ==============================================================================
# MODEL 1: POISSON WITH INFORMATIVE GAMMA PRIOR
# ==============================================================================

set.seed(2)
k <- 10; n <- 10
offspring_parameter <- 1.1
z_large <- generate_sample(k, rpois, offspring_parameter)
z <- z_large[1:n]

alpha <- 1; beta <- 10; z0 <- 1
new_alpha <- alpha + z[1]
new_beta <- beta + z0

for (i in 2:length(z)) {
  new_alpha[i] <- new_alpha[i-1] + z[i] 
  new_beta[i]  <- new_beta[i-1] + z[i-1]
}

MCR_estimate <- new_alpha / new_beta
MAP_estimate <- (new_alpha - 1) / new_beta

df_est_m1 <- data.frame(Generation = 2:10, MAP = MAP_estimate[2:10], MQR = MCR_estimate[2:10])

# SAMPLING POSTERIOR DENSITY (Gamma)
x_grid <- seq(0, 2.5, length.out = 500)
df_dens_m1 <- data.frame(x = x_grid, Density = dgamma(x_grid, shape = alpha, rate = beta), Etapa = "Prior")

for (i in 2:(length(z)-1)) {
  df_dens_m1 <- rbind(df_dens_m1, data.frame(x = x_grid, Density = dgamma(x_grid, shape = new_alpha[i], rate = new_beta[i]), Etapa = "Intermediate"))
}
df_dens_m1 <- rbind(df_dens_m1, data.frame(x = x_grid, Density = dgamma(x_grid, shape = new_alpha[length(z)], rate = new_beta[length(z)]), Etapa = "Posterior"))

# PROBABILITY OF EXTINCTION
extinction_probability_poisson_model_1 <- function(N, z0, z, alpha, beta) {
  A <- beta + z0 + sum(z[1:(N-2)])
  B <- alpha + sum(z[1:(N-1)])
  C <- beta + z0 + sum(z[1:(N-1)])
  return((A / C)^B)
}

ep <- numeric(n + 1)
ep[1] <- 1
for (j in 2:(n+1)) ep[j] <- extinction_probability_poisson_model_1(j, z0, z, alpha, beta)

df_ext_m1 <- data.frame(Generation = 2:11, ProbExt = ep[2:11])
df_pop_m1 <- data.frame(Generation = 1:10, Poblacion = z[1:10])

# Multi-panel fig.
p1_m1 <- plot_estimates_gg(df_est_m1, exact_val = 1.1, title = "")
p2_m1 <- plot_posterior_gg(df_dens_m1, title = "")
p3_m1 <- plot_series_gg(df_ext_m1, "ProbExt", "Prob of Extinction", "", "#2b5c8f")
p4_m1 <- plot_series_gg(df_pop_m1, "Poblacion", "Number of individuals", "", "#67000d")

figura_m1_gamma <- (p1_m1 | p2_m1) / (p3_m1 | p4_m1)
print(figura_m1_gamma)

# ==============================================================================
# MODEL 1: POISSON WITH IMPROPER PRIOR (ALPHA=1, BETA=0)
# ==============================================================================

set.seed(2)
z <- generate_sample(10, rpois, 1.1)[1:10]
alpha <- 1; beta <- 0; z0 <- 1

new_alpha <- alpha + z[1]
new_beta <- beta + z0
for (i in 2:length(z)) {
  new_alpha[i] <- new_alpha[i-1] + z[i] 
  new_beta[i]  <- new_beta[i-1] + z[i-1]
}

MCR_estimate <- new_alpha / new_beta
MAP_estimate <- (new_alpha - 1) / new_beta

df_est_m1_imp <- data.frame(Generation = 2:10, MAP = MAP_estimate[2:10], MQR = MCR_estimate[2:10])

# Posterior densities
x_grid <- seq(0, 4, length.out = 500)
df_dens_m1_imp <- data.frame()
for (i in 2:(length(z)-1)) {
  df_dens_m1_imp <- rbind(df_dens_m1_imp, data.frame(x = x_grid, Density = dgamma(x_grid, shape = new_alpha[i], rate = new_beta[i]), Etapa = "Intermediate"))
}
df_dens_m1_imp <- rbind(df_dens_m1_imp, data.frame(x = x_grid, Density = dgamma(x_grid, shape = new_alpha[length(z)], rate = new_beta[length(z)]), Etapa = "Posterior"))

p1_m1_imp <- plot_estimates_gg(df_est_m1_imp, exact_val = 1.1, title = "")
p2_m1_imp <- plot_posterior_gg(df_dens_m1_imp, title = "")

figura_m1_impropia <- p1_m1_imp | p2_m1_imp
print(figura_m1_impropia)

# ==============================================================================
# MODEL 1: GEOMETRIC WITH BETA PRIOR
# ==============================================================================

set.seed(2)
z <- generate_sample(10, rgeom, 0.49)[1:10]
alpha <- 10; beta <- 2; z0 <- 1

new_alpha <- alpha + z[1]
new_beta <- beta + z0
for (i in 2:length(z)) {
  new_alpha[i] <- new_alpha[i-1] + z[i] 
  new_beta[i]  <- new_beta[i-1] + z[i-1]
}

MCR_geom <- new_alpha / (new_alpha + new_beta)
MAP_geom <- (new_alpha - 1) / (new_alpha + new_beta - 2)

df_est_geom <- data.frame(Generation = 2:10, MAP = MAP_geom[2:10], MQR = MCR_geom[2:10])

# Densities
x_grid <- seq(0, 1, length.out = 500)
df_dens_geom <- data.frame(x = x_grid, Density = dbeta(x_grid, alpha, beta), Etapa = "Prior")
for (i in 2:(length(z)-1)) {
  df_dens_geom <- rbind(df_dens_geom, data.frame(x = x_grid, Density = dbeta(x_grid, new_alpha[i], new_beta[i]), Etapa = "Intermediate"))
}
df_dens_geom <- rbind(df_dens_geom, data.frame(x = x_grid, Density = dbeta(x_grid, new_alpha[length(z)], new_beta[length(z)]), Etapa = "Posterior"))

# Extinction
extinction_probability_geometric_simulated_model_1 <- function(N, z0, z, alpha, beta) {
  theta <- (alpha - 1) / (alpha + beta - 2)
  size <- ifelse(N == 1, 0, z[N-1])
  return(dnbinom(0, size, theta))
}

ep_geom <- numeric(n + 1)
ep_geom[1] <- 1
for (j in 2:(n+1)) ep_geom[j] <- extinction_probability_geometric_simulated_model_1(j, z0, z, alpha, beta)

df_ext_geom <- data.frame(Generation = 2:11, ProbExt = ep_geom[2:11])
df_pop_geom <- data.frame(Generation = 1:10, Poblacion = z[1:10])

p1_geom <- plot_estimates_gg(df_est_geom, exact_val = 0.49, title = "")
p2_geom <- plot_posterior_gg(df_dens_geom, title = "")
p3_geom <- plot_series_gg(df_ext_geom, "ProbExt", "Prob. Extinción", "", "#006d2c")
p4_geom <- plot_series_gg(df_pop_geom, "Poblacion", "Casos (Z_g)", "", "#e6550d")

figura_geom_beta <- (p1_geom | p2_geom) / (p3_geom | p4_geom)
print(figura_geom_beta)

# ==============================================================================
# MODELO 2: POISSON WITH EVOLVING PARAMETER
# ==============================================================================

set.seed(2)
z <- generate_sample(10, rpois, 1.1)[1:10]
alpha <- 1; beta <- 10; z0 <- 1

new_alpha <- alpha + z[1]
new_beta <- beta + z0

for (i in 2:length(z)) {
  estimated_lambda <- (new_alpha[i-1] - 1) / new_beta[i-1]
  new_alpha[i] <- 2 * estimated_lambda + z[i] 
  new_beta[i]  <- new_beta[i-1] + 2
}

MCR_estimate <- new_alpha / new_beta
MAP_estimate <- (new_alpha - 1) / new_beta

df_est_m2 <- data.frame(Generation = 2:10, MAP = MAP_estimate[2:10], MQR = MCR_estimate[2:10])

x_grid <- seq(0, 2.5, length.out = 500)
df_dens_m2 <- data.frame(x = x_grid, Density = dgamma(x_grid, shape = alpha, rate = beta), Etapa = "Prior")
for (i in 2:(length(z)-1)) {
  df_dens_m2 <- rbind(df_dens_m2, data.frame(x = x_grid, Density = dgamma(x_grid, shape = new_alpha[i], rate = new_beta[i]), Etapa = "Intermediate"))
}
df_dens_m2 <- rbind(df_dens_m2, data.frame(x = x_grid, Density = dgamma(x_grid, shape = new_alpha[length(z)], rate = new_beta[length(z)]), Etapa = "Posterior"))

ep_m2 <- numeric(n)
ep_m2[1] <- 1
for (j in 2:n) ep_m2[j] <- dpois(0, z[j-1] * MAP_estimate[j])

df_ext_m2 <- data.frame(Generation = 1:10, ProbExt = ep_m2)
df_pop_m2 <- data.frame(Generation = 1:10, Poblacion = z[1:10])

p1_m2 <- plot_estimates_gg(df_est_m2, exact_val = 1.1, title = "")
p2_m2 <- plot_posterior_gg(df_dens_m2, title = "")
p3_m2 <- plot_series_gg(df_ext_m2, "ProbExt", "Prob. of Extinction", "", "#7a0177")
p4_m2 <- plot_series_gg(df_pop_m2, "Poblacion", "No. of individuals", "", "#004529")

figura_m2_evolutivo <- (p1_m2 | p2_m2) / (p3_m2 | p4_m2)
print(figura_m2_evolutivo)


# ==============================================================================
# MODELO 2: SAMPLING WITH DYNAMIC LAMBDAS FROM TABLE 1
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. SIMULATIONS AND COMPUTATIONS
# ------------------------------------------------------------------------------

set.seed(20)
lambdas_tabla1 <- c(5, 0.97, 1.15, 1, 1.05, 0.95, 0.93, 1.12, 1.15, 1.1)
n <- length(lambdas_tabla1)
z0 <- 1

# Generation of the branching process (Z_g)
z <- numeric(n)
z[1] <- rpois(1, lambdas_tabla1[1] * z0)
for (i in 2:n) {
  z[i] <- rpois(1, lambdas_tabla1[i] * z[i-1])
}

# Prior Gamma(1, 10)
alpha_0 <- 1
beta_0  <- 10

new_alpha <- numeric(n)
new_beta  <- numeric(n)

new_alpha[1] <- alpha_0 + z[1]
new_beta[1]  <- beta_0 + z0

for (i in 2:n) {
  estimated_lambda <- (new_alpha[i-1] - 1) / new_beta[i-1]
  new_alpha[i] <- 2 * estimated_lambda + z[i]
  new_beta[i]  <- new_beta[i-1] + 2
}

# MAP and MQR
MCR_estimate <- new_alpha / new_beta
MAP_estimate <- (new_alpha - 1) / new_beta

# Extinction prob.
ep <- numeric(n)
ep[1] <- dpois(0, z0 * MAP_estimate[1])
for (j in 2:n) {
  ep[j] <- dpois(0, z[j-1] * MAP_estimate[j])
}

# ------------------------------------------------------------------------------
# 2. DATA STRUCTURE FOR GGPLOT2
# ------------------------------------------------------------------------------

# Panel A: Estimations (Generations 2 to 10)
df_est_dyn <- data.frame(
  Generation = 2:n,
  MAP = MAP_estimate[2:n],
  MQR = MCR_estimate[2:n],
  Exact = lambdas_tabla1[2:n]
)

# Panel B: Evolution of posterior densities
x_grid <- seq(0, 2.5, length.out = 500)

df_dens_dyn <- data.frame(
  x = x_grid, 
  Density = dgamma(x_grid, shape = alpha_0, rate = beta_0), 
  Etapa = "Prior"
)

for (i in 2:(n - 1)) {
  df_dens_dyn <- rbind(
    df_dens_dyn, 
    data.frame(
      x = x_grid, 
      Density = dgamma(x_grid, shape = new_alpha[i], rate = new_beta[i]), 
      Etapa = "Intermediate"
    )
  )
}

df_dens_dyn <- rbind(
  df_dens_dyn, 
  data.frame(
    x = x_grid, 
    Density = dgamma(x_grid, shape = new_alpha[n], rate = new_beta[n]), 
    Etapa = "Posterior"
  )
)

df_dens_dyn$Etapa <- factor(
  df_dens_dyn$Etapa, 
  levels = c("Prior", "Intermediate", "Posterior")
)

# Panels C and D: Time series
df_ext_dyn <- data.frame(Generation = 1:n, ProbExt = ep)
df_pop_dyn <- data.frame(Generation = 1:n, Poblacion = z)

# ------------------------------------------------------------------------------
# 3. Panel generation
# ------------------------------------------------------------------------------

# Panel A: Estimation ofLambda
p1_dyn <- ggplot(df_est_dyn, aes(x = Generation)) +
  geom_line(aes(y = MAP, color = "MAP"), linewidth = 0.8) +
  geom_point(aes(y = MAP, color = "MAP"), size = 2) +
  geom_line(aes(y = MQR, color = "MQR"), linewidth = 0.8) +
  geom_point(aes(y = MQR, color = "MQR"), size = 2) +
  geom_line(aes(y = Exact, color = "Exact"), linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y = Exact, color = "Exact"), size = 2) +
  scale_x_continuous(breaks = 2:n) +
  scale_color_manual(
    name = NULL,
    values = c("MAP" = "#08519c", "MQR" = "#cb181d", "Exacto" = "#238b45")
  ) +
  labs(
    title = "", 
    x = "Generation", 
    y = "Parameter value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"), 
    panel.grid.minor = element_blank(), 
    legend.position = "top"
  )

# Panel B: Posterior densities
p2_dyn <- ggplot(df_dens_dyn, aes(x = x, y = Density, color = Etapa, group = Etapa)) +
  geom_line(aes(linewidth = Etapa, linetype = Etapa)) +
  scale_color_manual(
    name = "",
    values = c("Prior" = "#08519c", "Intermediate" = "gray60", "Posterior" = "#cb181d")
  ) +
  scale_linewidth_manual(
    name = "",
    values = c("Prior" = 1, "Intermediate" = 0.4, "Posterior" = 1)
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Prior" = "dashed", "Intermediate" = "solid", "Posterior" = "solid")
  ) +
  labs(
    title = "", 
    x = expression(lambda[i]), 
    y = "Density value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"), 
    panel.grid.minor = element_blank(), 
    legend.position = "top"
  )

# Panel C: Prob. of Extinction
p3_dyn <- ggplot(df_ext_dyn, aes(x = Generation, y = ProbExt)) +
  geom_line(color = "#7a0177", linewidth = 0.8) +
  geom_point(color = "#7a0177", size = 2) +
  scale_x_continuous(breaks = 1:n) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "", 
    x = "Generation", 
    y = "Prob. of Extinction"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"), 
    panel.grid.minor = element_blank()
  )

# Panel D: Population size
p4_dyn <- ggplot(df_pop_dyn, aes(x = Generation, y = Poblacion)) +
  geom_line(color = "#004529", linewidth = 0.8) +
  geom_point(color = "#004529", size = 2) +
  scale_x_continuous(breaks = 1:n) +
  labs(
    title = "", 
    x = "Generation (g)", 
    y = expression("No. of individuals (" * Z[g] * ")")
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"), 
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 4. Multi panel composition
# ------------------------------------------------------------------------------

figura_m2_tabla1_homogenea <- (p1_dyn | p2_dyn) / (p3_dyn | p4_dyn)
print(figura_m2_tabla1_homogenea)

