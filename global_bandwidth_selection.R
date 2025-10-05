library(ggplot2)
set.seed(26)

m <- function(x) {
  sin( (x/3 + 0.1)^(-1) )
}

gen_X <- function(n, alpha, beta) {
  rbeta(n, alpha, beta)
}

gen_Y <- function(X, sigma2 = 1) {
  m(X) + rnorm(length(X), 0, sqrt(sigma2))
}

block_indices <- function(X, N) {
  quantiles <- quantile(X, probs = seq(0, 1, length.out = N + 1), na.rm = TRUE)
  bin <- cut(X, breaks = quantiles, include.lowest = TRUE, right = FALSE)
  split(seq_along(X), bin)
}

fit_block_poly <- function(x, y) {
  df <- data.frame(y = y, x = x)
  fit <- lm(y ~ poly(x, 4, raw = TRUE), data = df)
  
  coeff <- coef(fit)  
  m_hat <- fitted(fit)
  
  beta2 <- coeff[3]
  beta3 <- coeff[4]
  beta4 <- coeff[5]
  m_hat2 <- 2*beta2 + 6*beta3*x + 12*beta4*x^2
  
  list(coef = coeff, m_hat = m_hat, m_hat2 = m_hat2, rss = sum((y - m_hat)^2))
}

estimate_theta_sigma <- function(X, Y, N) {
  idx_list <- block_indices(X, N)
  
  total_rss <- 0
  m2_sq_sum <- 0
  n <- length(X)
  
  for (i in idx_list) {
    if (length(i) < 5 || length(unique(X[i])) < 5) next
    
    x <- X[i]; y <- Y[i]
    fit <- fit_block_poly(x, y)
    total_rss <- total_rss + fit$rss
    m2_sq_sum <- m2_sq_sum + sum(fit$m_hat2^2)
  }
  
  theta22_hat <- m2_sq_sum / n
  df <- n - 5 * N
  sigma2_hat <- total_rss / df
  
  list(theta22 = theta22_hat, sigma2 = sigma2_hat, RSS = total_rss)
}

h_amise_from_estimates <- function(n, sigma2_hat, theta22_hat, interval_len = 1) {
  n^(-1/5) * ((35 * sigma2_hat * interval_len) / theta22_hat)^(1/5)
}

mallows_cp <- function(RSS_N, RSS_Nmax, n, N, Nmax) {
  denom <- RSS_Nmax / (n - 5 * Nmax)
  RSS_N / denom - (n - 10 * N)
}

choose_N_by_cp <- function(X, Y) {
  n <- length(X)
  Nmax <- max(min(floor(n / 20), 5), 1)
  Ns <- 1:Nmax
  
  est_max <- estimate_theta_sigma(X, Y, Nmax)
  RSS_max <- est_max$RSS
  
  tab <- data.frame(N = integer(), theta22 = double(), sigma2 = double(),
                    RSS = double(), Cp = double(), h = double())
  
  for (N in Ns) {
    est <- estimate_theta_sigma(X, Y, N)
    CpN <- mallows_cp(est$RSS, RSS_max, n, N, Nmax)
    hN  <- h_amise_from_estimates(n, est$sigma2, est$theta22, interval_len = 1)
    tab <- rbind(tab, data.frame(
      N = N, theta22 = est$theta22, sigma2 = est$sigma2, RSS = est$RSS, Cp = CpN, h = hN
    ))
  }
  
  i <- which.min(tab$Cp)
  list(Nopt = tab$N[i], Cp = tab$Cp[i], table = tab)
}

estimate_h_for_dataset <- function(n, alpha, beta, sigma2 = 1) {
  
  X <- gen_X(n, alpha, beta)
  Y <- gen_Y(X, sigma2)
  
  sel <- choose_N_by_cp(X, Y)
  return(list(
    X = X, Y = Y,
    table_N = sel$table,
    Nopt = sel$Nopt,
    h_at_Nopt = if (is.na(sel$Nopt)) NA_real_ else sel$table$h[sel$table$N == sel$Nopt]
  ))
}

run_simulation <- function(
    n_vals = c(200, 500, 1000),
    alpha_beta_list = list(c(5,2), c(2,2), c(1,2), c(0.5,2), c(0.1,2)),
    sigma2 = 1,
    R = 100
) {
  res_all <- list()
  counter <- 1
  
  for (n in n_vals) {
    for (ab in alpha_beta_list) {
      alpha <- ab[1]; beta <- ab[2]
      
      for (r in 1:R) {
        out <- estimate_h_for_dataset(n, alpha, beta, sigma2)
        tab <- out$table_N
        tab$rep <- r
        tab$n   <- n
        tab$alpha <- alpha
        tab$beta  <- beta
        tab$h_opt <- out$h_at_Nopt
        tab$Nopt  <- out$Nopt
        
        res_all[[counter]] <- tab
        counter <- counter + 1
      }
    }
  }
  
  do.call(rbind, res_all)
}

df <- run_simulation()




n_pick <- 500
alpha_pick <- 2
beta_pick <- 2
df_sub <- subset(df, n == n_pick & alpha == alpha_pick & beta == beta_pick)
df_mean <- aggregate(h ~ N, data = df_sub, FUN = mean, na.rm = TRUE)

p1 <- ggplot(df_mean, aes(x = N, y = h)) +
  geom_line(color = "steelblue", size = 0.8) +
  geom_point(color = "darkblue", size = 1.5) +
  labs(
    title = sprintf("Behavior of h_AMISE as a function of N (n = %d, α = %.1f, β = %.1f)", 
                    n_pick, alpha_pick, beta_pick),
    x = "# of blocks (N)",
    y = expression(h[AMISE])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("h_vs_N.png", p1, width = 7, height = 4.5, dpi = 300)



df_Nopt <- aggregate(Nopt ~ n, data = df, FUN = mean, na.rm = TRUE)

p2 <- ggplot(df_Nopt, aes(x = n, y = Nopt)) +
  geom_line(color = "steelblue", size = 0.8) +
  geom_point(color = "darkblue", size = 1.5) +
  labs(
    title = "Evolution of optimal block number Nopt with sample size n",
    x = "Sample size (n)",
    y = "Optimal number of blocks (Nopt)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Nopt_vs_n.png", p2, width = 7, height = 4.5, dpi = 300)



df_hopt <- subset(df, !is.na(h_opt) & is.finite(h_opt))
df_mean_shape <- aggregate(h_opt ~ alpha, data = df_hopt, FUN = mean, na.rm = TRUE)

p3 <- ggplot(df_mean_shape, aes(x = factor(alpha), y = h_opt)) +
  geom_col(fill = "steelblue", color = "black", width = 0.7) +
  labs(
    title = expression(bold("Effect of the shape parameter " * alpha * " on estimated " * h[AMISE])),
    x = expression(alpha),
    y = expression("Mean estimated " ~ h[AMISE])
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12)
  )

ggsave("shape_effect_alpha.png", p3, width = 7, height = 4.5, dpi = 300)






