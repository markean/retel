## 1. Load packages
options(warn = -1)
options(scipen = 999)
suppressMessages(library(doParallel)) # parallel backend
suppressMessages(library(doRNG)) # reproducible parallel loop
library(nloptr)
library(retel)


## 2. Parameters
# Sample size
n <- 5L
# Prior scale
s <- 1
# Tau
tau <- log(n)


## 3. Constants
# Simulation replications
S <- 1e+04L
# Prior location
l <- 0
# Grid wing length
w <- 5
# Optimization
opts <- list("algorithm" = "NLOPT_LD_LBFGS")


## 4. Functions
f <- function(x, par) {
  x - par
}
# AETEL
obj_aetel <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  mean(exp(g %*% l))
}
gr_obj_aetel <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  colMeans(as.vector(exp(g %*% l)) * g)
}
aetel <- function(theta) {
  g <- c(x - theta, -log(n) / 2 * mean(x - theta))
  out <- nloptr(
    x0 = 0, eval_f = obj_aetel, eval_grad_f = gr_obj_aetel, opts = opts, g = g
  )
  lambda <- out$solution
  -(n + 1) * log(obj_aetel(lambda, g)) + lambda * sum(g)
}
# Posterior density functions
etel_post <- function(theta) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) +
    etel(f, x, theta)
  exp(out)
}
aetel_post <- function(theta) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) + aetel(theta)
  exp(out)
}
retel_f_post <- function(theta, tau, mu, sigma) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) +
    retel(f, x, theta, mu, sigma, tau, type = "full")
  exp(out)
}
retel_r_post <- function(theta, tau, mu, sigma) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) +
    retel(f, x, theta, mu, sigma, tau, type = "reduced")
  exp(out)
}


## 5. Simulations
cat("< Simulation replications =", S, ">\n")
cat("< Length of grid =", 2 * w, ">\n")
cat("< Number of points in grid =", 1000L, ">\n")
cat("< tau =", tau, ">\n")
cat("< Sample size =", n, ">\n")
cat("< Scale parameter =", s, ">\n")
set.seed(847675)
cl <- makeCluster(24L)
registerDoParallel(cl)
result <- foreach(
  i = icount(S), .combine = "rbind", .inorder = FALSE,
  .packages = c("nloptr", "retel")
) %dorng% {
  # Sample theta from prior
  theta <- rlogis(1L, l, s)
  # Sample data
  x <- rnorm(n, mean = theta, sd = 1)
  # Grids
  grid_coarse <- seq(mean(x) - w / log(n), mean(x) + w / log(n),
    length.out = 1000L
  )
  grid <- seq(mean(x) - w, mean(x) + w, length.out = 10000L)
  # Index in the grid up to theta
  grid_index <- which.min(grid < theta) - 1L
  # Interpolate function from the coarser grid
  aetel_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) aetel_post(k), FUN.VALUE = numeric(1L))
  )
  retel_f_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) {
      retel_f_post(k, tau = tau, mu = mean(x) - k, sigma = 1)
    },
    FUN.VALUE = numeric(1L)
    )
  )
  retel_r_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) {
      retel_r_post(k, tau = tau, mu = mean(x) - k, sigma = 1)
    },
    FUN.VALUE = numeric(1L)
    )
  )
  # Normalize so that prob > 0
  aetel_post_density <- pmax(0, aetel_post_density_approx(grid))
  retel_f_post_density <- pmax(0, retel_f_post_density_approx(grid))
  retel_r_post_density <- pmax(0, retel_r_post_density_approx(grid))
  # Posterior probability up to theta
  H_aetel <-
    sum(aetel_post_density[seq_len(grid_index)]) / sum(aetel_post_density)
  H_retel_f <-
    sum(retel_f_post_density[seq_len(grid_index)]) / sum(retel_f_post_density)
  H_retel_r <-
    sum(retel_r_post_density[seq_len(grid_index)]) / sum(retel_r_post_density)
  # Special treatment for ETEL due to the convex hull constraint
  if (theta <= min(x)) {
    H_etel <- 0
  } else if (theta >= max(x)) {
    H_etel <- 1
  } else {
    grid_etel <- seq(min(x), max(x), length.out = 1000L)
    grid_etel_index <- which.min(grid_etel < theta) - 1L
    etel_post_density <-
      vapply(grid_etel, function(k) etel_post(k), FUN.VALUE = numeric(1L))
    H_etel <-
      sum(etel_post_density[seq_len(grid_etel_index)]) / sum(etel_post_density)
  }
  c(H_etel, H_aetel, H_retel_f, H_retel_r)
}
stopCluster(cl)


## 6. Results
colnames(result) <- c("etel", "aetel", "retel_f", "retel_r")
etel_ks <- ks.test(result[, "etel"], "punif")
aetel_ks <- ks.test(result[, "aetel"], "punif")
retel_f_ks <- ks.test(result[, "retel_f"], "punif")
retel_r_ks <- ks.test(result[, "retel_r"], "punif")
cat(
  "KS test p-values: \n",
  "ETEL:   ", etel_ks$p.value, "\n",
  "AETEL:  ", aetel_ks$p.value, "\n",
  "RETEL_f: ", retel_f_ks$p.value, "\n",
  "RETEL_r: ", retel_r_ks$p.value, "\n\n"
)
