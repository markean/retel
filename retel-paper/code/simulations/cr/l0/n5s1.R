## 1. Load packages
options(warn = -1)
options(scipen = 999)
suppressMessages(library(doParallel)) # parallel backend
suppressMessages(library(doRNG)) # reproducible parallel loop
library(nloptr)
library(retel)


## 2. Parameters
# Prior location
l <- 0
# Sample size
n <- 5L
# Prior scale
s <- 1


## 3. Constants
# Tau
tau <- log(n)
# Simulation replications
S <- 1e+04L
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


## 5. Simulation
cat("< Simulation replications =", S, ">\n")
cat("< Length of grid =", 2 * w, ">\n")
cat("< Number of points in grid =", 1000L, ">\n")
cat("< tau =", tau, ">\n")
cat("< Location parameter =", l, ">\n")
cat("< Sample size =", n, ">\n")
cat("< Scale parameter =", s, ">\n")
set.seed(847675)
cl <- makeCluster(24L)
registerDoParallel(cl)
result <- foreach(
  i = icount(S), .combine = "rbind", .inorder = FALSE,
  .packages = c("nloptr", "retel")
) %dorng% {
  # Sample data
  x <- rnorm(n, mean = 0, sd = 1)
  # Grids
  grid_coarse <- seq(min(mean(x), l) - w / log(n), max(mean(x), l) + w / log(n),
    length.out = 1000L
  )
  grid <- seq(min(mean(x), l) - w, max(mean(x), l) + w, length.out = 10000L)
  grid_etel <- seq(min(x), max(x), length.out = 1000L)
  etel_post_density <-
    vapply(grid_etel, function(k) etel_post(k), FUN.VALUE = numeric(1L))
  # Interpolate function from the coarser grid
  aetel_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) aetel_post(k), FUN.VALUE = numeric(1L))
  )
  retel_f_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) {
      retel_f_post(k, tau = tau, mu = mean(x) - k, sigma = 1)
    }, FUN.VALUE = numeric(1L))
  )
  retel_r_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) {
      retel_r_post(k, tau = tau, mu = mean(x) - k, sigma = 1)
    }, FUN.VALUE = numeric(1L))
  )
  # Normalize so that prob > 0
  aetel_post_density <- pmax(0, aetel_post_density_approx(grid))
  retel_f_post_density <- pmax(0, retel_f_post_density_approx(grid))
  retel_r_post_density <- pmax(0, retel_r_post_density_approx(grid))

  # 95% credible interval
  etel_post_cumsum <- cumsum(etel_post_density) / sum(etel_post_density)
  etel_ci_lower <- grid_etel[which(etel_post_cumsum >= 0.025)[1L]]
  etel_ci_upper <- grid_etel[which(etel_post_cumsum >= 0.975)[1L] - 1L]

  aetel_post_cumsum <- cumsum(aetel_post_density) / sum(aetel_post_density)
  aetel_ci_lower <- grid[which(aetel_post_cumsum >= 0.025)[1L]]
  aetel_ci_upper <- grid[which(aetel_post_cumsum >= 0.975)[1L] - 1L]

  retel_f_post_cumsum <-
    cumsum(retel_f_post_density) / sum(retel_f_post_density)
  retel_f_ci_lower <- grid[which(retel_f_post_cumsum >= 0.025)[1L]]
  retel_f_ci_upper <- grid[which(retel_f_post_cumsum >= 0.975)[1L] - 1L]

  retel_r_post_cumsum <-
    cumsum(retel_r_post_density) / sum(retel_r_post_density)
  retel_r_ci_lower <- grid[which(retel_r_post_cumsum >= 0.025)[1L]]
  retel_r_ci_upper <- grid[which(retel_r_post_cumsum >= 0.975)[1L] - 1L]

  # Results
  etel_ci_length <- etel_ci_upper - etel_ci_lower
  aetel_ci_length <- aetel_ci_upper - aetel_ci_lower
  retel_f_ci_length <- retel_f_ci_upper - retel_f_ci_lower
  retel_r_ci_length <- retel_r_ci_upper - retel_r_ci_lower

  etel_ci_cover <- (0 > etel_ci_lower) & (0 < etel_ci_upper)
  aetel_ci_cover <- (0 > aetel_ci_lower) & (0 < aetel_ci_upper)
  retel_f_ci_cover <- (0 > retel_f_ci_lower) & (0 < retel_f_ci_upper)
  retel_r_ci_cover <- (0 > retel_r_ci_lower) & (0 < retel_r_ci_upper)

  c(
    etel_ci_length, etel_ci_cover, aetel_ci_length, aetel_ci_cover,
    retel_f_ci_length, retel_f_ci_cover, retel_r_ci_length, retel_r_ci_cover
  )
}
stopCluster(cl)


## 5. Results
colnames(result) <- c(
  "etel_length", "etel_cover", "aetel_length", "aetel_cover",
  "retel_f_length", "retel_f_cover", "retel_r_length", "retel_r_cover"
)
cat("\nETEL coverage rate:    ", mean(result[, "etel_cover"]), "\n")
cat("ETEL interval length:   ",
  mean(result[, "etel_length"]), " (",
  round(sd(result[, "etel_length"]) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("AETEL coverage rate:   ", mean(result[, "aetel_cover"]), "\n")
cat("AETEL interval length:  ",
  mean(result[, "aetel_length"]), " (",
  round(sd(result[, "aetel_length"]) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("RETEL_f coverage rate:  ", mean(result[, "retel_f_cover"]), "\n")
cat("RETEL_f interval length: ",
  mean(result[, "retel_f_length"]), " (",
  round(sd(result[, "retel_f_length"]) / sqrt(S), 4L), ")", "\n",
  sep = ""
)
cat("RETEL_r coverage rate:  ", mean(result[, "retel_r_cover"]), "\n")
cat("RETEL_r interval length: ",
  mean(result[, "retel_r_length"]), " (",
  round(sd(result[, "retel_r_length"]) / sqrt(S), 4L), ")", "\n\n",
  sep = ""
)
