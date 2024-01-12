## 1. Load packages ----
options(warn = -1)
options(scipen = 999)
suppressMessages(library(doParallel)) # parallel backend
suppressMessages(library(doRNG)) # reproducible parallel loop
library(nloptr)


## 2. Functions ----
# ETEL
obj <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  mean(exp(g %*% l))
}
gr_obj <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  colMeans(as.vector(exp(g %*% l)) * g)
}
ETEL <- function(g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  n <- nrow(g)
  out <- nloptr(
    x0 = rep(0, ncol(g)), eval_f = obj, eval_grad_f = gr_obj, opts = opts, g = g
  )
  lambda <- out$solution
  -n * log(obj(lambda, g)) + as.numeric(lambda %*% colSums(g))
}

# AETEL
obj_AETEL <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  mean(exp(g %*% l))
}
gr_obj_AETEL <- function(l, g) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  colMeans(as.vector(exp(g %*% l)) * g)
}
AETEL <- function(theta) {
  g <- c(x - theta, -log(n) / 2 * mean(x - theta))
  out <- nloptr(
    x0 = 0, eval_f = obj_AETEL, eval_grad_f = gr_obj_AETEL, opts = opts, g = g
  )
  lambda <- out$solution
  -(n + 1) * log(obj_AETEL(lambda, g)) + lambda * sum(g)
}

# RETEL
d <- function(l, g, tau = 1) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  n <- nrow(g)
  sum(exp(g %*% l)) / (n + tau)
}
penalty <- function(l, g, tau, mu, sigma) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  n <- nrow(g)
  as.numeric((tau / (n + tau)) * exp(l %*% mu + 0.5 * (t(l) %*% sigma %*% l)))
}
obj2 <- function(l, g, tau, mu, sigma) {
  d(l, g, tau) + penalty(l, g, tau, mu, sigma)
}
gr_obj2 <- function(l, g, tau, mu, sigma) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  n <- nrow(g)
  colSums(as.vector(exp(g %*% l)) * g) / (n + tau) +
    (tau / (n + tau)) *
      as.numeric(exp(l %*% mu + 0.5 * (t(l) %*% sigma %*% l))) *
      as.vector((mu + sigma %*% l))
}
RETEL1 <- function(g, tau, mu, sigma) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  n <- nrow(g)
  out <- nloptr(
    x0 = rep(0, ncol(g)), eval_f = obj2, eval_grad_f = gr_obj2, opts = opts,
    g = g, tau = tau, mu = mu, sigma = sigma
  )
  l <- out$solution
  -(n + 1) * log(obj2(l, g, tau, mu, sigma)) +
    as.numeric(l %*% colSums(g)) + as.numeric(l %*% mu + (t(l) %*% sigma %*% l))
}
RETEL2 <- function(g, tau, mu, sigma) {
  if (isFALSE(is.matrix(g))) {
    g <- as.matrix(g)
  }
  n <- nrow(g)
  out <- nloptr(
    x0 = rep(0, ncol(g)), eval_f = obj2, eval_grad_f = gr_obj2, opts = opts,
    g = g, tau = tau, mu = mu, sigma = sigma
  )
  l <- out$solution
  -n * log(obj2(l, g, tau, mu, sigma)) + as.numeric(l %*% colSums(g))
}
# Posterior density functions
ETEL_post <- function(theta) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) + ETEL(x - theta)
  exp(out)
}
AETEL_post <- function(theta) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) + AETEL(theta)
  exp(out)
}
RETEL1_post <- function(theta, tau, mu, sigma) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) +
    RETEL1(x - theta, tau, mu, sigma)
  exp(out)
}
RETEL2_post <- function(theta, tau, mu, sigma) {
  out <- dlogis(theta, location = l, scale = s, log = TRUE) +
    RETEL2(x - theta, tau, mu, sigma)
  exp(out)
}


## 3. Parameters ----
# Simulation replications
S <- 1e4
# Grid wing length
w <- 5
# Grid points
ng <- 1e3
# Sample size
n <- 100
# Tau
tau <- log(n)
# Prior location
l <- 0
# Prior scale
s <- 5
# Optimization
opts <- list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-04)


## 4. Simulation ----
cat("\n< Simulation replications =", S, ">\n")
cat("< Length of grid =", 2 * w, ">\n")
cat("< Number of points in grid =", ng, ">\n")
cat("< tau =", tau, ">\n")
cat("< Sample size =", n, ">\n")
cat("< Scale parameter =", s, ">\n")
set.seed(847675)
cl <- makeCluster(24L)
registerDoParallel(cl)
result <- foreach(
  i = icount(S), .combine = "rbind", .inorder = F, .packages = c("nloptr")
) %dorng% {
  # Sample theta from prior
  theta <- rlogis(1L, l, s)
  # Sample data
  x <- rnorm(n, mean = theta, sd = 1)
  # Grids
  grid_coarse <- seq(mean(x) - w / log(n), mean(x) + w / log(n),
    length.out = ng
  )
  grid <- seq(mean(x) - w, mean(x) + w, length.out = 10 * ng)
  # Index in the grid up to theta
  grid_index <- which.min(grid < theta) - 1L
  # Interpolate function from the coarser grid
  AETEL_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) AETEL_post(k), FUN.VALUE = numeric(1L))
  )
  RETEL1_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) {
      RETEL1_post(k,
        tau = tau, mu = mean(x) - k, sigma = 1
      )
    },
    FUN.VALUE = numeric(1L)
    )
  )
  RETEL2_post_density_approx <- splinefun(
    grid_coarse,
    vapply(grid_coarse, function(k) {
      RETEL2_post(k,
        tau = tau, mu = mean(x) - k, sigma = 1
      )
    },
    FUN.VALUE = numeric(1L)
    )
  )
  # Normalize so that prob >0
  AETEL_post_density <- pmax(0, AETEL_post_density_approx(grid))
  RETEL1_post_density <- pmax(0, RETEL1_post_density_approx(grid))
  RETEL2_post_density <- pmax(0, RETEL2_post_density_approx(grid))
  # Posterior probability up to theta
  H_AETEL <-
    sum(AETEL_post_density[seq_len(grid_index)]) / sum(AETEL_post_density)
  H_RETEL1 <-
    sum(RETEL1_post_density[seq_len(grid_index)]) / sum(RETEL1_post_density)
  H_RETEL2 <-
    sum(RETEL2_post_density[seq_len(grid_index)]) / sum(RETEL2_post_density)
  # Special treatment for ETEL due to the convex hull constraint
  if (theta <= min(x)) {
    H_ETEL <- 0
  } else if (theta >= max(x)) {
    H_ETEL <- 1
  } else {
    grid_etel <- seq(min(x), max(x), length.out = ng)
    grid_etel_index <- which.min(grid_etel < theta) - 1L
    ETEL_post_density <-
      vapply(grid_etel, function(k) ETEL_post(k), FUN.VALUE = numeric(1L))
    H_ETEL <-
      sum(ETEL_post_density[seq_len(grid_etel_index)]) / sum(ETEL_post_density)
  }
  c(H_ETEL, H_AETEL, H_RETEL1, H_RETEL2)
}
stopCluster(cl)


## 5. Result ----
colnames(result) <- c("etel", "aetel", "retel1", "retel2")
etel_ks <- ks.test(result[, "etel"], "punif")
aetel_ks <- ks.test(result[, "aetel"], "punif")
retel1_ks <- ks.test(result[, "retel1"], "punif")
retel2_ks <- ks.test(result[, "retel2"], "punif")
cat(
  "KS test p-values: \n",
  "ETEL:   ", etel_ks$p.value, "\n",
  "AETEL:  ", aetel_ks$p.value, "\n",
  "RETEL1: ", retel1_ks$p.value, "\n",
  "RETEL2: ", retel2_ks$p.value, "\n\n"
)
path <- paste0("/home/kim.7302/simulation/mb_logn/", "n", n, "s", s, ".RData")
save(result, file = path)
