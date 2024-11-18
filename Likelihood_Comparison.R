# Purpose: To check analyitcal equations and MLE for comparing different 
# parameterizations of sex-composition data
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 8/15/24

# Functions ---------------------------------------------------------------

rdirmult <- function(theta, n, prob) {
  dir_mult_alpha <- theta * n * prob # variability of dirichlet controlled by sample size and theta
  if(n > 0) {
    sim_dir_mult_probs_tmp <- as.vector(extraDistr::rdirichlet(1, dir_mult_alpha))
    sim_dir_mult <- rmultinom(1, n, sim_dir_mult_probs_tmp) # second stage - multinomial
  }
  if(n == 0) sim_dir_mult <- matrix(0, nrow = length(prob), ncol = 1)
  return(sim_dir_mult)
}

# Set up (Binomial Multinomial) ------------------------------------------------------------------
set.seed(123)
nreps <- 1e5 # sims to do
nsamps <- 100
n <- 100 # sample size
nages <- 15 # number of ages

# Get proportions
pfw <- runif(nages) # female proportions
pmw <- runif(nages) # male propotions
propfw <- sum(pfw) / sum(pfw, pmw) # female ratios
propmw <- 1 - propfw # male ratios
pmw <- pmw/sum(pmw) # normalize (split)
pfw <- pfw/sum(pfw) # normalize (split)
pa <- c(propfw*pfw, (1-propfw)*pmw) # joint
pw <- c(pfw, pmw) # split

# Simulate via two stage sampling
nfw <- rbinom(nreps, n, propfw) # female samples
nmw <- n - nfw # male samples
Xfw <- sapply(nfw, function(x) rmultinom(1, x, pfw)) # female multinomial (split)
Xmw <- sapply(nmw, function(x) rmultinom(1, x, pmw)) # male multinomial (split)
Xw <- rbind(Xfw, Xmw) # combine split samples

# Simulate via one stage sampling
Xa <- rmultinom(nreps, n, pa) # joint
Xfa <- Xa[1:nages,] # females (joint)
XMa <- Xa[-(1:nages),] # males 

# Rename variables for use later
Xa_bm <- Xa
Xw_bm <- Xw

# Comparison (Binomial Multinomial) ---------------------------------------

# joint - Expectation
exp_a <- n * pa # analytical
plot(rowMeans(Xa))
lines(exp_a) # analytical
round(rowMeans(Xa) - exp_a, 3)

# joint - Variance
var_a <- n * pa * (1 - pa) # analytical
plot(apply(Xa, 1, var)) # simulated
lines(var_a) 
round(apply(Xa, 1, var) - var_a, 3)

# split - Expectation (Females)
exp_fw <- n * propfw * pfw
plot(apply(Xfw, 1, mean))
lines(exp_fw)
round(apply(Xfw, 1, mean) - exp_fw, 3)

# split - Variance (Females)
var_fw <- n * pfw * propfw * (1-propfw) + propfw^2 * (n * pfw * (1 - pfw)) 
plot(apply(Xfw, 1, var))
lines(var_fw)
round(apply(Xfw, 1, var) - var_fw, 3)

# split - Expectation (Males)
exp_mw <- n * propmw * pmw
plot(apply(Xmw, 1, mean))
lines(exp_mw)
round(apply(Xmw, 1, mean) - exp_mw, 3)

# split - Variance (Males)
var_mw <- propmw^2 * (n * pmw * (1 - pmw)) + propmw * (1-propmw) * n * pmw
plot(apply(Xmw, 1, var))
lines(var_mw)
mean(apply(Xmw, 1, var) - var_mw, 3)

# Compare variances
plot(apply(rbind(Xfw, Xmw), 1, var), apply(Xa, 1, var)); abline(0,1)

# Compare expected values and variances
plot(exp_a , c(exp_fw, exp_mw)); abline(0, 1) 
plot(var_a , c(var_fw, var_mw)); abline(0, 1)

# Set up (Binomial Dirichlet Multinomial) ---------------------------------
set.seed(123)
n <- 100 # samples to draw
nreps <- 1e5 # sims to do
nsamps <- 500 # samples to estiamte with
nages <- 15 # number of ages
nsexes <- 2 # number of sexes
theta <- 1 # dispersion parameter

# Get proportions
pfw <- runif(nages) # female proportions
pmw <- runif(nages) # male propotions
propfw <- sum(pfw) / sum(pfw, pmw) # female ratios
propmw <- 1 - propfw # male ratios
pmw <- pmw/sum(pmw) # normalize (split)
pfw <- pfw/sum(pfw) # normalize (split)
pa <- c(propfw*pfw, (1-propfw)*pmw) # joint
pw <- c(pfw, pmw) # split

# Simulate via two stage sampling
nfw <- rbinom(nreps, n, propfw) # female samples
nmw <- n - nfw # male samples
Xfw <- sapply(nfw, function(x) rdirmult(theta, x, pfw)) # female multinomial (split)
Xmw <- sapply(nmw, function(x) rdirmult(theta, x, pmw)) # male multinomial (split)
Xw <- rbind(Xfw, Xmw) # combine split samples

# Simulate via one stage sampling
Xa <- sapply(rep(n, nreps), function(x) rdirmult(theta, x, pa)) # joint
Xfa <- Xa[1:nages,] # females (joint)
XMa <- Xa[-(1:nages),] # males 

# Rename variables for use later
Xa_bdm <- Xa
Xw_bdm <- Xw

# Comparison (Binomial Dirichlet Multinomial) ---------------------------------------
# Dirichlet Multinomial joint parameters
alpha_i <- n * theta * pa
alpha_0 <- sum(alpha_i)

# joint - Expectation
exp_a <- n * (pa / sum(pa)) # expectration of a single dir mult draw
plot(rowMeans(Xa))
lines(exp_a) 
sum(exp_a - rowMeans(Xa))

# joint - Variance
var_a <- n * (alpha_i / alpha_0) * (1 - (alpha_i / alpha_0)) * ((n + alpha_0)/(1 + alpha_0))
plot(var_a)
lines(apply(Xa, 1, var))
mean(apply(Xa, 1, var) - var_a)

# split - Expectation
exp_fw <- n *  propfw * (propfw * pfw) / sum(propfw * pfw)
exp_mw <- n *  propmw * (propmw * pmw) / sum(propmw * pmw)
plot(apply(rbind(Xfw, Xmw), 1, mean))
lines(c(exp_fw, exp_mw))
mean(apply(rbind(Xfw, Xmw), 1, mean) - c(exp_fw, exp_mw))

# split - Variance
# dirichlet multinomial females
alpha_i_fw <- theta * n * propfw * pfw
alpha_0_fw <- sum(alpha_i_fw)
var_fw <- ((pfw) * (1 - (pfw)) * (1 / (1 + alpha_0_fw))) * 
  ((n*propfw * (1 - propfw)) + (n*propfw)^2 + (alpha_0_fw *n*propfw))  + 
  (pfw)^2 * n*propfw * (1-propfw)

# dirichlet multinomial males
alpha_i_mw <- theta * n * propmw * pmw
alpha_0_mw <- sum(alpha_i_mw)
var_mw <- ((pmw) * (1 - (pmw)) * (1 / (1 + alpha_0_mw))) * 
  ((n*propmw * (1 - propmw)) + (n*propmw)^2 + (alpha_0_mw *n*propmw))  + 
  (pmw)^2 * n*propmw * (1-propmw)

# Dirichlet Multinomial Variance (split)
plot(apply(rbind(Xfw, Xmw), 1, var))
lines(c(var_fw, var_mw))
mean(apply(rbind(Xfw, Xmw), 1, var) - c(var_fw, var_mw))

# Compare variances
plot(apply(rbind(Xfw, Xmw), 1, var), apply(Xa, 1, var)); abline(0,1, col = 'red', lwd = 3)


# Plot Variance Comparisons -----------------------------------------------
par(mfrow = c(2,3))

plot(apply(Xw_bm, 1, mean), apply(Xa_bm, 1, mean), xlab = "Expected Value (Split)", ylab = "Expected Value (Joint)",
     main = "A) Multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

plot(apply(Xw_bm, 1, var), apply(Xa_bm, 1, var), xlab = "Variance (Split)", ylab = "Variance (Joint)",
     main = "B) Multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

ecdf_Xa_bm <- ecdf(Xa_bm)
ecdf_Xw_bm <- ecdf(Xw_bm)
plot(ecdf_Xa_bm(1:15) - ecdf_Xw_bm(1:15), type = 'l', lwd = 4, col = "blue4", ylim = c(-0.003, 0.003),
     ylab = "Difference in ECDF", main = "C) Multinomial"); abline(0, 0, lty = 2, lwd = 2)

plot(apply(Xw_bdm, 1, mean), apply(Xw_bdm, 1, mean), xlab = "Expected Value (Split)", ylab = "Expected Value (Joint)",
     main = "D) Dirichlet-multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

plot(apply(Xw_bdm, 1, var), apply(Xa_bdm, 1, var), xlab = "Variance (Split)", ylab = "Variance (Joint)",
     main = "E) Dirichlet-multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

ecdf_Xa_bdm <- ecdf(Xa_bdm)
ecdf_Xw_bdm <- ecdf(Xw_bdm)
plot(ecdf_Xw_bdm(1:15) - ecdf_Xa_bdm(1:15), type = 'l', lwd = 4, col = "blue4", ylim = c(-0.003, 0.003),
     ylab = "Difference in ECDF", main = "F) Dirichlet-multinomial"); abline(0, 0, lty = 2, lwd = 2)

