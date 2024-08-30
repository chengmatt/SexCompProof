# Purpose: To check analyitcal equations and MLE for comparing different 
# parameterizations of sex-composition data
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 8/15/24

library(dgof)
# Functions ---------------------------------------------------------------

inv_logit <- function(opt_pars) {
  p <- c(exp(opt_pars), 1) 
  p <- p / sum(p) 
  return(p)
}

nLL_mult_across <- function(data, pars, nsamp) {
  nLL <- vector(length = nsamp, mode = "integer")
  pars <- exp(pars) / (1 + sum(exp(pars)))
  p <- c(pars, 1 - sum(pars))
  for (i in 1:nsamp) nLL[i] <- dmultinom(data[, i], sum(data[, i]), p, log = TRUE)
  return(-sum(nLL))
}

nLL_bin_mult_within <- function(data, pars, nsamp) {
  nLL <- vector(length = nsamp, mode = "integer")
  # Set up parameter vectors
  prop1 <- exp(pars[length(pars)]) / (1 + exp(pars[length(pars)])) # logit constrain
  p1w <- exp(pars[1:(nages - 1)]) / (1 + sum(exp(pars[1:(nages - 1)]))) # logit constrain
  p1w <- c(p1w, 1 - sum(p1w)) # proportions sex 1
  p2w <- exp(pars[-c(1:(nages - 1), length(pars))]) / (1 + sum(exp(pars[-c(1:(nages - 1), length(pars))]))) # logit constrain
  p2w <- c(p2w, 1 - sum(p2w)) # proportions sex 2
  for(i in 1:nsamp) {
    obs_n1 <- sum(data[1:nages,i]) # get observed proportion sex1
    obs_n2 <- sum(data[,i]) - obs_n1 # get oberseved proportion sex2
    nLL[i] <- dbinom(obs_n1, sum(data[,i]), prop1, log = TRUE) # get nLL from binomial
    nLL[i] <- nLL[i] + dmultinom(data[1:nages, i], obs_n1, p1w, log = TRUE)
    nLL[i] <- nLL[i] + dmultinom(data[-c(1:nages), i], obs_n2, p2w, log = TRUE)
  }
  return(-sum(nLL))
}

rdirmult <- function(theta, n, prob) {
  dir_mult_alpha <- theta * n * prob # variability of dirichlet controlled by sample size and theta
  if(n > 0) {
    sim_dir_mult_probs_tmp <- as.vector(extraDistr::rdirichlet(1, dir_mult_alpha))
    sim_dir_mult <- rmultinom(1, n, sim_dir_mult_probs_tmp) # second stage - multinomial
  }
  if(n == 0) sim_dir_mult <- matrix(0, nrow = length(prob), ncol = 1)
  return(sim_dir_mult)
}

nLL_ddirmult_across <- function(data, pars, nsamp) {
  nLL <- vector(length = nsamp, mode = "integer")
  dm_theta <- exp(pars[length(pars)]) # theta dispersion parameter
  p <- exp(pars[-length(pars)]) / (1 + sum(exp(pars[-length(pars)]))) # logit transform
  p <- c(p, 1 - sum(p)) # set up parameter vector for probabilities
  
  for(i in 1:nsamp) {
    alpha <- dm_theta * sum(data[,i]) * p # derive alpha parameter for a dirichlet multinomial
    nLL[i] <- extraDistr::ddirmnom(data[,i], size = sum(data[,i]), alpha = alpha, log = TRUE)
  } # end i
  return(-sum(nLL))
}

nLL_bin_dirmult_within <- function(data, pars, nsamp) {
  nLL <- vector(length = nsamp, mode = "integer")
  # Set up parameter vectors
  prob <- exp(pars[length(pars)]) / (1 + sum(exp(pars[length(pars)]))) # logit constrain binomial probabilities
  dm_theta <- exp(pars[length(pars) - 1]) # dispersion parameter
  # Logit transform second stage dirichlet multinomial (get parameters here)
  pw_list <- list()
  for(s in 1:2) {
    # get indexing for parameters
    if(s == 1) idx <- 1:(nages-1)
    else idx <- idx + (nages - 1)
    pw_tmp <- exp(pars[idx]) / (1 + sum(exp(pars[idx]))) # logit constrain
    pw_list[[s]] <- c(pw_tmp, 1 - sum(pw_tmp)) # proportions sex 
  } # end s loop
  
  for(i in 1:nsamp) {
    obs_n1 <- sum(data[1:nages,i]) # get observed proportion sex1
    obs_n2 <- sum(data[,i]) - obs_n1 # get oberseved proportion sex2
    alpha_1 <- dm_theta * obs_n1 * pw_list[[1]] # alpha parameter for sex 1
    alpha_2 <- dm_theta * obs_n2 * pw_list[[2]] # alpha parameter for sex 2
    nLL[i] <- dbinom(obs_n1, sum(data[,i]), prob, log = TRUE) # get nLL from binomial
    nLL[i] <- nLL[i] + extraDistr::ddirmnom(data[1:nages, i], obs_n1, alpha_1, log = TRUE)
    nLL[i] <- nLL[i] + extraDistr::ddirmnom(data[-c(1:nages), i], obs_n2, alpha_2, log = TRUE)
  }
  return(-sum(nLL))
}

# nLL_mult_mult_within <- function(data, pars, nsamp, nsexes) {
#   nLL <- vector(length = nsamp, mode = "integer")
#   # Set up parameter vectors
#   # logit constrain multinomial parameters (first stage multinomial)
#   prop <- exp(pars[length(pars):(length(pars) - nsexes + 2)]) / (1 + sum(exp(pars[length(pars):(length(pars) - nsexes + 2)])))
#   prop <- c(prop, 1 - sum(prop)) # proportions for first stage multinomial
#   
#   # Logit transform second stage multinomial (get parameters here)
#   pw_list <- list()
#   for(s in 1:nsexes) {
#     # get indexing for parameters
#     if(s == 1) idx <- 1:(nages-1)
#     else idx <- idx + (nages - 1)
#     pw_tmp <- exp(pars[idx]) / (1 + sum(exp(pars[idx]))) # logit constrain
#     pw_list[[s]] <- c(pw_tmp, 1 - sum(pw_tmp)) # proportions sex 
#   } # end s loop
#   
#   for(i in 1:nsamp) {
#     obs_nw_tmp_vec <- vector() # temp vector to store observed n
#     for(s in 1:nsexes) {
#       # get indexing for data
#       if(s == 1) idx <- 1:nages
#       else idx <- idx + nages
#       obs_nw_tmp <- sum(data[idx,i]) # get observed within n
#       obs_nw_tmp_vec <- c(obs_nw_tmp_vec, obs_nw_tmp)
#       nLL[i] <- nLL[i] + dmultinom(data[idx,i], obs_nw_tmp, pw_list[[s]], log = TRUE) # do multinomial wihtin draws
#     } # s loop
#     # do multinomial sex probabilities
#     nLL[i] <- nLL[i] + dmultinom(obs_nw_tmp_vec, size = sum(obs_nw_tmp_vec), prop, log = TRUE) 
#   } # i loop
#   return(-sum(nLL))
# }

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
pmw <- pmw/sum(pmw) # normalize (within)
pfw <- pfw/sum(pfw) # normalize (within)
pa <- c(propfw*pfw, (1-propfw)*pmw) # across
pw <- c(pfw, pmw) # within

# Simulate via two stage sampling
nfw <- rbinom(nreps, n, propfw) # female samples
nmw <- n - nfw # male samples
Xfw <- sapply(nfw, function(x) rmultinom(1, x, pfw)) # female multinomial (within)
Xmw <- sapply(nmw, function(x) rmultinom(1, x, pmw)) # male multinomial (within)
Xw <- rbind(Xfw, Xmw) # combine within samples

# Simulate via one stage sampling
Xa <- rmultinom(nreps, n, pa) # across
Xfa <- Xa[1:nages,] # females (across)
XMa <- Xa[-(1:nages),] # males 

dgof_bm <- dgof::ks.test(Xw, ecdf(Xa))

# Rename variables for use later
Xa_bm <- Xa
Xw_bm <- Xw

# Comparison (Binomial Multinomial) ---------------------------------------

# Across - Expectation
exp_a <- n * pa # analytical
plot(rowMeans(Xa))
lines(exp_a) # analytical
round(rowMeans(Xa) - exp_a, 3)

# Across - Variance
var_a <- n * pa * (1 - pa) # analytical
plot(apply(Xa, 1, var)) # simulated
lines(var_a) 
round(apply(Xa, 1, var) - var_a, 3)

# Within - Expectation (Females)
exp_fw <- n * propfw * pfw
plot(apply(Xfw, 1, mean))
lines(exp_fw)
round(apply(Xfw, 1, mean) - exp_fw, 3)

# Within - Variance (Females)
var_fw <- n * pfw * propfw * (1-propfw) + propfw^2 * (n * pfw * (1 - pfw)) 
plot(apply(Xfw, 1, var))
lines(var_fw)
round(apply(Xfw, 1, var) - var_fw, 3)

# Within - Expectation (Males)
exp_mw <- n * propmw * pmw
plot(apply(Xmw, 1, mean))
lines(exp_mw)
round(apply(Xmw, 1, mean) - exp_mw, 3)

# Within - Variance (Males)
var_mw <- propmw^2 * (n * pmw * (1 - pmw)) + propmw * (1-propmw) * n * pmw
plot(apply(Xmw, 1, var))
lines(var_mw)
mean(apply(Xmw, 1, var) - var_mw, 3)

# Compare variances
plot(apply(rbind(Xfw, Xmw), 1, var), apply(Xa, 1, var)); abline(0,1)

# MLE Comparison (Across with Across)
fit_a_a <- nlminb(start =  rep(0.1, (nages * 2 - 1)) / (nages * 2 - 1), 
                  objective = nLL_mult_across, 
                  data = Xa[,1:nsamps], 
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# Compare estimate of p
plot(pa)
lines(inv_logit(fit_a_a$par))
round(pa - inv_logit(fit_a_a$par), 3)

# MLE Comparison (Within with Across)
fit_w_a <- nlminb(start =  rep(0.1, (nages * 2 - 1)) / (nages * 2 - 1), 
                  objective = nLL_mult_across, 
                  data = Xw[,1:nsamps], 
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# Compare estimate of p to single multinomial across approach
plot(pa)
lines(inv_logit(fit_w_a$par))
round(pa - inv_logit(fit_w_a$par), 3)

data <- Xa
pars <- c(rep(0.1, (nages - 1)) / (nages - 1), 
          rep(0.1, (nages - 1)) / (nages - 1), 
          0.5)

# MLE Comparison (Across with Within)
fit_a_w <- nlminb(start = pars, 
                  objective = nLL_bin_mult_within,
                  data = Xa[,1:nsamps],
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# Compare obj
fit_a_w$objective - fit_a_a$objective

# Compare estimate of p (females)
plot(pfw)
lines(inv_logit(fit_a_w$par[1:(nages-1)]))

# Males
plot(pmw)
lines(inv_logit(fit_a_w$par[-c(1:(nages-1), length(fit_a_w$par))]))

# Proportion of females
propfw
exp(fit_a_w$par[length(fit_a_w$par)]) / (1 + exp(fit_a_w$par[length(fit_a_w$par)]))

# MLE Comparison (Within with Within)
fit_w_w <- nlminb(start = pars, 
                  objective = nLL_bin_mult_within,
                  data = Xw[,1:nsamps],
                  nsamp = nsamp,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# Compare estimate of p (females)
plot(pfw)
lines(inv_logit(fit_w_w$par[1:(nages-1)]))

# Males
plot(pmw)
lines(inv_logit(fit_w_w$par[-c(1:(nages-1), length(fit_w_w$par))]))

# Proportion of females
propfw
exp(fit_w_w$par[length(fit_w_w$par)]) / (1 + exp(fit_w_a$par[length(fit_w_w$par)]))

# Compare obj
fit_w_a$objective - fit_w_w$objective

# Compare expected values and variances
plot(exp_a , c(exp_fw, exp_mw)); abline(0, 1) 
plot(var_a , c(var_fw, var_mw)); abline(0, 1)

# Set up (Binomial Dirichlet Multinomial) ---------------------------------
set.seed(123)
n <- 100 # samples to draw
nreps <- 1e5 # sims to do
nsamps <- 100 # samples to estiamte with
nages <- 15 # number of ages
nsexes <- 2 # number of sexes
theta <- 1 # dispersion parameter

# Get proportions
pfw <- runif(nages) # female proportions
pmw <- runif(nages) # male propotions
propfw <- sum(pfw) / sum(pfw, pmw) # female ratios
propmw <- 1 - propfw # male ratios
pmw <- pmw/sum(pmw) # normalize (within)
pfw <- pfw/sum(pfw) # normalize (within)
pa <- c(propfw*pfw, (1-propfw)*pmw) # across
pw <- c(pfw, pmw) # within

# Simulate via two stage sampling
nfw <- rbinom(nreps, n, propfw) # female samples
nmw <- n - nfw # male samples
Xfw <- sapply(nfw, function(x) rdirmult(theta, x, pfw)) # female multinomial (within)
Xmw <- sapply(nmw, function(x) rdirmult(theta, x, pmw)) # male multinomial (within)
Xw <- rbind(Xfw, Xmw) # combine within samples

# Simulate via one stage sampling
Xa <- sapply(rep(n, nreps), function(x) rdirmult(theta, x, pa)) # across
Xfa <- Xa[1:nages,] # females (across)
XMa <- Xa[-(1:nages),] # males 

dgof_bdm <- dgof::ks.test(Xw, ecdf(Xa))

# Rename variables for use later
Xa_bdm <- Xa
Xw_bdm <- Xw

# Comparison (Binomial Dirichlet Multinomial) ---------------------------------------
# Dirichlet Multinomial Across parameters
alpha_i <- n * theta * pa
alpha_0 <- sum(alpha_i)

# Across - Expectation
exp_a <- n * (pa / sum(pa)) # expectration of a single dir mult draw
plot(rowMeans(Xa))
lines(exp_a) 
sum(exp_a - rowMeans(Xa))

# Across - Variance
var_a <- n * (alpha_i / alpha_0) * (1 - (alpha_i / alpha_0)) * ((n + alpha_0)/(1 + alpha_0))
plot(var_a)
lines(apply(Xa, 1, var))
mean(apply(Xa, 1, var) - var_a)

# Within - Expectation
exp_fw <- n *  propfw * (propfw * pfw) / sum(propfw * pfw)
exp_mw <- n *  propmw * (propmw * pmw) / sum(propmw * pmw)
plot(apply(rbind(Xfw, Xmw), 1, mean))
lines(c(exp_fw, exp_mw))
mean(apply(rbind(Xfw, Xmw), 1, mean) - c(exp_fw, exp_mw))

# Within - Variance
# dirichlet multinomial females
alpha_i_fw <- theta * n * propfw * pfw
alpha_0_fw <- sum(alpha_i_fw)
var_fw <- ((alpha_i_fw / alpha_0_fw) * (1 - (alpha_i_fw / alpha_0_fw)) * (1 / (1 + alpha_0_fw))) * 
  ((n*propfw * (1 - propfw)) + (n*propfw)^2 + (alpha_0_fw *n*propfw))  + 
  (pfw)^2 * n*propfw * (1-propfw)

# dirichlet multinomial males
alpha_i_mw <- theta * n * propmw * pmw
alpha_0_mw <- sum(alpha_i_mw)
var_mw <- ((alpha_i_mw / alpha_0_mw) * (1 - (alpha_i_mw / alpha_0_mw)) * (1 / (1 + alpha_0_mw))) * 
  ((n*propmw * (1 - propmw)) + (n*propmw)^2 + (alpha_0_mw *n*propmw))  + 
  (pmw)^2 * n*propmw * (1-propmw)

# Dirichlet Multinomial Variance (Within)
plot(apply(rbind(Xfw, Xmw), 1, var))
lines(c(var_fw, var_mw))
mean(apply(rbind(Xfw, Xmw), 1, var) - c(var_fw, var_mw))

# Compare variances
plot(apply(rbind(Xfw, Xmw), 1, var), apply(Xa, 1, var)); abline(0,1)

# MLE Comparison (Across with Across)
fit_a_a <- nlminb(start =  rep(0.1, nages * nsexes), 
                  objective = nLL_ddirmult_across, 
                  data = Xa[,1:nsamps], 
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# Compare proportions
plot(pa)
lines(inv_logit(fit_a_a$par[-length(fit_a_a$par)]))
mean(pa - inv_logit(fit_a_a$par[-length(fit_a_a$par)]))

theta
exp(fit_a_a$par[length(fit_a_a$par)]) # estimates this just right (about)

# MLE Comparison (Within with Across)
fit_w_a <- nlminb(start =  rep(0.1, nages * nsexes), 
                  objective = nLL_ddirmult_across, 
                  data = Xw[,1:nsamps], 
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))


# Compare proportions
plot(pa)
lines(inv_logit(fit_w_a$par[-length(fit_w_a$par)]))
sum(pa - inv_logit(fit_w_a$par[-length(fit_w_a$par)]))

theta
exp(fit_w_a$par[length(fit_w_a$par)]) # (underestimates dispersion when applying within
# to across - i.e., think its less variable then it actually should be)

# MLE Comparison (Across with Within)
fit_a_w <- nlminb(start =  rep(0.1, nages * nsexes), 
                  objective = nLL_bin_dirmult_within, 
                  data = Xa[,1:nsamps], 
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# sex females
plot(pfw)
lines(inv_logit(fit_a_w$par[1:(nages - 1)]))
sum(pfw - inv_logit(fit_a_w$par[1:(nages - 1)]))

# sex males
plot(pmw)
lines(inv_logit(fit_a_w$par[-c(1:(nages - 1), length(fit_a_w$par), length(fit_a_w$par) - 1)]))
sum(pmw - inv_logit(fit_a_w$par[-c(1:(nages - 1), length(fit_a_w$par), length(fit_a_w$par) - 1)]))

theta
exp(fit_a_w$par[length(fit_a_w$par) - 1]) # underestimates theta - more overdispersed than it actually should be

# MLE Comparison (Within with Within)
fit_w_w <- nlminb(start =  rep(0.1, nages * nsexes), 
                  objective = nLL_bin_dirmult_within, 
                  data = Xw[,1:nsamps], 
                  nsamp = nsamps,
                  control = list(iter.max = 1e5, eval.max = 1e5))

# sex females
plot(pfw)
lines(inv_logit(fit_w_w$par[1:(nages - 1)]))
sum(pfw - inv_logit(fit_w_w$par[1:(nages - 1)]))

# sex males
plot(pmw)
lines(inv_logit(fit_w_w$par[-c(1:(nages - 1), length(fit_w_w$par), length(fit_w_w$par) - 1)]))
sum(pmw - inv_logit(fit_w_w$par[-c(1:(nages - 1), length(fit_w_w$par), length(fit_w_w$par) - 1)]))

theta
exp(fit_w_w$par[length(fit_w_w$par) - 1]) # just right!

# Slightly Different objs
fit_w_w$objective
fit_w_a$objective

# Slightly Different objs
fit_a_w$objective
fit_a_a$objective

# Plot Variance Comparisons -----------------------------------------------
par(mfrow = c(2,3))

plot(apply(Xw_bm, 1, mean), apply(Xa_bm, 1, mean), xlab = "Expected Value (Within)", ylab = "Expected Value (Across)",
     main = "A) Multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

plot(apply(Xw_bm, 1, var), apply(Xa_bm, 1, var), xlab = "Variance (Within)", ylab = "Variance (Across)",
     main = "B) Multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

ecdf_Xa_bm <- ecdf(Xa_bm)
ecdf_Xw_bm <- ecdf(Xw_bm)
plot(ecdf_Xa_bm(1:25) - ecdf_Xw_bm(1:25), type = 'l', lwd = 4, col = "blue4",
     ylab = "Difference in ECDF", main = "C) Multinomial"); abline(0, 0, lty = 2, lwd = 2)
text(15, -0.0001, paste("p = ", round(dgof_bm$p.value, 3)))

plot(apply(Xw_bdm, 1, mean), apply(Xw_bdm, 1, mean), xlab = "Expected Value (Within)", ylab = "Expected Value (Across)",
     main = "D) Dirichlet-Multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

plot(apply(Xw_bdm, 1, var), apply(Xa_bdm, 1, var), xlab = "Variance (Within)", ylab = "Variance (Across)",
     main = "E) Dirichlet-Multinomial", pch = 19, col = "blue4", lwd = 3)
abline(0,1, lty = 2, lwd = 2)

ecdf_Xa_bdm <- ecdf(Xa_bdm)
ecdf_Xw_bdm <- ecdf(Xw_bdm)
plot(ecdf_Xw_bdm(1:85) - ecdf_Xa_bdm(1:85), type = 'l', lwd = 4, col = "blue4",
     ylab = "Difference in ECDF", main = "F) Dirichlet-Multinomial"); abline(0, 0, lty = 2, lwd = 2)
text(50, -0.001, "p < 0.05")


# Set up (Multinomial Multinomial) ---------------------------------------
# nreps <- 1e5 # sims to do
# n <- 100 # sample size
# nages <- 15 # number of ages
# nsexes <- 3 # number of sexes
# 
# # Get proportions
# p1w <- runif(nages) # sex 1 proportions
# p2w <- runif(nages) # sex 2 propotions
# p3w <- runif(nages) # sex 3 propotions
# prop1w <- sum(p1w) / sum(p1w, p2w, p3w) # prob for sex 1
# prop2w <- sum(p2w) / sum(p1w, p2w, p3w) # prob for sex 2
# prop3w <- sum(p3w) / sum(p1w, p2w, p3w) # prob for sex 3
# prop <- c(prop1w, prop2w, prop3w) # cocanetate probs for multinomial draw
# 
# # Normalize
# p1w <- p1w/sum(p1w) # (within)
# p2w <- p2w/sum(p2w) # (within)
# p3w <- p3w/sum(p3w) # (within)
# 
# # Vectorize all
# pa <- c(prop1w*p1w, prop2w*p2w, prop3w*p3w) # across
# pw <- c(p1w, p2w, p3w) # within
# 
# # Simulate via two stage sampling
# n_w <- rmultinom(nreps, n, prop) # samples for each sex
# X1w <- sapply(n_w[1,], function(x) rmultinom(1, x, p1w)) # sex 1 multinomial (within)
# X2w <- sapply(n_w[2,], function(x) rmultinom(1, x, p2w)) # sex 2 multinomial (within)
# X3w <- sapply(n_w[3,], function(x) rmultinom(1, x, p3w)) # sex 2 multinomial (within)
# Xw <- rbind(X1w, X2w, X3w) # combine within samples
# 
# # Simulate via one stage sampling
# Xa <- rmultinom(nreps, n, pa) # across
# X1a <- Xa[1:nages,] # sex 1 (across)
# X2a <- Xa[(nages + 1):(nages * 2),] # sex 2 across 
# X3a <- Xa[(nages*2 + 1):(nages * 3),] # sex 2 across 
# 
# dgof_mm <- dgof::ks.test(Xw, ecdf(Xa))
# 
# # Rename variables for use later
# Xa_mm <- Xa
# Xw_mm <- Xw

# # Comparison (Multinomial Multinomial) ---------------------------------------
# 
# # Across - Expectation
# exp_a <- n * pa # analytical
# plot(rowMeans(Xa))
# lines(exp_a) # analytical
# round(rowMeans(Xa) - exp_a, 3)
# 
# # Across - Variance
# var_a <- n * pa * (1 - pa) # analytical
# plot(apply(Xa, 1, var)) # simulated
# lines(var_a) 
# round(apply(Xa, 1, var) - var_a, 3)
# 
# # Within - Expectation (Sex1)
# exp_1w <- n * prop1w * p1w
# plot(apply(X1w, 1, mean))
# lines(exp_1w)
# round(apply(X1w, 1, mean) - exp_1w, 3)
# 
# # Within - Variance (Sex 1)
# var_1w <- prop1w^2 * (n * p1w * (1 - p1w)) + prop1w * (1-prop1w) * n * p1w
# plot(apply(X1w, 1, var))
# lines(var_1w)
# round(apply(X1w, 1, var) - var_1w, 3)
# 
# # Within - Expectation (Sex2)
# exp_2w <- n * prop2w * p2w
# plot(apply(X2w, 1, mean))
# lines(exp_2w)
# round(apply(X2w, 1, mean) - exp_2w, 3)
# 
# # Within - Variance (Sex 1)
# var_2w <- prop2w^2 * (n * p2w * (1 - p2w)) + prop2w * (1-prop2w) * n * p2w
# plot(apply(X2w, 1, var))
# lines(var_2w)
# round(apply(X2w, 1, var) - var_2w, 3)
# 
# # Within - Expectation (Sex3)
# exp_3w <- n * prop3w * p3w
# plot(apply(X3w, 1, mean))
# lines(exp_3w)
# round(apply(X3w, 1, mean) - exp_3w, 3)
# 
# # Within - Variance (Sex 1)
# var_3w <- prop3w^2 * (n * p3w * (1 - p3w)) + prop3w * (1-prop3w) * n * p3w
# plot(apply(X3w, 1, var))
# lines(var_3w)
# round(apply(X3w, 1, var) - var_3w, 3)
# 
# # Compare variances
# plot(apply(rbind(X1w, X2w, X3w), 1, var), apply(Xa, 1, var)); abline(0,1)
# 
# # MLE Comparison (Across with Across)
# fit_a_a <- nlminb(start =  rep(0.1, (nages * 3 - 1)) / (nages * 3 - 1), 
#                   objective = nLL_mult_across, 
#                   data = Xa[,1:nsamps], 
#                   nsamp = nsamps,
#                   control = list(iter.max = 1e5, eval.max = 1e5))
# 
# # Compare estimate of p
# plot(pa)
# lines(inv_logit(fit_a_a$par))
# round(pa - inv_logit(fit_a_a$par), 5)
# 
# # MLE Comparison (Within with Across)
# fit_w_a <- nlminb(start =  rep(0.1, (nages * 3 - 1)) / (nages * 3 - 1), 
#                   objective = nLL_mult_across, 
#                   data = Xw[,1:nsamps], 
#                   nsamp = nsamps,
#                   control = list(iter.max = 1e5, eval.max = 1e5))
# 
# # Compare estimate of p to single multinomial across approach
# plot(pa)
# lines(inv_logit(fit_w_a$par))
# round(pa - inv_logit(fit_w_a$par), 5)
# 
# # MLE Comparison (Within with Within)
# fit_w_w <- nlminb(start =  rep(0.1, 44), 
#                   objective = nLL_mult_mult_within, 
#                   data = Xw[,1:nsamps], 
#                   nsamp = nsamps,
#                   nsexes = nsexes,
#                   control = list(iter.max = 1e5, eval.max = 1e5))
# 
# # Compare these MLEs (same)
# fit_w_w$objective - fit_w_a$objective
# 
# # Compare estimate of p (sex1)
# plot(p1w)
# lines(inv_logit(fit_w_w$par[1:(nages-1)]))
# round(p1w - inv_logit(fit_w_w$par[1:(nages-1)]), 5)
# 
# # Compare estimate of p (sex2)
# plot(p2w)
# lines(inv_logit(fit_w_w$par[1:(nages-1) + (nages - 1)]))
# round(p2w - inv_logit(fit_w_w$par[1:(nages-1) + (nages - 1)]), 5)
# 
# # Compare estimate of p (sex3)
# plot(p3w)
# lines(inv_logit(fit_w_w$par[1:(nages-1) + (nages - 1) + (nages - 1)]))
# round(p2w - inv_logit(fit_w_w$par[1:(nages-1) + (nages - 1) + (nages - 1)]), 5)
# 
# # Compare probabilities of multinomial draws
# plot(prop)
# lines(inv_logit(fit_w_w$par[length(fit_w_w$par):(length(fit_w_w$par) - nsexes + 2)]))
# 
# # MLE Comparison (Across with Within)
# fit_a_w <- nlminb(start =  rep(0.1, 44), 
#                   objective = nLL_mult_mult_within, 
#                   data = Xa[,1:nsamps], 
#                   nsamp = nsamps,
#                   nsexes = nsexes,
#                   control = list(iter.max = 1e5, eval.max = 1e5))
# 
# # Compare these MLEs (same)
# fit_a_w$objective - fit_a_a$objective
# 
# # Compare estimate of p (sex1)
# plot(p1w)
# lines(inv_logit(fit_a_w$par[1:(nages-1)]))
# round(p1w - inv_logit(fit_a_w$par[1:(nages-1)]), 5)
# 
# # Compare estimate of p (sex2)
# plot(p2w)
# lines(inv_logit(fit_a_w$par[1:(nages-1) + (nages - 1)]))
# round(p2w - inv_logit(fit_a_w$par[1:(nages-1) + (nages - 1)]), 5)
# 
# # Compare estimate of p (sex3)
# plot(p3w)
# lines(inv_logit(fit_a_w$par[1:(nages-1) + (nages - 1) + (nages - 1)]))
# round(p2w - inv_logit(fit_a_w$par[1:(nages-1) + (nages - 1) + (nages - 1)]), 5)
# 
# # Compare probabilities of multinomial draws
# plot(prop)
# lines(inv_logit(fit_a_w$par[length(fit_a_w$par):(length(fit_a_w$par) - nsexes + 2)]))


# Set up (Multinomial Dirichlet-Multinomial) ---------------------------------------
# nreps <- 1e5 # sims to do
# n <- 100 # sample size
# nages <- 15 # number of ages
# nsexes <- 3 # number of sexes
# theta <- 1 # dispersion parameter
# 
# # Get proportions
# p1w <- runif(nages) # sex 1 proportions
# p2w <- runif(nages) # sex 2 propotions
# p3w <- runif(nages) # sex 3 propotions
# prop1w <- sum(p1w) / sum(p1w, p2w, p3w) # prob for sex 1
# prop2w <- sum(p2w) / sum(p1w, p2w, p3w) # prob for sex 2
# prop3w <- sum(p3w) / sum(p1w, p2w, p3w) # prob for sex 3
# prop <- c(prop1w, prop2w, prop3w) # cocanetate probs for multinomial draw
# 
# # Normalize
# p1w <- p1w/sum(p1w) # (within)
# p2w <- p2w/sum(p2w) # (within)
# p3w <- p3w/sum(p3w) # (within)
# 
# # Vectorize all
# pa <- c(prop1w*p1w, prop2w*p2w, prop3w*p3w) # across
# pw <- c(p1w, p2w, p3w) # within
# 
# # Simulate via two stage sampling
# n_w <- rmultinom(nreps, n, prop) # samples for each sex
# X1w <- sapply(n_w[1,], function(x) rdirmult(theta, x, p1w)) # sex 1 multinomial (within)
# X2w <- sapply(n_w[2,], function(x) rdirmult(theta, x, p2w)) # sex 2 multinomial (within)
# X3w <- sapply(n_w[3,], function(x) rdirmult(theta, x, p3w)) # sex 2 multinomial (within)
# Xw <- rbind(X1w, X2w, X3w) # combine within samples
# 
# # Simulate via one stage sampling
# Xa <-  sapply(rep(n, nreps), function(x) rdirmult(theta, x, pa))
# X1a <- Xa[1:nages,] # sex 1 (across)
# X2a <- Xa[(nages + 1):(nages * 2),] # sex 2 across 
# X3a <- Xa[(nages*2 + 1):(nages * 3),] # sex 2 across 
# 
# dgof_mdm <- dgof::ks.test(Xw, ecdf(Xa))
# 
# # Rename variables for use later
# Xa_mdm <- Xa
# Xw_mdm <- Xw
