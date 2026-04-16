################################################################################
# Intensity and Compensator for the Hawkes Process
################################################################################
# Intensity and compensator function calls for likelihoods/posteriors
intensity_at_time <- function(t, lambda0, rho, theta, w, times) {
  diffs <- t - times
  N_k <- colSums(outer(diffs, theta, FUN = "<"))
  
  return(lambda0 + rho*sum(w*N_k))
}

compensator_at_time <- function(t, lambda0, rho, theta, w, times) {
  diffs <- t - times
  min_matrix <- outer(diffs, theta, FUN = pmin)
  inner_sums <- colSums(min_matrix)
  return(lambda0*t + rho*sum(w*inner_sums))
}

# Helpers
logit <- function(x) {log(x/(1 - x))}
logit_inv <- function(x) {exp(x)/(1 + exp(x))}


################################################################################
# full conditional log-posteriors
################################################################################
##### lambda_0
logpost_lambda0 <- function(log_lambda0, rho, theta, w, nu, T_max, times) {
  lambda0 <- exp(log_lambda0)
  
  loglik_sum <- 0
  for (t_i in times) {
    loglik_sum = loglik_sum + log(intensity_at_time(t_i, lambda0, rho, theta, w, times[times < t_i]))
  }
  
  return(
    loglik_sum - lambda0*(T_max + nu)
    + log_lambda0  #Jacobian
  )
}

##### rho
logpost_rho <- function(logit_rho, lambda0, theta, w, r_1, r_2, T_max, times) {
  rho <- logit_inv(logit_rho)
  
  loglik_sum <- 0
  for (t_i in times) {
    loglik_sum = loglik_sum + log(intensity_at_time(t_i, lambda0, rho, theta, w, times[times < t_i]))
  }
  
  return(
    loglik_sum 
    - compensator_at_time(T_max, lambda0, rho, theta, w, times)
    + (r_1 - 1)*log(rho) + (r_2 - 1)*log(1 - rho)
    + log(rho) + log(1 - rho)  #Jacobian (logit)
  )
  
}

##### theta_k
logpost_theta_k <- function(log_theta_k, k, lambda0, rho, theta, w, phi, T_max, times) {
  theta_k <- exp(log_theta_k)
  theta[k] <- theta_k
  
  loglik_sum <- 0
  for (t_i in times) {
    loglik_sum = loglik_sum + log(intensity_at_time(t_i, lambda0, rho, theta, w, times[times < t_i]))
  }
  
  return(
    loglik_sum 
    - compensator_at_time(T_max, lambda0, rho, theta, w, times)
    - phi*theta_k
    + log_theta_k  #Jacobian
  )
}

##### v_k
logpost_v_k <- function(logit_v_k, k, lambda0, rho, theta, v, alpha, T_max, times) {
  v_k <- logit_inv(logit_v_k)
  v[k] <- v_k
  
  # reconstruct raw stick-break weights
  K <- length(theta)
  remaining <- cumprod(c(1, 1 - v))
  w_raw <- c(v, 1)*remaining
  w_raw[K] <- 1 - sum(w_raw[-K])
  # scaled to make mu(t)/rho a probability density
  w <- w_raw/sum(w_raw * theta)
  
  loglik_sum <- 0
  for (t_i in times) {
    loglik_sum = loglik_sum + log(intensity_at_time(t_i, lambda0, rho, theta, w, times[times < t_i]))
  }
  
  return(
    loglik_sum 
    - compensator_at_time(T_max, lambda0, rho, theta, w, times)
    + (alpha - 1)*log(1 - v_k)
    + log(v_k) + log(1 - v_k)  #Jacobian (logit)
  )
}

##### alpha
sample_alpha <- function(alpha, K, w_raw, a_1, a_2) {
  return(rgamma(1, shape = a_1 + K - 1, rate = a_2 - log(w_raw[K])))
}

##### phi
sample_phi <- function(theta, f_1, f_2) {
  return(rgamma(1, shape = f_1 + length(theta), rate = f_2 + sum(theta)))
}

################################################################################
# Metropolis Hastings step function
################################################################################
mh_step <- function(current, logpost_fun, proposal_sd, ...) {
  # Automatically get the name of the function passed in
  fun_name <- deparse(substitute(logpost_fun))
  
  # random walk sampler (in transformed space e.g. log, logit)
  proposal <- rnorm(1, current, proposal_sd)
  
  lp_prop <- logpost_fun(proposal, ...)
  lp_curr <- logpost_fun(current, ...)
  
  # Check for length zero or NA
  if (length(lp_prop) == 0 || is.na(lp_prop)) {
    stop(paste("Error in", fun_name, ": Proposal returned length 0 or NA. Value:", lp_prop))
  }
  if (length(lp_curr) == 0 || is.na(lp_curr)) {
    stop(paste("Error in", fun_name, ": Current returned length 0 or NA. Value:", lp_curr))
  }
  #log(acceptance ratio). the log-posterior functions will have necessary Jacobians
  log_acc <- lp_prop - lp_curr
  
  if (log(runif(1)) < log_acc) {
    return(list(value = proposal, accept = 1))
  } else {
    return(list(value = current, accept = 0))
  }
}

## Code for Metropolis-within-Gibbs sampler
run_sampler <- function(times, T_max, n_iter, init, 
                        prior_params, proposal_sds, progress = TRUE) {
  
  samples <- matrix(NA, n_iter, sum(lengths(init)))
  acceptance <- matrix(NA, n_iter, 6)
  colnames(samples) <- c("lambda0", "rho", 
                         paste0("theta", 1:(length(init$theta))), 
                         paste0("v", 1:(length(init$v))),
                         "alpha", "phi")
  colnames(acceptance) <- c("lambda0", "rho", "theta", "v", "alpha", "phi")
  
  # initialise
  lambda0 <- init$lambda0
  rho     <- init$rho
  theta   <- init$theta
  v       <- init$v
  alpha   <- init$alpha
  phi     <- init$phi
  
  # some preliminary calculations
  N <- length(times)
  N_T <- times[N]
  K <- length(init$theta)
  
  # Stick-break weights
  remaining <- cumprod(c(1, 1 - v))
  w_raw <- c(v, 1)*remaining
  w_raw[K] <- 1 - sum(w_raw[-K])
  # Scaled to make mu(t)/rho a probability density
  w <- w_raw/sum(w_raw*theta)
  
  if (progress) {
    pb <- txtProgressBar(min = 0,      
                         max = n_iter, 
                         style = 3,    
                         width = 50,   
                         char = "=")   
  }
  
  for (iter in 1:n_iter) {
    
    # lambda0 (Metropolis-within-Gibbs)
    res <- mh_step(log(lambda0), logpost_lambda0, proposal_sds$lambda0,
                   rho = rho, theta = theta, w = w, nu = prior_params$nu, 
                   T_max = T_max, times = times)
    lambda0 <- exp(res$value)
    acceptance[iter, "lambda0"] <- res$accept
    
    # rho (Metropolis-within-Gibbs)
    res <- mh_step(logit(rho), logpost_rho, proposal_sds$rho,
                   lambda0 = lambda0, theta = theta, w = w, 
                   r_1 = prior_params$r_1, r_2 = prior_params$r_2, 
                   T_max = T_max, times = times)
    rho <- logit_inv(res$value)
    acceptance[iter, "rho"] <- res$accept
    
    # Theta_k (Metropolis-within-Gibbs)
    k_theta <- sample(1:K, 1, prob = w_raw)
    res <- mh_step(log(theta[k_theta]), logpost_theta_k, proposal_sds$theta_k, 
                   k = k_theta, lambda0 = lambda0, rho = rho, theta = theta, 
                   w = w, phi = phi, T_max = T_max, times = times)
    theta[k_theta] <- exp(res$value)
    acceptance[iter, "theta"] <- res$accept
    
    # v_k (Metropolis-within-Gibbs)
    k_v <- sample(1:(K-1), 1, prob = w_raw[-K]/(sum(w_raw[-K])))
    res <- mh_step(logit(v[k_v]), logpost_v_k, proposal_sds$v_k, k = k_v, 
                   lambda0 = lambda0, rho = rho, theta = theta, v = v, 
                   alpha = alpha, T_max = T_max, times = times)
    v[k_v] <- logit_inv(res$value)
    acceptance[iter, "v"] <- res$accept
    
    # recalculate stick-break weights after sampling v_k
    remaining <- cumprod(c(1, 1 - v))
    w_raw <- c(v, 1)*remaining
    w_raw[K] <- 1 - sum(w_raw[-K])
    w <- w_raw/sum(w_raw*theta)
    
    # alpha (Gibbs step)
    alpha <- sample_alpha(alpha = alpha, K = K, w_raw = w_raw, 
                          a_1 = prior_params$a_1, a_2 = prior_params$f_2)
    acceptance[iter, "alpha"] <- 1
    
    # phi (Gibbs step)
    phi <- sample_phi(theta = theta, 
                      f_1 = prior_params$f_1, f_2 = prior_params$f_2)
    acceptance[iter, "phi"] <- 1
    
    samples[iter, ] <- c(lambda0, rho, theta, v, alpha, phi)
    if (progress) {
      setTxtProgressBar(pb, iter)
    }
  }
  
  if (progress) {
    close(pb)
  }
  
  acceptance_rates <- colMeans(acceptance)
  
  return(list(samples = samples, acceptance_rates = acceptance_rates))
}

