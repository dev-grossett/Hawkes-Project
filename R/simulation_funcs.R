#' Simulation of a Poisson process
#' 
#' Simulation of a (possibly inhomogeneous) Poisson process using the thinning 
#' algorithm in  Laub, Taimre, and Pollett (2021) Chapter 4 (Algorithm 1). It 
#' works by generating a homogeneous Poisson process with rate `M` on 
#' \eqn{[0,T]}, and then, in the inhomogeneous case, thinning this process using
#' the provided intensity function `intens`
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param intens A non-negative function (or numeric if simulating a homogeneous
#'   Poisson process) representing the intensity function, \eqn{\lambda (t)} of 
#'   the Poisson process.
#' @param M A non-negative numeric value - upper bound on `fun`. Ignored for 
#'   a homogeneous Poisson process.
#'
#' @return A vector of simulated arrival times from the process.
#' @importFrom stats rpois runif
#' @export
#'
#' @examples
#' # simulate an inhomogeneous Poisson process with Exp(2) intensity
#' t <- sim_pp(T = 5, intens = function(t) {2*exp(-2*t)}, M = 2)
#' # Simulate a homogeneous Poisson process with λ = 2
#' x <- sim_pp(T = 5, intens = 2)
sim_pp <- function(T, intens, M = NULL) {
    
  # Homogeneous Case 
  if (is.numeric(intens)) {
    stopifnot(
      "'T' must be positive" = is.numeric(T) && length(T) == 1 && T > 0,
      "'intens' must be a single non-negative number" = 
        length(intens) == 1 && !is.na(intens) && intens >= 0
    )
    N_T <- rpois(1, intens * T)
    return(sort(runif(N_T, 0, T)))
  }
  
  # Inhomogeneous Case Guardrails
  # TODO: errors for missing 'M', could instead approximate across a grid
  stopifnot(
    "'T' must be positive" = is.numeric(T) && length(T) == 1 && T > 0,
    "'intens' must be a function" = is.function(intens),
    "Inhomogeneous simulation requires a single numeric 'M' >= 0" = 
      is.numeric(M) && length(M) == 1 && !is.na(M) && M >= 0
  )
  
  # simulate homogeneous poisson process with rate M on [0, T]
  N_T <- rpois(1, M * T)
  if (N_T == 0) return(numeric(0)) # Handle edge case of 0 arrivals
  
  unthinned_arrival_times <- sort(runif(N_T, 0, T))
  
  # keep arrivals with some probability
  # intens() must be able to handle the vector 'unthinned_arrival_times'
  u <- runif(N_T, 0, M)
  keep <- u <= intens(unthinned_arrival_times)
  
  unthinned_arrival_times[keep]
}


#' Simulation of a Hawkes process by thinning
#' 
#' Simulation of a Hawkes process using Ogata's modified thinning algorithm, 
#' following Chapter 4 (Algorithm 2) of Laub, Taimre, and Pollett (2021). This 
#' is similar to the thinning method for an inhomogeneous Poisson process, but
#' in this case we have no a.s. asymptotic bound \eqn{M} for the conditional 
#' intensity \eqn{\lambda^*(t)}. Instead, we restrict the function space to be
#' only non-increasing functions of \eqn{t}, which allows us to simply use the 
#' value of `fn(t_i)` where \eqn{t_i} is the time just after an arrival, which 
#' is updated after each arrival
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param lambda A non-negative numeric value - the background arrival 
#'   component \eqn{\lambda} of the Hawkes conditional intensity
#' @param mu A non-negative, non-increasing function - the excitation 
#'   function/kernel \eqn{\mu(\cdot)} of the Hawkes conditional intensity. Must 
#'   be a vectorised function
#'
#' @return A vector of simulated arrival times from the process.
#' @importFrom stats rexp runif
#' @export
#'
#' @examples
#' # simulate a Hawkes process with background rate 1 and Exp(2) kernel
#' t <- sim_hp(T=5, lambda = 1, mu = function(t) {2*exp(-2*t)})
sim_hp <- function(T, lambda, mu) {
  
  # some guardrails
  stopifnot(
    "'T' must be a single numeric" = 
      is.numeric(T) && length(T) == 1 && !is.na(T),
    
    "'T' must be positive" 
      = T > 0,
    
    "'lambda' must be a single numeric" 
      = is.numeric(lambda) && length(lambda) == 1 && !is.na(lambda),
    
    "'lambda' must be non-negative" 
      = lambda >= 0,
    
    "'mu' must be a function" 
      = is.function(mu)
  )
  
  # Start with a sensible ~8KB buffer for vector (will double later if required)
  arrival_times <- numeric(1024)
  t <- 0
  count <- 0
  
  # initialise simulation objects
  arrival_times <- c()
  t <- 0
  while (t < T) {
    # Use only the filled portion for intensity
    current_arrivals <- if (count == 0) numeric(0) else arrival_times[1:count]
    
    # find new upper bound
    M <- lambda + sum(mu(t - arrival_times))
    if (is.na(M)) {
      stop("Intensity calculation returned NA. Check your 'mu' function.")
    }
    
    # generate next candidate point
    t <- t + rexp(1, M)
    if (t > T) break
    
    current_intensity <- lambda + sum(mu(t - arrival_times))
    # keep it with some probability
    if (runif(1, 0, M) <= current_intensity) {
      count <- count + 1
      
      # double the capacity of arrival_times if out of space
      if (count > length(arrival_times)) {
        arrival_times <- c(arrival_times, numeric(length(arrival_times)))
      }
      arrival_times[count] <- t
    }
  }
  
  # trim trailing zeros
  if (count == 0) numeric(0) else arrival_times[1:count]
}


#' Simulation of a Hawkes process by clusters
#' 
#' Simulation of a Hawkes process by clusters i.e. using the immigrant-birth 
#' representation of a Hawkes process. It follows Chapter 4 (Algorithm 3) of
#' Laub, Taimre, and Pollett (2021). For simplicity, we factor the usual 
#' Hawkes kernel as \eqn{\mu(t) = \eta g(t)} where 
#' \eqn{\eta = \int_0^{\infty} \mu(u) du} is the branching ratio, and \eqn{g(t)}
#' is a probability density function on \eqn{\mathbb{R}^+}. 
#' 
#' Simulation via an immigrant-birth process proceeds in two stages. The first 
#' stage consists of randomly generating immigrant arrivals according to a 
#' Poisson process with intensity \eqn{\lambda} on the interval \eqn{[0, T]}. The 
#' second stage then proceeds by randomly generating 
#' \eqn{\mathrm{Poisson}(\eta)} first-generation offspring for each immigrant 
#' generated, where conditional on an immigrant arriving at time \eqn{s}, the 
#' waiting times of its offspring are independent draws from the density 
#' \eqn{g(t)}. Thus, first-generation offspring arrival times are given by 
#' \deqn{s + E_j, \quad E_j \sim g(·)}
#' Each first-generation offspring in turn generates further descendants 
#' according to the same procedure, which produces the characteristic branching 
#' structure. The recursion continues until no further offspring are generated 
#' within \eqn{[0, T]}.
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param lambda A non-negative numeric value - the background arrival 
#'   component \eqn{\lambda} of the Hawkes conditional intensity
#' @param eta A non-negative numeric value - the branching ratio of the Hawkes
#'   process
#' @param family A character string "family" which corresponds to a distribution
#'   with a random generation function "rfamily". 
#' @param ... Additional arguments passed to the random generation function.
#'
#' @return A \code{data.frame} containing the simulated Hawkes process, 
#'   sorted by arrival time, with the following columns:
#' \itemize{
#'   \item \code{time}: The arrival time of the event.
#'   \item \code{gen}: The generation number (0 for immigrants, 1 for 
#'     direct offspring, etc.).
#'   \item \code{id}: A unique integer identifier for each event.
#'   \item \code{parent_id}: The \code{id} of the parent event (0 for 
#'     immigrants).
#' }
#' @importFrom stats rpois
#' @export
#' 
sim_hp_clus <- function(T, lambda, eta, family, ...) {
  # TODO: some exception handling where "rfamily" isn't a functions
  # Fetch the random generation function (e.g., "rexp" if family="exp")
  random_fn <- get(paste0("r", family), mode = "function")
  
  # Stage 1: Generate immigrants using sim_pp function
  # TODO: possibility for inghomogeneous background (not standard Hawkes formulation though)
  immigrant_times <- sim_pp(T, lambda)
  n_immigrants <- length(immigrant_times)
  
  if (n_immigrants == 0) return(data.frame())
  
  # Store results: time, generation, id, parent_id
  res <- data.frame(
    time = immigrant_times,
    gen = 0,
    id = seq_len(n_immigrants),
    parent_id = 0
  )
  
  # 'recent_gen_df are the arrivals that can still "give birth"
  recent_gen_df <- res
  next_id <- n_immigrants + 1
  
  # Stage 2: Recursive branching (Generation by Generation)
  while (nrow(recent_gen_df) > 0) {
    # number of offspring for each event in the recent generation
    num_offspring <- rpois(nrow(recent_gen_df), eta)
    
    # if no offspring in this entire generation, we are done
    if (sum(num_offspring) == 0) break
    
    # map children to their specific parents
    parent_indices <- rep(seq_len(nrow(recent_gen_df)), times = num_offspring)
    parents <- recent_gen_df[parent_indices, ]
    
    # generate child arrival times
    child_times <- parents$time + random_fn(nrow(parents), ...)
    
    # filter by T and build the next generation dataframe slice
    in_window <- child_times <= T
    if (!any(in_window)) break
    
    next_gen_df <- data.frame(
      time = child_times[in_window],
      gen = parents$gen[in_window] + 1,
      id = seq(next_id, length.out = sum(in_window)),
      parent_id = parents$id[in_window]
    )
    
    # update for next loop and storage
    res <- rbind(res, next_gen_df)
    recent_gen_df <- next_gen_df
    next_id <- next_id + sum(in_window)
  }
  
  # Return all arrivals sorted by time
  res[order(res$time), ]
}

