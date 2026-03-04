#' Simulation of a Poisson process
#' 
#' Simulation of a (possibly inhomogeneous) Poisson process using the thinning 
#' algorithm in  Laub, Taimre, and Pollett (2021) Chapter 4 (Algorithm 1). It 
#' works by generating a homogeneous Poisson process with rate `M` on 
#' \eqn{[0,T]}, and then, in the inhomogeneous case, thinning this process using
#' the provided intensity function `fn`
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param fun A non-negative function (or numeric if simulating a homogeneous
#'   Poisson process) representing the intensity function, \eqn{\lambda (t)} of 
#'   the Poisson process.
#' @param M A non-negative numeric value - upper bound on `fun`. Ignored for 
#'   a homogeneous Poisson process.
#'
#' @return A vector of simulated arrival times from the process.
#'
#' @export
#'
#' @examples
#' # simulate an inhomogeneous Poisson process with Exp(2) intensity
#' t <- simulate_pp(T = 5, fn = function(y\t) {2*exp(-2*t)}, M = 2)
#' # Simulate a homogeneous Poisson process with λ = 2
#' x <- simulate_pp(T = 5, fn = 2)
sim_pp <- function(T, fn, M = NULL) {
    
    # Homogeneous Case 
    if (is.numeric(fn)) {
      stopifnot(
        "'T' must be positive" = is.numeric(T) && length(T) == 1 && T > 0,
        "'fn' must be a single non-negative number" = 
          length(fn) == 1 && !is.na(fn) && fn >= 0
      )
      N_T <- rpois(1, fn * T)
      return(sort(runif(N_T, 0, T)))
    }
    
    # Inhomogeneous Case Guardrails
    stopifnot(
      "'T' must be positive" = is.numeric(T) && length(T) == 1 && T > 0,
      "'fn' must be a function" = is.function(fn),
      "Inhomogeneous simulation requires a single numeric 'M' >= 0" = 
        is.numeric(M) && length(M) == 1 && !is.na(M) && M >= 0
    )
    
    # simulate homogeneous poisson process with rate M on [0, T]
    N_T <- rpois(1, M * T)
    if (N_T == 0) return(numeric(0)) # Handle edge case of 0 arrivals
    
    unthinned_arrival_times <- sort(runif(N_T, 0, T))
    
    # keep arrivals with some probability
    # fn() must be able to handle the vector 'unthinned_arrival_times'
    u <- runif(N_T, 0, M)
    keep <- u <= fn(unthinned_arrival_times)
    
    unthinned_arrival_times[keep]
}


#' Simulation of a Hawkes process by thinning
#' 
#' Simulation of a Hawkes process using Ogata's modified thinning algorithm, 
#' following Chapter 4 (Algorithm 2) of Laub, Taimre, and Pollett (2021). This 
#' is similar to the thinning method for an inhomogeneous Poisson process, but
#' in this case we have no a.s. asymptotic bound M for the conditional intensity 
#' \eqn{\lambda^*(t)}. Instead, we restrict the function space to be only 
#' non-increasing functions of $t$, which allows us to simply use the value of
#' `fn(t_i)` where $t_i$ is the time just after an arrival, which is updated 
#' with each arrival
#'
#' @param T A non-negative numeric value - the end of the interval \eqn{[0,T]}.
#' @param lambda A non-negative numeric value - the background arrival 
#'   component \eqn{\lambda} of the Hawkes conditional intensity
#' @param mu A non-negative, non-increasing function - the excitation 
#'   function/kernel \eqn{\mu(\cdot)} of the Hawkes conditional intensity. Must 
#'   be a vectorised function
#'
#' @return A vector of simulated arrival times from the process.
#'
#' @export
#'
#' @examples
#' # simulate a Hawkes process with background rate 1 and Exp(2) kernel
#' t <- simulate_hp(T=5, lambda = 1, mu = function(t) {2*exp(-2*t)})
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