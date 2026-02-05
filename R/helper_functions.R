# source files

# *** This function creates item-level matrices for true score, errors and observed scores, 
# to simulate data under specific conditions ***
get_TEX_matrices <- function(n_items_in, pop_reliability, covTiTj, lbE, ubE){
  
  if (F) {
    n_items_in = 10
    pop_reliability = .90
    d = 1
    lbE = 1
    ubE = 3
  }
  
  d <- covTiTj
  
  # initialize
  covTT <- 1 * d
  Tcovmat <- matrix(covTT, n_items_in, n_items_in); diag(Tcovmat) <- 1
  
  Ecovmat <- diag(n_items_in)
  diag(Ecovmat) <- runif(n_items_in, lbE, ubE)
  
  # rescale error variances to reach reliablity
  varT <- sum(Tcovmat)
  varE <- sum(Ecovmat)
  es <- (varT - pop_reliability * varT)/(pop_reliability*varE) 
  Ecovmat <- Ecovmat * es
  
  # standardize
  Xcovmat <- Tcovmat + Ecovmat
  varX <- sum(Xcovmat)
  Xcovmat <- Xcovmat / varX
  Tcovmat <- Tcovmat / varX
  Ecovmat <- Ecovmat / varX
  
  # compuate stats and checks
  varX <- sum(Xcovmat)
  varT <- sum(Tcovmat)
  varE <- sum(Ecovmat)
  
  pop_SEM <- sqrt(varE)
  pop_rho <- sum(Tcovmat)/sum(Xcovmat)
  pop_alpha <-  n_items_in/(n_items_in-1) * (1 - sum(diag(Xcovmat))/varX)
  pop_SEE <- pop_SEM * sqrt(pop_rho)
  
  return(list(sigmaT = Tcovmat, sigmaE = Ecovmat, sigmaX = Xcovmat, 
              pop_rho = pop_rho, pop_alpha = pop_alpha, 
              pop_SEM = pop_SEM, pop_SEE= pop_SEE ))
  
}



# *** This function creates data sets using the double stochastic model and computes
# sample statistics for reliability.

simulate_dataset <- function(n_persons, CTTpop, method = "single", format = "wide"){

  stop("Under construction")  
  
  if (F) {
    n_persons <- 100
    CTTpop <- CTT_population
  }
  
  if (method ==  "double"){ 
    Ts <- mvrnorm(n_persons, mu = rep(0, nrow(CTTpop$sigmaX)), Sigma = CTTpop$sigmaT)
    Es <- mvrnorm(n_persons, mu = rep(0, nrow(CTTpop$sigmaX)), Sigma = CTTpop$sigmaE)    
    X <- Ts + Es
  }
  
  
  if (method ==  "single"){
    X <- mvrnorm(n_persons, mu = rep(0, nrow(CTTpop$sigmaX)), Sigma = CTTpop$sigmaX)    
  }

  #stop("incorrect method specified")    
}




# Simulate item scores according to the graded response model
#'
#' @param nN Numeric; Sample size.
#' @param nJ Numeric; Number of items.
#' @param nM Numeric; Number of categories minus 1 per item.
#' @param abil Numeric; vector of length nN (default: N(0,1)).
#' @param itempars Numeric; matrix of size nJ x (nM + 1):
#'   - column 1: discrimination parameters
#'   - columns 2:(nM+1): category threshold parameters
#' @param itemfixed Logical; if TRUE, the same item parameters are generated across runs
#'
#' @return A list containing:
#' * \code{true.mat}: matrix of size nN x nJ: true model implied item score s 
#' * \code{isc.mat}: matrix of size nN x nJ: observed item scores 
#'
simulateGRM <- function(nN = 1000,
                        nJ = 10,
                        nM = 1,
                        abil = NULL,
                        itempars = NULL,
                        itemfixed = TRUE) {
  
  ## ---------------------------
  ## Input checks
  ## ---------------------------
  stopifnot(nN > 0, nJ > 0, nM >= 1)
  
  if (!is.null(abil) && !is.vector(abil)) {
    stop("`abil` must be a numeric vector.")
  }
  
  if (!is.null(itempars) && !is.matrix(itempars)) {
    stop("`itempars` must be a matrix.")
  }
  
  
  ## ---------------------------
  ## Item parameters
  ## ---------------------------
  
  if (is.null(itempars)) {
    itempars <- matrix(NA_real_, nrow = nJ, ncol = nM + 1)
    
    if (isTRUE(itemfixed)) {
      old_seed <- .Random.seed
      set.seed(50226)
      on.exit({ .Random.seed <<- old_seed }, add = TRUE)
    }
    
    # Discrimination parameters
    itempars[, 1] <- runif(nJ, 1, 3)
    
    # Thresholds: equally spaced intervals over [-2, 2]
    bounds <- seq(-2, 2, length.out = nM + 1)
    for (m in seq_len(nM)) {
      itempars[, m + 1] <- runif(nJ, bounds[m], bounds[m + 1])
    }
  }
  
  
  ## ---------------------------
  ## Ability parameters
  ## ---------------------------
  if (is.null(abil)) {
    set.seed(NULL)  # resets RNG using current time / entropy
    abil <- rnorm(nN)
  }
  
  ## ---------------------------
  ## Output objects
  ## ---------------------------
  true.mat <- matrix(NA_real_, nrow = nN, ncol = nJ)
  isc.mat  <- matrix(NA_integer_, nrow = nN, ncol = nJ)
  
  ## ---------------------------
  ## GRM simulation
  ## ---------------------------
  abil.mat <- matrix(abil, nrow = nN, ncol = nM)
  
  for (j in seq_len(nJ)) {
    
    a <- itempars[j, 1]
    b <- itempars[j, -1]
    
    if (nM == 1) {
      # Binary GRM (2PL model)
      eta <- a * (abil - b)
      p <- plogis(eta)
      
      isc.mat[, j]  <- as.integer(runif(nN) < p)
      true.mat[, j] <- p
      
    } else {
      # Cumulative logits for graded response model
      eta.mat <- a * sweep(abil.mat, 2, b, "-")
      p.mat <- plogis(eta.mat)
      
      # Category = number of thresholds exceeded
      u <- runif(nN)
      isc.mat[, j] <- rowSums(p.mat > u)
      
      # Expected score
      true.mat[, j] <- rowSums(p.mat)
    }
  }
  
  ## ---------------------------
  ## Return
  ## ---------------------------
  list(
    true = true.mat,
    isc  = isc.mat
  )
}
