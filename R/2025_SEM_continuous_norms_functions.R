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
  
  
  #stop("incorrect method specificied")    
}

simulateGRM <- function (abil = NULL, itempars = NULL) 
{
  if (is.null(itempars)) 
    stop("Item parameters not specified.")
  if (is.null(abil)) 
    stop("Person parameters not specified.")
  if (!is.matrix(itempars)) 
    stop("Item parameters should be provided as a matrix.")
  if (!is.vector(abil)) 
    stop("Person parameters should be provided as a vector.")
  nN <- length(abil)
  nJ <- dim(itempars)[1]
  nM <- dim(itempars)[2] - 1
  abil.mat <- matrix(abil, nN, nM)
  isc.mat <- matrix(NA, nN, nJ)
  true.mat <- matrix(NA, nN, nJ)
  prob <- matrix(NA, nN, nJ)
  for (.j in 1:nJ) {
    if (nM == 1) {
      logit <- itempars[.j, 1] * (abil - itempars[.j, 2])
      true.vec <- exp(logit)/(1 + exp(logit))
      isc.vec <- 1 * (true.vec - runif(nN, 0, 1) >= 0)
      true.mat[, .j] <- true.vec
      isc.mat[, .j] <- isc.vec
    }
    else {
      prob.mat <- itempars[.j, 1] * (sweep(abil.mat, 2, 
                                           itempars[.j, -1], "-"))
      prob.mat <- exp(prob.mat)/(1 + exp(prob.mat))
      isc.vec <- apply(sweep(prob.mat, 1, runif(nN, 0, 
                                                1), "-") >= 0, 1, sum)
      true.vec <- apply(prob.mat, 1, sum)
      isc.mat[, .j] <- isc.vec
      true.mat[, .j] <- true.vec
    }
  }
  return(list(true = true.mat, isc = isc.mat))
}
