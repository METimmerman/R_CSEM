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
    set.seed(50226)
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


#' Calculate model-implied composite SEM (cSEM) from an lme object
#'
#' This function computes the model-implied residual variance based on a
#' linear mixed-effects model (`lme`) with either a single `varExp`
#' variance structure or a combined `varComb` of two `varExp` functions,
#' and calculates the composite standard error of measurement (cSEM)
#' for a specified number of items.
#'
#' @param object An `lme` object (from `nlme::lme`) that includes a
#'   variance function (`varExp` or `varComb(varExp(...), varExp(...))`).
#' @param X_eval Numeric vector of predictor values for the first variance component.
#' @param X2_eval Optional numeric vector of predictor values for the second variance component.
#'   Required only if `object` uses `varComb` with two `varExp`s.
#' @param nJ Integer. Number of items used to compute the composite SEM.
#'
#' @return Numeric vector of cSEM values corresponding to each element of
#'   the input predictor(s).
#'
#' @details
#' The model-implied residual variance is computed as:
#' \deqn{Var(\epsilon) = \sigma^2 * \prod_i \exp(2 * \gamma_i * X_i)}
#' where `gamma_i` are the variance function coefficients and `X_i` are
#' the corresponding predictor values (`X_eval`, `X2_eval`).
#' The cSEM is then \eqn{cSEM = \sqrt{Var(\epsilon) * nJ}} assuming equal
#' residual variance across items.
#'
#' @examples
#' # Single gamma example
#' X_eval <- seq(0, 1, length.out = 5)
#' eval_resid_var_flex(model, X_eval, nJ = 10)
#'
#' # Two gamma example
#' X_eval <- seq(0, 1, length.out = 5)
#' X2_eval <- seq(-1, 1, length.out = 5)
#' eval_resid_var(model, X_eval, X2_eval, nJ = 10)
eval_resid_var <- function(object, X_eval, X2_eval = NULL, nJ = 10) {
  
  # 1. Residual standard deviation
  sigma_hat <- object$sigma
  
  # 2. Variance function coefficients (gamma)
  var_coefs <- coef(object$modelStruct$varStruct, unconstrained = FALSE)
  n_gamma <- length(var_coefs)
  
  # 3. Compute model-implied residual variance
  if(n_gamma == 1) {
    # Single variance function
    resid_var <- sigma_hat^2 * exp(2 * var_coefs[1] * X_eval)
    
  } else if(n_gamma == 2) {
    # Two variance functions combined
    if(is.null(X2_eval)) stop("X2_eval must be provided for a two-component variance function")
    if(length(X_eval) != length(X2_eval)) stop("X_eval and X2_eval must be the same length")
    
    resid_var <- sigma_hat^2 * exp(2 * var_coefs[1] * X_eval + 2 * var_coefs[2] * X2_eval)
    
  } else {
    stop("Variance function with more than 2 components not supported")
  }
  
  # 4. Composite SEM
  cSEM_2 <- sqrt(resid_var * nJ)
  
  return(cSEM_2)
}

#' Plot normalized level-1 residuals versus fitted values from an lme object
#'
#' @param object An `lme` object (from `nlme::lme`) that includes a
#'   variance function (`varExp` or `varComb(varExp(...), varExp(...))`).
#' @param Xplus Numeric vector of predictor values used for coloring points.
#'
#' @return A ggplot object showing residuals versus fitted values.
plot_resid_var <- function(object, Xplus, title = "Residuals vs Fitted") {
  
  fv  <- fitted(object, level = 1)          # conditional fitted values
  res <- resid(object, type = "normalized") # level-1 normalized residuals
  
  if (length(fv) != length(res) || length(res) != length(Xplus)) {
    stop("fv, res, and Xplus must have the same length")
  }
  
  df <- data.frame(
    fitted   = fv,
    residual = res,
    Xplus    = Xplus
  )
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = fitted, y = residual, color = Xplus)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      x = "Fitted values (including random intercepts)",
      y = "Normalized residuals",
      title = title
    )
  
  return(p)
}


