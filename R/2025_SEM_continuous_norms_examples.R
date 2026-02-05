# Estimating the SEM using multilevel modeling


# libraries and source files
  library(MASS)
  library(tidyverse)
  library(lme4)
  library(nlme)
  
  
# Read dedicated function (adjust the path):  
#  source("~/02 Onderzoek/Reliability of continuous norms/R codes/2025_SEM_continuous_norms_functions.R", echo=TRUE)

  library(here)  
  here::i_am("R/2025_SEM_continuous_norms_functions.R")
  source(
    here::here("R", "2025_SEM_continuous_norms_functions.R"),
    echo = TRUE
  )
  
# Deel I: Simulations using CTT philosophy --------------------------------


# I used this part to explore the role of the compound symmetry assumption in
# estimating the SEM using the random intercept model with fixed item effects
  

# Define the CTT population (Note: this routine creates variance-covariance matrices for the true scores, 
# error scores, and total scores given specific constraints such as desired reliability. I used this routine 
# to easily simulate data under different error-variance structures, while knowing the population values of the estimands)
  CTT_pop <- get_TEX_matrices(
    n_items_in =  10, 
    pop_reliability =  0.7, 
    covTiTj = .4, # covariantie tussen de item true scores 
    lbE = 1,  # ondergrens item error varianties (voor herschaling) 
    ubE = 5   # ondergrens item error varianties (voor herschaling) 
  )

  
# Simulate continuous CTT data (directly drawn from multivariate normal distribution)    
  n_persons <- 1000
  X <- mvrnorm(n_persons, mu = rep(0, nrow(CTT_pop$sigmaX)), Sigma = CTT_pop$sigmaX)    

# Compute sample statistics
  alpha <- psychometric::alpha(X)
  SEMsq <- sum(var(X)) * (1-alpha)
  c(alpha, SEMsq)

  
# Reshape data set to long format (needed for estimating random intercept models)
  Xwide <- as.data.frame(X)
  colnames(Xwide) <- paste0("I", 1:ncol(X))
  Xwide$Person <- rownames(X) %||% 1:nrow(X)
  
  Xlong <- pivot_longer(
    Xwide,
    cols = -Person,
    names_to = "Item",
    values_to = "Score"
  )
  head(Xlong)

  
# Estimating the SEM (unconditional)
  
  n_items <- (ncol(Xwide) - 1)
  
  # via lmer (from lme4)    ##MET: What model exactly??! 
  model <- lmer(Score ~ Item + (1 | Person), data = Xlong, REML = TRUE)
  summary(model)
  res_term <- as.data.frame(VarCorr(model))[2, 4]
  SEM_lme <- res_term * n_items
  all.equal(SEM_lme, SEMsq)

  
 
  # via nlme compound symmetry imposed
  model_cs <- lme(
    fixed = Score ~ Item,
    random = ~1 | Person,
    correlation = corCompSymm(form = ~1 | Person),
    data = Xlong
  )
  res_term <- as.numeric(VarCorr(model_cs)[2, 1])
  SEM_lme <- res_term * n_items
  all.equal(SEM_lme, SEMsq)
  
  
  
  # via nlme fully multivariate
  model_fmv <- lme(
    fixed = Score ~ Item,
    random = ~1 | Person,
    correlation = corSymm(form = ~1 | Person),
    data = Xlong
  )
  res_term <- as.numeric(VarCorr(model_fmv)[2, 1])
  SEM_lme <- res_term * n_items
  all.equal(SEM_lme, SEMsq)
  
  anova(model_cs, model_fmv)


      
# Conditional Approaches --------------------------------------------------

# Below you find some tries for estimating the SEMs conditional on sum scores.
# For this, I generated data under the 2PLM, otherwise conditioning should be based on binning. 
  
  
# Conditional approach (via IRT)
  
  # Define item parameters
  n_items <- 10
  ipars <- matrix(NA, n_items, 2)
  ipars[, 1] <- runif(n_items, 1, 3)
  ipars[, 2] <- runif(n_items, -2, 2)
  

# Explore the population values:
  
  # Approximate population values using a very large simulated sample (assuming standard normal theta)
  X <- simulateGRM(rnorm(1000000), ipars)
  appr_pop_sqSEM <- mean((rowSums(X$isc) - rowSums(X$true))^2)
  sqrt(appr_pop_sqSEM)
  
  
  # Condition SEM (Take care; see remarks other conceptual paper!!!!!)
  Xplus <- rowSums(X$isc)
  cSEMs <- matrix(NA, 11, 1)
  for (x in 0:10){
    sel <- Xplus == x
    cSEMs[x + 1, 1] <- mean((x - rowSums(X$true)[sel])^2)
  }
  plot(0:10, cSEMs[, 1], type="b", pch = 16)


# Apply the methodology to a smaller simulated data setje.     
  
# simulate data setje using 2PLM/GRM for the examples (wide format)    
  n_persons <- 150
  X <- simulateGRM(rnorm(n_persons), ipars)$isc
  Xwide <- as.data.frame(X)
  colnames(Xwide) <- paste0("I", 1:ncol(X))
  Xwide$Person <- rownames(X) %||% 1:nrow(X)
  
  # Compute sample statistics
  alpha <- psychometric::alpha(X)
  SEMsq <- sum(var(X)) * (1-alpha)
  c(alpha, SEMsq, appr_pop_sqSEM)
  c(alpha, sqrt(SEMsq), sqrt(appr_pop_sqSEM))
  

# Estimating the unconditional SEM for IRT generated data using different set-ups

  Xlong <- pivot_longer(
    Xwide,
    cols = -Person,
    names_to = "Item",
    values_to = "Score"
  )
  head(Xlong)
  
  
  # via lmer (from lme4) - random intercept model, persons random, items fixed
  model <- lmer(Score ~ Item + (1 | Person), data = Xlong, REML = TRUE)
  res_term <- as.data.frame(VarCorr(model))[2, 4]
  SEM_lme <- res_term * n_items
  all.equal(SEM_lme, SEMsq)  # numeric laten zien dat SEM via Alpha = SEM via lmer
  
  
  # via nlme compound symmetry imposed
  model_cs <- lme(
    fixed = Score ~ Item,
    random = ~1 | Person,
    correlation = corCompSymm(form = ~1 | Person),
    data = Xlong
  )
  res_term <- as.numeric(VarCorr(model_cs)[2, 1])
  SEM_lme_cs <- res_term * n_items
  all.equal(SEM_lme_cs, SEMsq)
  
  
  # via nlme fully multivariate (only for comparison; becomes easily slow)
  model_fmv <- lme(
    fixed = Score ~ Item,
    random = ~1 | Person,
    correlation = corSymm(form = ~1 | Person),
    data = Xlong
  )
  res_term <- as.numeric(VarCorr(model_fmv)[2, 1])
  SEM_lme_fmv <- res_term * n_items
  all.equal(SEM_lme_fmv, SEMsq)
  
  anova(model_cs, model_fmv)
  
  
# Estimating the SEM conditional on the sum score (compound symmetry only)

  # extend the data file with the summed scores and its square
  Xplus <- rowSums(Xwide[, 1:n_items])
  #Xplus <- as.factor(Xplus)
  Xwide  <- Xwide %>% mutate(Xplus = Xplus, sqXplus = Xplus^2)
  head(Xwide)
  
  Xlong <- pivot_longer(
    Xwide,
    cols = 1:n_items,
    names_to = "Item",
    values_to = "Score"
  )
  head(Xlong)
  
  
# Try it out for a selected sum score group:  
  
  # Select sum score group
  Xlong_select  <- Xlong %>% filter(Xplus == 6)
  
  # fit the model in the selected group (notice the warning messages; make)
  model_cs <- lme(
    fixed = Score ~ Item,
    random = ~ 1 | Person,
    correlation = corCompSymm(form = ~1 | Person),
    data = Xlong_select
  )
  summary(model_cs)
  


# Use summed score as covariate for residual varaiance
  Xlong$Xplus <- as.factor(Xlong$Xplus)
  model_cs <- lme(
    fixed = Score ~ Item,
    random = ~ 1 | Person,
    #correlation = corCompSymm(form = ~1 | Person),
    weights = varIdent(form = ~ 1| Xplus),
    data = Xlong
  )
  summary(model_cs)
  
  
# explore the results (code is rather messy because the objects created by lme have a weird structure)
  foo <- summary(model_cs)
  goo <- foo$modelStruct$varStruct
  coef(goo)
  goo <- c(1, 1 + as.numeric(foo$modelStruct$varStruct))
  
  # Extract estimated residual variances per observation
  base_var <- sigma(model_cs)^2 
  cSEMS <- sqrt(base_var * goo^2 * n_items)
  io <- c (      9 ,       4,        6,        7,        8,        2,        5,        1,        0,        3,       10)  # order of cSEMS in the output (must fix it)
  plot(io, cSEMS, pch = 16, ylim=c(0, 2))
  abline(h = sqrt(SEMsq))
  
  
  
  # Use summed score as covariate for residual varaiance (via exponential funciton; for some reason doesn'work)
  Xlong$Xplus <- as.numeric(Xlong$Xplus)
  
  model_cs <- lme(
    fixed = Score ~ Item,
    random = ~ 1 | Person,
    #correlation = corCompSymm(form = ~1 | Person),
    weights = varExp(form = ~ Xplus), # it seems you can only have one covariate predictor...
    data = Xlong
  )
  summary(model_cs)
  
  base_var <- sigma(model_cs)^2 
  sqrt(exp(0:10 * -0.04299932) *  base_var)
  
  
  ## heteroskedastic residual model. aanbevolen door ChatGPT obv mijn aanwijzingen.
  library(nlme)
  
  ##Model mét random intercept (kan op rand-oplossing uitkomen)
  m1 <- lme(
    y ~ 1,
    random = ~ 1 | subject,
    weights = varExp(form = ~ X),
    data = dat,
    control = lmeControl(msMaxIter = 200)
  )
  
  ##Model zonder random intercept (τ² ≡ 0)
  m0 <- gls(
    y ~ 1,
    weights = varExp(form = ~ X),
    data = dat
  )
  
  