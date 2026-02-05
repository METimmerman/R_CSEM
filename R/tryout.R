# Fitted values (marginal / inclusief random effects)
#fv <- fitted(model, level = 0)  # fixed-effects only
fv <- fitted(model, level = 1)  # inclusief random intercepts

# Residuals (standaard = conditional)
res <- resid(model, type = "normalized")  # of "pearson"

plot(fv, res,
     xlab = "Fitted values - incl. random intercepts",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "blue", lty = 2)

library(ggplot2)
ggplot(Xlong, aes(x = fv, y = res, color = Xplus)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted")



##Let's go Bayesian

library(brms)
brm_model <- brm(
  bf(Score ~ Item + (1|Person), sigma ~ s(Xplus)),
  data = Xlong,
  family = gaussian()
)

# Samenvatting van posterior
summary(brm_model)

# Rhat moet ~1 zijn voor alle parameters
# Neff = effectieve sample size

library(bayesplot)

pp_check(brm_model)  # standaard histogram / density overlay
pp_check(brm_model, type = "bars")  # discrete outcomes

library(ggplot2)

res <- residuals(brm_model, type = "pearson")[, "Estimate"]
fitted_vals <- fitted(brm_model)[, "Estimate"]  # posterior mean

df <- data.frame(fitted = fitted_vals, residuals = res, Xplus = Xlong$Xplus)

ggplot(df, aes(x = fitted, y = residuals, color = Xplus)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted")


X_grid <- data.frame(
  Xplus = seq(min(Xlong$Xplus), max(Xlong$Xplus), length.out = 100)
)


##Plot
X_grid <- data.frame(
  Xplus = seq(min(Xlong$Xplus), max(Xlong$Xplus), length.out = 100),
  Item  = Xlong$Item[1]   # any existing level is fine
)

sigma_grid <- posterior_epred(
  brm_model,
  dpar = "sigma",
  newdata = X_grid,
  re_formula = NA
)


sigma_mean  <- apply(sigma_grid, 2, mean)
sigma_lower <- apply(sigma_grid, 2, quantile, 0.025)
sigma_upper <- apply(sigma_grid, 2, quantile, 0.975)

plot(
  X_grid$Xplus, sigma_mean,
  type = "l",
  ylim = range(c(sigma_lower, sigma_upper)),
  xlab = "Xplus",
  ylab = expression(sigma),
  main = expression("Residual SD as a smooth function of Xplus")
)
lines(X_grid$Xplus, sigma_lower, lty = 2)
lines(X_grid$Xplus, sigma_upper, lty = 2)

