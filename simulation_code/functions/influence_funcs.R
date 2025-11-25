### Empirical influence functions and OS estimators ###

## Plug-In Estimators  ##

# Pearson plug-in can be obtained via cor() or cor.test()

# Simplified plug-in version
corr_simp = function(muhat, y){
  return(c(simp_hat = sqrt(var(muhat)/var(y))))
}


## Influence Functions ##
# Functions for empirical influence function (returns vector)

# Fully derived from either Pearson or Simplified parameter
if_full = function(muhat, y){
  # components
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  # influence function
  return((2*y*(muhat - my) - muhat^2 + my^2)/(2*sqrt(vmuhat*vy)) -
           (sqrt(vmuhat)*(y - my)^2/(2*(vy)^(3/2))))
  }


# Mix of influence function components for mux components and empirical for y components (starting from simplified)
# returns vector
# changed 3/14/25
if_mix_s = function(muhat, y){
  # components
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  
  # influence funcion
  return((2*y*(muhat - my) - muhat^2 + my^2 - vmuhat)/(2*sqrt(vmuhat*vy)))
}

# IF(logit(rho_0))
if_logit = function(muhat, y){
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  # fixed 3/5/25
  return((sqrt(vy)/sqrt(vmuhat) + sqrt(vy)/(sqrt(vy) - sqrt(vmuhat))) * if_full(muhat, y))
}

# IF(logit(rho_0^2))
# fixed 10/8/24
if_logit_sq = function(muhat, y){
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  
  return(((y-my)^2*(vy-vmuhat) - vy*(y-muhat)^2)/(vmuhat*(vy-vmuhat)))
}


# IF of Fisher transformed parameter
if_fisher = function(muhat, y){
  return((1 / (1 - cor(muhat, y)^2)) * if_full(muhat, y))
}

# IFs reparameterized in terms of rho_0 (psi)
# psi_hat is the plug-in estimate
if_full_psi = function(muhat, y, psi_hat){
  my = mean(y)
  vy = var(y)
  vmuhat=var(muhat)
  return((2*y*(muhat-my) - muhat^2 + my^2)/(2*psi_hat*vy) - (psi_hat*(y-my)^2)/(2*vy))
  
}

if_mixed_psi = function(muhat, y, psi_hat){
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  return((2*y*(muhat-my) - muhat^2 + my^2)/(2*psi_hat*vy) - psi_hat/2)
  
}

if_logit_psi = function(muhat, y, psi_hat){
  vy = var(y)
  vmuhat = var(muhat)
  (-1/(2*psi_hat^2*(1-psi_hat)))*((y-muhat)^2 - (1-psi_hat^2)*(y-mean(y))^2)/var(y)
  
}

if_logit_sq_psi = function(muhat, y, psi_hat){
  my = mean(y)
  vy = var(y)
  vmuhat = var(muhat)
  covymuhat = cov(y, muhat)
  
  p = covymuhat/sqrt(vmuhat*vy)
  (2*p)/(p^2*(1-p^2)) * if_full(muhat, y)
  
}

if_fisher_psi = function(muhat, y, psi_hat){
  0.5*((y-muhat)^2/(var(y)*(1-psi_hat^2)) - 1)
}

# Estimating Equation/conditional variance framework (parameterized in terms of theta = E[V(Y|X)]/V(Y))
# From Edward's notes: corr(Y, mu) = sqrt(1-theta) where theta = E[V(Y|X)]/V(Y)
# theta_hat = Pn((Y-muhat)^2)/Pn((Y-Ybar))
theta_est = function(muhat, y){
  return(mean((y-muhat)^2)/mean((y-mean(y))^2))
}

# for one-step estimator based on this parameterization
if_theta = function(muhat, y){
  plug_in = sqrt(1-mean((y-muhat)^2)/mean((y-mean(y))^2))
  (-1/(2*plug_in*var(y)))*((y-muhat)^2 - (1-plug_in^2)*(y-mean(y))^2)
}

# IF(logit(sqrt(1-theta))) based on Edward's notes
if_logit_theta = function(muhat, y){
  theta_hat = theta_est(muhat, y)
  return(1/(2*(1-theta_hat)*(1 - sqrt(1 - theta_hat))) * (((y-muhat)^2 - theta_hat*(y-mean(y))^2)/var(y)))
}

# IF(logit(1-theta)) from Edward's notes
if_logit_sq_theta = function(muhat, y){
  theta_hat = theta_est(muhat, y)
  return(-1/(theta_hat*(1-theta_hat)) * ((y-muhat)^2 - theta_hat*(y-mean(y))^2)/var(y))
}

# IF(atanh(sqrt(1-theta)))
if_fisher_theta = function(muhat, y){
  theta_hat = theta_est(muhat, y)
  return((-1/(4*sqrt(1-theta_hat)))*((1/(1+sqrt(1-theta_hat))) + 1/(1-sqrt(1-theta_hat)))*((y-muhat)^2/var(y) - (y-mean(y))^2*theta_hat/var(y)))
}


## New Estimators ##
# Some estimators are computed directly in the main "cf_corr_if" function

# logit transformed estimator (given plug in value and IF vector)
logit_pe = function(muhat, y, plug_in, if_emp){
  if(plug_in < 0) plug_in = 0
  if(plug_in > 1) plug_in = 1
  expit(logit(plug_in) + mean(if_emp))
}

# squared logit transformed estimator (given plug in value and IF vector)
# changed 10/8/24 to truncate for plug in outside of [-1, 1]
logit_sq_pe = function(muhat, y, plug_in, if_emp){
  if(plug_in < -1) plug_in = -1
  if(plug_in > 1) plug_in = 1
  sqrt(expit(logit(plug_in^2) + mean(if_emp)))
}

# logit point estimate and confidence interval
logit_ci = function(muhat, y, plug_in, if_emp, cl = 0.975){
  pe = logit_pe(muhat, y, plug_in, if_emp)
  bounds = expit(logit(pe) + c(1, -1)*qnorm(cl)*sqrt(var(if_emp)/length(muhat)))
  return(lCI = bounds[1], pe = pe, uCI = bounds[2])
}

# squared logit point estimate and confidence interval
logit_sq_ci = function(muhat, y, plug_in, if_emp, cl = 0.975){
  pe = logit_sq_pe(muhat, y, plug_in, if_emp)
  bounds = sqrt(expit(logit(pe^2) + c(1, -1)*qnorm(cl)*sqrt(var(if_emp)/length(muhat))))
  return(lCI = bounds[1], pe = pe, uCI = bounds[2])
}

# fisher transformed estimator (given plug in value and IF vector)
fisher_pe = function(muhat, y, plug_in, if_emp){
  if(plug_in < -1) plug_in = -1
  if(plug_in > 1) plug_in = 1
  tanh(atanh(plug_in) + mean(if_emp))
}

# fisher transformed point estimate and confidence interval
fisher_ci = function(muhat, y, plug_in, if_emp, cl = 0.975){
  pe = fisher_pe(muhat, y, plug_in, if_emp)
  bounds = tanh(atanh(pe) + c(1, -1)*qnorm(cl)*sqrt(var(if_emp)/length(muhat)))
  return(lCI = bounds[1], pe = pe, uCI = bounds[2])
}

# estimating equation estimator (CT) based on sqrt(1-theta) and IF of theta
ee_theta = function(muhat, y){
  sqrt(1-theta_est(muhat, y))
}
