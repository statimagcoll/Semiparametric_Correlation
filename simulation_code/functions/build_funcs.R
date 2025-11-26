### This file contains the functions used to set up the simulations and control the value of MAPA ###

### Create needed directories ###
check_dirs = function(path, rhos, ns){
  # Loop through each x value
  for (r in rhos) {
    # Create the path for the "rho_x" directory
    rho_dir <- file.path(path, paste0("rho_", r*100))
    
    # Check if the "rho_x" directory exists, if not, create it
    if (!dir.exists(rho_dir)) {
      dir.create(rho_dir, recursive = TRUE)
    }
    
    # Loop through each X value
    for (n in ns) {
      # Create the path for the "nX" directory within the "rho_x" directory
      n_dir <- file.path(rho_dir, paste0("n", n))
      
      # Check if the "nX" directory exists, if not, create it
      if (!dir.exists(n_dir)) {
        dir.create(n_dir, recursive = TRUE)
      }
    }
  }
}

### Functions to estimate the build model ###

# Split data
# prop = proportion of data to use for the build/test model
split_model = function(dt, prop = 0.2){
  # get number of participants to put in build/test model
  n =  round(nrow(dt) * prop)
  # sample the build/test participants
  inds = sample(nrow(dt), size = n)
  # separate data into build/test and eval/train data
  d1 = dt[inds,]
  d2 = dt[-inds,]
  return(list(d1 = d1, d2 = d2))
}


# Functions for numerical integration
int_f = function(x, dist){
  x*dlogspline(x, dist)
} 

int_fsq = function(x, dist){
  x^2*dlogspline(x, dist)
}

int_fthird = function(x, dist){
  x^3*dlogspline(x, dist)
}

int_ffourth = function(x, dist){
  x^4*dlogspline(x, dist)
}

# Error distributions
runif_target = function(n, rho_true, v_xi){
  runif(n, min = -sqrt(3*(1/rho_true^2 - 1)*v_xi), max = sqrt(3*(1/rho_true^2 - 1)*v_xi))
}

v_uniform <- function(rho_true, v_xi){
  1/12 * (2*sqrt(3*(1/rho_true^2 - 1)*v_xi))^2
}


# function to get correlation for ridge regression
# The limit is the OLS fit
cor_mu = function(X, xi, lambda, v_res, intercept = T){
  #browser()
  
  # add intercept 
  if(intercept){
    X = cbind(1, X)
  }
  
  XtX = t(X) %*% X
  XtX_linv = solve(XtX + diag(lambda, nrow = ncol(X)))
  H =  X %*% XtX_linv %*% t(X)
  
  mu = H %*% (as.matrix(xi))
  rho_0 = (mean(xi*mu)-mean(xi)*mean(mu))/sqrt((var(xi)+v_res)*var(mu))
  
  # correlations
  c(cor_mu_xi = cor(mu, xi), cor_mu_Y = rho_0, lambda = lambda)
}

# to get mu values in the population
get_mu = function(X, xi, lambda, intercept = T){
  if(intercept){
    X = cbind(1, X)
  }
  XtX = t(X) %*% X
  XtX_linv = solve(XtX + diag(lambda, nrow = ncol(X)))
  H =  X %*% XtX_linv %*% t(X)
  H %*% (as.matrix(xi))
}

### Old error distributions
## Edited 2/27/25 to center at 0
# function to return a logspline dist of residuals with variance corresponding to a desired rho_0 parameter
# res = residuals from build model
# mu = mu(X) values in the eval dataset
# target_rho = desired value for rho_0 [0,1]
logspline_targetrho = function(res, mu, target_rho){
  # get initial logspline distribution of the residuals
  initial_logspline = logspline(res)
  ex_res = integrate(int_f, lower = -Inf, upper = Inf, dist = initial_logspline)$value
  exsq_res =  integrate(int_fsq, lower = -Inf, upper = Inf, dist = initial_logspline)$value
  # get current variance of the distribution
  v_res = exsq_res -  ex_res^2
  # variance of mu(X)
  v_mu = var(mu)
  # what the variance of the error distribution should be for the desired rho
  new_v = v_mu/(target_rho^2) - v_mu
  # scale residuals
  
  # edited 4/9/25 - had shifted after scaling before (still took off the mean later so it shouldn't have caused a major issue)
  res = res - mean(res)
  
  m = sqrt(new_v/v_res)
  m_res = res*m
  
  # create new logspline object based on scaled residuals
  new_ls = logspline(m_res)
  return(new_ls)
}

# to try normal distribution
normal_targetrho = function(res, mu, target_rho){
  #browser()
  # get current variance of the distribution
  v_res = var(res)
  # variance of mu(X)
  v_mu = var(mu)
  # what the variance of the error distribution should be for the desired rho
  new_v = v_mu/(target_rho^2) - v_mu
  
  # hack - generate large number from this normal distribution and then use logspline so we don't have to change the other code
  s = rnorm(1000000, 0, sd = sqrt(new_v))
  
  # create new logspline object based on scaled residuals
  new_ls = logspline(s)
  return(new_ls)
  #return(v = new_v)
}


# Y is generated as muX + sample from logspline distribution - mean(logspline distribution)

# Get ground truth rho from a build model and logspline residual distribution
# mu = predicted y values from build model in build data (aka true muX)
# ls_obj = logspline object from residuals
get_rho_ls = function(mu, ls_obj, ...){
  # calculate variance of the residuals
  ex_res = integrate(int_f, lower = -Inf, upper = Inf, dist = ls_obj)
  exsq_res =  integrate(int_fsq, lower = -Inf, upper = Inf, dist = ls_obj)
  v_res = exsq_res$value - ex_res$value^2
  rho = sqrt(var(mu)/(var(mu) + v_res))
  # return rho and mean of the distribution
  return(list(mean_res = ex_res$value, rho = rho, v_res = v_res))
}



# higher moments of the error distribution
get_error_third = function(dist){
  integrate(int_fthird, lower = -Inf, upper = Inf, dist = dist)$value
}

get_error_fourth = function(dist){
  integrate(int_ffourth, lower = -Inf, upper = Inf, dist = dist)$value
}


# To not generate completely new error for each setting of rho_true for rho_2, just scale the error
scale_error_rho_2 = function(error_original, v_res_original, v_mu, new_rho_true){
  sqrt((v_mu/new_rho_true^2 - v_mu)/v_res_original)*error_original
} 


