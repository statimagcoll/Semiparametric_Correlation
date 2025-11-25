### Convenience functions for simulations ###
# logit and expit functions
logit = function(x){
  log(x/(1-x))
}

expit = function(x){
  1/(1+exp(-x))
}

## Functions to get predictions (based on model type) in a way that works nicely with other written functions ##

get_pred = function(mod, dt, ...){
  UseMethod("get_pred")
}

get_pred.default = function(mod, dt, ...){
  predict(mod, dt)
}

get_pred.ranger = function(mod, dt, ...){
  predict(mod, dt)$predictions
}

get_pred.glmnet = function(mod, dt, lambda, ...){
  predict(mod, dt, s = lambda)
}

get_pred.cv.glmnet = function(mod, dt, ...){
  dt = as.matrix(dt)
  # use best cvm error lambda
  muhat = predict(mod, dt, s = mod$lambda.min, ...)
  if (length(unique(muhat)) == 1){ # if using the largest lambda, all the predictions will be the same
    muhat = predict(mod, dt, mod$lambda[2]) # use the second largest lambda
  }
  return(muhat)
}

# tuned Ranger function
get_pred.WrappedModel = function(mod, dt, ...){
  dt = as.data.frame(dt)
  #colnames(dt) = paste0("x.",colnames(dt))
  muhat = predict(mod, newdata = dt)$data$response
  return(muhat)
}


# "alpha" is a defined class that basically just returns the predictions as already computed (used for rho_2 predictions and another type of simulation I had been doing)
get_pred.alpha = function(mod, ...){
  mod
}

# function just to get the rho_2 predictions to work with the overall estimation function (cf_corr_if below) I had already written
rho_2_func = function(preds, ...){
  class(preds) = "alpha"
  preds
}
