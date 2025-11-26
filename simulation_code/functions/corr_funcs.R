### Functions to compute the OS estimator and compare to Pearson in simulations ###

# User-friendly function to compute one-step estimator without having to explicitly compute the influence function
# Inputs: y, yhat, fold (vectors), confidence level (scalar)
# Outputs: List with: [1] vector of cor est plus confidence bounds, [2] data.frame with y, yhat, fold, and influence function estimates
# Only compute squared logit, leave other transformations for more manual construction for now
# fold argument should 
corOS <- function(yhat, y, fold, cl = 0.95){
  # compute influence function estimates
  out_df <- data.frame(y = y, yhat = yhat, fold = fold) %>%
    group_by(fold) %>%
    mutate(influence = if_logit_sq(yhat, y))
  
  # compute point estimate and bounds
  fold_summary <-  out_df %>%
    group_by(fold) %>%
    summarise(mean_if = mean(influence),
              pearson = cor(yhat, y),
              fold_cor_os = sqrt(expit(logit(pearson^2) + mean_if)),
              .groups = "drop")
  
  # final point estimate
  os <- mean(fold_summary$fold_cor_os)
  
  # CI bounds
  lb <- sqrt(expit(logit(os^2) - qnorm(cl + (1-cl)/2)*sqrt(var(out_df$influence)/length(y))))
  ub <- sqrt(expit(logit(os^2) + qnorm(cl + (1-cl)/2)*sqrt(var(out_df$influence)/length(y))))
  
  return(list(results = c(estimate = os, lb = lb, ub = ub), data = out_df, cl = cl))
    
}

# Function to compute one-step
# muhat: predicted y values
# y: actual y
# transform: estimator transformation (squared logit ("q") by default). Options: "none", "no" for no transform; "l" or "logit" for logit; "q" or "squared logit" for squared logit, "z", "f" or "fisher" for Fisher
cor_os = function(muhat, y, transform = c("q")){
  p = cor(muhat, y)
  if (transform %in% c("none", "no")){
    r = p + mean(if_full(muhat, y))
  }
  if (transform %in% c("l", "logit")){
    p = ifelse(p < 0, 0, ifelse(p > 1, 1, p))
    r = expit(logit(p) + mean(if_logit(muhat, y)))
  }
  if (transform %in% c("q", "squared logit")){
    p = ifelse(p < -1, -1, ifelse(p > 1, 1, p))
    r = sqrt(expit(logit(p^2) + mean(if_logit_sq(muhat, y)))) # fixed from if_full to if_logit_sq 08/12/2025
  }
  if (transform %in% c("z", "f", "fisher")){
    p = ifelse(p < -1, -1, ifelse(p > 1, 1, p))
    r = tanh(atanh(p) + mean(if_fisher(muhat, y)))
  }
  
  return(r)
}

# Function to compute difference between two one-step estimators (os2 - os1) (with confidence interval)
# os2: one-step estimate to subtract from
# os1: one-step estimate to subtract
# influence2: vector of influence functions for parameter 2 (corresponding to os2) for entire sample
# influence1: vector of influence functions for parameter 1 (corresponding to os1) for entire sample
# transform: estimator transformation (squared logit ("q") by default). Options: "none", "no" for no transform; "l" or "logit" for logit; "q" or "squared logit" for squared logit, "z", "f" or "fisher" for Fisher
# cl: confidence level
cor_os_diff = function(os2, os1, influence2, influence1, transform = "q", cl = 0.95){
  if(length(influence1) != length(influence2)){
    stop("Influence function estimates must come from the same sample")
  }
  n = length(influence1)
  
  # point estimate: difference
  point_est = os2 - os1
  
  # variance: v(os1) + v(os2) - 2*cov(os1, os2)
  # computed as sum of the empirical variances of the influence functions minus 2 times the empirical covariance of the influence function estimates
  v1 = var(influence1)
  v2 = var(influence2)
  cov12 = cov(influence1, influence2)
  if (transform %in% c("none", "no")){
    var_factor = function(rho_hat) 1
  }
  if (transform %in% c("l", "logit")){
    var_factor = function(rho_hat){
      rho_hat*(1-rho_hat)
    }
  }
  if (transform %in% c("q", "squared logit")){
    var_factor = function(rho_hat){
      0.5*rho_hat*(1-rho_hat^2)
    }
  }
  if (transform %in% c("z", "f", "fisher")){
    var_factor = function(rho_hat){
      (1-rho_hat^2)
    }
  }
  
  v_diff = v1*var_factor(os1)^2 + v2*var_factor(os2)^2 - 2*cov12*var_factor(os1)*var_factor(os2)
  
  lower = point_est - qnorm(cl+(1-cl)/2)*sqrt(v_diff/n)
  upper = point_est + qnorm(cl+(1-cl)/2)*sqrt(v_diff/n)
  
  return(c("estimate" = point_est, "lCI" = lower, "uCI" = upper))
}

# Function to compute ratio of Cohen's f^2 via two one-step (squared logit) estimators (with confidence interval)
# os2: one-step estimate (numerator)
# os1: one-step estimate (denominator)
# influence2: vector of influence functions for parameter 2 (corresponding to os2) for entire sample
# influence1: vector of influence functions for parameter 1 (corresponding to os1) for entire sample
# cl: confidence level
cor_os_ratio = function(os2, os1, influence2, influence1, cl = 0.95){
  #browser()
  if(length(influence1) != length(influence2)){
    stop("Influence function estimates must come from the same sample")
  }
  n = length(influence1)
  
  # point estimate: exponeniated log ratio
  point_est = exp(logit(os2^2)-logit(os1^2))
  
  # variance: v(os1) + v(os2) - 2*cov(os1, os2)
  # computed as sum of the empirical variances of the influence functions minus 2 times the empirical covariance of the influence function estimates
  v1 = var(influence1)
  v2 = var(influence2)
  cov12 = cov(influence1, influence2)

  v_diff = v1 + v2 - 2*cov12
  
  lower = exp(log(point_est) - qnorm(cl+(1-cl)/2)*sqrt(v_diff/n))
  upper = exp(log(point_est) + qnorm(cl+(1-cl)/2)*sqrt(v_diff/n))
  
  return(c("estimate" = point_est, "lCI" = lower, "uCI" = upper))
}



library(boot)

# Function to add the influence functions based on requested types
add_influence_functions <- function(out, muhat, y, types, p_hat, s_hat) {
  #browser()
  if (any(c("pf", "sf", "ct") %in% types)) out$if_full_hat <- if_full(muhat, y)
  if (any(c("pm", "sm") %in% types)) out$if_mixed_hat <- if_mix_s(muhat, y)
  if (any(c("pl", "sl") %in% types)) out$if_logit_hat <- if_logit(muhat, y)
  if (any(c("pq", "sq") %in% types)) out$if_logit_sq_hat <- if_logit_sq(muhat, y)
  if (any(c("pz", "sz") %in% types)) out$if_fisher_hat <- if_fisher(muhat, y)
  if ("cl" %in% types) out$if_logit_theta_hat <- if_logit_theta(muhat, y)
  if ("cq" %in% types) out$if_logit_sq_theta_hat <- if_logit_sq_theta(muhat, y)
  if ("cz" %in% types) out$if_fisher_theta_hat <- if_fisher_theta(muhat, y)
  if ("pfr" %in% types) out$if_full_rep_hat <- if_full_psi(muhat, y, p_hat)
  if ("pmr" %in% types) out$if_mixed_rep_hat <- if_mixed_psi(muhat, y, p_hat)
  if ("plr" %in% types) out$if_logit_p_rep_hat <- if_logit_psi(muhat, y, p_hat)
  if ("slr" %in% types) out$if_logit_s_rep_hat <- if_logit_psi(muhat, y, s_hat)
  if ("pqr" %in% types) out$if_logit_sq_p_rep_hat <- if_logit_sq_psi(muhat, y, p_hat)
  if ("sqr" %in% types) out$if_logit_sq_s_rep_hat <- if_logit_sq_psi(muhat, y, s_hat)
  if ("pzr" %in% types) out$if_fisher_p_rep_hat <- if_fisher_psi(muhat, y, p_hat)
  if ("szr" %in% types) out$if_fisher_s_rep_hat <- if_fisher_psi(muhat, y, s_hat)
  if ("co" %in% types) out$if_theta_hat <- if_theta(muhat, y)
  
  return(out)
}

# Function to compute estimate based on plug in and computing influence function
add_instant_estimators <- function(out, muhat, y, types, p_hat, s_hat) {
  #browser()
  if (any(substr(types, 1, 1) == "p")) out$pi_p <- p_hat
  if (any(substr(types, 1, 1) == "s")) out$pi_s <- s_hat
  if (any(c("ct", "cl", "cq", "cz") %in% types)) out$theta_hat <- theta_est(muhat, y)
  
  if ("pf" %in% types) out$pf <- out$pi_p[1] + mean(out$if_full_hat)
  if ("pm" %in% types) out$pm <- out$pi_p[1] + mean(out$if_mixed_hat)
  if ("sf" %in% types) out$sf <- out$pi_s[1] + mean(out$if_full_hat)
  if ("sm" %in% types) out$sm <- out$pi_s[1] + mean(out$if_mixed_hat)
  if ("ct" %in% types) out$ct <- ifelse(out$theta_hat[1] < 1, sqrt(1 - out$theta_hat[1]), 0)
  if ("pl" %in% types) out$pl <- logit_pe(muhat, y, out$pi_p[1], out$if_logit_hat)
  if ("sl" %in% types) out$sl <- logit_pe(muhat, y, out$pi_s[1], out$if_logit_hat)
  if ("cl" %in% types) out$cl <- ifelse(out$theta_hat[1] < 1, logit_pe(muhat, y, sqrt(1 - out$theta_hat[1]), out$if_logit_theta_hat), 0)
  if ("pq" %in% types) out$pq <- logit_sq_pe(muhat, y, out$pi_p[1], out$if_logit_sq_hat)
  if ("sq" %in% types) out$sq <- logit_sq_pe(muhat, y, out$pi_s[1], out$if_logit_sq_hat)
  if ("cq" %in% types) out$cq <- ifelse(out$theta_hat[1] < 1, logit_sq_pe(muhat, y, sqrt(1 - out$theta_hat[1]), out$if_logit_sq_theta_hat), 0)
  if ("pz" %in% types) out$pz <- fisher_pe(muhat, y, out$pi_p[1], out$if_fisher_hat)
  if ("sz" %in% types) out$sz <- fisher_pe(muhat, y, out$pi_s[1], out$if_fisher_hat)
  if ("cz" %in% types) out$cz <- ifelse(out$theta_hat[1] < 1, fisher_pe(muhat, y, sqrt(1 - out$theta_hat[1]), out$if_fisher_theta_hat), 0)
  if ("pfr" %in% types) out$pfr <- out$pi_p[1] + mean(out$if_full_rep_hat)
  if ("pmr" %in% types) out$pmr <- out$pi_p[1] + mean(out$if_mixed_rep_hat)
  if ("plr" %in% types) out$plr <- logit_pe(muhat, y, p_hat, out$if_logit_p_rep_hat)
  if ("slr" %in% types) out$slr <- logit_pe(muhat, y, s_hat, out$if_logit_s_rep_hat)
  if ("pqr" %in% types) out$pqr <- logit_sq_pe(muhat, y, p_hat, out$if_logit_sq_p_rep_hat)
  if ("sqr" %in% types) out$sqr <- logit_sq_pe(muhat, y, s_hat, out$if_logit_sq_s_rep_hat)
  if ("pzr" %in% types) out$pzr <- fisher_pe(muhat, y, p_hat, out$if_fisher_p_rep_hat)
  if ("szr" %in% types) out$szr <- fisher_pe(muhat, y, p_hat, out$if_fisher_s_rep_hat)
  if ("co" %in% types) out$co <- out$ct[1] + mean(out$if_theta_hat)
  
  return(out)
}

# Function to add average of estimates across folds for instant estimators
add_instant_final_estimates <- function(out, res, types) {
  #browser()
  for (type in types) {
    if (type %in% c("pf", "pm", "sf", "sm", "ct", "pl", "sl", "cl", "pq", "sq", "cq", "pz", "sz", "cz", "pfr", "pmr", "plr", "slr", "pqr", "sqr", "pzr", "szr", "co")) {
      out[[paste0(type, "_inst")]] <- mean(res[[type]], na.rm = TRUE)
    }
  }
  return(out)
}

# Function to add Hold estimates after cross fitting
# need to work more on this function to use the right transformations xx
add_hold_final_estimates <- function(out, res, types) {
  #browser()
  # get plug-in hold estimates
  if (any(substr(types, 1, 1) == "p")) {
    pi_p <- cor(res$muhat, res$y)
  }
  if (any(substr(types, 1, 1) == "s")) {
    pi_s <- unname(corr_simp(res$muhat, res$y))
  }
  if (any(c("ct", "cl", "cq", "cz") %in% types)) {
    theta_hat <- theta_est(res$muhat, res$y)
  }
  
  # compute the IF-based estimator using the results across all folds
  for (type in types) {
    if (type %in% c("pf", "pm", "sf", "sm", "pl", "sl", "pq", "sq", "pz", "sz", "pfr", "pmr", "plr", "slr", "pqr", "sqr", "pzr", "szr")) {
      # get the right label for the corresponding influence function
      tlab = get_tlab(type)
      plug_in = ifelse(substr(type,1,1) == "p", pi_p, pi_s)
      # calculate one-step estimate
      os = plug_in + mean(res[[paste0("if_", tlab, "_hat")]])
      # if transformation is needed
      t = substr(type, 2,2)
      if(t == "l") os =  logit_pe(res$muhat, res$y, plug_in, res[[paste0("if_", tlab, "_hat")]]) # expit(logit(plug_in) + mean(res[[paste0("if_", tlab, "_hat")]]))
      if(t == "q") os = logit_sq_pe(res$muhat, res$y, plug_in, res[[paste0("if_", tlab, "_hat")]]) #sqrt(expit(logit(plug_in^2) + mean(res[[paste0("if_", tlab, "_hat")]])))
      if(t == "z") os =  fisher_pe(res$muhat, res$y, plug_in, res[[paste0("if_", tlab, "_hat")]]) #tanh(atanh(plug_in) + mean(res[[paste0("if_", tlab, "_hat")]]))
      out[[paste0(type, "_hold")]] <- os
    } 
    ## EE based estimators
    if(type == "ct") out$ct_hold <- ifelse(theta_hat < 1, sqrt(1-theta_hat), 0)
    if(type=="cl") out$cl_hold <- ifelse(theta_hat < 1, logit_pe(res$muhat, res$y, sqrt(1-theta_hat), res$if_logit_theta_hat), 0)
    if(type=="cq") out$cq_hold <- ifelse(theta_hat < 1, logit_sq_pe(res$muhat, res$y, sqrt(1 - theta_hat), res$if_logit_sq_theta_hat), 0)
    if(type=="cz") out$cz_hold <- ifelse(theta_hat < 1, fisher_pe(res$muhat, res$y, sqrt(1 - theta_hat), res$if_fisher_theta_hat), 0)
    if(type=="co") out$co_hold = out$ct[1] + mean(res$if_theta_hat) # check this one xx
  }
  return(out)
}

# Add final estimates to output
compute_final_estimates <- function(res, types, cf_method, cl, n, ci, nboot) {
  #browser()
  out <- data.frame(p_hold = cor.test(res$muhat, res$y, conf.level = cl)$estimate)
  out$lCI_p_hold <- cor.test(res$muhat, res$y, conf.level = cl)$conf.int[1]
  out$uCI_p_hold <- cor.test(res$muhat, res$y, conf.level = cl)$conf.int[2]
  
  # instant pearson plug-in
  p_avg <- mean(res$pi_p, na.rm = TRUE)
  out$p_inst <- p_avg
  out$lCI_p_inst <- tanh(atanh(p_avg) - qnorm(cl+(1-cl)/2) * (1 / sqrt(n - 3)))
  out$uCI_p_inst <- tanh(atanh(p_avg) + qnorm(cl+(1-cl)/2) * (1 / sqrt(n - 3)))
  
  # simplified plug-ins for reference (don't have a varaiance method)
  out$s_hold = corr_simp(res$muhat, res$y)
  out$s_inst = ifelse("instant" %in% cf_method, mean(res$pi_s, na.rm = T), out$s_hold)
                         
  # add point estimates to output                       
  if ("instant" %in% cf_method | "inst" %in% cf_method) {
    out <- add_instant_final_estimates(out, res, types)
  }
  
  if ("hold" %in% cf_method) {
    out <- add_hold_final_estimates(out, res, types)
  }
                                                
  if (ci) {
    out <- add_confidence_intervals(out, res, types, cf_method, cl, n, nboot)
  }
  
  return(out)
}

# Function to add the confidence intervals to the results
add_confidence_intervals <- function(out, res, types, cf_method, cl, n, nboot) {
  #browser()
  for (type in types) {
    # instant/hold CIs (both use hold variance - not a clear theory for what the instant variance would be)
    if (type %in% c("pf", "pm", "sf", "sm", "ct", "pl", "sl", "cl", "pq", "sq", "cq", "pz", "sz", "cz", "pfr", "pmr", "plr", "slr", "pqr", "sqr", "pzr", "szr", "co")) {
      # get empirical variance of the influence functions (for instant/hold CIs)
      tlab = get_tlab(type)
      v <- var(res[[paste0("if_", tlab, "_hat")]])
      
      # get relevant transformation
      t = substr(type, 2,2)
      
      for(cf in cf_method){
        suff = substr(cf, 1,4)
        # add to the point estimate
        if(!(t %in% c("l","q","z"))){
          out[[paste0("lCI_", type, "_", suff)]] <- out[[paste0(sub("_.", "",type), "_", suff)]] - qnorm(cl+(1-cl)/2) * sqrt(v / n)
          out[[paste0("uCI_", type, "_", suff)]] <- out[[paste0(sub("_.", "",type), "_", suff)]] + qnorm(cl+(1-cl)/2) * sqrt(v / n)
        } 
        if(t == "l"){
          out[[paste0("lCI_", type, "_", suff)]] <- expit(logit(out[[paste0(sub("_.", "",type), "_", suff)]]) - qnorm(cl+(1-cl)/2) * sqrt(v / n))
          out[[paste0("uCI_", type, "_", suff)]] <- expit(logit(out[[paste0(sub("_.", "",type), "_", suff)]]) + qnorm(cl+(1-cl)/2) * sqrt(v / n))
        }
        if(t=="q"){
          out[[paste0("lCI_", type, "_", suff)]] <- sqrt(expit(logit(out[[paste0(sub("_.", "",type), "_", suff)]]^2) - qnorm(cl+(1-cl)/2) * sqrt(v / n)))
          out[[paste0("uCI_", type, "_", suff)]] <- sqrt(expit(logit(out[[paste0(sub("_.", "",type), "_", suff)]]^2) + qnorm(cl+(1-cl)/2) * sqrt(v / n)))
        }
        if(t == "z"){
          out[[paste0("lCI_", type, "_", suff)]] <- tanh(atanh(out[[paste0(sub("_.", "",type), "_", suff)]]) - qnorm(cl+(1-cl)/2) * sqrt(v / n))
          out[[paste0("uCI_", type, "_", suff)]] <- tanh(atanh(out[[paste0(sub("_.", "",type), "_", suff)]]) + qnorm(cl+(1-cl)/2) * sqrt(v / n))
        }
      }
      
      
    }
    
    # Multiplier/NP bootstrap CIs
    # Only applicable for one-step estimators
    if (type %in% c("pf", "pm", "sf", "sm", "pl", "sl", "pq", "sq", "pz", "sz", "pfr", "pmr", "plr", "slr", "pqr", "sqr", "pzr", "szr")) {
      # get the right label for the corresponding influence function
      tlab = get_tlab(type)
      
      est_se = sd(res[[paste0("if_", tlab, "_hat")]])*(1/sqrt(n))
      
      # multiplier bootstrap
      ## xx changed 3/4/25 to use both inst and hold as point estimate instead of just inst (need to recheck the literature)
      # xx changed 3/17/25 to work on correcting the quantiles/procedure, adding percentile-t bootstrap
      # boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_inst")]], transform = t)
      # out[[paste0("lCI_", type, "_MBinst")]] <- quantile(boot_obj$t, (1 - cl) / 2)
      # out[[paste0("uCI_", type, "_MBinst")]] <- quantile(boot_obj$t, cl + (1 - cl) / 2)
      
      
      boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_hold")]], transform = t)
      # out[[paste0("lCI_", type, "_MBhold")]] <- quantile(boot_obj$t, (1 - cl) / 2)
      # out[[paste0("uCI_", type, "_MBhold")]] <- quantile(boot_obj$t, cl + (1 - cl) / 2)
      # out[[paste0("lCI_", type, "_MBhold")]] <- out[[paste0(type, "_hold")]] - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2)  # centered at one-step
      # out[[paste0("uCI_", type, "_MBhold")]] <- out[[paste0(type, "_hold")]] - est_se*quantile(boot_obj$t, (1 - cl) / 2)
      out[[paste0("lCI_", type, "_MBhold")]] <- 2*out[[paste0(type, "_hold")]] - quantile(boot_obj$t, cl + (1 - cl) / 2)  # centered at one-step
      out[[paste0("uCI_", type, "_MBhold")]] <- 2*out[[paste0(type, "_hold")]] - quantile(boot_obj$t, (1 - cl) / 2)
      if("instant" %in% cf_method){
        boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_inst")]], transform = t)
        out[[paste0("lCI_", type, "_MBinst")]] <- 2*out[[paste0(type, "_inst")]] - quantile(boot_obj$t, cl + (1 - cl) / 2)  # centered at one-step
        out[[paste0("uCI_", type, "_MBinst")]] <- 2*out[[paste0(type, "_inst")]] - quantile(boot_obj$t, (1 - cl) / 2)
      }
      
      # nonparametric bootstrap
      # xx changed 3/17/25 - do NOT center the IFs before bootstrapping, want the mean to be the mean of the IF, not 0
      #boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]] - mean(res[[paste0("if_", tlab, "_hat")]]), statistic = bootstrap_stat_np, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_inst")]],  transform = t)
      # probably still need to change the quantiles
      if("instant" %in% cf_method){
      boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat_np, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_inst")]],  transform = t)
      # out[[paste0("lCI_", type, "_NPBinst")]] <- quantile(boot_obj$t, (1 - cl) / 2)
      # out[[paste0("uCI_", type, "_NPBinst")]] <- quantile(boot_obj$t, cl + (1 - cl) / 2)
      # one-step based quantiles
      out[[paste0("lCI_", type, "_NPBinst")]] <- 2*out[[paste0(type, "_inst")]] - quantile(boot_obj$t, cl + (1 - cl) / 2)
      out[[paste0("uCI_", type, "_NPBinst")]] <- 2*out[[paste0(type, "_inst")]] - quantile(boot_obj$t, (1 - cl) / 2)}
      
      boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat_np, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_hold")]],  transform = t)
      #boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]] - mean(res[[paste0("if_", tlab, "_hat")]]), statistic = bootstrap_stat_np, R = nboot, plug_in = out[[paste0(substr(type,1,1), "_hold")]],  transform = t)
      # out[[paste0("lCI_", type, "_NPBhold")]] <- quantile(boot_obj$t, (1 - cl) / 2)
      # out[[paste0("uCI_", type, "_NPBhold")]] <- quantile(boot_obj$t, cl + (1 - cl) / 2)
      out[[paste0("lCI_", type, "_NPBhold")]] <- 2*out[[paste0(type, "_hold")]] - quantile(boot_obj$t, cl + (1 - cl) / 2)
      out[[paste0("uCI_", type, "_NPBhold")]] <- 2*out[[paste0(type, "_hold")]] - quantile(boot_obj$t, (1 - cl) / 2)
      
      # percentile-t bootstrap (based on Hall textbook)
      boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat_pt, R = nboot, transform = t)
      if(t %in% c("f", "m")){
        out[[paste0("lCI_", type, "_PTBhold")]] <- out[[paste0(type, "_hold")]] - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2)  # centered at one-step
        out[[paste0("uCI_", type, "_PTBhold")]] <- out[[paste0(type, "_hold")]] - est_se*quantile(boot_obj$t, (1 - cl) / 2)
      }
      if(t == "l"){
        out[[paste0("lCI_", type, "_PTBhold")]] <- expit(logit(out[[paste0(type, "_hold")]]) - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2))  # centered at one-step
        out[[paste0("uCI_", type, "_PTBhold")]] <- expit(logit(out[[paste0(type, "_hold")]]) - est_se*quantile(boot_obj$t, (1 - cl) / 2))
      }
      if(t == "q"){
        out[[paste0("lCI_", type, "_PTBhold")]] <- sqrt(expit(logit(out[[paste0(type, "_hold")]]^2) - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2)))  # centered at one-step
        out[[paste0("uCI_", type, "_PTBhold")]] <- sqrt(expit(logit(out[[paste0(type, "_hold")]]^2) - est_se*quantile(boot_obj$t, (1 - cl) / 2)))
      }
      if(t == "z"){
        out[[paste0("lCI_", type, "_PTBhold")]] <- tanh(atanh(out[[paste0(type, "_hold")]]) - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2))  # centered at one-step
        out[[paste0("uCI_", type, "_PTBhold")]] <- tanh(atanh(out[[paste0(type, "_hold")]]) - est_se*quantile(boot_obj$t, (1 - cl) / 2))
      }
      
      if("instant" %in% cf_method){
        #browser()
        boot_obj <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat_pt, R = nboot, transform = t)
        # # double bootstrap
        # if(type %in% c("pq", "pqr")){
        #   boot_obj_double <- boot(data = res[[paste0("if_", tlab, "_hat")]], statistic = bootstrap_stat_pt_double, R = nboot, transform = t)
        # }
        if(t %in% c("f", "m")){
          out[[paste0("lCI_", type, "_PTBinst")]] <- out[[paste0(type, "_inst")]] - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2)  # centered at one-step
          out[[paste0("uCI_", type, "_PTBinst")]] <- out[[paste0(type, "_inst")]] - est_se*quantile(boot_obj$t, (1 - cl) / 2)
          }
        if(t == "l"){
          out[[paste0("lCI_", type, "_PTBinst")]] <- expit(logit(out[[paste0(type, "_inst")]]) - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2))  # centered at one-step
          out[[paste0("uCI_", type, "_PTBinst")]] <- expit(logit(out[[paste0(type, "_inst")]]) - est_se*quantile(boot_obj$t, (1 - cl) / 2))
         }
        if(t == "q"){
          out[[paste0("lCI_", type, "_PTBinst")]] <- sqrt(expit(logit(out[[paste0(type, "_inst")]]^2) - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2)))  # centered at one-step
          out[[paste0("uCI_", type, "_PTBinst")]] <- sqrt(expit(logit(out[[paste0(type, "_inst")]]^2) - est_se*quantile(boot_obj$t, (1 - cl) / 2)))
        #   if(type == "pq"){
        #   out[[paste0("lCI_", type, "_DPTBinst")]] <- sqrt(expit(logit(out[[paste0(type, "_inst")]]^2) - est_se*quantile(boot_obj_double$t, cl + (1 - cl) / 2)))  # centered at one-step
        #   out[[paste0("uCI_", type, "_DPTBinst")]] <- sqrt(expit(logit(out[[paste0(type, "_inst")]]^2) - est_se*quantile(boot_obj_double$t, (1 - cl) / 2)))
        # }
        }
        if(t == "z"){
          out[[paste0("lCI_", type, "_PTBinst")]] <- tanh(atanh(out[[paste0(type, "_inst")]]) - est_se*quantile(boot_obj$t, cl + (1 - cl) / 2))  # centered at one-step
          out[[paste0("uCI_", type, "_PTBinst")]] <- tanh(atanh(out[[paste0(type, "_inst")]]) - est_se*quantile(boot_obj$t, (1 - cl) / 2))
          }
      }
      
    }
  }

  return(out)
}

get_tlab = function(type){
  #browser()
  # get the right label for the corresponding influence function
  t = substr(type, 2,2)
  tlab = ifelse(t %in% c("f", "t"), "full", 
                ifelse(t == "m", "mixed", 
                       ifelse(t=="l", "logit", 
                              ifelse(t=="q", "logit_sq", 
                                     ifelse(t=="z", "fisher", 
                                            ifelse(t=="o", "theta", NULL))))))
  # for EE
  tlab = ifelse(type == "cl", "logit_theta",
                ifelse(type == "cq", "logit_sq_theta",
                       ifelse(type=="cz", "fisher_theta", tlab)))
  
  # for reparameterized version
  if(substr(type,nchar(type),nchar(type)) == "r"){
    if(t %in% c("f", "m")) tlab <- paste0(tlab, "_rep")
    else tlab <- paste0(tlab, "_", substr(type, 1, 1), "_rep")
  }
  return(tlab)
}

# model: object with corresponding model class (for determining proper method)
# dt: data.frame, data
# K: number of folds for cross-fitting
# yname: character, name of dependent variable
# model_args: list of arguments (with names indicating argument names)
# types: which estimators to compute (two or three characters)
# cf_method: cross-fitting method, either/both "instant" to compute estimates within each fold and then average, or "hold" to compute predicted values and then compute estimate with all observations
# ci: logical, default TRUE. Whether to calculate variance and confidence interval
# clevel: numeric, default 0.95. Confidence level for confidence interval
# nboot: numeric, number of bootstrap replicates for the multiplier bootstrap
cf_corr_if = function(model, dt, K, yname, model_args = list(), types,
                      cf_method = "hold", ci = TRUE,
                      cl = 0.95, nboot = 1000, ...){
  UseMethod("cf_corr_if")}

# form: optional formula specification. One of `formula` or `xname`/`yname` must be specified
# xname: character name of matrix or variable for X
# model_func: function or non-empty character string naming the ML function
cf_corr_if.default = function(model, dt, K, form = NULL, yname = NULL, xname = NULL,
                              model_args = list(), types,
                              cf_method = "hold", ci = TRUE,
                              cl = 0.95, model_func = lm, nboot = 1000, XtXinv = NULL, ...){
  #browser()
  n = nrow(dt)
  
  # randomly split data into K folds
  split_ind = suppressWarnings(split(1:n, 1:K))
  
  # in each fold:
  res = lapply(1:K, function(i){
    #browser()
    
    # define reference (i.e., training) and inference (i.e., test) groups
    ref_dt = dt[-split_ind[[i]],]
    inf_dt = dt[split_ind[[i]],]
    
    y = as.vector(inf_dt[,yname])
    out = data.frame(y = y) # initialize results
    
    # set model arguments
    # use either formula or set x and y arguments
    if(is.null(form)){
      model_args$y = ref_dt[,yname]
      model_args$x = as.matrix(ref_dt[,xname])
    } else{
      model_args$data = ref_dt
      model_args$formula = form}
    
    # for rho_2, not refitting a model, just using this function to grab estimators of correlation
    if(identical(model_func, rho_2_func)){
      model_args$preds = inf_dt$rho_2_preds
    }
    
    # build the reference model for this split
    ref_model = do.call(model_func, model_args)
    
    # predict on this split
    if(is.null(xname)){
      muhat = get_pred(ref_model, inf_dt, ...)
    } else{
      muhat = get_pred(ref_model, inf_dt[,xname], ...)
    }
    # get predicted values for inference group
    muhat = out$muhat = as.vector(muhat)
    
    # plug in estimators
    p_hat = cor(muhat, y) # pearson
    s_hat = unname(corr_simp(muhat, y)) # simplified
    
    # get needed influence function estimates for inference group based on types
    out <- add_influence_functions(out, muhat, y, types, p_hat, s_hat)
    
    if ("instant" %in% cf_method | "inst" %in% cf_method) {
      out <- add_instant_estimators(out, muhat, y, types, p_hat, s_hat)
    }
    
    out$pi_p = p_hat
    out$pi_s = s_hat
    
    # for rho_2
    # add fold index
    out$fold = i
    out$xi = inf_dt$fakeYMean
    # add lambda information
    if(identical(model_func, glmnet::cv.glmnet)) out$lambda =ref_model$lambda.min
    return(out)
  })
  # combine results from all folds
  res = do.call(rbind, res)

  out <- compute_final_estimates(res, types, cf_method, cl, n, ci, nboot)
  
  
  
  # add lambda values for ridge regression
  if(identical(model_func, glmnet::cv.glmnet)) out = cbind(out, unique(res$lambda))
  
  
  
  # return(out)
  return(list(estimates = out, obs_data = res))
}



# multiplier bootstrap function
# takes original data (vector of influence functions), and inds (resampled indices), plug-in estimate, optional transform character
# EDITED 3/17/25: changed to follow article more closely
# still trying to figure out
bootstrap_stat = function(dt, inds, plug_in, transform = NULL){
  # Draw multipliers from N(0, 1)
  multipliers = rnorm(length(inds)) 
  
  # center at the one-step (this will only be valid for pf/pm)
  
  # Deal with transformations (logit, sq logit, Fisher)
  if(transform %in% c("l", "q", "z")){
    # Truncate plug-in as needed
    if(plug_in < 0 & transform == "l") plug_in = 0
    if(plug_in < -1 & transform %in% c("q", "z")) plug_in = -1
    if(plug_in > 1) plug_in = 1

    # calculate one-step estimate
    if(transform == "l") os = expit(logit(plug_in) + mean(dt) + mean(multipliers*(dt[inds] - mean(dt))))
    if(transform == "q") os = sqrt(expit(logit(plug_in^2) + mean(dt) + mean(multipliers*(dt[inds] - mean(dt)))))
    if(transform == "z") os = tanh(atanh(plug_in) + mean(dt) + mean(multipliers*(dt[inds] - mean(dt))))

  } else{
    # calculate one-step estimate
    os = plug_in + mean(dt) +  mean(multipliers*(dt[inds] - mean(dt)))
  }
  return(os)
}



bootstrap_stat_np = function(dt, inds, plug_in, transform = NULL){
  # Deal with transformations (logit, sq logit, Fisher)
  if(transform %in% c("l", "q", "z")){
    # Truncate plug-in as needed
    if(plug_in < 0 & transform == "l") plug_in = 0
    if(plug_in < -1 & transform %in% c("q", "z")) plug_in = -1
    if(plug_in > 1) plug_in = 1
    
    # calculate one-step estimate
    if(transform == "l") os = expit(logit(plug_in) + mean(dt[inds]))
    if(transform == "q") os = sqrt(expit(logit(plug_in^2) + mean(dt[inds])))
    if(transform == "z") os = tanh(atanh(plug_in) + mean(dt[inds]))
    
  } else{
    # calculate one-step estimate
    os = plug_in + mean(dt[inds])
  }
  
  return(os)
}

# percentile-t bootstrap
bootstrap_stat_pt = function(dt, inds, transform = NULL){
  
  # transform to t-stat (plug-ins will cancel)
  sqrt(length(dt))*(mean(dt[inds]) - mean(dt))/sd(dt[inds])
  
}

# interested in seeing if double bootstrap works better
bootstrap_stat_pt_double = function(dt, inds, transform = NULL){
  # initial bootstrap sample
  # transform to t-stat (plug-ins will cancel)
  t_b = sqrt(length(dt))*(mean(dt[inds]) - mean(dt))/sd(dt[inds])
  
  # second bootstrap - sample from the bootstrap sample (50 replicates)
  boot_inner = boot::boot(data = dt[inds], statistic = bootstrap_stat_pt, R = 50, transform = transform)
  
  # bias correct t_b (taken from Chatgpt, need to find reference)
  t_b = 2*t_b - mean(boot_inner$t)
  t_b
}

