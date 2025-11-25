### Functions for compiling simulation results ###


# function to load in the results for a specific path
get_results = function(path, ns, nsim, r1name, ns_10){
  if(missing(ns_10)){
    ns_10 = ns
  }
  # load in all results 
  t_res_25 = list()
  for(n in ns){
    for (sim in 1:nsim){
      if(file.exists((paste0(path, "rho_25/n",n,"/simResults_",sim,".rds")))){
        t_res_25[[paste0("n",n)]][[sim]] = readRDS(paste0(path, "rho_25/n",n,"/simResults_",sim,".rds"))
      }
      else{
        t_res_25[[paste0("n",n)]][[sim]] = NULL
      }
    }}
  
  t_res_50 = list()
  for(n in ns){
    for (sim in 1:nsim){
      if(file.exists((paste0(path, "rho_50/n",n,"/simResults_",sim,".rds")))){
        t_res_50[[paste0("n",n)]][[sim]] = readRDS(paste0(path, "rho_50/n",n,"/simResults_",sim,".rds"))
      }
      else{
        t_res_50[[paste0("n",n)]][[sim]] = NULL
      }
    }}
  
  
  t_res_75 = list()
  for(n in ns){
    for (sim in 1:nsim){
      if(file.exists((paste0(path, "rho_75/n",n,"/simResults_",sim,".rds")))){
        t_res_75[[paste0("n",n)]][[sim]] = readRDS(paste0(path, "rho_75/n",n,"/simResults_",sim,".rds"))
      }
      else{
        t_res_75[[paste0("n",n)]][[sim]] = NULL
      }
    }}
  
  
  # rho_0 = 0.1 (to check when very low true correlation)
  t_res_10 = list()
  for(n in ns_10){
    for (sim in 1:nsim){
      if(file.exists((paste0(path, "rho_10/n",n,"/simResults_",sim,".rds")))){
        t_res_10[[paste0("n",n)]][[sim]] = readRDS(paste0(path, "rho_10/n",n,"/simResults_",sim,".rds"))
      }
      else{
        t_res_10[[paste0("n",n)]][[sim]] = NULL
      }
    }}
  
  ## parameter values
  # rho_0
  rho_0_25 = readRDS(paste0(path, "rho_25/rho_0.rds"))[2]
  rho_0_50 = readRDS(paste0(path, "rho_50/rho_0.rds"))[2]
  rho_0_75 = readRDS(paste0(path, "rho_75/rho_0.rds"))[2]
  rho_0_10 = readRDS(paste0(path, "rho_10/rho_0.rds"))[2]
  
  # rho_1
  rho_1_25 = calc_rho_1(t_res_25, ns, r1name)
  rho_1_50 = calc_rho_1(t_res_50, ns, r1name)
  rho_1_75 = calc_rho_1(t_res_75, ns, r1name)
  rho_1_10 = calc_rho_1(t_res_10, ns_10, r1name)
  
  # rho_2
  rho_2_25 = calc_rho_2(t_res_25, ns, "rho_2")
  rho_2_50 = calc_rho_2(t_res_50, ns, "rho_2")
  rho_2_75 = calc_rho_2(t_res_75, ns, "rho_2")
  rho_2_10 = calc_rho_2(t_res_10, ns_10, "rho_2")
  
  # estimates
  allResults_25 = get_if_res(t_res_25, ns)
  allResults_rho_2_25 = get_if_res_rho_2(t_res_25, ns)
  
  allResults_50 = get_if_res(t_res_50, ns)
  allResults_rho_2_50 = get_if_res_rho_2(t_res_50, ns)
  
  allResults_75 = get_if_res(t_res_75, ns)
  allResults_rho_2_75 = get_if_res_rho_2(t_res_75, ns)
  
  allResults_10 = get_if_res(t_res_10, ns_10)
  allResults_rho_2_10 = get_if_res_rho_2(t_res_10, ns_10)
  
  # return everything in a list
  return(list(r10 = list(rho_0 = rho_0_10, rho_1 = rho_1_10, rho_2 = rho_2_10, allResults = allResults_10, allResults_rho_2 = allResults_rho_2_10),
              r25 = list(rho_0 = rho_0_25, rho_1 = rho_1_25, rho_2 = rho_2_25, allResults = allResults_25, allResults_rho_2 = allResults_rho_2_25),
              r50 = list(rho_0 = rho_0_50, rho_1 = rho_1_50, rho_2 = rho_2_50, allResults = allResults_50, allResults_rho_2 = allResults_rho_2_50),
              r75 = list(rho_0 = rho_0_75, rho_1 = rho_1_75, rho_2 = rho_2_75, allResults = allResults_75, allResults_rho_2 = allResults_rho_2_75)))
  
}

# Truncate estimators to [0,1]
# last column in simulations is n
truncate_cols = function(df, num_end_cols=1, max_val = 1, min_val = 0){
  df[ , 1:(ncol(df)-num_end_cols)] <- lapply(df[ , 1:(ncol(df)-num_end_cols)], function(col) {
    if (is.numeric(col)) {
      pmax(pmin(col, max_val), min_val)
    } else {
      col
    }
  })
  df
}


## Calculating parameters from simulation results
# rho_0 (MAPA)
calc_rho_0_s = function(res, ns){
  rho_0_s_sim = map_depth(res, 2, "rho_0")
  rho_0_s_sim = map_depth(rho_0_s_sim, 2, "rho_0_s")
  rho_0_s_sim = as.data.frame(do.call(cbind, rho_0_s_sim))
  rho_0_s_sim = sapply(1:ncol(rho_0_s_sim), function(x) unlist(rho_0_s_sim[,x]))
  if("list" %in% class(rho_0_s_sim)){ # list indicates some of the sims didn't run (i.e. ran out of memory - oops)
    rho_0_s_res = unlist(lapply(rho_0_s_sim, mean))
  } else{
    rho_0_s_res = colMeans(rho_0_s_sim)
  }
  names(rho_0_s_res) = ns
  return(rho_0_s_res)
}

# MAPA
calc_rho_0_p = function(res, ns){
  rho_0_p_sim = map_depth(res, 2, "rho_0")
  rho_0_p_sim = map_depth(rho_0_p_sim, 2, "rho_0_p")
  rho_0_p_sim = as.data.frame(do.call(cbind, rho_0_p_sim))
  rho_0_p_sim = sapply(1:ncol(rho_0_p_sim), function(x) unlist(rho_0_p_sim[,x]))
  if("list" %in% class(rho_0_p_sim)){
    rho_0_p_res = unlist(lapply(rho_0_p_sim, mean))
  } else{
    rho_0_p_res = colMeans(rho_0_p_sim)
  }
  names(rho_0_p_res) = ns
  return(rho_0_p_res)
}

# rho_1 (not described in paper)
calc_rho_1 = function(res, ns, col_name){
  rho_1_sim = map_depth(res, 2, col_name)
  rho_1_sim = as.data.frame(do.call(cbind, rho_1_sim))
  rho_1_sim = sapply(1:ncol(rho_1_sim), function(x) unlist(rho_1_sim[,x]))
  if("list" %in% class(rho_1_sim)){
    rho_1 = unlist(lapply(rho_1_sim, mean))
  } else{
    rho_1= colMeans(rho_1_sim)
  }
  names(rho_1) = ns
  return(rho_1)
}

# rho_2 (TPA)
calc_rho_2 = function(res, ns, col_name){
  rho_2_sim = map_depth(res, 2, col_name)
  rho_2_sim = as.data.frame(do.call(cbind, rho_2_sim))
  rho_2_sim = sapply(1:ncol(rho_2_sim), function(x) unlist(rho_2_sim[,x]))
  if("list" %in% class(rho_2_sim)){
    rho_2 = unlist(lapply(rho_2_sim, mean))
  } else{
    rho_2= colMeans(rho_2_sim)
  }
  names(rho_2) = ns
  return(rho_2)
}

# TPA
calc_rho_2_simY = function(res, ns){
  rho_2_sim = map_depth(res, 2, "rho_2_simY")
  rho_2_sim = as.data.frame(do.call(cbind, rho_2_sim))
  rho_2_sim = sapply(1:ncol(rho_2_sim), function(x) unlist(rho_2_sim[,x]))
  if("list" %in% class(rho_2_sim)){
    rho_2 = unlist(lapply(rho_2_sim, mean))
  } else{
    rho_2= colMeans(rho_2_sim)
  }
  names(rho_2) = ns
  return(rho_2)
}

# calculate permutation bounds
calc_perm_bounds = function(res){
  #browser()
  # only for n = 2000
  res = res[["n2000"]]
  
  lower = map_depth(res, 1, "perm_lower")
  lower = as.data.frame(do.call(rbind, lower))
  upper = map_depth(res, 1, "perm_upper")
  upper = as.data.frame(do.call(rbind, upper))
  q95 = map_depth(res, 1, "perm_95")
  q95 = as.data.frame(do.call(rbind, q95))
  res = cbind(lower, upper, q95)

 
  return(res)
}

## Get influence function estimator results
get_if_res = function(res, ns, res_name = "res", delete_lambda = F, nsim=200){
  #browser()
  if_res = list()
  if_res_n = purrr::map_depth(res, 2, res_name)
  for(n in ns){
    if(delete_lambda){
      for(i in 1:nsim){
        # if_res_n[[paste0("n",n)]][[i]] = if_res_n[[paste0("n",n)]][[i]][1,(1:ncol(if_res_n[[paste0("n",n)]][[i]]) - 1)]
        if_res_n[[paste0("n",n)]][[i]] = if_res_n[[paste0("n",n)]][[i]][1,]
      }
    }
    mat = as.data.frame(do.call(rbind, if_res_n[[paste0("n",n)]]))
    mat$n = n
    if_res[[paste0('n', n)]] = mat
  }
  
  allResults = do.call(rbind, if_res)
  return(allResults)
}



get_if_res_rho_2 = function(res, ns){
  if_res = list()
  if_res_n = purrr::map_depth(res, 2, "res_rho_2")
  for(n in ns){
    mat = as.data.frame(do.call(rbind, if_res_n[[paste0("n",n)]]))
    mat$n = n
    if_res[[paste0('n', n)]] = mat
  }
  
  allResults = do.call(rbind, if_res)
  return(allResults)
}

# get coverage probability for plotting
get_cp <- function(data, types, rho) {
  # Check inputs
  stopifnot(is.data.frame(data))
  stopifnot(is.character(types) && length(types) > 0)
  stopifnot(is.numeric(rho))
  if (!"n" %in% names(data)) stop("Column 'n' (sample size) must exist in the data.")
  
  # Extract unique sample sizes
  unique_n <- sort(unique(data$n))
  
  # Handle scalar rho: Expand it to match all unique 'n'
  if (length(rho) == 1) {
    rho_lookup <- setNames(rep(rho, length(unique_n)), unique_n)
  } else if (length(rho) == length(unique_n)) {
    # Handle vector rho: Create lookup table for unique 'n'
    rho_lookup <- setNames(rho, unique_n)
  } else {
    stop("The length of 'rho' must be 1 (scalar) or equal to the number of unique values in 'n'.")
  }
  
  # Initialize an empty list to store results
  cp_list <- list()
  
  # Loop over each type and compute coverage probability
  for (type in types) {
    # Use regex to match columns for lower and upper bounds: "lCI_" and "uCI_"
    lower_cols <- grep(paste0("^lCI_", type, "_"), names(data), value = TRUE)
    upper_cols <- gsub("^lCI_", "uCI_", lower_cols)  # Infer upper columns dynamically
    
    # If no columns found for the type, issue a warning
    if (length(lower_cols) == 0) {
      warning(paste("No confidence interval columns found for type:", type))
      next
    }
    
    # Compute coverage for each confidence interval type under the current estimator type
    for (i in seq_along(lower_cols)) {
      lower_col <- lower_cols[i]
      upper_col <- upper_cols[i]
      
      # Verify that the corresponding upper column exists
      if (!(upper_col %in% names(data))) {
        warning(paste("Upper column", upper_col, "not found. Skipping", lower_col))
        next
      }
      
      # Map rho to rows based on the 'n' column
      rho_values <- rho_lookup[as.character(data$n)]
      
      # Check if rho is within the interval [lower, upper]
      coverage <- (data[[lower_col]] <= rho_values) & (data[[upper_col]] >= rho_values)
      
      # Store coverage probability as a data.frame
      cp_list[[paste0("cp_", type, "_", sub(paste0("^lCI_", type, "_"), "", lower_col))]] <- data.frame(
        n = data$n,
        est = type,
        ci = sub(paste0("^lCI_", type, "_"), "", lower_col),
        coverage = coverage
      )
    }
  }
  
  # Combine all results into a single data.frame
  cp_combined <- do.call(rbind, cp_list)
  
  # Calculate the mean coverage probability grouped by n, estimator, and confidence interval type
  coverage_prob_by_n <- cp_combined %>%
    group_by(n, est, ci) %>%
    summarise(coverage_probability = mean(coverage, na.rm = TRUE), .groups = "drop")
  
  return(coverage_prob_by_n)
}


## Construct tables of bias and coverage
# point estimates
get_est_tab_inst = function(allRes, types){
  inst_est_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$p_inst))))
  for (t in types){
    inst_est_tab = cbind(inst_est_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0(t, "_inst")]))))
  }
  colnames(inst_est_tab) = c("p", all_types)
  return(inst_est_tab)
}

get_est_tab_hold = function(allRes, types){
  hold_est_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$p_hold))))
  for (t in types){
    hold_est_tab = cbind(hold_est_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0(t, "_hold")]))))
  }
  colnames(hold_est_tab) = c("p", all_types)
  return(hold_est_tab)
}

# mean bias
get_bias_tab_inst = function(allRes, types, rho, ns = NULL){
  if(length(rho)==1){
    inst_bias_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$p_inst - rho))))
    for (t in types){
      inst_bias_tab = cbind(inst_bias_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0(t, "_inst")] - rho))))
    }
  } else{
    inst_bias_tab = as.data.frame(matrix(nrow = length(ns)))
    for (t in c("p", types)){
      col_out_inst = c()
      for(n in ns){
        col_out_inst = rbind(col_out_inst, mean(allRes[which(allRes$n == n),paste0(t, "_inst")] - rho[which(ns == n)]))
      }
      inst_bias_tab = cbind(inst_bias_tab, col_out_inst)
    }
    inst_bias_tab = inst_bias_tab[,-1]
  }
  
  colnames(inst_bias_tab) = c("p", types)
  return(inst_bias_tab)
}

get_bias_tab_hold = function(allRes, types, rho, ns = NULL){
  if(length(rho)==1){
    hold_bias_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$p_hold - rho))))
    for (t in types){
      hold_bias_tab = cbind(hold_bias_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0(t, "_hold")] - rho))))
    }
  } else{
    hold_bias_tab = as.data.frame(matrix(nrow = length(ns)))
    for (t in c("p", types)){
      col_out_hold = c()
      for(n in ns){
        col_out_hold = rbind(col_out_hold, mean(allRes[which(allRes$n == n),paste0(t, "_hold")] - rho[which(ns == n)]))
      }
      hold_bias_tab = cbind(hold_bias_tab, col_out_hold)
    }
    hold_bias_tab = hold_bias_tab[,-1]
  }
  
  colnames(hold_bias_tab) = c("p", types)
  return(hold_bias_tab)
}

# median
get_mbias_tab_hold = function(allRes, types, rho, ns = NULL){
  if(length(rho)==1){
    hold_bias_tab = data.frame(p = c(by(allRes, allRes$n, function(x) median(x$p_hold - rho))))
    for (t in types){
      hold_bias_tab = cbind(hold_bias_tab, c(by(allRes, allRes$n, function(x) median(x[,paste0(t, "_hold")] - rho))))
    }
  } else{
    hold_bias_tab = as.data.frame(matrix(nrow = length(ns)))
    for (t in c("p", types)){
      col_out_hold = c()
      for(n in ns){
        col_out_hold = rbind(col_out_hold, median(allRes[which(allRes$n == n),paste0(t, "_hold")] - rho[which(ns == n)]))
      }
      hold_bias_tab = cbind(hold_bias_tab, col_out_hold)
    }
    hold_bias_tab = hold_bias_tab[,-1]
  }
  
  colnames(hold_bias_tab) = c("p", types)
  return(hold_bias_tab)
}

# sd
get_sd_tab_inst = function(allRes, types){
  inst_sd_tab = data.frame(p = c(by(allRes, allRes$n, function(x) sd(x$p_inst))))
  for (t in types){
    inst_sd_tab = cbind(inst_sd_tab, c(by(allRes, allRes$n, function(x) sd(x[,paste0(t, "_inst")]))))
  }
  colnames(inst_sd_tab) = c("p", types)
  return(inst_sd_tab)
}

get_sd_tab_hold = function(allRes, types){
  hold_sd_tab = data.frame(p = c(by(allRes, allRes$n, function(x) sd(x$p_hold))))
  for (t in types){
    hold_sd_tab = cbind(hold_sd_tab, c(by(allRes, allRes$n, function(x) sd(x[,paste0(t, "_hold")]))))
  }
  colnames(hold_sd_tab) = c("p", types)
  return(hold_sd_tab)
}

# coverage probability
get_cp_tab_inst = function(allRes, types, rho, ns = NULL){
  if(length(rho)==1){
    inst_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho & x$uCI_p_inst>rho))))
    for (t in types){
      inst_cp_tab = cbind(inst_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_inst")] < rho & x[,paste0("uCI_",t, "_inst")] > rho))))
    }
  } else{
    inst_cp_tab = as.data.frame(matrix(nrow = length(ns)))
    for (t in c("p", types)){
      col_out_inst = c()
      for(n in ns){
        col_out_inst = rbind(col_out_inst, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_inst")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_inst")]>rho[which(ns==n)]))
      }
      inst_cp_tab = cbind(inst_cp_tab, col_out_inst)
    }
    inst_cp_tab = inst_cp_tab[,-1]
  }
  
  colnames(inst_cp_tab) = c("p", types)
  return(inst_cp_tab)
}

get_cp_tab_hold = function(allRes, types, rho, ns = NULL){
  if(length(rho)==1){
    hold_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho & x$uCI_p_hold>rho))))
    for (t in types){
      hold_cp_tab = cbind(hold_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_hold")] < rho & x[,paste0("uCI_",t, "_hold")] > rho))))
    }
  } else{
    hold_cp_tab = as.data.frame(matrix(nrow = length(ns)))
    for (t in c("p", types)){
      col_out_hold = c()
      for(n in ns){
        col_out_hold = rbind(col_out_hold, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_hold")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_hold")]>rho[which(ns==n)]))
      }
      hold_cp_tab = cbind(hold_cp_tab, col_out_hold)
    }
    hold_cp_tab = hold_cp_tab[,-1]
  }
  
  colnames(hold_cp_tab) = c("p", types)
  return(hold_cp_tab)
}


tab_bias = function(bias_tab, summary_types, cf_type = "Instant"){
  tab = bias_tab[,c("p", summary_types)]
  kable(tab, caption = paste0("Mean Bias, ",  cf_type), digits = 4) %>% kable_styling()
}

tab_est = function(est_tab, summary_types, cf_type = "Instant"){
  tab = est_tab[,c("p", summary_types)]
  kable(tab, caption = paste0("Mean Point Estimates, ",  cf_type), digits = 4) %>% kable_styling()
}

tab_sd = function(sd_tab, summary_types, cf_type = "Instant"){
  tab = sd_tab[,c("p", summary_types)]
  kable(tab, caption = paste0("Standard Deviation, ",  cf_type), digits = 4) %>% kable_styling()
}

tab_cp = function(cp_tab, summary_types, cf_type = "Instant"){
  tab = cp_tab[,c("p", summary_types)]
  kable(tab, caption = paste0("Coverage Probability, ",  cf_type), digits = 4) %>% kable_styling()
}

# MB and cross-fit covariances coverage
get_cp_tab_MB = function(allRes, types, rho, ns = NULL){
  #browser()
  # use instant Pearson as reference
  # changed to hold
  if(length(rho)==1){
    MB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho & x$uCI_p_hold>rho))))
    for (t in types){
      MB_cp_tab = cbind(MB_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_MB")] < rho & x[,paste0("uCI_",t, "_MB")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    MB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho[which(ns==mean(x$n))] & x$uCI_p_hold>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_MB = c()
      for(n in ns){
        col_out_MB = rbind(col_out_MB, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_MB")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_MB")]>rho[which(ns==n)]))
      }
      MB_cp_tab = cbind(MB_cp_tab, col_out_MB)
    }
  }
  
  colnames(MB_cp_tab) = c("p", types)
  return(MB_cp_tab)
}

get_cp_tab_MBinst = function(allRes, types, rho, ns = NULL){
  #browser()
  # use instant Pearson as reference
  # changed to hold
  if(length(rho)==1){
    MBinst_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho & x$uCI_p_inst>rho))))
    for (t in types){
      MBinst_cp_tab = cbind(MBinst_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_MBinst")] < rho & x[,paste0("uCI_",t, "_MBinst")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    MBinst_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho[which(ns==mean(x$n))] & x$uCI_p_inst>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_MBinst = c()
      for(n in ns){
        col_out_MBinst = rbind(col_out_MBinst, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_MBinst")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_MBinst")]>rho[which(ns==n)]))
      }
      MBinst_cp_tab = cbind(MBinst_cp_tab, col_out_MBinst)
    }
  }
  
  colnames(MBinst_cp_tab) = c("p", types)
  return(MBinst_cp_tab)
}

get_cp_tab_MBhold = function(allRes, types, rho, ns = NULL){
  #browser()
  # use holdant Pearson as reference
  # changed to hold
  if(length(rho)==1){
    MBhold_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho & x$uCI_p_hold>rho))))
    for (t in types){
      MBhold_cp_tab = cbind(MBhold_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_MBhold")] < rho & x[,paste0("uCI_",t, "_MBhold")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    MBhold_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho[which(ns==mean(x$n))] & x$uCI_p_hold>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_MBhold = c()
      for(n in ns){
        col_out_MBhold = rbind(col_out_MBhold, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_MBhold")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_MBhold")]>rho[which(ns==n)]))
      }
      MBhold_cp_tab = cbind(MBhold_cp_tab, col_out_MBhold)
    }
  }
  
  colnames(MBhold_cp_tab) = c("p", types)
  return(MBhold_cp_tab)
}

# NPB and cross-fit covariances coverage
get_cp_tab_NPB = function(allRes, types, rho, ns = NULL){
  #browser()
  # use instant Pearson as reference
  if(length(rho)==1){
    NPB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho & x$uCI_p_inst>rho))))
    for (t in types){
      NPB_cp_tab = cbind(NPB_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_NPB")] < rho & x[,paste0("uCI_",t, "_NPB")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    NPB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho[which(ns==mean(x$n))] & x$uCI_p_inst>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_NPB = c()
      for(n in ns){
        col_out_NPB = rbind(col_out_NPB, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_NPB")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_NPB")]>rho[which(ns==n)]))
      }
      NPB_cp_tab = cbind(NPB_cp_tab, col_out_NPB)
    }
  }
  
  colnames(NPB_cp_tab) = c("p", types)
  return(NPB_cp_tab)
}

get_cp_tab_NPBinst = function(allRes, types, rho, ns = NULL){
  #browser()
  # use instant Pearson as reference
  if(length(rho)==1){
    NPBinst_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho & x$uCI_p_inst>rho))))
    for (t in types){
      NPBinst_cp_tab = cbind(NPBinst_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_NPBinst")] < rho & x[,paste0("uCI_",t, "_NPBinst")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    NPBinst_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho[which(ns==mean(x$n))] & x$uCI_p_inst>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_NPBinst = c()
      for(n in ns){
        col_out_NPBinst = rbind(col_out_NPBinst, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_NPBinst")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_NPBinst")]>rho[which(ns==n)]))
      }
      NPBinst_cp_tab = cbind(NPBinst_cp_tab, col_out_NPBinst)
    }
  }
  
  colnames(NPBinst_cp_tab) = c("p", types)
  return(NPBinst_cp_tab)
}

get_cp_tab_NPBhold = function(allRes, types, rho, ns = NULL){
  #browser()
  # use holdant Pearson as reference
  if(length(rho)==1){
    NPBhold_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho & x$uCI_p_hold>rho))))
    for (t in types){
      NPBhold_cp_tab = cbind(NPBhold_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_NPBhold")] < rho & x[,paste0("uCI_",t, "_NPBhold")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    NPBhold_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_hold< rho[which(ns==mean(x$n))] & x$uCI_p_hold>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_NPBhold = c()
      for(n in ns){
        col_out_NPBhold = rbind(col_out_NPBhold, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_NPBhold")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_NPBhold")]>rho[which(ns==n)]))
      }
      NPBhold_cp_tab = cbind(NPBhold_cp_tab, col_out_NPBhold)
    }
  }
  
  colnames(NPBhold_cp_tab) = c("p", types)
  return(NPBhold_cp_tab)
}


get_cp_tab_MB_c = function(allRes, types, rho, ns = NULL){
  #browser()
  # use instant Pearson as reference
  if(length(rho)==1){
    MB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho & x$uCI_p_inst>rho))))
    for (t in types){
      MB_cp_tab = cbind(MB_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_MB_c")] < rho & x[,paste0("uCI_",t, "_MB_c")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    MB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho[which(ns==mean(x$n))] & x$uCI_p_inst>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_MB = c()
      for(n in ns){
        col_out_MB = rbind(col_out_MB, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_MB_c")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_MB_c")]>rho[which(ns==n)]))
      }
      MB_cp_tab = cbind(MB_cp_tab, col_out_MB)
    }
  }
  
  colnames(MB_cp_tab) = c("p", types)
  return(MB_cp_tab)
}

get_cp_tab_NPB_c = function(allRes, types, rho, ns = NULL){
  #browser()
  # use instant Pearson as reference
  if(length(rho)==1){
    NPB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho & x$uCI_p_inst>rho))))
    for (t in types){
      NPB_cp_tab = cbind(NPB_cp_tab, c(by(allRes, allRes$n, function(x) mean(x[,paste0("lCI_",t, "_NPB_c")] < rho & x[,paste0("uCI_",t, "_NPB_c")] > rho))))
    }
  } else{
    # 2/13/25: fixed bug in getting the right rho for p
    NPB_cp_tab = data.frame(p = c(by(allRes, allRes$n, function(x) mean(x$lCI_p_inst< rho[which(ns==mean(x$n))] & x$uCI_p_inst>rho[which(ns==mean(x$n))]))))
    for (t in types){
      col_out_NPB = c()
      for(n in ns){
        col_out_NPB = rbind(col_out_NPB, mean(allRes[which(allRes$n == n),paste0("lCI_",t, "_NPB_c")]<rho[which(ns==n)] &  allRes[which(allRes$n == n),paste0("uCI_",t, "_NPB_c")]>rho[which(ns==n)]))
      }
      NPB_cp_tab = cbind(NPB_cp_tab, col_out_NPB)
    }
  }
  
  colnames(NPB_cp_tab) = c("p", types)
  return(NPB_cp_tab)
}


