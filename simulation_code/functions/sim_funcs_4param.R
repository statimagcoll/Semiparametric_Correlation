## Functions to run simulations when xi =/= mu ###


# Main simulation function
# dt: Eval dataset
# model_func: ML model function to fit
# mu_name: character, name of true muX column in the Eval data (i.e., predictions generated from the build model)
# ls_dist: logspline error distribution for generating Y
# res_mean: numeric, mean of the logspline distribution
# res_var: numeric, variance of the logspline distribution
# rho_2_test_dat: subset of the Eval data used for resampling just for rho_2 estimates
# muhat_tr_mod: model that was fit on the subset of Eval data used to train the model for rho_2 estimates
# ns: vector of sample sizes
# nsim: number of simulation replicates
# K: number of folds for cross-fitting
# form: optional formula for ML model (either form or yname/xname should be specified)
# yname: character, name of generated Y
# xname: character, name of features in the data (i.e., "img" for the matrix of imaging data)
# model_args: list, other arguments to pass to the ML function
# pred_args: list, other arguments to pass to the get_pred() function that gets the predictions from the ML model
# types: character, which estimator types to compute (e.g., "pf", "sqr")
# cf_method: "inst", "hold", or both, which method to use for cross-fitting
# ci: logical, whether to compute confidence intervals
# cl: numeric, confidence level
# ncores: number of cores to use for pbmclapply
# save_path: character, folder path to save results in (each simulation replicate is saved separately)
sim_IF_4param = function(dt, model_func, xi_name, mu_name, ls_dist, res_mean, res_var, rho_2_test_dat, muhat_tr_mod,
                  ns, nsim, K = 5, form = NULL, yname = NULL, xname = NULL,
                  model_args = list(), pred_args = NULL, types, cf_method = "instant", ci = TRUE,
                  cl = 0.95, ncores = 8, save_path, XtXinv = NULL, permute = FALSE, rho_true, ...){
  
  for(n in ns){
    simResults = pbmclapply(1:nsim, function(sim){
      #browser()
      
      # simulation data
      seval = dt[sample(nrow(dt), n, replace = TRUE),]
      # generate fake age
      #seval$fakeY = seval[,xi_name] + rlogspline(nrow(seval), ls_dist) - res_mean
      seval$fakeY = seval[,xi_name] + runif_target(nrow(seval), rho_true, var(dt$fakeYMean))
      
      ## rho_true:
      # another way of getting the parameter
      # sqrt(var(xi(X))/var(Y))
      cor_true = sqrt(var(seval[,xi_name])/var(seval$fakeY)) # simplified
      cor_true_p = cor(seval[,xi_name], seval$fakeY) # pearson
      
      # rho_0 calculated outside of sim
      # check within sim
      cor_0_p = cor(seval[,mu_name], seval$fakeY) # within 2 SEs, seem confident that is accurate 3/19/25
      
      # rho_1
      split_ind_r1 = suppressWarnings(split(1:n, 1:K))
      
      # in each fold:
      r1 = sapply(1:K, function(i){
        #browser()

        # define reference (i.e., training)
        ref_dt = seval[-split_ind_r1[[i]],]
        # set model arguments
        # use either formula and data or set x and y arguments (depending on ML algo function arguments)
        if(is.null(form)){
          model_args$y = ref_dt$fakeY
          model_args$x = as.matrix(ref_dt[,xname])
        } else{
          model_args$data = ref_dt
          model_args$formula = form}

        sim_mod = do.call(model_func, model_args)
        # predict on full eval data for rho_1
        rho_1_preds = get_pred(sim_mod,as.matrix(dt$img), pred_args)

        cor_1 = (mean(rho_1_preds*dt[,xi_name]) - mean(rho_1_preds)*mean(dt[,xi_name]))/
          sqrt(var(rho_1_preds)*(var(dt[,xi_name])+res_var))
        return(cor_1)
      })

      #for one split
      cor_1_onesplit = r1[1]

      # for avg
      cor_1_avg = mean(r1)

      
      
      ## rho_2:
      # predict on sample of Eval_te using given model (data not used to fit model)
      seval_rho_2 = rho_2_test_dat[sample(nrow(rho_2_test_dat), n, replace = TRUE),]
      seval_rho_2$fakeY = seval_rho_2[,xi_name] + runif_target(nrow(seval_rho_2), rho_true, var(dt$fakeYMean))
      seval_rho_2$rho_2_preds = get_pred(muhat_tr_mod, as.matrix(seval_rho_2$img))
      # fixed 10/14/24 - had slightly wrong formula here
      cor_2 = (mean(seval_rho_2$rho_2_preds*seval_rho_2[,xi_name]) - mean(seval_rho_2$rho_2_preds)*mean(seval_rho_2[,xi_name]))/
        sqrt(var(seval_rho_2$rho_2_preds)*(var(seval_rho_2[,xi_name])+v_uniform(rho_true, var(dt$fakeYMean))))
      cor_2_simY = cor(seval_rho_2$rho_2_preds, seval_rho_2$fakeY)
      # need to get estimators without refitting any models
      res_rho_2 = cf_corr_if(muhat_tr_mod, dt = seval_rho_2, K = K, model_func = rho_2_func,
                             yname = "fakeY", xname = "img", types = types, cf_method = c("hold"), XtXinv = XtXinv)

      # cross-fit one-step estimates
      res = cf_corr_if(muhat_tr_mod, dt = seval, K = 5, model_func = model_func,
                       model_args = model_args, form = form, xname = xname,
                       yname = "fakeY", types = types, cf_method = cf_method,
                       pred_args = list(), ci = ci, cl=cl, nboot = 1000, XtXinv = XtXinv)

      # for figure 1, want a permutation test example
      # want the the upper and lower 95% bound of the null distribution - then take average across sims
      if(permute & n == 2000){
        permuted_cor = numeric(1000)
        for (p in 1:1000){
          # shuffle Y values
          seval$fakeY_permute = sample(seval$fakeY)

          # perform the cross-fitting
          perm_res = sapply(1:K, function(i){
            #browser()

            # define reference (i.e., training)
            ref_dt = seval[-split_ind_r1[[i]],]
            # set model arguments

            model_args$y = ref_dt$fakeY_permute
            model_args$x = as.matrix(ref_dt[,xname])

            sim_mod = do.call(model_func, model_args)
            # predict on held out data
            perm_preds = get_pred(sim_mod,as.matrix(seval[split_ind_r1[[i]],]$img), pred_args)
            # return Pearson correlation
            return(cor(perm_preds, seval[split_ind_r1[[i]],]$fakeY_permute))
          })
          # take average of the 5 folds
          permuted_cor[p] = mean(perm_res)
        }
        # return quantiles of the null distribution
        perm_lower = quantile(permuted_cor, 0.025)
        perm_upper = quantile(permuted_cor, 0.975)
        perm_95 = quantile(permuted_cor, 0.95)
        # needed elements in output: estimates of rho_0, rho_1, and rho_2, if estimates for rho_0/rho_1 and rho_2
        saveRDS(list(rho_true = c(rho_true_s = cor_true, rho_true_p = cor_true_p), rho_0 = cor_0_p, rho_2 = cor_2, rho_2_simY = cor_2_simY,   res = res, res_rho_2 = res_rho_2, perm_lower = perm_lower, perm_upper = perm_upper, perm_95 = perm_95),paste0(save_path,n,"/simResults_",sim,".rds"))
      } else{
        # needed elements in output: estimates of rho_0, rho_1, and rho_2, if estimates for rho_0/rho_1 and rho_2
        saveRDS(list(rho_true = c(rho_true_s = cor_true, rho_true_p = cor_true_p), rho_0 = cor_0_p, rho_2 = cor_2, rho_2_simY = cor_2_simY, res = res, res_rho_2 = res_rho_2),paste0(save_path,n,"/simResults_",sim,".rds"))
      }
      
      }
    , mc.cores=ncores, mc.set.seed = TRUE
    )
  }
}


