# set of functions to work out diganostics and model performance measures
## observed vs predicted for binary response, use 2Dhist to give density and jitter to visulise uncertianty
obs_v_pred_bin = function(obs, lp_samp, num_resamp, jitter_am, labels = labs(list(title = 'title', x = 'x_lab', y = 'y_lab'))){

  obs_jit = jitter(obs, amount = jitter_am)
  # now need to expand each obs so there are num_resamp predictions for each obs, the -1 dummy data is to trick the image function
  # into coloring the whole mat
  rs_row = sample.int(dim(lp_samp)[1], num_resamp)

  pred_obs = data.frame(ob = rep(obs_jit, each = num_resamp),
    pred_lp = as.vector(apply(lp_samp, MARGIN = 2, FUN = function(x) x[rs_row])))
  pred_obs$p_pred = 1 / (1 + exp(-pred_obs$pred_lp))

  p <- ggplot(pred_obs, aes(p_pred, ob))
  h3 <- p + stat_bin2d() + theme(legend.position = "none") + xlim(0, 1) + ylim(0, 1) + labels
  h3

}

## wrapper on rstan extract to flatten but not permute the output of extract
extract_flat = function(rstan_ob, ...){

  ext = rstan::extract(rstan_ob, permuted = FALSE, inc_warmup = FALSE, ...)
  ext_flat = sapply(1:dim(ext)[3], FUN = function(x) as.numeric(ext[, , x]))
  colnames(ext_flat) <- colnames(ext[1, , ])

  return(ext_flat)
}

## calculate the logLik from a set of postieriers of the linear predictor of a model
# post_mat is a num_samp by num_obs matrix of the postierier on the logit scale (plogis() to turn to prob)
logLik_bin = function(post_mat, obs){

  # turn the obs 1 into 0's and 0's into 1
  obs_flip = 1 - obs

  pred_prob = plogis(post_mat)

  return(log(apply(abs(obs_flip - t(pred_prob)), MARGIN = 1, FUN = mean)))

}

logLik_cont = function(post_mu, post_sd, obs){

  lik = numeric(length(obs))

  for(i in seq_along(obs)){
    lik[i] = mean(dnorm(obs[i], post_mu[,i], post_sd))
  }

  return(log(lik))

}
