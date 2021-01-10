## Lifetime demographic events for individual performance
library(rstan)
library(parallel)

# project an existing population forward one timestep
# note n_t is a distribution over the domaine Z
one_step = function(n_t, sur_year, sur_height, sur_G, grow_year, grow_height, grow_G, growth_sigma,
  cent_const_sur, Z, dz){

  # calculate survival for each z
  sur_lp = sur_year + sur_height * (Z - cent_const_sur) + sur_G
  S = 1 / (1 + exp(-sur_lp))

  # do the growth
  mean_growth = grow_year + grow_height * Z + grow_G
  mus = expand.grid(Z, mean_growth)
  height_dists = matrix(dnorm(mus[, 1], mus[, 2], growth_sigma),
    nrow = length(Z), byrow = TRUE)

  # combine the two
  n_t1 = colSums(n_t * S * height_dists) * dz

  return(n_t1)

}
ind_projection = function(n0, sur_year, sur_height, sur_G, grow_year, grow_height, grow_G, growth_sigma,
  cent_const_sur, Z, dz){

  time_hor = length(sur_year)
  pop_t = matrix(NA, nrow = time_hor, ncol = length(Z))

  # starting point
  i = 1
  pop_t[i, ] = one_step(n0, sur_year[i], sur_height, sur_G,
    grow_year[i], grow_height[i], grow_G, growth_sigma, cent_const_sur, Z, dz)
  # projection through time
  for(i in 2:time_hor){

    pop_t[i, ] = one_step(pop_t[i - 1,], sur_year[i], sur_height, sur_G,
      grow_year[i], grow_height[i], grow_G, growth_sigma, cent_const_sur, Z, dz)

  }

  return(pop_t)

}

# probability of survuival
pr_sur = function(proj, dz){

  return(rowSums(proj) * dz)

}

# probability of flowering
pr_rep_z = function(rep_year, rep_height, rep_G, cent_const, Z){

  time_hor = length(rep_year)
  res = matrix(NA, nrow = time_hor, ncol = length(Z))

  for(i in 1:time_hor){

    lp = rep_year[i] + rep_height * (Z - cent_const) + rep_G
    res[i, ] = 1 / (1 + exp(-lp))

  }

  return(res)

}

pr_rep = function(rep_Z_T, pop_proj, sur_proj, dz){

  return(rowSums((pop_proj / sur_proj) * rep_Z_T) * dz)

}

# number lifetime reproductive events
lifetime_Rn = function(rep_Z_T, pop_proj, dz){

  return(sum(rep_Z_T * pop_proj) * dz)

}

# age of first reproduction
age_R = function(rep_T){

  time_hor = length(rep_T)
  CDF_rep = numeric(time_hor)

  for(i in 1:time_hor){

    CDF_rep = 1 - prod(1 - rep_T[1:i])

  }

  return(CDF_rep)

}

# expected fecundity at each size in Z, take the modal number of seeds
fec_mode = function(int_f, beta_f, sigma_f, cent_const, Z){

  return(exp((int_f + beta_f * (Z - cent_const)) - sigma_f))

}

# calc study period fruit production using the expectation in IPM framework
R_E = function(n0, sur_year, sur_height, sur_G, grow_year, grow_height, grow_G, growth_sigma,
  rep_year, rep_height, rep_G, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep, cent_const_sur,
  Z, dz){

  #fruit in each year
  fruit = numeric(length(rep_year))

  # get the modal number of seeds predicted for each size
  fec_Z = fec_mode(int_f, beta_f, sigma_f, cent_const_fec, Z)

  # get the probability of flowering in each year for each height z
  rep_mat = pr_rep_z(rep_year, rep_height, rep_G, cent_const_rep, Z)

  # get the first year, since we know an indivudal is alive in the first year
  fruit[1] = sum(n0 * rep_mat[1, ] * fec_Z) * dz

  # project the population foward
  pop_pro = ind_projection(n0, sur_year, sur_height, sur_G, grow_year, grow_height, grow_G, growth_sigma,
    cent_const_sur, Z, dz)

  # calc expected fruit number in each year
  for(j in 1:length(sur_year)){

    fruit[j + 1] = sum(pop_pro[j, ] * rep_mat[j + 1, ] * fec_Z) * dz

  }

  return(sum(fruit))

}

# simulate an individual over period j0 to J to get number of fruit produced over that period
# int_sur, int_gr, int_rep, beta_gr are all vectors with one estimate per year simulated over
# -1 (except for int_rep since reproduciton happens in the firt year
R_sim = function(int_sur, beta_sur, SPP_sur, int_gr, beta_gr, SPP_gr, sigma_gr, int_rep,
  beta_rep, SPP_rep, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep, cent_const_sur, z0){

  # number of years to simulate over
  TT = length(int_rep)
  # hold the number fruit in each year
  Rt = numeric(TT)

  # prob reproduction
  lp_R = int_rep[1] + beta_rep * (z0 - cent_const_rep) + SPP_rep
  R = 1 / (1 + exp(-lp_R))

  # first year fruit
  Rt[1] = rbinom(1, 1, R) * rlnorm(1, int_f +  beta_f * (z0 - cent_const_fec), sigma_f)

  # step through time
  z = z0
  for(j in 2:TT){

    # get next size
    z = rnorm(1, int_gr[j - 1] + beta_gr[j - 1] * z + SPP_gr , sigma_gr)

    # get survival
    lp_S = int_sur[j - 1] + beta_sur * (z - cent_const_sur) + SPP_sur
    S = 1 / (1 + exp(-lp_S))
    S_bin = rbinom(1, 1, S)
    #stop the sim if individual dies
    if(S_bin == 0) break

    # get reproduction
    lp_R = int_rep[j] + beta_rep * (z - cent_const_rep) + SPP_rep
    R = 1 / (1 + exp(-lp_R))

    Rt[j] = rbinom(1, 1, R) * rlnorm(1, int_f +  beta_f * (z - cent_const_fec), sigma_f)

  }

  return(sum(Rt))

}

# asymototic simulation, gives R0 for a given individual assuming an average year
# and an average intial size, so the only difference comes from the spatial effect.
# calculated in the IPM framework
R_asy = function(rep_int, rep_slope, cent_rep, rep_G, int_f, slope_f, sigma_f, cent_f,
	sur_int, sur_slope, cent_sur, sur_G, grow_int, grow_slope, grow_G, growth_sigma,
	n0, Z){

	# make the fecundity kernel, vector of seeds produced by each size in Z
	# prob of flowering
    lp = rep_int + rep_slope * (Z - cent_rep) + rep_G
    prob_flower = 1 / (1 + exp(-lp))
	# size effect on seed production by the probability of flowering
	Fec = fec_mode(int_f, slope_f, sigma_f, cent_f, Z) * prob_flower

	# make the growth-survival kernel, matrix of |Z| by |Z|
	# calulualte surival for each size in Z
	sur_lp = sur_int + sur_slope * (Z - cent_sur) + sur_G
  	S = 1 / (1 + exp(-sur_lp))

	# do the growth
	mean_growth = grow_int + grow_slope * Z + grow_G
	mu = expand.grid(Z, mean_growth)
	G = matrix(dnorm(mu[, 1], mu[, 2], growth_sigma),
	    nrow = length(Z), byrow = TRUE)
	# make growth conditional on survival
	P = G * S

	# construct the identiy matrix
	I = diag(length(Z))

	# multiply fecundity by Neumann series by intial population eq. 3.1.8 from IPM book
	return((Fec %*% solve(I - P) %*% n0))

}

# function to take the sampled parameter values from the stan objects and produce a
# distribution of number of fruit
R_samp_space = function(sur_stan, gr_stan, rep_stan, fec_stan, num_resample, z0_mu, z0_sigma,
  cent_const_fec, cent_const_rep, cent_const_sur, loc_ind, loc_ind_list){

  # map loc_ind for each vital rate to the right index using the list
  sur_loc_ind = loc_ind_list[[loc_ind]]$sur_ind
  rep_loc_ind = loc_ind_list[[loc_ind]]$rep_ind
  gr_loc_ind = loc_ind_list[[loc_ind]]$gr_ind

  # extract the parameters from the stan objects as dataframes
  sur_pars = data.frame(extract_flat(sur_stan, pars = c('h_ef', 'year_int', paste0('spp[', sur_loc_ind, ']'))))
  gr_pars = data.frame(extract_flat(gr_stan, pars = c('b0', 'gr_rate', paste0('spp[', gr_loc_ind, ']'), 'sigma')))
  rep_pars = data.frame(extract_flat(rep_stan, pars = c('h_ef', 'year_int', paste0('spp[', rep_loc_ind, ']'))))
  fec_pars = data.frame(extract_flat(fec_stan, pars = c('h_ef', 'site_int[2]', 'sigma'))) # site 2 is site being estimated for

  # get the indicies for the rows to num_resample
  inds = sample.int(length(sur_pars[, 1]), num_resample)

  RJJ = numeric(num_resample)

  for(i in 1:num_resample){

    # get the parameters
    int_sur = as.numeric(sur_pars[inds[i], 2:6])
    beta_sur = as.numeric(sur_pars[inds[i], 1])
    SPP_sur = as.numeric(sur_pars[inds[i], 7])

    int_gr = as.numeric(gr_pars[inds[i], 1:5])
    beta_gr = as.numeric(gr_pars[inds[i], 6:10])
    SPP_gr = as.numeric(gr_pars[inds[i], 11])
    sigma_gr = as.numeric(gr_pars[inds[i], 12])

    int_rep = as.numeric(rep_pars[inds[i], 2:7])
    beta_rep = as.numeric(rep_pars[inds[i], 1])
    SPP_rep = as.numeric(rep_pars[inds[i], 8])

    int_f = as.numeric(fec_pars[inds[i], 2])
    beta_f = as.numeric(fec_pars[inds[i], 1])
    sigma_f = as.numeric(fec_pars[inds[i], 3])

    # draw the intial size from a distribution becuase we want to show spatial
    # structure, and we want that structure to be independent of spatial structure in
    # starting size
    z0 = rlnorm(1, z0_mu, z0_sigma)

    RJJ[i] = R_sim(int_sur, beta_sur, SPP_sur, int_gr, beta_gr, SPP_gr, sigma_gr, int_rep,
      beta_rep, SPP_rep, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep, cent_const_sur, z0)
  }

  return(RJJ)

}

R_samp_nospace = function(sur_stan, gr_stan, rep_stan, fec_stan, num_resample, z0_mu, z0_sigma,
  cent_const_fec, cent_const_rep, cent_const_sur){

  # extract the parameters from the stan objects as dataframes
  sur_pars = data.frame(extract_flat(sur_stan, pars = c('h_ef', 'year_int')))
  gr_pars = data.frame(extract_flat(gr_stan, pars = c('b0', 'gr_rate', 'sigma')))
  rep_pars = data.frame(extract_flat(rep_stan, pars = c('h_ef', 'year_int')))
  fec_pars = data.frame(extract_flat(fec_stan, pars = c('h_ef', 'site_int[2]', 'sigma'))) # site 2 is the site being estimated for

  # get the indicies for the rows to num_resample
  inds = sample.int(length(sur_pars[, 1]), num_resample)

  RJJ = numeric(num_resample)

  for(i in 1:num_resample){


    # get the parameters
    int_sur = as.numeric(sur_pars[inds[i], 2:6])
    beta_sur = as.numeric(sur_pars[inds[i], 1])
    SPP_sur = 0

    int_gr = as.numeric(gr_pars[inds[i], 1:5])
    beta_gr = as.numeric(gr_pars[inds[i], 6:10])
    SPP_gr = 0
    sigma_gr = as.numeric(gr_pars[inds[i], 11])

    int_rep = as.numeric(rep_pars[inds[i], 2:7])
    beta_rep = as.numeric(rep_pars[inds[i], 1])
    SPP_rep = 0

    int_f = as.numeric(fec_pars[inds[i], 2])
    beta_f = as.numeric(fec_pars[inds[i], 1])
    sigma_f = as.numeric(fec_pars[inds[i], 3])

    # draw the intial size from a distribution becuase we want to show spatial
    # structure, and we want that structure to be independent of spatial structure in
    # starting size
    z0 = rlnorm(1, z0_mu, z0_sigma)

    RJJ[i] = R_sim(int_sur, beta_sur, SPP_sur, int_gr, beta_gr, SPP_gr, sigma_gr, int_rep,
      beta_rep, SPP_rep, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep, cent_const_sur, z0)
  }

 return(RJJ)

}

# function to resample the expectations to include parameter uncertianty
R_E_samp = function(sur_stan, gr_stan, rep_stan, fec_stan, Z, dz, z0_mu, z0_sigma,
  cent_const_fec, cent_const_rep, cent_const_sur, loc_ind, loc_ind_list, num_resample){

  # map loc_ind for each vital rate to the right index using the list
  sur_loc_ind = loc_ind_list[[loc_ind]]$sur_ind
  rep_loc_ind = loc_ind_list[[loc_ind]]$rep_ind
  gr_loc_ind = loc_ind_list[[loc_ind]]$gr_ind

  # extract the parameters from the stan objects as dataframes
  sur_pars = data.frame(extract_flat(sur_stan, pars = c('h_ef', 'year_int', paste0('spp[', sur_loc_ind, ']'))))
  gr_pars = data.frame(extract_flat(gr_stan, pars = c('b0', 'gr_rate', paste0('spp[', gr_loc_ind, ']'), 'sigma')))
  rep_pars = data.frame(extract_flat(rep_stan, pars = c('h_ef', 'year_int', paste0('spp[', rep_loc_ind, ']'))))
  fec_pars = data.frame(extract_flat(fec_stan, pars = c('h_ef', 'site_int[2]', 'sigma'))) # site 2 is site being estimated for

  # get the indicies for the rows to num_resample
  inds = sample.int(length(sur_pars[, 1]), num_resample)

  RJJ = numeric(num_resample)
  sur_G = numeric(num_resample)
  rep_G = numeric(num_resample)
  gr_G = numeric(num_resample)

  for(i in 1:num_resample){

    # get the parameters
    int_sur = as.numeric(sur_pars[inds[i], 2:6])
    beta_sur = as.numeric(sur_pars[inds[i], 1])
    SPP_sur = as.numeric(sur_pars[inds[i], 7])

    int_gr = as.numeric(gr_pars[inds[i], 1:5])
    beta_gr = as.numeric(gr_pars[inds[i], 6:10])
    SPP_gr = as.numeric(gr_pars[inds[i], 11])
    sigma_gr = as.numeric(gr_pars[inds[i], 12])

    int_rep = as.numeric(rep_pars[inds[i], 2:7])
    beta_rep = as.numeric(rep_pars[inds[i], 1])
    SPP_rep = as.numeric(rep_pars[inds[i], 8])

    int_f = as.numeric(fec_pars[inds[i], 2])
    beta_f = as.numeric(fec_pars[inds[i], 1])
    sigma_f = as.numeric(fec_pars[inds[i], 3])

    # draw the intial size from a distribution becuase we want to show spatial
    # structure, and we want that structure to be independent of spatial structure in
    # starting size
    n0 = dlnorm(Z, z0_mu, z0_sigma)

    RJJ[i] = R_E(n0, int_sur, beta_sur, SPP_sur, int_gr, beta_gr, SPP_gr, sigma_gr,
      int_rep, beta_rep, SPP_rep, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep,
      cent_const_sur, Z, dz)

    sur_G[i] = SPP_sur
    rep_G[i] = SPP_rep
    gr_G[i] = SPP_gr

  }

  return(list(RJJ = RJJ, sur_spp = sur_G, rep_spp = rep_G, gr_spp = gr_G))

}

# fruit number under the null hypothesis that vital rates at each locatin are independent
R_E_null = function(sur_stan, gr_stan, rep_stan, fec_stan, Z, dz, z0_mu, z0_sigma,
  cent_const_fec, cent_const_rep, cent_const_sur, loc_ind_list, num_resample, num_samp_loc){


  sur_pars = data.frame(extract_flat(sur_stan, pars = c('h_ef', 'year_int', paste0('spp[', 1, ']'))))
  # get the indicies for the rows to num_resample
  inds = sample.int(length(sur_pars[, 1]), num_resample, replace = TRUE)

  RJJ = numeric(num_resample)
  sur_G = numeric(num_resample)
  rep_G = numeric(num_resample)
  gr_G = numeric(num_resample)

  max_ind = length(loc_ind_list)

  suff_sur_inds = sample.int(max_ind, num_resample, replace = TRUE)
  suff_rep_inds = sample.int(max_ind, num_resample, replace = TRUE)
  suff_gr_inds = sample.int(max_ind, num_resample, replace = TRUE)

  for(i in 1:num_resample){

    if((i %% 500) == 0) print(i / num_resample)

    # map loc_ind for each vital rate to the right index using the list
    sur_loc_ind = loc_ind_list[[ suff_sur_inds[i] ]]$sur_ind
    rep_loc_ind = loc_ind_list[[ suff_rep_inds[i] ]]$rep_ind
    gr_loc_ind = loc_ind_list[[ suff_gr_inds[i] ]]$gr_ind

    # extract the parameters from the stan objects as dataframes
    sur_pars = data.frame(extract_flat(sur_stan, pars = c('h_ef', 'year_int', paste0('spp[', sur_loc_ind, ']'))))
    gr_pars = data.frame(extract_flat(gr_stan, pars = c('b0', 'gr_rate', paste0('spp[', gr_loc_ind, ']'), 'sigma')))
    rep_pars = data.frame(extract_flat(rep_stan, pars = c('h_ef', 'year_int', paste0('spp[', rep_loc_ind, ']'))))
    fec_pars = data.frame(extract_flat(fec_stan, pars = c('h_ef', 'site_int[2]', 'sigma'))) # site 2 is site being estimated for

    # get the parameters
    int_sur = as.numeric(sur_pars[inds[i], 2:6])
    beta_sur = as.numeric(sur_pars[inds[i], 1])
    SPP_sur = as.numeric(sur_pars[inds[i], 7])

    int_gr = as.numeric(gr_pars[inds[i], 1:5])
    beta_gr = as.numeric(gr_pars[inds[i], 6:10])
    SPP_gr = as.numeric(gr_pars[inds[i], 11])
    sigma_gr = as.numeric(gr_pars[inds[i], 12])

    int_rep = as.numeric(rep_pars[inds[i], 2:7])
    beta_rep = as.numeric(rep_pars[inds[i], 1])
    SPP_rep = as.numeric(rep_pars[inds[i], 8])

    int_f = as.numeric(fec_pars[inds[i], 2])
    beta_f = as.numeric(fec_pars[inds[i], 1])
    sigma_f = as.numeric(fec_pars[inds[i], 3])

    # draw the intial size from a distribution becuase we want to show spatial
    # structure, and we want that structure to be independent of spatial structure in
    # starting size
    n0 = dlnorm(Z, z0_mu, z0_sigma)

    RJJ[i] = R_E(n0, int_sur, beta_sur, SPP_sur, int_gr, beta_gr, SPP_gr, sigma_gr,
      int_rep, beta_rep, SPP_rep, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep,
      cent_const_sur, Z, dz)

    sur_G[i] = SPP_sur
    rep_G[i] = SPP_rep
    gr_G[i] = SPP_gr

  }

  # return a structure like the observed data so can be analysied in the same way
  return(list(RJJ = matrix(RJJ, ncol = num_samp_loc, byrow = TRUE),
    sur_spp = matrix(sur_G, ncol = num_samp_loc, byrow = TRUE),
    rep_spp = matrix(rep_G, ncol = num_samp_loc, byrow = TRUE),
    gr_spp = matrix(gr_G, ncol = num_samp_loc, byrow = TRUE)))

}


# slightly differnt version of the randomization where only one shuffel is done, used so that by location summaries
# like means, are comparable to R0 generated from the fitted vital rate models

R_E_null_struct = function(sur_stan, gr_stan, rep_stan, fec_stan, Z, dz, z0_mu, z0_sigma,
  cent_const_fec, cent_const_rep, cent_const_sur, loc_ind_list, num_resample){


  num_locs = length(loc_ind_list)

  RJJ = matrix(NA, ncol = num_resample, nrow = num_locs )
  sur_G = matrix(NA, ncol = num_resample, nrow = num_locs )
  rep_G = matrix(NA, ncol = num_resample, nrow = num_locs )
  gr_G = matrix(NA, ncol = num_resample, nrow = num_locs )

  suff_sur_inds = sample.int(num_locs, num_locs, replace = FALSE)
  suff_rep_inds = sample.int(num_locs, num_locs, replace = FALSE)
  suff_gr_inds = sample.int(num_locs, num_locs, replace = FALSE)

  for(j in 1:num_locs){

    if((j %% 50) == 0) print(j / num_locs)

    # map loc_ind for each vital rate to the right index using the list
    sur_loc_ind = loc_ind_list[[ suff_sur_inds[j] ]]$sur_ind
    rep_loc_ind = loc_ind_list[[ suff_rep_inds[j] ]]$rep_ind
    gr_loc_ind = loc_ind_list[[ suff_gr_inds[j] ]]$gr_ind

    # extract the parameters from the stan objects as dataframes
    sur_pars = data.frame(extract_flat(sur_stan, pars = c('h_ef', 'year_int', paste0('spp[', sur_loc_ind, ']'))))
    gr_pars = data.frame(extract_flat(gr_stan, pars = c('b0', 'gr_rate', paste0('spp[', gr_loc_ind, ']'), 'sigma')))
    rep_pars = data.frame(extract_flat(rep_stan, pars = c('h_ef', 'year_int', paste0('spp[', rep_loc_ind, ']'))))
    fec_pars = data.frame(extract_flat(fec_stan, pars = c('h_ef', 'site_int[2]', 'sigma'))) # site 2 is site being estimated for

    # get the indicies for the rows to num_resample
    inds = sample.int(length(sur_pars[, 1]), num_resample, replace = FALSE)

    for(i in 1:num_resample){

      # get the parameters
      int_sur = as.numeric(sur_pars[inds[i], 2:6])
      beta_sur = as.numeric(sur_pars[inds[i], 1])
      SPP_sur = as.numeric(sur_pars[inds[i], 7])

      int_gr = as.numeric(gr_pars[inds[i], 1:5])
      beta_gr = as.numeric(gr_pars[inds[i], 6:10])
      SPP_gr = as.numeric(gr_pars[inds[i], 11])
      sigma_gr = as.numeric(gr_pars[inds[i], 12])

      int_rep = as.numeric(rep_pars[inds[i], 2:7])
      beta_rep = as.numeric(rep_pars[inds[i], 1])
      SPP_rep = as.numeric(rep_pars[inds[i], 8])

      int_f = as.numeric(fec_pars[inds[i], 2])
      beta_f = as.numeric(fec_pars[inds[i], 1])
      sigma_f = as.numeric(fec_pars[inds[i], 3])

      # draw the intial size from a distribution becuase we want to show spatial
      # structure, and we want that structure to be independent of spatial structure in
      # starting size
      n0 = dlnorm(Z, z0_mu, z0_sigma)

      RJJ[j, i] = R_E(n0, int_sur, beta_sur, SPP_sur, int_gr, beta_gr, SPP_gr, sigma_gr,
	int_rep, beta_rep, SPP_rep, int_f, beta_f, sigma_f, cent_const_fec, cent_const_rep,
	cent_const_sur, Z, dz)

      sur_G[j, i] = SPP_sur
      rep_G[j, i] = SPP_rep
      gr_G[j, i] = SPP_gr

    }

  }

  # return a structure like the observed data so can be analysied in the same way
  return(list(RJJ = RJJ, sur_spp = sur_G, rep_spp = rep_G, gr_spp = gr_G))

}


# little helper function for calculate shift in dist
dist_shift = function(dist_1, dist_2, X, dx){
  CDF1 = ecdf(dist_1)
  CDF2 = ecdf(dist_2)

  return(sum(CDF1(X)) * dx - sum(CDF2(X)) * dx)

}

rand_shift_wrapper = function(pl, sur_stan, gr_stan, rep_stan, fec_stan, loc_ind_list){

  print(pl$rand_ind)

  null_hold = R_E_null_struct(sur_stan, gr_stan, rep_stan, fec_stan, pl$Z, pl$dz, pl$z0_mu, pl$z0_sigma,
    pl$cent_const_fec, pl$cent_const_rep, pl$cent_const_sur, loc_ind_list, pl$num_resample)

  null_dist = apply(null_hold$RJJ, MARGIN = 1, FUN = median)

  return(dist_shift(pl$obs_Rf, null_dist, pl$X, pl$dx))

}

# randoimisation of the shift in dist of fruit number over space between observed and structured null
space_RJJ_shift_rand = function(sur_stan, gr_stan, rep_stan, fec_stan, Z, dz, z0_mu,
  z0_sigma, cent_const_fec, cent_const_rep, cent_const_sur, loc_ind_list, num_resample,
  n_cores, n_rands, obs_Rf, X, dx){

  dum_vec = 1:n_rands

  par_list = list()
  for(i in 1:n_rands){

    par_list[[i]] = list(Z = Z, dz = dz, z0_mu = z0_mu, z0_sigma = z0_sigma, cent_const_fec = cent_const_fec,
      cent_const_rep = cent_const_rep, cent_const_sur = cent_const_sur, num_resample = num_resample,
      obs_Rf = obs_Rf, X = X, dx = dx, rand_ind = i)

  }

  print(par_list)

  out = mclapply(par_list, FUN = rand_shift_wrapper, sur_stan = sur_stan, gr_stan = gr_stan,
    rep_stan = rep_stan, fec_stan = fec_stan, loc_ind_list = loc_ind_list, mc.cores = n_cores)

    return(out)

}

#########################################################################################################################
# DO I need these??


# function to fit a model to vital rates vs R0. Each R0 and vr is sampled from a postiroer and or a randimisation
# with gap structured random effects
lmer_fitter = function(form, res_pred_list, ranef_IDs, ...){

  out = list()
  num_locs = dim(res_pred_list$RJJ)[1]
  num_samp_loc = dim(res_pred_list$RJJ)[2]

  for(i in 1:num_samp_loc){

    if((i %% 10) == 0) print(i / num_samp_loc)

    # build a data frame with the preedictors and response, along with gapID
    dat = data.frame(RJJ = res_pred_list$RJJ[, i], sur_spp = res_pred_list$sur_spp[, i],
      rep_spp = res_pred_list$rep_spp[, i], gr_spp = res_pred_list$gr_spp[, i], gapID = ranef_IDs)

    out[[i]] = lmer(form, data = dat, ...)

  }

  return(out)

}
