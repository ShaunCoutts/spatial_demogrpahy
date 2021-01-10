## Nearest Knot version of the SPP model. With  a careful knot placement stratergy
## This means I do not have to the back interpolation step. This is computationaly
## expensive and eats RAM. This is pretty close to
## the vinailla version of a spatial errors model, but is much faster to run
## since the distance matrix is smaller and easier to invert or decompose, and 
## faster to build.

require(brms)
library(rstan)
library(here) # used to construct a path to you local version of the scripts 

hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
} 


# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

source('dist_neigh_setup.R')

# get data
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactor = FALSE)
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

# start with survival
sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]
# neighbour data we can use all data, even the first years observed height 
neigh_dat = sur_dat

# take out data from the first year observed since nothing can be observed dead in the first year
first_year = sapply(seq_along(sur_dat$year), FUN = function(x){

  min_year_group = min(sur_dat$year[sur_dat$uID == sur_dat$uID[x]])
  return(ifelse(sur_dat$year[x] == min_year_group, FALSE, TRUE))

})
sur_dat = sur_dat[first_year, ]

# recode year and patch to act as indicies 
sur_dat$year_num = (sur_dat$year - (max(sur_dat$year))) + max(abs(sur_dat$year - (max(sur_dat$year)))) + 1    
sur_dat$gapID = sapply(strsplit(sur_dat$uID, split = ':', fixed = TRUE), FUN = function(x) x[1])
inds = unique(sur_dat$gapID)
sur_dat$gapID_num = NA
for(i in 1:length(inds)) sur_dat$gapID_num[sur_dat$gapID %in% inds[i]] = i

# height: this requires some assumptions because 1) individuals change size over the course of a year, 
# 2) we do not have the final size for indivudals that died. 3) we have indivduals that were only found when they were very large.
# Best we can do here is assue that mortality occurs just after census, so that we can use previous height as a predictor,
# and in the case of indivduals in their first year we say they were 0 
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]
sur_dat$height_mod = sur_dat$height_prev - mean(sur_dat$height_prev)

# create knots at each data point, cluster some less than 50cm apart
# results in 575 knots. 50% of knots have another knot within 61cm and 75% have
# another knot within 0.82cm. 50% of knots have 8 or more neigbour knots within 2m
# and 75% have 5 or more neighbouring knots within 2m.
knot_cl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5) 

dm_knots = distance_matrix(knot_cl[,1], knot_cl[,2])
dm_knots_obs = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knot_cl[, 1], 
  y2 = knot_cl[, 2])
# All observations have a knot within 56cm, and 75% have a knot within 23cm
  
# some quick ckecks on the knot pattern
#dim(dm_knots)
#int_knot_dist = apply(dm_knots, MARGIN = 1, FUN = function(x) sort(x, FALSE)[2])
#summary(int_knot_dist)
#knot_NN = apply(dm_knots, MARGIN = 1, FUN = function(x) sum(x < 2))
#summary(knot_NN)

#knot_ob_dist_sum = apply(dm_knots_obs, MARGIN = 1, FUN = min)
#summary(knot_ob_dist_sum)
#hist(knot_ob_dist_sum)
#ob_knot_NN = apply(dm_knots_obs, MARGIN = 1, FUN = function(x) sum(x < 2))
#summary(ob_knot_NN)

# get index of nearerst knot for each data point
near_knot = apply(dm_knots_obs, MARGIN = 1, FUN = function(x) which(x == min(x)))

# prepare the data for stan
dd = data.frame(y = sur_dat$sur,
  year = sur_dat$year,
  height = sur_dat$height_mod)

stan_model2 = "
// We recommend generating the data with the 'make_standata' function. 
functions {
 
  //constructs a correlation matrix
  matrix cov_var_maker(matrix dist_mat, int mat_size, real dd_spp, real sigma_spp){
  
    matrix[mat_size, mat_size] C;
    
    for(k in 1:(mat_size - 1)){
      for(j in (k + 1):mat_size){
	
	C[k, j] = sigma_spp * exp(-(dd_spp * dist_mat[k, j])); 
	C[j, k] = C[k, j];
	
      }
    }
    //fill diagonals
    for(k in 1:mat_size){
      C[k, k] = sigma_spp;
    }
    
    return C;
    
  }

} 
data { 
  int<lower=1> N;  // total number of observations
  int<lower=1> N_knot; //total number of knots
  int sur[N];  // response variable 
  vector[N] height_cen;  // population-level design matrix 
  // data for group-level effects of ID 1 
  int<lower=1> year[N]; 
  int<lower=1> N_years;  
  int<lower=1> M_1;
  matrix[N_knot, N_knot] knot_dist; //distance matrix between observations
  int<lower=1, upper=N_knot> nn_knot_ind[N]; // index of the nearerst knot for each data point
} 
transformed data { 
} 
parameters { 
  real h_ef;  // population-level effects 
  vector[N_years] z_1[M_1];  // unscaled group-level effects
  
  //Predictive Proccess parameters 
  real<lower=0.00001> sigma_spp;
  real<lower=0.1, upper=1000> dd_spp_inv;
  vector[N_knot] spp_raw;
  
} 
transformed parameters { 
  // spatial error
  matrix[N_knot, N_knot] C;
  real<lower=0> dd_spp;
  vector[N_knot] spp; // spatial random effect, draw from MVnorm
  //year effect
  vector[N_years] year_int;
  //linear predcitor
  vector[N] eta; 
  
  //spatial error 
  dd_spp = inv(dd_spp_inv);
  C = cov_var_maker(knot_dist, N_knot, dd_spp, sigma_spp);
  
  spp = cholesky_decompose(C) * spp_raw; 
 
  //year effect
  year_int = 20.0 * (z_1[1]);// prior of N(0, 20), sd= 20 on logit scale is effectivly from 0-1 prob on real scale.
  
  //linear predictor fixed effects 
  eta = height_cen * h_ef; 
  // add both random effect of year and space
  for (n in 1:N) { 
    eta[n] = eta[n] + year_int[year[n]] + spp[nn_knot_ind[n]];
  } 
} 
model {
  // prior specifications 
  sigma_spp ~ cauchy(0, 5);
  dd_spp_inv ~ cauchy(0, 5);
  z_1[1] ~ normal(0, 1);
  
  //estimated spatial effect 
  spp_raw ~ normal(0, 1);
  
  // likelihood contribution 
  sur ~ bernoulli_logit(eta); 
} 
generated quantities { 
}"              

input = make_standata(formula = y ~ 0 + height + (1|year) ,
  data=dd, family="bernoulli")

stan_dat = list(N = input$N, N_knot = dim(dm_knots)[1], sur = input$Y,  
  height_cen = dd$height, N_years = input$N_1, year = input$J_1, knot_dist = dm_knots, M_1 = input$M_1,
  nn_knot_ind = near_knot)

SPP_NK_sur_stan <- stan(model_code = stan_model2, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.99), 
  pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv', 'spp', 'eta')) 

save(SPP_NK_sur_stan , file = 'SPP_NK_sur_stan.Rdata')
#load('SPP_NK_sur_stan.Rdata')

# diagnostic plots
pdf('SPP_NK_sur_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_NK_sur_stan, inc_warmup = FALSE, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = FALSE, pars = paste0('spp[', 1:100, ']'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = FALSE, pars = paste0('spp[', 400:500, ']'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP_NK_sur_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP_NK_sur_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'), binwidth = 0.05)
  hist_treedepth(SPP_NK_sur_stan)
  pairs(SPP_NK_sur_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
dev.off()
