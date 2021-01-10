# Spatial growth model, 
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
gr_dat = vr_loc_gen[!is.na(vr_loc_gen$height), c('year', 'height', 'height_prev', 'X', 'Y')]
gr_dat = gr_dat[!is.na(gr_dat$height_prev), ]

# create knots at each data point, cluster some less than 50cm apart
# results in 575 knots. 50% of knots have another knot within 61cm and 75% have
# another knot within 0.82cm. 50% of knots have 8 or more neigbour knots within 2m
# and 75% have 5 or more neighbouring knots within 2m.
knot_cl = knot_coords2(dat_x = gr_dat$X, dat_y = gr_dat$Y, min_dist = 0.5)

dm_knots = distance_matrix(knot_cl[,1], knot_cl[,2])
dm_knots_obs = distance_matrix2(x1 = gr_dat$X, y1 = gr_dat$Y, x2 = knot_cl[, 1],
    y2 = knot_cl[, 2])
# All observations have a knot within 56cm, and 75% have a knot within 23cm, all looks very nice

# some quick ckecks on the knot pattern
# dim(dm_knots)
# int_knot_dist = apply(dm_knots, MARGIN = 1, FUN = function(x) sort(x, FALSE)[2])
# summary(int_knot_dist)
# knot_NN = apply(dm_knots, MARGIN = 1, FUN = function(x) sum(x < 2))
# summary(knot_NN)
#
# knot_ob_dist_sum = apply(dm_knots_obs, MARGIN = 1, FUN = min)
# summary(knot_ob_dist_sum)
# hist(knot_ob_dist_sum)
# ob_knot_NN = apply(dm_knots_obs, MARGIN = 1, FUN = function(x) sum(x < 2))
# summary(ob_knot_NN)

# get index of nearerst knot for each data point
near_knot = apply(dm_knots_obs, MARGIN = 1, FUN = function(x) which(x == min(x)))

# prepare the data for stan
dd = data.frame(y = gr_dat$height,
  year = gr_dat$year,
  height_prev = gr_dat$height_prev)

stan_mod = "
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
  real height[N];  // response variable
  vector[N] height_prev;
  // data for year effect on growth rate
  int<lower=1> year[N];
  int<lower=1> N_years;
  // spatial predictive process data
  int<lower=1> N_knot; //total number of knots
  matrix[N_knot, N_knot] knot_dist; //distance matrix between observations
  int<lower=1, upper=N_knot> nn_knot_ind[N]; // index of the nearerst knot for each data point
}
transformed data {
  //might log things in here
}
parameters {
  real<lower=0> sigma; // standard deviaiton overall on height
  vector[N_years] b0; // year specifc intercept
  vector[N_years] gr_rate; // year specifc slope

  //Predictive Proccess parameters
  real<lower=0.00001> sigma_spp;
  real<lower=0.1> dd_spp_inv;
  vector[N_knot] spp_raw;
}
transformed parameters {
  // spatial error
  matrix[N_knot, N_knot] C;
  real<lower=0> dd_spp;
  vector[N_knot] spp; // spatial random effect, draw from MVnorm
  //linear predcitor
  vector[N] mu;

  //spatial error
  dd_spp = inv(dd_spp_inv);
  C = cov_var_maker(knot_dist, N_knot, dd_spp, sigma_spp);

  spp = cholesky_decompose(C) * spp_raw;

  //linear predictor fixed effects
  for (n in 1:N) {
    mu[n] = b0[year[n]] + height_prev[n] * gr_rate[year[n]] + spp[nn_knot_ind[n]];
  }
}
model {
  // prior specifications
  sigma ~ cauchy(0, 20);
  sigma_spp ~ cauchy(0, 5); 
  dd_spp_inv ~ cauchy(0, 5);

  //estimated spatial effect
  spp_raw ~ normal(0, 1);

  // likelihood contribution
  height ~ normal(mu, sigma);
}
generated quantities {
}"

input = make_standata(formula = y ~ 0 + height_prev + (1|year) ,
  data=dd)

stan_dat = list(N = input$N, height = input$Y, height_prev = dd$height_prev,
  N_years = input$N_1, year = input$J_1, N_knot = dim(dm_knots)[1], knot_dist = dm_knots,
  nn_knot_ind = near_knot)

SPP_gr_stan <- stan(model_code = stan_mod, data = stan_dat, iter = 3000, warmup = 1000,
  cores = 4, chains = 4, control = list(adapt_delta = 0.99),
  pars = c('b0', 'gr_rate', 'sigma', 'dd_spp_inv', 'sigma_spp', 'spp', 'mu'))

save(SPP_gr_stan , file = 'SPP_NK_gr_stan.Rdata')
#load('SPP_NK_gr_stan.Rdata')

# create a set of diagnostics
pdf('SPP_NK_gr_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_gr_stan, inc_warmup = FALSE, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP_gr_stan, inc_warmup = FALSE, pars = c(paste0('spp[', 1:50, ']')))
  stan_trace(SPP_gr_stan, inc_warmup = FALSE, pars = c(paste0('spp[', 400:450, ']')))
  stan_trace(SPP_gr_stan, inc_warmup = TRUE, pars = paste0('mu[', 1:50, ']'))
  stan_trace(SPP_gr_stan, inc_warmup = TRUE, pars = paste0('mu[', 1000:1050, ']'))
  stan_rhat(SPP_gr_stan, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP_gr_stan, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv'))
  hist_treedepth(SPP_gr_stan)
  pairs(SPP_gr_stan, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv'))
dev.off()