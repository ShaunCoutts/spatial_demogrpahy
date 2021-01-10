# seed per fruit prediction
# get data (in spss)

library(foreign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(brms)
library(rstan)
library(here) # used to construct a path to you local version of the scripts 

# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

burn_dat = read.csv('hcdem_ABS_2015_checked.csv')

## try a model with better data to estimate number of fruits (assume seed production is proportional to fruit production.
# first get a clean set of data and put into long form, only take years from 1994-2003 so that most indivudals have repeated
# meaures for them so individual level effects can be fitted
fruit_ht_dat = select(burn_dat, year, burn_yr, site, gap, tag, ht94, ht95, ht96, ht97, ht98, ht99, ht00,
  ht01, ht02, ht03)
fruit_rep_dat = select(burn_dat, year, site, gap, tag, rep94, rep95, rep96, rep97, rep98, rep99, rep00, rep01, rep02,
  rep03)

ht_long = gather(fruit_ht_dat, ht_lab, height, ht94:ht03)
rep_long = gather(fruit_rep_dat, rep_lab, fruit_num, rep94:rep03)

# drop NA's
ht_long = ht_long[ht_long$height != '#NULL!', ]
rep_long = rep_long[rep_long$fruit_num != '#NULL!', ]

# make the year labels numeric
ht_long$m_year = sapply(ht_long$ht_lab, FUN = function(x){
  a = as.numeric(strsplit(x, split = 'ht')[[1]][2])
  if(a >= 50) return(1900 + a) else return(2000 + a)
})

rep_long$m_year = sapply(rep_long$rep_lab, FUN = function(x){
  a = as.numeric(strsplit(x, split = 'rep')[[1]][2])
  if(a >= 50) return(1900 + a) else return(2000 + a)
})

# add some columns for ID
ht_long = mutate(ht_long, ID = paste0(site, ':', gap, ':', tag),
  join_ID = paste0(m_year, ':', site, ':', gap, ':', tag),
  height = as.numeric(height))

rep_long = mutate(rep_long, ID = paste0(site, ':', gap, ':', tag),
  join_ID = paste0(m_year, ':', site, ':', gap, ':', tag),
  fruit_num = as.numeric(fruit_num))

# merg the data frames
fruit_dat = inner_join(ht_long, rep_long, by = 'join_ID')

fruit_dat = select(fruit_dat, burn_yr = burn_yr, site = site.x, ID = ID.x,
  join_ID = join_ID, m_year = m_year.x, height = height, fruit_num = fruit_num)
# only take the rows with fruit number greater than 1, as these are counts of fruits
fruit_dat = filter(fruit_dat, fruit_num >= 2)
# get time since fire for each observation
fruit_dat = mutate(fruit_dat, time_since_fire = m_year - burn_yr, site = as.character(site))

# explore the data a bit
plot(fruit_dat[, c('m_year', 'height', 'time_since_fire', 'fruit_num')])

site_group = group_by(fruit_dat, ID)
fruit_ID_sum = summarise(site_group, n_obs = n())
hist(fruit_ID_sum$n_obs)

# Seeds are integer number and it looks like the maybe the variance increases with
# predicted numner so try a glm using poission, possibly try a negative binomial
# as well if over dispersed. This will make it easier to relate the model
# estimates to the real scale, where they are needed.
fruit_glm = "
// We recommend generating the data with the 'make_standata' function.
functions {
}
data {
	int<lower=1> N;  // total number of observations
	int<lower=0> fruit_num[N];  // response variable
	vector[N] height;  // height

	// data for group-level effects of site
	int<lower=1> site[N];
	int<lower=1> N_sites;

}
transformed data {
	vector[N] height_cen;
	height_cen = height - mean(height);
}
parameters {
	real h_ef;  // height effects
	vector[N_sites] site_int;  // unscaled group-level intercept

}
transformed parameters {
	//linear predcitor
	vector[N] lp;
	vector[N] lambda;

	for (n in 1:N) {
	  lp[n] = site_int[site[n]] + height_cen[n] * h_ef;
	}

	// use the log link-functions
	lambda = exp(lp);
}
model{

	// likelihood contribution
	fruit_num ~ poisson(lambda);
}
generated quantities {
}"

input = make_standata(formula = fruit_num ~ height + (1|site),
  data = fruit_dat, family = "normal")

stan_dat = list(N = input$N, fruit_num = input$Y, height = input$X[, 2],
  N_sites = input$N_1, site = input$J_1)

fruit_num_glm <- stan(model_code = fruit_glm, data = stan_dat, iter = 3000, warmup = 1000,
  cores = 4, chains = 4, control = list(adapt_delta = 0.90),
  pars = c('h_ef', 'site_int', 'lambda'))

traceplot(fruit_num_glm, pars = c('h_ef', 'site_int'))
stan_ess(fruit_num_glm, pars = c('h_ef', 'site_int'))
pairs(fruit_num_glm, pars = c('h_ef', 'site_int'))

save(fruit_num_glm , file = 'fruit_num_glm.Rdata')
#load('fruit_num_glm.Rdata')

# fitted vs predicted plot
str(fruit_num_glm)
fruit_pred = apply(extract(fruit_num_glm, pars = 'lambda')$lambda, MARGIN = 2, FUN = mean)
plot(fruit_dat$fruit_num, fruit_pred)
abline(0, 1)
# residual plot
plot(fruit_dat$fruit_num, fruit_pred - fruit_dat$fruit_num)

# Under predicts the really high fruit numbers, above about 2000 fruits/plant,
# but there is only a handful of plants with more than 2000 fruits, varaiance looks to
# scale pretty well with mean
