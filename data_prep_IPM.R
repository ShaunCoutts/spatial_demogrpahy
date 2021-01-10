# make the cohort fecundity and senesitivity of that fecundity to pertibation in
# intercept of vital rate functions
library(dplyr)
library(tidyr)
library(rstan)
library(here) # used to construct a path to you local version of the scripts 

# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

source('dist_neigh_setup.R')
source('model_perf_output_helper.R')

################################################################################
# data to get the centering constants and to match up the knot locations so spatial
# effect is propely alinged between vital rates
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
vr_loc_gen = filter(vr_loc_gen, !is.na(vr_loc_gen$X))

# for each vital rate I tool 8000 stan samples to get a decent number of effective
# samples, thin each vital rate samples so it is not so cumbersome to work with to get 1000
# independent samples
n.thin = 1000

# for each vital rate get the data used to fit the model, centering constant if used
# X and Y coords of each data location, knot index, SPP samples for that knot, coordinate
# for that knot, year averaged intercepts for each vital rate and slopes on height

# REPRODUCTION
rep_dat = select(vr_loc_gen, uID, uLoc, rep, year, height, X, Y) %>%
	filter(!is.na(height), rep <= 1)

# find mean height
rep_mean_height = mean(rep_dat$height)

# get the knot coordinates and knot index for each location to reference the
# stan samples
knot_rep = knot_coords2(dat_x = rep_dat$X, dat_y = rep_dat$Y, min_dist = 0.5)
rep_dat = mutate(rep_dat, knot_ind = apply(distance_matrix2(x1 = X, y1 = Y,
	x2 = knot_rep[, 1], y2 = knot_rep[, 2]), MARGIN = 1,
	FUN = function(x) which(x == min(x))))

# get the stan samples
ob_name = load('SPP_NK_rep_stan_ap.Rdata')
# get the SPP
rep_spp = extract_flat(SPP_NK_rep_stan, pars = c('spp'))

# thin samples
rep_thin = sample.int(dim(rep_spp)[1], n.thin)

# connect the spp sub_samples to the knot locations
rep_SPP = cbind(data.frame(X = knot_rep[, 1], Y = knot_rep[, 2],
		spp_lab = paste0('spp[', 1:dim(knot_rep)[1], ']')),
	setNames(data.frame(t(rep_spp[rep_thin, ])), paste0('s', 1:n.thin)))

# pull the year effects, note year 1 is 2002
rep_year_effect = data.frame(extract_flat(SPP_NK_rep_stan, pars = c('year_int'))[rep_thin, ])
# take the mean of the last 4 years 2004 - 2007 to get the intercept of a "typical" year
rep_int = apply(rep_year_effect[, 3:6], MARGIN = 1, FUN = mean)

# get the slope parameter
rep_slope = extract_flat(SPP_NK_rep_stan, pars = c('h_ef'))[rep_thin, ]

# SURVIVAL #####################################################################
sur_dat = select(vr_loc_gen, uID, uLoc, year,sur, height, height_prev, X, Y) %>%
	filter(!is.na(sur))

# take out data from the first year observed since nothing can be observed dead in the first year
first_year = sapply(seq_along(sur_dat$year), FUN = function(x){

  min_year_group = min(sur_dat$year[sur_dat$uID == sur_dat$uID[x]])
  return(ifelse(sur_dat$year[x] == min_year_group, FALSE, TRUE))

})
sur_dat = sur_dat[first_year, ]

# height: this requires some assumptions because 1) individuals change size over the course of a year,
# 2) we do not have the final size for indivudals that died. 3) we have indivduals that were only found when they were very large.
# Best we can do here is assue that mortality occurs just after census, so that we can use previous height as a predictor,
# and in the case of indivduals in their first year we say they were 0
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]

sur_mean_height = mean(sur_dat$height_prev)

# get knot locations for survival data
knot_sur = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5)
sur_dat = mutate(sur_dat, knot_ind = apply(distance_matrix2(x1 = X, y1 = Y,
	x2 = knot_sur[, 1], y2 = knot_sur[, 2]), MARGIN = 1,
	FUN = function(x) which(x == min(x))))

# get the stan sampled opject
# survival
ob_name = load('SPP_NK_sur_stan.Rdata')
sur_spp = extract_flat(SPP_NK_sur_stan, pars = c('spp'))

# thin samples
sur_thin = sample.int(dim(sur_spp)[1], n.thin)

# connect the spp sub_samples to the knot locations
sur_SPP = cbind(data.frame(X = knot_sur[, 1], Y = knot_sur[, 2],
		spp_lab = paste0('spp[', 1:dim(knot_sur)[1], ']')),
	setNames(data.frame(t(sur_spp[sur_thin, ])), paste0('s', 1:n.thin)))

# pull the year effects, note year 1 is 2002
sur_year_effect = data.frame(extract_flat(SPP_NK_sur_stan, pars = c('year_int'))[sur_thin, ])
# take the mean of the last 4 years 2004 - 2007 to get the intercept of a "typical" year
sur_int = apply(sur_year_effect[, 2:5], MARGIN = 1, FUN = mean)

# get the slope parameter
sur_slope = extract_flat(SPP_NK_sur_stan, pars = c('h_ef'))[sur_thin, ]

# GROWTH #######################################################################
gr_dat = select(vr_loc_gen, year, height, height_prev, X, Y) %>%
	filter(!is.na(height), !is.na(height_prev))

# get the knot points
knot_gr = knot_coords2(dat_x = gr_dat$X, dat_y = gr_dat$Y, min_dist = 0.5)
gr_dat = mutate(gr_dat, knot_ind = apply(distance_matrix2(x1 = X, y1 = Y,
	x2 = knot_gr[, 1], y2 = knot_gr[, 2]), MARGIN = 1,
	FUN = function(x) which(x == min(x))))

# get samples from stan for model paramerters
ob_name = load('SPP_NK_gr_stan.Rdata')

gr_spp = extract_flat(SPP_gr_stan, pars = c('spp'))

gr_thin = sample.int(dim(gr_spp)[1], n.thin)

gr_SPP = cbind(data.frame(X = knot_gr[, 1], Y = knot_gr[, 2],
		spp_lab = paste0('spp[', 1:dim(knot_gr)[1], ']')),
	setNames(data.frame(t(gr_spp[gr_thin, ])), paste0('s', 1:n.thin)))

gr_year_effect = data.frame(extract_flat(SPP_gr_stan, pars = c('b0'))[gr_thin, ])
# take the mean of the last 4 years 2004 - 2007 to get the intercept of a "typical" year
gr_int = apply(gr_year_effect[, 2:5], MARGIN = 1, FUN = mean)

gr_slope_year = data.frame(extract_flat(SPP_gr_stan, pars = c('gr_rate'))[gr_thin, ])
# take the mean of the last 4 years 2004 - 2007 to get the intercept of a "typical" year
gr_slope = apply(gr_slope_year[, 2:5], MARGIN = 1, FUN = mean)
# sigma of error distribution to help make the growth kernel
gr_sigma = extract_flat(SPP_gr_stan, pars = c('sigma'))[gr_thin]

# FRUIT PRODUCTIONS ############################################################
burn_dat = read.csv('hcdem_ABS_2015_checked.csv')

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
fruit_dat = inner_join(ht_long, rep_long, by = 'join_ID') %>%
	select(burn_yr = burn_yr, site = site.x, ID = ID.x,
		join_ID = join_ID, m_year = m_year.x, height = height, fruit_num = fruit_num) %>%
	filter(fruit_num >= 2) %>%
	mutate(time_since_fire = m_year - burn_yr, site = as.character(site))

mean_height_fruit = mean(fruit_dat$height)
site_num = 2

ob_name = load('fruit_num_glm.Rdata')

fruit_int = extract_flat(fruit_num_glm, pars = c('site_int'))[, site_num]
fruit_thin = sample.int(length(fruit_int), n.thin)
fruit_int = fruit_int[fruit_thin]
fruit_slope = extract_flat(fruit_num_glm, pars = c('h_ef'))[fruit_thin]

################################################################################
# create a key that can unify the SPP dataframes so that each knot has one location
# and lots of samples for each vital rate.

# build a common set of knot locations by taking vital rate with least number of
# knots (gr_spp), then finding all the knot locatons in the other 2 that are within
# 25 cm of the knot locations in gr_spp, and add those indexes to a key

gr_sur_dist = apply(distance_matrix2(x1 = gr_SPP$X, y1 = gr_SPP$Y, x2 = sur_SPP$X, y2 = sur_SPP$Y),
	MARGIN = 1, FUN = min)

gr_rep_dist = apply(distance_matrix2(x1 = gr_SPP$X, y1 = gr_SPP$Y, x2 = rep_SPP$X, y2 = rep_SPP$Y),
	MARGIN = 1, FUN = min)

# find the locations in gr_SPP that are 25cm of a location in sur_SPP and rep_SPP
# add a set of unified X and Y to each vr_SPP that gives the co-ordinate or NA
# if it is not in all three sets
gr_SPP = mutate(gr_SPP,
	unify_X = ifelse(gr_sur_dist < 0.25 & gr_rep_dist < 0.25, X, NA),
	unify_Y = ifelse(gr_sur_dist < 0.25 & gr_rep_dist < 0.25, Y, NA)) %>%
	filter(!is.na(unify_X)) %>%
	mutate(uni_spp = paste0('uspp_', seq_along(unify_X)))

# now need to mesh these to sur_SPP and rep_SPP, do this in for loop as a bit easier to
# follow. Becuase there are more knot points in both sur_SPP and rep_SPP some of the
# unify locations are close to more than one knot, take the first, as they will
# all be similar
sur_SPP$unify_X = NA
sur_SPP$unify_Y = NA
sur_SPP$uni_spp = NA
rep_SPP$unify_X = NA
rep_SPP$unify_Y = NA
rep_SPP$uni_spp = NA

for(i in seq_along(gr_SPP$unify_X)){

	#get distance to every knot in sur_SPP and find index of minimum
	sur_dist = distance_matrix2(x1 = sur_SPP$X, y1 = sur_SPP$Y,
		x2 = gr_SPP$unify_X[i], y2 = gr_SPP$unify_Y[i])
	# check distance less than 25cm
	if(min(sur_dist) < 0.25){

		sur_ind = which(sur_dist == min(sur_dist))[1]
		# put the unfied co-ordinates at that index
		sur_SPP$unify_X[sur_ind] = gr_SPP$unify_X[i]
		sur_SPP$unify_Y[sur_ind] = gr_SPP$unify_Y[i]
		sur_SPP$uni_spp[sur_ind] = gr_SPP$uni_spp[i]
	}

	#get distance to every knot in rep_SPP and find index of minimum
	rep_dist = distance_matrix2(x1 = rep_SPP$X, y1 = rep_SPP$Y,
		x2 = gr_SPP$unify_X[i], y2 = gr_SPP$unify_Y[i])
	# check distance less than 25cm
	if(min(rep_dist) < 0.25){

		rep_ind = which(rep_dist == min(rep_dist))[1]
		# put the unfied co-ordinates at that index
		rep_SPP$unify_X[rep_ind] = gr_SPP$unify_X[i]
		rep_SPP$unify_Y[rep_ind] = gr_SPP$unify_Y[i]
		rep_SPP$uni_spp[rep_ind] = gr_SPP$uni_spp[i]
	}
}

gr_SPP[1:25, c(1:4, (dim(gr_SPP)[2] - 2):dim(gr_SPP)[2])]
sur_SPP[1:25, c(1:4, (dim(sur_SPP)[2] - 2):dim(sur_SPP)[2])]
rep_SPP[1:25, c(1:4, (dim(rep_SPP)[2] - 2):dim(rep_SPP)[2])]

################################################################################
# make some data structures to save all this stuff for the IPMing
gr_samps = list(SPP = gr_SPP, slope = gr_slope, int = gr_int, sigma = gr_sigma,
	year_int = gr_year_effect, year_slope = gr_slope_year, cent_const = NA)
save(gr_samps, file = 'gr_IPM_samp.Rdata')

sur_samps = list(SPP = filter(sur_SPP, !is.na(unify_X)), slope = sur_slope, int = sur_int, sigma = NA,
	year_int = sur_year_effect, year_slope = NA, cent_const = sur_mean_height)
save(sur_samps, file = 'sur_IPM_samp.Rdata')

rep_samps = list(SPP = filter(rep_SPP, !is.na(unify_X)), slope = rep_slope, int = rep_int, sigma = NA,
	year_int = rep_year_effect, year_slope = NA, cent_const = rep_mean_height)
save(rep_samps, file = 'rep_IPM_samp.Rdata')

fruit_samps = list(SPP = NA, slope = fruit_slope, int = fruit_int, sigma = NA,
	year_int = NA, year_slope = NA, cent_const = mean_height_fruit)
save(fruit_samps, file = 'fruit_IPM_samp.Rdata')

# note on log scale be cuase the distribution is well described by the log-normal
gr_dat_firt = filter(gr_dat, year == 2003)
n0_pars = list(z0_mean = mean(log(gr_dat_firt$height_prev)),
	z0_sd = sd(log(gr_dat_firt$height_prev)))
save(n0_pars, file = 'n0_pars.Rdata')
################################################################################
