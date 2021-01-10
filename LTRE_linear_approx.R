# LTRE to unpick the effect of vital rate variation over time and space on
# an integrative measure of reproduction R0.

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(parallel)
library(here) # used to construct a path to you local version of the scripts 

# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

source('LTER_space_time.R')
source('IPM_life_metrics.R')
################################################################################
LTRE_wrapper_ln = function(par_ob, Th = Inf, delta, rep_cent, cent_f, sur_cent, Z,  U, n0){

	print(par_ob$samp_ID)

	df = LTRE_space_time_ln(par_ob, cent_rep = rep_cent, cent_f = cent_f,
		cent_sur = sur_cent, Z = Z, U = U, n0 = n0,
		Th = Th, Delta = delta)

	df$samp_ID = par_ob$samp_ID

	return(df)
}


################################################################################
# get the data objects for the samples
gr_name = load('gr_IPM_samp.Rdata')
sur_name = load('sur_IPM_samp.Rdata')
rep_name = load('rep_IPM_samp.Rdata')
fruit_name = load('fruit_IPM_samp.Rdata')
n0_name = load('n0_pars.Rdata')

# make the size mesh note the max size observed is 80cm from a pretty big data
# set of fruit numbers, so take 85cm the max height plants can achive
Z = 0:120
U = 80
# get the intial off-spring size distribution
n0 = dlnorm(Z, n0_pars$z0_mean, n0_pars$z0_sd)

# break the samples of each parameter up into a list so can pass to mclapply
pars_list = make_par_list(rep_samps, sur_samps, gr_samps, fruit_samps)

vr_ob = make_par_ave(rep_samps, sur_samps, gr_samps, fruit_samps)
cent_rep = rep_samps$cent_const
cent_f = fruit_samps$cent_const
cent_sur = sur_samps$cent_const
Th = Inf
Delta = delta


# make a dataframe with varaition over time and space in each vital rate and R0
# start with averages over samples
rep_year = apply(rep_samps$year_int[, 2:6], MARGIN = 2, FUN = mean)
rep_slope = mean(rep_samps$slope)
rep_cent = rep_samps$cent_const
int_f = mean(fruit_samps$int)
slope_f = mean(fruit_samps$slope)
cent_f = fruit_samps$cent_const
sur_year = apply(sur_samps$year_int, MARGIN = 2, FUN = mean)
sur_slope = mean(sur_samps$slope)
sur_cent = sur_samps$cent_const
grow_year = apply(gr_samps$year_int, MARGIN = 2, FUN = mean)
grow_slope = apply(gr_samps$year_slope, MARGIN = 2, FUN = mean)
grow_sigma = mean(gr_samps$sigma)

loc_IDs = gr_samps$SPP$uni_spp
rep_loc = apply(select(rep_samps$SPP, s1:s1000), MARGIN = 1, FUN = mean)
sur_loc = apply(select(sur_samps$SPP, s1:s1000), MARGIN = 1, FUN = mean)
gr_loc = apply(select(gr_samps$SPP, s1:s1000), MARGIN = 1, FUN = mean)

# look at the Kernels
P = make_P(sur_int = sur_year[1], sur_slope = sur_slope, sur_G = sur_loc[1],
		cent_sur = sur_cent, grow_int = grow_year[1], grow_slope = grow_slope[1],
		grow_G = gr_loc[1], grow_sigma = grow_sigma, Z, U)

colnames(P) = paste0('z', Z)

df_P = as_data_frame(P) %>%
	mutate(z_t1 = Z) %>%
	gather(z_t0_lab, P_den, z0:z120) %>%
	mutate(z_t0 = as.numeric(gsub('z', '', z_t0_lab)))

ggplot(df_P, aes(x = z_t0, y = z_t1, fill = P_den)) + geom_tile() # looks as it should

Fec = make_Fec(rep_int = rep_year[1], rep_slope = rep_slope, rep_G = rep_loc[1],
	cent_rep = rep_cent, int_f = int_f, slope_f = slope_f, cent_f = cent_f, Z, U)

plot(Z, Fec)

R0_asy(P, Fec, n0)

R0_sens_ays(rep_int = rep_year[1], rep_slope = rep_slope, rep_G = rep_loc[1],
	cent_rep = rep_cent, int_f = int_f, slope_f = slope_f, cent_f = cent_f,
	sur_int = sur_year[1], sur_slope = sur_slope, sur_G = sur_loc[1], cent_sur = sur_cent,
	grow_int = grow_year[1], grow_slope = grow_slope[1], grow_G = gr_loc[1],
	grow_sigma = grow_sigma, Z, U, n0, Delta = delta, vital_rate = 'grow')

test = LTRE_space_time(pars_list[[1]], cent_rep, cent_f, cent_sur, Z, U, n0,
	Th = Inf, Delta = delta)

pdf(file = 'pred_vs_obs_R0.pdf')
	ggplot(test, aes(x = R0, y = R0_pred)) + geom_point() + geom_abline(intercept = 0, slope = 1)
dev.off()

# look at the under and over prediction
hist(test$R0_pred / test$R0)
# prortionally the biggest over predictions are small values of R0
hist(test$R0_pred)
hist(test$R0)



test2 = LTRE_wrapper_ln(pars_list[[1]], Th = Inf, delta = delta,
	rep_cent = rep_cent, cent_f = cent_f, sur_cent = sur_cent, Z = Z,  U = U,
	n0 = n0)

#################################################################################
# make the LTRE for each sampled parameter set
LTRE_full_ln_list = mclapply(pars_list, FUN = LTRE_wrapper_ln, Th = Inf, delta = delta,
	rep_cent = rep_cent, cent_f = cent_f, sur_cent = sur_cent, Z = Z,  U = U,
	n0 = n0, mc.cores = 3)

LTRE_df_ln = bind_rows(LTRE_full_ln_list)

save(LTRE_df_ln, file = 'LTRE_df_ln_full_samp.Rdata')

################################################################################
# make a dataframe that is R0 over space, in the average year

test = R0_year_ave(pars_list[[3]], cent_rep = rep_cent, cent_f = cent_f, cent_sur = sur_cent,
	Z = Z, U = U, n0 = n0, Th = Inf)

R0_time_ave = mclapply(pars_list, FUN = R0_year_ave, cent_rep = rep_cent,
	cent_f = cent_f, cent_sur = sur_cent, Z = Z, U = U, n0 = n0, Th = Inf,
	mc.cores = 3)

R0_time_ave_df = bind_rows(R0_time_ave)

save(R0_time_ave_df, file = 'R0_time_ave_samp.Rdata')

#########################################################################################
# do a rondomisation to show how covariance in vital rates affects estimates of R0
R0_shuffel = function(par_set, cent_rep = rep_cent, cent_f = cent_f, cent_sur = sur_cent,
	Z = Z, U = U, n0 = n0, Th = Inf, num_rands = 3, rand_suff = ''){

	# observed R0 for that sample
	R0_samp = list()
	R0_samp[[1]] = R0_year_ave(par_set, cent_rep = rep_cent, cent_f = cent_f, cent_sur = sur_cent,
		Z = Z, U = U, n0 = n0, Th = Inf) %>%
		mutate(rand_ID = 'obs')

	# shuffel location effects 
	for(i in 2:(num_rands + 1)){

		par_set$sur_G = base::sample(par_set$sur_G, size = length(par_set$sur_G))
		par_set$rep_G = base::sample(par_set$rep_G, size = length(par_set$rep_G))
		par_set$grow_G = base::sample(par_set$grow_G, size = length(par_set$grow_G))

		R0_samp[[i]] = R0_year_ave(par_set, cent_rep = rep_cent, cent_f = cent_f, 
			cent_sur = sur_cent, Z = Z, U = U, n0 = n0, Th = Inf) %>%
			mutate(rand_ID = paste0('r', (i - 1), rand_suff))
	}

	return(bind_rows(R0_samp))
	
}

test = R0_shuffel(pars_list[[1]], cent_rep = rep_cent, cent_f = cent_f, 
	cent_sur = sur_cent, Z = Z, U = U, n0 = n0, Th = Inf, num_rands = 10, 
	rand_suff = '')

R0_null_test = mclapply(pars_list, FUN = R0_shuffel, cent_rep = rep_cent,
	cent_f = cent_f, cent_sur = sur_cent, Z = Z, U = U, n0 = n0, Th = Inf,
	num_rands = 100, rand_suff = 'b', mc.cores = 6)

R0_null_df = bind_rows(R0_null_test)

write.csv(R0_null_df, file = 'R0_null_101_to_200_df.csv')

#########################################################################################
# Now I need a metric to say how different each randomised distribution 
# over space is from the observed distribution. Look at shifts in quantiles 
# as a metric becuase this can be used to see both shifts and narrowing, even  
# in the tails

# read in the null distibutions 
R0_null_df = read.csv('R0_null_1_to_10_df.csv', header = TRUE, stringsAsFactors = FALSE) 

# get the shifts for the different quantiles 
R0_sum_df = group_by(R0_null_df, samp_ID, rand_ID) %>%
	summarise( 
		R0_001 = quantile(R0, probs = 0.01),
		R0_025 = quantile(R0, probs = 0.025),
		R0_05 = quantile(R0, probs = 0.05),
		R0_50 = quantile(R0, probs = 0.5),
		R0_95 = quantile(R0, probs = 0.95),
		R0_975 = quantile(R0, probs = 0.975),
		R0_99 = quantile(R0, probs = 0.99),
		R0_mean = mean(R0))

# split the data frame into observed and null dist
R0_obs_df = filter(R0_sum_df, rand_ID == 'obs') %>%
	select(-rand_ID)
R0_rand_df = filter(R0_sum_df, rand_ID != 'obs')

# join them so the obs are repeated for each null realisation
# then get the relative shifts, null / obs means left shift of null are <1 
shift_df = left_join(R0_rand_df, R0_obs_df, by = 'samp_ID') %>%
	mutate(
		shift_001 = R0_001.x / R0_001.y,
		shift_025 = R0_025.x / R0_025.y,
		shift_05 = R0_05.x / R0_05.y,
		shift_50 = R0_50.x / R0_50.y,
		shift_95 = R0_95.x / R0_95.y,
		shift_975 = R0_975.x / R0_975.y,
		shift_99 = R0_99.x / R0_99.y,
		shift_mean = R0_mean.x / R0_mean.y) %>%
	select(samp_ID, rand_ID, matches("^shift_"))

# plot the histograms for 4 of the measures
rand_plt_df1 = select(ungroup(shift_df), shift_001, shift_05, shift_50, shift_mean, shift_95, shift_99) %>%
	gather(key = 'quant', value = 'R0_shift', shift_001:shift_99) %>%
	mutate(quant_lab = factor(recode(quant, 
		shift_001 = '1% quantile',
		shift_05 = '5% quantile',
		shift_50 = 'median',
		shift_mean = 'mean',
		shift_95 = '95% quantile',
		shift_99 = '99% quantile'), 
		levels = c('1% quantile', 'median', '99% quantile', '5% quantile', 'mean', '95% quantile')))
gc()

write.csv(rand_plt_df1, file = 'rand_plt_df_1_to_10.csv')

# second lot of randomisations 
R0_null_df = NULL
gc()
R0_null_df = read.csv('R0_null_11_to_100_df.csv', header = TRUE, stringsAsFactors = TRUE) 
gc()

# get the shifts for the different quantiles 
R0_sum_df = group_by(R0_null_df, samp_ID, rand_ID) %>%
	summarise( 
		R0_001 = quantile(R0, probs = 0.01),
		R0_025 = quantile(R0, probs = 0.025),
		R0_05 = quantile(R0, probs = 0.05),
		R0_50 = quantile(R0, probs = 0.5),
		R0_95 = quantile(R0, probs = 0.95),
		R0_975 = quantile(R0, probs = 0.975),
		R0_99 = quantile(R0, probs = 0.99),
		R0_mean = mean(R0))
gc()

# split the data frame into observed and null dist
R0_obs_df = filter(R0_sum_df, rand_ID == 'obs') %>%
	select(-rand_ID)
R0_rand_df = filter(R0_sum_df, rand_ID != 'obs')
gc()

# join them so the obs are repeated for each null realisation
# then get the relative shifts, null / obs means left shift of null are <1 
shift_df = left_join(R0_rand_df, R0_obs_df, by = 'samp_ID') %>%
	mutate(
		shift_001 = R0_001.x / R0_001.y,
		shift_025 = R0_025.x / R0_025.y,
		shift_05 = R0_05.x / R0_05.y,
		shift_50 = R0_50.x / R0_50.y,
		shift_95 = R0_95.x / R0_95.y,
		shift_975 = R0_975.x / R0_975.y,
		shift_99 = R0_99.x / R0_99.y,
		shift_mean = R0_mean.x / R0_mean.y) %>%
	select(samp_ID, rand_ID, matches("^shift_"))
gc()

# plot the histograms for 4 of the measures
rand_plt_df2 = select(ungroup(shift_df), shift_001, shift_05, shift_50, shift_mean, shift_95, shift_99) %>%
	gather(key = 'quant', value = 'R0_shift', shift_001:shift_99) %>%
	mutate(quant_lab = factor(recode(quant, 
		shift_001 = '1% quantile',
		shift_05 = '5% quantile',
		shift_50 = 'median',
		shift_mean = 'mean',
		shift_95 = '95% quantile',
		shift_99 = '99% quantile'), 
		levels = c('1% quantile', 'median', '99% quantile', '5% quantile', 'mean', '95% quantile')))
gc()

write.csv(rand_plt_df2, file = 'rand_plt_df_11_to_100.csv')

# third lot of randomisations 
R0_null_df = NULL 
R0_sum_df = NULL
gc()
R0_null_df = read.csv('R0_null_101_to_200_df.csv', header = TRUE, stringsAsFactors = TRUE) 
gc()

# get the shifts for the different quantiles 
R0_sum_df = group_by(R0_null_df, samp_ID, rand_ID) %>%
	summarise( 
		R0_001 = quantile(R0, probs = 0.01),
		R0_025 = quantile(R0, probs = 0.025),
		R0_05 = quantile(R0, probs = 0.05),
		R0_50 = quantile(R0, probs = 0.5),
		R0_95 = quantile(R0, probs = 0.95),
		R0_975 = quantile(R0, probs = 0.975),
		R0_99 = quantile(R0, probs = 0.99),
		R0_mean = mean(R0))
gc()

# split the data frame into observed and null dist
R0_obs_df = filter(R0_sum_df, rand_ID == 'obs') %>%
	select(-rand_ID)
R0_rand_df = filter(R0_sum_df, rand_ID != 'obs')
gc()

# join them so the obs are repeated for each null realisation
# then get the relative shifts, null / obs means left shift of null are <1 
shift_df = left_join(R0_rand_df, R0_obs_df, by = 'samp_ID') %>%
	mutate(
		shift_001 = R0_001.x / R0_001.y,
		shift_025 = R0_025.x / R0_025.y,
		shift_05 = R0_05.x / R0_05.y,
		shift_50 = R0_50.x / R0_50.y,
		shift_95 = R0_95.x / R0_95.y,
		shift_975 = R0_975.x / R0_975.y,
		shift_99 = R0_99.x / R0_99.y,
		shift_mean = R0_mean.x / R0_mean.y) %>%
	select(samp_ID, rand_ID, matches("^shift_"))
gc()

# plot the histograms for 4 of the measures
rand_plt_df3 = select(ungroup(shift_df), shift_001, shift_05, shift_50, shift_mean, shift_95, shift_99) %>%
	gather(key = 'quant', value = 'R0_shift', shift_001:shift_99) %>%
	mutate(quant_lab = factor(recode(quant, 
		shift_001 = '1% quantile',
		shift_05 = '5% quantile',
		shift_50 = 'median',
		shift_mean = 'mean',
		shift_95 = '95% quantile',
		shift_99 = '99% quantile'), 
		levels = c('1% quantile', 'median', '99% quantile', '5% quantile', 'mean', '95% quantile')))
gc()

write.csv(rand_plt_df3, file = 'rand_plt_df_101_to_200.csv')

#########################################################################################
# take the saved summaries and bind together, making these is very memory intesive 
rand_plt_df = bind_rows(read.csv('rand_plt_df_1_to_10.csv', header = TRUE),
	read.csv('rand_plt_df_11_to_100.csv', header = TRUE),
	read.csv('rand_plt_df_101_to_200.csv', header = TRUE)) %>%
mutate(quant_lab = factor(recode(quant, 
	shift_001 = '1% quantile',
	shift_05 = '5% quantile',
	shift_50 = 'median',
	shift_mean = 'mean',
	shift_95 = '95% quantile',
	shift_99 = '99% quantile'), 
		levels = c('1% quantile', 'median', '99% quantile', '5% quantile', 'mean', '95% quantile')))

rand_R0_plt = ggplot(rand_plt_df, aes(x = R0_shift)) +
	geom_histogram(fill = grey(0.5), bins = 100) +
	labs(x = expression(R[0]~at~quantile~of~null~dist*" / "*R[0]~at~quantile~of~obs~dist)) +
	theme_grey() +
	theme(axis.title.y = element_blank(),
		axis.text.y = element_blank(),
		axis.text.x = element_text(size = 14),
		axis.title.x = element_text(size = 18),
		strip.text = element_text(size = 18),
		panel.grid.minor = element_blank()) +
	annotate('segment', x = 1, xend = 1, y = -Inf, yend = Inf) +
	facet_wrap(~ quant_lab, ncol = 3, scales = 'free_x') 

pdf(file = 'rand_test_dem_comp.pdf', width = 12, height  = 8)
	rand_R0_plt
dev.off()
#