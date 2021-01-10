# plotting the LTRE results

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(here) # used to construct a path to you local version of the scripts 

# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

###################################################################################
# Look at the plots in ln(R0) that improves the accuracy of the linear approximation
LTRE_ob_name = load('LTRE_df_ln_full_samp.Rdata')

LTRE_hist_ln_df = select(LTRE_df_ln, loc_ID, year, cont_T, cont_rep_T, cont_sur_T, cont_gr_T, cont_gr_s_T,
	 	cont_L, cont_rep_L, cont_sur_L, cont_gr_L, cont_LT, cont_rep_LT, cont_sur_LT, cont_gr_LT, cont_gr_s_LT,
		samp_ID) %>%
	rename(
			cont_grs_T = cont_gr_s_T,
			cont_grs_LT = cont_gr_s_LT
		) %>%
	gather(key = 'vital_rate_TL', value = 'cont', cont_T:cont_grs_LT) %>%
	mutate(TL = ifelse(grepl('_T', vital_rate_TL), 'time',
			ifelse(grepl('_LT', vital_rate_TL), 'location:time', 'location')),
		vital_rate = ifelse(grepl('rep', vital_rate_TL), 'Prob. flower',
			ifelse(grepl('sur', vital_rate_TL), 'Survival',
			ifelse(grepl('gr_', vital_rate_TL), 'Growth',
			ifelse(grepl('grs', vital_rate_TL), 'Growth Slope', 'Full'))))) %>%
	mutate(vital_rate = factor(vital_rate, levels = c('Full', 'Survival', 'Growth', 
		'Growth Slope', 'Prob. flower')))

LTRE_means_ln_L = filter(LTRE_hist_ln_df, TL == 'location') %>%
	group_by(loc_ID, year, vital_rate, TL) %>%
	summarise(cont = mean(cont))

LTRE_means_ln_T = filter(LTRE_hist_ln_df, TL == 'time') %>%
	group_by(loc_ID, year, vital_rate, TL) %>%
	summarise(cont = mean(cont))

cont_hist_ln = ggplot(LTRE_hist_ln_df, aes(x = cont)) +
	geom_histogram(bins = 25, aes(y = ..density..)) +
	geom_rug(data = LTRE_means_ln_L, sides = 'b', alpha = 0.2, colour = 'goldenrod') +
	geom_rug(data = LTRE_means_ln_T, sides = 'b', alpha = 0.2, colour = 'goldenrod') +
	labs(x = 'contribution to ln(R0)', y = 'density') +
	xlim(-6, 6) +
	theme_grey() +
	theme(panel.grid.minor = element_blank()) +
	facet_grid(vital_rate~TL)

pdf(file = 'cont_R0_ln_hist.pdf', width = 8, height = 12)
	cont_hist_ln
dev.off()

# split the data into time and location as better ploted in different ways
# Just location
LTRE_hist_ln_L = filter(LTRE_hist_ln_df, TL == 'location')

# make a box plot for the time contirbutions
LTRE_dat_ln_T = filter(LTRE_hist_ln_df, TL == 'time') %>%
	group_by(year, TL, vital_rate, samp_ID) %>%
	summarise(cont = mean(cont))

# box plots for interations, need mean for each year:location combination to take out estimation uncertianty 
LTRE_dat_ln_LT = filter(LTRE_hist_ln_df, TL == 'location:time')
# means over estimation uncertianty to show eccect of location seperate to sample variation
LTRE_means_ln_LT = group_by(LTRE_dat_ln_LT, loc_ID, year, vital_rate) %>%
	summarise(cont = mean(cont)) 
	
	%>%
	group_by(year, vital_rate) %>%
	summarise(
		med_cont = median(cont),
		min_cont = min(cont),
		max_cont = max(cont)
	)
# a few samples make the tails of violin plots very long trim each location in each year to the 95% credible interval
# find interval for each group
LTRE_dat_ln_LT_CI = group_by(LTRE_dat_ln_LT, loc_ID, year, vital_rate) %>%
	summarise(
		l95_CI = quantile(cont, probs = 0.025),
		u95_CI = quantile(cont, probs = 0.975)
	)

LTRE_dat_ln_LT_trim = left_join(LTRE_dat_ln_LT, LTRE_dat_ln_LT_CI, 
		by = c('loc_ID', 'year', 'vital_rate')) %>%  
	filter(cont >= l95_CI & cont <= u95_CI)

# Location plot constant over years variation from location (rug plot) and sample variaion (grey bars)
cont_hist_ln_L = ggplot(LTRE_hist_ln_L, aes(x = cont)) +
	geom_histogram(bins = 25, aes(y = ..density..)) +
	geom_rug(data = LTRE_means_ln_L, sides = 'b', alpha = 0.2, colour = 'goldenrod') +
	labs(x = bquote('contribution to ln('*E[F]*')')) +
	xlim(-6, 6) +
	theme_grey() +
	theme(panel.grid.minor = element_blank(),
		axis.title.y = element_blank(),
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.title.x = element_text(size = 18),
		axis.text.x = element_text(size = 16),
		strip.text = element_text(size = 18)) +
	facet_grid(vital_rate~TL)

pdf(file = 'cont_R0_ln_hist_L.pdf', width = 4, height = 8)
	cont_hist_ln_L
dev.off()

# time plot, constant over locations, variation from year (x-axis) and sample variation (violin plot) 
cont_viol_ln_T = ggplot(LTRE_dat_ln_T, aes(x = year, y = cont)) +
	geom_violin(fill = 'deepskyblue', colour = 'deepskyblue', alpha = 0.4) +
	labs(x = 'year', y = bquote('contribution to ln('*E[F]*')')) +
	scale_x_discrete(labels = c('2003', '2004', '2005', '2006', '2007')) +
	theme_grey() +
	theme(panel.grid.minor = element_blank(),
		axis.title = element_text(size = 18),
		axis.text = element_text(size = 16),
		strip.text = element_text(size = 18)) +
	facet_grid(vital_rate~TL)

# location:time interation, variation from location (box-plot), year (x-axis) and sample variation (violin plot) 
cont_viol_ln_LT = ggplot(LTRE_dat_ln_LT_trim, aes(x = year, y = cont)) +
	geom_violin(fill = 'deepskyblue', colour = 'deepskyblue', alpha = 0.4) +
	geom_point(data = LTRE_means_ln_LT, size = 0.5, colour = 'goldenrod', alpha = 0.5) +
	labs(x = 'year', y = bquote('contribution to ln('*E[F]*')')) +
	scale_x_discrete(labels = c('2003', '2004', '2005', '2006', '2007')) +
	theme_grey() +
	theme(panel.grid.minor = element_blank(),
		axis.title = element_text(size = 18),
		axis.text = element_text(size = 16),
		strip.text = element_text(size = 18)) +
	facet_grid(vital_rate~TL)

pdf(file = 'LTRE_space_time_cont_lnR0.pdf', width = 12, height = 9)
	plot_grid(cont_hist_ln_L, cont_viol_ln_T, cont_viol_ln_LT, ncol = 3, labels = c('a)', 'b)', 'c)'), 
		rel_widths = c(1, 1.2, 1.2))
dev.off()

################################################################################
# make a map of R0, using mean of each location, averaged across years
################################################################################
# load the dataframe with the R0 estimated for an average year
R0_ob_name = load('R0_time_ave_samp.Rdata')

# get the unified locations
gr_name = load('gr_IPM_samp.Rdata')
loc_xy = select(gr_samps$SPP, unify_Y, unify_X, uni_spp) %>%
	rename(loc_ID = uni_spp, X = unify_X, Y = unify_Y)

R0_df = group_by(R0_time_ave_df, loc_ID) %>%
	summarise(R0_med = median(R0), R0_men = mean(R0), R0_sd = sd(R0),
		R0_IQ = (quantile(R0, probs = 0.75) - quantile(R0, probs = 0.25)),
		ln_R0_med = median(log(R0)), ln_R0_men = mean(log(R0))) %>%
	left_join(loc_xy, by = 'loc_ID')

R0_med_map = ggplot(R0_df, aes(X, Y)) +
	geom_point(aes(color = R0_med), size = 1) + xlim(-20, 140) + ylim(-180, 40) +
	scale_color_viridis() +
  	theme(panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = grey(0.95)),
	    axis.text = element_blank(), axis.ticks = element_blank(),
		axis.title = element_blank(), axis.line = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    legend.key.size = unit(10, 'mm'), legend.text = element_text(size = 15),
	    legend.title = element_text(size = 20)) +
	labs(colour = bquote('median '*E[F]*'')) +
	annotate('segment', x = 80, xend = 130, y = 30, yend = 30, color = grey(0)) +
	annotate('text', x = 105, y = 40, label = '50 m', size  = 7) +
	annotate("rect", xmin = 10.5, xmax = 43, ymin = -58, ymax = 13, colour = grey(0), alpha = 0.0) +
	annotate("text", x = 10, y = 20, label = 'panel c', hjust = 0, size = 6) +
	coord_fixed(ratio = 1.375)

R0_med_hist = ggplot(R0_df, aes(x = R0_med / 1000)) +
	geom_histogram(bins = 25) +
	labs(x = bquote('median '*E[F]*' (x1000)')) +
	theme(axis.text.y = element_blank(), axis.title.x = element_text(size = 20),
		axis.ticks.y = element_blank(), axis.title.y = element_blank(),
		axis.line = element_blank(), plot.background = element_rect(fill = grey(1)))

# add zoomed in map 
EF_med_map_zoomed = ggplot(R0_df, aes(X, Y)) +
	geom_point(aes(color = R0_med), size = 1) + xlim(13.5,  40) + ylim(-50, 5.2) +
	scale_color_viridis() +
  	theme(panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = grey(0.95)),
	    axis.text = element_blank(), axis.ticks = element_blank(),
		axis.title = element_blank(), axis.line = element_blank(),
		legend.position = 'none') +
	annotate('segment', x = 30, xend = 35, y = 5, yend = 5, color = grey(0)) +
	annotate('text', x = 32.5, y = 2, label = '5 m', size = 8) + theme(legend.position = 'none') +
	coord_fixed(ratio = 2.075)

EF_map_hist = ggdraw() +
		draw_plot(R0_med_map, 0, 0, 1, 1) +
		draw_plot(R0_med_hist, 0.05, 0.02, 0.5, 0.35) +
		draw_plot_label(c('a)', 'b)'), c(0.02, 0.02), c(1, 0.38), size = 20)
dev.off()


pdf(file = 'R0_med_map_true_asp.pdf', width = 8, height = 10)
	  plot_grid(EF_map_hist, EF_med_map_zoomed, ncol = 2, rel_widths = c(1, 0.45),
	  	label_size = 20, labels  = c('', 'c)'))
dev.off()


##############################################################################################
# quick randomisation to test if the median absolute deviation of the distribtuions 
# of contributions are darwn from the same distribution 
LTRE_mad_random = function(vect1, vect2, n_samp = 100){

	#observed range
	mad1 = mad(vect1)
	mad2 = mad(vect2)

	# split into two randomized groups 
	joint_vect = c(vect1, vect2)

	rand_mad_dif = numeric(n_samp) 

	for(n in 1:n_samp){
	
		# split into two vectors with replacement 
		rand1 = sample(joint_vect, size = length(vect1), replace = TRUE) 
		rand2 = sample(joint_vect, size = length(vect2), replace = TRUE) 

		# get median absolute deviation of each group then get difference
		rand_mad1 = mad(rand1)

		rand_mad2 = mad(rand2)

		rand_mad_dif[n] = abs(rand_mad1 - rand_mad2)

	}
	
	obs_diff = abs(mad1 - mad2)

	pvalue_lower_tail = sum(rand_mad_dif >= obs_diff) / n_samp
	pvalue_upper_tail = sum(rand_mad_dif <= obs_diff) / n_samp
	
	out_df = list(mad1 = mad1, mad2 = mad2, obs_diff = obs_diff, rand_mad_diff = rand_mad_dif, 
		p_value = ifelse(pvalue_upper_tail < pvalue_lower_tail, pvalue_upper_tail, pvalue_lower_tail))

	return(out_df)
}

# get the sample data and get the means for each LTRE component 

LTRE_ob_name = load('LTRE_df_ln_full_samp.Rdata')

LTRE_hist_ln_df = select(LTRE_df_ln, loc_ID, year, cont_T, cont_rep_T, cont_sur_T, cont_gr_T, cont_gr_s_T,
	 	cont_L, cont_rep_L, cont_sur_L, cont_gr_L, cont_LT, cont_rep_LT, cont_sur_LT, cont_gr_LT, cont_gr_s_LT,
		samp_ID) %>%
	rename(
			cont_grs_T = cont_gr_s_T,
			cont_grs_LT = cont_gr_s_LT
		) %>%
	gather(key = 'vital_rate_TL', value = 'cont', cont_T:cont_grs_LT) %>%
	mutate(TL = ifelse(grepl('_T', vital_rate_TL), 'time',
			ifelse(grepl('_LT', vital_rate_TL), 'location:time', 'location')),
		vital_rate = ifelse(grepl('rep', vital_rate_TL), 'Prob. flower',
			ifelse(grepl('sur', vital_rate_TL), 'Survival',
			ifelse(grepl('gr_', vital_rate_TL), 'Growth',
			ifelse(grepl('grs', vital_rate_TL), 'Growth Slope', 'Full'))))) %>%
	mutate(vital_rate = factor(vital_rate, levels = c('Full', 'Survival', 'Growth', 
		'Growth Slope', 'Prob. flower')))


# get the means for each component for each vital rate 
LTRE_var_rand_df = group_by(LTRE_hist_ln_df, loc_ID, year, TL, vital_rate) %>%
	summarise(cont = mean(cont))

loc_cont_full = unique(filter(LTRE_var_rand_df, TL == 'location', vital_rate == 'Full')$cont)
loc_cont_sur = unique(filter(LTRE_var_rand_df, TL == 'location', vital_rate == 'Survival')$cont)
loc_cont_gr = unique(filter(LTRE_var_rand_df, TL == 'location', vital_rate == 'Growth')$cont)
loc_cont_rep = unique(filter(LTRE_var_rand_df, TL == 'location', vital_rate == 'Prob. flower')$cont)

time_cont_full = unique(filter(LTRE_var_rand_df, TL == 'time', vital_rate == 'Full')$cont)
time_cont_sur = unique(filter(LTRE_var_rand_df, TL == 'time', vital_rate == 'Survival')$cont)
time_cont_gr = unique(filter(LTRE_var_rand_df, TL == 'time', vital_rate == 'Growth')$cont)
time_cont_gr_slp = unique(filter(LTRE_var_rand_df, TL == 'time', vital_rate == 'Growth Slope')$cont)
time_cont_rep = unique(filter(LTRE_var_rand_df, TL == 'time', vital_rate == 'Prob. flower')$cont)

TL_cont_full = unique(filter(LTRE_var_rand_df, TL == 'location:time', vital_rate == 'Full')$cont)
TL_cont_sur = unique(filter(LTRE_var_rand_df, TL == 'location:time', vital_rate == 'Survival')$cont)
TL_cont_gr = unique(filter(LTRE_var_rand_df, TL == 'location:time', vital_rate == 'Growth')$cont)
TL_cont_gr_slp = unique(filter(LTRE_var_rand_df, TL == 'location:time', vital_rate == 'Growth Slope')$cont)
TL_cont_rep = unique(filter(LTRE_var_rand_df, TL == 'location:time', vital_rate == 'Prob. flower')$cont)

# obtain comparisions between each difference for the full contribution of each component to reduce number of comparisons 
# 1. location - time
# 2. location - interaction
# 3. time - interaction

loc_vs_time = LTRE_mad_random(loc_cont_full, time_cont_full, n_samp = 10000)
loc_vs_TL = LTRE_mad_random(loc_cont_full, TL_cont_full, n_samp = 10000)
time_vs_TL = LTRE_mad_random(time_cont_full, TL_cont_full, n_samp = 10000)


rand_df = data.frame(
	TL1 = c('location', 'location', 'time'),
	TL1_mad = c(mad(loc_cont_full), mad(loc_cont_full), mad(time_cont_full)),
	TL2 = c('time', 'location:time', 'location:time'),
	TL2_mad = c(mad(time_cont_full), mad(TL_cont_full), mad(TL_cont_full)),
	obs_diff= c(loc_vs_time$obs_diff, loc_vs_TL$obs_diff, time_vs_TL$obs_diff),
	p_value = c(loc_vs_time$p_value, loc_vs_TL$p_value, time_vs_TL$p_value) 
)


