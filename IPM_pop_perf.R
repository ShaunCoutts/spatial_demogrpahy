# Calculating and mapping R0 and sensitivity of R0
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(here) # used to construct a path to you local version of the scripts 

# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

################################################################################
# get the data objects for the samples
gr_name = load('gr_IPM_samp.Rdata')
sur_name = load('sur_IPM_samp.Rdata')
rep_name = load('rep_IPM_samp.Rdata')
fruit_name = load('fruit_IPM_samp.Rdata')
n0_name = load('n0_pars.Rdata')

# make the size mesh note the max size observed is 80cm from a pretty big data
# set of fruit numbers, so take 85cm the max height plants can achive
Z = 0:85

# get the intial off-spring size distribution
n0 = dlnorm(Z, n0_pars$z0_mean, n0_pars$z0_sd)

# find R0 for each location, set up a dataframe to hold these results
# age at which change in mean reproduction is assesed
rep_age = 1
R0_summary = R0_summary_df(rep_samps, sur_samps, gr_samps, rep_age, Z, n0)


#This takes a long time to produce so save the dataframe
write.csv(R0_summary, file = 'R0_summary.csv')

################################################################################
# reload in the results and start from here
R0_summary = read.csv('R0_summary.csv', header = TRUE, stringsAsFactors = FALSE)

# Map these points
R0_mean_map = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = R0_mean), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('a) '*R[0]*''), colour = '') +
	scale_colour_gradient2(space = "Lab") +
	annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
	annotate('text', x = 25, y =-140, label = '50 m')

R0_mean_map

# make a four panel mapping the R0 and the sensivities
rep_elas_map = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = elas_rep), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('d) elasticity of '*R[0]*' to rep: '*delta*'ln('*R[0]*')/'*delta*'ln('*beta[0]^r*')'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'darkorange', space = "Lab")

sur_elas_map = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = elas_sur), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('c) elasticity of '*R[0]*' to survival: '*delta*'ln('*R[0]*')/'*delta*'ln('*beta[0]^s*')'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'darkorange', space = "Lab")

gr_elas_map = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = elas_gr), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('b) elasticity of '*R[0]*' to growth: '*delta*'ln('*R[0]*')/'*delta*'ln('*beta[0]^g*')'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'darkorange', space = "Lab")

pdf(file = 'R0_and_elasticity.pdf', height = 8, width = 8)
	plot_grid(R0_mean_map, gr_elas_map, sur_elas_map, rep_elas_map, ncol = 2)
dev.off()

# Now the same map but lok at the elasticity of R0 to life expectancy and mean rep in firsrt year
sur_elas_le = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = elas_sur_le), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('c) elasticity of '*R[0]*' to sur:'*delta*'ln('*R[0]*')/'*delta*'ln(life ex)'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'darkorange', space = "Lab")

rep_elas_fl = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = elas_rep_fl), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('d) elasticity of '*R[0]*' to rep:'*delta*'ln('*R[0]*')/'*delta*'ln(mean rep)'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'darkorange', space = "Lab")

life_ex_map = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = life_expect), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('e) Life expectancy'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'magenta', space = "Lab")

mean_rep_map = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = mean_repA), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('f) mean repoduction first year'), colour = '') +
	scale_colour_gradient(low = 'white', high = 'magenta', space = "Lab")

sen_life_ex = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = sen_life_ex), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('g) sensitivity life ex to '*beta[0]^s*''), colour = '') +
	scale_colour_gradient(low = 'white', high = 'chartreuse', space = "Lab")

sen_mean_rep = ggplot(R0_summary, aes(uX, uY)) +
	geom_point(aes(color = sen_mean_rep), size = 1) + xlim(-20, 140) + ylim(-200, 50) +
  	theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
		panel.background = element_rect(fill = grey(0.9)), axis.text.x = element_blank(),
		axis.ticks = element_blank(), axis.title.x = element_blank(),
	    axis.title.y = element_blank(),axis.line.x = element_blank(),
		axis.line.y = element_blank(), axis.text.y = element_blank(),
		legend.position = c(0.95, 0.6), legend.justification = c(1,0),
	    plot.title = element_text(hjust = 0)) +
	labs(title = bquote('g) sensitivity mean rep to '*beta[0]^r*''), colour = '') +
	scale_colour_gradient(low = 'white', high = 'chartreuse', space = "Lab")


pdf(file = 'R0_and_elasticity_inter.pdf', height = 16, width = 8)
	plot_grid(R0_mean_map, gr_elas_map, sur_elas_le, rep_elas_fl,
		life_ex_map, mean_rep_map, sen_life_ex, sen_mean_rep, ncol = 2)
dev.off()

##################################################################################
# Explore the negative co-variance between R0 and elacticity a bit more, first just
# plot it with PCA
PCA_elas_dat = select(R0_summary, R0_mean, elas_rep:elas_gr)
plot(log(PCA_elas_dat))

life_met_PCA = prcomp(log(PCA_elas_dat), center = TRUE, scale. = TRUE)

summary(life_met_PCA)
biplot(life_met_PCA)
# Looks like 2 axis explain most of the variance (88%), and survival elast and R0 strong
# negative co-variance, R0 pretty orthoganal to the other two

PCA_int_dat =  select(R0_summary, R0_mean, elas_gr, elas_rep_fl, elas_sur_le)
plot(log(PCA_int_dat))
# pattern much stronger here
life_met_PCA = prcomp(log(PCA_int_dat), center = TRUE, scale. = TRUE)

summary(life_met_PCA)
biplot(life_met_PCA)
# here 2 axis explain 97% of variance, R0 strongly negativly covaries with elast_sr and elast_rep
# elast growth orthoganal.

# Also test if the covariance is jsut the result of constraint, if it is then within the
# distibutions of estimates for each plant (i.e. the stan samples) we should see correlation
rep_age = 1
con_test = covar_test(rep_samps, sur_samps, gr_samps, rep_age, Z, n0)

#write.csv(con_test, file = 'constraint_test.csv')
con_test = read.csv('constraint_test.csv', header = TRUE, stringsAsFactors = FALSE)

hist_dat = select(con_test, C_R0_sur:C_R0_gr) %>%
	gather(key = 'comp', value = 'corr', C_R0_sur:C_R0_gr) %>%
	mutate(comp_pretty = recode(comp, C_R0_sur = 'R0 vs elas R0 to sur',
		C_R0_sur_le = 'R0 vs elas R0 to life ex',
		C_R0_rep = 'R0 vs elas R0 to rep',
		C_R0_rep_fl = 'R0 vs elas R0 to mean rep',
		C_R0_gr = 'R0 vs elas R0 to growth'))

hist_plt = ggplot(hist_dat, aes(x = corr)) + geom_histogram() +
	xlim(-1, 1) + theme_minimal() +
	facet_wrap(~comp_pretty, ncol = 3)

hist_plt

plt_dat = select(con_test, R0_mean, C_R0_sur:C_R0_gr) %>%
	gather(key = 'comp', value = 'corr', C_R0_sur:C_R0_gr) %>%
	mutate(comp_pretty = recode(comp, C_R0_sur = 'R0 vs elas R0 to sur',
		C_R0_sur_le = 'R0 vs elas R0 to life ex',
		C_R0_rep = 'R0 vs elas R0 to rep',
		C_R0_rep_fl = 'R0 vs elas R0 to mean rep',
		C_R0_gr = 'R0 vs elas R0 to growth'))

cor_plt = ggplot(plt_dat, aes(x = R0_mean, y = corr)) +
	geom_point() +
	labs(x = 'R0 at each location', y = 'corrolation between R0 and elasticity within each location') +
	scale_x_log10() +
	theme_grey() +
	facet_wrap(~comp_pretty, ncol = 3)

pdf(file = 'constraint_test.pdf', width = 8, height = 5)
	cor_plt
dev.off()

#########################################################################################
# Null test of effect of co-variance of vital rates on population growth

# do this in parallel and save the results as it looks like it will take a while to run
pars_list = list()
pars_list[[1]] =list(sur_samp = sur_samps, rep_samp = rep_samps, gr_samp = gr_samps,
	age = rep_age, Z = Z, n0 = n0, rand_num = 'obs') 
# Do 500 randomisations 


null_R0 = null_R0_df(pars_list[[1]])

