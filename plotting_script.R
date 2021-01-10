# plotting script for individual performance over space

library(dplyr)
library(tidyr)
library(colorspace)
library(rstan)
library(arrayhelpers)
library(ggplot2)
library(cowplot)
library(viridis)
library(here) # used to construct a path to you local version of the scripts 


# build the file path form the local root, may need to modify if repo is saved in a different directory
# This will build /home/some_user_account/spatial_dem_perf
dir_path = here('spatial_dem_perf') 
setwd(dir_path)

vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
source('dist_neigh_setup.R')
source('model_perf_output_helper.R')

#need location data so drop the obervations that do not have a location
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X),]

#aggregate the data for plotting
plot_dat = group_by(vr_loc_gen, uID) %>%
	summarise(X = unique(X), Y = unique(Y), sur = min(sur),
	height = max(height, na.rm = TRUE), rep = sum(rep, na.rm = TRUE),
	num_years = unique(num_years))

# start with survival
sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]

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

plot_dat_sur = ddply(sur_dat, .(uID), summarise, X = unique(X), Y = unique(Y),
  sur = min(sur))

# plot the vital rates in space
pdf('spatial_vital_rates.pdf', height = 7.5, width = 11.25)
  par(mfrow = c(1, 3))
  # start with survival
  plot(x = plot_dat_sur$X[plot_dat_sur$sur == 1], y = plot_dat_sur$Y[plot_dat_sur$sur == 1],
    xlim = c(-20, 140), ylim = c(-345, 46), xlab = 'location (m)', ylab = 'loaction (m)',
    pch = 19, , bty = 'n', tck = 0.015, col = grey(0.7), main = 'survival')
  points(x = plot_dat_sur$X[plot_dat$sur == 0], y = plot_dat_sur$Y[plot_dat$sur == 0],
    pch = 19)
  points(knot_cl, col = 'red', pch = 3, cex = 0.25)
mtext('a)', side = 3, adj = 0)

  #height, color by max height
  col_height = sequential_hcl(max(round(plot_dat$height)), h = 260, c. = c(0, 100), l = c(10, 50))
  cols_plotted = col_height[round(plot_dat$height)]
  plot(x = plot_dat$X, y = plot_dat$Y, xlim = c(-20, 140), ylim = c(-345, 46),
    xlab = 'location (m)', ylab = '', pch = 19, , bty = 'n', tck = 0.015,
    col = cols_plotted, main = 'max height')
mtext('b)', side = 3, adj = 0)

  #reproduction, color by number of years reproductive
  col_height = sequential_hcl(max(plot_dat$rep) + 1, h = 260, c. = c(0, 100), l = c(10, 50))
  cols_plotted = col_height[plot_dat$rep + 1]
  plot(x = plot_dat$X, y = plot_dat$Y, xlim = c(-20, 140), ylim = c(-345, 46),
    xlab = 'location (m)', ylab = '', pch = 19, , bty = 'n', tck = 0.015,
    col = cols_plotted, main = 'reproduction')
mtext('c)', side = 3, adj = 0)
dev.off()


hist(plot_dat$num_years)  # not many indviduals with only one uear of data, and lots with 3, 4, and 5

# need to find a proxy for survival probability. A good candidate might be max height observed. It can be
# observed for each indivudual. To see if it is a good proxy take the all the idividuals that were observed >10cm in the
# first year and where their death was observed.

vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)

year1_height = vr_loc_gen[vr_loc_gen$age == 0, ]
hist(year1_height$height, breaks = seq(0, 50, 1))

first_height = tapply(seq_along(vr_loc_gen$uID), INDEX = vr_loc_gen$uID, FUN = function(x){
  heights = vr_loc_gen$height[x]
  ages = which(vr_loc_gen$age[x] == 0)
  return(c(ifelse(length(ages) > 0, heights[ages], NA), min(vr_loc_gen$sur[x])))
})

suit = sapply(first_height, FUN = function(x) x[1] < 10 & x[2] == 0)
suit = ifelse(is.na(suit), FALSE, suit)
uIDs = names(suit)
uIDs = uIDs[suit]

vr_full_lc = vr_loc_gen[vr_loc_gen$uID %in% uIDs,]
plot_dat = ddply(vr_full_lc, .(uID), summarise, max_age = max(age, na.rm = TRUE),
  max_height = max(height, na.rm = TRUE))

plot(x = plot_dat$max_age, y = log(plot_dat$max_height))
cor(plot_dat$max_age, plot_dat$max_height) # 0.825
h_mod = lm(log(max_height) ~ max_age + I(max_age ^ 2), data = plot_dat)
summary(h_mod)
#plot(h_mod)

new_dat = data.frame(max_age = seq(min(plot_dat$max_age), max(plot_dat$max_age), 0.1))
pred_age = predict(h_mod, newdata = new_dat, interval = 'prediction', level = 0.95)

pdf('age_vs_height.pdf', height = 8, width = 8)
  plot(x = plot_dat$max_age, y = log(plot_dat$max_height), type = 'n', bty = 'n', xlab = 'max age observed (years)',
    ylab = 'ln(max height observed (cm))', ylim = c(-1, 5))
  polygon(x = c(new_dat$max_age, rev(new_dat$max_age)), y = c(pred_age[,2], rev(pred_age[,3])),
    col = 'lightgoldenrod', border = NA)
  points(x = jitter(plot_dat$max_age, amount = 0.1), y = log(plot_dat$max_height), pch = 19, cex = 0.5)
  lines(x = new_dat$max_age, y = pred_age[,1])
dev.off()
##NOTE: used prediction interval not confidence interval

# high correlation, (R^2 = 0.79 for log model with quadratic term). Lots of spread but max height should act as
# okay porxy for age reached, it will at least discrimnate between first year, second year (to an extent)
# and then 3 and 4 together (line flattens out).

## observed versus predicted for the 4 survival models (no space), 25m grid, fine grid, combined 25m grid and fine grid.
# fitted versus predicted for 4 survival models
ob_name = load('SPP_NK_sur_stan.Rdata')
SPPfine_sur = SPP_NK_sur_stan

## split out the parameter samples and linear predictor samples
SPPfine_sur_pred = extract_flat(SPPfine_sur, pars = c('eta'))

GPPfine = obs_v_pred_bin(obs = sur_dat$sur, lp_samp = SPPfine_sur_pred, num_resamp = 1000, jitter_am = 0.2,
  labels = labs(list(title = 'GPP fine scale', x= 'predicted', y = '')))

# Make violin plot of year effects on prob scale.
GPPfine_sur_param = data.frame(extract_flat(SPPfine_sur, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv', 'sigma_spp')))

GPPfine_sur_year = gather(GPPfine_sur_param, year, effect_size, year_int.1.:year_int.5.)
df_year = rbind(df_year, data.frame(model = rep('GPP fine scale', dim(GPPfine_sur_year)[1]), year = GPPfine_sur_year$year,
  effect_size = GPPfine_sur_year$effect_size))
#GPP25fine_sur_year = gather(GPP25fine_sur_param, year, effect_size, year_int.1.:year_int.5.)
#df_year = rbind(df_year, data.frame(model = rep('GPP 25m res + GPP fine scale', dim(GPP25fine_sur_year)[1]), year = GPP25fine_sur_year$year,
#  effect_size = GPP25fine_sur_year$effect_size))

year_plt = ggplot(df_year, aes(year, effect_size)) +
  geom_violin(fill = hcl(h = 180), color = hcl(h = 180, c = 70)) +
  geom_abline(intercept = 0, slope = 0) + theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9))) +
  scale_x_discrete(labels = c('2003', '2004', '2005', '2006', '2007')) +
  facet_grid(.~model, scales="free", space="free")

year_plt


# plot the decay curves
decay_plotter = function(x, dd_inv){
  dd = 1 / dd_inv
  return(exp(-dd * x))
}

# data for decay plot
x_ax = seq(0, 50, 1)
num_resamp = 200
rs_row = sample.int(dim(GPP25_sur_param)[1], num_resamp)
GPP25_sur_dd_rs = GPP25_sur_param$dd_spp_inv[rs_row]
GPPfine_sur_dd_rs = GPPfine_sur_param$dd_spp_inv[rs_row]

df_dd_plot = data.frame(model = rep(c('GPP 25m res', 'GPP fine scale'), each = length(x_ax) * num_resamp),
  rep_lab = rep(rep(as.character(1:length(GPP25_sur_param$dd_spp_inv)), each = length(x_ax)), times = 2),
  dd_inv = c(rep(GPP25_sur_dd_rs, each = length(x_ax)), rep(GPPfine_sur_dd_rs, each = length(x_ax))),
  dist = rep(x_ax, times = 2 * num_resamp))
df_dd_plot$cor_curve = decay_plotter(df_dd_plot$dist, df_dd_plot$dd_inv)


dd_plt = ggplot(df_dd_plot, aes(dist, cor_curve, group = rep_lab)) +
  geom_line(color = hcl(h = 180, c = 70)) +
  theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9))) +
  facet_grid(model~., scales="free", space="free")

dd_plt

# try a color mat version
GPP25_dd_quant = quantile(GPP25_sur_param$dd_spp_inv, probs = c(0.025, 0.5, 0.975))
GPPfine_dd_quant = quantile(GPPfine_sur_param$dd_spp_inv, probs = c(0.025, 0.5, 0.975))

x_ax = seq(0, 50, 1)
df_dd_quant = data.frame(model = rep(c('GPP 25m res', 'GPP fine scale'), each = length(x_ax) * 2),
  dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)),
  quant_dd = c(rep(GPP25_dd_quant[1], length(x_ax)), rep(GPP25_dd_quant[3], length(x_ax)),
    rep(GPPfine_dd_quant[1], length(x_ax)), rep(GPPfine_dd_quant[3], length(x_ax))),
  med_dd = c(rep(GPP25_dd_quant[2], length(x_ax) * 2), rep(GPPfine_dd_quant[2], length(x_ax) * 2)))

df_dd_quant$quant95 = decay_plotter(df_dd_quant$dist, df_dd_quant$quant_dd)
df_dd_quant$quant50 = decay_plotter(df_dd_quant$dist, df_dd_quant$med_dd)

dd_quant_plt = ggplot(df_dd_quant, aes(dist, quant95)) +
  geom_polygon(aes(fill = hcl(h = 180), color = hcl(h = 180)), alpha = 0.5) +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_line(colour = grey(1)),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none") +
  geom_line(aes(dist, quant50)) + labs(x = 'distance', y = 'correlation') +
  facet_grid(model~., scales="free", space="free")

dd_quant_plt


########################################################################################################################################
MAY DELETE THE ABOVE
############################################################################################################################
## plot to look at performance of models and likelihood, mapped also by year
# get the data sets used to fit the model
sur_dat = read.csv('sur_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
rep_dat = read.csv('rep_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
gr_dat = read.csv('gr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)

# Make a little table to give the numbers of indivdual samples for each year in each vital rates
sur_vr = select(sur_dat, X, Y, year) %>%
	mutate(vital_rate = 'sur')
rep_vr = select(rep_dat, X, Y, year) %>%
	mutate(vital_rate = 'rep')
gr_vr = select(gr_dat, X, Y, year) %>%
	mutate(vital_rate = 'gr')

vr_dat_all = bind_rows(sur_vr, rep_vr, gr_vr) %>%
	group_by(vital_rate, year) %>%
	summarise(n_ob = n())

#	vital_rate  year  n_ob
#      <chr>      <int> <int>
#    1 gr          2003   642
#    2 gr          2004  1273
#    3 gr          2005   784
#    4 gr          2006   482
#    5 gr          2007   360
#    6 rep         2002   711
#    7 rep         2003  1471
#    8 rep         2004  1475
#    9 rep         2005   801
#   10 rep         2006   487
#   11 rep         2007   362
#   12 sur         2003   710
#   13 sur         2004  1468
#   14 sur         2005   997
#   15 sur         2006   799
#   16 sur         2007   486

# Also set up a table to show the number of individuals in each span 
# 2003-2003, 2003-2004, 2003-2005, 2003-2006, 2003-2007
# 2004-2004, 2004-2005, 2004-2006, 2004-2007
# 2005-2005, 2005-2006, 2005-2007
# 2006-2006, 2006-2007
# 2007-2007 
sur_span = select(sur_dat, uID, X, Y, year) %>%
  mutate(vital_rate = 'sur') %>% 
  group_by(vital_rate, uID) %>%
  summarise(
    n_ob = n(),
    X = mean(X),
    Y = mean(Y),
    first_year = min(year),
    last_year = max(year)
  ) %>%
  mutate(span = paste0(first_year, '-', last_year))

rep_span = select(rep_dat, uID, X, Y, year) %>%
  mutate(vital_rate = 'rep') %>% 
  group_by(vital_rate, uID) %>%
  summarise(
    n_ob = n(),
    X = mean(X),
    Y = mean(Y),
    first_year = min(year),
    last_year = max(year)
  ) %>%
  mutate(span = paste0(first_year, '-', last_year))

gr_span = select(gr_dat, uID, X, Y, year) %>%
  mutate(vital_rate = 'gr') %>% 
  group_by(vital_rate, uID) %>%
  summarise(
    n_ob = n(),
    X = mean(X),
    Y = mean(Y),
    first_year = min(year),
    last_year = max(year)
  ) %>%
  mutate(span = paste0(first_year, '-', last_year))


# SURVIVAL
# frequency tables for each vital rate
sur_tab = table(sur_span$first_year, sur_span$last_year)  

# map the survial samples in each year 
sur_span_map = ggplot(sur_span, aes(x = X, y = Y)) +
  geom_point(size = 0.5) +
  theme_light() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  facet_wrap(~span) +
  coord_fixed(ratio = 1.2)

pdf('sur_observation_by_span.pdf', height = 12, width = 12)
  sur_span_map
dev.off()

# GROWTH
gr_tab = table(gr_span$first_year, gr_span$last_year)  

# map the growth samples in each year 
gr_span_map = ggplot(gr_span, aes(x = X, y = Y)) +
  geom_point(size = 0.5) +
  theme_light() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  facet_wrap(~span) +
  coord_fixed(ratio = 1.2)

pdf('gr_observation_by_span.pdf', height = 12, width = 12)
  gr_span_map
dev.off()

# PROB. FLOWER
rep_tab = table(rep_span$first_year, rep_span$last_year)  

# map the prob flower samples in each year 
rep_span_map = ggplot(rep_span, aes(x = X, y = Y)) +
  geom_point(size = 0.5) +
  theme_light() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  facet_wrap(~span) +
  coord_fixed(ratio = 1.2)

pdf('rep_observation_by_span.pdf', height = 12, width = 12)
  rep_span_map
dev.off()


# build up obs vs predicted for spatial and non-spatial model for survival rep and growth
# get the predictions that we need
# survival
ob_name = load('SPP_NK_sur_stan.Rdata')
SPP_sur = SPP_NK_sur_stan
# reproduction
ob_name = load('SPP_NK_rep_stan_ap.Rdata')
SPP_rep = SPP_NK_rep_stan
# growth
ob_name = load('SPP_NK_gr_stan.Rdata')
SPP_gr = SPP_gr_stan

SPP_sur_pred = extract_flat(SPP_sur, pars = c('eta'))
SPP_rep_pred = extract_flat(SPP_rep, pars = c('eta'))
SPP_gr_pred = extract_flat(SPP_gr, pars = c('mu'))
SPP_gr_sigma = extract_flat(SPP_gr, pars = c('sigma'))

pred_obs_sur = data.frame(
  ob = jitter(sur_dat$sur, amount = 0.1),
  pred = as.vector(apply(SPP_sur_pred, MARGIN = 2, FUN = function(x) plogis(mean(x)))))

p <- ggplot(pred_obs_sur, aes(ob, pred))
pvo_sur <- p + geom_point() +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none", plot.title = element_text(size = 20), axis.title = element_text(size = 20),
    strip.background = element_blank(), strip.text = element_blank()) +
  annotate('text', x = 0.5, y = 0.80, label = paste0('logLik ',
    c(round(sum(logLik_bin(post_mat = NS_sur_pred, obs = sur_dat$sur)), 2),
    round(sum(logLik_bin(post_mat = SPP_sur_pred, obs = sur_dat$sur)), 2))), size = 5) +
  labs(title = 'survival', x = 'observed', y = 'predicted') + ylim(0, 1)

pred_obs_rep = data.frame(
  ob = jitter(rep_dat$rep, amount = 0.1),
  pred = as.vector(apply(SPP_rep_pred, MARGIN = 2, FUN = function(x) plogis(mean(x)))))

p <- ggplot(pred_obs_rep, aes(ob, pred))
pvo_rep <- p + geom_point() +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none", plot.title = element_text(size = 20), axis.title = element_text(size = 20),
    strip.background = element_blank(), strip.text = element_blank(),
    axis.title.y = element_blank()) +
  annotate('text', x = 0.5, y = 0.80, label = paste0('logLik ',
    c(round(sum(logLik_bin(post_mat = NS_rep_pred, obs = rep_dat$rep)), 2),
    round(sum(logLik_bin(post_mat = SPP_rep_pred, obs = rep_dat$rep)), 2))), size = 5) +
  labs(title = 'reproduction', x = 'observed', y = 'predicted') + ylim(0, 1) 

pred_obs_gr = data.frame(
  ob = jitter(gr_dat$height, amount = 0.1),
  pred = as.vector(apply(SPP_gr_pred, MARGIN = 2, FUN = function(x) mean(x))))

p <- ggplot(pred_obs_gr, aes(ob, pred))
pvo_gr <- p + geom_point() +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none", plot.title = element_text(size = 20), axis.title = element_text(size = 20),
    axis.title.y = element_blank(), strip.text = element_text(size = 20)) +
  annotate('text', x = 60, y = 10, label = paste0('logLik ',
    c(round(sum(logLik_cont(post_mu = NS_gr_pred, obs = gr_dat$height, post_sd = NS_gr_sigma)), 2),
    round(sum(logLik_cont(post_mu = SPP_gr_pred, obs = gr_dat$height, post_sd = SPP_gr_sigma)), 2))), size = 5) +
  labs(title = 'growth', x = 'observed', y = 'predicted') + geom_abline(intercept = 0, slope = 1) + ylim(0, 65) 

pdf('prev_vs_obs.pdf', height = 10, width = 15)
  grid.arrange(pvo_sur, pvo_rep, pvo_gr, nrow = 1, ncol = 3)
dev.off()

# residual plot gr
gr_df = data.frame(obs = gr_dat$height, height_prev = gr_dat$height_prev, 
    pred = as.vector(apply(SPP_gr_pred, MARGIN = 2, FUN = function(x) mean(x)))) %>%
  mutate(residual = pred - obs)

gr_res_plt = ggplot(gr_df, aes(x = pred, y = residual)) +
  geom_point() +
  geom_abline(slope = 0, intercept = 0)

jpeg(file = "gr_residual.jpg", res = 100)
  gr_res_plt
dev.off()

mean_obs = mean(gr_df$obs)
gr_r2 = 1 - (sum(gr_df$residual ^ 2) / sum((mean_obs - gr_df$obs) ^ 2)) 

# map of likliehood for each point for each year for both non-spatial and spatial models
# survival first
sur_pred_map_df = data.frame(X = sur_dat$X, Y = sur_dat$Y, year = sur_dat$year,
  pred = apply(plogis(SPP_sur_pred), MARGIN = 2, FUN = mean))

sur_pmap = ggplot(sur_pred_map_df, aes(X, Y)) + geom_point(aes(color = pred), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) +
  labs(color = 'prob surv') +
  annotate('segment', x = 0, xend = 50, y =-250, yend = -250, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 25, y =-240, label = c('50m', '', '', '', '')) +
  facet_grid(.~year, scales="free", space="free")

sur_lik_map_df = data.frame(X = sur_dat$X, Y = sur_dat$Y, year = sur_dat$year,
  logLik = logLik_bin(post_mat = SPP_sur_pred, obs = sur_dat$sur))

sur_lmap = ggplot(sur_lik_map_df, aes(X, Y)) + geom_point(aes(colour = exp(logLik)), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) +
  labs(color = 'likliehood') + facet_grid(.~year, scales="free", space="free")

pdf('sur_prev_lik_map.pdf', height = 8, width = 13)
  grid.arrange(sur_pmap, sur_lmap, ncol = 1, nrow = 2)
dev.off()

# make a zomed in verison see how patches look closer in
sur_pmap_zoom = sur_pmap + xlim(30,  45) + ylim(-55, -20) +
  annotate('segment', x = 31, xend = 36, y =-21, yend = -21, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 33.5, y =-20, label = c('5m', '', '', '', ''))

sur_lmap_zoom = sur_lmap + xlim(30,  45) + ylim(-55, -20)

pdf('sur_prev_lik_zoom_map.pdf', height = 8, width = 13)
  grid.arrange(sur_pmap_zoom, sur_lmap_zoom, ncol = 1, nrow = 2)
dev.off()

# reproduction
rep_pred_map_df = data.frame(X = rep_dat$X, Y = rep_dat$Y, year = rep_dat$year,
  pred = apply(plogis(SPP_rep_pred), MARGIN = 2, FUN = mean))

rep_pmap = ggplot(rep_pred_map_df, aes(X, Y)) + geom_point(aes(color = pred), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) +
  labs(color = 'prob repo') +
  annotate('segment', x = 0, xend = 50, y =-250, yend = -250, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 25, y =-240, label = c('50m', '', '', '', '', '')) +
  facet_grid(.~year, scales="free", space="free")

rep_lik_map_df = data.frame(X = rep_dat$X, Y = rep_dat$Y, year = rep_dat$year,
  logLik = logLik_bin(post_mat = SPP_rep_pred, obs = rep_dat$rep))

rep_lmap = ggplot(rep_lik_map_df, aes(X, Y)) + geom_point(aes(colour = exp(logLik)), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) +
  labs(color = 'likliehood') + facet_grid(.~year, scales="free", space="free")

pdf('rep_pred_lik_map.pdf', height = 8, width = 13)
  grid.arrange(rep_pmap, rep_lmap, ncol = 1, nrow = 2)
dev.off()

# make a zoomed in verison see how patches look closer in
rep_pmap_zoom = rep_pmap + xlim(30,  45) + ylim(-55, -20) +
  annotate('segment', x = 31, xend = 36, y =-21, yend = -21, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 33.5, y =-20, label = c('5m', '', '', '', '', ''))

rep_lmap_zoom = rep_lmap + xlim(30,  45) + ylim(-55, -20)

pdf('rep_pred_lik_zoom_map.pdf', height = 8, width = 13)
  grid.arrange(rep_pmap_zoom, rep_lmap_zoom, ncol = 1, nrow = 2)
dev.off()

# growth
gr_pred_map_df = data.frame(X = gr_dat$X, Y = gr_dat$Y, year = gr_dat$year,
  pred = apply(SPP_gr_pred, MARGIN = 2, FUN = mean))

gr_pmap = ggplot(gr_pred_map_df, aes(X, Y)) + geom_point(aes(color = pred), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) +
  labs(color = 'height (cm)') +
  annotate('segment', x = 0, xend = 50, y =-250, yend = -250, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 25, y =-240, label = c('50m', '', '', '', '')) +
  facet_grid(.~year, scales="free", space="free")

gr_lik_map_df = data.frame(X = gr_dat$X, Y = gr_dat$Y, year = gr_dat$year,
  logLik = logLik_cont(post_mu = SPP_gr_pred, obs = gr_dat$height, post_sd = SPP_gr_sigma))

gr_lmap = ggplot(gr_lik_map_df, aes(X, Y)) + geom_point(aes(colour = exp(logLik)), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) +
  labs(color = 'likliehood') + facet_grid(.~year, scales="free", space="free")

pdf('gr_pred_lik_map.pdf', height = 8, width = 13)
  grid.arrange(gr_pmap, gr_lmap, ncol = 1, nrow = 2)
dev.off()

# make a zomed in verison see how patches look closer in
gr_pmap_zoom = gr_pmap + xlim(30,  45) + ylim(-55, -20) +
  annotate('segment', x = 31, xend = 36, y =-21, yend = -21, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 33.5, y =-20, label = c('5m', '', '', '', ''))

gr_lmap_zoom = gr_lmap + xlim(30,  45) + ylim(-55, -20)

pdf('gr_pred_lik_zoom_map.pdf', height = 8, width = 13)
  grid.arrange(gr_pmap_zoom, gr_lmap_zoom, ncol = 1, nrow = 2)
dev.off()

## todo mapp the spp term for each vital rate, year invariant so only need 3 plots, all need a differnt z scale so cant facet
# spatial effect plot
# pull the spp from the models
sur_spp = apply(extract_flat(SPP_sur, pars = c('spp')), MARGIN = 2, FUN = mean)
sur_kl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5)
# the height and year effects to plot the actual probability for an average sized indivdual in and average year
# IMPORTANT: recall that the models are fit to centered hieght.
sur_ave_year = mean(apply(extract_flat(SPP_sur, pars = c('year_int')), MARGIN = 2, FUN = mean))

rep_spp = apply(extract_flat(SPP_rep, pars = c('spp')), MARGIN = 2, FUN = mean)
rep_kl = knot_coords2(dat_x = rep_dat$X, dat_y = rep_dat$Y, min_dist = 0.5)
# the height and year effects to plot the actual probability for an average sized indivdual in and average year
# IMPORTANT: recall that the models are fit to centered hieght.
rep_ave_year = mean(apply(extract_flat(SPP_rep, pars = c('year_int')), MARGIN = 2, FUN = mean))

gr_spp = apply(extract_flat(SPP_gr, pars = c('spp')), MARGIN = 2, FUN = mean)
gr_kl = knot_coords2(dat_x = gr_dat$X, dat_y = gr_dat$Y, min_dist = 0.5)
# the height and year effects to plot the actual probability for an average sized indivdual in and average year
gr_ave_year = mean(apply(extract_flat(SPP_gr, pars = c('b0')), MARGIN = 2, FUN = mean))
gr_ave_ef = mean(apply(extract_flat(SPP_gr, pars = c('gr_rate')), MARGIN = 2, FUN = mean))
gr_ave_height = mean(gr_dat$height_prev)

sur_dd = extract_flat(SPP_sur, pars = c('dd_spp_inv'))
rep_dd = extract_flat(SPP_rep, pars = c('dd_spp_inv'))
gr_dd = extract_flat(SPP_gr, pars = c('dd_spp_inv'))


# spp on survial big map
sur_spp_map_df = data.frame(X = sur_kl[,1], Y = sur_kl[,2], spp = sur_spp,
  spp_prob_sur = plogis(sur_ave_year + sur_spp))

sur_spp_map = ggplot(sur_spp_map_df, aes(X, Y)) + geom_point(aes(color = spp_prob_sur), size = 0.5) +
  xlim(-20, 140) + ylim(-180, 40) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), plot.title = element_text(size = 20),
    legend.position = c(0.95, 0.95), legend.justification = c(1, 1),
    legend.title = element_text(size = 20), legend.text = element_text(size = 17)) +
  labs(title = 'Survival', colour = 'survival\n(prob.)') + scale_color_viridis() +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m', size = 8) +
  coord_fixed(ratio = 1.375)

# reproduction big map
rep_spp_map_df = data.frame(X = rep_kl[,1], Y = rep_kl[,2], spp = rep_spp,
  spp_prob_rep = plogis(rep_ave_year + rep_spp))

rep_spp_map = ggplot(rep_spp_map_df, aes(X, Y)) + 
  geom_point(aes(color = spp_prob_rep), size = 0.5) + 
  xlim(-20, 140) + ylim(-180, 40) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), plot.title = element_text(size = 20),
    legend.position = c(0.95, 0.95), legend.justification = c(1, 1),
    legend.title = element_text(size = 20), legend.text = element_text(size = 17)) +
  labs(title = 'Prob of flowering', colour = 'prob.\nflower') + scale_color_viridis() +
  annotate("rect", xmin = 10.5, xmax = 43, ymin = -58, ymax = 13, colour = grey(0), alpha = 0.0) +
  annotate("text", x = 10, y = 20, label = 'panels d,e,f', hjust = 0, size = 6) +
  coord_fixed(ratio = 1.375)

# growth big map # we use the ave height for survival data, which is 32cm,
# we usde that value so the legends are consistant across the text
gr_spp_map_df = data.frame(X = gr_kl[,1], Y = gr_kl[,2], spp = gr_spp,
  spp_ave_scale = gr_ave_year + gr_ave_ef * 32 + gr_spp)

gr_spp_map = ggplot(gr_spp_map_df, aes(X, Y)) + geom_point(aes(color = spp_ave_scale), size = 0.5) +
	xlim(-20, 140) + ylim(-180, 40) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(),
    axis.text.y = element_blank(), plot.title = element_text(size = 20),
    legend.position = c(1, 0.95), legend.justification = c(1,1), 
    legend.title = element_text(size = 20), legend.text = element_text(size = 17)) +
  labs(title = 'Growth', color = "height t+1 \n (cm)") + scale_color_viridis() +
  coord_fixed(ratio = 1.375)

# zoomed in areas
sur_spp_zoom_map = sur_spp_map  + geom_point(aes(color = spp_prob_sur), size = 2) +
	xlim(13.5,  40) + ylim(-50, 5.2) + labs(title = '') +
	annotate('segment', x = 30, xend = 35, y = 5, yend = 5, color = grey(0)) +
	annotate('text', x = 32.5, y = 2, label = '5 m', size = 8) + theme(legend.position = 'none') +
  coord_fixed(ratio = 2.075)
rep_spp_zoom_map = rep_spp_map + geom_point(aes(color = spp_prob_rep), size = 2) +
	xlim(13.5,  40) + ylim(-50, 5.2) + labs(title = '') + theme(legend.position = 'none') +
  coord_fixed(ratio = 2.075)
gr_spp_zoom_map = gr_spp_map + geom_point(aes(color = spp_ave_scale), size = 2) +
	xlim(13.5,  40) + ylim(-50, 5.2) + labs(title = '') + theme(legend.position = 'none') +
  annotate("segment", x = 22, xend = 24.5, y = 3.5, yend = 3.5, colour = grey(0)) +
  annotate("segment", x = 24.5, xend = 24.5, y = 1, yend = 5.2, colour = grey(0)) +
  annotate("text", x = 25, y = 3.5, label = 'example\ngap', hjust = 0, size = 5) +
  coord_fixed(ratio = 2.075)

# make two plots, so the zoomed in ones are rectangular
patch_plt = plot_grid(sur_spp_map, gr_spp_map, rep_spp_map,
  ncol = 3, label_size = 18, labels  = paste0(letters[1:3], ')'))

zoom_plt = plot_grid(sur_spp_zoom_map, gr_spp_zoom_map, rep_spp_zoom_map,
  ncol = 3, label_size = 18, labels  = paste0(letters[4:6], ')'))

dev.off()

pdf('SPP_plot_true_aspect.pdf', height = 15, width = 10)
  plot_grid(patch_plt, zoom_plt, nrow = 2, rel_heights = c(1, 1.3))
dev.off()

# split the distance decay plots

# show the distance decay curves
# plot the decay curves
decay_plotter = function(x, dd_inv){
  dd = 1 / dd_inv
  return(exp(-dd * x))
}

sur_dd_quant = quantile(sur_dd, probs = c(0.025, 0.5, 0.975))

x_ax = seq(0, 40, 1)
sur_dd_df = data.frame(dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)),
    quant_dd = c(rep(sur_dd_quant[1], length(x_ax)), rep(sur_dd_quant[3], length(x_ax))),
    med_dd = rep(sur_dd_quant[2], length(x_ax) * 2)) %>%
  mutate(vr = 'Survival',
    quant95 = decay_plotter(dist, quant_dd), 
    quant50 = decay_plotter(dist, med_dd))

rep_dd_quant = quantile(rep_dd, probs = c(0.025, 0.5, 0.975))
rep_dd_df = data.frame(dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)),
    quant_dd = c(rep(rep_dd_quant[1], length(x_ax)), rep(rep_dd_quant[3], length(x_ax))),
    med_dd = rep(rep_dd_quant[2], length(x_ax) * 2)) %>%
  mutate(vr = 'Prob. flowering',
    quant95 = decay_plotter(dist, quant_dd), 
    quant50 = decay_plotter(dist, med_dd))

gr_dd_quant = quantile(gr_dd, probs = c(0.025, 0.5, 0.975))
gr_dd_df = data.frame(dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)),
    quant_dd = c(rep(gr_dd_quant[1], length(x_ax)), rep(gr_dd_quant[3], length(x_ax))),
    med_dd = rep(gr_dd_quant[2], length(x_ax) * 2))%>%
  mutate(vr = 'Growth',
    quant95 = decay_plotter(dist, quant_dd), 
    quant50 = decay_plotter(dist, med_dd))

dd_plt_dat = bind_rows(sur_dd_df, rep_dd_df, gr_dd_df)

dd_plt = ggplot(dd_plt_dat, aes(dist, quant95)) +
  geom_polygon(aes(fill = hcl(h = 180), color = hcl(h = 180)), alpha = 0.5) +
  theme_grey() +
  theme(legend.position = "none", 
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()) +
  geom_line(aes(dist, quant50)) + 
  labs(x = 'distance (m)', y = 'correlation') +
  facet_grid(.~factor(vr, levels = c('Survival', 'Growth', 'Prob. flowering')))

pdf(file = 'dist_decay_plt.pdf', height = 5, width = 15)
  dd_plt
dev.off()

## parameter plots
# make violine plots of the coeffiecents
SPP_sur_param = data.frame(extract_flat(SPP_sur, pars = c('h_ef', 'year_int', 'sigma_spp')))
SPP_rep_param = data.frame(extract_flat(SPP_rep, pars = c('h_ef', 'year_int', 'sigma_spp')))


SPP_sur_long = mutate(SPP_sur_param, h_ef_mh = h_ef * quantile(sur_dat$height_mod, probs = 0.75, na.rm = TRUE)) %>%
	gather(param, effect_size, h_ef:h_ef_mh)

#make the slope height effect of an individual of mean height to put on same scale as intercepts
SPP_sur_long[SPP_sur_long$param == 'h_ef', 'effect_size'] = SPP_sur_long[SPP_sur_long$param == 'h_ef', 'effect_size'] *
  quantile(sur_dat$height_mod, probs = 0.75)
SPP_sur_long$model = 'SPP'
SPP_sur_long$vrate = 'Survival'
#only keep the 99% central part of the distribution
SPP_sur_long = SPP_sur_long[unlist(tapply(SPP_sur_long$effect_size, INDEX = SPP_sur_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

SPP_rep_long = mutate(SPP_rep_param, h_ef_mh = h_ef * quantile(sur_dat$height_mod, probs = 0.75, na.rm = TRUE)) %>%
	gather(param, effect_size, h_ef:h_ef_mh)

# recode years so they start at year 0 for reproduction
for(i in 1:6){
  SPP_rep_long[SPP_rep_long$param == paste0('year_int.', i, '.'), 'param'] = paste0('year_int.', i - 1, '.')
}

# make growth seperate since there is a different set of predictors and the effect sizes are on a different scale
SPP_gr_param = data.frame(extract_flat(SPP_gr, pars = c('b0', 'gr_rate', 'sigma_spp', 'sigma')))
#make the slope height effect of an individual of mean height to put on same scale as intercepts
SPP_gr_param[, 6:10] = SPP_gr_param[, 6:10] * mean(gr_dat$height_prev)

SPP_gr_long = gather(SPP_gr_param, param, effect_size, b0.1.:sigma) %>%
	mutate(effect_mh = ifelse(grepl('gr_rate', param),
		effect_size * quantile(sur_dat$height, probs = 0.75, na.rm = TRUE), effect_size))

SPP_gr_long = SPP_gr_long[unlist(tapply(SPP_gr_long$effect_size, INDEX = SPP_gr_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

# make a version of the above but only for the spatial model
par_s_plt = ggplot(SPP_sur_long , aes(param, effect_size)) + geom_abline(intercept = 0, slope = 0) +
  geom_violin(fill = 'blue', color = 'blue', alpha = 0.2) +
  theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9)),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 20), axis.title.y = element_blank(),
    axis.text.x = element_text(size = 18), axis.text.y = element_blank(),
	axis.line.y = element_blank(), axis.ticks.y = element_blank(),
	plot.title = element_text(size = 20)) +
  labs(x = '\n\nParameter', y = ' ', title = 'Survival', size = 20) +
  scale_x_discrete(labels = c(bquote(''*beta^s*''), bquote(''*beta["03"]^s*''),
  	bquote(''*beta["04"]^s*''), bquote(''*beta["05"]^s*''), bquote(''*beta["06"]^s*''),
	bquote(''*beta["07"]^s*'')),
    limits = c('h_ef_mh', 'year_int.1.', 'year_int.2.', 'year_int.3.', 'year_int.4.', 'year_int.5.')) +
	scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(-0.4, 6.3))

par_r_plt = ggplot(SPP_rep_long , aes(param, effect_size)) + geom_abline(intercept = 0, slope = 0) +
  geom_violin(fill = 'blue', color = 'blue', alpha = 0.2) +
  theme(
	    panel.grid.major.y = element_line(colour = grey(1)),
	    panel.grid.major.x = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.background = element_rect(fill = grey(0.9)),
	    strip.background = element_blank(),
	    strip.text.y = element_blank(),
	    strip.text.x = element_text(size = 20),
	    axis.title = element_text(size = 20),
		axis.text = element_text(size = 18),
		plot.title = element_text(size = 20)) +
	labs(x = ' \n ', y = 'Effect size', title = 'Probability of flowering', size = 20) +
	scale_x_discrete(labels = c(bquote(''*beta^r*''), bquote(''*beta["02"]^r*''), bquote(''*beta["03"]^r*''),
    bquote(''*beta["04"]^r*''), bquote(''*beta["05"]^r*''), bquote(''*beta["06"]^r*''),
	bquote(''*beta["07"]^r*'')),
		limits = c('h_ef_mh', 'year_int.0.', 'year_int.1.', 'year_int.2.', 'year_int.3.', 'year_int.4.', 'year_int.5.')) +
	scale_y_continuous(breaks = c(0, 2, 4, 6), limits = c(-0.4, 6.3))

par_gr_plt = ggplot(SPP_gr_long, aes(param, effect_size)) + geom_abline(intercept = 0, slope = 0) +
  geom_violin(fill = 'blue', color = 'blue', alpha = 0.2) +
  labs(x = ' ', y = '', title = 'Growth', size = 20)+
  theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9)),
    strip.background = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
	plot.title = element_text(size = 20)) +
  scale_x_discrete(labels = c(bquote(''*beta["3"]^g*''), bquote(''*beta["4"]^g*''),
    bquote(''*beta["5"]^g*''), bquote(''*beta["6"]^g*''), bquote(''*beta["7"]^g*''),
	bquote(''*beta["03"]^g*''), bquote(''*beta["04"]^g*''), bquote(''*beta["05"]^g*''),
	bquote(''*beta["06"]^g*''), bquote(''*beta["07"]^g*'')),
	    limits = c('gr_rate.1.', 'gr_rate.2.', 'gr_rate.3.', 'gr_rate.4.', 'gr_rate.5.',
			'b0.1.', 'b0.2.', 'b0.3.', 'b0.4.', 'b0.5.')) + labs(x = ' ', y = ' ')

# all three plots toghter
comb_plt = plot_grid(par_r_plt, par_s_plt, par_gr_plt, align = 'h', ncol = 3)

# make some data for the boxes and annotation
ann_box = data.frame(x = c(0.035, 0.09, 0.35, 0.405, 0.7, 0.845),
	x_right = c(0.085, 0.32, 0.4, 0.65, 0.84, 0.99), y = 0.065, y_top = 0.24)

pdf('pars_vr_plot_spp_only.pdf', height = 6, width = 16)
	ggdraw() +
	geom_rect(data = ann_box, aes(xmin = x, xmax = x_right, ymin = y, ymax = y_top),
		fill = grey(0.97), colour = grey(0.4)) +
	geom_text(data = ann_box, aes(x = (x_right + x) / 2, y = (y_top + y) / 2), vjust = 0.8,
 		label = c('height\neffect', 'year\neffects', 'height\neffect', 'year\neffects',
			'height\neffects', 'year\neffects'), size = 6) +
	draw_plot(comb_plt)
dev.off()

# make a plot of within vs between patch relationship of the estimated growth vs survival
# first step is get the data from stan out put tp get the GPP

################################################################################
# Make some plots to show how the locations relate to the data and each other
# get data
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactor = FALSE)
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

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

min_l2l = apply(dm_knots, MARGIN = 1, FUN = function(x) min(x[x > 0]))
min_l2d = apply(dm_knots_obs, MARGIN = 1, FUN = min)

plt_df_l2l = data.frame(min_dist = min_l2l,
	dist_mat = rep('location to location', length(min_l2l)))
plt_df_l2d = data.frame(min_dist = min_l2d,
	dist_mat = rep('location to observation', length(min_l2d)))

plt_l2l = ggplot(plt_df_l2l, aes(x = min_dist)) + geom_histogram(aes(y = ..density..)) +
	theme_grey() +
	scale_x_log10()+
	theme(axis.title = element_blank(),
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank(),
		panel.grid.minor = element_blank()) +
	facet_grid(.~ dist_mat, scales = 'free')

plt_l2d = ggplot(plt_df_l2d, aes(x = min_dist)) + geom_histogram(aes(y = ..density..)) +
	theme_grey() +
	theme(axis.title = element_blank(),
		axis.text.y = element_blank(),
		axis.ticks.y = element_blank(),
		panel.grid.minor = element_blank()) +
	facet_grid(.~ dist_mat, scales = 'free')

plt_pair = plot_grid(plt_l2l, plt_l2d, align = 'v', labels = c('a)', 'b)'))

pdf(file = 'loc_dist_hists.pdf', height = 4, width = 8)
	ggdraw(add_sub(plt_pair, 'distance to nearest location (m)',
		vpadding = grid::unit(0.5, 'lines'), y = 6, x = 0.5, vjust = 9.5))
dev.off()







#
