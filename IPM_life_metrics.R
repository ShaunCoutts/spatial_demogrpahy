# calculate some demogrpahic metrics, in particular R0, and sensitivity of these
# metrics. Need to do this sample wise so set these functions to take and
# return vectors of sampled parameters.

# expected fecundity at each size in Z, modelled using a Poisson glm with a
# log link-function, so the linear predictor is log(expected number), exp this to
# get expected number. Looks quiet nice in terms of expected seeds when compaired to
# data

# fec_mode = function(int_f, beta_f, cent_const, Z){

#   return(exp(int_f + beta_f * (Z - cent_const)))

# }

# use code from the IPM book Chapter 3 to calculate R0, mean life expectancy
# and mean age of reproduction for a single sample, to save recalculation of
# N and S and Fec, return as a named list
life_cycle_metrics = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	age, Z, n0){

	# make the fecundity kernel, vector of seeds produced by each size in Z
	# prob of flowering
    lp = rep_int + rep_slope * (Z - cent_rep) + rep_G
    prob_flower = 1 / (1 + exp(-lp))

	# size effect on seed production by the probability of flowering
	Fec = matrix(fec_mode(int_f, slope_f, cent_f, Z) * prob_flower, nrow = 1)

	# make the growth-survival kernel, matrix of |Z| by |Z|
	# calulualte surival for each size in Z
	sur_lp = sur_int + sur_slope * (Z - cent_sur) + sur_G
  	S = 1 / (1 + exp(-sur_lp))

	# do the growth
	mean_growth = grow_int + grow_slope * Z + grow_G
	mu = expand.grid(Z, mean_growth)
	G = matrix(dnorm(mu[, 1], mu[, 2], grow_sigma),
	    nrow = length(Z), byrow = TRUE)

	# make growth conditional on survival, transpose to make P from col to row
	# so the multiplication with Fec works
	P = t(G * S)

	# construct the identiy matrix
	I = diag(length(Z))

	# Make the Neumann operator
	N = solve(I - P)

	R = Fec %*% N

	#Define the e vector which we will use for summing down columns
	# don't think I need since my Fec kernel already a vector
	e = matrix(1, nrow = 1, ncol = dim(P)[1])

	R0 = sum(R * n0)
	mean_life_expect = sum((e %*% N) * n0)

	# get mean seed production a given age
	Pa = P
	if(age > 1){

		for(a in 2:age){

			Pa = Pa %*% P

		}

	}

	# survival to each age used to turn reproduction to per capita
	la = sum((e %*% Pa) * n0)
	mean_rep_age = sum((Fec %*% Pa) * n0) / la

	return(list(R0 = R0, life_expect = mean_life_expect, mean_rep = mean_rep_age))

}

# Faster version that just does R0
R0_asym = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	age, Z, n0){

	# make the fecundity kernel, vector of seeds produced by each size in Z
	# prob of flowering
    lp = rep_int + rep_slope * (Z - cent_rep) + rep_G
    prob_flower = 1 / (1 + exp(-lp))

	# size effect on seed production by the probability of flowering
	Fec = matrix(fec_mode(int_f, slope_f, cent_f, Z) * prob_flower, nrow = 1)

	# make the growth-survival kernel, matrix of |Z| by |Z|
	# calulualte surival for each size in Z
	sur_lp = sur_int + sur_slope * (Z - cent_sur) + sur_G
  	S = 1 / (1 + exp(-sur_lp))

	# do the growth
	mean_growth = grow_int + grow_slope * Z + grow_G
	mu = expand.grid(Z, mean_growth)
	G = matrix(dnorm(mu[, 1], mu[, 2], grow_sigma),
	    nrow = length(Z), byrow = TRUE)

	# make growth conditional on survival, transpose to make P from col to row
	# so the multiplication with Fec works
	P = t(G * S)

	# construct the identiy matrix
	I = diag(length(Z))

	# Make the Neumann operator
	N = solve(I - P)

	R = Fec %*% N

	# Define the e vector which we will use for summing down columns
	# don't think I need since my Fec kernel already a vector
	# e = matrix(1, nrow = 1, ncol = dim(P)[1])

	R0 = sum(R * n0)

	return(R0)

}

# Sample version to return a dataframe of life cycle metrics incorporating uncertianty 
# in parameter estimation
life_cycle_samp = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	age, Z, n0){

	life_metrics_samps = list()

	for(i in seq_along(rep_int)){

		life_metrics_samps[[i]] = life_cycle_metrics(rep_int[i], rep_slope[i], rep_G[i], cent_rep, int_f[i],
			slope_f[i], cent_f, sur_int[i], sur_slope[i], sur_G[i], cent_sur,
			grow_int[i], grow_slope[i], grow_G[i], grow_sigma[i], age, Z, n0)

	}

	return(bind_rows(life_metrics_samps))

}

# Sample version to return a vector of R0 incorporating uncertianty in parameter
# estimation
R0_asym_samp = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	age, Z, n0){

	R0_samps = numeric()

	for(i in seq_along(rep_int)){

		R0_samps[i] = R0_asym(rep_int[i], rep_slope[i], rep_G[i], cent_rep, int_f[i],
			slope_f[i], cent_f, sur_int[i], sur_slope[i], sur_G[i], cent_sur,
			grow_int[i], grow_slope[i], grow_G[i], grow_sigma[i], age, Z, n0)

	}

	return(R0_samps)

}

################################################################################
# Use the finite difference method to find sensitivity to intercept in each
# vital rate (IPM book pp. 219)
# Error for the difference
delta = .Machine$double.eps ^ (1/3)

R0_sens = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	age, Z, n0, Delta = delta, vital_rate){

	sense = numeric(length(rep_int))

	if(!vital_rate %in% c('rep', 'sur', 'grow')){
		stop("vital_rate needs to be one of [rep, sur, grow]")
	}

	if(vital_rate == 'rep'){

		rep_int_up = rep_int + Delta
		sur_int_up = sur_int
		grow_int_up = grow_int

		rep_int_down = rep_int - Delta
		sur_int_down = sur_int
		grow_int_down = grow_int

		pr = rep_int

	}

	if(vital_rate == 'sur'){

		rep_int_up = rep_int
		sur_int_up = sur_int + Delta
		grow_int_up = grow_int

		rep_int_down = rep_int
		sur_int_down = sur_int - Delta
		grow_int_down = grow_int

		pr = sur_int

	}

	if(vital_rate == 'grow'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int + Delta

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int - Delta

		pr = grow_int

	}

	mets_up = life_cycle_samp(rep_int_up, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
		sur_int_up, sur_slope, sur_G, cent_sur, grow_int_up, grow_slope, grow_G, grow_sigma,
		age, Z, n0)

	mets_down = life_cycle_samp(rep_int_down, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
		sur_int_down, sur_slope, sur_G, cent_sur, grow_int_down, grow_slope, grow_G, grow_sigma,
		age, Z, n0)

	R0_sen = (mets_up$R0 - mets_down$R0) / (2 * delta)

	delta_ln_R0 = log(mets_up$R0) - log(mets_down$R0)

	if(vital_rate == 'rep'){

		# elasticity to intercpet
		delta_ln = log(rep_int_up) - log(rep_int_down)
		R0_elas = delta_ln_R0 / delta_ln

		# elacticity to mean reproduction in first year
		delta_ln_pr = log(mets_up$mean_rep) - log(mets_down$mean_rep)
		R0_elas_pr = delta_ln_R0 / delta_ln_pr

		sen_inter = (mets_up$mean_rep - mets_down$mean_rep) / (2 * Delta)

	}

	if(vital_rate == 'sur'){

		# elasticity to intercpet
		delta_ln = log(sur_int_up) - log(sur_int_down)
		R0_elas = delta_ln_R0 / delta_ln

		# elacticity to life expectancy
		delta_ln_pr = log(mets_up$life_expect) - log(mets_down$life_expect)
		R0_elas_pr = (log(mets_up$R0) - log(mets_down$R0)) / delta_ln_pr

		sen_inter = (mets_up$life_expect - mets_down$life_expect) / (2 * Delta)

	}

	if(vital_rate == 'grow'){

		# elasticity to intercpet
		delta_ln = log(grow_int_up) - log(grow_int_down)
		R0_elas = delta_ln_R0 / delta_ln

		R0_elas_pr = NA
		sen_inter = NA

	}

	return(list(R0_sen = R0_sen, R0_elas = R0_elas, R0_elas_pr = R0_elas_pr,
			interm_sen = sen_inter))

}

################################################################################
# put this all toghter and return a data frame with summaries of the R0 and sensitivity
# for each location
R0_summary_df = function(rep_samps, sur_samps, gr_samps, age, Z, n0){

	IPM_summary = data.frame(uni_spp = gr_samps$SPP$uni_spp, uX = gr_samps$SPP$unify_X,
		uY = gr_samps$SPP$unify_Y, R0_mean = NA, R0_med = NA, R0_lq95 = NA, R0_uq95 = NA,
		life_expect = NA, life_expect_lq95 = NA, life_expect_uq95 = NA,
		mean_repA = NA, mean_rep_lq95 = NA, mean_rep_uq95 = NA,
		sen_rep = NA, sen_sur = NA, sen_gr = NA, sen_life_ex = NA, sen_mean_rep = NA,
		elas_rep = NA, elas_sur = NA, elas_gr = NA, elas_rep_fl = NA, elas_sur_le = NA)

	# go through every location and get the distirbution of R0
	for(loc in seq_along(IPM_summary$uni_spp)){

		# get the SPP samples
		sur_SPP = filter(sur_samps$SPP, uni_spp == IPM_summary$uni_spp[loc]) %>%
			select(s1:s1000) %>% as.numeric()
		rep_SPP = filter(rep_samps$SPP, uni_spp == IPM_summary$uni_spp[loc]) %>%
			select(s1:s1000) %>% as.numeric()
		gr_SPP = filter(gr_samps$SPP, uni_spp == IPM_summary$uni_spp[loc]) %>%
			select(s1:s1000) %>% as.numeric()

		life_met_samples = life_cycle_samp(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age = age,
			Z = Z,
			n0 = n0)

		IPM_summary$R0_mean[loc] = mean(life_met_samples$R0)
		IPM_summary$R0_med[loc] = median(life_met_samples$R0)
		IPM_summary$R0_lq95[loc] = quantile(life_met_samples$R0, probs = 0.025)
		IPM_summary$R0_uq95[loc] = quantile(life_met_samples$R0, probs = 0.975)

		IPM_summary$life_expect[loc] = mean(life_met_samples$life_expect)
		IPM_summary$life_expect_lq95[loc] = quantile(life_met_samples$life_expect, probs = 0.025)
		IPM_summary$life_expect_uq95[loc] = quantile(life_met_samples$life_expect, probs = 0.975)

		IPM_summary$mean_repA[loc] = mean(life_met_samples$mean_rep)
		IPM_summary$mean_rep_lq95[loc] = quantile(life_met_samples$mean_rep, probs = 0.025)
		IPM_summary$mean_rep_uq95[loc] = quantile(life_met_samples$mean_rep, probs = 0.975)

		# get the sensitivity to intercept of each vital rate
		rep_sense_samp = R0_sens(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age, Z, n0, vital_rate = 'rep')

		IPM_summary$sen_rep[loc] = mean(rep_sense_samp$R0_sen)
		IPM_summary$elas_rep[loc] = mean(rep_sense_samp$R0_elas)
		IPM_summary$elas_rep_fl[loc] = mean(rep_sense_samp$R0_elas_pr)
		IPM_summary$sen_mean_rep[loc] = mean(rep_sense_samp$interm_sen)

		sur_sense_samp = R0_sens(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age, Z, n0, vital_rate = 'sur')

		IPM_summary$sen_sur[loc] = mean(sur_sense_samp$R0_sen)
		IPM_summary$elas_sur[loc] = mean(sur_sense_samp$R0_elas)
		IPM_summary$elas_sur_le[loc] = mean(sur_sense_samp$R0_elas_pr)
		IPM_summary$sen_life_ex[loc] = mean(sur_sense_samp$interm_sen)

		gr_sense_samp = R0_sens(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age, Z, n0, vital_rate = 'grow')

		IPM_summary$sen_gr[loc] = mean(gr_sense_samp$R0_sen)
		IPM_summary$elas_gr[loc] = mean(gr_sense_samp$R0_elas)

		print(loc)

	}

	return(IPM_summary)

}

################################################################################
# I want to test the constriant between R0 and the elasticities. We have multipule
# estimates of R0 and the elasticites for each plant. If the pattern is constraint
# we should see it at the indvidual plant scale as well. But if it is a real pattern
# in the plant life histories then should be no correlation within each plant since the
# samples are error in the estimate of each parameter, for each plant.
covar_test = function(rep_samps, sur_samps, gr_samps, age, Z, n0){

	IPM_summary = data.frame(uni_spp = gr_samps$SPP$uni_spp, uX = gr_samps$SPP$unify_X,
		uY = gr_samps$SPP$unify_Y, R0_mean = NA, life_expect = NA, mean_repA = NA,
		elas_rep = NA, elas_sur = NA, elas_gr = NA, elas_rep_fl = NA, elas_sur_le = NA,
		C_R0_sur = NA, C_R0_sur_le = NA, C_R0_rep = NA, C_R0_rep_fl = NA, C_R0_gr = NA,
		C_sur_rep = NA, C_sur_gr = NA, C_rep_gr = NA, C_surle_repfl = NA,
		C_surle_gr = NA, C_repfl_gr = NA)

	# go through every location and get the distirbution of R0
	for(loc in seq_along(IPM_summary$uni_spp)){

		# get the SPP samples
		sur_SPP = filter(sur_samps$SPP, uni_spp == IPM_summary$uni_spp[loc]) %>%
			select(s1:s1000) %>% as.numeric()
		rep_SPP = filter(rep_samps$SPP, uni_spp == IPM_summary$uni_spp[loc]) %>%
			select(s1:s1000) %>% as.numeric()
		gr_SPP = filter(gr_samps$SPP, uni_spp == IPM_summary$uni_spp[loc]) %>%
			select(s1:s1000) %>% as.numeric()

		life_met_samples = life_cycle_samp(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age = age, Z = Z, n0 = n0)

		IPM_summary$R0_mean[loc] = mean(life_met_samples$R0, na.rm = TRUE)
		IPM_summary$life_expect[loc] = mean(life_met_samples$life_expect, na.rm = TRUE)
		IPM_summary$mean_repA[loc] = mean(life_met_samples$mean_rep, na.rm = TRUE)

		# get the sensitivity to intercept of each vital rate
		rep_sense_samp = R0_sens(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age, Z, n0, vital_rate = 'rep')

		IPM_summary$elas_rep[loc] = mean(rep_sense_samp$R0_elas, na.rm = TRUE)
		IPM_summary$elas_rep_fl[loc] = mean(rep_sense_samp$R0_elas_pr, na.rm = TRUE)
		IPM_summary$C_R0_rep[loc] = cor(log(life_met_samples$R0), log(rep_sense_samp$R0_elas))
		IPM_summary$C_R0_rep_fl[loc] = cor(log(life_met_samples$R0), log(rep_sense_samp$R0_elas_pr))

		sur_sense_samp = R0_sens(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age, Z, n0, vital_rate = 'sur')

		IPM_summary$elas_sur[loc] = mean(sur_sense_samp$R0_elas, na.rm = TRUE)
		IPM_summary$elas_sur_le[loc] = mean(sur_sense_samp$R0_elas_pr, na.rm = TRUE)
		IPM_summary$C_R0_sur[loc] = cor(log(life_met_samples$R0), log(sur_sense_samp$R0_elas))
		IPM_summary$C_R0_sur_le[loc] = cor(log(life_met_samples$R0), log(sur_sense_samp$R0_elas_pr))
		IPM_summary$C_sur_rep[loc] = cor(log(rep_sense_samp$R0_elas), log(sur_sense_samp$R0_elas))
		IPM_summary$C_surle_repfl[loc] = cor(log(rep_sense_samp$R0_elas_pr), log(sur_sense_samp$R0_elas_pr))

		gr_sense_samp = R0_sens(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age, Z, n0, vital_rate = 'grow')

		IPM_summary$elas_gr[loc] = mean(gr_sense_samp$R0_elas, na.rm = TRUE)
		IPM_summary$C_R0_gr[loc] = cor(log(life_met_samples$R0), log(gr_sense_samp$R0_elas))
		IPM_summary$C_sur_gr[loc] = cor(log(gr_sense_samp$R0_elas), log(sur_sense_samp$R0_elas))
		IPM_summary$C_rep_gr[loc] = cor(log(gr_sense_samp$R0_elas), log(rep_sense_samp$R0_elas))
		IPM_summary$C_surle_gr[loc] = cor(log(gr_sense_samp$R0_elas), log(sur_sense_samp$R0_elas_pr))
		IPM_summary$C_repfl_gr[loc] = cor(log(gr_sense_samp$R0_elas), log(rep_sense_samp$R0_elas_pr))

		print(loc)

	}

	return(IPM_summary)
}

#####################################################################################
# sample R0 calc to do a quicker null model test, only use 500 samples for speed up
null_R0_df = function(par_list){
	
	# unpack the parameter list
	sur_samps = par_list$sur_samp 
	rep_samps = par_list$rep_samp
	gr_samps = par_list$gr_samp 
	age = par_list$age
	Z = par_list$Z
	n0 = par_list$n0 
	rand_num = par_list$rand_num

	print(rand_num)

	null_summary = data.frame(uni_spp = gr_samps$SPP$uni_spp, uX = gr_samps$SPP$unify_X,
		uY = gr_samps$SPP$unify_Y, R0_mean = NA, rand_run = rand_num)

	# randomise the order of the locations, 
	# special case for the observed, order maintained
	if(rand_num == 'obs'){

		sur_loc = seq_along(null_summary$uni_spp)
		rep_loc = seq_along(null_summary$uni_spp)
		gr_loc = seq_along(null_summary$uni_spp)
	
	}else{ #not observed randomise

		sur_loc = sample.int(n = length(null_summary$uni_spp))
		rep_loc = sample.int(n = length(null_summary$uni_spp))
		gr_loc = sample.int(n = length(null_summary$uni_spp))

	}

	# go through every location and get the distirbution of R0
	for(loc in seq_along(null_summary$uni_spp)){

		# get the SPP samples
		sur_SPP = filter(sur_samps$SPP, uni_spp == null_summary$uni_spp[sur_loc[loc]]) %>%
			select(s1:s500) %>% as.numeric()
		rep_SPP = filter(rep_samps$SPP, uni_spp == null_summary$uni_spp[rep_loc[loc]]) %>%
			select(s1:s500) %>% as.numeric()
		gr_SPP = filter(gr_samps$SPP, uni_spp == null_summary$uni_spp[gr_loc[loc]]) %>%
			select(s1:s500) %>% as.numeric()

		R0_sampels = R0_asym_samp(rep_int = rep_samps$int,
			rep_slope = rep_samps$slope,
			rep_G = rep_SPP,
			cent_rep = rep_samps$cent_const,
			int_f = fruit_samps$int,
			slope_f = fruit_samps$slope,
			cent_f = fruit_samps$cent_const,
			sur_int = sur_samps$int,
			sur_slope = sur_samps$slope,
			sur_G = sur_SPP,
			cent_sur = sur_samps$cent_const,
			grow_int = gr_samps$int,
			grow_slope = gr_samps$slope,
			grow_G = gr_SPP,
			grow_sigma = gr_samps$sigma,
			age = age, Z = Z, n0 = n0)

		null_summary$R0_mean[loc] = mean(R0_sampels, na.rm = TRUE)
	
	}

	return(null_summary)

}
