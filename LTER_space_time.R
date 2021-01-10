# Do a LTRE following Caswell 1989 to decompose the effects of space and time
# on R0.

################################################################################
# functions to make the survival growth kernels and the Fecundity kernels
# add large Upper size class to fixe the eviction and extrapolation problems.
make_P = function(sur_int, sur_slope, sur_G, cent_sur, grow_int,
	grow_slope, grow_G, grow_sigma, Z, U){

	# make the growth-survival kernel, matrix of |Z| by |Z|
	# calulualte surival for each size in Z
	sur_lp = sur_int + sur_slope * (pmin(U, Z) - cent_sur) + sur_G
  	S = 1 / (1 + exp(-sur_lp))

	# do the growth
	mean_growth = grow_int + grow_slope * Z + grow_G
	mu = expand.grid(Z, pmin(U, mean_growth))
	G = matrix(dnorm(mu[, 1], mu[, 2], grow_sigma),
	    nrow = length(Z), byrow = TRUE)

	# make growth conditional on survival, transpose to make P from col to row
	# so the multiplication with Fec works
	return(t(G * S))

}

make_Fec = function(rep_int, rep_slope, rep_G, cent_rep,
	int_f, slope_f, cent_f, Z, U){

	# make the fecundity kernel, vector of seeds produced by each size in Z
	# prob of flowering
    lp = rep_int + rep_slope * (pmin(U, Z) - cent_rep) + rep_G
    prob_flower = 1 / (1 + exp(-lp))

	# size effect on seed production by the probability of flowering
	return(matrix(fec_mode(int_f, slope_f, cent_f, Z, U) * prob_flower, nrow = 1))

}

################################################################################

fec_mode = function(int_f, beta_f, cent_const, Z, U){

  return(exp(int_f + beta_f * (pmin(U, Z) - cent_const)))

}

# Define two functions to calculate R0, both asymtpoic and over a fixed time,
# taking a survial growth kernel P and fecundity kernel Fec. Kernels are
# transition from size col -> to size row
R0_asy = function(P, Fec, n0){

	# construct the identiy matrix
	I = diag(length(Z))

	# Make the Neumann operator
	N = solve(I - P)

	return(sum((Fec %*% N) * n0))

}

R0_T = function(P, Fec, Th, n0){

	# iterate through years
	Pa = P
	R0_t = sum((Fec %*% Pa) * n0)

	if(Th > 1){

		for(t in 2:Th){

			Pa = Pa %*% P
			R0_t = R0_t + sum((Fec %*% Pa) * n0)

		}

	}

	return(R0_t)

}

################################################################################
# I also need the sesntivities to get effect of each parameter
# Use the finite difference method to find sensitivity to intercept in each
# vital rate (IPM book pp. 219)
# Error for the difference
delta = .Machine$double.eps ^ (1/3)

R0_sens_ays = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	Z, U, n0, Delta = delta, vital_rate){

	if(!vital_rate %in% c('rep', 'sur', 'grow', 'grow_s')){
		stop("vital_rate needs to be one of [rep, sur, grow, grow_s]")
	}

	if(vital_rate == 'rep'){

		rep_int_up = rep_int + Delta
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int - Delta
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'sur'){

		rep_int_up = rep_int
		sur_int_up = sur_int + Delta
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int - Delta
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int + Delta
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int - Delta
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow_s'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope + Delta

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope - Delta

	}

	P_up = make_P(sur_int_up, sur_slope, sur_G, cent_sur, grow_int_up, gr_slope_up,
		grow_G, grow_sigma, Z, U)
	P_down = make_P(sur_int_down, sur_slope, sur_G, cent_sur, grow_int_down, gr_slope_down,
		grow_G, grow_sigma, Z, U)

	Fec_up = make_Fec(rep_int_up, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)
	Fec_down = make_Fec(rep_int_down, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)


	R0_up = R0_asy(P_up, Fec_up, n0)

	R0_down = R0_asy(P_down, Fec_down, n0)

	return((R0_up - R0_down) / (2 * Delta))

}

# version using time limted R0
R0_sens_T = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	Th, Z, U, n0, Delta = delta, vital_rate){

	sense = numeric(length(rep_int))

	if(!vital_rate %in% c('rep', 'sur', 'grow', 'grow_s')){
		stop("vital_rate needs to be one of [rep, sur, grow, grow_s]")
	}

	if(vital_rate == 'rep'){

		rep_int_up = rep_int + Delta
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int - Delta
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'sur'){

		rep_int_up = rep_int
		sur_int_up = sur_int + Delta
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int - Delta
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int + Delta
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int - Delta
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow_s'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope + Delta

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope - Delta

	}


	P_up = make_P(sur_int_up, sur_slope, sur_G, cent_sur, grow_int_up, gr_slope_up,
		grow_G, grow_sigma, Z, U)
	P_down = make_P(sur_int_down, sur_slope, sur_G, cent_sur, grow_int_down, gr_slope_down,
		grow_G, grow_sigma, Z, U)

	Fec_up = make_Fec(rep_int_up, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)
	Fec_down = make_Fec(rep_int_down, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)

	R0_up = R0_T(P_up, Fec_up, Th, n0)

	R0_down = R0_T(P_down, Fec_down, Th, n0)

	return((R0_up - R0_down) / (2 * Delta))

}


R0_sens_ays_ln = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	Z, U, n0, Delta = delta, vital_rate){

	if(!vital_rate %in% c('rep', 'sur', 'grow', 'grow_s')){
		stop("vital_rate needs to be one of [rep, sur, grow, grow_s]")
	}

	if(vital_rate == 'rep'){

		rep_int_up = rep_int + Delta
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int - Delta
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'sur'){

		rep_int_up = rep_int
		sur_int_up = sur_int + Delta
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int - Delta
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int + Delta
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int - Delta
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow_s'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope + Delta

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope - Delta

	}

	P_up = make_P(sur_int_up, sur_slope, sur_G, cent_sur, grow_int_up, gr_slope_up,
		grow_G, grow_sigma, Z, U)
	P_down = make_P(sur_int_down, sur_slope, sur_G, cent_sur, grow_int_down, gr_slope_down,
		grow_G, grow_sigma, Z, U)

	Fec_up = make_Fec(rep_int_up, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)
	Fec_down = make_Fec(rep_int_down, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)


	R0_up = log(R0_asy(P_up, Fec_up, n0))

	R0_down = log(R0_asy(P_down, Fec_down, n0))

	return((R0_up - R0_down) / (2 * Delta))

}

# version using time limted R0
R0_sens_T_ln = function(rep_int, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f,
	sur_int, sur_slope, sur_G, cent_sur, grow_int, grow_slope, grow_G, grow_sigma,
	Th, Z, U, n0, Delta = delta, vital_rate){

	sense = numeric(length(rep_int))

	if(!vital_rate %in% c('rep', 'sur', 'grow', 'grow_s')){
		stop("vital_rate needs to be one of [rep, sur, grow, grow_s]")
	}

	if(vital_rate == 'rep'){

		rep_int_up = rep_int + Delta
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int - Delta
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'sur'){

		rep_int_up = rep_int
		sur_int_up = sur_int + Delta
		grow_int_up = grow_int
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int - Delta
		grow_int_down = grow_int
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int + Delta
		gr_slope_up = grow_slope

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int - Delta
		gr_slope_down = grow_slope

	}

	if(vital_rate == 'grow_s'){

		rep_int_up = rep_int
		sur_int_up = sur_int
		grow_int_up = grow_int
		gr_slope_up = grow_slope + Delta

		rep_int_down = rep_int
		sur_int_down = sur_int
		grow_int_down = grow_int
		gr_slope_down = grow_slope - Delta

	}


	P_up = make_P(sur_int_up, sur_slope, sur_G, cent_sur, grow_int_up, gr_slope_up,
		grow_G, grow_sigma, Z, U)
	P_down = make_P(sur_int_down, sur_slope, sur_G, cent_sur, grow_int_down, gr_slope_down,
		grow_G, grow_sigma, Z, U)

	Fec_up = make_Fec(rep_int_up, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)
	Fec_down = make_Fec(rep_int_down, rep_slope, rep_G, cent_rep, int_f, slope_f, cent_f, Z, U)

	R0_up = log(R0_T(P_up, Fec_up, Th, n0))

	R0_down = log(R0_T(P_down, Fec_down, Th, n0))

	return((R0_up - R0_down) / (2 * Delta))

}

#################################################################################
# LTRE
# use Th as a indicator, if Inf use the asymtopic function otherwise it is time horizon

LTRE_space_time = function(vr_ob, cent_rep, cent_f, cent_sur, Z, U, n0,
	Th = Inf, Delta = delta){

	# I need some arrays to hold all the kernels for each location/year, will be
	# 4D first 2 dims are the from z to z matricies, and the next 2 are space and time

	# 1: parameter values for each year and location
	param_array = array(NA,
		dim = c(length(vr_ob$grow_G), length(vr_ob$grow_year), 4),
		dimnames = list(L = NULL, year = NULL, pr = c('rep_int', 'sur_int', 'gr_int', 'gr_slope')))

	# fill the array
	for(l in 1:dim(param_array)[1]){
		for(y in 1:dim(param_array)[2]){

			param_array[l, y, 'rep_int'] = vr_ob$rep_year[y] + vr_ob$rep_G[l]
			param_array[l, y, 'sur_int'] = vr_ob$sur_year[y] + vr_ob$sur_G[l]
			param_array[l, y, 'gr_int'] = vr_ob$grow_year[y] + vr_ob$grow_G[l]
			param_array[l, y, 'gr_slope'] = vr_ob$grow_slope[y]
		}
	}

	# 2. averaged parameters
	rep_int_ave = mean(param_array[, , 'rep_int'])
	sur_int_ave = mean(param_array[, , 'sur_int'])
	gr_int_ave = mean(param_array[, , 'gr_int'])
	gr_slope_ave = mean(param_array[, , 'gr_slope'])

	# 3. K under parameters averaged over locations and year
	P_ave = make_P(sur_int = sur_int_ave, sur_slope = vr_ob$sur_slope, sur_G = 0, # rep_G here is subsummed in to rep_int[l,y]
		cent_sur = cent_sur, grow_int = gr_int_ave, grow_slope = gr_slope_ave, grow_G = 0,
		grow_sigma = vr_ob$grow_sigma, Z, U)

	F_ave = make_Fec(rep_int = rep_int_ave, rep_slope = vr_ob$rep_slope, rep_G = 0, # rep_G here is subsummed in to rep_int[l,y]
		cent_rep = cent_rep, int_f = vr_ob$int_f, slope_f = vr_ob$slope_f,
		cent_f = cent_f, Z = Z, U)

	# R0 | mean kernel, used for comparison
	R0_ave = ifelse(Th == Inf, R0_asy(P_ave, F_ave, n0), R0_T(P_ave, F_ave, Th, n0))

	# get R0_lt and the means over time and location and senesitivites to parameters
	R0 = param_array[, , 1]
	R0[,] = NA
	sen_rep = param_array[, , 1]
	sen_rep[,] = NA
	sen_sur = param_array[, , 1]
	sen_sur[,] = NA
	sen_gr = param_array[, , 1]
	sen_gr[,] = NA
	sen_gr_s = param_array[, , 1]
	sen_gr_s[,] = NA

	# time averaged intercepts
	sen_rep_l = param_array[, 1, 1]
	sen_rep_l[] = NA
	sen_sur_l = param_array[, 1, 1]
	sen_sur_l[] = NA
	sen_gr_l = param_array[, 1, 1]
	sen_gr_l[] = NA

	rep_int_l = c()
	sur_int_l = c()
	gr_int_l = c()

	# location averaged intercepts
	sen_rep_t = param_array[1, , 1]
	sen_rep_t[] = NA
	sen_sur_t = param_array[1, , 1]
	sen_sur_t[] = NA
	sen_gr_t = param_array[1, , 1]
	sen_gr_t = NA
	sen_gr_s_t = param_array[1, , 1]
	sen_gr_s_t[] = NA

	rep_int_t = c()
	sur_int_t = c()
	gr_int_t = c()
	gr_slope_t = c()

	for(t in 1:length(sen_rep_t)){
		rep_int_t[t] = mean(param_array[, t, 'rep_int'])
		sur_int_t[t] = mean(param_array[, t, 'sur_int'])
		gr_int_t[t] = mean(param_array[, t, 'gr_int'])
		gr_slope_t[t] = mean(param_array[, t, 'gr_slope'])

		# Get the location averaged sensitivites
		sen_rep_t[t] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'rep'),
			R0_sens_T((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'rep'))

		sen_sur_t[t] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'sur'),
			R0_sens_T((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'sur'))

		sen_gr_t[t] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'grow'),
			R0_sens_T((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'grow'))

		sen_gr_s_t[t] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'grow_s'),
			R0_sens_T((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'grow_s'))

	}

	# fill each slice of this array with the Kernels
	for(l in 1:dim(param_array)[1]){

		# take the average sensitvites over time
		# start by take means across time of the paramerters
		rep_int_l[l] = mean(param_array[l, , 'rep_int'])
		sur_int_l[l] = mean(param_array[l, , 'sur_int'])
		gr_int_l[l] = mean(param_array[l, , 'gr_int'])
		gr_slope_l = mean(param_array[l, , 'gr_slope']) # slope does not vary over location

		# find the sensitivies to parameter at mid point between overall average
		# and time average
		# Sensitivity to R0 under the kernel built at mid-point of parameters
		sen_rep_l[l] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'rep'),
			R0_sens_T((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'rep'))

		# Get the time averaged sensitivites
		sen_sur_l[l] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'sur'),
			R0_sens_T((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'sur'))

		sen_gr_l[l] = ifelse(Th == Inf,
			R0_sens_ays((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'grow'),
			R0_sens_T((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'grow'))

		for(y in 1:dim(param_array)[2]){

			P_lt = make_P(sur_int = param_array[l, y, 'sur_int'],
				sur_slope = vr_ob$sur_slope, sur_G = 0, # rep_G here is subsummed in to rep_int[l,y]
				cent_sur = cent_sur, grow_int = param_array[l, y, 'gr_int'],
				grow_slope = param_array[l, y, 'gr_slope'], grow_G = 0,
				grow_sigma = vr_ob$grow_sigma, Z, U)

			F_lt = make_Fec(rep_int = param_array[l, y, 'rep_int'],
				rep_slope = vr_ob$rep_slope, rep_G = 0, # rep_G here is subsummed in to rep_int[l,y]
				cent_rep = cent_rep, int_f = vr_ob$int_f, slope_f = vr_ob$slope_f,
				cent_f = cent_f, Z = Z, U)

			R0[l, y] = ifelse(Th == Inf, R0_asy(P_lt, F_lt, n0),
				R0_T(P_lt, F_lt, Th, n0))

			# get the interaction senesitivity
			sen_rep[l, y] = ifelse(Th == Inf,
				R0_sens_ays((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'rep'),
				R0_sens_T((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'rep'))

			sen_sur[l, y] = ifelse(Th == Inf,
				R0_sens_ays((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'sur'),
				R0_sens_T((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'sur'))

			sen_gr[l, y] = ifelse(Th == Inf,
				R0_sens_ays((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'grow'),
				R0_sens_T((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'grow'))

			sen_gr_s[l, y] = ifelse(Th == Inf,
				R0_sens_ays((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'grow_s'),
				R0_sens_T((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'grow_s'))

		}
	}

	# get all the contributions and wrap it all up into a dataframe
	df = make_cont_df(param_array, rep_int_t, rep_int_l, rep_int_ave,
		sur_int_t, sur_int_l, sur_int_ave, gr_int_t, gr_int_l, gr_int_ave,
		gr_slope_t, gr_slope_l, gr_slope_ave, R0, sen_rep_l, sen_sur_l, sen_gr_l,
		sen_rep_t, sen_sur_t, sen_gr_t, sen_gr_s_t, sen_rep, sen_sur,
		sen_gr, sen_gr_s, R0_ave, loc_IDs = names(vr_ob$grow_G))

	return(df)

}

LTRE_space_time_ln = function(vr_ob, cent_rep, cent_f, cent_sur, Z, U, n0,
	Th = Inf, Delta = delta){

	# I need some arrays to hold all the kernels for each location/year, will be
	# 4D first 2 dims are the from z to z matricies, and the next 2 are space and time

	# 1: parameter values for each year and location
	param_array = array(NA,
		dim = c(length(vr_ob$grow_G), length(vr_ob$grow_year), 4),
		dimnames = list(L = NULL, year = NULL, pr = c('rep_int', 'sur_int', 'gr_int', 'gr_slope')))

	# fill the array
	for(l in 1:dim(param_array)[1]){
		for(y in 1:dim(param_array)[2]){

			param_array[l, y, 'rep_int'] = vr_ob$rep_year[y] + vr_ob$rep_G[l]
			param_array[l, y, 'sur_int'] = vr_ob$sur_year[y] + vr_ob$sur_G[l]
			param_array[l, y, 'gr_int'] = vr_ob$grow_year[y] + vr_ob$grow_G[l]
			param_array[l, y, 'gr_slope'] = vr_ob$grow_slope[y]
		}
	}

	# 2. averaged parameters
	rep_int_ave = mean(param_array[, , 'rep_int'])
	sur_int_ave = mean(param_array[, , 'sur_int'])
	gr_int_ave = mean(param_array[, , 'gr_int'])
	gr_slope_ave = mean(param_array[, , 'gr_slope'])

	# 3. K under parameters averaged over locations and year
	P_ave = make_P(sur_int = sur_int_ave, sur_slope = vr_ob$sur_slope, sur_G = 0, # rep_G here is subsummed in to rep_int[l,y]
		cent_sur = cent_sur, grow_int = gr_int_ave, grow_slope = gr_slope_ave, grow_G = 0,
		grow_sigma = vr_ob$grow_sigma, Z, U)

	F_ave = make_Fec(rep_int = rep_int_ave, rep_slope = vr_ob$rep_slope, rep_G = 0, # rep_G here is subsummed in to rep_int[l,y]
		cent_rep = cent_rep, int_f = vr_ob$int_f, slope_f = vr_ob$slope_f,
		cent_f = cent_f, Z = Z, U)

	# R0 | mean kernel, used for comparison
	R0_ave = log(ifelse(Th == Inf, R0_asy(P_ave, F_ave, n0), R0_T(P_ave, F_ave, Th, n0)))

	# get R0_lt and the means over time and location and senesitivites to parameters
	R0 = param_array[, , 1]
	R0[,] = NA
	sen_rep = param_array[, , 1]
	sen_rep[,] = NA
	sen_sur = param_array[, , 1]
	sen_sur[,] = NA
	sen_gr = param_array[, , 1]
	sen_gr[,] = NA
	sen_gr_s = param_array[, , 1]
	sen_gr_s[,] = NA

	# time averaged intercepts
	sen_rep_l = param_array[, 1, 1]
	sen_rep_l[] = NA
	sen_sur_l = param_array[, 1, 1]
	sen_sur_l[] = NA
	sen_gr_l = param_array[, 1, 1]
	sen_gr_l[] = NA

	rep_int_l = c()
	sur_int_l = c()
	gr_int_l = c()

	# location averaged intercepts
	sen_rep_t = param_array[1, , 1]
	sen_rep_t[] = NA
	sen_sur_t = param_array[1, , 1]
	sen_sur_t[] = NA
	sen_gr_t = param_array[1, , 1]
	sen_gr_t = NA
	sen_gr_s_t = param_array[1, , 1]
	sen_gr_s_t[] = NA

	rep_int_t = c()
	sur_int_t = c()
	gr_int_t = c()
	gr_slope_t = c()

	for(t in 1:length(sen_rep_t)){
		rep_int_t[t] = mean(param_array[, t, 'rep_int'])
		sur_int_t[t] = mean(param_array[, t, 'sur_int'])
		gr_int_t[t] = mean(param_array[, t, 'gr_int'])
		gr_slope_t[t] = mean(param_array[, t, 'gr_slope'])

		# Get the location averaged sensitivites
		sen_rep_t[t] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'rep'),
			R0_sens_T_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'rep'))

		sen_sur_t[t] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'sur'),
			R0_sens_T_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'sur'))

		sen_gr_t[t] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'grow'),
			R0_sens_T_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'grow'))

		sen_gr_s_t[t] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'grow_s'),
			R0_sens_T_ln((rep_int_ave + rep_int_t[t]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_t[t]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_t[t]) / 2, (gr_slope_ave + gr_slope_t[t]) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'grow_s'))

	}

	# fill each slice of this array with the Kernels
	for(l in 1:dim(param_array)[1]){

		# take the average sensitvites over time
		# start by take means across time of the paramerters
		rep_int_l[l] = mean(param_array[l, , 'rep_int'])
		sur_int_l[l] = mean(param_array[l, , 'sur_int'])
		gr_int_l[l] = mean(param_array[l, , 'gr_int'])
		gr_slope_l = mean(param_array[l, , 'gr_slope']) # slope does not vary over location

		# find the sensitivies to parameter at mid point between overall average
		# and time average
		# Sensitivity to R0 under the kernel built at mid-point of parameters
		sen_rep_l[l] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'rep'),
			R0_sens_T_ln((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'rep'))

		# Get the time averaged sensitivites
		sen_sur_l[l] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'sur'),
			R0_sens_T_ln((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'sur'))

		sen_gr_l[l] = ifelse(Th == Inf,
			R0_sens_ays_ln((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Z, U, n0, Delta = delta, 'grow'),
			R0_sens_T_ln((rep_int_ave + rep_int_l[l]) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
				vr_ob$int_f, vr_ob$slope_f, cent_f,
				(sur_int_ave + sur_int_l[l]) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
				(gr_int_ave + gr_int_l[l]) / 2, (gr_slope_ave + gr_slope_l) / 2, grow_G = 0, vr_ob$grow_sigma,
				Th, Z, U, n0, Delta = delta, 'grow'))

		for(y in 1:dim(param_array)[2]){

			P_lt = make_P(sur_int = param_array[l, y, 'sur_int'],
				sur_slope = vr_ob$sur_slope, sur_G = 0, # rep_G here is subsummed in to rep_int[l,y]
				cent_sur = cent_sur, grow_int = param_array[l, y, 'gr_int'],
				grow_slope = param_array[l, y, 'gr_slope'], grow_G = 0,
				grow_sigma = vr_ob$grow_sigma, Z, U)

			F_lt = make_Fec(rep_int = param_array[l, y, 'rep_int'],
				rep_slope = vr_ob$rep_slope, rep_G = 0, # rep_G here is subsummed in to rep_int[l,y]
				cent_rep = cent_rep, int_f = vr_ob$int_f, slope_f = vr_ob$slope_f,
				cent_f = cent_f, Z = Z, U)

			R0[l, y] = log(ifelse(Th == Inf, R0_asy(P_lt, F_lt, n0),
				R0_T(P_lt, F_lt, Th, n0)))

			# get the interaction senesitivity
			sen_rep[l, y] = ifelse(Th == Inf,
				R0_sens_ays_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'rep'),
				R0_sens_T_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'rep'))

			sen_sur[l, y] = ifelse(Th == Inf,
				R0_sens_ays_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'sur'),
				R0_sens_T_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'sur'))

			sen_gr[l, y] = ifelse(Th == Inf,
				R0_sens_ays_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'grow'),
				R0_sens_T_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'grow'))

			sen_gr_s[l, y] = ifelse(Th == Inf,
				R0_sens_ays_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Z, U, n0, Delta = delta, 'grow_s'),
				R0_sens_T_ln((rep_int_ave + param_array[l, y, 'rep_int']) / 2 , vr_ob$rep_slope, rep_G = 0, cent_rep,
					vr_ob$int_f, vr_ob$slope_f, cent_f,
					(sur_int_ave + param_array[l, y, 'sur_int']) / 2, vr_ob$sur_slope, sur_G = 0, cent_sur,
					(gr_int_ave + param_array[l, y, 'gr_int']) / 2,
					(gr_slope_ave + param_array[l, y, 'gr_slope']) / 2, grow_G = 0, vr_ob$grow_sigma,
					Th, Z, U, n0, Delta = delta, 'grow_s'))

		}
	}

	# get all the contributions and wrap it all up into a dataframe
	df = make_cont_df(param_array, rep_int_t, rep_int_l, rep_int_ave,
		sur_int_t, sur_int_l, sur_int_ave, gr_int_t, gr_int_l, gr_int_ave,
		gr_slope_t, gr_slope_l, gr_slope_ave, R0, sen_rep_l, sen_sur_l, sen_gr_l,
		sen_rep_t, sen_sur_t, sen_gr_t, sen_gr_s_t, sen_rep, sen_sur,
		sen_gr, sen_gr_s, R0_ave, loc_IDs = names(vr_ob$grow_G))

	return(df)

}

# make  list of objects to pass to the LTRE function, to allow easy mclapply work
make_par_list = function(rep_samps, sur_samps, gr_samps, fruit_samps){

	# make the template to fill
	pr_temp = list(rep_year = c(NA, NA, NA, NA, NA),
		rep_slope = NA, rep_G = numeric(length(gr_samps$SPP[, 'uni_spp'])),
		int_f = NA, slope_f = NA,
		sur_year = c(NA, NA, NA, NA, NA),
		sur_slope = NA, sur_G = numeric(length(gr_samps$SPP[, 'uni_spp'])),
		grow_year = c(NA, NA, NA, NA, NA),
		grow_slope = c(NA, NA, NA, NA, NA),
		grow_G = numeric(length(gr_samps$SPP[, 'uni_spp'])), grow_sigma = NA)

	names(pr_temp$rep_G) = gr_samps$SPP[, 'uni_spp']
	names(pr_temp$sur_G) = gr_samps$SPP[, 'uni_spp']
	names(pr_temp$grow_G) = gr_samps$SPP[, 'uni_spp']

	names(pr_temp$rep_year) = c('y03', 'y04', 'y05', 'y06', 'y07')
	names(pr_temp$sur_year) = c('y03', 'y04', 'y05', 'y06', 'y07')
	names(pr_temp$grow_year) = c('y03', 'y04', 'y05', 'y06', 'y07')
	names(pr_temp$grow_slope) = c('y03', 'y04', 'y05', 'y06', 'y07')

	pr_list = list()

	for(i in seq_along(gr_samps$sigma)){

		pr_list[[i]] = pr_temp

		# fill the fields
		pr_list[[i]]$rep_year[1:5] = as.numeric(rep_samps$year_int[i, 2:6])
		pr_list[[i]]$rep_slope = rep_samps$slope[i]
		pr_list[[i]]$rep_G[seq_along(pr_temp$rep_G)] = rep_samps$SPP[, 3 + i]

		pr_list[[i]]$sur_year[1:5] = as.numeric(sur_samps$year_int[i, 1:5])
		pr_list[[i]]$sur_slope = sur_samps$slope[i]
		pr_list[[i]]$sur_G[seq_along(pr_temp$sur_G)] = sur_samps$SPP[, 3 + i]

		pr_list[[i]]$grow_year[1:5] = as.numeric(gr_samps$year_int[i, 1:5])
		pr_list[[i]]$grow_slope[1:5] = as.numeric(gr_samps$year_slope[i, 1:5])
		pr_list[[i]]$grow_G[seq_along(pr_temp$grow_G)] = gr_samps$SPP[, 3 + i]
		pr_list[[i]]$grow_sigma = gr_samps$sigma[i]

		pr_list[[i]]$int_f = fruit_samps$int[i]
		pr_list[[i]]$slope_f = fruit_samps$slope[i]

		pr_list[[i]]$samp_ID = i

	}

	return(pr_list)

}

# make a parameter object with mean values for each value
make_par_ave = function(rep_samps, sur_samps, gr_samps, fruit_samps){

	# make the template to fill
	pr_temp = list(rep_year = c(NA, NA, NA, NA, NA),
		rep_slope = NA, rep_G = numeric(length(gr_samps$SPP[, 'uni_spp'])),
		int_f = NA, slope_f = NA,
		sur_year = c(NA, NA, NA, NA, NA),
		sur_slope = NA, sur_G = numeric(length(gr_samps$SPP[, 'uni_spp'])),
		grow_year = c(NA, NA, NA, NA, NA),
		grow_slope = c(NA, NA, NA, NA, NA),
		grow_G = numeric(length(gr_samps$SPP[, 'uni_spp'])), grow_sigma = NA)

	names(pr_temp$rep_G) = gr_samps$SPP[, 'uni_spp']
	names(pr_temp$sur_G) = gr_samps$SPP[, 'uni_spp']
	names(pr_temp$grow_G) = gr_samps$SPP[, 'uni_spp']

	names(pr_temp$rep_year) = c('y03', 'y04', 'y05', 'y06', 'y07')
	names(pr_temp$sur_year) = c('y03', 'y04', 'y05', 'y06', 'y07')
	names(pr_temp$grow_year) = c('y03', 'y04', 'y05', 'y06', 'y07')
	names(pr_temp$grow_slope) = c('y03', 'y04', 'y05', 'y06', 'y07')

	pr_temp$rep_year[1:5] = as.numeric(apply(rep_samps$year_int[, 2:6], MARGIN = 2, FUN = mean))
	pr_temp$rep_slope = mean(rep_samps$slope)
	pr_temp$rep_G[seq_along(pr_temp$rep_G)] = apply(select(rep_samps$SPP, s1:s1000), MARGIN = 1, FUN = mean)

	pr_temp$sur_year[1:5] = as.numeric(apply(sur_samps$year_int[, 1:5], MARGIN = 2, FUN = mean))
	pr_temp$sur_slope = mean(sur_samps$slope)
	pr_temp$sur_G[seq_along(pr_temp$sur_G)] = apply(select(sur_samps$SPP, s1:s1000), MARGIN = 1, FUN = mean)

	pr_temp$grow_year[1:5] = as.numeric(apply(gr_samps$year_int[, 1:5], MARGIN = 2, FUN = mean))
	pr_temp$grow_slope = as.numeric(apply(gr_samps$year_slope[, 1:5], MARGIN = 2, FUN = mean))
	pr_temp$grow_G[seq_along(pr_temp$grow_G)] = apply(select(gr_samps$SPP, s1:s1000), MARGIN = 1, FUN = mean)
	pr_temp$grow_sigma = mean(gr_samps$sigma)

	pr_temp$int_f = mean(fruit_samps$int)
	pr_temp$slope_f = mean(fruit_samps$slope)

	return(pr_temp)

}

# wrape all the contributions into a dataframe
make_cont_df = function(param_array, rep_int_t, rep_int_l, rep_int_ave,
	sur_int_t, sur_int_l, sur_int_ave, gr_int_t, gr_int_l, gr_int_ave,
	gr_slope_t, gr_slope_l, gr_slope_ave, R0_lt, sen_rep_l, sen_sur_l, sen_gr_l,
	sen_rep_t, sen_sur_t, sen_gr_t, sen_gr_s_t, sen_rep_lt, sen_sur_lt,
	sen_gr_lt, sen_gr_s_lt, R0_ave, loc_IDs){

	# start by making a dataframe of all the parameters and R0s
	rep_R0 = param_array[, , 'rep_int']
	colnames(rep_R0) = paste0('y0', 3:7)
	rep_R0_df = as_data_frame(rep_R0) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'rep_int', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	sur_R0 = param_array[, , 'sur_int']
	colnames(sur_R0) = paste0('y0', 3:7)
	sur_R0_df = as_data_frame(sur_R0) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'sur_int', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	gr_R0 = param_array[, , 'gr_int']
	colnames(gr_R0) = paste0('y0', 3:7)
	gr_R0_df = as_data_frame(gr_R0) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'gr_int', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	gr_s_R0 = param_array[, , 'gr_slope']
	colnames(gr_s_R0) = paste0('y0', 3:7)
	gr_s_R0_df = as_data_frame(gr_s_R0) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'gr_slope', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	R0_mat = R0_lt
	colnames(R0_mat) = paste0('y0', 3:7)
	R0_df = as_data_frame(R0_mat) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'R0', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	#join these all up
	out_df = inner_join(rep_R0_df, select(sur_R0_df, row_ID, sur_int), by = 'row_ID') %>%
		inner_join(select(gr_R0_df, row_ID, gr_int), by = 'row_ID') %>%
		inner_join(select(gr_s_R0_df, row_ID, gr_slope), by = 'row_ID') %>%
		inner_join(select(R0_df, row_ID, R0), by = 'row_ID')

	# join the sesnitivites for each vital rate at each location and then each time
	out_df = inner_join(out_df,
				data.frame(loc_ID = loc_IDs, sen_rep_L = sen_rep_l, sen_sur_L = sen_sur_l,
					sen_gr_L = sen_gr_l, sen_gr_s_L = 0), by = 'loc_ID') %>%
			inner_join(data.frame(year = paste0('y0', 3:7), sen_rep_T = sen_rep_t, sen_sur_T = sen_sur_t,
					sen_gr_T = sen_gr_t, sen_gr_s_T = sen_gr_s_t), by = 'year')

	# join the interaction senitivites to the data frame
	sen_rep_mat = sen_rep_lt
	colnames(sen_rep_mat) = paste0('y0', 3:7)
	sen_rep_df = as_data_frame(sen_rep_mat) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'sen_rep_LT', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	sen_sur_mat = sen_sur_lt
	colnames(sen_sur_mat) = paste0('y0', 3:7)
	sen_sur_df = as_data_frame(sen_sur_mat) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'sen_sur_LT', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	sen_gr_mat = sen_gr_lt
	colnames(sen_gr_mat) = paste0('y0', 3:7)
	sen_gr_df = as_data_frame(sen_gr_mat) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'sen_gr_LT', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	sen_gr_s_mat = sen_gr_s_lt
	colnames(sen_gr_s_mat) = paste0('y0', 3:7)
	sen_gr_s_df = as_data_frame(sen_gr_s_mat) %>% mutate(loc_ID = loc_IDs) %>%
		gather(key = 'year', value = 'sen_gr_s_LT', y03:y07) %>%
		mutate(row_ID = paste0(loc_ID, ':', year))

	out_df = inner_join(out_df, select(sen_rep_df, row_ID, sen_rep_LT), by = 'row_ID') %>%
		inner_join(select(sen_sur_df, row_ID, sen_sur_LT), by = 'row_ID') %>%
		inner_join(select(sen_gr_df, row_ID, sen_gr_LT), by = 'row_ID') %>%
		inner_join(select(sen_gr_s_df, row_ID, sen_gr_s_LT), by = 'row_ID')

	# get the averaged parameter values to work out the contributions
	out_df = inner_join(out_df, data.frame(year = paste0('y0', 3:7), rep_int_T = rep_int_t), by = 'year') %>%
		inner_join(data.frame(loc_ID = loc_IDs, rep_int_L = rep_int_l), by = 'loc_ID') %>%
		mutate(rep_int_A = rep_int_ave) %>%
		inner_join(data.frame(year = paste0('y0', 3:7), sur_int_T = sur_int_t), by = 'year') %>%
		inner_join(data.frame(loc_ID = loc_IDs, sur_int_L = sur_int_l), by = 'loc_ID') %>%
		mutate(sur_int_A = sur_int_ave) %>%
		inner_join(data.frame(year = paste0('y0', 3:7), gr_int_T = gr_int_t), by = 'year') %>%
		inner_join(data.frame(loc_ID = loc_IDs, gr_int_L = gr_int_l), by = 'loc_ID') %>%
		mutate(gr_int_A = gr_int_ave) %>%
		inner_join(data.frame(year = paste0('y0', 3:7), gr_slope_T = gr_slope_t), by = 'year') %>%
		inner_join(data.frame(loc_ID = loc_IDs, gr_slope_L = gr_slope_l), by = 'loc_ID') %>%
		mutate(gr_slope_A = gr_slope_ave)

	# now make the actual contributions
	out_df = mutate(out_df, cont_rep_T = (rep_int_T - rep_int_A) * sen_rep_T,
			cont_sur_T = (sur_int_T - sur_int_A) * sen_sur_T,
			cont_gr_T = (gr_int_T - gr_int_A) * sen_gr_T,
			cont_gr_s_T = (gr_slope_T - gr_slope_A) * sen_gr_s_T,
			cont_rep_L = (rep_int_L - rep_int_A) * sen_rep_L,
			cont_sur_L = (sur_int_L - sur_int_A) * sen_sur_L,
			cont_gr_L = (gr_int_L - gr_int_A) * sen_gr_L,
			cont_gr_s_L = (gr_slope_L - gr_slope_A) * sen_gr_s_L,
			cont_rep_LT = ((rep_int - rep_int_A) * sen_rep_LT) - cont_rep_L - cont_rep_T,
			cont_sur_LT = ((sur_int - sur_int_A) * sen_sur_LT) - cont_sur_L - cont_sur_T,
			cont_gr_LT = ((gr_int - gr_int_A) * sen_gr_LT) - cont_gr_L - cont_gr_T,
			cont_gr_s_LT = ((gr_slope - gr_slope_A) * sen_gr_s_LT) - cont_gr_s_L - cont_gr_s_T,
			cont_L = cont_rep_L + cont_sur_L + cont_gr_L + cont_gr_s_L,
			cont_T = cont_rep_T + cont_sur_T + cont_gr_T + cont_gr_s_T,
			cont_LT = cont_rep_LT + cont_sur_LT + cont_gr_LT + cont_gr_s_LT,
			R0_ave = R0_ave, R0_pred = R0_ave + cont_L + cont_T + cont_LT)

	return(out_df)

}

################################################################################
# Make a data frame of changes in each vital rate and R0
make_R0_df = function(rep_year, rep_slope, rep_loc, rep_cent, int_f, slope_f, cent_f,
	sur_year, sur_slope, sur_loc, sur_cent, grow_year, grow_slope, grow_loc, grow_sigma,
	Z, n0, Th = Inf, loc_IDs){

	time_df = data.frame(year = c('y03', 'y04', 'y05', 'y06', 'y07'),
		rep_year = rep_year, rep_slope = rep_slope, rep_cent = rep_cent,
		fruit_int = int_f, fruit_slope = slope_f, cent_f, sur_year = sur_year,
		sur_slope = sur_slope, sur_cent = sur_cent, grow_year = grow_year,
		grow_slope = grow_slope)

	# repeate this for each location to make the mesh up
	time_list = list()
	for(l in seq_along(loc_IDs)){

		time_list[[l]] = time_df
		time_list[[l]]$loc_ID = loc_IDs[l]

	}

	time_df_spread = bind_rows(time_list)

	loc_df = data.frame(loc_ID = loc_IDs, rep_loc = rep_loc, sur_loc = sur_loc, gr_loc = sur_loc)

	df = inner_join(loc_df, time_df_spread, by = 'loc_ID') %>%
		mutate(R0 = NA, rep_int = rep_year + rep_loc,
			sur_int = sur_year + sur_loc, gr_int = grow_year + gr_loc)

	df$R0 = NA

	if(Th == Inf){
		for(i in seq_along(df$loc_ID)){

			P = make_P(df$sur_int[i], df$sur_slope[i], sur_G = 0, cent_sur = sur_cent,
					df$gr_int[i], df$grow_slope[i], grow_G = 0, grow_sigma, Z)
			Fec = make_Fec(df$rep_int[i], df$rep_slope[i], rep_G = 0, cent_rep = rep_cent,
				int_f, slope_f, cent_f, Z)

			df$R0[i] = R0_asy(P, Fec, n0)

		}
	}else{
		for(i in seq_along(df$loc_ID)){

			P = make_P(df$sur_int[i], df$sur_slope[i], sur_G = 0, cent_sur = sur_cent,
					df$gr_int[i], df$grow_slope[i], grow_G = 0, grow_sigma, Z)
			Fec = make_Fec(df$rep_int[i], df$rep_slope[i], rep_G = 0, cent_rep = rep_cent,
				int_f, slope_f, cent_f, Z)

			df$R0[i] = R0_T(P, Fec, Th, n0)

		}
	}

	return(df)

}

################################################################################
# cut down function to find R0 in the average year
R0_year_ave = function(vr_ob, cent_rep, cent_f, cent_sur, Z, U, n0, Th = Inf){

	print(vr_ob$samp_ID)
	# find mean across years for each vital rate
	mean_rep = mean(vr_ob$rep_year)
	mean_sur = mean(vr_ob$sur_year)
	mean_gr = mean(vr_ob$grow_year)
	mean_grs = mean(vr_ob$grow_slope)

	# get the unified location IDs for each estiamted R0
	loc_IDs = names(vr_ob$grow_G)
	R0_ave = numeric(length(vr_ob$grow_G))

	# for each location
	for(l in seq_along(vr_ob$grow_G)){

		P = make_P(sur_int = mean_sur, sur_slope = vr_ob$sur_slope, sur_G = vr_ob$sur_G[l],
			cent_sur = cent_sur, grow_int = mean_gr, grow_slope = mean_grs, grow_G = vr_ob$grow_G[l],
			grow_sigma = vr_ob$grow_sigma, Z, U)

		Fec = make_Fec(rep_int = mean_rep, rep_slope = vr_ob$rep_slope, rep_G = vr_ob$rep_G[l],
			cent_rep = cent_rep, int_f = vr_ob$int_f, slope_f = vr_ob$slope_f,
			cent_f = cent_f, Z = Z, U)

		# R0
		R0_ave[l] = ifelse(Th == Inf, R0_asy(P, Fec, n0), R0_T(P, Fec, Th, n0))

	}

	return(data.frame(samp_ID = vr_ob$samp_ID, loc_ID = loc_IDs, R0 = R0_ave))

}



#
