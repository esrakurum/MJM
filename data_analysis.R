
MJM_fit = function(data, xvar.names, zvar.names, yvar.name, timevar.name, lambda, nknots, initial.values, n.boot)
{
##Load the required packages	


	subject_id = data$subject_id
	facility_id = data$facility_id
#facility ids:

	fac_id = unique(facility_id)  

##Determine the number of facilities:

	n = length(fac_id) 


#Determine the number of subjects within each facility:
	n_i = rep(0, n)
	for(i in 1:n)
		n_i[i] = length(unique(subject_id[facility_id == fac_id[i]]))


##Determine the number of observations for each subject within each facility

	n_ij1 = list()

	for(i in 1:n)
		n_ij1[[i]] = count(subject_id[facility_id == fac_id[i]])$freq 

###Prepare the data for analysis: everything needs to be in data-frame format.


	T = list() ###A data frame to pick up the survival times for each facility, for instance T[[1]] has the survival times for all subjects within the first facility-

	X = list() # A data frame to pick up the predictors for the longitudinal outcome for each facility. For instance X[[i]] has the observation times for all subjects within the ith facility and in X[[i]], each subject is a row, each column is an observation.

	Z = list() # A data frame to pick up the predictors for the survival outcome for each facility. For instance Z[[i]] has the observation times for all subjects within the ith facility and in Z[[i]], each subject is a row, each column is an observation.

	Y = list()   ##A data frame to pick up the longitudinal outcome for each facility.

	delta_ij = list() ## A data frame to pick up survival status for each facility, 1 is the event happened, 0 is event did not happen.

###Subject-level predictors
	time = data[timevar.name]
	X_data = cbind(time, data[xvar.names])


###Facility-level predictors

	Z_data1 = data[zvar.names]
	Z_data = unique(cbind(subject_id, facility_id, Z_data1))
	
	
###Create the data frames:

	surv_data = data.frame(subject_id, survival = data$survival, time_surv = data$time_surv, facility_id)

	Y_data = data[yvar.name]
	
	for(i in 1:n)
	{
		T[[i]] = unique.matrix(surv_data[facility_id == fac_id[i], ])$time_surv
		X[[i]] = as.matrix(X_data[facility_id == fac_id[i], ])
		Z[[i]] = as.matrix(Z_data[Z_data[,2] == fac_id[i], 3:dim(Z_data)[2]]) 	

	 
		Y[[i]] = Y_data[facility_id == fac_id[i],]
		delta_ij[[i]] = unique.matrix(surv_data[facility_id == fac_id[i], ])$survival
	}

	facility_id = facility_id
	subject_id = subject_id

	
###Number of parameters

	p = dim(X_data)[2]  #dimension of beta, subject-specific covariates, excluding the intercept, zeta is p-1, since it does not include the time_long

	q = dim(Z_data)[2]-2 #dimension of eta and psi, facility-related covariates
  


#### Initial values for the parameters


### Longitudinal component:


	t_betas = initial.values[1:(p+1)] 
	t_psi = initial.values[((p+2):(p + q + 1))]


###Variances of the random effects:

 
	sigma_b_2 = initial.values [(p+q+p+q+2)]
	sigma_xi_2 = initial.values [(p+q+p+q+3)]



###Survival component: 

	t_zetas = initial.values[(p+q+2): (p+q+p)]
	t_etas = initial.values[(p+q+p+1): (p+q+p+q)]

	t_alpha = initial.values[(p+q+p+q+1)]

	####Load the necessary functions for the EM algorithm

	#source('functions.R', local = FALSE, chdir = "TRUE")


	####Function for the EM algorithm:

  	source('data_fitting.R', local = FALSE, chdir = "TRUE")
	

	environment(JM_est) <- environment()
	trt = JM_est()
################
## Estimated parameters and random effects
	final_parameters = trt$final_parameters
	likelihood_se = sqrt(trt$Variances)
	final_xis = trt$final_xi
	final_bs = trt$final_b

	names_par = cbind("Intercept", "Time")
	for (i in 1:(p-1))
		names_par <- cbind(names_par, paste0('X_', i))

	for (i in 1:q)
		names_par <- cbind(names_par, paste0('Z_', i))
		
	names_par = c(names_par, names_par[,3:dim(names_par)[2]], "alpha", "lambda", "Subject level", "Facility level")
	
	names(final_parameters) = names_par
	
	###Bootstrap variances:
	
	X_old = X
	Z_old = Z
	Y_old = Y

	T_old = T


	delta_old = delta_ij

	nis_old = n_i 
	nij_old = n_ij1 
	
	final_par_boot = list()
	
	if(n.boot!=0)
	{
		print("Bootstrap variance calculation...")

	
		for(l in 1:n.boot)
		{
		
			boot.ids = sample(1:n, n, replace = TRUE) #new sample's facility_ids

			new_X = list()
			new_Y = list()
			new_Z = list()
			new_T = list()
			new_delta = list()

			new_b = list()
			new_xi_i = list()

			new_nis = list()
			new_nijs = list()
			new_X = X_old[boot.ids]
			new_Z = Z_old[boot.ids]
			new_Y = Y_old[boot.ids]

			new_T = T_old[boot.ids]
			new_delta = delta_old[boot.ids]

			new_nis = nis_old[boot.ids]
			new_nijs = nij_old[boot.ids]

			X = new_X
			Z = new_Z
			Y = new_Y
			T = new_T
			delta_ij = new_delta
			n_i = new_nis
			n_ij1 = new_nijs
		
##R code to obtain estimates using this sample
			source('data_fitting_boot.R', local = FALSE)
			environment(JM_boot) <- environment()
			final_par_boot[[l]] = JM_boot()


		}
	
	
		std_samples = matrix(0, nrow = n.boot, ncol = length(initial.values))
		for(l in 1:n.boot)
			std_samples[l,] = final_par_boot[[l]][[1]]
	

		boot_se = sqrt(colVars(std_samples))

		result = list(estimates = final_parameters, LikSE = likelihood_se, BootSE = boot_se, subjectREst = final_bs, facilityREst = final_xis)
	
		result_1 = cbind(Estimates = (final_parameters), Bootstrap.std.error = (boot_se))
	
		cat(noquote("Longitudinal Outcome Results\n"))
		print(result_1[1:(p+q+1),])
		cat(noquote("Survival Outcome Results\n"))
		print(result_1[(p+q+2):(p+q+p+q+2),])
		cat(noquote("Variance of Random Effects\n"))
		print(result_1[(p+q+p+q+3):(p+q+p+q+4),])
		
	}
	
	else	
	{
		
		result = list(estimates = final_parameters, LikSE = likelihood_se, BootSE = NULL, subjectREst = final_bs, facilityREst = final_xis)

		result_1 = cbind(Estimates = (final_parameters), Likelihood.std.error = (likelihood_se))
		cat(noquote("Longitudinal Outcome Results\n"))
		print(round(result_1[1:(p+q+1),], 4))
		cat(noquote("Survival Outcome Results\n"))
		print(result_1[(p+q+2):(p+q+p+q+2),])
		cat(noquote("Variance of Random Effects\n"))
		print(result_1[(p+q+p+q+3):(p+q+p+q+4),])
		
	}
	
	###Plot the estimated baseline hazard function
	
	times = seq(0.05, 1, length = 200)
	
	thi = max(data$time_surv) + 0.05 
	tlo = 0			
	nseg = nknots
	bdeg = 2
	
	tau_est = final_parameters[(p+q+p+q+2):(length(final_parameters)-2)]
	est_baseline = exp(bbase1(times, tlo, thi, nseg, bdeg) %*% tau_est)
	quartz()
	par(cex.axis = 1.2)
	par(cex.lab = 1.2)
	par(mai = c(1.2, 1.2, 0.5, 0.5))
	plot(times, est_baseline, type="l", lty=1, ylim = c(min(est_baseline)-0.05, max(est_baseline)+0.05), ylab = expression(hat(h)[0](t)), xlab = "t", xaxs = "i", yaxs = "i")
	
	return(result)
}


###Load the following packages:

	install.packages("pracma")
	install.packages("lme4")
	install.packages("survival")
	install.packages("plyr")
	install.packages("Rfast")


	require(pracma)
	require(lme4)
	require(survival)
	require(Rfast)
	require(plyr)


###Run for the sample data:

timevar.name = "time"
xvar.names = c("X1", "X2")
zvar.names = c("Z1", "Z2")
yvar.name = "Y"

lambda = 5
nknots = 8
load('sample_data.rdata')

###Starting values:

initial.values = c(1.053, 0.49368, 0.18701, -1.53714, 2.54604, 3.13039, 0.128901, -0.158899, 0.316630, 0.917720, 0.5269, 0.9756735, 0.7672)

result_MJM = MJM_fit(data = data,  xvar.names, zvar.names, yvar.name, timevar.name, lambda, nknots, initial.values, n.boot = 50)


initials_long = glmer(Y ~ time + X1 + X2 + Z1 + Z2 +  (1 | subject_id) + (1 | facility_id), family = binomial, nAGQ = 1, data = data)

initials_surv = coxph(Surv(time_surv, survival) ~  X1 + X2 + Z1 + Z2, data = data)


getME(initials_long,"theta")^2



