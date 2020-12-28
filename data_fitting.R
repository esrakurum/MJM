##############################
######
###### EM algorithm
######
##############################

## Assign initial values obtained from traditional models:
JM_est = function()
{
  
  beta_c = t(t_betas)
  psi_c = t(t_psi)
  zeta_c = t(t_zetas)
  eta_c = t(t_etas)
  
  alpha_c = t_alpha
  
  
  sigma_b_2_c = sigma_b_2
  sigma_xi_2_c = sigma_xi_2
  
##P-spline for baseline hazard
  thi = 1 + 0.05 
  tlo = 0			
  nseg = nknots
  bdeg = 2

 
##difference penalty
  DD <- diag(nseg+2)   
  DdD = crossprod(diff(DD, diff = 2))

  ##FOR nseg = 18
#  t_tau = c(-0.110000000, -0.098684211, -0.300000000, -0.076052632, -0.064736842, -0.053421053, -0.042105263, -0.030789474, -0.019473684, -0.008157895,  0.003157895,  0.014473684,  0.025789474,
#0.037105263,  0.048421053,  0.059736842,  0.071052632,  0.082368421,  0.093684211,  0.105000000) ###Hypothetical true value

  ##for nseg = 8
  t_tau = c(-0.41159607, -0.15659665, -0.07723866, -0.02593982,  0.01267109,  0.04391009,  0.07028668,  0.09319634,  0.11349951,  0.13176619)
  tau_c = t_tau 

##Initial values for the facility- and subject-level random effects
  
  xi_i_c = rep(0, n)  
  
  
  b_c = list()
  for(i in 1:n)
    b_c[[i]] = rep(0, n_i[i]) 
  
  pars = matrix(NA, nrow = 400, ncol = length(c(t_betas, t_psi, t_zetas,  t_etas, t_alpha, t_tau)))
  
  
##Initial likelihood value is calculated, this will be used to stop the EM
  
  
  environment(likelihood) <- environment()
  
  old_likelihoodi = rep(0, n)
  for(i in 1:n)
  { 
    b_i_c = b_c[[i]] # subject-specific random effect
    old_likelihoodi[i] = likelihood(i)
  }
  old_likelihood = sum(old_likelihoodi)	
  
  EMdiff = 11
  iterE = 0
  
  

environment(A_ij) <- environment()
environment(B_ij) <- environment()
environment(C_ij) <- environment()
environment(W_s) <- environment()
environment(D_ij) <- environment()
environment(F_ij) <- environment()
environment(logL) <- environment()
environment(logL2) <- environment()


while(EMdiff>=2 & iterE<=20) ##If the iteration number is greater than 20, stop and also if the two consecutive likelihoods are more than 1% different, carry on. 
{
  meanb.C = c()
  varb.C = c()
	varxi = rep(0, n)
	iterE = iterE+1
	varb = list()    
  corbxi =  list()
  woexp_b = c()
	woexp_xi = c()
	rij.C = c()
	####################################
	print( paste('E-step, iteration: ', iterE, sep='') ) 
	####################################
	
	
########		
######## LAPLACE
########
	Itertrack <- c()
	
	for(i in 1:n) 
	{
		#print(i)
		X_i = X[[i]] 
		Z_i = Z[[i]]
		Y_i = Y[[i]]
		T_i = T[[i]]
		b_c[[i]] = rep(0, n_i[i])
		b_i_c = b_c[[i]]
		delta_i = delta_ij[[i]]
		n_ij = n_ij1[[i]]
		iter1 = 0
		facdiff = 1	
		
		xi_i_c[i] = 0  #initial value for facility-level random effects
	
		while(facdiff >= 0.05 & iter1 <= 100)
		{	
			old_log = likelihood(i)
			iter1 = iter1+1
				
			all_As = rep(0, n_i[i])
			all_Bs = rep(0, n_i[i])
			all_Cs = rep(0, n_i[i])
			all_Ds = rep(0, n_i[i])						
			Gradient = rep(0, (n_i[i]+1))
			Hessian = matrix(0, nrow = (n_i[i]+1), ncol = (n_i[i]+1))
			for(j in 1:n_i[i])
			{
				all_As[j] = A_ij(j)
				Gradient[j] = (-all_As[j] + b_i_c[j]/sigma_b_2_c)
				all_Bs[j] =  B_ij(j)
				all_Cs[j] = C_ij(j)
				all_Ds[j] = D_ij(j)
				Hessian[j, j] =  (all_Bs[j] + 1/sigma_b_2_c)
				Hessian[j, (n_i[i]+1)] = Hessian[(n_i[i]+1), j] = all_Bs[j]
			}
		
			Gradient[(n_i[i] + 1)] =  (-sum(all_As) + xi_i_c[i]/sigma_xi_2_c)
			Hessian[(n_i[i] +1), (n_i[i] +1)] = (sum(all_Bs) + 1/sigma_xi_2_c)
			A = all_Bs - 1/sigma_b_2_c
			D =  sum(all_Bs) - 1/sigma_xi_2_c
			E = D - sum(all_Bs^2/A)
		
			FF = all_Bs/A
			
			invHess = FF %*% t(FF) / E + diag(1/A)
      invHess = cbind(invHess, -FF/E)
      invHess = (-1) * rbind(invHess, c((-FF/E),1/E))
      update =  invHess %*% Gradient

			bi_old = b_i_c
			xii_old = xi_i_c
			
			for(lp in 1:10)
			{
				s1 = 0.5^(lp-1)
				randoms_new = as.vector(c(bi_old, xii_old[i])) - s1 * as.vector(update)
				b_i_c = randoms_new[1:n_i[i]]
				xi_i_c[i] = randoms_new[(n_i[i]+1)]
				new_log = likelihood(i)	
				if(new_log<old_log) break
			}
				
			uij0 = b_i_c
			b_c[[i]] = uij0
			ui0 = xi_i_c[i]

			facdiff = max(abs(s1 * update))
 				
		}
		# Non zero entries in Aij from Supplement page 2
		C.1 = (-1) * all_Cs
		# Non zero entries in Bij from Supplement page 2
		C.2 = (-1) * all_Ds
		
		# Combine the correction term for gamma
		C.end <- c(C.1,sum(C.1))
		C.end.TTT <- invHess %*% C.end
		C.TTT <- matrix(0, n_i[i]+1,n_i[i]+1)
		for(subj in 1:n_i[i]){
		  # Construct columns of Aij
		  C.subvec <- rep(0, (n_i[i] + 1))
		  C.subvec[(n_i[i] + 1)] <- C.1[subj]
		  C.subvec[subj] <- C.1[subj]
		  # Elements in derivative of Sigma wrt c from Supplement (2)
		  C.TTT[subj,] <- invHess %*% C.subvec
		}
		C.comb <- NULL
		C.var.comb <- NULL
		C.cov.comb <- NULL
		# Correction terms for gamma 
		subjj <- n_i[i] + 1
		C.mat1 <- diag(C.TTT[,subjj]) # Derivative of Sigma wrt c from Supplement (2)
		C.mat1[(n_i[i]+1),] <- C.TTT[,subjj]
		C.mat1[,(n_i[i]+1)] <- C.TTT[,subjj]
		C.mat1[(n_i[i]+1),(n_i[i]+1)] <- C.end.TTT[subjj]
		# Correction term for gamma posterior mean, trace(V) in main paper (3)
		C.comb[subjj] <- sum(diag(invHess %*% C.mat1))
		# Correction term 1 for gamma posterior variance, trace(VV^T) in main paper (3)
		C.var.comb[subjj] <- sum(diag(invHess %*% C.mat1 %*% C.mat1 %*% invHess))
		C.matF <- C.mat1
		# Correction terms for bij
		for(subjj in 1:n_i[i]){
		  C.mat1 <- diag(C.TTT[,subjj]) # Derivative of Sigma wrt c from Supplement (2)
		  C.mat1[n_i[i]+1,] <- C.TTT[,subjj]
		  C.mat1[,n_i[i]+1] <- C.TTT[,subjj]
		  C.mat1[n_i[i]+1,n_i[i]+1] <- C.end.TTT[subjj]
		  #Correction term for bij posterior mean, trace(V) in main paper (3)
		  C.comb[subjj] <- sum(diag(invHess %*% C.mat1))
		  #Correction term 1 for bij posterior variance, trace(VV^T) in main paper (3)
		  C.var.comb[subjj] <- sum(diag(invHess %*% C.mat1 %*% C.mat1 %*% invHess))
		  #Correction term 1 for posterior covariance, trace(VV^T) in main paper (3)
		  C.cov.comb[subjj] <- sum(diag(invHess %*% C.mat1 %*% C.matF %*% invHess))
		}
		#Correction term 2 for posterior variance and covariance
		C.var.mat <- list()
		C.var.mat.2 <- list()
		C.var.TTT.2 <- matrix(0, n_i[i]+1,n_i[i]+1)
		for(subj2 in 1:n_i[i]){
		  # Second derivative of Sigma wrt ui for computing the first term in the second derivative of Sigma wrt c from Supplement (2)
		  C.var.TTT <- matrix(0, n_i[i]+1,n_i[i]+1)
		  C.var.TTT[subj2, subj2] <- C.2[subj2]
		  C.var.TTT[subj2, n_i[i] + 1] <- C.2[subj2]
		  C.var.TTT[n_i[i] + 1, subj2] <- C.2[subj2]
		  C.var.TTT[n_i[i] + 1, n_i[i] + 1] <- C.2[subj2]
		  C.var.mat[[subj2]] <- C.var.TTT
		  # First derivative of Sigma wrt ui for computing the second term in the second derivative of Sigma wrt c from Supplement (2)
		  C.var.TTT[subj2, subj2] <- C.1[subj2]
		  C.var.TTT[subj2, n_i[i] + 1] <- C.1[subj2]
		  C.var.TTT[n_i[i] + 1, subj2] <- C.1[subj2]
		  C.var.TTT[n_i[i] + 1, n_i[i] + 1] <- C.1[subj2]
		  C.var.mat.2[[subj2]] <- C.var.TTT
		  C.var.TTT.2 <- C.var.TTT.2 + C.var.TTT
		}
		C.var.mat.2[[n_i[i] + 1]] <- C.var.TTT.2
		# Compute the first term in the second derivative of Sigma wrt c from Supplement (2)
		C.var.2 <- NULL 
		C.cov.2 <- NULL
		for(subj2 in 1:(n_i[i] + 1)){
		  ttt <- (invHess[subj2,] + invHess[subj2, n_i[i]+1])^2
		  C.var.mat2 <- matrix(0, n_i[i]+1,n_i[i]+1)
		  C.cov.mat2 <- matrix(0, n_i[i]+1,n_i[i]+1)
		  for(subjj2 in 1:n_i[i]){
		    C.var.mat2 <- C.var.mat2 + C.var.mat[[subjj2]] * ttt[subjj2]
		    ttt2 <- (invHess[subj2,subjj2] + invHess[subj2, n_i[i] + 1]) * (invHess[n_i[i]+1, subjj] + invHess[n_i[i]+1,n_i[i]+1])
		    C.cov.mat2 <- C.cov.mat2 + C.var.mat[[subjj2]] * ttt2
		  }
		  # First term in correction term 2 for posterior variance 
		  C.var.2[subj2] <- sum(diag(invHess %*% C.var.mat2))
		  C.cov.2[subj2] <- sum(diag(invHess %*% C.cov.mat2))
		}
		# Compute the second term in the second derivative of Sigma wrt c from Supplement (2)
		C.var.2.2 <- NULL
		C.cov.2.2 <- NULL
		for(subj2 in 1:(n_i[i] + 1)){
		  J.mat <- -invHess %*% C.var.mat.2[[subj2]] %*% invHess %*% invHess
		  C.var.mat2.2 <- matrix(0, n_i[i]+1,n_i[i]+1)
		  C.cov.mat2.2 <- matrix(0, n_i[i]+1,n_i[i]+1)
		  for(subjj2 in 1:(n_i[i]+1)){
		    C.var.mat2.2 <- C.var.mat2.2 + C.var.mat.2[[subjj2]] * J.mat[subj2, subjj2]
		    C.cov.mat2.2 <- C.cov.mat2.2 + C.var.mat.2[[subjj2]] * J.mat[n_i[i] + 1, subjj2]
		  }
		  # Second term in correction term 2 for posterior variance 
		  C.var.2.2[subj2] <- sum(diag(invHess %*% C.var.mat2.2))
		  C.cov.2.2[subj2] <- sum(diag(invHess %*% C.cov.mat2.2))
		}
		# Combine the correction terms for posterior variance and covariance
		C.var.comb.2 <- C.var.2 + C.var.2.2
		C.cov.comb.2 <- (C.cov.2 + C.cov.2.2)[1:n_i[i]]
		# Apply correction to posterior mean as in main paper (3)
		uij0.C <- uij0 - 0.5 * C.comb[1:n_i[i]]
		ui0.C <- ui0 - 0.5 * C.comb[n_i[i] + 1]
		# Apply correction to posterior variance and covariance as in main paper (3)
		# Posterior variance of bij
		vij0.C <- diag(invHess)[1:n_i[i]] + 0.5 * C.var.comb[1:n_i[i]] - 0.5 * C.var.comb.2[1:n_i[i]]
		# Posterior variance of gamma
		vi0.C <- invHess[(n_i[i]+1), (n_i[i]+1)] + 0.5 * C.var.comb[n_i[i] + 1] - 0.5 * C.var.comb.2[n_i[i] + 1]
		# Posterior covariance
		rij0.C <- invHess[(1:n_i[i]), (n_i[i]+1)] + 0.5 * C.cov.comb - 0.5 * C.cov.comb.2
		##### End of fully exponential Laplace approximation #####
		
		xi_i_c[i] = ui0.C
		varxi[i] = vi0.C 
		b_c[[i]] = uij0.C
		varb[[i]] = vij0.C    
		corbxi[[i]] = rij0.C
		meanb.C <- c(meanb.C, uij0.C)
		varb.C <- c(varb.C, vij0.C)   
		rij.C <- c(rij.C, rij0.C)
		woexp_b = c(woexp_b, uij0.C)	
		woexp_xi = c(woexp_xi, ui0.C)
		
	}
	
##### End of E-step
	
	
##### M-step starts:
	
	print( paste('M-step, iteration: ', iterE, sep='') ) 
	
#####  Estimation of variances for the random effects    
	
   	sigma_xi_2_c = mean(xi_i_c^2 + varxi)   	
  	
    sigma_b_2_c = mean(meanb.C^2 + varb.C)

#####  NEWTON-RAPHSON Algorithm is used in the M-step to obtain the estimates of the rest of the parameters
    
	theta.start = c(beta_c, psi_c, zeta_c,  eta_c, alpha_c, tau_c)
	

	result_1 = try(optim(theta.start, fn = logL, sigma_b_2_c = sigma_b_2_c, sigma_xi_2_c = sigma_xi_2_c, hessian = FALSE, method = "BFGS", control = list(trace=6)), silent = TRUE)
	
	if (class(result_1) == "try-error")
	{
	  EMdiff = 0
	  new_likelihood = -100
	}    
	else
	{
	  
	  new_likelihood = result_1$value
	  pars_new = result_1$par	
	  
	  beta_c = t(pars_new[1:(p+1)])
	  
	  psi_c = t(pars_new[(p+2): (p+q+1)])
	  
	  zeta_c = t(pars_new[ (p+q+2): (p+q+p)])
	  
	  eta_c = t(pars_new[(p+q+p+1): (p+q+p+q)])
	  
	  alpha_c = pars_new[(p+q+p+q+1)]
	  
	  tau_c = pars_new[(p+q+p+q+2):length(pars_new)] 
	  
	  EMdiff = abs((abs(old_likelihood)-new_likelihood)/old_likelihood)*100
	  old_likelihood = new_likelihood
	}

}


print( paste('EM algorithm converged', sep='') ) 

print( paste('Likelihood-based variance calculation...', sep='') ) 

# thetas = c(pars_new, sigma_b_2_c, sigma_xi_2_c)
# 
# Hess = hessian(logL2, thetas)
# 
# ##Variances for the variance of random effects
# 
# score.sigmaxi = .5 * ((xi_i_c^2 + varxi)/sigma_xi_2_c^2 - 1/sigma_xi_2_c)
# 
# score.sigmab = rep(0, n)
# for(i in 1:n)
#   score.sigmab[i] =  sum((b_c[[i]]^2 + varb[[i]])/sigma_b_2_c^2 - 1/sigma_b_2_c)
# 
# score.sigmab = .5 * score.sigmab
# score.sig = cbind(score.sigmab,score.sigmaxi)
# inf.sig = t(score.sig) %*% score.sig
# 
# est_var_xi  = solve(inf.sig)[2,2]
# est_var_b = solve(inf.sig)[1,1]
# 
# ##Put all variances together
# 
# Variances = c(diag(solve(Hess))[1:(p+q+p+q+1)], est_var_b, est_var_xi)

Variances = 0

final_parameters = c(as.vector(pars_new), sigma_b_2_c, sigma_xi_2_c)

final_b = b_c
final_xi = xi_i_c

list(final_parameters = final_parameters, final_xi = final_xi, final_b = final_b, final_likelihood = new_likelihood, Variances = Variances)

}