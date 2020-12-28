
tpower <- function(x, t, p)
  # Function for truncated p-th power function
  (x - t) ^ p * (x > t)


 # Function for B-spline basis

bbase1 <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B 
	}


invLogit = function(x)
  return(exp(x) / (1 + exp(x)))

wk <- c(0.0229353220105292, 0.0630920926299785, 0.10479001032225, 
        0.140653259715526, 0.169004726639268, 0.190350578064785, 
        0.204432940075299, 0.209482141084728, 0.204432940075299, 
        0.190350578064785, 0.169004726639268, 0.140653259715526, 
        0.10479001032225, 0.0630920926299785, 0.0229353220105292)
sk <- c(-0.991455371120813, -0.949107912342758, -0.864864423359769, 
        -0.741531185599394, -0.586087235467691, -0.405845151377397, 
        -0.207784955007899, 0, 0.207784955007899, 0.405845151377397, 
        0.586087235467691, 0.741531185599394, 0.864864423359769, 
        0.949107912342758, 0.991455371120813)


likelihood = function(i)
{
	old_logl = 0
	X_i = X[[i]] 
	Z_i = Z[[i]]
	Y_i = Y[[i]]
	T_i = T[[i]]
	delta_i = delta_ij[[i]]
	n_ij = n_ij1[[i]]
	penalty = lambda * (t(tau_c) %*% DdD %*% tau_c)[1, 1]
	
	for(j in 1:n_i[i])
	{
		
		gk = 0.5 * (T_i[j] * sk + T_i[j])
		
		new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
		mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
		Mm = invLogit(mu_ij)
				
				
		first_part = delta_i[j] * (bbase1(T_i[j], tlo, thi, nseg, bdeg) %*% tau_c  + X_i[(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c) + alpha_c * Mm)
		
		out = 0
				
		for(tty in 1:15)
		{
				
			mu_ij = beta_c[1] + beta_c[2]*gk[tty] + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
			Mm = invLogit(mu_ij)
					
			out = wk[tty] * exp(bbase1(gk[tty], tlo, thi, nseg, bdeg) %*% tau_c + alpha_c * Mm) + out
			
		}
		second_part = exp(X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c)) * out * (T_i[j]/2)
		if(n_ij[j]==1)
		mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
		else
		mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
		Mm = invLogit(mu_ij)
		third_part = sum(mu_ij * Y_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j])] + log(1 - Mm)) - 0.5 * (b_i_c[j]/sigma_b_2_c) - 0.5*log(2 * pi * sigma_b_2_c)
		old_logl = old_logl + first_part - second_part + third_part 
	}
	old_logl = old_logl - (0.5) * penalty - 0.5 * (xi_i_c[i]/sigma_xi_2_c) - 0.5 * log(2 * pi * sigma_xi_2_c)
	 
current_l = (-1)*old_logl
return(current_l)
}

##For sim:
logL = function(theta, sigma_b_2_c, sigma_xi_2_c)
{


	beta_c = t(c(theta[[1]], theta[[2]], theta[[3]], theta[[4]]))
	psi_c = t(c(theta[[5]], theta[[6]]))

 	zeta_c = t(c(theta[[7]], theta[[8]]))
	eta_c = t(c(theta[[9]], theta[[10]]))
	
	alpha_c = theta[[11]]
	tau_c = theta[12:length(theta)]
	old_logl = 0
	penalty = lambda * (t(tau_c) %*% DdD %*% tau_c)[1, 1]


	for(i in 1:n) 
	{#print(i)
		X_i = X[[i]] 
		Z_i = Z[[i]]
		Y_i = Y[[i]]
		T_i = T[[i]]

		b_i_c = b_c[[i]]
		

		delta_i = delta_ij[[i]]
		n_ij = n_ij1[[i]]

		for(j in 1:n_i[i])
		{
			gk = 0.5 * (T_i[j] * sk + T_i[j])

			#print(j)
			key = (varxi[i]  +  varb[[i]][j] + 2 * corbxi[[i]][j])/2
			new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
			mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
			Mm = invLogit(mu_ij)
				
				
			first_part = delta_i[j] * (bbase1(T_i[j], tlo, thi, nseg, bdeg) %*% tau_c  + X_i[(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c) + alpha_c * Mm)
				
			B_first = delta_i[j] * alpha_c * Mm * (1-Mm) * (1 - 2 * Mm)
			
			out = 0
			
			for(tty in 1:15)
			{
				
				mu_ij = beta_c[1] + beta_c[2]*gk[tty] + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
				Mm = invLogit(mu_ij)
					
				out = wk[tty] * exp(bbase1(gk[tty], tlo, thi, nseg, bdeg) %*% tau_c + alpha_c * Mm) + out
			
			}
			second_part = exp(X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c)) * out * (T_i[j]/2)
			if(n_ij[j]==1)
				mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
			else
				mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
			Mm = invLogit(mu_ij)
			B_third = sum(Mm * (1 - Mm))
			third_part = sum(mu_ij * Y_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j])] + log(1 - Mm)) - 0.5 * ((b_i_c[j]^2 +varb[[i]][j])/sigma_b_2_c)  - 0.5*log(2 * pi * sigma_b_2_c) 
			key_part = key * 0.5 * (B_first - second_part - B_third)
			old_logl = old_logl + first_part - second_part + third_part 
				}
	
			old_logl = old_logl - (0.5) * penalty - 0.5 * ((xi_i_c[i]^2 + varxi[i])/sigma_xi_2_c) - 0.5 * log(2 * pi * sigma_xi_2_c)	 
	 }
	 
current_l = (-1)*old_logl
#print(paste("current_l", current_l))
return(current_l)
}


logL2 = function(thetas)
{

	theta = thetas[1:(length(thetas)-2)]
	sigma_b_2_c = thetas[(length(thetas)-1)]
	sigma_xi_2_c = thetas[length(thetas)]
		
	beta_c = t(c(theta[1], theta[2], theta[3], theta[4]))
	psi_c = t(c(theta[5], theta[6]))

 	zeta_c = t(c(theta[7], theta[8]))
	eta_c = t(c(theta[9], theta[10]))
	
	alpha_c = theta[11]
	tau_c = theta[12:length(theta)]
	old_logl = 0
	penalty = lambda * (t(tau_c) %*% DdD %*% tau_c)[1, 1]


	for(i in 1:n) 
	{#print(i)
		X_i = X[[i]] 
		Z_i = Z[[i]]
		Y_i = Y[[i]]
		T_i = T[[i]]

		b_i_c = b_c[[i]]
		delta_i = delta_ij[[i]]
		n_ij = n_ij1[[i]]

		for(j in 1:n_i[i])
		{
			#print(j)
			key = (varxi[i]  +  varb[[i]][j] + 2 * corbxi[[i]][j])/2
			new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
			mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
			Mm = invLogit(mu_ij)
				
				
			first_part = delta_i[j] * (tau_c %*% t(bbase1(T_i[j], tlo, thi, nseg, bdeg))  + X_i[(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c) + alpha_c * Mm)
				
			B_first = 	delta_i[j] * alpha_c * Mm * (1-Mm) * (1 - 2 * Mm)
			int_part1 = function(s)
			{
				mu_ij = beta_c[1] + beta_c[2]*s + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
				Mm = invLogit(mu_ij)
					
				out = exp((tau_c %*% t(bbase1(s, tlo, thi, nseg, bdeg)))[1, 1] * X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c) + alpha_c * Mm) 
				return(out)
			}

		if(T_i[j]>=0.0000001)
		second_part = suppressWarnings(gauss_kronrod(int_part1, 0.0000001, T_i[j])$value)
		else 
		second_part = suppressWarnings(gauss_kronrod(int_part1, 0, T_i[j])$value)
						#print(second_part)
			if(n_ij[j]==1)
				mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
			else
				mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
			Mm = invLogit(mu_ij)
			B_third = sum(Mm * (1 - Mm))
			third_part = sum(mu_ij * Y_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j])] + log(1 - Mm)) - 0.5 * ((b_i_c[j]^2 +varb[[i]][j])/sigma_b_2_c)  - 0.5*log(2 * pi * sigma_b_2_c) 
			key_part = key * 0.5 * (B_first - second_part - B_third)
			old_logl = old_logl + first_part - second_part + third_part 
				}
	
			old_logl = old_logl - (0.5) * penalty - 0.5 * ((xi_i_c[i]^2 + varxi[i])/sigma_xi_2_c) - 0.5 * log(2 * pi * sigma_xi_2_c)	 
	 }
	 
current_l = (-1)*old_logl
#print(paste("current_l", current_l))
return(current_l)
}




A_ij = function(j) 
{
	new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
	mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
	Mm = invLogit(mu_ij)
	first_part = delta_i[j] * alpha_c * Mm * (1 - Mm)
		
	gk = 0.5 * (T_i[j] * sk + T_i[j])
	
		
	out = 0
	for(tty in 1:15)
	{
				
		mu_ij = beta_c[1] + beta_c[2]*gk[tty] + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
		Mm = invLogit(mu_ij)
					
		out = wk[tty] * exp(bbase1(gk[tty], tlo, thi, nseg, bdeg) %*% tau_c) * alpha_c * Mm * (1 - Mm) * exp(alpha_c * Mm) + out 

	}
	
	second_part = exp(X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c)) * out * (T_i[j]/2)

	if(n_ij[j]==1)
	mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	else
	mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	Mm = invLogit(mu_ij)
	third_part = sum(Y_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j])] - Mm)
	result = first_part - second_part + third_part
	return(result)
}
		
B_ij = function(j)
{ 
	gk = 0.5 * (T_i[j] * sk + T_i[j])

	new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
	mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
	Mm = invLogit(mu_ij)

	first_part = delta_i[j] * alpha_c * Mm * (1-Mm) * (1 - 2 * Mm)
			
	out = 0
	for(tty in 1:15)
	{
				
		mu_ij = beta_c[1] + beta_c[2]*gk[tty] + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
		Mm = invLogit(mu_ij)
		out = wk[tty] * exp(bbase1(gk[tty], tlo, thi, nseg, bdeg) %*% tau_c) * alpha_c * Mm * (1 - Mm) * exp(alpha_c * Mm)  * (alpha_c * Mm * (1 - Mm) + (1 - 2* Mm))  + out 

	}
	second_part = exp(X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c)) * out * (T_i[j]/2)
	
	if(n_ij[j]==1)
		mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	else
		mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
		Mm = invLogit(mu_ij)
		third_part = sum(Mm * (1 - Mm))
		result = first_part - second_part - third_part
		return(result)
	}



C_ij = function(j)
{ 
	gk = 0.5 * (T_i[j] * sk + T_i[j])

	new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
	mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
	Mm = invLogit(mu_ij)

	first_part = delta_i[j] * alpha_c * Mm * (1-Mm) *  ((1 - 2 * Mm)^2 - 2* Mm * (1 - Mm))
			
	out = 0
	for(tty in 1:15)
	{	
		mu_ij = beta_c[1] + beta_c[2]*gk[tty] + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
		Mm = invLogit(mu_ij)
		out = wk[tty] * exp(bbase1(gk[tty], tlo, thi, nseg, bdeg) %*% tau_c) * alpha_c * Mm * (1 - Mm) * exp(alpha_c * Mm)  * ((alpha_c * Mm * (1 - Mm) + (1 - 2*Mm))^2 + Mm * (1 - Mm) * (alpha_c * (1 - Mm) - alpha_c * Mm -2)) + out

	}
	second_part = exp(X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c)) * out * (T_i[j]/2)
		
	if(n_ij[j]==1)
		mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	else
		mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	Mm = invLogit(mu_ij)
	third_part = sum(Mm * (1 - Mm) * (1 - 2 * Mm))
	result = first_part - second_part - third_part
	return(result)
}


W_s = function(s)
{
		mu_ij = beta_c[1] + beta_c[2]*s + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
		Mm = invLogit(mu_ij)
		P_ij = alpha_c * Mm * (1 - Mm)		
		out =  (P_ij + (1 - 2*Mm))^3 + Mm * (1 - Mm) * (alpha_c * (1 - Mm) - alpha_c * Mm - 2) * (3 * P_ij + 4 * (1 - 2*Mm)) - 2 * alpha_c * Mm^2 * (1 - Mm)^2 
		return(out)
}

D_ij = function(j)
{ 
	new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
	mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
	Mm = invLogit(mu_ij)

	first_part = delta_i[j] * alpha_c * Mm * (1-Mm) *  ((1 - Mm)^3 - 11 * Mm * (1-Mm)^2 + 11 * Mm^2 * (1 - Mm) - Mm^3)
	
	gk = 0.5 * (T_i[j] * sk + T_i[j])
		
	out = 0
	for(tty in 1:15)
	{	
		mu_ij = beta_c[1] + beta_c[2]*gk[tty] + beta_c[3:(p+1)] %*% X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
		out = wk[tty] * exp(bbase1(gk[tty], tlo, thi, nseg, bdeg) %*% tau_c)* alpha_c * Mm * (1 - Mm) * exp(alpha_c * Mm) * W_s(gk[tty]) + out

	}
	second_part = exp(X[[i]][(sum(n_ij[0:(j-1)])+1), 2:p] %*% t(zeta_c) + Z_i[j, ] %*% t(eta_c)) * out * (T_i[j]/2)
	
	
	if(n_ij[j]==1)
		mu_ij = c(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	else
		mu_ij = cbind(1, X_i[(sum(n_ij[0:(j-1)])+1):sum(n_ij[0:j]), ]) %*% t(beta_c) + rep(Z_i[j, ] %*% t(psi_c),  n_ij[j]) + rep(b_i_c[j], n_ij[j]) + rep(xi_i_c[i], n_ij[j])		
	Mm = invLogit(mu_ij)
	third_part = sum(Mm * (1 - Mm) * ((1 - 2 * Mm)^2 - 2 * Mm * (1 - Mm)))
	result = first_part - second_part - third_part
	return(result)
}
		
F_ij = function(j)
{
	new_X = c(1, T_i[j], X_i[(sum(n_ij[0:(j-1)])+1), 2:p])
	mu_ij = new_X %*% t(beta_c) + Z_i[j, ] %*% t(psi_c) + b_i_c[j] + xi_i_c[i]
	Mm = invLogit(mu_ij)
	
	P_ij = alpha_c * Mm * (1 - Mm)

	result = delta_i[j] * P_ij * (( (1 - 2 * Mm)^2 - 2* Mm * (1-Mm)) * (key) + 1)
	return(result)
}
