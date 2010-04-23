## MANOVA functions

# 4.8.09 - also effects mu_c, mu_g and mu_cg sampled from 1 to S to cope with 2-level 3-way data

sampleMu_c = function(z,case,gender,state,mu_g,mu_s,mu_cg,mu_sc,mu_sg,mu_scg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_s)[2] # HMM-tilojen lkm.
	S = max(state)
	
	# ajasta riippumattomien parametrien vaikutus
	mu_c0 = rowSums(matrix(z[,(case==1)],nrow=K)) - sum(case==1&gender==1)*(mu_g+mu_cg)
	# Poistetaan tilariippuvien parametrien vaikutus.
	if (S>0) {
		for (s in 1:S) { # 's' goes from 2 to 'S' because time-point 1 does not have time-specific effect. - 3.4.09
			mu_c0 = mu_c0 - ( sum(state==s&case==1)*mu_s[,s] + sum(state==s&case==1)*mu_sc[,s] + sum(state==s&case==1&gender==1)*mu_sg[,s] + sum(state==s&case==1&gender==1)*mu_scg[,s] )
		}
	}
	Ns = sum(case==1)
	mu_c0 = mu_c0 / (Ns+1)
	mu_c = mvrnorm(mu=mu_c0,Sigma=diag(1/(Ns+1),nrow=K))
	
	return(mu_c)

}


sampleMu_g = function(z,case,gender,state,mu_c,mu_s,mu_cg,mu_sc,mu_sg,mu_scg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_s)[2] # HMM-tilojen lkm.
	S = max(state)

	# ajasta riippumattomien parametrien vaikutus
	mu_g0 = rowSums(matrix(z[,(gender==1)],nrow=K)) - sum(case==1&gender==1)*(mu_c+mu_cg)
	# Poistetaan tilariippuvien parametrien vaikutus.
	if (S>0) {
		for (s in 1:S) { # 's' goes from 2 to 'S' because time-point 1 does not have time-specific effect. - 3.4.09
			mu_g0 = mu_g0 - ( sum(state==s&gender==1)*mu_s[,s] + sum(state==s&case==1&gender==1)*mu_sc[,s] + sum(state==s&gender==1)*mu_sg[,s] + sum(state==s&case==1&gender==1)*mu_scg[,s] )
		}
	}
	Ns = sum(gender==1)
	mu_g0 = mu_g0 / (Ns+1)
	mu_g = mvrnorm(mu=mu_g0,Sigma=diag(1/(Ns+1),nrow=K))

	return(mu_g)

}

sampleMu_cg = function(z,case,gender,state,mu_c,mu_g,mu_s,mu_sc,mu_sg,mu_scg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_s)[2] # HMM-tilojen lkm.
	S = max(state)

	# ajasta riippumattomien parametrien vaikutus
	mu_cg0 = rowSums(matrix(z[,(case==1&gender==1)],nrow=K)) - sum(case==1&gender==1)*(mu_c+mu_g)
	# Poistetaan tilariippuvien parametrien vaikutus.
	if (S>0) {
		for (s in 1:S) { # 's' goes from 2 to 'S' because time-point 1 does not have time-specific effect. - 3.4.09
			mu_cg0 = mu_cg0 - sum(state==s&case==1&gender==1)*( mu_s[,s]+mu_sc[,s]+mu_sg[,s]+mu_scg[,s] )
		}
	}
	Ns = sum(case==1&gender==1)
	mu_cg0 = mu_cg0 / (Ns+1)
	mu_cg = mvrnorm(mu=mu_cg0,Sigma=diag(1/(Ns+1),nrow=K))

	return(mu_cg)

}

sampleMu_s = function(z,case,gender,state,mu_c,mu_g,mu_cg,mu_sc,mu_sg,mu_scg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_sc)[2] # HMM-tilojen lkm.
	S = max(state)

	if (S>0) {
		mu_s = matrix(nrow=K,ncol=S)
		mu_s[,1] = 0 # No time point 1-specific effect. -3.4.09
		for (s in 1:S) { # 's' goes from 2 to 'S' because time-point 1 does not have time-specific effect. - 3.4.09
			Ns = sum(state==s) # havaintojen lkm. HMM-tilassa 's'
			mu_s0 = ( rowSums(matrix(z[,(state==s)],nrow=K)) - ( sum(state==s&case==1)*mu_c + sum(state==s&gender==1)*mu_g + sum(state==s&case==1&gender==1)*mu_cg + sum(state==s&case==1)*mu_sc[,s] + sum(state==s&gender==1)*mu_sg[,s] + sum(state==s&case==1&gender==1)*mu_scg[,s] ) ) / (Ns+1)
			# Oletetaan kullekin 'mu_s':n klusterille sama varianssi.
			Sgm = diag(1/(Ns+1),nrow=K)
			mu_s[,s] = mvrnorm(mu=mu_s0,Sigma=Sgm)
		}
	} else {
		mu_s = matrix(0,nrow=K,ncol=1)
	}

	return(mu_s)

}

sampleMu_sc = function(z,case,gender,state,mu_c,mu_g,mu_s,mu_cg,mu_sg,mu_scg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_s)[2] # HMM-tilojen lkm.
	S = max(state)

	if (S>0) {
		mu_sc = matrix(nrow=K,ncol=S)
# 		mu_sc[,1] = 0 # No time point 1-specific effect. -3.4.09
# 		for (s in 2:S) { # No time point 1-specific effect. -3.4.09
		for (s in 1:S) { # Experiment 3.4.09: Sample mu_s for 1:S, mu_sc for 2:S, and no mu_c.
			Ns = sum(state==s&case==1) # havaintojen lkm. HMM-tilassa 's'
			mu_sc0 = ( rowSums(matrix(z[,(state==s&case==1)],nrow=K)) - (sum(state==s&case==1)*mu_c + sum(state==s&case==1&gender==1)*mu_g + sum(state==s&case==1)*mu_s[,s] + sum(state==s&case==1&gender==1)*mu_cg + sum(state==s&case==1&gender==1)*mu_sg[,s] + sum(state==s&case==1&gender==1)*mu_scg[,s] ) ) / (Ns+1)
			# Oletetaan kullekin 'mu_s':n klusterille sama varianssi.
			Sgm = diag(1/(Ns+1),nrow=K)
			mu_sc[,s] = mvrnorm(mu=mu_sc0,Sigma=Sgm)
		}
	} else {
		mu_sc = matrix(0,nrow=K,ncol=1)
	}

	return(mu_sc)

}

sampleMu_sg = function(z,case,gender,state,mu_c,mu_g,mu_s,mu_cg,mu_sc,mu_scg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_s)[2] # HMM-tilojen lkm.
	S = max(state)

	if (S>0) {
		mu_sg = matrix(nrow=K,ncol=S)
# 		mu_sg[,1] = 0 # No time point 1-specific effect. -3.4.09
# 		for (s in 2:S) { # No time point 1-specific effect. -3.4.09
		for (s in 1:S) { # Experiment 3.4.09: Sample mu_s for 1:S, mu_sc for 2:S, and no mu_c.
			Ns = sum(state==s&gender==1) # havaintojen lkm. HMM-tilassa 's'
			mu_sg0 = ( rowSums(matrix(z[,(state==s&gender==1)],nrow=K)) - ( 
							sum(state==s&case==1&gender==1)*mu_c + sum(state==s&gender==1)*mu_g + sum(state==s&gender==1)*mu_s[,s] + sum(state==s&case==1&gender==1)*mu_cg + sum(state==s&case==1&gender==1)*mu_sc[,s] + sum(state==s&case==1&gender==1)*mu_scg[,s] ) ) / (Ns+1)
			# Oletetaan kullekin 'mu_s':n klusterille sama varianssi.
			Sgm = diag(1/(Ns+1),nrow=K)
			mu_sg[,s] = mvrnorm(mu=mu_sg0,Sigma=Sgm)
		}
	} else {
		mu_sg = matrix(0,nrow=K,ncol=1)
	}

	return(mu_sg)

}

sampleMu_scg = function(z,case,gender,state,mu_c,mu_g,mu_s,mu_cg,mu_sc,mu_sg) {

	K = nrow(z) # klusterien lkm.
	#S = dim(mu_s)[2] # HMM-tilojen lkm.
	S = max(state)

	if (S>0) {
		mu_scg = matrix(nrow=K,ncol=S)
# 		mu_scg[,1] = 0 # No time point 1-specific effect. -3.4.09
# 		for (s in 2:S) { # No time point 1-specific effect. -3.4.09
		for (s in 1:S) { # Experiment 3.4.09: Sample mu_s for 1:S, mu_sc for 2:S, and no mu_c.
			Ns = sum(state==s&case==1&gender==1) # havaintojen lkm. HMM-tilassa 's'
			mu_scg0 = ( rowSums(matrix(z[,(state==s&case==1&gender==1)],nrow=K)) - sum(state==s&case==1&gender==1) * (mu_c+mu_g+mu_s[,s]+mu_cg+mu_sc[,s]+mu_sg[,s]) ) / (Ns+1)
			# Oletetaan kullekin 'mu_s':n klusterille sama varianssi.
			Sgm = diag(1/(Ns+1),nrow=K)
			mu_scg[,s] = mvrnorm(mu=mu_scg0,Sigma=Sgm)
		}
	} else {
		mu_scg = matrix(0,nrow=K,ncol=1)
	}

	return(mu_scg)

}
