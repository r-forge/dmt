# Functions of the dimension reduction part

## Function sampleFromPriorDR
# This is used in initializing the variables of the dimensionality reduction part.
# Input:
# X,Y - datasets
# xlatDim,ylatDim - number of latent variables of the datasets X and Y, respectively
# N - number of samples (i.e. nrow(X)==nrow(Y))
# Output:
# model - a list of variables 'xlat', 'ylat', 'SigmaX', 'SigmaY', 'thetaX', 'thetaY', 'Vx', 'Vy', 'varsX', 'varsY', 'muX' and 'muY'
# Tommi Suvitaival 16.4.09

sampleFromPriorDR = function(data,priors,xlatDim,ylatDim,N) {

	# Initialize residual covariance matrices of datasets 'X' and 'Y'.
	SigmaX <- sd(t(data$X))
	SigmaY <- sd(t(data$Y))
	
	# Initialize dataset-specific latent variables 'xlat' and 'ylat'.
	xlat = array(rnorm(xlatDim*ncol(data$X)),dim=c(xlatDim,ncol(data$X)))
	ylat = array(rnorm(ylatDim*ncol(data$Y)),dim=c(ylatDim,ncol(data$Y)))
	
	# Initialize dataset-specific clustering matrices 'Vx' and 'Vy'.
	thetaX <- array(1/xlatDim,dim=c(xlatDim,1)) # priors for the clusters
	thetaY <- array(1/ylatDim,dim=c(ylatDim,1))
	Vx <- array(0,dim=c(nrow(data$X),xlatDim)) # clustering matrices
	Vy <- array(0,dim=c(nrow(data$Y),ylatDim))
	init_clus <- floor(runif(length(SigmaX), min=1, max=xlatDim+0.99))
	for (i in 1:nrow(data$X)) {
		if (init_clus[i]>0) {
			Vx[i,init_clus[i]] = 1;
		}
	}
	init_clus <- floor(runif(length(SigmaY), min=1, max=ylatDim+0.99))
	for (i in 1:nrow(data$Y)) {
		if (init_clus[i]>0) {
			Vy[i,init_clus[i]] = 1;
		}
	}
	rm(init_clus)

	model = list()
	model$xlat = xlat
	model$ylat = ylat
	model$SigmaX = SigmaX
	model$SigmaY = SigmaY
	model$thetaX = thetaX
	model$thetaY = thetaY
	model$Vx = Vx
	model$Vy = Vy
	# Priors of the following dimension of data-specific variables are expected to be global.
	model$varsX = priors$varsX0
	model$varsY = priors$varsY0
	model$muX = priors$muX0
	model$muY = priors$muY0
	
	return(model)

}

## sampleMuDR
# This function was renamed from "sample_mus" on 16.4.09.
# Uusi funktio, 28.1.09.
# New prior for x_hat in form of var0. -Tommi 9.3.09

sampleMuDR <- function(x,Sigma,V,vars,xlat,mu0,var0) {

    n=ncol(x)
    m=nrow(x)
    K=nrow(xlat)
    library('mvtnorm')
    #n0=n
	 var0 = var0^2 # Take square of prior -12.3.09
    Siginv <- solve(diag(Sigma))
	 var0Inv = solve(diag(var0))
    #Cov_inv<-n*Siginv+n0*diag(1,m)
    Cov_inv<-n*Siginv+var0Inv

    Cov<-solve(Cov_inv)
    V<-matrix((vars),m,K)*V

    #x_hat <- (Cov)%*%(Siginv%*%rowMeans(x-V%*%xlat)*n+n0*diag(1,m)%*%mu0 )
    x_hat <- (Cov)%*%(Siginv%*%rowMeans(x-V%*%xlat)*n+var0Inv%*%mu0 )
	#print(dim(x_hat));print(dim(Cov))
    #Cov<-solve((V%*%t(V)+Sigma )*n+eye(m))

    mu<-mvrnorm(1,x_hat,Cov)
    mu

}

## sampleV
# 16.4.09 - This function was renamed from "sample_V".

sampleV <- function(x,z,mu,Sigma,theta,ccasta=TRUE,vars) {

	factorr = 1
	
	V <- array()
	n = ncol(x) # näytteiden lkm
	m = nrow(x) # x:n eli yhden näytteen dimensio
	K = nrow(z) # klusterien lkm
	posterior <- array(dim=c(K,1))
	#if (!ccasta) {
	#	z[K,] <- 0
	#}

	for (i in 1:m) {
			summ <- array(0,c(K,1))
		
		#Calculating the sum is a bit complicated
		for (k in 1:K) {
			for (j in 1:n) {
				summ[k] = summ[k]+max(-300,log(dnorm(x[i,j],mu[i]+vars[i]*z[k,j],sqrt(Sigma[i])) ))
			}
		}
		
		evidence <- safesum(summ+log(theta))

		# Bayes: posterior = prior * likelihood / evidence
		for (k in 1:K) {
			posterior[k] <- sum(c(summ[k],theta[k],-evidence))
		}
		posterior <- exp(posterior)
		# rmultinom(n,size,prob)
		choise = rmultinom(1, 1,t(posterior))
		if (i==1) {
			V_new <- choise
		} else {
			V_new <- cbind(V_new,choise)
		}
	}

# VERSION:
# Random sampling of one metabolite into empty cluster.
	if (arvonta) {
		for (k in 1:K) {
			if (sum(V_new[k,])==0) {
				arvonta = floor(runif(1)*m+1)
			
				V_new[,arvonta] <- 0
				V_new[k,arvonta] <- 1
			}
		}
	}
	
	V <- V_new*factorr
	V <- t(V)
	
}

## sampleSigma
# 17.4.09 - Input parameters 'case', 'gender' and 'K' were removed.
# 16.4.09 - This function was renamed from "sample_Sigmas".

sampleSigma = function(x,mu,V,z,vars) {

	K = ncol(V)
	Sigma <- array(dim=c(K,1))
	n = ncol(x)
	m = nrow(x)
	x<-x-matrix((mu),m,n)
	for (i in 1:m) {
		V[i,] <- V[i,]*vars[i]
	}
	
	temppi = (x-V%*%z)
	
	sigma_hat = NA # initialize the variable
	for (i in 1:m) {
		# Computation of variance does not work unless the vector is transposed. -Tommi 10.2.09
		sigma_hat = var(temppi[i,])
		#This is traditional chi2
		# Jäännösvarianssi, jota ei voida selittää piilomuuttujilla
		Sigma[i] = sigma_hat*n/rchisq(1, n)
		#13.5 Try Janne's suggestion
		#Sigma[i]=n/rchisq(1, n)
		#13.5 Try cutoff
		#Sigma[i]=min(1,sigma_hat*n/rchisq(1, n))
		#15.5 Use prior for variance
		#n0=2000	# data
		#v0=0.001	# priori
		
		#Sigma[i]=(sigma_hat*n+v0*n0)/rchisq(1, (n+n0))
		
		
	}
	#print("sigma_hat");print(sigma_hat)
	#print(Sigma)
	#Sigma[]<-0.01
	#There might have been a bug since this was inside the for loop
	#print(Sigma)
	sample_Sigmas <- Sigma

}

## sampleVars4
# Tämä funktio on käytössä.
# var0 on skaalaparametrin empiirinen priori.
# 17.4.09 - Function was renamed from "sample_vars4". Parameter "jama" was removed from the input. It is a global parameter now.
# Parameter 'n0' added: number of samples used in estimation of empirical prior var0. -Tommi 9.3.09
sampleVars4<-function(x,mu,V,z,Sigma,var0,n0) {

	n=ncol(x) # hiirten lkm. x aikapisteiden lkm.
	m=nrow(x)
	K=nrow(z)
	vars<-Sigma
	vars[] <- 0
	var0 = var0^2 # Take square of prior -12.3.09
	meaan<-Sigma
	meaan[]<-0
	meaan3<-meaan
	for (j in 1:m) {
		klusteri<-which(V[j,]>0)
		if (klusteri==K & jama==1) {
			Vb=0
		} else {
			Vb<-1/sum(z[klusteri,]*z[klusteri,])
		}
		#8.10.2008 This is a strange attempt trying to make var(z)=1
		#Vb=1
		#17.9.2008 truncate the scale parameter such that it can't be negative since that would ruin the model. This might be a bit questionable
		meaan3[j]<-max(Vb*sum(z[klusteri,]*(x[j,]-mu[j])),0)
		vars[j]=(meaan3[j]*n+var0[j]*n0)/rchisq(1,(n+n0))
		#vars[j]=(meaan3[j]*(n)+var0[j]*(20))/rchisq(1, (n+(20)))
		#vars[j]=(meaan3[j]*(n)+var0[j]*(0))/rchisq(1, (n+0))
	}
	vars<-sqrt(vars)
	vars

}

## sampleXlatsWithPriorDR
# This is the sampling function that is to be used in the DRCCA. Effects need to be added into the sampling.
# 17.4.09 - Renamed from "sample_xlats_with_prior".

sampleXlatsWithPriorDR <- function(x,mu,Sigma,V,z,prior_mu,prior_cov,alku,vars) {

	n=ncol(x)
	m=nrow(x)
	K=nrow(z)
	xlatnew<-array(dim=dim(z))
	
	varr2<-solve(t(V)%*%solve(diag(Sigma))%*%V+solve(prior_cov))
	meaan2<-varr2%*%(t(V)%*%solve(diag(Sigma))%*%x+solve(prior_cov)%*%prior_mu)
	xlatnew<-meaan2+t(mvrnorm(n = ncol(meaan2),t(t(rep(0,K))),varr2, tol = 1e-6, empirical = FALSE))

	if (jama) { # Set latent variables of "jämä" cluster to zero.
		xlatnew[K,] = 0
	}
	
	like = 0
	for (n in 1:N) { # Compute likelihood of the new variable values -24.3.10
		like = like + dmvnorm(x=xlatnew[,n], mean=z[,n], sigma=prior_cov, log=T)
	}
	
	return(xLat=xlatnew, like=like)

}

## sampleXlats2
# 30.4.09 - Input parameter "z" was removed as it is not used.
# 17.4.09 - Handling of possible "jämä" cluster was moved inside this function.
# 16.4.09 - This function was renamed from "sample_xlats2".

sampleXlats2 = function(x,case,gender,state,mu,mu_c,mu_g,mu_s,mu_cg,mu_sc,mu_sg,mu_scg,Sigma,V,vars) {

	N = ncol(x) # näytteiden lkm.
	M = nrow(x) # metaboliittien lkm.
	K = ncol(V) # latent variable dimension (number of clusters)
	x<-x-matrix((mu),M,N)
	if (!is.null(vars)) { # Scale the projection matrix 'V' row-wise by vector 'vars'.
		for (m in 1:M) {
			V[m,] <- V[m,]*vars[m]
		}
	}
	xlatnew <- array(dim=c(K,N))
	if (is.matrix(Sigma)) { # 'Sigma' is already a matrix.
		invSigma = solve(Sigma)
	} else { # 'Sigma' is a vector which then is expanded to a diagonal matrix.
		invSigma = diag(1/Sigma)
	}
	covvi = diag(1,K)
	varr = solve(t(V)%*%invSigma%*%V+covvi)

	# effects that are independent of covariate 'state'
	# The products were changed to use 'outer()'. -Tommi 2.2.09
	effects = outer(mu_c,case) + outer(mu_g,gender) + outer(mu_cg,(case&gender)*1)
# 	mean_n = t(V)%*%invSigma%*%x + outer(mu_c,case) + outer(mu_g,gender) + outer(mu_cg,(case&gender)*1)
	# effects that are dependent of the covariate 'state'
	if (sum(state==0)==0) {
		# Create a vector of ones used in extending other vectors to matrices. 'outer(k_ones,x)' copies vector 'x' as 'k_ones' number of rows in a matrix. 'outer()' is more efficient than operator '%*%' - Tommi 30.1.09
		k_ones = rep(1,K)
# 		mean_n = mean_n + mu_s[,state] + mu_sc[,state]*(k_ones%o%case) + mu_sg[,state]*(k_ones%o%gender) + mu_scg[,state]*(k_ones%o%(case&gender))
		effects = effects + mu_s[,state] + mu_sc[,state]*(k_ones%o%case) + mu_sg[,state]*(k_ones%o%gender) + mu_scg[,state]*(k_ones%o%(case&gender))
	}
# 	mean_n = t(V)%*%invSigma%*%x

# 	mean_n =  varr %*% mean_n
	mean_n =  varr %*% (crossprod(V,invSigma)%*%x+effects)
	# Arvotaan uudet 'x_lat'-arvot normaalijakaumasta kullekin näytteelle
	# yllä lasketuilla odotusarvoilla.
	for (n in 1:N) {
		xlatnew[,n] = rmvnorm(n=1,mean_n[,n],varr)
	}

	if (jama) { # Set latent variables of "jämä" cluster to zero.
		xlatnew[K,] = 0
	}
	
	like = 0
	for (k in 1:K) { # Compute likelihood of the new variable values -24.3.10
		like = like + sum(dnorm(xlatnew[k,], mean=effects[k,]), log=T)
	}
	
	
	return(list(z=xlatnew, like=like))

}

## safesum

safesum <- function(x) { 

	xmax <- max(x);
	safesum<-xmax+log(sum(exp(x-xmax)))

}