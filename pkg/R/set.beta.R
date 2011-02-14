
# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     






set.beta <-
function (M, W, phi) {
	# assuming isotropic marginal covariance
	M%*%t(W)/phi
}


set.beta.fullcov <- function (M, W, phi.inv) {
	# assuming full marginal covariance
	# as in section 4.1 EM algorithm from BachJordan probCCA paper
	#phi.inv%*%W%*%M
	M%*%t(W)%*%phi.inv
}


