# Compare two alternative ways to calculate CCA solution
# in terms of speed and accuracy

# simCCA.optimize3 is far worse than calc.pcca!!!

library(dmt)

#source("../R/pcca.R")
#source("../R/internals.R")
#source("../R/simCCA.optimize3.R")
#source("../R/optimize.parameters.R")
#source("../R/costs.R")
#source("../R/get.W.R")
#source("../R/M.set.R")
#source("../R/phi.updates.R")

####################################################

# Two data sets, different dimensions

xdim <- 8
ydim <- 5
cors1 <- c()
cors2 <- c()
for (zdim in seq(1, min(xdim, ydim), 1)) {

  print(zdim)

  toy <- generate.toydata(N = 100, zDim = zdim, xDim = xdim, yDim = ydim, marginal.covariances = "full")
  X <- toy$X
  Y <- toy$Y
  zDimension <- zdim
  
  ############################################################

  res <- calc.pcca(X, Y, zDimension)

  covX.estimated <- res$W$X%*%t(res$W$X)
  covX.true <- toy$Wx%*%t(toy$Wx)
  
  covY.estimated <- res$W$Y%*%t(res$W$Y)
  covY.true <- toy$Wy%*%t(toy$Wy)  
  
  phiX.estimated <- res$phi$X
  phiX.true <- toy$Bx%*%t(toy$Bx)

  phiY.estimated <- res$phi$Y
  phiY.true <- toy$By%*%t(toy$By)
  
  corsx <- cor(as.vector(covX.estimated), as.vector(covX.true))
  corsy <- cor(as.vector(covY.estimated), as.vector(covY.true))
  
  cormx <- cor(as.vector(phiX.estimated), as.vector(phiX.true))
  cormy <- cor(as.vector(phiY.estimated), as.vector(phiY.true))  

  cors1 <- rbind(cors1, c(corsx, corsy, cormx, cormy))

  ##############################################################
  
  res <- simCCA.optimize3(X, Y, zDimension, epsilon = 1e-3, priors = NULL, marginalCovariances = "full")

  covX.estimated <- res$W$X%*%t(res$W$X)
    covX.true <- toy$Wx%*%t(toy$Wx)
    
      covY.estimated <- res$W$Y%*%t(res$W$Y)
        covY.true <- toy$Wy%*%t(toy$Wy)
	
	  phiX.estimated <- res$phi$X
	    phiX.true <- toy$Bx%*%t(toy$Bx)
	    
	      phiY.estimated <- res$phi$Y
	        phiY.true <- toy$By%*%t(toy$By)
		
		  corsx <- cor(as.vector(covX.estimated), as.vector(covX.true))
		    corsy <- cor(as.vector(covY.estimated), as.vector(covY.true))
		    
		      cormx <- cor(as.vector(phiX.estimated), as.vector(phiX.true))
		        cormy <- cor(as.vector(phiY.estimated), as.vector(phiY.true))
			
			  cors2 <- rbind(cors2, c(corsx, corsy, cormx, cormy))
			  
			  


}

colnames(cors1) <- c("wx", "wy", "phix", "phiy")
colnames(cors2) <- c("wx", "wy", "phix", "phiy")

print(cors1)
print(cors2)


xdim <- 8
ydim <- 5
zdim <- 2 
toy <- generate.toydata(N = 100, zDim = zdim, xDim = xdim, yDim = ydim, marginal.covariances = "full")
X <- toy$X
Y <- toy$Y
t1 <- system.time(for (zDimension in 1:2) {res <- calc.pcca(X, Y, zDimension)})
t2 <- system.time(for (zDimension in 1:2) {res <- simCCA.optimize3(X, Y, zDimension, epsilon = 1e-3, priors = NULL, marginalCovariances = "full")})

print(t1)
print(t2)
