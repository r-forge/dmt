library(dmt)

cors.list <- list()

####################################################

# Two data sets, nonmatched dimensions case

N <- 100
xdim <- 8
ydim <- 5
cors <- c()
for (zdim in seq(1, min(xdim, ydim), 1)) {

  print(zdim)

  toy <- generate.toydata(N = N, zDim = zdim, xDim = xdim, yDim = ydim, marginal.covariances = "isotropic")

  res <- ppca(toy$X, toy$Y, zDimension = zdim)

  covX.estimated <- res@W$X%*%t(res@W$X)
  covX.true <- toy$Wx%*%t(toy$Wx)
  
  covY.estimated <- res@W$Y%*%t(res@W$Y)
  covY.true <- toy$Wy%*%t(toy$Wy)  
  
  phiX.estimated <- res@phi$X
  phiX.true <- toy$Bx%*%t(toy$Bx)

  phiY.estimated <- res@phi$Y
  phiY.true <- toy$By%*%t(toy$By)

  corsx <- cor(as.vector(covX.estimated), as.vector(covX.true))
  corsy <- cor(as.vector(covY.estimated), as.vector(covY.true))
  
  cormx <- cor(as.vector(phiX.estimated), as.vector(phiX.true))
  cormy <- cor(as.vector(phiY.estimated), as.vector(phiY.true))  

  cors <- rbind(cors, c(corsx, corsy, cormx, cormy))

}

colnames(cors) <- c("wx", "wy", "phix", "phiy")

print(cors)

cors.list[["two data sets"]] <- cors

#######################################################

# One data set

xdim <- 10
cors <- c()
for (zdim in seq(1, xdim, 2)) {

  print(zdim)

  toy <- generate.toydata(N = N, zDim = zdim, xDim = xdim, yDim = xdim, marginal.covariances = "isotropic")

  res <- ppca(toy$X, zDimension = zdim)

  covX.estimated <- res@W$total%*%t(res@W$total)
  covX.true <- toy$Wx%*%t(toy$Wx)
  
  phiX.estimated <- res@phi$total
  phiX.true <- toy$Bx%*%t(toy$Bx)

  corsx <- cor(as.vector(covX.estimated), as.vector(covX.true))  
  cormx <- cor(as.vector(phiX.estimated), as.vector(phiX.true))

  cors <- rbind(cors, c(corsx, cormx))

}

colnames(cors) <- c("wx", "phix")
print(cors)

cors.list[["one data set"]] <- cors

#################################################################

# TODO:
# regularization missing -> give option to set priors also for ppca function
# tests with higher noise levels (now quite moderate)

# NOTE: in ppca and pfa (maybe also pcca?) marginal matrices correlate
# well but are typically estimated in about (not exactly) double the scale of the truth?

# Latent variables
# Z <- z.expectation(res, X, Y)

pdf("pic/ppca.pdf")
plot(as.vector(covX.estimated), as.vector(covX.true)); abline(0,1)
dev.off()
