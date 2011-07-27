library(dmt)

cors.list <- list()

id <- "pfa.nnw"

####################################################

# Single data set
N <- 100
xdim <- 10 # features
cors <- c()
for (zdim in seq(1, 5, 2)) {

  print(zdim)

  toy <- generate.toydata(N = N, zDim = zdim, xDim = xdim, yDim = xdim, marginal.covariances = "diagonal", priors = list(W = 1e-3))

  res <- pfa(toy$X, zDimension = zdim, priors = list(W = 1e-3))

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

pdf("pic/pfa.nnw.pdf")
plot(as.vector(covX.estimated), as.vector(covX.true)); abline(0,1)
dev.off()
