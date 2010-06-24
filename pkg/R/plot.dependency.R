plot.dependency <- function (X, Y, type) {
  
  if (type == "sw") {

    models <- dependency.screen(X, Y, verbose = FALSE)

    plot(locs, getScore(models), type = 'l', main = "Shared signal", ylab = "Score", cex.main = 2, cex.lab = 1.5, cex.axis =1.8 , xlab = "Chromosomal location (Mb)")
  } else {
    message("Under construction.") # FIXME
  
  }

}