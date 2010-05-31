
# this program will run sharedVar.R and specificVar.R for a given set of
# dimensions which will be used to project the data back for regCCA and PCA.

# it will also have an option to plot the variations for two programs

#input 1. datasets <- a list of data matrices
#input 2. regcca <- output of regCCA
#input 2. dimVector <- a vectors of integers, dimensions which will be used for
                       # for back projection
#input 3. plot = FALSE , if it is TRUE, a plot will be generated

"plotVar" <- 
function(datasets,regcca, dimVector, plot = FALSE)
{


        mat <- datasets        

        regcca <- regcca

        dims <- dimVector # list of values for dimensions

        m <- length(dims)

   ### pwiseVar.R

         pwcca <- c(rep(0,m)) # sum of pwise variation for each dims[i]

         pwpca <- c(rep(0,m)) # sum of pwise variation for each dims[i]

         for(i in 1:m)
         {

           dummy <- sharedVar(mat,regcca,dims[i],pca = TRUE)

           pwcca[i] <- sum(dummy$cc)/sum(dummy$oo) 

           pwpca[i] <- sum(dummy$pc)/sum(dummy$oo)
         }



   ### specificVar.R 

         withcca <- c(rep(0,m)) # sum of within data variation for each dims[i]

         withpca <- c(rep(0,m)) # sum of within data variation for each dims[i]

        for(i in 1:m)
         {

           dummy <- specificVar(mat,regcca,dims[i], pca = TRUE)

            withcca[i] <- sum(dummy$cc)


            withpca[i] <- sum(dummy$pc)
         }





      


 
   ###plotting the graphs 

  if(plot)
  {
   pw <- rbind(pwcca,pwpca)
   wn <- rbind(withcca,withpca)

   split.screen(c(2,1))
   plot(pw[1,],t='l',xlab='Dimensionality of the projection', ylab='Shared variation', main='Shared variation, drCCA vs PCA(red)', ylim=c(min(pw),max(pw)))
   lines(pw[2,], col ='red')

     
   screen(2)
   plot(wn[1,],t='l',xlab='Dimensionality of the projection', ylab='Data-specific variation', main='Data-specific variation, drCCA vs PCA(red)', ylim=c(min(wn),max(wn)))
   lines(wn[2,], col='red')

   close.screen(all = TRUE)
  }

 return(list(pw_cca = pwcca, pw_pca = pwpca, within_cca = withcca, within_pca = withpca))

   

}