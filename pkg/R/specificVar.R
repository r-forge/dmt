
# this program is like sharedVar.R but it will calculate the variances
# with in each data sets. This finds the within data variation captured by
# cca and pca projections.

#input 1. datasets : list of data sets
#input 2. regcca : output of regCCA
#input 3. dim : number of dimensions to be used for back projection


#input 4. pca = FALSE : default parameter

# output : a list of 2 vectors containing within data variation
#          for cca projected data and for pca projected data

"specificVar" <- 
function(datasets, regcca,dim, pca = FALSE)
{


#########################################################
#####subroutines start heremi/doc/highlights/drafts2006/


      reverse <- function(proj,dir,dim)
      {

        proj <- proj #projected data
        dir <- dir # direction matrix

        dim <- dim # number of dimension to use

              if(dim > ncol(proj)) stop("Number of input dimension is greater than the number of dimensions in the projected data")



        p <- proj[,1:dim]

        d <- dir[,1:dim]

        inv_d <- solve(t(d)%*% d)%*% t(d)
             
        x <- p  %*% inv_d # same as dimension of proj
      
       return(x)
       }

#subroutine for pca

    withinpca <- function(mat,dim)
    {

     mat <- mat # list of data matrices

     dim <- dim # number of dimensions of projected data to be used to
              # regenerate original data sets

     m <- length(mat)


     fea <- c(rep(0,m)) # numebr of dimensions in each data matrix
     for(i in 1:m)
         {

           fea[i] <- ncol(mat[[i]])
         }

        #concatenate all matrices

          x <- concatenate(mat)


      #performing pca#

         
          pcax <- prcomp(x) #pca of x


          z <- pcax$x  # z is pca rotated data
          eig_z <- pcax$rotation # pca directions
 
      # reverse the pca projection, z[,1:dim]

          x_z <- reverse(z,eig_z,dim)

      # separating all matrices


             arr_xz <- array(list(),m) #storing matrices form x_z

             c <- 0
             for(i in 1:m)
             {

               t <- x_z[,(c+1):(c+fea[i])] #ith data set

               #storing matrices in a array
                arr_xz[[i]] <- t

                c <- c + fea[i]
              }
 
         
       # calculating within data  variation

            pc_var <- c(rep(0,m))

           for(i in 1:m)
           {
       
            var <- cov(arr_xz[[i]])

            dummy <- svd(var)

            pc_var[[i]] <- sum(dummy$d)

          }

           return(pc_var)

        }

#########################################################
# MAIN PROGRAM STARTS HERE



     mat <- datasets  # list of matrices

     cc <- regcca # output of regCCA
     eigV <- cc$eigvecs #, whitened projections directions
     proj <- cc$proj #, projected data from regCCA
     white <- cc$white #, whitening matrices from regCCA
     
    dim <- dim # number of dimensions to be used for back projection


     m <- length(mat)


           covm <- array(list(),m) #covariance matrices
           fea <- c(rep(0,m)) # number of dimensions in each data matrix

           for(i in 1:m)
           {
             
            mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})
           
            covm[[i]]<- cov(mat[[i]])  #covariance matrix

            fea[i] <- ncol(mat[[i]])
           }


        ##inverse whitening for projection directions
       
          for(i in 1:m)
          {
            eigV[[i]] <- solve(white[[i]]) %*% eigV[[i]]
          }

       ## and concatenating all eigV's rowwise

           eigfull <- eigV[[1]]
           for(i in 2:m)
           {
             eigfull <- rbind(eigfull,eigV[[i]])
           }

      # reverse process for regcca proj data for a given 'dim'

         fullrev_cca <- reverse(proj,eigfull,dim)

      # separate all matrices

          rev_cca <- array(list(),m) # list of regenerated data sets from
                                     # proj for given 'dim
          
             c <- 0
             for(i in 1:m)
             {

               t <- fullrev_cca[,(c+1):(c+fea[i])] #ith data set

               t <- t %*% solve(white[[i]])

               #storing matrices in a array
                rev_cca[[i]] <- t

                c <- c + fea[i]
              }

         # within data variation for regCCA projected data

           cc_var <- c(rep(0,m)) #variance for each data set

                 for(i in 1:m)
                 {

                   var <- cov(rev_cca[[i]])

                   dummy <- svd(var)

                   cc_var[[i]] <- sum(dummy$d)
                 }

         # within data variance for original data sets

           orig <- c(rep(0,m)) #for each data set

                 for(i in 1:m)
                 {
                    y <- svd(covm[[i]])

                    orig[i] <- sum(y$d)
                 }

         # mean of within data variance for  drCCA

           mcc <- mean(cc_var/orig)


          if(pca)
           {
      pc_var <- withinpca(mat,dim)

      mpc <- mean(pc_var/orig)

     result <- list(cc= sum(cc_var)/sum(orig), pc = sum(pc_var)/sum(orig),mcca = mcc, mpca=mpc)  

     return(result)
   }else{ return(cc =sum(cc_var)/sum(orig), mcca = mcc)}


}

