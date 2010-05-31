   
#It will be used to calculate cca and pca and then to reverse the 
# the projected data back to original one using first few projected
# components; then calculating pairwise mutual information and comparing
# the result for drcca and PCA

#input 1. file: a file containing nanmes of the matrices
#inout 2. regcca_output : output of regCCA
#input 3. dimension : number of dimensions in projected data to be used to
#                     regenarate original data sets
#input 4. pca : if pca = TRUE , only then the pca variances will be calculated
#output : returns a list of three matrices, one for original data sets, one for cca and other for pca
#         entries are pairwise mutual information



"sharedVar" <- 
function(datasets,regcca,dimension,pca=FALSE)
{

        mat <- datasets # list of data matrices

        cc <- regcca  # projection direction from regCCA

        #print(length(cc))
        #print(dim(cc[[3]]))
        eigV <- cc$eigvecs #, whitened projections directions
        proj <- cc$proj #, projected data from regCCA
        white <- cc$white #, whitening matrices from regCCA
        xmean <- cc$meanvec #, array of columnwise mean for each data matrix


        dim <- dimension #dimensions to be used for back projection


        m <- length(mat)

#############################################################
#####subroutines start here


#subroutine 1
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

### subroutine 2  ###

 xvar <- function(index,array)
 {

       ind <- index

       arr <- array

       len <- length(arr)

        var <- c(rep(0,(len))) # to store pairwise mi with each index

        for(i in (ind+1):len)
        {

           dummy <- t(arr[[ind]]) %*% arr[[i]]

           scov <- svd(dummy)

           var[i] <- sum(scov$d)

           #den1 <- svd(cov(arr[[ind]]))
           #d1 <- sum(den1$d)
           #den2 <- svd(cov(arr[[i]]))
           #d2 <- sum(den2$d)
           #var[i] <- sum(scov$d)/sqrt(d1*d2)
           
        }

       return(var)

 }

##pca subroutine 3

pwisepca <- function(mat,dim)
{

    mat <- mat # array of matrices
 
    dim <- dim 

    m <- length(mat)

    fea <- c(rep(0,m)) # numebr of dimensions in each data matrix
    for(i in 1:m)
         {

            mat[[i]] <- mat[[i]]
           #mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})

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
 
         
       # calculating pairwise variation

        cross_pca <- matrix(0,(m-1),m)

       for(i in 1:(m-1))
       {
       
       cross_pca[i,] <- xvar(i,arr_xz)

       #print(offdiag_pca[[i]])
       }

       return(cross_pca)

}


################################################################




          
           fea <- c(rep(0,m)) # number of dimensions in each data matrix
           for(i in 1:m)
           {
             
            mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})
           

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


       # reverse process for each data set for a given 'dim'

         fullrev_cca <- reverse(proj,eigfull,dim)

       # separate all matrices


          rev_cca <- array(list(),m) # list of regenerated adat sets from
                                     # y for given 'dim
          
             c <- 0
             for(i in 1:m)
             {

               t <- fullrev_cca[,(c+1):(c+fea[i])] #ith data set

               t <- t %*% solve(white[[i]])

               #storing matrices in a array
                rev_cca[[i]] <- t

                #print(rev_cca[[i]][1,1:4])

                c <- c + fea[i]
              }

      
    ########## CALCULATING PAIRWISE CROSS VARIATION ##########

    # pairwise cross variation for drcca, rev_cca

      cross_cca <- matrix(0,(m-1),m) 

      for(i in 1:(m-1))
      {

       cross_cca[i,] <- xvar(i,rev_cca) # var of x1%*%t(xj)

       #print(offdiag_cca[[i]])
      }

     #pairwise cross variation for original data sets

      orig <- matrix(0,(m-1),m)
      for(i in 1:(m-1))
      {

       orig[i,] <- xvar(i,mat)

      }

     # calculating the mean of pairwise variation

      mean_cc <- mean(cross_cca[cross_cca !=0]) # for drCCA
        


      mean_oo <- mean(orig[orig !=0]) # for original

           mcc <- mean_cc/mean_oo  #normalized mean for drCCA
           


   if(pca)
   {
      cross_pca <- pwisepca(mat,dim)

      mean_pc <- mean(cross_pca[cross_pca !=0]) # for PCA
      mpc <- mean_pc/mean_oo  #normalized mean for PCA

      result <- list(oo = orig, cc= cross_cca, pc = cross_pca, mcca =mcc, mpca=mpc)            

      return(result)
   }else{ return(list(oo = orig, cc =cross_cca, mcca =mcc))}

}



  



