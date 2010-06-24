
##########################################################################
# Regularized CCA
# We will drop those eigen values which contributes less than a threshold 
# in the frobenius norm of the matrix
#########################################################################

###INPUTS###

#input 1. list of  data matrices
#input 2. regularization parameter (0 to 1), default is 0

# regularization <- a value p  will drop smaller eigen values 
# contributing less than p% of total frobenius norm or 1-norm
# "A value of zero mean no regularization chosen"

###OUTPUT VALUES###

# A list containing three components
# 1. x$eigval <- canonical correlations 
# 2. x$eigvecs <- canonical variates
# 3. x$proj <- Projected data on canonical variates
# 4. x$meanvec <- array of columnwise means for each data matrix
# 5. x$white <- array of whitening matrices for each data matrix

 "regCCA" <- 
 function(datasets,reg = 0)
{
 
    mat <- datasets #list of data matrices
  
    m <- length(mat)
  
    #print('total number of files\n')
    #print(m);

    para <- reg # between 0 and 1

    if(para < 0 || para > 1)
    stop(" Value of 'reg' should be greater than zero and less than one")


##################################################
####SUBROUTINES####


drop <- function(mat,parameter)
{
      mat <- as.matrix(mat)

      p <- parameter

      fac <- svd(mat)

      len <- length(fac$d)
 
      ord <- order(fac$d) #sort in increasing order,ie smallest eig val first
   
      ord_eig <- fac$d[ord] #order eig vals in increasing order


     # cumulative percentage of Frobenius norm
     
     portion <- cumsum(ord_eig^2)/sum(ord_eig^2)

     #portion <- cumsum(ord_eig)/sum(ord_eig)


     ind <- which(portion < p) #ind r those who contribute less than p%


     #drop eig vecs corrsponding to small eig vals dropped
       num <- len -length(ind)
  
      fac$v <- as.matrix(fac$v[,1:num])

      fac$u <- as.matrix(fac$u[,1:num])

   mat_res <- (fac$v %*% sqrt(diag(1/fac$d[1:num], nrow = length(fac$d[1:num])))%*% t(fac$u))

     w <- list(mat_res = mat_res, num = num);

     return(w);
}



########SUBROUTINES END HERE
##################################################
    # mat is pointer to array of matrices

       covm <- array(list(),m) #array of covariance matrices

       fea <- c(rep(0,m)) # numebr of dimensions in each data matrix

      mean_m <- array(list(),m) # array of columnwise mean vector for each mat

        for(i in 1:m)
         {

          mean_m[[i]] <- apply(mat[[i]],2,mean)

           #subtracting columnwise mean
           mat[[i]] <- apply(mat[[i]],2,function(x){x - mean(x)})

           
           covm[[i]]<- cov(mat[[i]]); #covariance matrix

           fea[i] <- ncol(mat[[i]])
         }


# whitening of data...Regularized or Normal

     whiten_mat <- array(list(),m) # array of matrices after whitening

     cov_whiten <- array(list(),m) # covariance of whitened matrices

      white <- array(list(),m) #array of matrices,which does whitening

     d <- c(rep(0,m)) # to store the number of dimensions we kept in each
                      # data matrix after regularization

    # approximation of covariance matrices for regularization


      if(para > 0)
      {
             
               for(i in 1:m)
             {

               dummy <- drop(covm[[i]], para)

              whiten_mat[[i]] <- mat[[i]] %*% dummy$mat_res 
              #  print(det(cov(whiten_mat[[i]])))
              d[i] <- dummy$num #dimensions kept for whitening

              white[[i]] <- dummy$mat_res #storing whitening matrix
              }
      }else
      {
 
           for(i in 1:m)
           {
            whiten_mat[[i]] <- mat[[i]] %*% SqrtInvMat(covm[[i]])

            d[i] <- ncol(mat[[i]]) #all dimensions, no whitening

            white[[i]] <- SqrtInvMat(covm[[i]])

            #print('dim of sqrt_mat inverse mat')
            # print(det(cov(whiten_mat[[i]])))
            #print(dim(sqrt_mat_inv(covm[[i]])))
           }
      }


            if(sum(d) ==  sum(fea))
            {
                 print('Normal gCCA')
            }else{print('Regularized gCCA')}

            #print(fea)
            #print(d)
           # print(sum(d))

 
    #concactenating the whitened matrices

           a <- concatenate(whiten_mat) #internal function
   
           #print('concatenate')
           #print(dim(a))


   # z is covariance matrix of concatenated data

           z <- cov(a)
       
           eig <- svd(z) #use svd, eigen is not so good here

        # projected data

           proj_data <- a %*% eig$v


##-------------------------------------------##

     # whitening matrix %*% eig$v

       
     eig_wh <- array(list(),m) # whitened eig vecs for eaxh matrix

      if(para > 0)
      {
             
               j <- 0

               for(i in 1:m)
             {

               dummy <- drop(covm[[i]], para)

              eig_wh[[i]] <- dummy$mat_res %*% eig$v[(j+1) : (j+fea[i]),]

              j <- j + fea[i]
              #print(dim(eig_wh[[i]]))

              }
      }else
      {
 
           j <- 0
           for(i in 1:m)
           {
            
            eig_wh[[i]] <- SqrtInvMat(covm[[i]]) %*% eig$v[(j+1) : (j+fea[i]),]

              j <- j + fea[i]
              
            }
      }


       
        #fullwh <- eig_wh[[1]] #rbind all eig_wh
        #if(m >1)
        #{
           #for(i in 2:m)
           #{
            # fullwh <- rbind(fullwh,eig_wh[[i]])           
 
           #}
        #}

        #fullmat <- concatenate(mat)

        #wproj <- fullmat %*% fullwh
          
            


##-------------------------------------------##





     # make list of eigen values, eigen vectors and projected data

      gencorr <- list(eigval = eig$d,eigvecs = eig_wh, proj = proj_data, meanvec = mean_m, white = white)
      
      return(gencorr)

   

}

