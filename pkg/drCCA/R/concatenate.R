
#internal subroutine extra

"concatenate" <- 
function(datasets)
{

   mat <- datasets #list of data sets

  m <- length(mat)

  com <- mat[[1]]

  if(m > 1)
  {
       for(i in 2:m)
       {
         com <- cbind(com,mat[[i]])
       }
  }


  return(com)
}
