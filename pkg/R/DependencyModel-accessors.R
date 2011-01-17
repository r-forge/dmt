setMethod("getW","DependencyModel", 
  function(model) { 
    return(model@W) 
  } 
) 

setMethod("getPhi","DependencyModel", 
  function(model) { 
    return(model@phi) 
  } 
) 

setMethod("getScore","DependencyModel", 
  function(model) { 
    return(model@score) 
  } 
) 

setMethod("getParams","DependencyModel", 
  function(model) { 
    return(model@params) 
  } 
) 

setMethod("getModelMethod","DependencyModel", 
  function(model) { 
    return(model@method) 
  } 
) 

setMethod("getWindowSize","DependencyModel", 
  function(model) { 
    if(is.null(getW(model)$X)){
      return(dim(getW(model))[1])
    }
    else {
      return(c(dim(getW(model)$X)[1],dim(getW(model)$Y)[1]))
    }
  } 
) 
