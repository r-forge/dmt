# TODO: move application-specific functions to other classes, keep only model-specific stuff


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

#setMethod("getLoc","DependencyModel", 
#	function(model) { 
#		return(model@loc) 
#	} 
#) 

setMethod("getFeatureName","DependencyModel", 
          function(model) { 
            return(model@featureName) 
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
		return(model@windowSize) 
	} 
) 

