#setReplaceMethod(f="setLoc", signature("DependencyModel"),#
#	definition=(function(model,value) {#
#		model@loc <- value
#		return(model)
#	}
#))


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

