setMethod(f="[[", signature("DependencyScreenModels"),
			definition=(function(x,i,j,drop) {
				return(x@models[[i]])
			} 
))

setReplaceMethod(f="[[",signature("DependencyScreenModels"),
				definition=(function(x,i,j,value) {
					x@models[[i]] <- value
 					return(x)
 				}
))
	
setMethod("getScore","DependencyScreenModels", 
	function(model) {
		scores <- vector()
		for (i in seq_along(model@models)) {
				scores[i] <- getScore(model[[i]])
		}
		return(scores) 
	} 
)			

#setMethod("getLoc","DependencyScreenModels", 
#	function(model) {
#		locs <- vector()
#		for (i in seq_along(model@models)) {
#			locs[i] <- getLoc(model[[i]])
#		}
#		return(locs) 
#	} 
#)			

setMethod("getFeatureName","DependencyScreenModels", 
	function(model) {
          featureNames <- vector()
          for (i in seq_along(model@models)){
            featureNames[i] <- getFeatureName(model[[i]])
          }
          return(featureNames) 
	} 
)	

setMethod("getModelMethod","DependencyScreenModels", 
	function(model) { 
		return(model@method) 	
	} 
)

setMethod("getParams","DependencyScreenModels", 
	function(model) { 
		return(model@params) 
	} 
)

setMethod("getModelNumbers","DependencyScreenModels", 
	function(model) { 
		return(length(model@models)) 
	} 
)

setMethod("getWindowSize","DependencyScreenModels", 
	function(model) { 
		return(model@windowSize) 
	} 
)

setMethod("isEmpty","DependencyScreenModels",
	function(model) {
		return(length(model@models) == 0)
	}
)

#setMethod("topFeatures","DependencyScreenModels",
#	function(model, num = 1) {
#          scores <- getScore(model)
#          data <- data.frame(scores)
#          # order dataframe and take num names of genes with highest scores 
#          return(as.character(data[order(scores,decreasing=TRUE),]$genes[1:num]))	
#	}
#          
#)

setMethod("topModels","DependencyScreenModels",
	function(model,num = 1) {
		scores <- getScore(model)
		indices <- 1:length(scores)
		data <- data.frame(scores,indices)
		# Order dataframe
		data <- data[order(scores,decreasing=TRUE),]
		returnList = list()
		for(i in 1:num){
			returnList = c(returnList,model[[data$indices[i]]])
		}
		return(returnList)
	}
)

#setMethod("orderGenes","DependencyScreenModels",
# function(model){
#
#    scores <- getScore(model)
#    genes <- getGeneName(model)
#    data <- data.frame(scores,genes)
#    data <- data[order(scores,decreasing=TRUE),]
#    return(as.character(data$genes))
#  }	
#)

#setMethod("findModel","DependencyScreenModels",
#  function(model, name){
#    index = which(getGeneName(model) == name)
#    if (length(index) > 0)
#      return(model[[which(getGeneName(model) == name)]])
#    stop("No model found")
#  }
#)
    
