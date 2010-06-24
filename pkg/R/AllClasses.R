setClass("DependencyModel", representation(W = "list", phi = "list", score = "numeric", loc = "numeric",
		featureName = "character", windowSize = "numeric", method = "character", params = "list"))


setClass("DependencyScreenModels", representation(models = "list", windowSize = "numeric", method = "character", params = "list"))


#setClass("SparseDependencyModel", contains = "DependencyModel")


