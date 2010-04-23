rm(list=ls()) # Empty the workspace from other variables.
options(error=recover) # debug mode

path = "/share/mi/workspace/tsuvitai/mi/CCA/package/package-100326/" # This is the root folder of the package. The code should be located in folder "sourcecode" and the results will be saved into sub-folders of the folder "results".
runId = "twowayCCA-100326-gen-demo" # Keep this ID unique for each simulation. The program will save the results into a folder named after this.

NburnIn = 1000 # length of the burn-in; remember to ensure sufficient burn-in
Niterfinal <- 1000 # number of iterations; remember to ensure sufficient number of samples from the posterior

dataset = 0 # 0/1: use generated data (0) or real data (1)
takeLog = FALSE # TRUE/FALSE: do/do not perform a log-transform on real data

nXlat <- 3 # dimensionality of the latent variable of the first data set in learning
nYlat <- 3 # dimensionality of the latent variable of the second data set in learning

doPlotting = TRUE # TRUE/FALSE: do/do not perform plotting of results

source(paste(path,"sourcecode/loadSource.R",sep="")) # Load the function that loads other source code. Do not change this line.
loadSource(path=path) # Load the source code. Do not change this line.

if (dataset==1) { ## Load and process real data. Write your own data loading procedure here!
	data = list() # 'data' is a list-structure containing the observation matrices 'X' and 'Y', and covariate vectors 'a', 'b' and 'c'. The number of columns in matrices 'X' and 'Y' must agree, as the samples are assumed to be paired. Do not edit this line.
	data$X = as.matrix(...) # 'data' is a MxN matrix (features 'm' as rows, samples 'n' as columns). Write your own data loading procedure here. 
	data$Y = as.matrix(...) # 'data' is a MxN matrix (features 'm' as rows, samples 'n' as columns). Write your own data loading procedure here. 
	data$covariates = list()
	# Load covariate values into vectors 'a', 'b' and 'c'. Set non-used covariates to zero vectors.
	data$covariates$a = rep(0,ncol(data)) # 'a' is a two-level covariate: allowed values are 0 and 1.
	data$covariates$b = rep(0,ncol(data)) # 'b' is a two-level covariate: allowed values are 0 and 1.
	data$covariates$c = rep(0,ncol(data)) # 'c' is a multi-level covariate: allowed values are 1, 2, 3, ... If the covariate is not used, set all values to 0.
	 # Optional character vectors containing the names of the features (e.g. 'c("Name1", "Name2", "Name3")'). The length must be equal to the number of features in the corresponding data matrix or 'NULL'.
	data$Xnames = NULL # can be edited
	data$Ynames = NULL # can be edited

} else if (dataset==0) { ## Generate simulated two-way data from the model
	# The following lines generate a simulated two-way, two-source, data set. Parameters can be changed by the user, unless otherwise stated.
	
	## Generated data parameters
	N <- 100 # number of samples in the generated data; can be edited
	nX = 200 # number of features in view X; can be edited
	nY = 210 # number of features in view Y; can be edited

	genParams = list() # Do not edit this line.
	genParams$Xnoise = 1 # noise of the data set 'X'; can be edited
	genParams$Ynoise = genParams$Xnoise # noise of the data set 'Y'; can be edited
	
	## Projection matrix 'W'
	# Proper CCA data generation: data sets X and Y have a shared latent variable 'z'.
	fixedW = TRUE # TRUE/FALSE: do/do not use fixed projection from z to xlat and ylat on generated data; can be edited
	if (fixedW) { # Use fixed projection 'W' from 'z' to 'xlat' and 'ylat'.
		Wgen = list() # Do not edit this line.
		Wgen$X = matrix(c(1,0,-1, 1,0,0, 0,0,0),ncol=3,byrow=F) # projection of 'z' to 'xlat'; can be edited
		Wgen$Y = matrix(c(1,0,-1, 0,0,0, 0,1,0),ncol=3,byrow=F) # projection of 'z' to 'ylat'; can be edited
	} else { # Use random projection.
		Wgen = NULL # Do not edit this line.
	}
	
	## Covariates
	# Covariate vectors can be edited by the user but the vectors must equal to the number of samples.
	# The following code generates two-way data from two separate sources.
	treatments = list() # Do not edit this line.
	treatments$a = rep(c(0,1),length.out=N) # 'a' is a two-level covariate: allowed values are 0 and 1. Can be edited.
	treatments$b = rep(c(0,0,1,1),length.out=N) # 'b' is a two-level covariate: allowed values are 0 and 1. Can be edited.
	
	## Covariate effects
	# Two-way effects of latent variable 'z'. Effects are ordered such that three consecutive numbers correspond to effects 'a', 'b' and '(ab)' of one component in latent variable 'z'. In other words, the groups of three values correspond to shared effects (first three), X-specific effects (middle three), and Y-specific effects (last three). Warning: effects here are in wrong order compared to 'sampleEffects'!
	effects = array(c(2,0,0, 0,0,2, 0,2,0),dim=c(3,3)) # Can be edited.
	
	## Generate data
	data = gener2wayCCAdata(N_samples=N, nZ=3, xdim=nXlat, ydim=nYlat, xxdim=nX, yydim=nY, W=Wgen, treatments=treatments, effects=effects, params=genParams) # Do not edit this line.

} ## End of simulated data generation.

## Run the simulation.

posterior = multiWayCCA(data=data, nXlat=nXlat, nYlat=nYlat, takeLog=takeLog, maxBurnIn=NburnIn, Niterfinal=Niterfinal, path=path, runId=runId) # Do not edit this line.