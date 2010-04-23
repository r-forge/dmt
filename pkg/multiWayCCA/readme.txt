Multivariate multi-way analysis of multi-source data
Ilkka Huopaniemi, Tommi Suvitaival, Janne Nikkilä, Matej Orešič, and Samuel Kaski. Multivariate multi-way analysis of multi source data. Bioinformatics, To appear. (ISMB 2010).

Documentation for the R implementation 'multiWayCCA'.

26.03.2010

Tommi Suvitaival, tommi.suvitaival@tkk.fi
Ilkka Huopaniemi, ilkka.huopaniemi@tkk.fi

Aalto University
Department of Computer and Information Science


This package includes:
-source code for the R implementation of multi-way, multi-view analysis model
-a script for carrying out the analysis

Description for the model is available at http://www.cis.hut.fi/projects/mi/software/multiWayCCA/ .

The implementation is based on Arto Klami's Bayesian CCA implementation in R. Inference is done with Gibbs sampling.
R is available at http://www.r-project.org/ .


SCRIPTS

3-way, 2-source, analysis
-script file multiWayCCA-script-100326.R


Steps in preparing the package for new experiments
1. write a procedure for loading your data
2. set sampling parameters to suit your interests

Simulated data

The package includes a routine for generating simulated normally distributed data. The user will have to define, what covariate effects are generated. This is done by editing the matrix "effects". The effects are generated on latent variables 'z', where the first component is shared between all data sets and the rest are data set-specific, one per each data set. The first column of the matrix "effects" includes the effects of the shared latent component. Rest of the columns are the effects of the data set-specific latent components (2 for 2-source).

The latent variables 'z' are projected onto data set-specific latent variables 'xlat' and 'ylat'. Note that these are a different thing than the data set-specific components of the latent variable 'z'. The number of components in these data set'specific latent variables is selectable by the user. These latent variables are projected onto higher dimensions to generate the observed data X and Y.

Real data

The user should write a routine for loading your data. The function provided in this package is not universal and works only for the type of data files that were used in the development of this method. To write the routine, edit the lines after "Load and process real data". You should be apple to supply the method with the data matrices, the number of features in data matrices, the covariates ("a", "b", "c") and the names of the features.

The covariates are called "a", "b" and "c". Covariates "a" and "b" have been implemented to have two levels, 0 and 1. Covariate "c" has been implemented to have S levels but more than two levels have not been tested in the 3-way setting.

In the script, the data is Z-normalized (zero mean, unit variance) according to the base-level (the population, where the covariates are zero).

Inference of the model is achieved with standard Gibbs sampling. The implementation does not include convergence control even though the code has some hints about that. Thus, the user should pay attention to allowing long enough time for the burn-in.



CONTENTS OF A SCRIPT

First, the user will need to define the file path ('path') of the package. The package requires the source code to be located in directory 'path/sourcecode/'. When path is successfully defined, required source code files are loaded. The package also requires some additional R libraries, such as 'mvtnorm'.

PARAMETERS

Parameters for running the analysis are defined after the source files. The user can modify these parameters to e.g. set the number of latent components in the model.

'path' defines the root folder to which results from individual analysis are saved. Also the script file is copied into the results directory from folder 'path/scripts'.
'runId' is a unique name for the analysis. Results of the analysis are saved into a folder named according to this variable. To avoid deletion of previous results, the user should rename this variable every time he runs an analysis. The script file is assumed to be named 'runid.R' and the user should follow this convention, as then the script file will be correctly copied to the results directory for later reference.

'NburnIn' is the length of the burn-in phase for the MCMC chain
'Niterfinal' is the number of Gibbs samples to be drawn after the burn-in

'dataset' is a selector for the data type to be analyzed. Value 0 will lead to generation of data from the model. Value 1 will lead to loading of external data. Typically, value 1 is in use.
'takeLog' decides whether a logarighm will be taken of the data before the actual analysis. It is recommended if the raw data is not normally distributed. This parameter applies to real data only.

'nXlat' and 'nYlat' are the numbers of view-specific latent components. The data features will be clustered in to a number of clusters defined by these values.

'doPlotting' decides whether the results are saved both numerically and visually (TRUE) or whether only numerical samples are saved (FALSE).

DATA

The user should write a new data import function. Additionally, the covariates 'a', 'b' and 'c' need to be served. The function should return the data as a list with entries 'X' and 'Y' corresponding to view observation matrices, and vectors 'a', 'b' and 'c' corresponding to the three covariates.

Number of columns ('N') is required to be the same for all views. Also covariate vectors are required to be of length 'N'. Number of features (rows) in views can vary.

When generated data is determined to be used (by variable 'dataset'), dimensionality of the data needs to be defined.
'N' is the number of samples (observations). 'nX' and 'nY' are the numbers of features in views.
'genParams$Xnoise' is the noise level in first view. Similarly for other views.
'Wgen$X' is the projection matrix from "higher level" latent variable 'z' to view-specific latent variable 'xLat'. This applies only if variable 'fixedW' is set to 'TRUE'. A column of the matrix corresponds to weights from one component of 'z'. Thus, when a column is set to zero, the projection will be completely independent of the corresponding component of 'z'. Row 'i' of the matrix corresponds to a projection weights from latent variable 'z' to component 'i' of 'xLat'. Similarly for other views.
'treatments' is a list of covariates ('case', 'gender', 'state') that can be defined by the user. Covariates, whose effects are not estimated, can be set to zero, still conserving the correct length of the vector.
'effects' is a matrix of covariates effects to be generated. Rows of the matrix correspond to different covariate effects and are in order 'a', 'b', 'ab'. Columns correspond to components of latent variable 'z', i.e. into which component the effect is generated. If non-zero values of the matrix mean that an effect of that size is generated to the corresponding population.
The data is generated using functions 'gener2wayCCAdata' for 2-way, 2-view, experiments.

SAMPLING

First stage of sampling is the burn-in, which ensures that the sampler converges to the proper posterior distribution. The user should pay attention to allow for long enough burn-in. Burn is followed by additional sampling from the posterior distribution. Posterior samples are saved into list 'posterior'. Gibbs samples are saved into a .Rdata file named after 'runId' and placed into the folder 'path'.

For generated data experiment, view-specific clusterings are compared to the real clusterings that were used to generate the data. The clusters in the posterior samples are ordered according to this comparison.

If parameter 'doPlotting' is set to value 'TRUE', visual plotting of the marginal posterior distributions is done.



FURTHER REMARKS

If the user is interested in less-than-three-way analysis, additional covariates should be set to zero such that they are of correct length.
