
Pattern component modelling toolbox
===================================

Jörn Diedrichsen, Spencer Arbuckle, and Atusushi Yokoi

# Overview 
The pattern-component modelling (PCM) toolbox provides function for model fitting, comparison, and some basic visualization. Users can use these provided functions to test various hypotheses about how the brain represents features of external world or actions. This readme page just briefly summarizes the contents of the toolbox and short comments. The complete-version of this manual accompanying with some example use of the codes and full-details of mathematical derivations is provided [here](https://github.com/jdiedrichsen/pcm_toolbox/blob/master/pcm_toolbox_manual.pdf). Also, some scientific background, discussions, and relationship between other relevant approaches can be found [here](http://biorxiv.org/content/early/2017/03/25/120584) and [there](http://biorxiv.org/content/early/2017/04/02/071472). 

What the toolbox does *not* provide are functions to extract the required data from the first-level GLM or raw data, search-light code or code for ROI definition. We have omitted these function as they strongly depend on the analysis package used for basic imaging analysis. Some useful tools for the extraction of multivariate data from the standard first-level GLM, please see the [RSA-toolbox](https://github.com/rsagroup/rsatoolbox) and [Surfing toolbox](https://github.com/nno/surfing). 

The functions provided in the toolbox can be categorized into different categories:  

### Basic likelihood and optimization

These are the functions that perform the core statistical functions of the toolbox.  

| 	Function 			    | Comment  
|:-----------------------------------|:-----------------------------
|  `pcm_likelihood` |  Likelihood of a single data set under a model
|  `pcm_likelihoodIndivid` | pcm_likelihood with optional random or fixed block effect
|  `pcm_likelihoodGroup`   | Likelihood of a group data set under a model
|  `pcm_NR`				| Newton Raphson optimisation 
|  `pcm_minimize`			| Conjugate gradient descent 

### Model Evaluation
These functions are higher level functions that perform fitting and crossvalidation of either individual data set or group data sets.  

| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_fitModelIndivid`        | Fits G and noise parameter to individual data
| `pcm_fitModelIndividCrossval`| Within-subject crossvalidation of models
| `pcm_fitModelGroup`          | Fit common G to all subjects
| `pcm_fitModelGroupCrossval`  | Between-subject crossvalidation of models  
| `pcm_setUpFit`         |  Generally prepares models and data for fitting 
| `pcm_knockModels`     | Inference on components using simple knock-in knock-out likelihoods
| `pcm_componentPosterior`     | Inference on components by model averaging 

### Visualization functions
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_classicalMDS`           | Multidimensional scaling  
| `pcm_estimateU`         | Estimates voxel-patterns under model M
| `pcm_estimateW`         | Estimates  voxel-feature weights under  model M
| `pcm_plotModelLikelihood`    | Displays marginal likelihood for different models

### Model building 
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_constructModelFamily`     | Makes a family of models from different components
| `pcm_buildModelFromFeatures`     | Makes a component model from featuresets


### Utility functions
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_addModelComp`                  | Adds a model component to model M  
| `pcm_blockdiag`                  | Makes a blockdiagonal matrix  
| `pcm_calculateG`                  | Determines G for models   
| `pcm_checkderiv`                  | Checks derivate of a nonlinear model 
| `pcm_diagonalize`                  | Decomposes G into A * A'  
| `pcm_estGcrossval`              | Cross-validated estimate of G 
| `pcm_generateData`                  | Generates data simulations from model M  
| `pcm_getStartingval`               | Provides a starting value estimate for model M  
| `pcm_indicatorMatrix`       | Generates indicator matrices 
| `pcm_vararginoptions`	    | Deals with variable input options 
| `pcm_getUserOptions`	    | Deals with variable input options (struct-based, faster) 
| `pcm_makePD`    | Moves a matrix to closest semi-postive definite Matrix 
| `pcm_optimalAlgorithm`    | Recommendation for the best optimisation algorithm
| `pcm_prepFreeModel`         |  Sets up free model (freechol)


### Recipes 
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_recipe_finger`          | Example of a fixed and component models 
| `pcm_recipe_correlation`     | Example of feature model   
| `pcm_recipe_nonlinear`       | Example for non-linear model   

### Currently not maintained / under developement
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
|  `pcm_NR_diag`			| Newton Raphson for diagonalized component models
|  `pcm_NR_free`			| Newton Raphson for a free model 
|  `pcm_EM`				| Expectation-Maximization 