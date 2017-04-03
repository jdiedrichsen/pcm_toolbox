
Pattern component modelling toolbox
===================================

Jörn Diedrichsen, Spencer Arbuckle, and Atusushi Yokoi

# Overview 
The pattern-component modelling (PCM) toolbox provides function for model fitting, comparison, and some basic visualization. Users can use these provided functions to test various hypotheses about how the brain represents features of external world or actions. This readme page just briefly summarizes the contents of the toolbox and short comments. The complete-version of this manual accompanying with some example use of the codes and full-details of mathematical derivations is provided [here](/url/to/complete_manual). Also, some scientific background, discussions, and relationship between other relevant approaches can be found [here](http://biorxiv.org/content/early/2017/03/25/120584) and [there](http://biorxiv.org/content/early/2017/04/02/071472). 

What the toolbox does *not* provide are functions to extract the required data from the first-level GLM or raw data, search-light code or code for ROI definition. We have omitted these function as they strongly depend on the analysis package used for basic imaging analysis. Some useful tools for the extraction of multivariate data from the standard first-level GLM, please see the [RSA-toolbox](https://github.com/rsagroup/rsatoolbox) and [Surfing toolbox](https://github.com/nno/surfing). 

The functions provided in the toolbox can be categorized into different categories:  

## Basic likelihood and optimization

These are the functions that perform the core statistical functions of the toolbox.  

| 	Function 			    | Comment  
|:--------------------------|:-----------------------------
|	`pcm_likelihoodIndivid` | Likelihood of a single data set under a model
|   `pcm_likelihoodGroup`   | Likelihood of a group data set under a model
|   `pcm_NR`				| Newton Raphson optimisation 
|   `pcm_NR_diag`			| Newton Raphson for diagonalized models (faster)
|   `pcm_NR_free`			| Newton Raphson for a free model 
|   `pcm_EM`				| Expectation-Maximization 
|   `pcm_minimize`			| Conjugate gradient descent 

## Model Evaluation
These functions are higher level functions that perform fitting and crossvalidation of either individual data set or group data sets.  

| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_fitModelIndivid`        | Fits G and noise parameter to individual data
| `pcm_fitModelIndividCrossval`| Within-subject crossvalidation of models 
| `pcm_fitModelGroup`          | Fit common G to all subjects, using individual noise and scale parameter
| `pcm_fitModelGroupCrossval`  | Between-subject crossvalidation of models  


## Utility functions
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_checkderiv`             | Checks derivate of a nonlinear model 
| `pcm_estGcrossval`           | Crosscvalidated estimate of G 
| `pcm_indicatorMatrix`        | Generates indicator matrices 
| `pcm_vararginoptions`		   | Handles optional input arguments

## Visualization functions
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_classicalMDS`           | Multidimensional scaling  
| `pcm_plotModelLikelihood`    | Displays marginal likelihood (results from `pcm_fitModel` functions)
| `pcm_plotFittedG`            | Displays second-moment matrices with fitted parameters

## Recipes 
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_recipe_finger`          | Example of fixed ad component models 
| `pcm_recipe_correlation`     | Example of feature model  
| `pcm_recipe_nonlinear`       | Example for non-linear model   
