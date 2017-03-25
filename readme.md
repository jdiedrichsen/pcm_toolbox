
Pattern component modelling toolbox
===================================

Jörn Diedrichsen, Spencer Arbuckle, and Atusushi Yokoi

# How to use this manual
This manual provides an introduction to how to use the Pattern component modelling (PCM) toolbox. The theory behind this approach is laid out in an accompanying paper (REF) - but the main ideas are sketched out here in the introduction. We then provide an overview over the toolbox functions, and explain the different steps of model specification, model estimation, visualisation, and model comparison following real examples presented in the paper. The toolbox comes with a few example "recipes", which we hope will be useful to provide a full example of usage. Finally, the last section contains some of the mathematical details readers that would like to understand all algorithmic details.  


# Introduction 
The study of brain representations aims to illuminate the relationship between complex patterns of activity occurring in the brain and "things in the world" - be it objects, actions, or abstract concepts. By understanding internal syntax of brain representations, and especially how the structure of representations changes across different brain regions, we ultimately hope to gain insight into the way the brain processes information.

Central to the definition of representation is the concept of decoding \citep{RN3623}. A feature (i.e. a variable that describes some aspect of the "things in the world") that can be decoded from the ongoing neural activity in a region is said to be represented there. For example, a feature could be the direction of a movement, the orientation and location of a visual stimulus, or the semantic meaning of a word. Of course, if we allow the decoder to be arbitrarily complex, we would use the term representation in the most general sense. For example, using a computer vision algorithm, one may be able to identify objects based on activity in primary visual cortex. However, we may not conclude necessarily that object identity is represented in V1 - at least not explicitly. Therefore, it makes sense to restrict our definition of an explicit representation to features that can be linearly decoded by a single neuron from some population activity \citep{RN3627, RN3626, RN3625, RN3573}.

While decoding approaches are very popular in the study of multi-voxel activity patterns \citep{RN2832, RN2839, RN2909}, they are not the most useful tool when we aim to make inferences about the nature of brain representations. The fact that we can decode feature X well from region A does not imply that the representation in A is well characterized by feature X - there may be many other features that better determine the activity patterns in this region. 

![Figure 1](Figures/Figure_1.png)
**Figure 1.** 
*Decoding, encoding and representational models. (**A**) The matrix of activity data consists of rows of activity patterns for each condition or of columns of activity profiles for each voxel (or more generally, measurement channel). The data can be used to decode specific features that describe the experimental conditions (decoding). Alternatively, a set of features can be used to predict the activity data (encoding). Representational models work at the level of a sufficient statistics (the second moment) of the activity profiles. Models are formulated in this space and possibly combined and changed using higher-order model parameters ($\theta$). (**B**) The activity profiles of different voxels are plotted as points in the space of the experimental conditions. Features in encoding models are vectors that describe the overall distribution of the activity profiles. (**C**) The distribution can also be directly described using a multivariate normal distribution (PCM).(**D**) Representational similarity analysis (RSA) provides an alternative view by plotting the activity patterns in the space defined by different voxel activities. The distances between activity patterns serves here as the sufficient statistic, which is fully defined by the second moment matrix.*

Encoding models, on the other hand, characterize how well we can explain the activities in a specific region using a sets of features. The activity profile of each voxel (here shown as columns in the activity data matrix), is modeled as the linear combination of a set of features (\hyperref[fig1]{Fig. 1a}). We will use the term voxels interchangeably with the more general term measurement channel, which could, depending on the measurement modality, refer to a single neuron, an electrode, or sensor. Each voxel has its own set of parameters ($\mathbf{W}$) that determine the weight of each feature. This can visualized by plotting the activity profile of each voxel into the space spanned by the experimental conditions (\hyperref[fig1]{Fig. 1b}). Each dot refers to the activity profile of a channel (here a voxel), indicating how strongly the voxel is activated by each condition. Estimating the weights is equivalent to a projection of each of the activity profiles onto the feature vectors. The quality of the model can then be evaluated by determining how well unseen activity data can be predicted. When estimating the weights, encoding models often use some form of regularization, which essentially imposes a prior on the feature weights. This prior is an important component of the model. It determines a predicted distribution of the activity profiles \citep{RN3573}. An encoding model that matches the real distribution of activity profiles best will show the best prediction performance.

The interpretational problem for encoding models is that for each feature set that predicts the data well, there is an infinite number of other (rotated) features sets that describe the same distribution of activity profiles and, hence, predict the data equally well. The argument may be made that to understand brain representations, we should not think about specific features that are encoded, but rather about the distribution of activity profiles. This can be justified by considering a read-out neuron that receives input from a population of neurons. From the standpoint of this neuron, it does not matter which neuron has which activity profile (as long as it can adjust input weights), and which features were chosen to describe these activity profiles - all that matters is what information can read out from the code. Thus, from this perspective it may be argued that the formulation of specific feature sets and the fitting of feature weights for each voxel are unnecessary distractions. 

Therefore, our approach of pattern component modeling (PCM) abstracts from specific activity patterns. This is done by summarizing the data using a suitable summary statistic (\hyperref[fig1]{Fig. 1a}), that describes the shape of the activity profile distribution (\hyperref[fig1]{Fig. 1c}). This critical characteristic of the distribution is the covariance matrix of the activity profile distribution or - more generally -  the second moment. The second moment determines how well we can linearly decode any feature from the data. If, for example, activity measured for two experimental conditions is highly correlated in all voxels, then the difference between these two conditions will be very difficult to decode. If however, the activities are uncorrelated, then decoding will be very easy. Thus, the second moment is a central statistical quantity that determines the representational content of the brain activity patterns of an area \citep{RN3573}.

Similarly, a representational model is formulated in PCM not by its specific feature set, but by its predicted second moment matrix. If two feature sets have the same second moment matrix , then the two models are equivalent. Thus, PCM makes hidden equivalences between encoding models explicit. To evaluate models, PCM simply compares the likelihood of the data under the distribution predicted by the model. To do so, we rely on an generative model of brain activity data, which fully specifies the distribution and relationship between the random variables. Specifically, true activity profiles are assumed to have a multivariate Gaussian distribution and the noise is also assumed to be Gaussian, with known covariance structure. Having a fully-specified generative model allows us to calculate the likelihood of data under the model, averaged over all possible values of the feature weights. This results in the so-called model evidence, which can be used to compare different models directly, even if they have different numbers of features. In summarizing the data using a sufficient statistic, PCM is closely linked to representation similarity analysis (RSA), which characterizes the second moment of the activity profiles in terms of the distances between activity patterns (\hyperref[fig1]{Fig. 1d}, also see \citealt{RN3573}).

By removing the requirement to fit and cross-validate individual voxel weights, PCM enables the user to concentrate on a different kind of free parameter, namely model parameters that determine the shape of the distribution of activity profiles. From the perspective of encoding models, these would be hyper-parameters that change the form of the feature or regression matrix. For example, we can fit the distribution of activity profiles using a weighted combination of 3 different feature sets (\hyperref[fig1]{Fig. 1a}). Such component models (see section \ref{componentmodels}) are extremely useful if we hypothesize that a region cares about different groups of features (i.e.\ colour, size, orientation), but we do not know how strongly each feature is represented. In encoding models, this would be equivalent to providing a separate scaling factor to different parts of the feature matrix. Most encoding models, however, use a single model feature matrix, making them equivalent to a fixed PCM model.

In this manual we will show how to use the PCM toolbox to estimate and compare flexible representational models. We will present the fundamentals of the generative approach taken in PCM and outline different ways in which flexible representational models with free parameters can be specified. We will then discuss methods for model fitting and for model evaluation. We will also walk in detail through three illustrative examples from our work on finger representations in primary sensory and motor cortices, also providing recipe code for the examples presented in the paper. 

# Overview 
the toolbox provides function for model fitting, comparison, and some basic visualization. What the toolbox does *not* provide are functions to extract the required data from the first-level GLM or raw data, search-light code or code for ROI definition. We have omitted these function as they strongly depend on the analysis package used for basic imaging analysis. Some useful tools for the extraction of multivariate data from the standard first-level GLM, please see the RSA-toolbox and Surfing toolbox. 

The functions provided in the toolbox can be categorized into different categories:  

### Basic likelihood and optimization

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

### Model Evaluation
These functions are higher level functions that perform fitting and crossvalidation of either individual data set or group data sets.  

| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_fitModelIndivid`        | Fits G and noise parameter to individual data
| `pcm_fitModelIndividCrossval`| Within-subject crossvalidation of models 
| `pcm_fitModelGroup`          | Fit common G to all subjects, using individual noise and scale parameter
| `pcm_fitModelGroupCrossval`  | Between-subject crossvalidation of models  


### Utility functions
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_checkderiv`             | Checks derivate of a nonlinear model 
| `pcm_estGcrossval`           | Crosscvalidated estimate of G 
| `pcm_indicatorMatrix`        | Generates indicator matrices 
| `pcm_vararginoptions`		   | 

### Visualization functions
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_classicalMDS`           | Multidimensional scaling  

### Recipes 
| 	Function 			       | Comment  
|:-----------------------------|:-----------------------------
| `pcm_recipe_finger`          | Example of fixed ad component models 
| `pcm_recipe_correlation`     | Example of feature model  
| `pcm_recipe_nonlinear`       | Example for non-linear model   
 
# Model specification 

## Generative model 

Central to PCM is a generative model of the measured brain activity data $\mathbf{Y}$, a matrix of N x P activity measurements, referring to N time points (or trials) and P voxels. The data can refer to the minimally preprocessed raw activity data, or to already deconvolved activity estimates, such as those obtained as beta weights from a first-level time series model. $\mathbf{U}$ is the matrix of true activity patterns (a number of conditions x number of voxels matrix) and $\mathbf{Z}$ the design matrix. Also influencing the data are effects of no interest $\mathbf{B}$ and noise:


\\[
\mathbf{Y} = \mathbf{ZU+XB}+\epsilon\\
\mathbf{u}_{p}  \sim N(\mathbf{0},\mathbf{G})\\
\epsilon_p \sim N(\mathbf{0},\mathbf{S}\sigma^{2})
\\]

There are a five assumptions in this generative model. First, the activity profiles ( $\mathbf{u}_p,$ columns of \textbf{U}) are considered to be a random variable drawn from a normal distribution. Representational models therefore do not specify the exact activity profiles of specific voxels, but simply the characteristics of the distribution from which they originate. Said differently, PCM is not interested in which voxel has which activity profiles - it ignores their spatial arrangement. This makes sense considering that activity patterns can vary widely across different participants \citep{RN3415} and do not directly impact what can be decoded from a region. For this, only the distribution of these activity profiles in this region is considered.

The second assumption is that the mean of the activity profiles (across voxels) is the same across conditions, and that it is modeled using the effects of no interests . Therefore, we most often model in  $\mathbf{X}$ the mean of each voxel across conditions. While one could also artificially remove the mean of each condition across voxels \citep{RN3565}, this approach would remove differences that, from the persepctive of decoding and representation, are highly meaningful \citep{RN3573}. 

The third assumption is that the activity profiles come from a multivariate Gaussian distribution. This is likely the most controversial assumption, but it is motivated by a few reasons: First, for fMRI data the multi-variate Gaussian is often a relatively appropriate description, especially if the mean of each voxel across conditions has been removed by the model. Secondly, the definition causes us to focus on the mean and covariance matrix, $\mathbf{G}$, as sufficient statistics, as these completely determine the Gaussian. Thus, even if the true distribution of the activity profiles is better described by a non-Gaussian distribution, the focus on the second moment is sensible as it characterizes the linear decodability of any feature of the stimuli.

Fourthly, the model assumes that different voxels are independent from each other. If we used raw data, this assumption would be clear violated, given the strong spatial correlation of noise processes in fMRI data. To reduce these dependencies we typically uses spatially pre-whitened data, which is divided by a estimate of the spatial covariance matrix \citep{RN3565, RN3573}. One complication here is that spatial pre-whitening usually does not remove spatial dependencies completely, given the estimation error in the spatial covariance matrix  \citep{RN3543}. 

Finally, we assume that the noise of each voxel is Gaussian with a temporal covariance that is known up to a constant term $\sigma^{2}$. Given the many additive influences of various noise sources on fMRI signals, Gaussianity of the noise is, by the central limit theorem, most likely a very reasonable assumption, which is commonly made in the fMRI literature. The original formulation of PCM used a model which assumed that the noise is also temporally independent and identically distributed (i.i.d.) across different trials, i.e. $\mathbf{S} = \mathbf{I}$ . However, as pointed out recently \citep{RN3638}, this assumption is often violated in non-random experimental designs with strong biasing consequences for estimates of the covariance matrix. If this is violated, we can either assume that we have a valid estimate of the true covariance structure of the noise ($S$), or we can model different parts of the noise structure (see section \ref{noise}).

When we fit a PCM model, we are not trying to estimate specific values of the the estimates of the true activity patterns $\mathbf{U}$. This is a difference to encoding approaches, in which we would estimate the values of $\mathbf{U}$ by estimating the feature weights $\mathbf{W}$. In PCM, we want to assess how likely the data is under any possible value of $\mathbf{U}$, as specified by the prior distribution. Thus we wish to calculate the marginal likelihood

\\[
p\left(\mathbf{Y}|\theta\right)=\int p\left(\mathbf{Y}|\mathbf{U},\theta\right) p\left(\mathbf{U}|\theta\right) d\mathbf{U}.
\\]

This is the likelihood that is maximized in PCM in respect to the model parameters $\theta$. For more details, see mathematical and algorithmic details.  

## Model types 
### Fixed models

In fixed models, the second moment matrix $\mathbf{G}$ is exactly predicted by the model. The simplest and most common example is the Null model, which states that $\mathbf{G} = \mathbf{0}$. This is equivalent to assuming that there is no difference between the activity patterns measured under any of the conditions. The Null-model is useful if we want to test whether there are any differences between experimental conditions.

Fixed models also occur when the representational structure can be predicted from some independent data. An example for this is shown in section \ref{example1}, where we predict the structure of finger representations directly from the correlational structure of finger movements in every-day life \citep{RN3415}. Importantly, fixed models only predict the the second moment matrix up to a proportional constant. The width of the distribution will vary with the overall signal-to-noise-level (assuming we use pre-whitened data). Thus, when evaluating fixed models we allow the predicted second moment matrix to be scaled by an arbitrary positive constant.

#### Example

An empirical example to for a fixed representational model comes from Ejaz et al (2015). Here the representational structure of 5 finger movements was compared to the representational structure predicted by the way the muscles are activated during finger movements (Muscle model), or by the covariance structure of natural movements of the 5 fingers. That is the predicted second moment matrix is derived from data completely independent of our imaging data.

Models are stored in structures, with the field `type` indicating the model type. To define a fixed model, we simple need to load the predicted second moment matrix and define a model structure as follows (see `pcm_recipe_finger`): 

```{matlab}
M.type = 'fixed’; % Type set to fixed
M.numGparams = 0; % Number of parameters
M.Gc = Model(1).G_cent; % This is the predicted second moment matrix from the behavioural data
M.name = 'muscle’;% This is the name of the
```

When evaluating the likelihood of a data set under the prediction, the pcm toolbox still needs to estimate the scaling factor and the noise variance, so even in the case of fixed models, an iterative maximization of the likelihood is required (see below).

### Component models

A more flexible model is to express the second moment matrix as a linear combination of different components. For example, the representational structure of activity patterns in the human object recognition system in inferior temporal cortex can be compared to the response of a convolutional neural network that is shown the same stimuli \citep{RN3544}. Each layer of the network predicts a specific structure of the second moment matrix and therefore constitutes a fixed model. However, the real representational structure seems to be best described by a mixture of multiple layers. In this case, the overall predicted second moment matrix is a linear sum of the weighted components matrices:


\\[
\label{componentModel}
\mathbf{G}= \sum_{h}{\exp(\theta_{h})\mathbf{G}_{h}}.
\\]


The weights for each component need to be positive - allowing negative weights would not guarantee that the overall second moment matrix would be positive definite. Therefore we use the exponential of the weighing parameter here, such that we can use unconstrained optimization to estimate the parameters.

For fast optimization of the likelihood, we require the derivate of the second moment matrix in respect to each of the parameters. Thus derivative can then be used to calculate the derivative of the log-likelihood in respect to the parameters (see section **4.3. Derivative of the log-likelihood**). In the case of linear component models this is easy to obtain.

\\[
\frac{\partial G}{\partial {\theta }_{h}}=\exp(\theta_{h}) {\bf{G}}_{h} 
\\]

#### Example

In the example `pcm_recipe_finger`, we have two fixed models, the Muscle and the natural statistics model. One question that arises in the paper is whether the real observed structure is better fit my a linear combination of the natural statistics and the muscle activity structure. So we can define a third model, which allows any arbitrary mixture between the two type.

```{matlab}
M.type = ‘component’;
M.numGparams = 2;
M.Gc(:,:,1) = Model(1).G_cent;
M.Gc(:,:,2) = Model(2).G_cent;
M.name = 'muscle + usage';
```

### Feature models
A representational model can be also formulated in terms of the features that are thought to be encoded in the voxels. Features are hypothetical tuning functions, i.e.\ models of what activation profiles of single neurons could look like. Examples of features would be Gabor elements for lower-level vision models \citep{RN3098}, elements with cosine tuning functions for different movement directions for models of motor areas \citep{RN2960}, and semantic features for association areas \citep{RN3566}. The actual activity profiles of each voxel are a weighted combination of the feature matrix $\mathbf{u}_p = \mathbf{M} \mathbf{w}_p$. The predicted second moment matrix of the activity profiles is then $\mathbf{G} = \mathbf{MM}^{T}$, assuming that all features are equally strongly and independently encoded, i.e.\ $E \left(\mathbf{w}_p\mathbf{w}_p^{T} \right)=\mathbf{I}$. A feature model can now be flexibly parametrized by expressing the feature matrix as a weighted sum of linear components.

\\[
\mathbf{M}= \sum_{h} \theta_h \mathbf{M}_{h}
\\]

Each parameter $\theta_h$ determines how strong the corresponding set of features is represented across the population of voxels. Note that this parameter is different from the actual feature weights $\mathbf{W}$.  Under this model, the second moment matrix becomes

\\[
	\mathbf{G}=\mathbf{UU}^{T}/P=\frac{1}{P}\sum_{h}\theta_{h}^{2}\mathbf{M}_{h}\mathbf{M}_{h}^{T}+\sum_{i}\sum_{j}\theta_{i}\theta_{j}\mathbf{M}_{i}\mathbf{M}_{j}^{T}.
\\]

From the last expression we can see that, if features that belong to different components are independent of each other, i.e. $\mathbf{M}_{i} \mathbf{M}_{j} = \mathbf{0}$, then a feature model is equivalent to a component model with $\mathbf{G}_h = \mathbf{M}_{h}\mathbf{M}_{h}^{T}$.  The only technical difference is that we use the square of the parameter $\theta_h$, rather than its exponential, to enforce non-negativity. Thus, component models assume that the different features underlying each component are encoded independently in the population of voxels - i.e.\ knowing something about the tuning to feature of component A does not tell you anything about the tuning to a feature of component B. If this cannot be assumed, then the representational model is better formulated as a feature model. 

By the product rule for matrix derivatives, we have
\\[
\frac{{\partial {\bf{G}}}}{{\partial {\theta _h}}} = {{\bf{M}}_h}{\bf{M}}{\left( \bf{\theta} \right)^T} + {\bf{M}}\left( \theta \right){\bf{M}}_h^T
\\]

#### Example 
Provide from correlation recipe 


### Nonlinear models

The most flexible way of defining a representational model is to express the second moment matrix as a non-linear (matrix valued) function of the parameters, $\mathbf{G}=F\left(\theta\right)$. While often a representational model can be expressed as a component or feature model, sometimes this is not possible. One example is a representational model in which the width of the tuning curve (or the width of the population receptive field) is a free parameter \citep{RN3558}. Such parameters would influence the features, and hence also the second-moment matrix in a non-linear way. Computationally, such non-linear models are not much more difficult to estimate than component or feature models, assuming that one can analytically derive the matrix derivatives $\partial \mathbf{G} / \partial \theta_{h}$. 

For this, the user needs to define a function that takes the parameters as an input and returns **G** the partial derivatives of **G** in respect to each of these parameters. The derivates are returned as a (KxKxH) tensor, where H is the number of parameters. 

```
\[G,dGdtheta\]=fcn(theta,data,…)
```
Note that this function is repeatedly called by the optimization routine and needs to execute fast. That is, any computation that does not depend on the current value of \theta should be performed outside the function and then passed to it.

#### Example 

SPENCER: PROVIDE FROM RECIPE

### Free models 

The most flexible representational model is the free model, in which the predicted second moment matrix is unconstrained. Thus, when we estimate this model, we would simply derive the maximum-likelihood estimate of the second-moment matrix. This can be useful for a number of reasons. First, we may want an estimate of the second moment matrix to derive the corrected correlation between different patterns, which is less influenced by noise than the simple correlation estimate \citep{RN3638, RN3033}. Furthermore, we may want to estimate the likelihood of the data under a free model to obtain a noise ceiling - i.e.\ an estimate of how well the best model should fit the data (see section \ref{noiseceilings}).

In estimating an unconstrained $\mathbf{G}$, it is important to ensure that the estimate will still be a positive definite matrix. For this purpose, we express the second moment as the square of an upper-triangular matrix, $\mathbf{G} = \mathbf{AA}^{T}$ \citep{RN3638, RN3033}. The parameters are then simply all the upper-triangular entries of $\mathbf{A}$.

#### Example 

QUICK EXAMPLE CODE 

## Noise models 

The noise is assumed to come from a multivariate normal distribution with covariance matrix $\mathbf{S}\sigma^{2}$. What is a reasonable noise structure to assume? First, the data can usually be assumed to be independent across imaging runs. If the data are regression estimates from a first-level model, and if the design of the experiment is balanced, then it is usually also permissible to make the assumption that the noise is independent within each imaging run $\mathbf{S}=\mathbf{I}$, \citep{RN3033}. Usually, however, the regression coefficients from a single imaging run show positive correlations with each other. This is due to the fact that the regression weights measure the activation during a condition as compared to a resting baseline, and the resting baseline is common to all conditions within the run \citep{RN3033}. To account for this, one can model the mean activation (across conditions) for each voxel with a separate fixed effect for each run. This effectively accounts for any uniform correlation.

Usually, assuming equal correlations of the activation estimates within a run is only a rough approximation to the real co-varince structure. A better estimate can be obtained by using an estimate derived from the design matrix and the estimated temporal autocorrelation of the raw signal. As pointed out recently \citep{RN3638}, the particular design can have substantial influence on the estimation of the second moment matrix. This is especially evident in cases where the design is such that the trial sequence is not random, but has an invariant structure (where trials of one condition are often to follow trials of another specific condition). The accuracy of our approximation hinges critically on the quality of our estimate of the temporal auto-covariance structure of the true noise. Note that it has been recently demonstrated that especially for high sampling rates, a simple autoregressive model of the noise is insufficient \citep{RN3550}. 

The last option is to estimate the covariance structure of the noise from the data itself. This can be achieved by introducing random effects into the generative model equation in section \ref{generativemodel}, which account for the covariance structure across the data. One example used here is to assume that the data are independent within each imaging run, but share an unknown covariance within each run, which is then estimated as a part of the covariance matrix \citep{RN3033}. While this approach is similar to just removing the run mean from the data as a fixed effect (see above) it is a good strategy if we actually want to model the difference of each activation pattern against the resting baseline. When treating the mean activation pattern in each run as a random effect, the algorithm finds a compromise between how much of the shared pattern in each run to ascribe to the random run-to-run fluctuations, and how much to ascribe to a stable mean activation.

ADD EXAMPLE AND CODE SNIPPETS FOR THE DIFFERENT OPTIONS

# Model fitting and comparison 
Details of the different optimization routines that maximize the likelihood can be found in the appendix. 

## Fitting to individual data sets

Function pcm_fitModelIndivid.m implements this step.

####  Example 
pcm_recipe_correlation.m (line 46)

## Fitting to individual data sets with cross-validation across partitions 

Function pcm_fitModelIndividCrossval.m implements these steps.

#### Example 

##  Fitting to group data sets

Function pcm_fitModelGroup.m implements this step.

#### Example

pcm_recipe_fixed.m (line 71)

pcm_recipe_nonlinear.m (line 128)

## Fitting to group data sets with cross-validation across participants 

Function pcm_fitModelGroupCrossval.m implements these steps.

#### Example 


# Data visualisation 

### Plotting Second Moment matrices 
 
Details of the different optimization routines that maximize the likelihood can be found in the appendix. 

### Multidimensional scaling 



### Plotting model evidence 
 


# Mathematical and Algorithmic details

## 4.1. Likelihood 

In this section, we derive the likelihood in the case that there are no fixed effects. In this case the distribution of the data would be

\\[
{\bf{y}} \sim N \left(0,{\bf{V}} \right)\\ {\bf{V}}=\bf{ZGZ^{T}+S}\sigma^{2}_{\epsilon}
\\]

To calculate the likelihood, let us consider at the level of the single voxel, namely , , and . Then, the likelihood for i-th single voxel is Then the likelihood over all voxels, assuming that the voxels are independent (e.g. effectively pre-whitened) is

\\[
p \left( {\bf{Y}}|{\bf{V}} \right)= \prod^{P}_{i=1} (2π)^{-\frac{N}{2}} |{\bf{V}}|^{-\frac{1}{2}} exp \left( -\frac{1}{2}{\bf{y}}_i^T {\bf{V}}^{-1} {\bf{y}}_i \right)
\\]

When we take the logarithm of this expression, the product over the individual Gaussian probabilities becomes a sum. and the exponential disappears:

\\[
L=\mathrm{ln}\left(p\left(\bf{Y}|V\right)\right) =\sum _{i=1}^{P}p\left(\bf{y}_{i}\right)\\ =\sum _{i=1}^{P}\left\[-\frac{N}{2}\mathrm{ln}\left(2\pi \right)-\frac{1}{2}\mathrm{ln}\left(|\bf{V}|\right)-\frac{1}{2}{\bf{y}}_{i}^{T}{\bf{V}}^{-1}{\bf{y}}_{i}\right\]\\ =-\frac{NP}{2}\mathrm{ln}\left(2\pi \right)-\frac{P}{2}\mathrm{ln}\left(|\bf{V}|\right)-\frac{1}{2}\sum _{i=1}^{P}{\bf{y}}_{i}^{T}{\bf{V}}^{-1}{\bf{y}}_{i}\\ =-\frac{NP}{2}\mathrm{ln}\left(2\pi \right)-\frac{P}{2}\mathrm{ln}\left(|\bf{V}|\right)-\frac{1}{2}trace\left({\bf{Y}}^{T}{\bf{V}}^{-1}\bf{Y}\right)
\\]

Using the trace trick, which allows $\mathrm{trace}\left(\bf{ABC}\right) = \mathrm{trace}\left(\bf{BCA}\right)$, we can obtain a form of the likelihood that does only depend on the second moment of the data, $\bf{YY}^{T}$ ,as a sufficient statistics:

\\[
L =-\frac{NP}{2}\mathrm{ln}\left(2\pi \right)-\frac{P}{2}\mathrm{ln}\left(|\bf{V}|\right)-\frac{1}{2}trace\left({\bf{Y}\bf{Y}}^{T}{\bf{V}}^{-1}\right)
\\]

## 4.2. Restricted likelihood 

In the presence of fixed effects (usually effects of no interest), we have the problem that the estimation of these fixed effects depends iterativly on the current estimate of \bf{V} and hence on the estimates of the second moment matrix and the noise covariance.

\\[
\bf{\hat{B}}} = {\left( {{{\bf{X}}^T}{{\bf{V}}^{ - 1}}{\bf{X}}} \right)^{ - 1}}{{\bf{X}}^T}{{\bf{V}}^{ - 1}}{\bf{Y}
\\]

Under the assumption of fixed effects, the distribution of the data is

\\[
{\bf{y_i}} \sim N \left(\bf{Xb_i},{\bf{V}} \right)
\\]

To compute the likelihood we need to remove these fixed effects from the data, using the residual forming matrix

\\[
{\bf{R}} = \bf{X}{\left( {{{\bf{X}}^T}{{\bf{V}}^{ - 1}}{\bf{X}}} \right)^{ - 1}}{{\bf{X}}^T}{{\bf{V}}^{ - 1}}\\ {\bf{r_i}} = \bf{Ry_i}
\\]

For the optimization of the random effects we therefore also need to take into account the uncertainty in the random effects estimates. Together this leads to a modified likelihood - the restricted likelihood that we which to optimize.

\\[
L_{ReML} =-\frac{NP}{2}\mathrm{ln}\left(2\pi \right)-\frac{P}{2}\mathrm{ln}\left(|\bf{V}|\right)-\frac{1}{2}trace\left({\bf{Y}\bf{Y}}^{T}\bf{R}^{T}{\bf{V}}^{-1}\bf{R}\right)-\frac{P}{2}\mathrm{ln}|\bf{X}^{T}\bf{V}^{-1}\bf{X}|
\\[

Note that the third term can be simplified by noting that

\\[
\bf{R}^{T}{\bf{V}}^{-1}\bf{R} = \bf{V}^{-1} - \bf{V}^{-1}\bf{X} (\bf{X}{\bf{V}}^{-1}\bf{X})^{-1}\bf{X}^{T}\bf{V}^{-1}=\bf{V}^{-1}\bf{R}=\bf{V}_{R}^{-1}
\\]

## 4.3. Derivatives of the log-likelihood 
Next, we find the derivatives of *L* with respect to each hyper parameter $ \theta_{i}$, which influence G. Also we need to estimate the hyper parameters that describe the noise, at least the noise parameter $\sigma_{\epsilon}^{2}$. To take these derivatives we need to use two general rules of taking derivatives of matrices (or determinants) of matrices:

\\[
\frac{{\partial \ln \left( {\bf{V}} \right)}}{{\partial {\theta _i}}} = tr\left( {{{\bf{V}}^{ - 1}}\frac{{\partial {\bf{V}}}}{{\partial {\theta _i}}}} \right)
\\]

\\[
\frac{{\partial {{\bf{V}}^{ - 1}}}}{{\partial {\theta _i}}} = {{\bf{V}}^{ - 1}}\left( {\frac{{\partial {\bf{V}}}}{{\partial {\theta _i}}}} \right){{\bf{V}}^{ - 1}}
\\]

Therefore the derivative of the log-likelihood in Eq. 5. in respect to each parameter is given by:

\\[
\frac{{\partial {L_{ML}}}}{{\partial {\theta _i}}} = - \frac{P}{2}tr\left( {{{\bf{V}}^{ - 1}}\frac{{\partial {\bf{V}}}}{{\partial {\theta _i}}}} \right) + \frac{1}{2}tr\left( {{{\bf{V}}^{ - 1}}\frac{{\partial {\bf{V}}}}{{\partial {\theta _i}}}{{\bf{V}}^{ - 1}}{\bf{Y}}{{\bf{Y}}^T}} \right)
\\]

## 4.4. Derivatives of the restricted log-likelihood 
## 4.5. Acceleration of matrix inversion 
To speed up the computation of the inverse of the variance-covariance matrix, we can exploit its special structure


\\[
 \begin{array}{c} {{\bf{V}}^{ - 1}} = {\left( {s{\bf{ZG}}{{\bf{Z}}^T} + {\bf{S}}\sigma _\varepsilon ^2} \right)^{ - 1}}\\ = {{\bf{S}}^{ - 1}}\sigma _\varepsilon ^{ - 2} - {{\bf{S}}^{ - 1}}{\bf{Z}}\sigma _\varepsilon ^{ - 2}{\left( {{s^{ - 1}}{G^{ - 1}} + {{\bf{Z}}^T}{{\bf{S}}^{ - 1}}{\bf{Z}}\sigma _\varepsilon ^{ - 2}} \right)^{ - 1}}{{\bf{Z}}^T}{{\bf{S}}^{ - 1}}\sigma _\varepsilon ^{ - 2}\\ = \left( {{{\bf{S}}^{ - 1}} - {{\bf{S}}^{ - 1}}{\bf{Z}}{{\left( {{s^{ - 1}}{G^{ - 1}}\sigma _\varepsilon ^2 + {{\bf{Z}}^T}{{\bf{S}}^{ - 1}}{\bf{Z}}} \right)}^{ - 1}}{{\bf{Z}}^T}{{\bf{S}}^{ - 1}}} \right)/\sigma _\varepsilon ^2 \end{array}
\\]
