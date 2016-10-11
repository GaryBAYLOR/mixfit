# mixfit
an R package for mixture model fitting, model comparison and model selection for both raw data and grouped data
## Motivation
This project is motivated by my Preliminary Presentation Project (PPP), in which I studied a bootstrap method to select the number of components for Weibull mixture models fitted to grouped data. During this project I wrote R code for Weibull mixture model fitting, plotting and model selection by bootstrap. Even though the code I have is only for Weibull mixture model, there are a lot of similarities between Weibull mixture models and the ones using other distribution families. Being able to expand the existing code to include other common mixture models will be very meaningful.
## Goals
There are already some R packages out there on CRAN doing similar things, including `mixtools`, `mixdist`, `mclust` etc. All of them serve the purpose of working with mixture model, each has different emphases: `mixtools` fits and selects mixture models for raw data; `mixdist` fits mixture model with grouped and conditional data; `mclust` mainly focuses on doing classification by mixture models.Our goal is to create a package that takes advantages of existing packages (e.g. computation speed of `mclust`) and at the same time improve them in areas where we think is critical for users (e.g. the plotting system). In general, we want the package to be 
* versatile: being able to handle different mixture models, and for both raw data and grouped data
* computationally efficient: being able to run comparatively fast, especially in model selction by resampling methods
* good graphical output: use both base R plotting system and ggplot2 plotting system

## Methods
We use EM algorithm and extended EM algorithm, combined with other algorithms like Newton-Raphson algorithm.
