#####################################################
## beta regression with fixed precision parameter ####
#####################################################

source("mbrbetareg.R")
library(betareg)
data("FoodExpenditure", package = "betareg")
attach(FoodExpenditure)


## Fit the GLM using maximum likelihood using R package betareg 

fe_mle <- betareg(I(food/income) ~ income + persons, data = FoodExpenditure,type="ML")
coef(fe_mle)

## Mean bias-reduced fit usind the R package betareg

fe_br <- betareg(I(food/income) ~ income + persons, data = FoodExpenditure,type="BR")
coef(fe_br)

## Median bias-reduced fit 

par <- coef(fe_mle)                  ## initial values for the vector of parameters
y<- food/income                      ## response variable
X<- model.matrix(fe_mle)             ## design matrix involving covariates of the response mean
Z <- matrix(rep(1,nrow(X)),nrow(X),1)## design matrix invoving covariates of the response precision
eps <- 1e-06                         ## positive convergence tolerance 
maxit <- 100                         ## maximal number of IWLS iterations
link <-"logit"                       ##  link function for the mean 
phi_link <- "identity"               ## link function for the precision parameter
mlembr=mbrbetareg.fit(par=par,y=y,X=X,Z=Z,eps=eps,maxit=maxit,link=link,phi_link=phi_link)
mlembr$par

