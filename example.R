#####################################################
## beta regression with fixed precision parameter ####
#####################################################

source("mbrbetareg.R")


##################################################
###       Food expenditure dataset             ###
### available on betareg R package             ###     
###  no covariates on the precision parameter  ###
##################################################

data("FoodExpenditure", package = "betareg")
attach(FoodExpenditure)


## maximum likelihood fit using R package betareg 
fe_mle <- betareg(I(food/income) ~ income + persons, data = FoodExpenditure,type="ML")
summary(fe_mle)

## Mean bias-reduced fit usind the R package betareg
fe_meanBR <- betareg(I(food/income) ~ income + persons, data = FoodExpenditure,type="BR")
summary(fe_meanBR)

## Median bias-reduced fit usind the R function  mbrbetareg 
fe_medianBR <- mbrbetareg(I(food/income) ~ income + persons, data = FoodExpenditure,type="medianBR")
summary(fe_medianBR)


#####################################################
###      Reading Skills  dataset               ######
###   available on betareg R package           ######
### covariates on the precision parameters     ######
#####################################################

data("ReadingSkills", package = "betareg")

## maximum likelihood fit using R package betareg 
rs_f <- accuracy ~ dyslexia * iq | dyslexia * iq
rs_ML<- betareg(rs_f, data = ReadingSkills, type = "ML")
summary(rs_ML)

## Mean bias-reduced fit usind the R package betareg
rs_meanBR<- betareg(rs_f, data = ReadingSkills, type = "BR")
summary(rs_meanBR)

## Median bias-reduced fit usind the R function  mbrbetareg
rs_medianBR<- mbrbetareg(rs_f, data = ReadingSkills, type = "medianBR")
summary(rs_medianBR)
