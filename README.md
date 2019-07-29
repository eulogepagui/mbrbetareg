# Welcome to the ongoing work of the library mbrbetareg

# [library mbrbetareg]

* * *

## Set up
Here are the steps to follow in order to obtain the  
beta regression parameter estimates based on the median modified score function
proposed in Kenne Pagui et al. (2017):

**1-** Download the repository mbrbetareg on your own computer

**2-** Set the working directory on the terminal where there is the mbrbetareg repository

**3-** Compile the C code by typing the following command line on the terminal:

*R CMD SHLIB mod1.c*

**4-**  Open the file examples.R for the numerical example on beta regression

* * *

## Authors

The mbrbetareg implementation has been written by Euloge Clovis Kenne Pagui (kenne@stat.unipd.it) and Nicola Sartori (sartori@stat.unipd.it) From the University of Padova.

## References

Kenne Pagui, E. C., A. Salvan, and N. Sartori (2019). Efficient implementation of median bias reduction (in preparation).

Kenne Pagui, E. C., A. Salvan, and N. Sartori (2017). Median bias
reduction of maximum likelihood estimates. Biometrika, 104,923-938.
