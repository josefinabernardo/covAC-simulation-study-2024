# Prepare environment
rm(list = ls())

# Load in packages
library(MASS)
library(OpenMx)
library(geepack)
library(tidyverse)

# Function to calculate power from non-centrality parameter
ncp_power_function <- function(ncp, df, alpha = .0001){
  critical_chi2 = qchisq(alpha, df, lower = FALSE)
  if(abs(ncp) < .0001){
    ncp = 0
    }
  power = pchisq(critical_chi2, df, ncp, lower = FALSE)
  return(power)
}

# Function to simulate data
dolan_simulation_function <- function(nrep = 500, # Number of repetitions
             alpha = .05, # Alpha for power
             cmethod = 'independence', # Gee Covariance Structure
             seed = NA, # Set a seed if desired
             standprs = FALSE, # Standardize the PRS
             nmz = 3999, # Sample size monozygotic twins
             ndz = 4001, # Sample size dizygotic twins
             a = sqrt(c(.5, .6)), # Additive genetic path coefficient
             c = sqrt(c(.2, .1)), # Shared environmental path coefficient
             e = sqrt(c(.3, .3)), # Unique environmental path coefficient
             ct = sqrt(c(0, .01, .0025)), # Parent genotype to child phenotype
             si = sqrt(c(0, .01, .0025)) # Twin 1 genotype to twin 2 phenotype
    ){
      nset = length(a) * length(c) * length(e) * length(ct) * length(si)
      nset
}

dolan_simulation_function(a = sqrt(c(.5, .6)), c = sqrt(c(.2, .1)),
                          e = sqrt(c(.3, .3)), ct = sqrt(c(0, .01, .0025)))


# ISSUES
# Why do we round down the NCP?
# Why is there another type of sibling interaction model?