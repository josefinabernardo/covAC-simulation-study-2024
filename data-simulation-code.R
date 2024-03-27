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
             ct = sqrt(c(0, .01, .0025)), # Cultural Transmission - Parent genotype to child phenotype
             si = sqrt(c(0, .01, .0025)), # Sibling Interaction - Sibling 1 genotype to sibling 2 phenotype
             nloci = 100, # Number of diallelic loci
             npgsloci = c(2, 5, 10, 15) # Number of loci comprising the PGS
    ){
  
      #Create all possible parameter combinations
      param_combinations <- expand.grid(a = a, c = c, e = e, ct = ct, si = si)
      
      # Only retain distinct combinations where A + C + E = 1
      filtered_combinations <- param_combinations %>%
        filter(a^2 + c^2 + e^2 == 1) %>%
        distinct()
      
      # Number of settings we iterate through
      n_set = nrow(filtered_combinations)
      
      # R2 of the polygenic scores
      p_pgs = npgsloci/nloci
      
      # Initiate counters
      counter_within = 0  # counts sets within PGS setting
      counter_overall = 0 # counts sets overall
      
      # Determine numner of rows for data frames
      n_rows = n_set * length(p_pgs)
      
      # Pre-allocate data frames with the appropriate dimensions
      final_estimates <- data.frame(matrix(NA, nrow = n_rows, ncol = 19))
      final_power <- data.frame(matrix(NA, nrow = n_rows, ncol = 19))
      
      # Create 
      setkeep = matrix(NA, n_set, 30)   # to keep settings 
      reskeep = matrix(NA, n_set, 30)   # to keep results each data set
      mxkeep = matrix(NA, n_set, 30) # openmx results
      
      colnames(setkeep) = c(c('nmz','ndz','a','c','e','g','b','x','prs','A'), rep(NA, 20))
      
      # Print number of settings to the user
      print(paste('The factorial design has', n_set, 'setting(s).'))
}

dolan_simulation_function(a = sqrt(c(.5, .6)), c = sqrt(c(.2, .1)),
                          e = sqrt(c(.3, .3)), ct = sqrt(c(0, .01, .0025)))


# ISSUES

# Why do we round down the NCP?
# Why is there another type of sibling interaction model?

# Do some lintering and code-formatting.
# Why does the function say it's recursive when I start it without arguments?
# Have a closer look at the number of rows and columns for final_power & final_estimates.
# Make sure I can get out mx and gee regression results at the same time