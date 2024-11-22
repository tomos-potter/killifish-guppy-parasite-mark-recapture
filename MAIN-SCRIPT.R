# The MAIN SCRIPT to run to reproduce the analyses in the manuscript:

# “Novel parasite invasion destabilizes species coexistence” 

# by By Tomos Potter*, Ryan S. Mohammed, Joshua Goldberg, 
# David N Reznick, Joseph Travis, and Ronald D Bassar

#=========================================================================================================
# Load required packages
#=========================================================================================================

# data wrangling
library(plyr) 
library(tidyr)
library(data.table) 
library(reshape2)
# R interface for running mark-recapture analysis with Program MARK
library(RMark)      
# model fitting and estimating marginal means for analysis of parasite prevalance
library(glmmTMB)
library(emmeans)
# for building matrix projection models
library(popbio)
# for making figures
library(ggplot2)
library(ggimage)
library(ggthemes)
library(patchwork)

#=========================================================================================================
# Load the data
#=========================================================================================================

# killifish mark-recapture histories data file
H_killi <- read.csv('./data/killifish-mr-data.csv')
H_killi$capture_date <- as.Date(H_killi$capture_date)
H_killi$processing_date <- as.Date(H_killi$processing_date)
# guppy mark-recapture histories data file
H_guppy <- read.csv('./data/guppy-mr-data.csv')
H_guppy$capture_date <- as.Date(H_guppy$capture_date)
H_guppy$processing_date <- as.Date(H_guppy$processing_date)
# stream area data
stream_data <- fread('./data/stream-data.csv')
# dissection data
dissection_data <- fread('./data/dissection-data.csv')

#=========================================================================================================
# Run the mark-recapture analyses
#=========================================================================================================

# Run the POPAN mark-recapture analyses for killifish
# this script processes the data, runs model selection and averaging
# outputs model results, AIC tables, and parameter values
# N.B. takes about 5 minutes to run on my machine
source("./sub-scripts/POPAN-killifish.R")

# Run the POPAN mark-recapture analyses for guppies
# as for killifish
# N.B. takes about a minute to run
source("./sub-scripts/POPAN-guppies.R")

# Run the Multistrata (size-class structured) mark-recapture analysis for killifish
# N.B. this takes about 40 minutes to run on my machine
source("./sub-scripts/Multistrata-killifish.R")

#=========================================================================================================
# Analyses of parasite prevalence
#=========================================================================================================

# Run analysis of parasite prevalance, based on visual diagnoses from the mark-recpature data
# AND based on disection data from 5 other river systems
# N.B. this runs very quickly
source("./sub-scripts/Parasite-prevalance.R")

#=========================================================================================================
# Run MPMS
#=========================================================================================================

# Build the matrix projection models for killifish that estimate
# population growth rate (lambda) etc
# N.B. this takes ~2hrs on my machine
# because there are 10k iterations of each of 12 MPMs to build and analyse
source("./sub-scripts/MPMs.R")

#=========================================================================================================
# Figures
#=========================================================================================================

# Draw all the figures from the manuscript
source("./sub-scripts/Figures.R")

# END OF SCRIPT