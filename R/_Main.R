# Project: heterogeneity in labor mobility
# Crechet (2022)

rm(list = ls())
shell("cls")

# directory path [fill in for replication]
mypath = ""

# working directory
setwd(paste(mypath,"/R",sep = ""))

# Matlab path
matlab_path = paste(mypath,"/Matlab",sep = "") 

# libraries
library(matlib)
library(zoo)

# 1. Compute time-aggregation adjusted aggregate flows 
source("aggregate_flows.R")

# 2. Age profiles of transition rates 
source("age_profiles.R")

# 3. Unemployment duration profiles
source("uduration_profiles.R")

# 4. Job tenure profiles 
source("tenure_profiles.R")
