# Tenure profiles

###
###

# load dataset
ten <- read.csv("Data/CPS_flows_tenure.csv")

# sample
my_data = subset( ten, ten$tenure >= 0 & ten$tenure <= 25 )

## regressions

# tenure coeff vector size
K1 = length( unique( my_data$tenure ) )


## estimations

# EU
results = lm( eu ~ 0 + factor(tenure) + factor(date), data = my_data )
coeffs = results$coefficients

# compute profile
eu = mean(coeffs[ (K1+1):length(coeffs) ] ) + coeffs[1:K1]

# EE 
results = lm( ee ~ 0 + factor(tenure) + factor(date), data = my_data )
coeffs = results$coefficients
ee = mean(coeffs[ (K1+1):length(coeffs) ] ) + coeffs[1:K1]

# print
print("Tenure profiles computation: done")

###
###


## stack in new dataset
tenure = unique( my_data$tenure )
data_profiles = data.frame( tenure, eu, ee)

# save in R directory
write.csv(data_profiles,"Data/CPS_tenure_profiles.csv")

# Matlab export
file_path = paste(matlab_path, "/Data/CPS_tenure_profiles.csv", sep = "")
write.csv(data_profiles,file_path)





