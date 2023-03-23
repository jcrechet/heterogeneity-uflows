# Compute unemployment duration profiles

##

# load dataset
udur <- read.csv("Data/CPS_flows_uduration.csv")

# quarterly date
udur$date = as.yearqtr(udur$date)

# subset data for regression
my_data = subset( udur, udur$duration >= 0 & udur$duration <= 24 & udur$date > "1995 Q4" )

# size of vector duration parameters
K1 = length( unique( my_data$duration ) )


## estimate duration effects

# regression
results = lm( ue ~ 0 + factor( duration ) + factor( date ), data = my_data  )

# coefficients
coeffs = results$coefficients

# compute duration profile
ue = coeffs[1:K1] + mean( coeffs[(K1+1):length(coeffs)] )

# print
print("Unemployment duration profile computation: done")


###
###


# stack in dataset
uduration = unique( my_data$duration )
data_profiles = data.frame(uduration, ue)

# save in R directory
write.csv(data_profiles,"Data/CPS_udur_profiles.csv")

# Matlab export
file_path = paste(matlab_path, "/Data/CPS_udur_profiles.csv", sep = "")
write.csv(data_profiles,file_path)
