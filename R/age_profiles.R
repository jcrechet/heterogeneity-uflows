# load dataset
CPS_flows_age <- read.csv("Data/CPS_flows_age.csv")

# date variable converted in R date format
CPS_flows_age$date = as.yearqtr(CPS_flows_age$date)

# generate new variables
CPS_flows_age$ue_s = NaN
CPS_flows_age$eu_s = NaN
CPS_flows_age$ei_s = NaN
CPS_flows_age$ui_s = NaN
CPS_flows_age$ie_s = NaN
CPS_flows_age$iu_s = NaN

#####
#####

## Time aggregation adjustment

# Date range
date_range = unique( subset( CPS_flows_age$date, CPS_flows_age$date >= "1995 Q4" ) )

# Age range
age_range = 18:64

# loop over time range
for (tt in date_range) { 
  for (aa in age_range) {
    
    # get vector
    tmp <- subset( CPS_flows_age, date == tt & age == aa )
    
    # transition matrices
    
    # row 1
    IE = tmp$ie
    IU = tmp$iu
    II = 1 - IE - IE
    
    # row 2
    UE = tmp$ue
    UI = tmp$ui
    UU = 1 - UE - UE
    
    # row 3
    EU = tmp$eu
    EI = tmp$ei 
    EE = 1 - EU - EI
    
    # matrix with unadjusted transition probabilities
    data = c( II, IU, IE, UI, UU, UE, EI, EU, EE )
    Pi_hat = matrix( unlist( data ), nrow = 3, ncol = 3, byrow = TRUE)
    
    # diagonalisation
    ev <- eigen( Pi_hat )
    llambda = ev$values  # eigenvalues
    pp = ev$vectors      # eigenvectors
    
    # check if eigenvalues are real
    if (is.complex( llambda ) == 'FALSE') {
      
      # check if eigenvalues are positive
      if ( all(llambda >= 0) ) {
        
        # compute log diag of eigenvalues
        log_llambda = log( llambda )
        Llambda = diag( log_llambda )
        
        # compute instantaneous rates
        instant_rates = pp %*% Llambda %*% inv( pp )
        
        # compute the adjusted monthly transition probability 
        Pi = 0*Pi_hat
        
        # construct the adjusted matrix
        for (i in 1:3) { for (j in 1:3) { if (i!=j) {
          Pi[i,j] = 1 - exp( - instant_rates[i,j] )
        } } }
        
        # stack into dataset
        # from I
        CPS_flows_age$iu_s[(CPS_flows_age$age == aa) & (CPS_flows_age$date == tt)] <-  Pi[1,2]
        CPS_flows_age$ie_s[(CPS_flows_age$age == aa) & (CPS_flows_age$date == tt)] <-  Pi[1,3]
        
        # from U
        CPS_flows_age$ui_s[(CPS_flows_age$age == aa) & (CPS_flows_age$date == tt)] <-  Pi[2,1]
        CPS_flows_age$ue_s[(CPS_flows_age$age == aa) & (CPS_flows_age$date == tt)] <-  Pi[2,3]
        
        # from E
        CPS_flows_age$ei_s[(CPS_flows_age$age == aa) & (CPS_flows_age$date == tt)] <-  Pi[3,1]
        CPS_flows_age$eu_s[(CPS_flows_age$age == aa) & (CPS_flows_age$date == tt)] <-  Pi[3,2]
        
      } # end of condition if eig positive
    } # end of condition if eig real
    
    if ( is.complex( llambda ) ) {
      
      print("eigenvalues are not real") 
      print("date, age: ") 
      print(tt) 
      print(aa)
      
    }
    
    if ( is.complex( llambda ) == FALSE ) {
      
      if ( any(llambda < 0) ) {
        
        print("eigenvalues are real but negative") 
        print("date, age: ") 
        print(tt) 
        print(aa)
        
      }
    }
    
  } # end loop age 
} # end loop date

# save adjusted flows
# save in R directory
write.csv(CPS_flows_age,"Data/CPS_flows_age.csv")


####
####

# Compute age profile (with time fixed-effect) 

# subset dataset for regression
my_data = subset( CPS_flows_age, CPS_flows_age$age >= 18 & CPS_flows_age$age <= 64 & CPS_flows_age$date > "1995 Q4" )

# size of vector age param
K1 = length( unique( my_data$age ) )

# UE rate
results = lm( ue_s ~ 0 + factor( age ) + factor( date ), data = my_data )

# coefficients
coeffs = results$coefficients

# age profile
ue = coeffs[1:K1] + mean( coeffs[(K1+1):length(coeffs)] )

# EU rate
results = lm( eu_s ~ 0 + factor( age ) + factor( date ), data = my_data )
coeffs = results$coefficients
eu = coeffs[1:K1] + mean( coeffs[(K1+1):length(coeffs)] )

# EE rate
results = lm( ee ~ 0 + factor( age ) + factor( date ), data = my_data )
coeffs = results$coefficients
ee = coeffs[1:K1] + mean( coeffs[(K1+1):length(coeffs)] )

# E pop ratio
results = lm( e ~ 0 + factor( age ) + factor( date ), data = my_data )
coeffs = results$coefficients
e = coeffs[1:K1] + mean( coeffs[(K1+1):length(coeffs)] )

# U pop ratio
results = lm( u ~ 0 + factor( age ) + factor( date ), data = my_data )
coeffs = results$coefficients
u = coeffs[1:K1] + mean( coeffs[(K1+1):length(coeffs)] )


# print
print("Age profiles computation: done")

####
####

# Age profile dataset
age = unique( my_data$age )
data_profiles = data.frame(age, e, u, ue, eu, ee)

# save in R directory
write.csv(data_profiles, "Data/CPS_age_profiles.csv")

# Matlab export
file_path = paste(matlab_path, "/Data/CPS_age_profiles.csv", sep = "")
write.csv(data_profiles,file_path)


