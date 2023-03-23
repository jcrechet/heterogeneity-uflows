# Compute time-aggregation adjusted aggregate transitions

##

# load dataset
CPS_flows <- read.csv("Data/CPS_flows.csv")

# date variable converted in R format
CPS_flows$date = as.yearqtr(CPS_flows$date)

##

CPS_flows$ue_s = NaN
CPS_flows$eu_s = NaN
CPS_flows$ei_s = NaN

CPS_flows$ui_s = NaN
CPS_flows$ie_s = NaN
CPS_flows$iu_s = NaN

#####
#####


## Time aggregation adjustment

# Date range
date_range = subset( CPS_flows$date, CPS_flows$date >= "1995 Q4" )

# loop over time range
for (tt in date_range) {
  
    # get vector
    tmp <- subset( CPS_flows, date == tt )
    
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
    
    # compute log diag of eigenvalues
    llambda = log( llambda )
    Llambda = diag( llambda )
    
    # compute instantaneous rates
    instant_rates = pp %*% Llambda %*% inv( pp )
    
    # compute the adjusted monthly transition probability 
    Pi = 0*Pi_hat
    
    # construct the adjusted matrix
    for (i in 1:3) {
      for (j in 1:3) {
        
        if (i!=j) {
          Pi[i,j] = 1 - exp( - instant_rates[i,j] )
        }
      }
    }
    
# stack into dataset
# from I
CPS_flows$iu_s[(CPS_flows$date == tt)] <-  Pi[1,2]
CPS_flows$ie_s[(CPS_flows$date == tt)] <-  Pi[1,3]
    
# from U
CPS_flows$ui_s[(CPS_flows$date == tt)] <-  Pi[2,1]
CPS_flows$ue_s[(CPS_flows$date == tt)] <-  Pi[2,3]
    
# from E
CPS_flows$ei_s[(CPS_flows$date == tt)] <-  Pi[3,1]
CPS_flows$eu_s[(CPS_flows$date == tt)] <-  Pi[3,2]
    
}

# print
print("Aggregates: adjusted probabilities computation done")

# save
write.csv(CPS_flows,"Data/CPS_flows.csv")
