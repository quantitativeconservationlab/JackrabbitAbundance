# Prep workspace
rm( list = ls() )
library(nimble)
library(dplyr)

# Load data
load( "./JckModel_Data.Rdata" )

# Reference:
#   i = sampling site (grid cell of 1 km)
#   j = replicated visits in a sampling occasion, i.e., survey nights.
#   T = sampling occasion (year)
#   t = number of nights within years
#   k = covariates for abundance
#   k2 = covariates for detection

# Model
Code <- nimbleCode({
  
  #### Priors
  
  # Priors for beta distribution to estimate mean Jackrabbit availability across locations and surveys. Mean availability then divides the estimated abundance (N[i, t]) to get density for a given cell.
  sh1~dunif(0, 15)
  sh2~dunif(0, 15)
  
  # Priors for covariates for detection (k2). 
  ## rd speed, temp, wind, moonphase
  for (k2 in 1:4) {
    pB[k2] ~ dnorm(0, 0.2) 
  }
  
  # Priors for t = year. 
  for (t in 1:T){
    
    # Prior for Intercept abundance
    lamB[1, t] ~ dnorm(0, 0.1)
    
    # Priors for covariates for abundance 
    ## k = 3 vegetation types, 2 distances (+ 1 intercept). 
    for (k in 2:4){ 
      lamB[k, t] ~ dnorm(0, 0.2) 
    }
    
    # Prior for intercept detection
    mean.pA[t] ~ dbeta(4,4)
    pA[t] <- log( mean.pA[t] / (1 - mean.pA[t]) )
    
    # Define formula to estimate Abundance for each year at each grid cell (Ni), with expected value N;i.
    for (i in 1:nsites){
      
       log(lambda[i, t]) <- lamB[1, t] + 
        lamB[2, t]*Annual[i,t] + 
        lamB[3, t]*Peren[i,t] +
        lamB[4, t]*Shrub[i,t] +
        log(SampArea[i] )
      N[i, t]~dpois(lambda[i, t])
    }
  }
  
  # Estimation of detection probability for each year
  for (i in 1:nsites){
    
  # Binomial detection probability (one observer)
  for (j in 1:t1){ 
    
    # For 2022, 8 surveys. 
    theta[i, j, 1] ~ dbeta(sh1, sh2)
    A[i, j, 1]~dbin(theta[i, j, 1], N[i, 1])
    # Probability of detection.
    logit(p[i, j, 1]) <- pA[1] + 
      pB[1]*Speed[i] + 
      pB[2]*Temp[j, 1] + 
      pB[3]*Moon[j, 1] + 
      pB[4]*Wind[j, 1]
    # Number of animals observed, depends on availability and detection probability. 
    # c because model is collapsed
    c[i, j] ~ dbin(p[i, j, 1], A[i, j, 1])
  }
  
    # Multinomial detection probability for 2023

    # For 2023 - 6 surveys (j)
    for (j in 1:t2){ 
      # Mean availability of Jackrabbits in each cell (i) and survey (j)
      theta[i, j, 2]~dbeta(sh1, sh2)
      # Individuals available for detection depends on mean availability and abundance.
      A[i, j, 2]~dbin(theta[i, j, 2], N[i, 2])
      # probability of detection. 
      logit(p[i, j, 2]) <- pA[2] + 
        pB[1]*Speed[i] + 
        pB[2]*Temp[j, 2] + 
        pB[3]*Moon[j, 2] + 
        pB[4]*Wind[j, 2]
      # Number of animals observed, depends on availability and detection probability
      n[i, j] ~ dbin(  p[i, j, 2], A[i, j, 2] )
      # We have two observers. 
      pi1[i, j, 1] <- p[i, j, 2]
      pi1[i, j, 2] <- (1-p[i, j, 2])*p[i, j, 2]
      y[i, j, 1:2] ~ dmulti(size=n[i, j], pi1[i, j, 1:2])
    }

    # Multinomial detection probability for 2024
    
        # For 2024 - 8 surveys (j)
      for (j in 1:t3){ 
        # Mean availability of Jackrabbits in each cell (i) and survey (j)
        theta[i, j, 3]~dbeta(sh1, sh2)
        # Individuals available for detection depends on mean availability and abundance.
        A[i, j, 3]~dbin(theta[i, j, 3], N[i, 3])
        # probability of detection. 
        logit(p[i, j, 3]) <- pA[3] + 
          pB[1]*Speed[i] + 
          pB[2]*Temp[j, 3] + 
          pB[3]*Moon[j, 3] + 
          pB[4]*Wind[j, 3]
        # Number of animals observed, depends on availability and detection probability
        m[i, j]~dbin(p[i, j, 3], A[i, j, 3])
        # We have two observers. 
        pi2[i, j, 1] <- p[i, j, 3] 
        pi2[i, j, 2] <- (1-p[i, j, 3])*p[i, j, 3]
        # Changed y by the corresponding detection matrix (z)
        z[i, j, 1:2] ~ dmulti(size = m[i, j], pi2[i, j, 1:2])
      }
    } 
})

# Inits --------
# dim: n_sites, n_surveys, n_years
nsites <- 200
#n_surveys <- 8
n_years <- 3
det_predictors <- 4
t1 <- 8
t2 <- 6
t3 <- 8
T <- 3

# Notes: Data is a list where c = all obs for 2022, and n = all obs for 2022
AIn <- array( NA, dim=c(nsites, t1, T) )
for (i in 1:nsites){
  # taking the count and adding 1 so they are more than are seen
  # number available for detection
  AIn[i,1:t1,1]<-Data$c[i, ]+1 
  AIn[i,1:t2,2]<-Data$n[i, ]+1
  AIn[i,1:t3,3]<-Data$m[i, ]+1
  }

# You cannot have an NA as initial values
AIn[is.na(AIn)]<-1
# just for the first year, 2022 that we don't have two teams
cIN<-matrix(NA, nsites, t1) 
cIN[is.na(Data$c)]<-0
# 2023
nIN<-matrix(NA, nsites, t2) 
nIN[is.na(Data$n)]<-0
# 2024
mIN<-matrix(NA, nsites, t3) 
mIN[is.na(Data$m)]<-0

# Doble observers for 2023
yIN <- array(NA, dim = c(nsites, t2, 2))
yIN[is.na(Data$y)]<-0

# Doble observers for 2024
zIN <- array(NA, dim = c(nsites, t3, 2))
zIN[is.na(Data$z)]<-0

N=apply(AIn, c(1, 3), max)*2

Inits=list(
  sh1=10, 
  sh2=10,
  # initial values for theta
  theta=array(.5, dim=c(nsites, t3, T)), 
  #intercept around 1 and coefficients as zero. Dim= predictors, and years. Initial values inside.
  lamB=matrix(c(1,0,0,0, 1,0,0,0, 1,0,0,0 ), 4, T),
  # Matrix for the detection predictor. Dim= predictors, and years.
  pA=rep(0, T),
  pB=rep(0, 4), 
  A=AIn, 
  # Initial value for Nt. It has to be more than could be available.
  N=N, 
  # Actual detection probability
  p=array(0.5, dim=c(nsites, t3, T)), 
  # For year 2022
  c=cIN,
  n=nIN,
  m=mIN,
  y=yIN,
  z=zIN
  )

# Model details --------
mod <- nimbleModel(
  code = Code, 
  name = "cc", 
  data = Data,
  constants = Consts,
  inits = Inits)

modconf <- configureMCMC(
  mod, 
  monitors = c("sh1","sh2","pB","pA","lamB"), 
  useConjugacy = FALSE) 

RMCMC <- buildMCMC(modconf)
CMCMC <- compileNimble(RMCMC, mod)

start <- Sys.time()

samps <- runMCMC(CMCMC$RMCMC,
                 niter = 150000,
                 nburnin = 75000,
                 thin = 25,
                 nchains = 3,
                 samplesAsCodaMCMC = TRUE)

end <- Sys.time()

end - start # 13 hours
  
### Write out posterior samples...
save(samps, file="./Postsamps.RData") 