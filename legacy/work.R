
################################################################################################
################################################################################################
## Tor's original (full) model
## exploring the data on my own
################################################################################################
################################################################################################

# DATA:
# - N[1]: Number of individuals that are only seen in the last primary session
#         or the session they are censored (these individuals must be sorted
#         first in the data)
# - N[2]:     Total number of individuals
# - n.prim:   Number of primary sessions
# - dt[k]:    Length of interval k
# - tod[k,j]: Category of trapping session k,j; two categories (e.g., time-of-day)
# - first[i]: First primary session of individual i
# - K[i]:     Last primary session for individual i (allows censoring)
# - J[i,k]:   Last secondary session for individual i in primary session k 
#             (allows censoring)
# - gr[i]:    Group (sex) of individual i
# - R:        Number of traps
# - X[r,]:    Location of trap r = 1..R (coordinate)
# - H[i,j,k]: Trap index for capture of individual i in secondary session j
#             within primary session k. 1 = not captured, other values are r+1
# - Ones[i,j,k]: An array of all ones used in the "Bernoulli-trick" (see 
#                OpenBUGS user manual) in the observation likelihood. This saves
#                computation time as it is not necessary to compute the complete 
#                capture probabilities for every individual*trap*session (it is
#                sufficient to compute the g-variable for every
#                ind*trap*primary)
# - xlow[i]: Lower bound in uniform prior for first x-location of individual i
# - xupp[i]: Upper bound in uniform prior for first x-location of individual i
# - ylow[i]: Lower bound in uniform prior for first y-location of individual i
# - yupp[i]: Upper bound in uniform prior for first y-location of individual i
#
# PARAMETERS:
# - kappa, sigma and lambda: Space use and recapture probability parameters
#                            (eq. 5 in paper)
# - beta: additive effects on log(lambda):
#           - beta[1]: effect of tod==2
#           - beta[2]: effect of gr==2
# - Phi[gr[i]]:   Survival for group gr[i] per time-unit 
# - phi[gr[i],k]: Survival probability for individual i belonging to group (sex) 
#                 gr[i] over the k'th interval given that the individual is
#                 alive at the beginning of the interval.
# - dmean[gr[i]]: Mean dispersal distance (exponential distribution) given
#                 dispersal for group gr[i]
# - Psi[gr[i]]:   Probability of not emigrating for group gr[i] per time-unit
#                 (only for zero-inflated model, commented out)
# - psi[gr[i],k]: Probability of dispersal (1-zero-inflation probability) for
#                 group gr[i] in interval k (only for zero-inflated model, 
#                 commented out)
# - d.offset[gr[i]]: Minimum dispersal distance given dispersal for group gr[i]
#                    (only for zero-inflated model, commented out)                 
# STATE VARIABLES:
# - z[i,k]: Alive state of individual i in primary session k (1=alive, 0=dead)
# - S[i,,k]: Centre of activity of individual i in primary session k (x- and
#            y-coordinate)
# - u[i,k]: Dispersal state for individual i in interval k (1=dispersal,
#           0=no dispersal)(only for zero-inflated model, commented out)

volesData <- dget('~/github/scr_ergon/original_files/dryad/ErgonAndGardner2013.rdat')

names(volesData)
## [1] "N"      "K"      "R"      "J"      "tod"    "first"  "X"      "H"     
## [9] "n.prim" "dt"     "gr"

volesData$N     ## c(8, 158)... 8 individuals only in last primary season. 158 tot. ind.
volesData$n.prim    ## 4 primary seasons
volesData$dt      ## 0.7666667 0.7666667 0.7000000  length of primary seasons 1,2,3
volesData$tod
volesData$first  ## length = 158 = N[2], first primary season each ind. was seen
sum(volesData$first == 4)   ## 5 individuals first seen in 4th (final) primary season
volesData$K  ## length  = 158
head(volesData$J)
dim(volesData$J)
volesData$gr  ## length  = 158
volesData$R   ## R = 192 traps
head(volesData$X)    ## X trap locations:
dim(volesData$X)   ##  dim = 192 x 2, 192 traps each with (x,y) coord
dim(volesData$H)  ## dim = 158 x 5 x 4 = individuals x 2nd session x primary session
dimnames(volesData$H)



################################################################################################
################################################################################################
## saving a standardized form of the FULL voles dataset
## including constants, data, and initial values,
## for testing JAGS and NIMBLE 
################################################################################################
################################################################################################

### FUNCTION FOR GENERATING INITIAL VALUES
## May need to be adjusted based on what are plausible starting values
Inits = function(H,X,K){
    n = dim(H)[1]
    n.prim = dim(H)[3]
    mean.x = apply(H, c(1,3), function(i) mean(X[i-1,1], na.rm=T))
    mean.y = apply(H, c(1,3), function(i) mean(X[i-1,2], na.rm=T))
                                        # For initial values of dispersal distance:
    for(i in (n.prim-1):1){
        mean.x[,i][is.na(mean.x[,i])] = mean.x[,i+1][is.na(mean.x[,i])]
        mean.y[,i][is.na(mean.y[,i])] = mean.y[,i+1][is.na(mean.y[,i])]
    }
    dx = mean.x[,2:n.prim] - mean.x[,1:(n.prim-1)]
    dy = mean.y[,2:n.prim] - mean.y[,1:(n.prim-1)]
    d.emp = sqrt(dx^2 + dy^2)
    ch = apply(H,c(1,3), function(i) any(i!=1))
    first.last = apply(ch, 1, function(i) range(which(i)))
    z = ch
    theta = atan2(dy,dx)
    theta[is.na(theta)] = runif(sum(is.na(theta)), -pi, pi)
    d = d.emp
    d[is.na(d)] = 0.001
    S = array(NA, c(n,2,n.prim)) # For initial values og FIRST location
    for(i in 1:n){
        S[i,1,first.last[1,i]] = mean.x[i,first.last[1,i]]
        S[i,2,first.last[1,i]] = mean.y[i,first.last[1,i]]
        z[i, first.last[1,i]:first.last[2,i]] = 1   # 1 when known to be alive, 0 otherwise
        if(first.last[1,i] != 1){
            theta[i,1:(first.last[1,i]-1)] = NA
            z[i,1:(first.last[1,i]-1)] = NA
            d[i,1:(first.last[1,i]-1)] = NA
        }
        if(K[i]!=n.prim){ # Adding NA's for censored individuals
            theta[i,K[i]:(n.prim-1)] = NA
            z[i, (K[i]+1):n.prim] = NA # after being removed
            d[i, K[i]:(n.prim-1)] = NA
        }
    }
    list(
        kappa = runif(2, 1.9, 2.1),
        sigma = runif(2, 4, 11), # Do not set too high if kappa is low
        PL = runif(1,0.7,0.9),
        beta = runif(2,0.5,2),
        dmean = c(10,15) + runif(2,-3,3),
        Phi = runif(2,0.5,0.8),
        theta = theta,
        d = d,
        z = z,
        S = S
    )
}

volesData <- dget('~/github/scr_ergon/original_files/dryad/ErgonAndGardner2013.rdat')

## Generating initial values for 1 chain
set.seed(0)
inits <- Inits(volesData$H, volesData$X, volesData$K)

## Priors for first centre of activity:
## Assuming that first centre of activity is always within +/- maxd
## from the mean capture location in both x and y direction.
trap.dist = 7
X <- volesData$X
H <- volesData$H
maxd = 2*trap.dist
mean.x = apply(H, c(1,3), function(i) mean(X[i-1,1], na.rm=T))
mean.y = apply(H, c(1,3), function(i) mean(X[i-1,2], na.rm=T))
first.mean.x = apply(mean.x,1, function(i) i[min(which(!is.na(i)))])
first.mean.y = apply(mean.y,1, function(i) i[min(which(!is.na(i)))])
xlow = first.mean.x - maxd
xupp = first.mean.x + maxd
ylow = first.mean.y - maxd
yupp = first.mean.y + maxd

constants <- list(
    N = volesData$N,
    K = volesData$K,
    R = volesData$R,
    J = volesData$J,
    tod = volesData$tod,
    first = volesData$first,
    X = volesData$X,
    H = volesData$H,
    n.prim = volesData$n.prim,
    dt = volesData$dt,
    gr = volesData$gr,
    xlow = xlow,
    xupp = xupp,
    ylow = ylow,
    yupp = yupp
)

data <- list(Ones = array(1, dim(volesData$H)))

save(constants, data, inits, file = '~/github/scr_ergon/analysis/volesData.RData')


for(i in 1:dim(inits$S)[1]) {print('=========================='); print(i); print(inits$S[i,,])}

for(i in 1:dim(inits$S)[1]) {
    print('==========================')
    print('==========================')
    print('==========================')
    print(i)
    print('H:')
    print(constants$H[i,,])
    print('z:')
    print(inits$z[i,])
}




################################################################################################
################################################################################################
## generating a smaller (reduced) version of the voles dataset,
## including constants, data, and initial values,
## for testing other NIMBLE models
################################################################################################
################################################################################################



volesData <- dget('~/github/scr_ergon/original_files/dryad/ErgonAndGardner2013.rdat')
names(volesData)
## [1] "N"      "K"      "R"      "J"      "tod"    "first"  "X"      "H"     
## [9] "n.prim" "dt"     "gr"

ind1 <- c(1, 6)
ind2 <- c(25, 51, 110)
ind <- c(ind1, ind2)

volesData_reduced <- list()
volesData_reduced$N <- c(length(ind1), length(ind))
volesData_reduced$K <- volesData$K[ind]
volesData_reduced$R <- volesData$R
volesData_reduced$J <- volesData$J[ind,]
volesData_reduced$tod <- volesData$tod
volesData_reduced$first <- volesData$first[ind]
volesData_reduced$X <- volesData$X
volesData_reduced$H <- volesData$H[ind,,]
volesData_reduced$n.prim <- volesData$n.prim
volesData_reduced$dt <- volesData$dt
volesData_reduced$gr <- volesData$gr[ind]

identical(names(volesData), names(volesData_reduced))

## FUNCTION FOR GENERATING INITIAL VALUES
## May need to be adjusted based on what are plausible starting values
Inits = function(H,X,K){
    n = dim(H)[1]
    n.prim = dim(H)[3]
    mean.x = apply(H, c(1,3), function(i) mean(X[i-1,1], na.rm=T))
    mean.y = apply(H, c(1,3), function(i) mean(X[i-1,2], na.rm=T))
                                        # For initial values of dispersal distance:
    for(i in (n.prim-1):1){
        mean.x[,i][is.na(mean.x[,i])] = mean.x[,i+1][is.na(mean.x[,i])]
        mean.y[,i][is.na(mean.y[,i])] = mean.y[,i+1][is.na(mean.y[,i])]
    }
    dx = mean.x[,2:n.prim] - mean.x[,1:(n.prim-1)]
    dy = mean.y[,2:n.prim] - mean.y[,1:(n.prim-1)]
    d.emp = sqrt(dx^2 + dy^2)
    ch = apply(H,c(1,3), function(i) any(i!=1))
    first.last = apply(ch, 1, function(i) range(which(i)))
    z = ch
    theta = atan2(dy,dx)
    theta[is.na(theta)] = runif(sum(is.na(theta)), -pi, pi)
    d = d.emp
    d[is.na(d)] = 0.001
    S = array(NA, c(n,2,n.prim)) # For initial values og FIRST location
    for(i in 1:n){
        S[i,1,first.last[1,i]] = mean.x[i,first.last[1,i]]
        S[i,2,first.last[1,i]] = mean.y[i,first.last[1,i]]
        z[i, first.last[1,i]:first.last[2,i]] = 1   # 1 when known to be alive, 0 otherwise
        if(first.last[1,i] != 1){
            theta[i,1:(first.last[1,i]-1)] = NA
            z[i,1:(first.last[1,i]-1)] = NA
            d[i,1:(first.last[1,i]-1)] = NA
        }
        if(K[i]!=n.prim){ # Adding NA's for censored individuals
            theta[i,K[i]:(n.prim-1)] = NA
            z[i, (K[i]+1):n.prim] = NA # after being removed
            d[i, K[i]:(n.prim-1)] = NA
        }
    }
    list(
        kappa = runif(2, 1.9, 2.1),
        sigma = runif(2, 4, 11), # Do not set too high if kappa is low
        PL = runif(1,0.7,0.9),
        beta = runif(2,0.5,2),
        dmean = c(10,15) + runif(2,-3,3),
        Phi = runif(2,0.5,0.8),
        theta = theta,
        d = d,
        z = z,
        S = S
    )
}

## Priors for first centre of activity:
## Assuming that first centre of activity is always within +/- maxd
## from the mean capture location in both x and y direction.
trap.dist = 7
X <- volesData_reduced$X
H <- volesData_reduced$H
maxd = 2*trap.dist
mean.x = apply(H, c(1,3), function(i) mean(X[i-1,1], na.rm=T))
mean.y = apply(H, c(1,3), function(i) mean(X[i-1,2], na.rm=T))
first.mean.x = apply(mean.x,1, function(i) i[min(which(!is.na(i)))])
first.mean.y = apply(mean.y,1, function(i) i[min(which(!is.na(i)))])
xlow = first.mean.x - maxd
xupp = first.mean.x + maxd
ylow = first.mean.y - maxd
yupp = first.mean.y + maxd

set.seed(0)
inits_reduced <- Inits(volesData_reduced$H, volesData_reduced$X, volesData_reduced$K)

names(inits_reduced)
## [1] "kappa" "sigma" "PL"    "beta"  "dmean" "Phi"   "theta" "d"     "z"    "S"    

dim(inits_reduced$S)
for(i in 1:dim(inits_reduced$S)[1]) {print('=========================='); print(i); print(inits_reduced$S[i,,])}

constants_reduced <- list(
    N = volesData_reduced$N,
    K = volesData_reduced$K,
    R = volesData_reduced$R,
    J = volesData_reduced$J,
    tod = volesData_reduced$tod,
    first = volesData_reduced$first,
    X = volesData_reduced$X,
    H = volesData_reduced$H,
    n.prim = volesData_reduced$n.prim,
    dt = volesData_reduced$dt,
    gr = volesData_reduced$gr,
    xlow = xlow,
    xupp = xupp,
    ylow = ylow,
    yupp = yupp
)

constants_reduced$N

data_reduced <- list(Ones = array(1, dim(volesData_reduced$H)))

save(constants_reduced, data_reduced, inits_reduced, file = '~/github/scr_ergon/analysis/volesData_reduced.RData')






################################################################################################
################################################################################################
## saving and comparing MCMC results (NIMBLE, JAGS, NIMBLE_SCR)
## using the smaller (reduced) version of the voles dataset
################################################################################################
################################################################################################


load('~/github/scr_ergon/analysis/samples_reduced.RData')


save(samples_jags_reduced, runtime_jags_reduced, ess_jags_reduced, eff_jags_reduced,
     samples_nimble_reduced, runtime_nimble_reduced, ess_nimble_reduced, eff_nimble_reduced,
     samples_scr1_reduced, runtime_scr1_reduced, ess_scr1_reduced, eff_scr1_reduced,
     samples_scr2_reduced, runtime_scr2_reduced, ess_scr2_reduced, eff_scr2_reduced,
     file = '~/github/scr_ergon/analysis/samples_reduced.RData')



load('~/github/scr_ergon/analysis/samples_reduced.RData')

## rename JAGS samples columns names to match NIMBLE (with spaces)
colnames(samples_jags_reduced) <- colnames(samples_nimble_reduced)
identical(colnames(samples_jags_reduced), colnames(samples_nimble_reduced))

## change 'lambda0' into 'log.lambda0' to match original model
## change needs to be made in: scr1, scr2
lambdaInd1 <- which(colnames(samples_scr1_reduced) == 'lambda0')
lambdaInd2 <- which(colnames(samples_scr2_reduced) == 'lambda0')
colnames(samples_scr1_reduced)[lambdaInd1] <- 'log.lambda0'
colnames(samples_scr2_reduced)[lambdaInd2] <- 'log.lambda0'
samples_scr1_reduced[, lambdaInd1] <- log(samples_scr1_reduced[, lambdaInd1])
samples_scr2_reduced[, lambdaInd2] <- log(samples_scr2_reduced[, lambdaInd2])

## create the combined samplesList
samplesList <- list(nimble = samples_nimble_reduced, jags = samples_jags_reduced, SCR1 = samples_scr1_reduced, SCR2 = samples_scr2_reduced)

chainsPlot(samplesList, nrow=2, legend.location = 'topleft')

## NIMBLE model improvement factors over JAGS
nimble_eff_ratio <- eff_nimble_reduced / eff_jags_reduced
nimble_eff_ratio
##     Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
##   5.102269    4.348108    4.002090    6.471835    4.074611    4.453131 
##   kappa[1]    kappa[2] log.lambda0   phi[1, 1]   phi[2, 1]   phi[1, 2] 
##   3.953121    3.619372    4.117930    4.581857    4.624109    4.581857 
##  phi[2, 2]   phi[1, 3]   phi[2, 3]    sigma[1]    sigma[2] 
##   4.624109    4.556765    4.642408    4.130959    2.889778 
range(nimble_eff_ratio)
##  2.889778 6.471835
summary(nimble_eff_ratio)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  2.890   4.075   4.453   4.398   4.624   6.472 


## dSCR1 improvement over simply porting original model into NIMBLE
scr1_eff_ratio <- eff_scr1_reduced / eff_nimble_reduced
scr1_eff_ratio
##    Phi[1]    Phi[2]   beta[1]   beta[2]  dmean[1]  dmean[2]  kappa[1]  kappa[2] 
##  1.715805  1.811289  1.715067  2.059461  1.762678  2.468944  1.492654  2.394276 
##   lambda0 phi[1, 1] phi[2, 1] phi[1, 2] phi[2, 2] phi[1, 3] phi[2, 3]  sigma[1] 
##  2.061677  1.898029  1.714301  1.898029  1.714301  1.893465  1.710554  1.662081 
##  sigma[2] 
##  2.431405 
range(scr1_eff_ratio)
##  1.492654 2.468944
summary(scr1_eff_ratio)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.493   1.714   1.811   1.906   2.059   2.469 


## dSCR2 improvement over simply porting original model into NIMBLE
scr2_eff_ratio <- eff_scr2_reduced / eff_nimble_reduced
scr2_eff_ratio
##    Phi[1]    Phi[2]   beta[1]   beta[2]  dmean[1]  dmean[2]  kappa[1]  kappa[2] 
##  2.433030  3.055343  1.717993  2.928221  2.231794  2.691159  1.275310  3.393752 
##   lambda0 phi[1, 1] phi[2, 1] phi[1, 2] phi[2, 2] phi[1, 3] phi[2, 3]  sigma[1] 
##  1.436885  2.600000  2.899075  2.600000  2.899075  2.537439  2.895992  1.383146 
##  sigma[2] 
##  3.167960 
range(scr2_eff_ratio)
##  1.275310 3.393752
summary(scr2_eff_ratio)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.275   2.232   2.600   2.479   2.899   3.394 








################################################################################################
################################################################################################
## WRITING Tor's original (full) model
## code out to file: 'RD-SCR.Exp.txt'
################################################################################################
################################################################################################



setwd('~/github/scr_ergon/analysis')

sink("RD-SCR.Exp.txt")
cat("
model{
## PRIORS AND CONSTRAINTS: ##
#  Space use and recapture probability parameters:
for(sex in 1:2){
  kappa[sex] ~ dunif(0,50)
  sigma[sex] ~ dunif(0.1,20)
}
for(sex in 1:2){
  for(TOD in 1:2){
    lambda[TOD, sex] <- lambda0 * pow(beta[1],(TOD-1)) * pow(beta[2],(sex-1))
  }
}
PL ~ dunif(0.01,0.99)
lambda0 <- -log(1-PL) 
beta[1] ~ dunif(0.1,10)
beta[2] ~ dunif(0.1,10)
# Survival parameters:
for(sex in 1:2){
  Phi[sex] ~ dunif(0,1)
  for(k in 1:(n.prim-1)){
    phi[sex,k] <- pow(Phi[sex], dt[k])
  }
}
# Dispersal parameters:
for(sex in 1:2){
  dmean[sex] ~ dunif(0,100)
  dlambda[sex] <- 1/dmean[sex]
# For zero-inflated model:  
#  d.offset[sex] ~ dunif(5,100)
#  Psi0[sex] ~ dunif(0,1)
#  for(k in 1:(n.prim-1)){
#    psi[sex,k] <- 1 - pow(Psi0[sex], dt[k])
#  }
}
## MODEL: ##
# Loop over individuals that are only seen in the last primary session or the
# session they are censored
for(i in 1:N[1]){
  z[i,first[i]] ~ dbern(1)
  S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
  S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
  g[i,first[i],1] <- 0
  for(r in 1:R){ # trap
      D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
      g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
  }
  G[i,first[i]] <- sum(g[i,first[i],]) # Total trap exposure
  for(j in 1:J[i,first[i]]){
      P[i,j,first[i]] <- 1 - exp(-lambda[tod[first[i],j],gr[i]]*G[i,first[i]]) # Probability of being captured
      PPII[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
      Ones[i,j,first[i]] ~ dbern(PPII[i,first[i],j])
  }
}
# Loop over all other individuals
for(i in (N[1]+1):N[2]){
  z[i,first[i]] ~ dbern(1)
  S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
  S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
  # First primary session:
  g[i,first[i],1] <- 0
  for(r in 1:R){ # trap
      D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
      g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
  }
  G[i,first[i]] <- sum(g[i,first[i],]) # Total trap exposure
  for(j in 1:J[i,first[i]]){
      P[i,j,first[i]] <- 1 - exp(-lambda[tod[first[i],j],gr[i]]*G[i,first[i]]) # Probability of being captured
      PPII[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
      Ones[i,j,first[i]] ~ dbern(PPII[i,first[i],j])
  }
  ## Later primary sessions
  for(k in (first[i]+1):K[i]){ # primary session
    theta[i,k-1] ~ dunif(-3.141593,3.141593) # Prior for dispersal direction 
    z[i,k] ~ dbern(Palive[i,k-1])
    Palive[i,k-1] <- z[i,k-1]*phi[gr[i],k-1] # Pr(alive in primary session k) gr[i] = sex
    d[i,k-1] ~ dexp(dlambda[gr[i]])
    # For zero-inflated model, replace line above with:
    #u[i,k-1] ~ dbern(psi[gr[i],k-1])
    #dd[i,k-1] ~ dexp(dlambda[gr[i]])
    #d[i,k-1] <- u[i,k-1]*(d.offset[gr[i]] + dd[i,k-1])
    S[i,1,k] <- S[i,1,k-1] + d[i,k-1]*cos(theta[i,k-1])
    S[i,2,k] <- S[i,2,k-1] + d[i,k-1]*sin(theta[i,k-1])
    g[i,k,1] <- 0
    for(r in 1:R){ # trap
      D[i,r,k] <- sqrt(pow(S[i,1,k]-X[r,1],2) + pow(S[i,2,k]-X[r,2],2))  # Squared distance to trap
      g[i,k,r+1] <- exp(-pow(D[i,r,k]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
    }
    G[i,k] <- sum(g[i,k,]) # Total trap exposure
    for(j in 1:J[i,k]){
      P[i,j,k] <- (1 - exp(-lambda[tod[k,j],gr[i]]*G[i,k]))*z[i,k] # Probability of being captured
      PPII[i,k,j] <- step(H[i,j,k]-2)*(g[i,k,H[i,j,k]]/(G[i,k] + 0.000000001))*P[i,j,k] + (1-step(H[i,j,k]-2))*(1-P[i,j,k])
      Ones[i,j,k] ~ dbern(PPII[i,k,j])
    }
  }
}
}
",fill = TRUE)
sink()




################################################################################################
################################################################################################
## Tor's original (full) model
## try running the model in JAGS
################################################################################################
################################################################################################

# Fitting the model in JAGS (may take several hours, and may need to be run for
# for much longer to get good convergence)
library(rjags)

load('~/github/scr_ergon/analysis/volesData.RData')
JAGS.data = c(constants, data)

niter <- 5000

t1 <- Sys.time()
set.seed(0)
jags = jags.model(file="~/github/scr_ergon/analysis/RD-SCR.Exp.txt", data=JAGS.data, inits=inits, n.chains=1, n.adapt=1000)
##   Observed stochastic nodes: 2297
##   Unobserved stochastic nodes: 1595
##   Total graph size: 1609000
t2 <- Sys.time()
as.numeric(difftime(t2, t1, units = 'secs'))  ## 1056 seconds = 17.6 minutes

## Parameters to monitor
##params = c("kappa", "sigma", "lambda0", "beta", "Psi0", "psi", "dmean", "phi", "Phi", "d.offset")
params = c("kappa", "sigma", "lambda0", "beta", "dmean", "phi", "Phi")

t1 <- Sys.time()
set.seed(0)
out <- coda.samples(jags, variable.names = params, n.iter = niter, thin = 1)
t2 <- Sys.time()

runtime <- as.numeric(difftime(t2, t1, units = 'secs'))
runtime  ## don't know the runtime... ?

samples <- as.matrix(out)
colnames(samples)
dim(samples)   ## [1] 5000   17
samplesSummary(samples)
##                 Mean    Median    St.Dev.  95%CI_low  95%CI_upp
##Phi[1]      0.7492555 0.7506094 0.03553235 0.67582945  0.8147092
##Phi[2]      0.9068083 0.9122053 0.04290200 0.81185509  0.9753736
##beta[1]     1.8001976 1.7968295 0.12065501 1.57399407  2.0476714
##beta[2]     0.1674828 0.1607082 0.04555049 0.10413850  0.2722192
##dmean[1]    3.0940841 3.0667923 0.37639386 2.43861277  3.9191676
##dmean[2]    5.9242398 5.8379447 1.05542248 4.10852384  8.2451105
##kappa[1]    0.7711348 0.7695902 0.05360774 0.67596312  0.8791385
##kappa[2]    1.6661953 1.6256970 0.26270614 1.27524480  2.3279814
##log.lambda0 0.5134471 0.5057048 0.22204791 0.06776161  0.9313259
##phi[1,1]    0.8013000 0.8025721 0.02918040 0.74052852  0.8546114
##phi[2,1]    0.9275565 0.9319751 0.03378292 0.85231512  0.9810650
##phi[1,2]    0.8013000 0.8025721 0.02918040 0.74052852  0.8546114
##phi[2,2]    0.9275565 0.9319751 0.03378292 0.85231512  0.9810650
##phi[1,3]    0.8168416 0.8180688 0.02717243 0.76012669  0.8663669
##phi[2,3]    0.9335923 0.9377019 0.03108295 0.86424119  0.9826972
##sigma[1]    1.4521012 1.4443410 0.29247181 0.96507891  2.0735402
##sigma[2]    8.7847620 8.6889354 1.28222175 6.48480517 11.5170777

samplesPlot(samples)

library(coda)
apply(samples, 2, effectiveSize)

setwd('~/github/scr_ergon/analysis')
load('results_reduced.RData')
class(out_jags)
names(out_jags)
class(out_jags$jags)
names(out_jags$jags)
#[1] "summary"    "timing"     "efficiency"
dimnames(out_jags$jags$summary)
out_jags$jags$timing
dim(out_jags$jags$summary)
dimnames(out_jags$jags$summary)[[3]]
dimnames(out_jags$jags$summary)[[2]]
#[1] "mean"       "median"     "sd"         "CI95_low"   "CI95_upp"   "n"         
#[7] "ess"        "efficiency"
orig_jags_colnames <- colnames(samples)
updated_colnames <- gsub(',', ', ', orig_jags_colnames)
colnames(samples) <- updated_colnames
identical(sort(colnames(samples)), sort(dimnames(out_jags$jags$summary)[[3]]))
out_jags_colnames <- dimnames(out_jags$jags$summary)[[3]]
out_jags_colnames
## I changed the column names of JAGS samples to be correct (but wrong order)
out_jags$jags$summary[1,'mean',] <- apply(samples, 2, mean)[out_jags_colnames]
out_jags$jags$summary[1,'median',] <- apply(samples, 2, median)[out_jags_colnames]
out_jags$jags$summary[1,'sd',] <- apply(samples, 2, sd)[out_jags_colnames]
out_jags$jags$summary[1,'CI95_low',] <- apply(samples, 2, function(x) quantile(x, 0.025))[out_jags_colnames]
out_jags$jags$summary[1,'CI95_upp',] <- apply(samples, 2, function(x) quantile(x, 0.975))[out_jags_colnames]
out_jags$jags$summary[1,'n',] <- nrow(samples)
out_jags$jags$summary[1,'ess',] <- apply(samples, 2, effectiveSize)[out_jags_colnames]
out_jags$jags$summary[1,'efficiency',] <- apply(samples, 2, effectiveSize)[out_jags_colnames] / runtime
##
out_jags$jags$timing['jags'] <- runtime
out_jags$jags$timing
##
out_jags$jags$efficiency$min <- min(out_jags$jags$summary[1,'efficiency',])
out_jags$jags$efficiency$mean <- mean(out_jags$jags$summary[1,'efficiency',])
names(out_jags$jags$efficiency$min) <- 'jags'
names(out_jags$jags$efficiency$mean) <- 'jags'
out_jags$jags$efficiency
out_jags$jags$efficiency$min
out_jags$jags$efficiency$mean
##
save(out_jags, file = 'results_jags.RData')
##
1




################################################################################################
################################################################################################
## Tor's original (full) model
## running Tor's model in nimble
################################################################################################
################################################################################################


library(nimble)

code <- nimbleCode({
    ## PRIORS AND CONSTRAINTS: ##
    ##  Space use and recapture probability parameters:
    for(sex in 1:2){
        kappa[sex] ~ dunif(0,50)
        sigma[sex] ~ dunif(0.1,20)
    }
    for(sex in 1:2){
        for(TOD in 1:2){
            lambda[TOD, sex] <- exp(log.lambda0) * pow(beta[1],(TOD-1)) * pow(beta[2],(sex-1))
        }
    }
    PL ~ dunif(0.01,0.99)
    log.lambda0 <- log(-log(1-PL)) 
    beta[1] ~ dunif(0.1,10)
    beta[2] ~ dunif(0.1,10)
    ##
    ## Survival parameters:
    for(sex in 1:2){
        Phi[sex] ~ dunif(0,1)
        for(k in 1:(n.prim-1)){
            phi[sex,k] <- pow(Phi[sex], dt[k])
        }
    }
    ##
    ## Dispersal parameters:
    for(sex in 1:2){
        dmean[sex] ~ dunif(0,100)
        dlambda[sex] <- 1/dmean[sex]
        ## For zero-inflated model:  
        ##  d.offset[sex] ~ dunif(5,100)
        ##  Psi0[sex] ~ dunif(0,1)
        ##  for(k in 1:(n.prim-1)){
        ##    psi[sex,k] <- 1 - pow(Psi0[sex], dt[k])
        ##  }
    }
    ##
    ## MODEL: ##
    ##
    ## Loop over individuals that are only seen in the last primary session or the
    ## session they are censored
    for(i in 1:N[1]){
        z[i,first[i]] ~ dbern(1)
        S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
        S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
        g[i,first[i],1] <- 0
        for(r in 1:R){ # trap
            D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
            g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
        }
        G[i,first[i]] <- sum(g[i,first[i],1:(R+1)]) # Total trap exposure
        for(j in 1:J[i,first[i]]){
            P[i,j,first[i]] <- 1 - exp(-lambda[tod[first[i],j],gr[i]]*G[i,first[i]]) # Probability of being captured
            PPII[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
            Ones[i,j,first[i]] ~ dbern(PPII[i,first[i],j])
        }
    }
    ##
    ## Loop over all other individuals
    for(i in (N[1]+1):N[2]){
        z[i,first[i]] ~ dbern(1)
        S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
        S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
        ## First primary session:
        g[i,first[i],1] <- 0
        for(r in 1:R){ # trap
            D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
            g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
        }
        G[i,first[i]] <- sum(g[i,first[i],1:(R+1)]) # Total trap exposure
        for(j in 1:J[i,first[i]]){
            P[i,j,first[i]] <- 1 - exp(-lambda[tod[first[i],j],gr[i]]*G[i,first[i]]) # Probability of being captured
            PPII[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
            Ones[i,j,first[i]] ~ dbern(PPII[i,first[i],j])
        }
        ## Later primary sessions
        for(k in (first[i]+1):K[i]){ # primary session
            theta[i,k-1] ~ dunif(-3.141593,3.141593) # Prior for dispersal direction 
            z[i,k] ~ dbern(Palive[i,k-1])
            Palive[i,k-1] <- z[i,k-1]*phi[gr[i],k-1] # Pr(alive in primary session k) gr[i] = sex
            d[i,k-1] ~ dexp(dlambda[gr[i]])
            ## For zero-inflated model, replace line above with:
            ##u[i,k-1] ~ dbern(psi[gr[i],k-1])
            ##dd[i,k-1] ~ dexp(dlambda[gr[i]])
            ##d[i,k-1] <- u[i,k-1]*(d.offset[gr[i]] + dd[i,k-1])
            S[i,1,k] <- S[i,1,k-1] + d[i,k-1]*cos(theta[i,k-1])
            S[i,2,k] <- S[i,2,k-1] + d[i,k-1]*sin(theta[i,k-1])
            g[i,k,1] <- 0
            for(r in 1:R){ # trap
                D[i,r,k] <- sqrt(pow(S[i,1,k]-X[r,1],2) + pow(S[i,2,k]-X[r,2],2))  # Squared distance to trap
                g[i,k,r+1] <- exp(-pow(D[i,r,k]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
            }
            G[i,k] <- sum(g[i,k,1:(R+1)]) # Total trap exposure
            for(j in 1:J[i,k]){
                P[i,j,k] <- (1 - exp(-lambda[tod[k,j],gr[i]]*G[i,k]))*z[i,k] # Probability of being captured
                PPII[i,k,j] <- step(H[i,j,k]-2)*(g[i,k,H[i,j,k]]/(G[i,k] + 0.000000001))*P[i,j,k] + (1-step(H[i,j,k]-2))*(1-P[i,j,k])
                Ones[i,j,k] ~ dbern(PPII[i,k,j])
            }
        }
    }
})


load('~/github/scr_ergon/analysis/volesData.RData')

system.time(Rmodel <- nimbleModel(code, constants, data, inits, calculate = FALSE))  ## 17 seconds
system.time(lp <- Rmodel$calculate())    ## 224.148 seconds!
lp
## [1] -8735.789

monitors = c('kappa', 'sigma', 'log.lambda0', 'beta', 'dmean', 'phi', 'Phi')

system.time(conf <- configureMCMC(Rmodel, monitors = monitors))  ## 5 seconds
conf$printSamplers()
conf$printMonitors()
system.time(Rmcmc <- buildMCMC(conf))  ## 22 seconds

system.time(Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE))   ## 20 seconds
system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE))   ## 60 seconds

system.time(lp <- Cmodel$calculate())    ## 28 seconds
lp
## [1] -8735.789

set.seed(0)
system.time(samples <- runMCMC(Cmcmc, 5000)) ## 5,000 samples: 446.168 seconds = 7.5 minutes

set.seed(0)
system.time(samples <- runMCMC(Cmcmc, 20000)) ## 20,000 samples: 1560 sec = 26 min

colnames(samples)
dim(samples)   ## [1] 5000   17
samplesSummary(samples)   ## below is from 20,000 samples:
##                 Mean    Median    St.Dev.  95%CI_low  95%CI_upp
##Phi[1]      0.7485944 0.7500457 0.03477904  0.6781672  0.8125263
##Phi[2]      0.9089585 0.9136052 0.04483934  0.8106180  0.9832178
##beta[1]     1.8234433 1.8033933 0.21675514  1.5789238  2.0979951
##beta[2]     0.2001723 0.1715974 0.15807143  0.1060900  0.3586842
##dmean[1]    3.0175298 3.0013089 0.36859083  2.3604041  3.7462296
##dmean[2]    6.0216409 5.9320390 1.19858882  3.9454782  8.6231615
##kappa[1]    0.8059716 0.7864722 0.11979516  0.6994998  1.0989429
##kappa[2]    1.6948574 1.6685302 0.27577474  1.1802471  2.3078590
##log.lambda0 0.3825414 0.4418947 0.39301714 -0.4711747  0.8070829
##phi[1, 1]   0.8007644 0.8021100 0.02856812  0.7424916  0.8528553
##phi[2, 1]   0.9292255 0.9330715 0.03530050  0.8513193  0.9871082
##phi[1, 2]   0.8007644 0.8021100 0.02856812  0.7424916  0.8528553
##phi[2, 2]   0.9292255 0.9330715 0.03530050  0.8513193  0.9871082
##phi[1, 3]   0.8163448 0.8176386 0.02660398  0.7619663  0.8647413
##phi[2, 3]   0.9351215 0.9387091 0.03247749  0.8633192  0.9882226
##sigma[1]    1.6652135 1.5343069 0.76434439  1.0838367  3.3320468
##sigma[2]    8.8778716 8.9262063 1.37246401  5.8070067 11.4433757

samplesPlot(samples)

library(coda)
ess <- apply(samples, 2, effectiveSize)
ess    ## below is from 20,000 samples
##     Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
## 2936.15385   927.10859   213.75203    45.86233   354.44970   243.04328 
##   kappa[1]    kappa[2] log.lambda0   phi[1, 1]   phi[2, 1]   phi[1, 2] 
##   15.51700    60.00479    31.72400  2932.49267   932.01296  2932.49267 
##  phi[2, 2]   phi[1, 3]   phi[2, 3]    sigma[1]    sigma[2] 
##  932.01296  2931.45332   933.41105    18.27393    48.45210 

eff <- ess / 1560
eff
##     Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
##1.882149906 0.594300376 0.137020533 0.029398928 0.227211348 0.155796973 
##   kappa[1]    kappa[2] log.lambda0   phi[1, 1]   phi[2, 1]   phi[1, 2] 
##0.009946795 0.038464611 0.020335899 1.879802994 0.597444208 1.879802994 
##  phi[2, 2]   phi[1, 3]   phi[2, 3]    sigma[1]    sigma[2] 
##0.597444208 1.879136744 0.598340416 0.011714060 0.031059035 





################################################################################################
################################################################################################
## Running the **reduced** data mdoel (volesData_reduced.RData)
## try running the model in JAGS
################################################################################################
################################################################################################

library(rjags)

load('~/github/scr_ergon/analysis/volesData_reduced.RData')
JAGS.data_reduced = c(constants_reduced, data_reduced)

niter <- 20000

t1 <- Sys.time()
set.seed(0)
jags_reduced <- jags.model(file='~/github/scr_ergon/analysis/RD-SCR.Exp.txt',data=JAGS.data_reduced,inits=inits_reduced,n.chains=1,n.adapt=1000)
##   Observed stochastic nodes: 48
##   Unobserved stochastic nodes: 44
##   Total graph size: 32004
t2 <- Sys.time()
as.numeric(difftime(t2, t1, units = 'secs'))  ## 12.34008 seconds


## Parameters to monitor
##params = c("kappa", "sigma", "log.lambda0", "beta", "Psi0", "psi", "dmean", "phi", "Phi", "d.offset")
params = c("kappa", "sigma", "log.lambda0", "beta", "dmean", "phi", "Phi")

t1 <- Sys.time()
set.seed(0)
out <- coda.samples(jags_reduced, variable.names = params, n.iter = niter, thin = 1)
t2 <- Sys.time()
runtime <- as.numeric(difftime(t2, t1, units = 'secs'))
runtime_jags_reduced <- runtime
runtime  ## 114.9386 seconds

samples <- as.matrix(out)
samples_jags_reduced <- samples
colnames(samples)
dim(samples)   ## 20000    17
samplesSummary(samples)  ## below is from 20,000 samples:
##                   Mean     Median    St.Dev.   95%CI_low  95%CI_upp
## Phi[1]       0.3790627  0.3262907  0.2738549  0.01294151  0.9429426
## Phi[2]       0.7582324  0.8023950  0.1903157  0.30974723  0.9921208
## beta[1]      2.6650054  2.3422824  1.3779395  0.92715205  6.3638410
## beta[2]      1.6295805  1.0402492  1.6984543  0.23217831  7.2090489
## dmean[1]    56.6223201 59.5241375 28.1013126  3.68652325 98.3306019
## dmean[2]    12.3439204  7.9288125 13.4066977  0.82322708 52.4052764
## kappa[1]    27.1446638 27.6942945 13.5837899  3.52145783 48.8098959
## kappa[2]    16.7866230 10.9399319 15.0317325  1.27968599 47.7416397
## log.lambda0 -1.4268371 -1.3939283  0.6555171 -2.79260121 -0.2263236
## phi[1,1]     0.4503035  0.4237375  0.2674330  0.03568797  0.9559578
## phi[2,1]     0.8035038  0.8446905  0.1620376  0.40716681  0.9939537
## phi[1,2]     0.4503035  0.4237375  0.2674330  0.03568797  0.9559578
## phi[2,2]     0.8035038  0.8446905  0.1620376  0.40716681  0.9939537
## phi[1,3]     0.4751250  0.4565867  0.2634918  0.04768585  0.9597093
## phi[2,3]     0.8174581  0.8571795  0.1526649  0.44025592  0.9944780
## sigma[1]     8.5507477  8.0132494  2.3411788  5.52307573 14.8531014
## sigma[2]    10.4974808 10.2807156  2.8245350  4.77859208 16.4594130

samplesPlot(samples)

library(coda)
ess <- apply(samples, 2, effectiveSize)
ess_jags_reduced <- ess
ess
##     Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
##  2708.9363   6080.9797   2860.5828    222.3537   1330.7588   1035.6647 
##   kappa[1]    kappa[2] log.lambda0    phi[1,1]    phi[2,1]    phi[1,2] 
##  5878.1254    684.6472    805.6314   2635.5746   6041.5196   2635.5746 
##   phi[2,2]    phi[1,3]    phi[2,3]    sigma[1]    sigma[2] 
##  6041.5196   2662.0019   6032.0622    709.1884    265.8138 

eff <- ess / runtime
eff_jags_reduced <- eff
eff
##     Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
##  12.925774   29.015585   13.649360    1.060968    6.349757    4.941706 
##   kappa[1]    kappa[2] log.lambda0    phi[1,1]    phi[2,1]    phi[1,2] 
##  28.047659    3.266815    3.844095   12.575726   28.827299   12.575726 
##   phi[2,2]    phi[1,3]    phi[2,3]    sigma[1]    sigma[2] 
##  28.827299   12.701825   28.782173    3.383914    1.268339 


################################################################################################
################################################################################################
## Running the **reduced** data mdoel (volesData_reduced.RData)
## try running the model in NIMBLE,
## using Tor's originally written model code
################################################################################################
################################################################################################

library(nimble)

code <- nimbleCode({
    for(sex in 1:2){
        kappa[sex] ~ dunif(0,50)
        sigma[sex] ~ dunif(0.1,20)
    }
    for(sex in 1:2){
        for(TOD in 1:2){
            lambda[TOD, sex] <- exp(log.lambda0) * pow(beta[1],(TOD-1)) * pow(beta[2],(sex-1))
        }
    }
    PL ~ dunif(0.01,0.99)
    log.lambda0 <- log(-log(1-PL)) 
    beta[1] ~ dunif(0.1,10)
    beta[2] ~ dunif(0.1,10)
    for(sex in 1:2){
        Phi[sex] ~ dunif(0,1)
        for(k in 1:(n.prim-1)){
            phi[sex,k] <- pow(Phi[sex], dt[k])
        }
    }
    for(sex in 1:2){
        dmean[sex] ~ dunif(0,100)
        dlambda[sex] <- 1/dmean[sex]
    }
    for(i in 1:N[1]){
        z[i,first[i]] ~ dbern(1)
        S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
        S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
        g[i,first[i],1] <- 0
        for(r in 1:R){ # trap
            D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
            g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
        }
        G[i,first[i]] <- sum(g[i,first[i],1:(R+1)]) # Total trap exposure
        for(j in 1:J[i,first[i]]){
            P[i,j,first[i]] <- 1 - exp(-lambda[tod[first[i],j],gr[i]]*G[i,first[i]]) # Probability of being captured
            PPII[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
            Ones[i,j,first[i]] ~ dbern(PPII[i,first[i],j])
        }
    }
    for(i in (N[1]+1):N[2]){
        z[i,first[i]] ~ dbern(1)
        S[i,1,first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
        S[i,2,first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
        ## First primary session:
        g[i,first[i],1] <- 0
        for(r in 1:R){ # trap
            D[i,r,first[i]] <- sqrt(pow(S[i,1,first[i]]-X[r,1],2) + pow(S[i,2,first[i]]-X[r,2],2))
            g[i,first[i],r+1] <- exp(-pow(D[i,r,first[i]]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
        }
        G[i,first[i]] <- sum(g[i,first[i],1:(R+1)]) # Total trap exposure
        for(j in 1:J[i,first[i]]){
            P[i,j,first[i]] <- 1 - exp(-lambda[tod[first[i],j],gr[i]]*G[i,first[i]]) # Probability of being captured
            PPII[i,first[i],j] <- step(H[i,j,first[i]]-2)*(g[i,first[i],H[i,j,first[i]]]/(G[i,first[i]]+ 0.000000001))*P[i,j,first[i]] + (1-step(H[i,j,first[i]]-2))*(1-P[i,j,first[i]])
            Ones[i,j,first[i]] ~ dbern(PPII[i,first[i],j])
        }
        for(k in (first[i]+1):K[i]){ # primary session
            theta[i,k-1] ~ dunif(-3.141593,3.141593) # Prior for dispersal direction 
            z[i,k] ~ dbern(Palive[i,k-1])
            Palive[i,k-1] <- z[i,k-1]*phi[gr[i],k-1] # Pr(alive in primary session k) gr[i] = sex
            d[i,k-1] ~ dexp(dlambda[gr[i]])
            S[i,1,k] <- S[i,1,k-1] + d[i,k-1]*cos(theta[i,k-1])
            S[i,2,k] <- S[i,2,k-1] + d[i,k-1]*sin(theta[i,k-1])
            g[i,k,1] <- 0
            for(r in 1:R){ # trap
                D[i,r,k] <- sqrt(pow(S[i,1,k]-X[r,1],2) + pow(S[i,2,k]-X[r,2],2))  # Squared distance to trap
                g[i,k,r+1] <- exp(-pow(D[i,r,k]/sigma[gr[i]], kappa[gr[i]])) # Trap exposure
            }
            G[i,k] <- sum(g[i,k,1:(R+1)]) # Total trap exposure
            for(j in 1:J[i,k]){
                P[i,j,k] <- (1 - exp(-lambda[tod[k,j],gr[i]]*G[i,k]))*z[i,k] # Probability of being captured
                PPII[i,k,j] <- step(H[i,j,k]-2)*(g[i,k,H[i,j,k]]/(G[i,k] + 0.000000001))*P[i,j,k] + (1-step(H[i,j,k]-2))*(1-P[i,j,k])
                Ones[i,j,k] ~ dbern(PPII[i,k,j])
            }
        }
    }
})

load('~/github/scr_ergon/analysis/volesData_reduced.RData')

system.time(Rmodel <- nimbleModel(code, constants_reduced, data_reduced, inits_reduced, calculate = FALSE))  ## 3 seconds
system.time(lp <- Rmodel$calculate())    ## 1 second
lp  ##  -167.1596

monitors = c('kappa', 'sigma', 'log.lambda0', 'beta', 'dmean', 'phi', 'Phi')

system.time(conf <- configureMCMC(Rmodel, monitors = monitors))  ## 0.1 seconds
conf$printSamplers()
length(conf$getSamplers())   ## 44 un-observed stochastic nodes
conf$printMonitors()
system.time(Rmcmc <- buildMCMC(conf))  ## 0.5 seconds

system.time(Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE))   ## 21 seconds
system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE))   ## 6 seconds

system.time(lp <- Cmodel$calculate())    ## 0.5 seconds
lp    ##  -167.1596

niter <- 20000

set.seed(0)
t1 <- Sys.time()
samples <- runMCMC(Cmcmc, niter)
t2 <- Sys.time()
runtime <- as.numeric(difftime(t2, t1, units = 'secs'))
runtime_nimble_reduced <- runtime
runtime  ## 21 seconds

samples_nimble_reduced <- samples

colnames(samples)
dim(samples)   ##  20000    17
samplesSummary(samples)   ## below is from 20,000 samples:
##                   Mean     Median    St.Dev.   95%CI_low  95%CI_upp
## Phi[1]       0.3711870  0.3194464  0.2650928  0.01501541  0.9367876
## Phi[2]       0.7603304  0.8004303  0.1858904  0.31984946  0.9919560
## beta[1]      2.6935573  2.4286067  1.3440834  0.91917071  6.1636370
## beta[2]      1.4887397  0.9947450  1.4887032  0.22711649  6.0586816
## dmean[1]    56.9697239 59.1063781 27.3014817  6.12401795 97.5857572
## dmean[2]    13.2056796  8.2410685 14.8693100  0.92430397 61.4057467
## kappa[1]    26.9368238 27.3050553 13.6810872  2.99588117 48.7740175
## kappa[2]    15.8379543  9.9787899 14.4776512  1.26822402 47.8536396
## log.lambda0 -1.4027894 -1.3788776  0.6165128 -2.69328795 -0.2670376
## phi[1, 1]    0.4439772  0.4169063  0.2593910  0.03999556  0.9511701
## phi[2, 1]    0.8054820  0.8431044  0.1581477  0.41730954  0.9938271
## phi[1, 2]    0.4439772  0.4169063  0.2593910  0.03999556  0.9511701
## phi[2, 2]    0.8054820  0.8431044  0.1581477  0.41730954  0.9938271
## phi[1, 3]    0.4693403  0.4498612  0.2556385  0.05291466  0.9553198
## phi[2, 3]    0.8193718  0.8557097  0.1489742  0.45025852  0.9943624
## sigma[1]     8.5991051  8.0719104  2.2918016  5.59300864 14.7615577
## sigma[2]    10.8733746 10.5966919  2.9032874  4.93685901 17.2023690

samplesPlot(samples)

library(coda)
ess <- apply(samples, 2, effectiveSize)
ess_nimble_reduced <- ess
ess
##      Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
##  1438.93330  2752.65910  1191.84547   149.81328   564.50007   480.13475 
##    kappa[1]    kappa[2] log.lambda0   phi[1, 1]   phi[2, 1]   phi[1, 2] 
##  2419.12076   257.97557   345.37742  1257.17406  2908.39114  1257.17406 
##   phi[2, 2]   phi[1, 3]   phi[2, 3]    sigma[1]    sigma[2] 
## 2908.39114  1262.82611  2915.32983   304.99383    79.96877

eff <- ess / runtime
eff_nimble_reduced <- eff
eff
##     Phi[1]      Phi[2]     beta[1]     beta[2]    dmean[1]    dmean[2] 
##  65.950770  126.162893   54.625970    6.866407   25.872787   22.006063 
##   kappa[1]    kappa[2] log.lambda0   phi[1, 1]   phi[2, 1]   phi[1, 2] 
## 110.875798   11.823819   15.829717   57.620181  133.300575   57.620181 
##  phi[2, 2]   phi[1, 3]   phi[2, 3]    sigma[1]    sigma[2] 
## 133.300575   57.879232  133.618596   13.978812    3.665217 



################################################################################################
################################################################################################
## NEW NIMBLE MODEL #1: dSCR1, or scr1
## using a custom distribution to calculate likelihood (EXACTLY)
## using REDUCED voles dataset
################################################################################################
################################################################################################


library(nimble)

code_dSCR1 <- nimbleCode({
    ## space use and recapture probability parameters
    PL ~ dunif(0.01, 0.99)
    lambda0 <- -log(1-PL)
    for(sex in 1:2) {
        kappa[sex] ~ dunif(0,   50)
        sigma[sex] ~ dunif(0.1, 20)
        beta[sex]  ~ dunif(0.1, 10)    # misnomer: beta[1] is coeff of tod, beta[2] is coeff of sex
        for(TOD in 1:2) {
            lambda[TOD, sex] <- lambda0 * beta[1]^(TOD-1) * beta[2]^(sex-1)
        }
        ## survival parameters
        Phi[sex] ~ dunif(0, 1)
        for(k in 1:(nPrimary-1)) {
            phi[sex, k] <- Phi[sex]^dt[k]
        }
        ## dispersal parameters
        dmean[sex] ~ dunif(0, 100)
        dlambda[sex] <- 1/dmean[sex]
    }
    for(i in 1:nInd) {
        S[i, 1, first[i]] ~ dunif(xlow[i], xupp[i])  # initial center of activity (x)
        S[i, 2, first[i]] ~ dunif(ylow[i], yupp[i])  # initial center of activity (y)
        z[i, first[i]] <- 1
        for(k in first[i]:last[i]) {
            D[i, k, 1:R] <- sqrt((S[i, 1, k] - X[1:R, 1])^2 + (S[i, 2, k] - X[1:R, 2])^2)
            g[i, k, 1:R] <- exp(-(D[i, k, 1:R]/sigma[gr[i]])^kappa[gr[i]])  # trap exposure
            G[i, k] <- sum(g[i, k, 1:R])                                    # total trap exposure
        }
        for(k in first[i]:(last[i]-1)) {
            theta[i, k] ~ dunif(-3.141593, 3.141593)   # dispersal direction
            d[i, k] ~ dexp(dlambda[gr[i]])
            S[i, 1, k+1] <- S[i, 1, k] + d[i, k] * cos(theta[i, k])
            S[i, 2, k+1] <- S[i, 2, k] + d[i, k] * sin(theta[i, k])
            Palive[i, k] <- z[i, k] * phi[gr[i], k]
            z[i, k+1] ~ dbern(Palive[i, k])
        }
        ## likelihood
        H[i, 1:nSecondary, 1:nPrimary] ~ dSCR1(
            first = first[i], last = last[i], J = J[i,1:nPrimary],
            lambda = lambda[1:2,gr[i]], tod = tod[1:nPrimary,1:nSecondary],
            g = g[i,1:nPrimary,1:R], G = G[i,1:nPrimary], z = z[i,1:nPrimary])
    }
})

## define custom distribution
dSCR1 <- nimbleFunction(
    run = function(x = double(2),
        first = double(), last = double(), J = double(1),
        lambda = double(1), tod = double(2),
        g = double(2), G = double(1), z = double(1),
        log = double()) {
        lp <- 0
        for(k in first:last) {    # primary session
            for(j in 1:J[k]) {    # secondary session
                logPnoCapture <- -lambda[tod[k,j]] * G[k] * z[k]
                if(x[j,k] == 1) {    # not captured
                    lp <- lp + logPnoCapture
                } else {             # captured
                    lp <- lp + log(1-exp(logPnoCapture)) + log(g[k, x[j,k]-1]) - log(G[k])
                }
            }
        }
        returnType(double())
        if(log) return(lp) else return(exp(lp))
    }
)

rSCR1 <- nimbleFunction(
    run = function(n = integer(),
        first = double(), last = double(), J = double(1),
        lambda = double(1), tod = double(2),
        g = double(2), G = double(1), z = double(1)) {
        x <- array(1, c(dim(tod)[2], dim(tod)[1]))
        returnType(double(2))
        return(x)
    }
)

registerDistributions(list(
    dSCR1 = list(
        BUGSdist = 'dSCR1(first, last, J, lambda, tod, g, G, z)',
        types = c('value = double(2)', 'first = double()', 'last = double()', 'J = double(1)', 'lambda = double(1)', 'tod = double(2)', 'g = double(2)', 'G = double(1)', 'z = double(1)'),
        discrete = TRUE,
        mixedSizes = TRUE
    )
))


load('~/github/scr_ergon/analysis/volesData_reduced.RData')
constants_dSCR1 <- constants_reduced[c('R', 'J', 'tod', 'first', 'X',
                                       'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR1$last       <- constants_reduced$K
constants_dSCR1$nInd       <- dim(constants_reduced$H)[1]
constants_dSCR1$nPrimary   <- dim(constants_reduced$H)[3]
constants_dSCR1$nSecondary <- dim(constants_reduced$H)[2]
data_dSCR1 <- list(H=constants_reduced$H)
inits_dSCR1 <- inits_reduced


system.time(Rmodel <- nimbleModel(code_dSCR1, constants_dSCR1, data_dSCR1, inits_dSCR1, calculate = FALSE))  ## 1.9 seconds

system.time(lp <- Rmodel$calculate())    ## .05 second
lp  ##  -167.1596

## NOTE: 'log.lambda0' changed to 'lambda0'
monitors = c('kappa', 'sigma', 'lambda0', 'beta', 'dmean', 'phi', 'Phi')

system.time(conf <- configureMCMC(Rmodel, monitors = monitors))  ## 0.1 seconds
conf$printSamplers()
system.time(Rmcmc <- buildMCMC(conf))  ## 0.25 seconds

system.time(Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE))   ## 16 seconds
system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE))   ## 4.5 seconds

system.time(lp <- Cmodel$calculate())    ## 0.02 seconds
lp    ##  -167.1596

niter <- 20000

set.seed(0)
t1 <- Sys.time()
samples <- runMCMC(Cmcmc, niter)
t2 <- Sys.time()
runtime <- as.numeric(difftime(t2, t1, units = 'secs'))
runtime_scr1_reduced <- runtime
runtime  ## 12 seconds

samples_scr1_reduced <- samples

colnames(samples)
dim(samples)   ##  20000    17
samplesSummary(samples)   ## below is from 20,000 samples:
##                 Mean     Median    St.Dev.  95%CI_low  95%CI_upp
## Phi[1]     0.3811684  0.3305808  0.2756644 0.00981843  0.9430733
## Phi[2]     0.7704439  0.8142343  0.1836984 0.32821171  0.9940873
## beta[1]    2.7081528  2.3769716  1.3756124 0.96733862  6.3478021
## beta[2]    1.5291659  1.0696473  1.4673276 0.22759882  6.1198698
## dmean[1]  54.1283388 56.3582949 28.3435779 3.51223792 97.7677991
## dmean[2]  13.0356277  8.3865147 14.1630269 1.07660783 55.8173397
## kappa[1]  27.1158301 27.6161065 13.4763350 3.36387826 48.7495260
## kappa[2]  16.4770625 10.1651176 14.8558692 1.75930337 47.9956812
## lambda0    0.2763837  0.2305539  0.1828862 0.06631507  0.7484122
## phi[1, 1]  0.4519105  0.4280024  0.2695279 0.02887790  0.9560593
## phi[2, 1]  0.8138991  0.8542294  0.1561065 0.42564883  0.9954638
## phi[1, 2]  0.4519105  0.4280024  0.2695279 0.02887790  0.9560593
## phi[2, 2]  0.8138991  0.8542294  0.1561065 0.42564883  0.9954638
## phi[1, 3]  0.4765415  0.4607808  0.2657029 0.03930334  0.9598023
## phi[2, 3]  0.8272435  0.8660133  0.1470100 0.45846674  0.9958575
## sigma[1]   8.6082784  8.1031316  2.3216378 5.55651166 14.7026496
## sigma[2]  10.8578539 10.4370309  2.6921622 6.14002148 17.4127553


samplesPlot(samples)

library(coda)
ess <- apply(samples, 2, effectiveSize)
ess_scr1_reduced <- ess
ess
##    Phi[1]    Phi[2]   beta[1]   beta[2]  dmean[1]  dmean[2]  kappa[1]  kappa[2] 
## 1355.1139 2736.5759 1121.9367  169.3445  546.1404  650.6415 1981.9109  339.0159 
##   lambda0 phi[1, 1] phi[2, 1] phi[1, 2] phi[2, 2] phi[1, 3] phi[2, 3]  sigma[1] 
##  390.8246 1309.6811 2736.5735 1309.6811 2736.5735 1312.4055 2737.1066  278.2343 
##  sigma[2] 
##  106.7198 

eff <- ess / runtime
eff_scr1_reduced <- eff
eff
##     Phi[1]     Phi[2]    beta[1]    beta[2]   dmean[1]   dmean[2]   kappa[1] 
## 101.324187 204.618460  83.889129  12.662180  40.835855  48.649580 148.190869 
##   kappa[2]    lambda0  phi[1, 1]  phi[2, 1]  phi[1, 2]  phi[2, 2]  phi[1, 3] 
##  25.348801  29.222622  97.927099 204.618284  97.927099 204.618284  98.130805 
##  phi[2, 3]   sigma[1]   sigma[2] 
## 204.658142  20.804052   7.979622 





################################################################################################
################################################################################################
## NEW NIMBLE MODEL #2: dSCR2, or scr2
## using a custom distribution to calculate likelihood (EXACTLY),
## but also removing the unknown stochastic 'z' latent states
## (integrating over them in the likelihood calculation)
## using REDUCED voles dataset
################################################################################################
################################################################################################


library(nimble)

code_dSCR2 <- nimbleCode({
    ## space use and recapture probability parameters
    PL ~ dunif(0.01, 0.99)
    lambda0 <- -log(1-PL)
    for(sex in 1:2) {
        kappa[sex] ~ dunif(0,   50)
        sigma[sex] ~ dunif(0.1, 20)
        beta[sex]  ~ dunif(0.1, 10)    # misnomer: beta[1] is coeff of tod, beta[2] is coeff of sex
        for(TOD in 1:2) {
            lambda[TOD, sex] <- lambda0 * beta[1]^(TOD-1) * beta[2]^(sex-1)
        }
        ## survival parameters
        Phi[sex] ~ dunif(0, 1)
        for(k in 1:(nPrimary-1)) {
            phi[sex, k] <- Phi[sex]^dt[k]
        }
        ## dispersal parameters
        dmean[sex] ~ dunif(0, 100)
        dlambda[sex] <- 1/dmean[sex]
    }
    for(i in 1:nInd) {
        S[i, 1, first[i]] ~ dunif(xlow[i], xupp[i])  # initial center of activity (x)
        S[i, 2, first[i]] ~ dunif(ylow[i], yupp[i])  # initial center of activity (y)
        for(k in first[i]:last[i]) {
            D[i, k, 1:R] <- sqrt((S[i, 1, k] - X[1:R, 1])^2 + (S[i, 2, k] - X[1:R, 2])^2)
            g[i, k, 1:R] <- exp(-(D[i, k, 1:R]/sigma[gr[i]])^kappa[gr[i]])  # trap exposure
            G[i, k] <- sum(g[i, k, 1:R])                                    # total trap exposure
        }
        for(k in first[i]:(last[i]-1)) {
            theta[i, k] ~ dunif(-3.141593, 3.141593)   # dispersal direction
            d[i, k] ~ dexp(dlambda[gr[i]])
            S[i, 1, k+1] <- S[i, 1, k] + d[i, k] * cos(theta[i, k])
            S[i, 2, k+1] <- S[i, 2, k] + d[i, k] * sin(theta[i, k])
        }
        ## likelihood
        H[i, 1:nSecondary, 1:nPrimary] ~ dSCR2(
            first = first[i], last = last[i], J = J[i,1:nPrimary],
            lambda = lambda[1:2,gr[i]], tod = tod[1:nPrimary,1:nSecondary],
            g = g[i,1:nPrimary,1:R], G = G[i,1:nPrimary],
            z = z[i,1:nPrimary], phi = phi[gr[i],1:(nPrimary-1)])
    }
})

## define custom distribution
dSCR2 <- nimbleFunction(
    run = function(x = double(2),
        first = double(), last = double(), J = double(1),
        lambda = double(1), tod = double(2),
        g = double(2), G = double(1), z = double(1), phi = double(1),
        log = double()) {
        pAlive <- 1
        lp <- 0
        ## probability of surviving from k to (k+1): phi[k]
        for(k in first:last) {    # primary session
            if(z[k] == 0)   pAlive <- pAlive * phi[k-1]
            for(j in 1:J[k]) {    # secondary session
                ## PnoCapture|alive = exp(-lambda[tod[k,j]] * G[k])
                PnoCapture <- (1-pAlive) + (pAlive)*exp(-lambda[tod[k,j]] * G[k])
                if(x[j,k] == 1) {    # not captured
                    lp <- lp + log(PnoCapture)
                } else {             # captured
                    lp <- lp + log(1-PnoCapture) + log(g[k, x[j,k]-1]) - log(G[k])
                }
            }
        }
        returnType(double())
        if(log) return(lp) else return(exp(lp))
    }
)

rSCR2 <- nimbleFunction(
    run = function(n = integer(),
        first = double(), last = double(), J = double(1),
        lambda = double(1), tod = double(2),
        g = double(2), G = double(1), z = double(1), phi = double(1)) {
        x <- array(1, c(dim(tod)[2], dim(tod)[1]))
        returnType(double(2))
        return(x)
    }
)

registerDistributions(list(
    dSCR2 = list(
        BUGSdist = 'dSCR2(first, last, J, lambda, tod, g, G, z, phi)',
        types = c('value = double(2)', 'first = double()', 'last = double()', 'J = double(1)', 'lambda = double(1)', 'tod = double(2)', 'g = double(2)', 'G = double(1)', 'z = double(1)', 'phi = double(1)'),
        discrete = TRUE,
        mixedSizes = TRUE
    )
))


load('~/github/scr_ergon/analysis/volesData_reduced.RData')
constants_dSCR2 <- constants_reduced[c('R', 'J', 'tod', 'first', 'X',
                                       'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR2$z          <- inits_reduced$z
constants_dSCR2$last       <- constants_reduced$K
constants_dSCR2$nInd       <- dim(constants_reduced$H)[1]
constants_dSCR2$nPrimary   <- dim(constants_reduced$H)[3]
constants_dSCR2$nSecondary <- dim(constants_reduced$H)[2]
data_dSCR2 <- list(H=constants_reduced$H)
zInd <- which(names(inits_reduced) == 'z')
inits_dSCR2 <- inits_reduced[-zInd]


system.time(Rmodel <- nimbleModel(code_dSCR2, constants_dSCR2, data_dSCR2, inits_dSCR2, calculate = FALSE))  ## 1.8 seconds

system.time(lp <- Rmodel$calculate())    ## .05 second
lp  ##  -181.7042

## NOTE: 'log.lambda0' changed to 'lambda0'
monitors = c('kappa', 'sigma', 'lambda0', 'beta', 'dmean', 'phi', 'Phi')

system.time(conf <- configureMCMC(Rmodel, monitors = monitors))  ## 0.1 seconds
conf$printSamplers()
system.time(Rmcmc <- buildMCMC(conf))  ## 0.25 seconds

system.time(Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE))   ## 15 seconds
system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE))   ## 4.5 seconds

system.time(lp <- Cmodel$calculate())    ## 0.02 seconds
lp    ##  -181.7042

niter <- 20000

set.seed(0)
t1 <- Sys.time()
samples <- runMCMC(Cmcmc, niter)
t2 <- Sys.time()
runtime <- as.numeric(difftime(t2, t1, units = 'secs'))
runtime_scr2_reduced <- runtime
runtime  ## 12 seconds

samples_scr2_reduced <- samples

colnames(samples)
dim(samples)   ##  20000    17
samplesSummary(samples)   ## below is from 20,000 samples:
##                 Mean     Median    St.Dev.    95%CI_low  95%CI_upp
## Phi[1]     0.3816868  0.3105915  0.2962007  0.009537986  0.9580860
## Phi[2]     0.5057092  0.5097251  0.2844527  0.026597307  0.9720209
## beta[1]    2.6979009  2.3781258  1.3852046  0.947767291  6.3055778
## beta[2]    1.3161554  0.9673094  1.1038177  0.226164121  4.5117192
## dmean[1]  61.6021827 64.9642341 25.9022343 10.472399253 98.3608487
## dmean[2]  13.3751406  8.7214269 14.2178584  0.851043704 55.7681167
## kappa[1]  26.4572806 26.6914398 13.6594690  2.792589396 48.7093599
## kappa[2]  18.6381782 14.0538027 14.7056469  2.048842025 47.9521903
## lambda0    0.2720497  0.2314629  0.1752559  0.061681372  0.7151843
## phi[1, 1]  0.4484283  0.4080174  0.2893544  0.028243388  0.9677060
## phi[2, 1]  0.5720607  0.5965190  0.2682250  0.061997694  0.9784785
## phi[1, 2]  0.4484283  0.4080174  0.2893544  0.028243388  0.9677060
## phi[2, 2]  0.5720607  0.5965190  0.2682250  0.061997694  0.9784785
## phi[1, 3]  0.4719149  0.4410956  0.2852270  0.038514094  0.9704723
## phi[2, 3]  0.5942837  0.6239290  0.2613254  0.078956222  0.9803314
## sigma[1]   8.5912832  8.0562907  2.3264091  5.537397452 14.9186639
## sigma[2]  10.7168299 10.2269220  2.4561655  6.951261343 16.9326225


samplesPlot(samples)

library(coda)
ess <- apply(samples, 2, effectiveSize)
ess_scr2_reduced <- ess
ess
##    Phi[1]    Phi[2]   beta[1]   beta[2]  dmean[1]  dmean[2]  kappa[1]  kappa[2] 
## 1916.0076 4602.7947 1120.5996  240.0841  689.4887  707.1502 1688.4281  479.1460 
##   lambda0 phi[1, 1] phi[2, 1] phi[1, 2] phi[2, 2] phi[1, 3] phi[2, 3]  sigma[1] 
##  271.5970 1788.8666 4614.4655 1788.8666 4614.4655 1753.6719 4620.5551  230.8703 
##  sigma[2] 
##  138.6466 

eff <- ess / runtime
eff_scr2_reduced <- eff
eff
##    Phi[1]    Phi[2]   beta[1]   beta[2]  dmean[1]  dmean[2]  kappa[1]  kappa[2] 
## 160.46019 385.47096  93.84703  20.10636  57.74272  59.22182 141.40105  40.12711 
##   lambda0 phi[1, 1] phi[2, 1] phi[1, 2] phi[2, 2] phi[1, 3] phi[2, 3]  sigma[1] 
##  22.74548 149.81249 386.44836 149.81249 386.44836 146.86503 386.95834  19.33473 
##  sigma[2] 
##  11.61126 











system.time(Rmodel <- nimbleModel(code_dSCR3, constants_dSCR3, data_dSCR3, inits_dSCR3, calculate = FALSE))  ## 1.8 seconds

Rmodel$calculate()

system.time(lp <- Rmodel$calculate())    ## .05 second
lp  ##  -181.7042

## NOTE: 'log.lambda0' changed to 'lambda0'
monitors = c('kappa', 'sigma', 'lambda0', 'beta', 'dmean', 'phi', 'Phi')

system.time(conf <- configureMCMC(Rmodel, monitors = monitors))  ## 0.1 seconds
conf$printSamplers()
system.time(Rmcmc <- buildMCMC(conf))  ## 0.25 seconds

system.time(Cmodel <- compileNimble(Rmodel, showCompilerOutput = TRUE))   ## 15 seconds
system.time(Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = TRUE))   ## 4.5 seconds

system.time(lp <- Cmodel$calculate())    ## 0.02 seconds
lp    ##  -181.7042

niter <- 20000




