


################################################################################################
################################################################################################
## saving a standardized form of the FULL voles dataset
## including constants, data, and initial values,
## for NIMLBE and JAGS
################################################################################################
################################################################################################

1

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

volesData <- dget('original_files/dryad/ErgonAndGardner2013.rdat')

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


code <- nimbleCode({
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



## original monitors from Tor:
##monitors <- c('kappa', 'sigma', 'lambda0', 'beta', 'dmean', 'phi', 'Phi')
monitors <- c('kappa', 'sigma', 'lambda0', 'beta', 'dmean', 'Phi')



modelInfo <- list(code = code,
                  constants = constants,
                  data = data,
                  inits = inits,
                  monitors = monitors)

save(modelInfo, file = 'data/modelInfo_voles.RData')



################################################################################################
################################################################################################
## generating a smaller (reduced) version of the voles dataset,
## including constants, data, and initial values,
## for testing other NIMBLE models
################################################################################################
################################################################################################


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

##identical(names(volesData), names(volesData_reduced))



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

##names(inits_reduced)
#### [1] "kappa" "sigma" "PL"    "beta"  "dmean" "Phi"   "theta" "d"     "z"    "S"    
## 
##dim(inits_reduced$S)
##for(i in 1:dim(inits_reduced$S)[1]) {print('=========================='); print(i); print(inits_reduced$S[i,,])}

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

data_reduced <- list(Ones = array(1, dim(volesData_reduced$H)))




modelInfo <- list(code = code,
                  constants = constants_reduced,
                  data = data_reduced,
                  inits = inits_reduced,
                  monitors = monitors)


save(modelInfo, file = 'data/modelInfo_voles_reduced.RData')




################################################################################################
################################################################################################
## NEW NIMBLE MODEL #1: dSCR1, or scr1
## using a custom distribution to calculate likelihood (EXACTLY)
## using voles dataset
################################################################################################
################################################################################################





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



load('data/modelInfo_voles.RData')
constants <- modelInfo$constants
data <- modelInfo$data
inits <- modelInfo$inits

constants_dSCR1 <- constants[c('R', 'J', 'tod', 'first', 'X',
                                       'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR1$last       <- constants$K
constants_dSCR1$nInd       <- dim(constants$H)[1]
constants_dSCR1$nPrimary   <- dim(constants$H)[3]
constants_dSCR1$nSecondary <- dim(constants$H)[2]
data_dSCR1 <- list(H=constants$H)
inits_dSCR1 <- inits


modelInfo <- list(code = code_dSCR1,
                  constants = constants_dSCR1,
                  data = data_dSCR1,
                  inits = inits_dSCR1,
                  monitors = monitors)


save(modelInfo, file = 'data/modelInfo_volesSCR1.RData')



##out_dSCR1 <- compareMCMCs(modelInfo=modelInfo_dSCR1, monitors=monitors, niter=niter)[[1]]
##out_dSCR1 <- rename_MCMC_comparison_method('nimble', 'SCR1', out_dSCR1)
##outList$SCR1 <- out_dSCR1
##save(outList, file = saveFile)
##message('finished SCR1')





################################################################################################
################################################################################################
## NEW NIMBLE MODEL #1: dSCR1, or scr1
## using a custom distribution to calculate likelihood (EXACTLY)
## using REDUCED voles dataset
################################################################################################
################################################################################################


load('data/modelInfo_voles_reduced.RData')
constants <- modelInfo$constants
data <- modelInfo$data
inits <- modelInfo$inits

constants_dSCR1 <- constants[c('R', 'J', 'tod', 'first', 'X',
                                       'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR1$last       <- constants$K
constants_dSCR1$nInd       <- dim(constants$H)[1]
constants_dSCR1$nPrimary   <- dim(constants$H)[3]
constants_dSCR1$nSecondary <- dim(constants$H)[2]
data_dSCR1 <- list(H=constants$H)
inits_dSCR1 <- inits


modelInfo <- list(code = code_dSCR1,
                  constants = constants_dSCR1,
                  data = data_dSCR1,
                  inits = inits_dSCR1,
                  monitors = monitors)


save(modelInfo, file = 'data/modelInfo_volesSCR1_reduced.RData')





################################################################################################
################################################################################################
## NEW NIMBLE MODEL #2: dSCR2, or scr2
## using a custom distribution to calculate likelihood (EXACTLY),
## but also removing the unknown stochastic 'z' latent states
## (integrating over them in the likelihood calculation)
## using voles dataset
################################################################################################
################################################################################################




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
        pDead <- 0
        lp <- 0
        ## probability of surviving from k to (k+1): phi[k]
        for(k in first:last) {
            if(z[k] == 1) {    # known to be alive
                if(k > first)           # survived
                    lp <- lp + log(phi[k-1])
                for(j in 1:J[k]) {
                    pNoCaptureGivenAlive <- exp(-lambda[tod[k,j]] * G[k])
                    if(x[j,k] == 1) {   # not captured
                        lp <- lp + log(pNoCaptureGivenAlive)
                    } else {            # captured
                        lp <- lp + log(1-pNoCaptureGivenAlive) + log(g[k, x[j,k]-1]) - log(G[k])
                    }
                }
            } else {           # could be dead or alive
                pTheseNonSightings <- 1
                for(j in 1:J[k]) {
                    pNoCaptureGivenAlive <- exp(-lambda[tod[k,j]] * G[k])
                    pTheseNonSightings <- pTheseNonSightings * pNoCaptureGivenAlive
                }
                pAlive_new <- phi[k-1] * pAlive
                pDead_new <- (1-phi[k-1]) * pAlive + pDead
                L <- pAlive_new * pTheseNonSightings + pDead_new
                pAlive <- (pAlive_new * pTheseNonSightings) / L
                pDead <- pDead_new / L
                lp <- lp + log(L)
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


load('data/modelInfo_voles.RData')
constants <- modelInfo$constants
data <- modelInfo$data
inits <- modelInfo$inits


constants_dSCR2 <- constants[c('R', 'J', 'tod', 'first', 'X',
                               'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR2$z          <- inits$z
constants_dSCR2$last       <- constants$K
constants_dSCR2$nInd       <- dim(constants$H)[1]
constants_dSCR2$nPrimary   <- dim(constants$H)[3]
constants_dSCR2$nSecondary <- dim(constants$H)[2]
data_dSCR2 <- list(H=constants$H)
zInd <- which(names(inits) == 'z')
inits_dSCR2 <- inits[-zInd]


modelInfo <- list(code = code_dSCR2,
                  constants = constants_dSCR2,
                  data = data_dSCR2,
                  inits = inits_dSCR2,
                  monitors = monitors)


save(modelInfo, file = 'data/modelInfo_volesSCR2.RData')



################################################################################################
################################################################################################
## NEW NIMBLE MODEL #2: dSCR2, or scr2
## using a custom distribution to calculate likelihood (EXACTLY),
## but also removing the unknown stochastic 'z' latent states
## (integrating over them in the likelihood calculation)
## using REDUCED voles dataset
################################################################################################
################################################################################################


load('data/modelInfo_voles_reduced.RData')
constants <- modelInfo$constants
data <- modelInfo$data
inits <- modelInfo$inits


constants_dSCR2 <- constants[c('R', 'J', 'tod', 'first', 'X',
                               'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR2$z          <- inits$z
constants_dSCR2$last       <- constants$K
constants_dSCR2$nInd       <- dim(constants$H)[1]
constants_dSCR2$nPrimary   <- dim(constants$H)[3]
constants_dSCR2$nSecondary <- dim(constants$H)[2]
data_dSCR2 <- list(H=constants$H)
zInd <- which(names(inits) == 'z')
inits_dSCR2 <- inits[-zInd]


modelInfo <- list(code = code_dSCR2,
                  constants = constants_dSCR2,
                  data = data_dSCR2,
                  inits = inits_dSCR2,
                  monitors = monitors)


save(modelInfo, file = 'data/modelInfo_volesSCR2_reduced.RData')



























