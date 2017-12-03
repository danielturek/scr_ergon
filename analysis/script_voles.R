
setwd('~/github/scr_ergon/analysis')

##niter <- 3000
niter <- 20000

saveFile <- 'results_scr3.RData'
monitors <- c("kappa", "sigma", "lambda0", "beta", "dmean", "phi", "Phi")



load('volesData_reduced.RData')
constants <- constants_reduced
data <- data_reduced
inits <- inits_reduced


##load('volesData.RData')






message('starting NIMBLE')




if(Sys.info()['nodename'] == 'gandalf') library(nimble, lib.loc = '~/Documents/') else library(nimble)



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





##modelInfo <- list(code = code, constants = constants, data=data, inits = inits, name = 'jags')
##out_jags <- compareMCMCs(modelInfo=modelInfo, MCMCs = 'jags', monitors=monitors, niter=niter)
##save(out_jags, file = saveFile)
## 
##message('finished JAGS')


##modelInfo <- list(code = code, constants = constants, data=data, inits = inits, name = 'nimble')
##out_nimble <- compareMCMCs(modelInfo=modelInfo, monitors=monitors, niter=niter)
##save(out_jags, out_nimble, file = saveFile)
## 
##message('finished NIMBLE')




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


constants_dSCR1 <- constants[c('R', 'J', 'tod', 'first', 'X',
                                       'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR1$last       <- constants$K
constants_dSCR1$nInd       <- dim(constants$H)[1]
constants_dSCR1$nPrimary   <- dim(constants$H)[3]
constants_dSCR1$nSecondary <- dim(constants$H)[2]
data_dSCR1 <- list(H=constants$H)
inits_dSCR1 <- inits




modelInfo_dSCR1 <- list(code = code_dSCR1, constants = constants_dSCR1, data=data_dSCR1, inits = inits_dSCR1, name = 'SCR1')
out_dSCR1 <- compareMCMCs(modelInfo=modelInfo_dSCR1, monitors=monitors, niter=niter)
##save(out_jags, out_nimble, out_dSCR1, file = saveFile)
save(out_dSCR1, file = saveFile)

message('finished SCR1')






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






modelInfo_dSCR2 <- list(code = code_dSCR2, constants = constants_dSCR2, data=data_dSCR2, inits = inits_dSCR2, name = 'SCR2')
out_dSCR2 <- compareMCMCs(modelInfo=modelInfo_dSCR2, monitors=monitors, niter=niter)
##save(out_jags, out_nimble, out_dSCR1, out_dSCR2, file = saveFile)
save(out_dSCR1, out_dSCR2, file = saveFile)

message('finished SCR2')




code_dSCR3 <- nimbleCode({
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
##############theta[i, k] ~ dunif(-3.141593, 3.141593)   # dispersal direction
##############d[i, k] ~ dexp(dlambda[gr[i]])
##############S[i, 1, k+1] <- S[i, 1, k] + d[i, k] * cos(theta[i, k])
##############S[i, 2, k+1] <- S[i, 2, k] + d[i, k] * sin(theta[i, k])
            S[i, 1:2, k+1] ~ dS3(S = S[i, 1:2, k], lambda = dlambda[gr[i]])
        }
        ## likelihood
        H[i, 1:nSecondary, 1:nPrimary] ~ dSCR3(
            first = first[i], last = last[i], J = J[i,1:nPrimary],
            lambda = lambda[1:2,gr[i]], tod = tod[1:nPrimary,1:nSecondary],
            g = g[i,1:nPrimary,1:R], G = G[i,1:nPrimary],
            z = z[i,1:nPrimary], phi = phi[gr[i],1:(nPrimary-1)])
    }
})

dS3 <- nimbleFunction(
    run = function(x = double(1), S = double(1), lambda = double(), log = double()) {
        D <- sqrt((x[1]-S[1])^2 + (x[2]-S[2])^2)
        lp <- dexp(D, lambda, log = TRUE)
        returnType(double())
        return(lp)
    }
)

rS3 <- nimbleFunction(
    run = function(n = integer(), S = double(1), lambda = double()) {
        x <- c(1,1)
        returnType(double(1))
        return(x)
    }
)

registerDistributions(list(
    dS3 = list(
        BUGSdist = 'dS3(S, lambda)',
        types = c('value = double(1)', 'S = double(1)', 'lambda = double()'),
        discrete = TRUE,
        mixedSizes = TRUE
    )
))

## define custom distribution
dSCR3 <- nimbleFunction(
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

rSCR3 <- nimbleFunction(
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
    dSCR3 = list(
        BUGSdist = 'dSCR3(first, last, J, lambda, tod, g, G, z, phi)',
        types = c('value = double(2)', 'first = double()', 'last = double()', 'J = double(1)', 'lambda = double(1)', 'tod = double(2)', 'g = double(2)', 'G = double(1)', 'z = double(1)', 'phi = double(1)'),
        discrete = TRUE,
        mixedSizes = TRUE
    )
))


constants_dSCR3 <- constants[c('R', 'J', 'tod', 'first', 'X',
                               'dt', 'gr', 'xlow', 'xupp', 'ylow', 'yupp')]
constants_dSCR3$z          <- inits$z
constants_dSCR3$last       <- constants$K
constants_dSCR3$nInd       <- dim(constants$H)[1]
constants_dSCR3$nPrimary   <- dim(constants$H)[3]
constants_dSCR3$nSecondary <- dim(constants$H)[2]
data_dSCR3 <- list(H=constants$H)
zInd <- which(names(inits) == 'z')
thetaInd <- which(names(inits) == 'theta')
dInd <- which(names(inits) == 'd')
inits_dSCR3 <- inits[c(-zInd, -thetaInd, -dInd)]
Sinit <- inits_dSCR3$S
for(i in 1:constants_dSCR3$nInd) {
    if(constants_dSCR3$last[i] > constants_dSCR3$first[i]) {
        for(k in (constants_dSCR3$first[i]+1):constants_dSCR3$last[i]) {
            Sinit[i,1:2,k] <- Sinit[i,1:2,constants_dSCR3$first[i]]
        }
    }
}
inits_dSCR3$S <- Sinit



modelInfo_dSCR3 <- list(code = code_dSCR3, constants = constants_dSCR3, data=data_dSCR3, inits = inits_dSCR3, name = 'SCR3')
out_dSCR3 <- compareMCMCs(modelInfo=modelInfo_dSCR3, monitors=monitors, niter=niter)
##save(out_jags, out_nimble, out_dSCR1, out_dSCR2, out_dSCR3, file = saveFile)
save(out_dSCR1, out_dSCR2, out_dSCR3, file = saveFile)

message('finished SCR3')







if(FALSE) {

    ## making comparison pages for everything (in 'results_all.RData')
    ## jags, nimble, dSCR1, dSCR2
    setwd('~/github/scr_ergon/analysis')
    ls()
    load('results_all.RData')
    ls()
    library(nimble)
    ## rename results
    out_dSCR1[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR1', out_dSCR1[[1]])
    out_dSCR2[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR2', out_dSCR2[[1]])
    ## combine results
    results <- combine_MCMC_comparison_results(out_jags[[1]], out_nimble[[1]], out_dSCR1[[1]], out_dSCR2[[1]])
    ## make comparison pages
    make_MCMC_comparison_pages(results, dir = 'pages', pageComponents = list(timing = TRUE, efficiencySummary = FALSE, efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, efficiencyDetails = TRUE, posteriorSummary = TRUE))
    system('open pages/MCMCresults.html')
    save(out_jags, out_nimble, out_dSCR1, out_dSCR2, results, file = 'results_all.RData')

    ## without jags:
    setwd('~/github/scr_ergon/analysis')
    load('results_reduced2.RData')
    library(nimble)
    ## rename results
    out_dSCR1[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR1', out_dSCR1[[1]])
    out_dSCR2[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR2', out_dSCR2[[1]])
    ## combine results
    results <- combine_MCMC_comparison_results(out_nimble[[1]], out_dSCR1[[1]], out_dSCR2[[1]])
    ## make comparison pages
    make_MCMC_comparison_pages(results, dir = 'pages', pageComponents = list(timing = TRUE, efficiencySummary = FALSE, efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, efficiencyDetails = TRUE, posteriorSummary = TRUE))

    ## adding dSCR2 results (from 'results_reduced3.RData') to
    ## nimble and dSCR1 results (in 'results_reduced2.RData')
    setwd('~/github/scr_ergon/analysis')
    load('results_reduced3.RData')
    ls()
    ## remove everything from workspace
    rrr
    ls()
    load('results_reduced2.RData')
    ls()
    load('results_reduced3.RData')
    save(out_nimble, out_dSCR1, out_dSCR2, file = 'results_reduced2.RData')

    ## adding jags results (in 'results_jags')
    ## into nimble, dSCR1, and dSCR2 results (in 'results_reduced2.RData')
    setwd('~/github/scr_ergon/analysis')
    load('results_jags.RData')
    ls()
    load('results_reduced2.RData')
    ls()
    save(out_jags, out_nimble, out_dSCR1, out_dSCR2, file = 'results_all.RData')
    


    rrr
    setwd('~/github/scr_ergon/analysis')
    ls()
    load('results_noJags.RData')
    ls()
    out_nimble$nimble$timing
    out_dSCR1$SCR1$timing
    out_dSCR2$SCR2$timing
    out_dSCR2_save <- out_dSCR2
    rrr
    load('results_reduced.RData')   ## dSCR2 is wrong here
    ls()
    out_nimble$nimble$timing
    out_dSCR1$SCR1$timing
    out_dSCR2$SCR2$timing
    out_jags$jags$timing
    out_dSCR2 <- out_dSCR2_save
    save(out_jags, out_nimble, out_dSCR1, out_dSCR2, file = 'results_all.RData')


    ## looking at "average" improvement in Efficieny between:
    ## nimble vs. jags
    ## dSCRx vs. nimble
    setwd('~/github/scr_ergon/analysis')
    load('results_all.RData')
    unlist(out_nimble$nimble$efficiency)
    unlist(out_jags$jags$efficiency)
    unlist(out_dSCR1$SCR1$efficiency)
    unlist(out_dSCR2$SCR2$efficiency)
    unlist(out_nimble$nimble$efficiency) / unlist(out_jags$jags$efficiency)
    ## min.nimble mean.nimble 
    ## 4.175435    3.834809 
    unlist(out_dSCR2$SCR2$efficiency) / unlist(out_nimble$nimble$efficiency)
    ## min.SCR2 mean.SCR2 
    ## 1.723233  2.311017

    ## compare and look at dSCR1, dSCR2, and dSCR3
    setwd('~/github/scr_ergon/analysis')
    load('results_scr3.RData')
    ls()
    library(nimble)
    ## rename results
    out_dSCR1[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR1', out_dSCR1[[1]])
    out_dSCR2[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR2', out_dSCR2[[1]])
    out_dSCR3[[1]] <- rename_MCMC_comparison_method('nimble', 'SCR3', out_dSCR3[[1]])
    ## combine results
    results <- combine_MCMC_comparison_results(out_dSCR1[[1]], out_dSCR2[[1]], out_dSCR3[[1]])
    ## make comparison pages
    make_MCMC_comparison_pages(results, dir = 'pages_scr3', pageComponents = list(timing = TRUE, efficiencySummary = FALSE, efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, efficiencyDetails = TRUE, posteriorSummary = TRUE), control = list(mainPageName = 'scr3'))
    system('open pages_scr3/MCMCresults.html')
    


}
















