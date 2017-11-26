



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



