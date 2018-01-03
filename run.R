
if(Sys.info()['nodename'] == 'gandalf') library(nimble, lib.loc = '~/Documents/') else library(nimble)
library(coda)

system('rm data/*')
source('create_data_voles.R')

##reduced <- TRUE
reduced <- FALSE

runVoles <- TRUE
##runVoles <- FALSE

##niter <- 3000
##niter <- 10000
niter <- 20000




runComparison <- function(modelInfoFile, reduced, name, MCMCs, niter, MCMCdefs = list(), add = FALSE, saveFile, verbose = TRUE) {
    if(verbose) message(paste0('running ', MCMCs, ' on ', modelInfoFile, '...'))
    modelInfoFileToLoad <- modelInfoFile
    if(reduced) modelInfoFileToLoad <- paste0(modelInfoFileToLoad, '_reduced')
    modelInfoFileToLoad <- paste0('data/modelInfo_', modelInfoFileToLoad, '.RData')
    load(modelInfoFileToLoad)
    outList <- if(add) dget(saveFile) else list()
    if(length(MCMCs) > 1) stop('only one MCMC at a time, please')
    out <- compareMCMCs(modelInfo = modelInfo, MCMCs = MCMCs, MCMCdefs = MCMCdefs,
                        monitors = modelInfo$monitors, niter = niter)[[1]]
    out <- rename_MCMC_comparison_method(MCMCs, name, out)
    outList[[name]] <- out
    if(!missing(saveFile)) dput(outList, file = saveFile)
    if(verbose) message(paste0('finished running ', name, ' on ', modelInfoFile))
    return(invisible(outList))
}

makePages <- function(saveFile, dir, open = TRUE) {
    outList <- dget(saveFile)
    results <- do.call(combine_MCMC_comparison_results, unname(outList))
    pagesDir <- paste0('pages/', dir, '/')
    make_MCMC_comparison_pages(results, dir = pagesDir, pageComponents = list(timing = TRUE, efficiencySummary = FALSE, efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, efficiencyDetails = TRUE, posteriorSummary = TRUE))
    if(open) system(paste0('open ', pagesDir, 'MCMCresults.html'))
}

makeRmodel <- function(modelInfoFile, reduced) {
    modelInfoFileToLoad <- modelInfoFile
    if(reduced) modelInfoFileToLoad <- paste0(modelInfoFileToLoad, '_reduced')
    modelInfoFileToLoad <- paste0('data/modelInfo_', modelInfoFileToLoad, '.RData')
    load(modelInfoFileToLoad)
    Rmodel <- nimbleModel(code = modelInfo$code, constants = modelInfo$constants, data = modelInfo$data, inits = modelInfo$inits)
    return(Rmodel)
}

if(runVoles) {
    saveFile <- 'results/voles.rda'
    ##
    runComparison(modelInfoFile = 'voles', name = 'nimble', MCMCs = 'nimble', reduced = reduced, niter = niter, saveFile = saveFile)
    ##runComparison(modelInfoFile = 'voles', name = 'jags', MCMCs = 'jags', reduced = reduced, niter = niter, saveFile = saveFile, add = TRUE)
    runComparison(modelInfoFile = 'volesSCR1', name = 'SCR1', MCMCs = 'nimble', reduced = reduced, niter = niter, saveFile = saveFile, add = TRUE)
    runComparison(modelInfoFile = 'volesSCR2', name = 'SCR2', MCMCs = 'nimble', reduced = reduced, niter = niter, saveFile = saveFile, add = TRUE)
    ##
    runComparison(modelInfoFile = 'volesSCR2', name = 'dth', MCMCs = 'X', MCMCdefs = list(X = quote({
        conf <- configureMCMC(Rmodel)
        conf$removeSamplers(c('d', 'theta'))
        for(nn in Rmodel$expandNodeNames('theta'))
            conf$addSampler(c(nn, gsub('theta','d',nn)), 'RW_block', silent = TRUE)
        conf$printSamplers()
        conf
    })), reduced = reduced, niter = niter, saveFile = saveFile, add = TRUE)
    ##
    runComparison(modelInfoFile = 'volesSCR2', name = 'skb', MCMCs = 'X', MCMCdefs = list(X = quote({
        conf <- configureMCMC(Rmodel)
        conf$removeSamplers(c('sigma', 'kappa', 'beta'))
        conf$addSampler(c('sigma[1]', 'kappa[1]'), 'RW_block', silent = TRUE)
        conf$addSampler(c('sigma[2]', 'kappa[2]'), 'RW_block', silent = TRUE)
        conf$addSampler(c('beta[1]',  'beta[2]'),  'RW_block', silent = TRUE)
        conf$printSamplers()
        conf
    })), reduced = reduced, niter = niter, saveFile = saveFile, add = TRUE)
    ##
    runComparison(modelInfoFile = 'volesSCR2', name = 'SLC', MCMCs = 'X', MCMCdefs = list(X = quote({
        conf <- configureMCMC(Rmodel)
        nns <- c('sigma', 'kappa', 'beta', 'PL', 'Phi', 'dmean')
        conf$removeSamplers(nns)
        for(nn in Rmodel$expandNodeNames(nns))
            conf$addSampler(nn, 'slice')
        conf$printSamplers()
        conf
    })), reduced = reduced, niter = niter, saveFile = saveFile, add = TRUE)
    ##makePages(saveFile = saveFile, dir = 'voles', open = TRUE)
    makePages(saveFile = saveFile, dir = 'voles', open = FALSE)
}


if(FALSE) {
    modelInfoFile <- 'volesSCR2'
    reduced <- FALSE
    reduced <- TRUE
    Rmodel <- makeRmodel(modelInfoFile = modelInfoFile, reduced = reduced)
    conf <- configureMCMC(Rmodel)
    conf$printSamplers()
    Rmodel$getNodeNames(stochOnly = TRUE, includeData = FALSE)
}









