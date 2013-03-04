#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2)

# species involved in reaction
n_c <- 25
metabs <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(metabs) <- c("fu", "bar")
enzyme <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(enzyme) <- c("yfg1", "yfg2")
# Model parameters, and correspondence between parameters and species - i.e. kcat - species 1 matches the activity of the first enzyme catalyzing the reaction. km - species 1 means this the affinity of the first metabolite (column 1) for  the enzyme
assocConst <- data.frame(name = c("yfg1_kcat", "fu_km", "C"), type = c("kcat", "km", "ki"), specie = c(1, 1, 2), priorMean = NA, priorVar = NA)

assocConst$priorMean <- assocConst$priorVar <- c(1, 0.25, 0.3)
trueFlux <- (enzyme[,1] * assocConst$priorMean[1] * metabs[,1] / (metabs[,1] + assocConst$priorMean[2])) #simulate measured fluxes from the rxn equation form assuming that metabolites and enzymes are measured accurately
measuredFlux <- trueFlux*rlnorm(n_c, 0, 0.25) #simulate measured fluxes from the rxn equation form with added lognormal noise
  


# function predicting flux from provided paramters
#exampleformula
#enzyme activity
#substrate occupancy

rxnEqn <- as.formula(" ~ I(yfg1_kcat * yfg1 * fu / (fu + fu_km)) + 0")


## Genetic algorithm parameters ##
N = 1000 # the number of parameter sets competing for ultimate fit
mu = 0.05 # the lambda mutation rate across all parameters 

#initialize with parameters drawn from the parameters prior
gaConstants <- sapply(1:length(assocConst[,1]), function(initConstN){
  rnorm(N, assocConst[initConstN,colnames(assocConst) == "priorMean"], assocConst[initConstN,colnames(assocConst) == "priorVar"])
  })
colnames(gaConstants) <- assocConst$name


indParams <- t(t(rep(1, n_c))) %*% gaConstants[1,]; colnames(indParams) <- assocConst$name
indDat <- data.frame(metabs, enzyme, indParams)

predFlux <- model.matrix(rxnEqn, data = indDat)[,1] #evaluates the combination of data and parameters in indDat according to the reaction equation

# setting the mean of measured flux to the mean of the predicted flux is equivalent to scaling flux magnitude to account for kcat (when weights are not used this is true, otherwise the weighted mean needs to be determined)


# Determine how to propogate weights associated with measurement through to get a weight on predicted flux

log2(predFlux) - log2(measuredFlux)