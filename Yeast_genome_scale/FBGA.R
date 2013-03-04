#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2)

# species involved in reaction
n_c <- 25
metabs <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(metabs) <- c("fu", "bar")
enzyme <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(enzyme) <- c("yfg1", "yfg2")
# Model parameters, and correspondence between parameters and species - i.e. kcat - species 1 matches the activity of the first enzyme catalyzing the reaction. km - species 1 means this the affinity of the first metabolite (column 1) for  the enzyme
assocConst <- data.frame(name = c("yfg1_kcat", "fu_km", "C"), type = c("kcat", "km", "ki"), specie = c(1, 1, 2), priorMean = NA, priorSD = NA)

assocConst$priorSD <- assocConst$priorMean <- c(1, 0.25, 0.3)
trueFlux <- (enzyme[,1] * assocConst$priorMean[1] * metabs[,1] / (metabs[,1] + assocConst$priorMean[2])) #simulate measured fluxes from the rxn equation form assuming that metabolites and enzymes are measured accurately
measuredFlux <- trueFlux*rlnorm(n_c, 0, 0.25) #simulate measured fluxes from the rxn equation form with added lognormal noise

parNum <- length(assocConst[,1]) #how many parameters are there in the model


# function predicting flux from provided paramters
#exampleformula
#enzyme activity
#substrate occupancy

rxnEqn <- as.formula(" ~ I(yfg1_kcat * yfg1 * fu / (fu + fu_km)) + 0")



 ## Genetic algorithm parameters ##
N = 1000 # the number of parameter sets competing for ultimate fit
mu = 0.1 # the lambda mutation rate across all parameters 
generations <- 100
genMeanlogFit <- rep(NA, generations)

#initialize with parameters drawn from the parameters prior
gaConstants <- sapply(1:parNum, function(initConstN){
  rlnorm(N, assocConst[initConstN,colnames(assocConst) == "priorMean"], assocConst[initConstN,colnames(assocConst) == "priorSD"])
  })
colnames(gaConstants) <- assocConst$name

for(genN in 1:generations){
  ### generations of selection (implicitely includes drift - because sampling is proportional to fitness but N is finite), and mutation
  
  indPrior <- apply(sapply(1:parNum, function(parN){
    log(gaConstants[,parN]), assocConst$priorMean[parN]
  dlnorm(gaConstants[,parN], assocConst$priorMean[parN], assocConst$priorSD[parN])
  }), 1, prod) #presolve the prior probability of a model (of parameter values)

  ### Selection ###
  
  indLogFit <- sapply(1:N, indFitnessFxn)
  genMeanlogFit[genN] <- mean(indLogFit) #save the average log fitness of individuals

  rFit <- exp(indLogFit - max(indLogFit)) #fitness relative to most fit individual
  progeny <- rowSums(rmultinom(N, 1, rFit))
  progeny <- unlist(sapply(1:N, function(x){rep(x, progeny[x])}))
  
  gaConstants <- gaConstants[progeny,] #replace current population with selected progeny
  
  ### Mutation ###
  
  mutations <- rpois(N, mu)
  mutations <- sapply(mutations, function(x){min(c(x, length(assocConst[,1])))}) #floor to the maximum number of parameters
  mutInds <- c(1:N)[mutations != 0]
  
  newPar <- sapply(1:length(mutInds), function(mutInd){
    replpar <- sample(1:parNum, mutations[mutInd], replace = TRUE)
    indPars <- gaConstants[mutInd,]
    for(mut in replpar){
      indPars[mut] <- rlnorm(1, assocConst[mut,colnames(assocConst) == "priorMean"], assocConst[mut,colnames(assocConst) == "priorSD"])
      }
    indPars
    })
  
  gaConstants[mutInds,] <- t(newPar)
  
  if(genN/10 == floor(genN/10)){
    print(paste("Generation", genN, "complete.  The mean log fitness is", signif(genMeanlogFit[genN], 5), sep = " "))
    }
  
  }











indFitnessFxn <- function(i){
  indParams <- t(t(rep(1, n_c))) %*% gaConstants[i,]; colnames(indParams) <- assocConst$name
  indDat <- data.frame(metabs, enzyme, indParams)
  
  predFlux <- model.matrix(rxnEqn, data = indDat)[,1] #evaluates the combination of data and parameters in indDat according to the reaction equation
  
  # setting the mean of measured flux to the mean of the predicted flux is equivalent to scaling flux magnitude to account for kcat (when weights are not used this is true, otherwise the weighted mean needs to be determined)
  # Determine how to propogate weights associated with measurement through to get a weight on predicted flux
  
  lpredFluxcent <- log(predFlux) + (mean(log(measuredFlux)) - mean(log(predFlux))) #centered log flux
  
  SSE <- sum((lpredFluxcent - log(measuredFlux))^2) #determine the sum of squared error between predicted log flux and observed log flux
  var_fit <- SSE / (n_c - 1) #correct estimate of variance according to number of fitted parameters + 1
  sum(dnorm(lpredFluxcent, log(measuredFlux), sqrt(var_fit), log = TRUE)) + log(indPrior[i]) # bayes factor for this parameter set is returned
  
  }
