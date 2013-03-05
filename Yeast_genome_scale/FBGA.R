#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###



setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2)

# form a Km hyperparameter using the distribution of [met]/km estimates





# species involved in reaction
n_c <- 100
metabs <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(metabs) <- c("fu", "bar")
enzyme <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(enzyme) <- c("yfg1", "yfg2")
# Model parameters, and correspondence between parameters and species - i.e. kcat - species 1 matches the activity of the first enzyme catalyzing the reaction. km - species 1 means this the affinity of the first metabolite (column 1) for  the enzyme
assocConst <- data.frame(name = c("yfg1_kcat", "fu_km", "C"), type = c("kcat", "km", "ki"), specie = c(1, 1, 2), priorMean = NA, priorSD = NA) # all of these parameters are lognormally distributed, the parameterization that I am using for this is the parameter MLE (prior mean) and the log standard deviation (priorSD)
assocConst$priorSD <- assocConst$priorMean <- c(1, 0.25, 0.8)
assocConstTRUE <- assocConst; assocConstTRUE$priorMean <- c(1, 0.5, 0.6)

trueFlux <- (enzyme[,1] * assocConstTRUE$priorMean[1] * metabs[,1] / (metabs[,1] + assocConstTRUE$priorMean[2])) #simulate measured fluxes from the rxn equation form assuming that metabolites and enzymes are measured accurately
measuredFlux <- trueFlux*rlnorm(n_c, 0, 0.25) #simulate measured fluxes from the rxn equation form with added lognormal noise

parNum <- length(assocConst[,1]) #how many parameters are there in the model
plot(measuredFlux ~ metabs[,1])
plot(measuredFlux ~ enzyme[,1])


# function predicting flux from provided paramters
#exampleformula
#enzyme activity
#substrate occupancy

rxnEqn <- as.formula(" ~ I(yfg1_kcat * yfg1 * fu / (fu + fu_km)) + 0")



 ## Genetic algorithm parameters ##
N = 5000 # the number of parameter sets competing for ultimate fit
mu = 0.1*length(assocConst[,1]) # the lambda mutation rate across all parameters 
generations <- 100
genMeanlogFit <- rep(NA, generations)

#initialize with parameters drawn from the parameters prior
gaConstants <- sapply(1:parNum, function(initConstN){
  #rlnorm(N, assocConst[initConstN,colnames(assocConst) == "priorMean"], assocConst[initConstN,colnames(assocConst) == "priorSD"])
  exp(rnorm(N, log(assocConst[initConstN,colnames(assocConst) == "priorMean"]), assocConst[initConstN,colnames(assocConst) == "priorSD"]))
  })
colnames(gaConstants) <- assocConst$name

for(genN in 1:generations){
  ### generations of selection (implicitely includes drift - because sampling is proportional to fitness but N is finite), and mutation
  
  ### Mutation ###
  
  mutations <- rpois(N, mu)
  mutations <- sapply(mutations, function(x){min(c(x, length(assocConst[,1])))}) #floor to the maximum number of parameters
  mutInds <- c(1:N)[mutations != 0]
  
  newPar <- sapply(1:length(mutInds), function(mutInd){
    replpar <- sample(1:parNum, mutations[mutInd], replace = TRUE)
    indPars <- gaConstants[mutInd,]
    for(mut in replpar){
      indPars[mut] <- exp(rnorm(1, log(assocConst[mut,colnames(assocConst) == "priorMean"]), assocConst[mut,colnames(assocConst) == "priorSD"]))
        #rlnorm(1, assocConst[mut,colnames(assocConst) == "priorMean"], assocConst[mut,colnames(assocConst) == "priorSD"])
      }
    indPars
    })
  
  gaConstants[mutInds,] <- t(newPar)
  
  ######
  
  indPrior <- apply(sapply(1:parNum, function(parN){
    dnorm(log(gaConstants[,parN]), log(assocConst$priorMean[parN]),  assocConst$priorSD[parN])
    #dlnorm(gaConstants[,parN], assocConst$priorMean[parN], assocConst$priorSD[parN])
  }), 1, prod) #presolve the prior probability of a model (of parameter values)

  ### Selection ###
  
  indLogFit <- sapply(1:N, indFitnessFxn)
  
  rFit <- exp(indLogFit - max(indLogFit)) #fitness relative to most fit individual
  progeny <- rowSums(rmultinom(N, 1, rFit))
  progeny <- unlist(sapply(1:N, function(x){rep(x, progeny[x])}))
  
  genMeanlogFit[genN] <- mean(indLogFit[progeny]) #save the average log fitness of individuals

  gaConstants <- gaConstants[progeny,] #replace current population with selected progeny
  
  if(genN/10 == floor(genN/10)){
    print(paste("Generation", genN, "complete.  The mean log fitness is", signif(genMeanlogFit[genN], 5), sep = " "))
    }
  
  }

plot(genMeanlogFit, pch = 16, col = "RED", cex = 0.3)

finalFit <- t(sapply(1:N, indFitnessFxn, splitFit = TRUE))
colnames(finalFit) <- c("fit", "prior")

# visualize the parameter correlation matrix
library(colorRamps)
library(gplots)
parCorrMat <- cor(log(gaConstants))
#parSpearCorrMat <- cor(log(gaConstants), method = "spearman")
heatmap.2(parCorrMat, symm = TRUE, col = blue2yellow(1000), dendrogram = "row", symkey = TRUE)

# visualize parameter distributions colored by score

par_fit <- data.frame(gaConstants, finalFit)
par_fit_stack <- melt(par_fit, measure.vars = assocConst$name)
par_fit_stack$logVal = log(par_fit_stack$value)
#par_fit_stack2 <- melt(par_fit_stack, measure.vars = c("fit", "prior"))
#colnames(par_fit_stack2) <- c("Parameter", "ParValue", "BayesFactorComponent", "BFVal"); par_fit_stack2$BFVal <- exp(par_fit_stack2$BFVal)
#dist_plotter <- ggplot(par_fit_stack, aes(x =  logVal)) + facet_grid(variable ~ .) + geom_abline
#dist_plotter + geom_histogram(par_fit_stack, aes(x = logVal), binwidth = 0.01)

colnames(par_fit_stack) <- c("Fit", "Prior", "name", "Value", "lnValue")
maxCatVal <- sapply(assocConst$name, function(x){max(table(floor(par_fit_stack$lnValue[par_fit_stack$name == x] * 100)))})
assocConst$CImin <- log(assocConst$priorMean) - 2 * assocConst$priorSD; assocConst$CImax <- log(assocConst$priorMean) + 2 * assocConst$priorSD
assocConst$maxBinCount = maxCatVal


hist_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "none", 
  panel.grid.minor = element_blank(), legend.key.width = unit(6, "line"), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan")) 

dist_plotter <- ggplot() + facet_grid(name ~ ., scales = "free_y") + hist_theme
dist_plotter + geom_histogram(data = par_fit_stack, aes(x = lnValue), binwidth = 0.01) + geom_vline(data = assocConst, aes(xintercept = log(priorMean), colour = "RED", size = 2)) +
  geom_errorbarh(data = assocConst, aes(x = log(priorMean), xmin = CImin, xmax = CImax, y = maxCatVal/2, height = maxCatVal/10, size = 2, colour = "RED"))







#dist_plotter <- ggplot(par_fit_stack2, aes(x =  ParValue, fill = BFVal)) + facet_grid(BayesFactorComponent ~ Parameter)
#dist_plotter + geom_histogram() + scale_fill_gradient(name = "woot", low = "black", high = "firebrick1")

#scale_fill_gradient(name = "Counts", low = "black", high = "firebrick1", trans = "log", breaks = hex_breaks, labels = hex_breaks)

# for pairs of parameters exceeding a correlation cutoff plot a bivariate histogram
# find a good example first
corCO <- .2

hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "gray90"), legend.position = "top", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line")) 

parPoint <- sapply(rownames(parCorrMat), function(x){paste(x, colnames(parCorrMat), sep = "/")})
parAssoc <- parPoint[upper.tri(parCorrMat)][abs(parCorrMat[upper.tri(parCorrMat)]) > corCO]
if(length(parAssoc) != 0){
  
  for(parSet in parAssoc){
    bivariateDF <- data.frame(log(gaConstants[,c(1:parNum)[assocConst$name %in% strsplit(parSet, split = '/')[[1]]]]))
    colnames(bivariateDF) <- c("P1", "P2")
    print(ggplot(bivariateDF, aes(x = P1, y = P2)) + geom_hex() + scale_x_continuous(assocConst$name[assocConst$name %in% strsplit(parSet, split = '/')[[1]]][1], expand = c(0.02,0.02)) +
      scale_y_continuous(assocConst$name[assocConst$name %in% strsplit(parSet, split = '/')[[1]]][2], expand = c(0.02,0.02)) + hex_theme +
      scale_fill_gradient(name = "Counts", low = "black", high = "firebrick1"))
    }
   
     
  }






indFitnessFxn <- function(i, splitFit = FALSE){
  indParams <- t(t(rep(1, n_c))) %*% gaConstants[i,]; colnames(indParams) <- assocConst$name
  indDat <- data.frame(metabs, enzyme, indParams)
  
  predFlux <- model.matrix(rxnEqn, data = indDat)[,1] #evaluates the combination of data and parameters in indDat according to the reaction equation
  
  # setting the mean of measured flux to the mean of the predicted flux is equivalent to scaling flux magnitude to account for kcat (when weights are not used this is true, otherwise the weighted mean needs to be determined)
  # Determine how to propogate weights associated with measurement through to get a weight on predicted flux
  
  lpredFluxcent <- log(predFlux) + (mean(log(measuredFlux)) - mean(log(predFlux))) #centered log flux
  
  SSE <- sum((lpredFluxcent - log(measuredFlux))^2) #determine the sum of squared error between predicted log flux and observed log flux
  var_fit <- SSE / (n_c - 1) #correct estimate of variance according to number of fitted parameters + 1
  
  if(splitFit == FALSE){
  sum(dnorm(lpredFluxcent, log(measuredFlux), sqrt(var_fit), log = TRUE)) + log(indPrior[i]) # bayes factor for this parameter set is returned
    }else{
      c(sum(dnorm(lpredFluxcent, log(measuredFlux), sqrt(var_fit), log = TRUE)), log(indPrior[i]))
      }
  }





fullFitness <- function(gaConstants){
  
  
  
  
  
  
  
  }














