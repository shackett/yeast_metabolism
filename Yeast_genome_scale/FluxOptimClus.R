#qsub -l 1day -cwd -sync n Rscript FluxOptimClus.R runNum=0

setwd("/Genomics/grid/users/shackett/FBA/FluxParOptim/")

library(nnls)

options(stringsAsFactors = FALSE)

args <- commandArgs()
runNum = as.numeric(unlist(strsplit(args[grep("runNum", args)], "="))[2])

load("paramOptim.Rdata")
run_summary <- list() #### MCMC run output and formatted inputs

markov_pars <- list()
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 10000 #how many total markov samples are desired
markov_pars$burn_in <- 500 #how many initial samples should be skipped

markov_pars <- list()
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 100 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped


run_summary$markov_pars <- markov_pars


########### Functions ##########

par_draw <- function(updates){
  #### update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ####
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    } else if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    }
  }
  draw
}

lik_calc <- function(proposed_params){
  #### determine the likelihood of predicted flux as a function of metabolite abundance and kinetics parameters relative to actual flux ####
  
  par_stack <- rep(1, n_c) %*% t(proposed_params); colnames(par_stack) <- kineticPars$formulaName
  par_stack <- exp(par_stack)
  occupancy_vals <- data.frame(met_abund, par_stack)
  
  predOcc <- model.matrix(occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*enzyme_abund #occupany of enzymes * relative abundance of enzymes
  
  flux_fit <- nnls(enzyme_activity, flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  fit_resid_error <- sqrt(mean((flux_fit$resid - mean(flux_fit$resid))^2))
  
  sum(dnorm(flux, flux_fit$fitted, fit_resid_error, log = TRUE))
  
}

n_c <- length(rxnList_form[[1]]$flux)


############# Body ###########

for(rxN in 1:length(rxnList_form)){
  
  rxnSummary <- rxnList_form[[rxN]]
  
  occupancyEq <- as.formula(paste("~", sub(paste(paste("E", rxnSummary$rxnID, sep = "_"), " \\* ", paste("V", rxnSummary$rxnID, sep = "_"), sep = ""), "1", rxnSummary$rxnForm)[2], sep = " "))
  
  ### Create a data.frame describing the relevent parameters for the model ###
  kineticPars <- data.frame(rel_spec = c(rxnSummary$enzymeAbund[,1], colnames(rxnSummary$rxnMet)), 
                            SpeciesType = c(rep("Enzyme", times = length(rxnSummary$enzymeAbund[,1])), 
                                            rep("Metabolite", times = length(colnames(rxnSummary$rxnMet)))), 
                            modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  kineticPars$formulaName[kineticPars$SpeciesType == "Enzyme"] <- paste("E", rxnSummary$rxnID, sep = "_")
  kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
  kineticPars$commonName[kineticPars$SpeciesType == "Metabolite"] <- rxnSummary$metNames[kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]]
  kineticPars$commonName[kineticPars$SpeciesType == "Enzyme"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  kineticPars$formulaName[kineticPars$SpeciesType == "Metabolite"] <- paste("K", rxnSummary$rxnID, kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"], sep = "_")
  
  all_species <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, occupancyEq)) != 0}) | kineticPars$SpeciesType == "Enzyme",]
  
  kineticPars <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, occupancyEq)) != 0}),] #remove species which dont appear in the reaction equation
  kineticPars <- rbind(kineticPars, c("keq", "keq", NA, NA, paste("Keq", rxnSummary$rxnID, sep = ""), NA))
  
  ### Create a matrix containing the metabolites and enzymes 
  enzyme_abund <- t(rxnSummary$enzymeAbund[,cond_mapping$enzyme_reordering]); colnames(enzyme_abund) <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  met_abund <- rxnSummary$rxnMet
  met_abund <- met_abund[,colnames(met_abund) %in% kineticPars$rel_spec]
  
  if(length(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]) <= 1){
    met_abund <- data.frame(met_abund)
    colnames(met_abund) <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
    kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- !all(is.na(met_abund))
  }else{
    kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){(apply(is.na(met_abund), 2, sum) == 0)[names((apply(is.na(met_abund), 2, sum) == 0)) == x]}))
  }
  
  
  ### set missing data to invariant across conditions
  met_abund[,!as.logical(kineticPars$measured[kineticPars$rel_spec %in% colnames(met_abund)])] <- 0
  met_abund <- 2^met_abund
  colnames(met_abund) <- unname(sapply(colnames(met_abund), function(x){kineticPars$modelName[kineticPars$rel_spec == x]}))
  
  enzyme_abund <- 2^enzyme_abund
  
  flux <- rxnSummary$flux/median(abs(rxnSummary$flux[rxnSummary$flux != 0])) #flux, scaled to a prettier range
  
  #### write flux parameters to a list ####
  run_summary[[names(rxnList_form)[rxN]]]$metabolites <- met_abund
  run_summary[[names(rxnList_form)[rxN]]]$enzymes <- enzyme_abund
  run_summary[[names(rxnList_form)[rxN]]]$flux <- flux
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq <- occupancyEq
  run_summary[[names(rxnList_form)[rxN]]]$all_species <- all_species
  run_summary[[names(rxnList_form)[rxN]]]$rxnSummary <- rxnSummary
  
  
  
  kineticParPrior <- data.frame(distribution = rep(NA, times = length(kineticPars[,1])), par_1 = NA, par_2 = NA) #par1/2 of a uniform are the lower bound and upper bound; par1/2 of a normal are the mean and variance
  kineticParPrior$distribution <- "unif"
  kineticParPrior$par_1 <- -10; kineticParPrior$par_2 <- 10
  for(exp_param in kineticPars$modelName[!is.na(kineticPars$measured) & kineticPars$measured == TRUE]){
    kineticParPrior[kineticPars$modelName == exp_param & !is.na(kineticPars$modelName), c(2:3)] <- median(log(met_abund[,colnames(met_abund) == exp_param])) + c(-10,10)
  }#priors for measured metabolites (either in absolute or relative space) are drawn about the median
  
  
  lik_track <- NULL
  markov_par_vals <- matrix(NA, ncol = length(kineticPars[,1]), nrow = markov_pars$n_samples)
  colnames(markov_par_vals) <- kineticPars$formulaName
  
  current_pars <- rep(NA, times = length(kineticParPrior[,1]))
  current_pars <- par_draw(1:length(kineticParPrior[,1]))
  current_lik <- lik_calc(current_pars)
  
  proposed_params <- current_pars
  
  for(i in 1:(markov_pars$burn_in + markov_pars$n_samples*markov_pars$sample_freq)){
    for(j in 1:length(kineticPars[,1])){#loop over parameters values
      proposed_par <- par_draw(j)
      proposed_lik <- lik_calc(proposed_par)
      if(runif(1, 0, 1) < exp(proposed_lik - current_lik)){
        current_pars <- proposed_par
        current_lik <- proposed_lik
      }
    }
    
    if(i > markov_pars$burn_in){
      if((i - markov_pars$burn_in) %% markov_pars$sample_freq == 0){
        markov_par_vals[(i - markov_pars$burn_in)/markov_pars$sample_freq,] <- current_pars
        lik_track <- c(lik_track, current_lik)
      }
    }
  }
  
  #colnames(markov_par_vals) <- ifelse(kineticPars$SpeciesType == "keq", "keq", kineticPars$commonName)
  colnames(markov_par_vals) <- kineticPars$rel_spec
  
  kineticPars$formatted[kineticPars$SpeciesType != "keq"] <- unname(sapply(kineticPars$commonName[kineticPars$SpeciesType != "keq"], function(name_int){
    if(length(strsplit(name_int, split = "")[[1]]) >= 25){
      split_name <- strsplit(name_int, split = "")[[1]]
      split_pois <- c(1:length(split_name))[split_name %in% c(" ", "-")][which.min(abs(20 - c(1:length(split_name)))[split_name %in% c(" ", "-")])]
      split_name[split_pois] <- "\n"
      paste(split_name, collapse = "")
    }else{name_int}
  }))
  kineticPars$formatted[kineticPars$SpeciesType == "keq"] <- "keq"
  
  run_summary[[names(rxnList_form)[rxN]]]$kineticPars <- kineticPars
  run_summary[[names(rxnList_form)[rxN]]]$markovChain <- markov_par_vals
  run_summary[[names(rxnList_form)[rxN]]]$likelihood <- lik_track
  
}


save(run_summary, file = paste(c("paramSets/paramSet", runNum, ".Rdata"), collapse = ""))

