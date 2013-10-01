#qsub -l 1day -cwd -sync n Rscript FluxOptimClus.R runNum=$a_run chunk=$a_chunk

#setwd("/Genomics/grid/users/shackett/FBA/FluxParOptim/")
setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(nnls)

options(stringsAsFactors = FALSE)

args <- commandArgs()
runNum = as.numeric(unlist(strsplit(args[grep("runNum", args)], "="))[2])
chunkNum = as.numeric(unlist(strsplit(args[grep("chunk", args)], "="))[2])


run_summary <- list() #### MCMC run output and formatted inputs

#markov_pars <- list()
#markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
#markov_pars$n_samples <- 2000 #how many total markov samples are desired
#markov_pars$burn_in <- 500 #how many initial samples should be skipped

markov_pars <- list()
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 100 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped


run_summary$markov_pars <- markov_pars

load("paramOptim.Rdata")
rxnList_form <- rxnList_form[names(rxnList_form) %in% chunk_assignment$set[chunk_assignment$chunk == chunkNum]]

if(chunkNum %in% chunk_assignment$chunk[grep('metX', chunk_assignment$set)]){
  # If there are reactions proposing hypothetical regulators, load principal components
  load('flux_cache/metaboliteTables.RData')
  
  apply(metSVD$u, 2, mean)
  apply(metSVD$u, 2, sd)
  
  }

  
#cherryPicked <- c("r_0232-rm", "r_0277-rm", "r_0789-rm", "r_0484-rm", "r_0859-rm", "r_0941-rm",
# "r_0232-rm-t_0717-inh-noncomp", "r_0232-rm-t_0287-inh-noncomp", "r_0232-rm-t_0582-act-allo",  
# "r_0789-rm-t_0248-act-allo", "r_0859-rm-t_0296-act-allo", "r_0859-rm-t_0248-act-allo",  
# "r_0859-rm-t_0604-act-allo", "r_0859-rm-t_0231-act-allo", "r_0859-rm-t_0254-inh-noncomp",
# "r_0859-rm-t_0283-inh-noncomp", "r_0941-rm-t_0296-act-allo", "r_0941-rm-t_0283-inh-noncomp",
# "r_0859-rm-t_0446-inh-noncomp", "r_0859-rm-t_0446-inh-comp", "r_0859-rm-t_0446-inh-uncomp")

  
  

########### Functions ##########

par_draw <- function(updates){
  #### update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ####
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    } else if(kineticParPrior$distribution[par_n] == "norm"){
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
  
  predOcc <- model.matrix(l_occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*2^enzyme_abund #occupany of enzymes * relative abundance of enzymes
  
  # fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  # flux objective is set as the average of the minimal and maximal allowable flux flowing through the reaction at the optimal solution
  
  flux_fit <- nnls(enzyme_activity, (flux$FVAmax + flux$FVAmin)/2) 
  
  # calculate SD of fitted flux
  # propagate variance from enzyme + metabolite measurements into predicted flux using partial derivatives of reaction form
  # A) partial derivatives by condition
  # -@-@ form of partial derivatives: D_full_kinetic_form
  # -@-@ metabolite abundance and current parameter value
  # A) covariance matrix by condition - residual correlation * SDa * SDb
  # -@-@ residual correlation calculated across all conditions
  # -@-@ standard deviation calculated by condition
  
  nnlsCoef <- t(t(rep(1, n_c)))  %*% flux_fit$x; colnames(nnlsCoef) <- all_species$formulaName[all_species$SpeciesType == "Enzyme"]
  
  all_components <- data.frame(occupancy_vals, enzyme_abund, nnlsCoef)

  # partial derivatives of each measured specie in a condition
  comp_partials <- matrix(NA, nrow = n_c, ncol = length(D_full_kinetic_form))
  colnames(comp_partials) <- names(D_full_kinetic_form)
  
  for(j in 1:ncol(comp_partials)){
    comp_partials[,j] <- with(all_components ,eval(D_full_kinetic_form[[j]]))
  }
  
  # calculate the fitted standard deviation after first finding the by-condition residual covariance matrix
  
  flux_SD <- rep(NA, n_c)
  for(i in 1:n_c){
    sampleCov <- species_corr * t(t(species_SD[i,])) %*% species_SD[i,]
    flux_SD[i] <- sqrt(t(comp_partials[i,]) %*% sampleCov %*% t(t(comp_partials[i,])))
    }
  
  # evaluate the relative density a gaussian centered about fitted flux with the SD calculated above
  
  # if p(x = Xmax) - p(x = Xmin) != 0 (for numerical reasons)
  # Di = p(x - Xmax) - p(x = Xmin) / Xmax - Xmin
  
  lik <- (pnorm(flux$FVAmax, flux_fit$fitted, flux_SD) - pnorm(flux$FVAmin, flux_fit$fitted, flux_SD))/(flux$FVAmax - flux$FVAmin)
  
  log_cumsum <- data.frame(high_max = pnorm(flux$FVAmax[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = F),
    high_min = pnorm(flux$FVAmin[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = F),
    low_max = pnorm(flux$FVAmax[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = T),
    low_min = pnorm(flux$FVAmax[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = T))
  
  # For minute densities: calculate them in log space
  # To find log[p(x) / xmax - xmin] -> log[p(x)] - log(xmax - xmin)
  # -@-@ Find log[p(x)] -> log[p(max)] + log(1 - exp(log[p(max)] - log[p(min)]))
  # -@-@ subtract log[FVAmax - FVAmin]
  
  logLik = log(lik)
  
  logLik[lik == 0] <- apply(data.frame(RT = log_cumsum$high_min + log(1 - exp(log_cumsum$high_max - log_cumsum$high_min)) - log(flux$FVAmax[lik == 0] - flux$FVAmin[lik == 0]), 
             LT = log_cumsum$low_max + log(1 - exp(log_cumsum$low_min - log_cumsum$low_max)) - log(flux$FVAmax[lik == 0] - flux$FVAmin[lik == 0])), 1, max)
  
  if(any(logLik == "-Inf")){
    die
    }else{
     sum(logLik)
      }
 
  }


+############# Body ###########

for(rxN in 1:length(rxnList_form)){
  
  rxnSummary <- rxnList_form[[rxN]]
  n_c <- nrow(rxnSummary$flux)
  occupancyEq <- rxnSummary$rxnForm # a parametric form relating metabolites and constants to fraction of maximal activity
  l_occupancyEq <- as.formula(paste(gsub('([^_])(t_)', '\\12^\\2', occupancyEq), collapse = "")) # same equation using naturally using log2 data
  
  
  ### Create a data.frame describing the relevent parameters for the model ###
  kineticPars <- data.frame(rel_spec = c(rownames(rxnSummary$enzymeComplexes), colnames(rxnSummary$rxnMet)), 
  SpeciesType = c(rep("Enzyme", times = nrow(rxnSummary$enzymeComplexes)), rep("Metabolite", times = ncol(rxnSummary$rxnMet))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  kineticPars$formulaName[kineticPars$SpeciesType == "Enzyme"] <- paste(paste("E", rxnSummary$rxnID, sep = "_"), kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"], sep = "_")
  kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
  kineticPars$commonName[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){rxnSummary$metNames[names(rxnSummary$metNames) == x]}))
  kineticPars$commonName[kineticPars$SpeciesType == "Enzyme"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  kineticPars$formulaName[kineticPars$SpeciesType == "Metabolite"] <- paste("K", rxnSummary$rxnID, kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"], sep = "_")
  
  all_species <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, occupancyEq)) != 0}) | kineticPars$SpeciesType == "Enzyme",]
  kineticPars <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, occupancyEq)) != 0}),] #remove species which dont appear in the reaction equation
  
  kineticPars <- rbind(kineticPars, c("keq", "keq", NA, NA, paste("Keq", rxnSummary$rxnID, sep = ""), NA))
  
  ### Create a matrix containing the metabolites and enzymes 
  enzyme_abund <- t(rxnSummary$enzymeComplexes); colnames(enzyme_abund) <- all_species$rel_spec[all_species$SpeciesType == "Enzyme"]
  met_abund <- rxnSummary$rxnMet
  met_abund <- met_abund[,colnames(met_abund) %in% kineticPars$rel_spec]
  
  if(length(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]) <= 1){
    met_abund <- data.frame(met_abund)
    colnames(met_abund) <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
    kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- !all(is.na(met_abund))
    }else{
      kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){(apply(is.na(met_abund), 2, sum) == 0)[names((apply(is.na(met_abund), 2, sum) == 0)) == x]}))
      }
  kineticPars$measured <- as.logical(kineticPars$measured)
  
  met_abund[,!kineticPars$measured[kineticPars$rel_spec %in% colnames(met_abund)]] <- 0 # set missing data (unascertained) to invariant across conditions
  
  # expression combining the log-occupancy equation and scaling of enzyme abundance by activity
  
  KcatEs <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = sapply(all_species$rel_spec[all_species$SpeciesType == "Enzyme"], function(x){paste("2^", x, sep = "")}), Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
  KcatExpression <- paste('(', paste(KcatEs, collapse = " + "), ')', sep = "")
  Kcatpaste <- paste('I(', KcatExpression, '*')
  
  full_kinetic_form <- as.formula(gsub('(I\\()', Kcatpaste, l_occupancyEq))
  # find the partial derivatives of the kinetic form for each reaction specie
  eq <- eval(parse(text = paste('expression(',gsub('I\\(', '\\(', as.character(full_kinetic_form)[2]),')')))
  
  D_full_kinetic_form <- list()
  for(spec in c(kineticPars$rel_spec[kineticPars$measured & !is.na(kineticPars$measured)], all_species$rel_spec[all_species$SpeciesType == "Enzyme"])){
    D_full_kinetic_form[[spec]] <- D(eq, spec)
    }
      
  
  #colnames(met_abund) <- unname(sapply(colnames(met_abund), function(x){kineticPars$modelName[kineticPars$rel_spec == x]}))
  #met_abund <- 2^met_abund
  #enzyme_abund <- 2^enzyme_abund
  
  flux <- rxnSummary$flux/median(abs(rxnSummary$flux$standardQP[rxnSummary$flux$standardQP != 0])) #flux, scaled to a prettier range
  
  #If FVA min flux is greater than max flux, switch them (and print a warning).
  
  if(sum(!(flux$FVAmax >= flux$FVAmin)) != 0){
    print("maximum flux is less than minimum flux")
    }
  
  flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmin', 'FVAmax')] <- flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmax', 'FVAmin')]
  
  #If bounds are exactly equal, then introduce a minute spread so a range can be calculated ###
  
  flux$FVAmax[flux$FVAmax == flux$FVAmin] <- flux$FVAmax[flux$FVAmax == flux$FVAmin] + flux$FVAmax[flux$FVAmax == flux$FVAmin]*10^-4
  
  
  
  species_SD <- rxnSummary$all_species_SD
  species_corr <- rxnSummary$all_species_corr
  
  #### write flux parameters to a list ####
  run_summary[[names(rxnList_form)[rxN]]]$metabolites <- met_abund
  run_summary[[names(rxnList_form)[rxN]]]$enzymes <- enzyme_abund
  run_summary[[names(rxnList_form)[rxN]]]$flux <- flux
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq$linear <- occupancyEq
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq$log <- l_occupancyEq
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq$full <- full_kinetic_form
  run_summary[[names(rxnList_form)[rxN]]]$all_species <- all_species
  run_summary[[names(rxnList_form)[rxN]]]$rxnSummary <- rxnSummary
  run_summary[[names(rxnList_form)[rxN]]]$specSD <- species_SD
  run_summary[[names(rxnList_form)[rxN]]]$specCorr <- species_corr
  
  
  
  kineticParPrior <- data.frame(distribution = rep(NA, times = length(kineticPars[,1])), par_1 = NA, par_2 = NA) #par1/2 of a uniform are the lower bound and upper bound; par1/2 of a normal are the mean and variance
  kineticParPrior$distribution <- "unif"
  kineticParPrior$par_1 <- -10; kineticParPrior$par_2 <- 10
  for(exp_param in kineticPars$modelName[!is.na(kineticPars$measured) & kineticPars$measured == TRUE]){
    kineticParPrior[kineticPars$modelName == exp_param & !is.na(kineticPars$modelName), c(2:3)] <- median(met_abund[,colnames(met_abund) == exp_param]) + c(-10,10)
    }#priors for measured metabolites (either in absolute or relative space) are drawn about the median
  run_summary[[names(rxnList_form)[rxN]]]$kineticParPrior <- data.frame(kineticPars, kineticParPrior)
  
  
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
  
save(run_summary, file = paste(c("paramSets/paramSet", "C", chunkNum, "R", runNum, ".Rdata"), collapse = ""))
