# Works with R version >= 2.15

#qsub -l 1day -cwd -sync n Rscript FluxOptimClus.R runNum=$a_run chunk=$a_chunk

#setwd("/Genomics/grid/users/shackett/FBA/FluxParOptim/")
setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(nnls)

options(stringsAsFactors = FALSE)

args <- commandArgs()
runNum = as.numeric(unlist(strsplit(args[grep("runNum", args)], "="))[2])
chunkNum = as.numeric(unlist(strsplit(args[grep("chunk", args)], "="))[2])

print(paste("CHUNK NUMBER: ", chunkNum, ", RUN NUMBER: ", runNum, sep = ""))

rurun_summary <- list() #### MCMC run output and formatted inputs

#markov_pars <- list()
#markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
#markov_pars$n_samples <- 10000 #how many total markov samples are desired
#markov_pars$burn_in <- 500 #how many initial samples should be skipped

markov_pars <- list()
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 10 #how many total markov samples are desired
markov_pars$burn_in <- 0 #how many initial samples should be skipped


run_summary$markov_pars <- markov_pars

load("paramOptim.Rdata")
rxnList_form <- rxnList_form[names(rxnList_form) %in% chunk_assignment$set[chunk_assignment$chunk == chunkNum]]

if(chunkNum %in% chunk_assignment$chunk[grep('metX', chunk_assignment$set)]){
  # If there are reactions proposing hypothetical regulators, load principal components
  load('flux_cache/metaboliteTables.RData')
  npc <- ncol(metSVD$v)
  
  PC_loading_pars <- data.frame(mean = apply(metSVD$u, 2, mean), sd = apply(metSVD$u, 2, sd)) # metabolomic loadings of principal components
  
}



##@##@##@##@##@##@##@##@##@##@##@##@##@
############# Functions ###############
##@##@##@##@##@##@##@##@##@##@##@##@##@


par_draw <- function(updates){
  ### Update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ###
  ### Parameters are all returned in log-space (base e) ###
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    } else if(kineticParPrior$distribution[par_n] == "norm"){
      draw[par_n] <- rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
    } else if(kineticParPrior$distribution[par_n] == "SpSl"){
      draw[par_n] <- ifelse(rbinom(1, 1, kineticParPrior$par_3[par_n]) == 0, 0, rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n]))
    } else {
      print("invalid distribution")
      die
    }
  }
  draw
}



lik_calc <- function(proposed_params){
  ### Determine the likelihood of predicted flux as a function of metabolite abundance and kinetics parameters relative to actual flux ###
  
  par_stack <- rep(1, n_c) %*% t(proposed_params); colnames(par_stack) <- kineticPars$formulaName
  par_stack <- par_stack[,!(kineticPars$SpeciesType %in% "PCL")]
  
  par_stack <- 2^par_stack
  occupancy_vals <- data.frame(met_abund, par_stack)
  
  predOcc <- model.matrix(rxnEquations[["l_occupancyEq"]], data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
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
  
  logLik = log(lik)
  
  if(any(logLik == "-Inf")){
    
    log_cumsum <- data.frame(high_max = pnorm(flux$FVAmax[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = F),
                             high_min = pnorm(flux$FVAmin[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = F),
                             low_max = pnorm(flux$FVAmax[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = T),
                             low_min = pnorm(flux$FVAmax[lik == 0], flux_fit$fitted[lik == 0], flux_SD[lik == 0], log = T, lower.tail = T))
    
    RT <- ifelse(log_cumsum$high_max < log_cumsum$low_min, T, F) # is flux FVA >> par -> TRUE
    
    log_cumsum$high_min[RT][log_cumsum$high_max[RT] == log_cumsum$high_min[RT]] <- log_cumsum$high_min[RT][log_cumsum$high_max[RT] == log_cumsum$high_min[RT]] + 10^-4
    log_cumsum$low_max[!RT][log_cumsum$low_min[!RT] == log_cumsum$low_max[!RT]] <- log_cumsum$low_max[!RT][log_cumsum$low_min[!RT] == log_cumsum$low_max[!RT]] + 10^-4
    
    # For minute densities: calculate them in log space
    # To find log[p(x) / xmax - xmin] -> log[p(x)] - log(xmax - xmin)
    # -@-@ Find log[p(x)] -> log[p(max)] + log(1 - exp(log[p(max)] - log[p(min)]))
    # -@-@ subtract log[FVAmax - FVAmin]
    
    logLik[lik == 0] <- apply(data.frame(RT = log_cumsum$high_min + log(1 - exp(log_cumsum$high_max - log_cumsum$high_min)) - log(flux$FVAmax[lik == 0] - flux$FVAmin[lik == 0]), 
                                         LT = log_cumsum$low_max + log(1 - exp(log_cumsum$low_min - log_cumsum$low_max)) - log(flux$FVAmax[lik == 0] - flux$FVAmin[lik == 0])), 1, max)
    
  }
  
  
  if(any(logLik == "-Inf")){ # if there is a parmeter set that is so infeasible that it rounds to -Inf in log-space!
    -Inf
  }else if(any(flux_SD == 0)){ # if SD collapses to zero, because kcat = 0 (e.g. flux is always negative and prediction is always positive)
    -Inf
  }else{
    sum(logLik)
  }
  
}



metX_calc <- function(proposed_params, kineticPars, treatmentPCs){
  ### Using current values of principal component loadings and the fixed eigenvalues and principal components,
  ### calculate the log-abunance of the hypothetical metabolite
  
  c(proposed_params[kineticPars$SpeciesType == "PCL"] %*% diag(metSVD$d[1:npc]) %*% t(treatmentPCs))
  
}




##@##@##@##@##@##@##@##@##@##@##@##@##@
############### Body ##################
##@##@##@##@##@##@##@##@##@##@##@##@##@


#for(rxN in c(1:length(rxnList_form))[names(rxnList_form) %in% c("r_1088-rm", "r_0250-rm-t_0674-inh-comp")]){
for(rxN in 1:length(rxnList_form)){
  
  t_start = proc.time()[3]
  print(paste(names(rxnList_form)[rxN], "started", sep = " "))
  
  rxnSummary <- rxnList_form[[rxN]]
  n_c <- nrow(rxnSummary$flux)
  
  rxnEquations <- list()
  rxnEquations[["occupancyEq_list"]] <- paste(deparse(as.list(rxnSummary$rxnForm)[[2]]), collapse = "") # a parametric form relating metabolites and constants to fraction of maximal activity
  rxnEquations[["occupancyEq_list"]] <- gsub('[ ]+', ' ', rxnEquations[["occupancyEq_list"]])
  rxnEquations[["occupancyEq_list"]] <- gsub('\\^1234', '^h_allo', rxnEquations[["occupancyEq_list"]]) # 1234 was a standin for a hill coefficient that isn't pre-specified as a numbewr
  
  rxnEquations[["l_occupancyEq_list"]] <- rxnEquations[["occupancyEq_list"]]
  rxnEquations[["l_occupancyEq_list"]] <- gsub('([^_])(t_)', '\\12^\\2', rxnEquations[["l_occupancyEq_list"]])
  rxnEquations[["l_occupancyEq"]] <- as.formula(paste("~", rxnEquations[["l_occupancyEq_list"]], collapse = "")) # same equation using naturally using log2 data
  
  ### Create a data.frame describing the relevent parameters for the model ###
  kineticPars <- data.frame(rel_spec = c(rownames(rxnSummary$enzymeComplexes), colnames(rxnSummary$rxnMet)), 
                            SpeciesType = c(rep("Enzyme", times = nrow(rxnSummary$enzymeComplexes)), rep("Metabolite", times = ncol(rxnSummary$rxnMet))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  kineticPars$formulaName[kineticPars$SpeciesType == "Enzyme"] <- paste(paste("E", rxnSummary$rxnID, sep = "_"), kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"], sep = "_")
  kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
  kineticPars$commonName[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){rxnSummary$metNames[names(rxnSummary$metNames) == x]}))
  kineticPars$commonName[kineticPars$SpeciesType == "Enzyme"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  kineticPars$formulaName[kineticPars$SpeciesType == "Metabolite"] <- paste("K", rxnSummary$rxnID, kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"], sep = "_")
  
  all_species <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, rxnEquations[["occupancyEq_list"]])) != 0}) | kineticPars$SpeciesType == "Enzyme",]
  kineticPars <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, rxnEquations[["occupancyEq_list"]])) != 0}),] #remove species which dont appear in the reaction equation
  
  kineticPars <- rbind(kineticPars, c("keq", "keq", NA, NA, paste("Keq", rxnSummary$rxnID, sep = ""), NA))
  if(length(grep('\\^h', rxnEquations[["occupancyEq_list"]])) != 0){ # an unspecified hill coefficient was found
    kineticPars <- rbind(kineticPars, c("h_allo", "hillCoefficient", NA, NA, "h_allo", NA))
  }
  if("t_metX" %in% all_species$rel_spec){
    # if searching for a hypothetical allosteric modifier informed by the metabolomic principal components is desiered
    # principal component loadings are sampled in log-space from N(mean(loading), sd(loading)) - they are approximately log-normal
    # treating each conditions principal components (v) and the eigenvalues (d) as fixed, sample loadings one-by-one, treating each as a seperate parameter w.r.t. metropolis optimization
    
    treatmentPCs <- metSVD$v[rownames(metSVD$v) %in% rownames(rxnSummary$rxnMet),]
    
    kineticPars <- rbind(kineticPars, data.frame(rel_spec = paste("PCL", 1:npc, sep = "_"), SpeciesType = "PCL", modelName = NA, commonName = NA, formulaName = NA, measured = NA))
    
  }
  
  
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
  Kcatpaste <- paste('I(', KcatExpression, ' * ', sep = "")
  
  rxnEquations[["full_kinetic_form_list"]] <- rxnEquations[["l_occupancyEq_list"]]
  rxnEquations[["full_kinetic_form_list"]] <- gsub('(I\\()', Kcatpaste, rxnEquations[["full_kinetic_form_list"]])
  rxnEquations[["full_kinetic_form_list"]] <- sub('I\\(', '\\(', rxnEquations[["full_kinetic_form_list"]])
  
  # find the partial derivatives of the kinetic form for each reaction specie
  eq <- eval(parse(text = paste('expression(',rxnEquations[["full_kinetic_form_list"]],')')))
  
  D_full_kinetic_form <- list()
  for(spec in c(kineticPars$rel_spec[kineticPars$measured & !is.na(kineticPars$measured)], all_species$rel_spec[all_species$SpeciesType == "Enzyme"])){
    D_full_kinetic_form[[spec]] <- D(eq, spec)
  }
  rxnEquations[["kinetic_form_partials"]] <- D_full_kinetic_form
  
  
  flux <- rxnSummary$flux/median(abs(rxnSummary$flux$standardQP[rxnSummary$flux$standardQP != 0])) #flux, scaled to a prettier range
  
  # If FVA min flux is greater than max flux, switch them (and print a warning).
  
  if(sum(!(flux$FVAmax >= flux$FVAmin)) != 0){
    print("maximum flux is less than minimum flux")
  }
  flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmin', 'FVAmax')] <- flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmax', 'FVAmin')]
  
  # If bounds are exactly equal, then introduce a minute spread so a range can be calculated ###
  
  flux$FVAmax[flux$FVAmax == flux$FVAmin] <- flux$FVAmax[flux$FVAmax == flux$FVAmin] + flux$FVAmax[flux$FVAmax == flux$FVAmin]*10^-4
  
  # Metabolite SD and correlation specified so that they can be passed to v(mets, par)
  
  species_SD <- rxnSummary$all_species_SD
  species_corr <- rxnSummary$all_species_corr
  
  
  kineticParPrior <- data.frame(distribution = rep(NA, times = length(kineticPars[,1])), par_1 = NA, par_2 = NA, par_3 = NA) 
  
  # Options for these parameters are:
  # -@-@ unif: uniform in log-space: par_1 = lowerbound, par_2 = upperbound. draw in log2 space and exponentiate back to linear space
  # -@-@ norm: lognormal: in log2 space draw a value from mean = par_1, sd = par_2
  # -@-@ SpSl: spike and slab (In log2 space: the spike is a point mass at zero with p = par_3, the slab is a normal with mean = par_1 and sd = par_2)
  
  # Specify prior for michaelis constants
  kineticParPrior$distribution[kineticPars$SpeciesType %in% c("Metabolite", "keq")] <- "unif"
  kineticParPrior$par_1[kineticParPrior$distribution == "unif"] <- -10; kineticParPrior$par_2[kineticParPrior$distribution == "unif"] <- 10 # default value to 2^-10:2^10
  
  for(exp_param in kineticPars$modelName[!is.na(kineticPars$measured) & kineticPars$measured == TRUE]){
    kineticParPrior[kineticPars$modelName == exp_param & !is.na(kineticPars$modelName), c(2:3)] <- median(met_abund[,colnames(met_abund) == exp_param]) + c(-10,10)
  }#priors for measured metabolites (either in absolute or relative space) are drawn about the median
  
  # Prior for keq is centered around sum log(sub) - sum log(prod) - this deals with some species being in absolute space and others being absolute measurements
  kineticParPrior[kineticPars$SpeciesType == "keq", c(2:3)] <- median(rowSums(met_abund * c(rep(1,n_c)) %*% t(rxnSummary$rxnFormData$Stoi))) + c(-10,10)
  
  # Specify prior for hill constants
  kineticParPrior$distribution[kineticPars$SpeciesType == "hillCoefficient"] <- "SpSl"
  kineticParPrior$par_1[kineticPars$SpeciesType == "hillCoefficient"] <- 0
  kineticParPrior$par_2[kineticPars$SpeciesType == "hillCoefficient"] <- 0.5
  kineticParPrior$par_3[kineticPars$SpeciesType == "hillCoefficient"] <- 0.5
  
  if("t_metX" %in% all_species$rel_spec){
    # Specify prior for principal component loadings
    kineticParPrior$distribution[kineticPars$SpeciesType == "PCL"] <- "norm"
    kineticParPrior$par_1[kineticPars$SpeciesType == "PCL"] <- PC_loading_pars$mean
    kineticParPrior$par_2[kineticPars$SpeciesType == "PCL"] <- PC_loading_pars$sd
  }
  
  
  #### Optimize l(par|X) using Metropolis-Hastings MCMC ####
  
  
  ### Initialize parameters & setup tracking of likelihood and parameters ###
  
  
  lik_track <- NULL
  markov_par_vals <- matrix(NA, ncol = nrow(kineticParPrior), nrow = markov_pars$n_samples)
  colnames(markov_par_vals) <- kineticPars$formulaName
  
  current_pars <- rep(NA, times = nrow(kineticParPrior))
  current_pars <- par_draw(1:nrow(kineticParPrior))
  if("t_metX" %in% all_species$rel_spec){met_abund$t_metX <- metX_calc(current_pars, kineticPars, treatmentPCs)}
  
  current_lik <- lik_calc(current_pars)
  
  proposed_params <- current_pars
  
  ### Generate markov chain ###
  
  for(i in 1:(markov_pars$burn_in + markov_pars$n_samples*markov_pars$sample_freq)){
    for(j in 1:nrow(kineticPars)){#loop over parameters values
      proposed_par <- par_draw(j)
      if("t_metX" %in% all_species$rel_spec){met_abund$t_metX <- metX_calc(proposed_par, kineticPars, treatmentPCs)}
      proposed_lik <- lik_calc(proposed_par)
      if(runif(1, 0, 1) < exp(proposed_lik - current_lik) | (proposed_lik == current_lik & proposed_lik == -Inf)){
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
  
  ## Save MC information
  run_summary[[names(rxnList_form)[rxN]]]$kineticPars <- kineticPars
  run_summary[[names(rxnList_form)[rxN]]]$markovChain <- markov_par_vals
  run_summary[[names(rxnList_form)[rxN]]]$likelihood <- lik_track
  ##  Save flux and rxn information
  run_summary[[names(rxnList_form)[rxN]]]$metabolites <- met_abund
  run_summary[[names(rxnList_form)[rxN]]]$enzymes <- enzyme_abund
  run_summary[[names(rxnList_form)[rxN]]]$flux <- flux
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq <- rxnEquations
  run_summary[[names(rxnList_form)[rxN]]]$all_species <- all_species
  run_summary[[names(rxnList_form)[rxN]]]$rxnSummary <- rxnSummary
  run_summary[[names(rxnList_form)[rxN]]]$specSD <- species_SD
  run_summary[[names(rxnList_form)[rxN]]]$specCorr <- species_corr
  ## Save priors
  run_summary[[names(rxnList_form)[rxN]]]$kineticParPrior <- data.frame(kineticPars, kineticParPrior)
  
  print(paste(names(rxnList_form)[rxN], "finished in ", round(proc.time()[3] - t_start, 0), " seconds", sep = " "))
  
  
}

save(run_summary, file = paste(c("paramSets/paramSet", "C", chunkNum, "R", runNum, ".Rdata"), collapse = ""))
