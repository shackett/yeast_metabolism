# Works with R version >= 2.15

#qsub -l 1day -cwd -sync n Rscript FluxOptimClus.R runNum=$a_run chunk=$a_chunk

library(nnls)

options(stringsAsFactors = FALSE)


run_location <- "local" # either local or cluster
if(run_location == "local"){
  
  setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")
  
}else if(run_location == "cluster"){
  
  setwd("/Genomics/grid/users/shackett/FBA/FluxParOptim/")
  
  .libPaths(append(.libPaths(), "/Genomics/grid/users/shackett/R/x86_64-unknown-linux-gnu-library/2.15"))
  
  args <- commandArgs()
  runNum = as.numeric(unlist(strsplit(args[grep("runNum", args)], "="))[2])
  chunkNum = as.numeric(unlist(strsplit(args[grep("chunk", args)], "="))[2])
  
  print(paste("CHUNK NUMBER: ", chunkNum, ", RUN NUMBER: ", runNum, sep = ""))
  
}else{
  stop("provide a valid value to run_location") 
}

run_summary <- list() #### MCMC run output and formatted inputs
run_summary$markov_pars <- markov_pars


load("paramOptim.Rdata")

if(run_location == "cluster"){
  rxnList_form <- rxnList_form[names(rxnList_form) %in% chunk_assignment$set[chunk_assignment$chunk == chunkNum]]
}

if(chunkNum %in% chunk_assignment$chunk[grep('metX', chunk_assignment$set)] | run_location == "local"){
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
      stop("invalid distribution")
    }
  }
  draw
}


lik_calc_fittedSD <- function(proposed_params){

  #### determine the likelihood of predicted flux as a function of metabolite abundance and kinetics parameters relative to actual flux ####
  
  par_stack <- rep(1, n_c) %*% t(proposed_params); colnames(par_stack) <- kineticPars$formulaName
  par_stack <- par_stack[,!(kineticPars$SpeciesType %in% "PCL")]
  
  par_stack <- 2^par_stack
  occupancy_vals <- data.frame(met_abund, par_stack)
  
  if(!(kinetically_differing_isoenzymes)){
    predOcc <- eval(rxnEquations[["l_occupancyExpression"]], occupancy_vals) #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
    enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*2^enzyme_abund #occupany of enzymes * relative abundance of enzymes
  }else{
    enzyme_activity <- NULL
    for(isoenzyme in names(rxnSummary$rxnForm)){
      predOcc <- eval(rxnEquations[["l_occupancyExpression"]][[isoenzyme]], occupancy_vals)
      enzyme_activity <- cbind(enzyme_activity, predOcc %*% t(rep(1, sum(occEqtn_complex_match$occEqtn == isoenzyme)))*2^enzyme_abund[,colnames(enzyme_abund) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]])
    }
  }
  
  # fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  # flux objective is set as the average of the minimal and maximal allowable flux flowing through the reaction at the optimal solution
  
  flux_fit <- nnls(enzyme_activity, (flux$FVAmax + flux$FVAmin)/2) 
  
  # setting the variance from residual mean-squared error
  
  fit_resid_error <- sqrt(mean((flux_fit$resid - mean(flux_fit$resid))^2))*sqrt(n_c/(n_c-1))
  
  lik <- (pnorm(flux$FVAmax, flux_fit$fitted, fit_resid_error) - pnorm(flux$FVAmin, flux_fit$fitted, fit_resid_error))/(flux$FVAmax - flux$FVAmin)
  
  return(sum(log(lik)))
  
  }
  
  
lik_calc_propagatedSD <- function(proposed_params){
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


shmatch <- function(x,y){
  ### A quick (less efficient) replacement for data.table::chmatch to avoid its dependencies
  sapply(x, function(z){which(y == z)[1]})
  }


##@##@##@##@##@##@##@##@##@##@##@##@##@
############### Body ##################
##@##@##@##@##@##@##@##@##@##@##@##@##@

# test cases

rxTests <- c(which(names(rxnList_form) == "r_0005-rm-t_metX-inh-uncomp_ultra"), # ultra-sensitive allostery
             which(names(rxnList_form) == "r_0042_Y_F_inhibition_isoenzymeSpecific"), # isoenzyme specific regulation
              which(names(rxnList_form) == "r_0042_E4P_enzyme_specific_affinity_test2"), # isoenzyme specific kinetics w.r.t. substrate
             which(names(rxnList_form) == "r_0148-rm-t_0234-inh-uncomp"),# auto-regulation
              which(names(rxnList_form) == "r_0962-rm-t_0468-act-mm_ultra"),
             which(names(rxnList_form) == "r_0250-rm-t_0219-inh-uncomp"),
             which(names(rxnList_form) == "r_0211-rm-t_0234-inh-noncomp"),
             which(names(rxnList_form) == "r_0250-cc")) 

library(dplyr)

NLS_rxns <- data.frame(rxnID = names(rxnList_form), row = 1:length(rxnList_form)) %>% tbl_df() %>%
  filter(!grepl('allMetabolites$', rxnID)) %>%
  filter(!grepl('metX', rxnID)) %>%
  filter(!grepl('forward', rxnID))

NLS_fit_list <- list()

#rxN <- 21
#rxN <- 20
for(rxN in NLS_rxns$row){
#for(rxN in rxTests){  
  
  rxnSummary <- rxnList_form[[rxN]]
  n_c <- nrow(rxnSummary$flux)
  
  ##### Generate non-linear part(s) of kinetic form ####
  
  # if isoenzymes differ w.r.t. kinetics or regulation then their occupancy equation are stored as different elements of a list
  # and enzyme concentrations need to be paired with these seperate equations
  
  kinetically_differing_isoenzymes <- any(names(rxnSummary$rxnForm) %in% rownames(rxnSummary$enzymeComplexes))
  
  rxnEquations <- list()
  
  if(kinetically_differing_isoenzymes){
    for(isoenzyme in names(rxnSummary$rxnForm)){
      
      rxnEquations$occupancyEq_list[[isoenzyme]] <- paste(deparse(as.list(rxnSummary$rxnForm[names(rxnSummary$rxnForm) == isoenzyme])[[1]][[2]]), collapse = "") # a parametric form relating metabolites and constants to fraction of maximal activity
      rxnEquations$occupancyEq_list[[isoenzyme]] <- gsub('[ ]+', ' ', rxnEquations$occupancyEq_list[[isoenzyme]])
      
      rxnEquations$l_occupancyEq_list[[isoenzyme]] <- rxnEquations$occupancyEq_list[[isoenzyme]]
      rxnEquations$l_occupancyEq_list[[isoenzyme]] <- gsub('([^_])(t_)', '\\12^\\2', rxnEquations$l_occupancyEq_list[[isoenzyme]])
      rxnEquations$l_occupancyEq_list[[isoenzyme]] <- gsub('([0-9]+)\\^(t_[0-9]{4})\\^([0-9]+)', '(\\1\\^\\2)\\^\\3', rxnEquations$l_occupancyEq_list[[isoenzyme]]) # correct 2^X^2 -> (2^X)^2
      
      # create an expression for each isoenzyme
      rxnEquations[["l_occupancyExpression"]][[isoenzyme]] <- parse(text = sub(' \\+ 0$', '', sub('^I', '', rxnEquations$l_occupancyEq_list[[isoenzyme]]))) 
    }
  }else{
    # the typical case - single enzyme or kinetically equivalent isoenzymes
    rxnEquations[["occupancyEq_list"]] <- paste(deparse(as.list(rxnSummary$rxnForm)[[2]]), collapse = "") # a parametric form relating metabolites and constants to fraction of maximal activity
    rxnEquations[["occupancyEq_list"]] <- gsub('[ ]+', ' ', rxnEquations[["occupancyEq_list"]])
    
    rxnEquations[["l_occupancyEq_list"]] <- rxnEquations[["occupancyEq_list"]]
    rxnEquations[["l_occupancyEq_list"]] <- gsub('([^_])(t_)', '\\12^\\2', rxnEquations[["l_occupancyEq_list"]])
    rxnEquations[["l_occupancyEq_list"]] <- gsub('([0-9]+)\\^(t_[0-9]{4})\\^([0-9]+)', '(\\1\\^\\2)\\^\\3', rxnEquations[["l_occupancyEq_list"]]) # correct 2^X^2 -> (2^X)^2
    
    # create an expression for each isoenzyme
    rxnEquations[["l_occupancyExpression"]] <- parse(text = sub(' \\+ 0$', '', sub('^I', '', rxnEquations$l_occupancyEq_list))) 
  }
  
 
  #### Describe the relevent kinetic parameters ####
  
  if(kinetically_differing_isoenzymes & length(grep('.[0-9]$', names(rxnSummary$rxnForm))) != 0){
    # A single enzyme is applied to multiple expressions (ex. partial constitutive activity)
    kineticPars <- data.frame(rel_spec = c(names(rxnSummary$rxnForm), rxnSummary$rxnFormData$SubstrateID), 
                              SpeciesType = c(rep("Enzyme", times = length(rxnSummary$rxnForm)), rep("Metabolite", times = nrow(rxnSummary$rxnFormData))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  }else{
    # One enzyme, one Kcat
    kineticPars <- data.frame(rel_spec = c(rownames(rxnSummary$enzymeComplexes), rxnSummary$rxnFormData$SubstrateID), 
                            SpeciesType = c(rep("Enzyme", times = nrow(rxnSummary$enzymeComplexes)), rep("Metabolite", times = nrow(rxnSummary$rxnFormData))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  }
  kineticPars$formulaName[kineticPars$SpeciesType == "Enzyme"] <- paste(paste("E", rxnSummary$rxnID, sep = "_"), kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"], sep = "_")
  kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
  kineticPars$commonName[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){rxnSummary$metNames[names(rxnSummary$metNames) == x]}))
  kineticPars$commonName[kineticPars$SpeciesType == "Enzyme"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  kineticPars$formulaName[kineticPars$SpeciesType == "Metabolite"] <-  rxnSummary$rxnFormData$affinityParameter
  
  all_species <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, rxnEquations[["occupancyEq_list"]])) != 0}) | kineticPars$SpeciesType == "Enzyme",]
  kineticPars <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, rxnEquations[["occupancyEq_list"]])) != 0}),] #remove species which dont appear in the reaction equation
  
  kineticPars <- rbind(kineticPars, c("keq", "keq", NA, NA, paste("Keq", rxnSummary$rxnID, sep = ""), NA))
  if(any(rxnSummary$rxnFormData$Hill == 0)){ # an unspecified hill coefficient was found - introduce hill parameter
    kineticPars <- rbind(kineticPars, c(rxnSummary$rxnFormData$SubstrateID[rxnSummary$rxnFormData$Hill == 0], "hillCoefficient", NA, NA, sub('^K', 'KH', rxnSummary$rxnFormData$affinityParameter[rxnSummary$rxnFormData$Hill == 0]), NA))
  }
  
  if("t_metX" %in% all_species$rel_spec){
    # if searching for a hypothetical allosteric modifier informed by the metabolomic principal components is desiered
    # principal component loadings are sampled in log-space from N(mean(loading), sd(loading)) - they are approximately log-normal
    # treating each conditions principal components (v) and the eigenvalues (d) as fixed, sample loadings one-by-one, treating each as a seperate parameter w.r.t. metropolis optimization
    
    treatmentPCs <- metSVD$v[rownames(metSVD$v) %in% rownames(rxnSummary$rxnMet),]
    
    kineticPars <- rbind(kineticPars, data.frame(rel_spec = paste("PCL", 1:npc, sep = "_"), SpeciesType = "PCL", modelName = NA, commonName = NA, formulaName = NA, measured = NA))
    
  }
  
  #### Setup metabolite concentration, enzyme abundance and flux carried ####
  
  ### Create a matrix containing the metabolites and enzymes
  if(kinetically_differing_isoenzymes & length(grep('.[0-9]$', names(rxnSummary$rxnForm))) != 0){
    # A single enzyme is applied to multiple expressions (ex. partial constitutive activity)
    enzyme_abund <- t(rxnSummary$enzymeComplexes)[,shmatch(gsub('.[0-9]$', '', names(rxnSummary$rxnForm)), rownames(rxnSummary$enzymeComplexes))]
    colnames(enzyme_abund) <- names(rxnSummary$rxnForm)
    
    species_SD <- rxnSummary$all_species_SD
    remapping_indices <- shmatch(gsub('.[0-9]$', '', names(rxnSummary$rxnForm)), colnames(species_SD)[colnames(species_SD) %in% rownames(rxnSummary$enzymeComplexes)])
    remapping_table <- data.frame(remap = names(rxnSummary$rxnForm), enzyme = names(remapping_indices))
    renamed_enzyme_SD <- species_SD[,colnames(species_SD) %in% rownames(rxnSummary$enzymeComplexes)][,remapping_indices]
    colnames(renamed_enzyme_SD) <- names(rxnSummary$rxnForm)
    species_SD <- cbind(species_SD[,!(colnames(species_SD) %in% rownames(rxnSummary$enzymeComplexes))], renamed_enzyme_SD)
    
    species_corr <- matrix(0, nrow = ncol(species_SD), ncol = ncol(species_SD))
    rownames(species_corr) <- colnames(species_corr) <- colnames(species_SD)
    species_corr[!(colnames(species_SD) %in% names(rxnSummary$rxnForm)),!(colnames(species_SD) %in% names(rxnSummary$rxnForm))] <- rxnSummary$all_species_corr[!(rownames(rxnSummary$all_species_corr) %in% rownames(rxnSummary$enzymeComplexes)),!(rownames(rxnSummary$all_species_corr) %in% rownames(rxnSummary$enzymeComplexes))]
    
    remapped_enzyme_corr <- species_corr[(colnames(species_SD) %in% names(rxnSummary$rxnForm)),(colnames(species_SD) %in% names(rxnSummary$rxnForm))]
    for(i in 1:nrow(remapped_enzyme_corr)){
      for(j in 1:ncol(remapped_enzyme_corr)){
        remapped_enzyme_corr[i,j] <- ifelse(remapping_table$enzyme[remapping_table$remap == colnames(remapped_enzyme_corr)[j]] == remapping_table$enzyme[remapping_table$remap == rownames(remapped_enzyme_corr)[i]], 1, 0)
      }
    }
    
    species_corr[rownames(species_corr) %in% rownames(remapped_enzyme_corr), colnames(species_corr) %in% colnames(remapped_enzyme_corr)] <- remapped_enzyme_corr    
    
  }else{
    enzyme_abund <- t(rxnSummary$enzymeComplexes); colnames(enzyme_abund) <- all_species$rel_spec[all_species$SpeciesType == "Enzyme"]
    species_SD <- rxnSummary$all_species_SD
    species_corr <- rxnSummary$all_species_corr
    
  }
  
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
  
  met_abund[,!kineticPars$measured[shmatch(colnames(met_abund), kineticPars$modelName)]] <- 0 # set missing data (unascertained) to invariant across conditions
  
  flux <- rxnSummary$flux/median(abs(rxnSummary$flux$standardQP[rxnSummary$flux$standardQP != 0])) #flux, scaled to a prettier range
  
  # If FVA min flux is greater than max flux, switch them (and print a warning).
  
  if(sum(!(flux$FVAmax >= flux$FVAmin)) != 0){
    print(paste("maximum flux is less than minimum flux for", sum(!(flux$FVAmax >= flux$FVAmin)), "conditions"))
  }
  flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmin', 'FVAmax')] <- flux[!(flux$FVAmax >= flux$FVAmin),c('FVAmax', 'FVAmin')]
  
  # If bounds are exactly equal, then introduce a minute spread so a range can be calculated ###
  
  flux$FVAmax[flux$FVAmax == flux$FVAmin] <- flux$FVAmax[flux$FVAmax == flux$FVAmin] + abs(flux$FVAmax[flux$FVAmax == flux$FVAmin]*10^-4)
  
  #### Combine enzyme(s) with non-linear portion of the kinetic form to generate final kinetic form ####
  
  # expression combining the log-occupancy equation and scaling of enzyme abundance by activity
  # determined w.r.t. metabolites/enzymes in log-space and linear-space
  # if there are multiple isoenzymes which differ kinetically then their expressions are generated seperately and then concatenated
  # a single equation is returned
  
  if(kinetically_differing_isoenzymes){
    # if isoenzymes differ then their occupancy equations are not distributed e.g. Vmax1[O1] + Vmax2[O2] rather than (Vmax1 + Vmax2)*O
    KcatEs_log <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = sapply(all_species$rel_spec[all_species$SpeciesType == "Enzyme"], function(x){paste("2^", x, sep = "")}), Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    KcatEs_linear <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = all_species$rel_spec[all_species$SpeciesType == "Enzyme"], Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    
    if(!all(names(rxnEquations$occupancyEq_list) %in% c(names(KcatEs_log), "other"))){
      stop(paste("isoenzyme name does not match available complexes for reaction", names(rxnList_form)[rxN],"\nCheck the field \"enzymeInvolved\" in manual_ComplexRegulation.tsv"))
      }
    
    # a complex is matched to an enzyme either by being directly specified or otherwise being placed with "other" enzymes 
    occEqtn_complex_match <- data.frame(complex = names(KcatEs_log), occEqtn = NA)
    occEqtn_complex_match$occEqtn[occEqtn_complex_match$complex %in% names(rxnEquations$occupancyEq_list)] <- occEqtn_complex_match$complex[occEqtn_complex_match$complex %in% names(rxnEquations$occupancyEq_list)]
    occEqtn_complex_match$occEqtn[is.na(occEqtn_complex_match$occEqtn)] <- "other"
    
    KcatExpressions_linear <- NULL
    KcatExpressions_log <- NULL
    for(isoenzyme in unique(occEqtn_complex_match$occEqtn)){
      # log
      KcatExpression <- paste('(', paste(KcatEs_log[names(KcatEs_log) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]], collapse = " + "), ')', sep = "")
      isoenzymeO <- rxnEquations$l_occupancyEq_list[names(rxnEquations$occupancyEq_list) == isoenzyme]
      isoenzymeO <- sub('^I', paste(KcatExpression, '*', sep = ''), isoenzymeO)
      isoenzymeO <- sub(' \\+ 0$', '', isoenzymeO)
      
      KcatExpressions_log <- c(KcatExpressions_log, isoenzymeO)
      
      # linear
      KcatExpression <- paste('(', paste(KcatEs_linear[names(KcatEs_linear) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]], collapse = " + "), ')', sep = "")
      isoenzymeO <- rxnEquations$occupancyEq_list[names(rxnEquations$occupancyEq_list) == isoenzyme]
      isoenzymeO <- sub('^I', paste(KcatExpression, '*', sep = ''), isoenzymeO)
      isoenzymeO <- sub(' \\+ 0$', '', isoenzymeO)
      
      KcatExpressions_linear <- c(KcatExpressions_linear, isoenzymeO)
    }
    
    rxnEquations[["full_kinetic_form_list"]] <- paste(unlist(KcatExpressions_log), collapse = " + ")
    rxnEquations[["elasticity_calc"]] <- paste(unlist(KcatExpressions_linear), collapse = " + ")
    
  }else{
    
    # (Vmax1 + Vmax2)*O
    
    KcatEs <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = sapply(all_species$rel_spec[all_species$SpeciesType == "Enzyme"], function(x){paste("2^", x, sep = "")}), Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    KcatExpression <- paste('(', paste(KcatEs, collapse = " + "), ')', sep = "")
    Kcatpaste <- paste('I(', KcatExpression, ' * ', sep = "")
    
    rxnEquations[["full_kinetic_form_list"]] <- rxnEquations[["l_occupancyEq_list"]]
    rxnEquations[["full_kinetic_form_list"]] <- gsub('(I\\()', Kcatpaste, rxnEquations[["full_kinetic_form_list"]])
    rxnEquations[["full_kinetic_form_list"]] <- sub('I\\(', '\\(', rxnEquations[["full_kinetic_form_list"]])
    
    # save the same expression in linear space, so that it can be used later on to look at elasticitity
    
    KcatEs <- mapply(function(E, Kcat){paste(Kcat, E, sep = " * ")}, E = all_species$rel_spec[all_species$SpeciesType == "Enzyme"], Kcat = all_species$formulaName[all_species$SpeciesType == "Enzyme"])
    KcatExpression <- paste('(', paste(KcatEs, collapse = " + "), ')', sep = "")
    Kcatpaste <- paste('I(', KcatExpression, ' * ', sep = "")
    
    rxnEquations[["elasticity_calc"]] <- rxnEquations[["occupancyEq_list"]]
    rxnEquations[["elasticity_calc"]] <- gsub('(I\\()', Kcatpaste, rxnEquations[["elasticity_calc"]])
    rxnEquations[["elasticity_calc"]] <- sub('I\\(', '\\(', rxnEquations[["elasticity_calc"]])
    
  }
  
  # find the partial derivatives of the kinetic form for each reaction specie
  eq <- eval(parse(text = paste('expression(',rxnEquations[["full_kinetic_form_list"]],')')))
  
  D_full_kinetic_form <- list()
  for(spec in c(kineticPars$rel_spec[(kineticPars$measured & !is.na(kineticPars$measured)) | kineticPars$modelName == "t_metX" & !is.na(kineticPars$modelName)], all_species$rel_spec[all_species$SpeciesType == "Enzyme"])){
    D_full_kinetic_form[[spec]] <- D(eq, spec)
  }
  rxnEquations[["kinetic_form_partials"]] <- D_full_kinetic_form
  D_full_kinetic_form <- D_full_kinetic_form[names(D_full_kinetic_form) %in% c(kineticPars$rel_spec[(kineticPars$measured & !is.na(kineticPars$measured))], all_species$rel_spec[all_species$SpeciesType == "Enzyme"])] # remove un-measured species, as these will have no variance
  
  #### Generate a prior for each non-linear kinetic parameter (i.e. not kcat) ####
  
  kineticParPrior <- data.frame(distribution = rep(NA, times = nrow(kineticPars)), par_1 = NA, par_2 = NA, par_3 = NA) 
  
  # Options for these parameters are:
  # -@-@ unif: uniform in log-space: par_1 = lowerbound, par_2 = upperbound. draw in log2 space and exponentiate back to linear space
  # -@-@ norm: lognormal: in log2 space draw a value from mean = par_1, sd = par_2
  # -@-@ SpSl: spike and slab (In log2 space: the spike is a point mass at zero with p = par_3, the slab is a normal with mean = par_1 and sd = par_2)
  
  # Specify prior for michaelis constants
  kineticParPrior$distribution[kineticPars$SpeciesType %in% c("Metabolite", "keq")] <- "unif"
  kineticParPrior$par_1[kineticParPrior$distribution == "unif"] <- -15; kineticParPrior$par_2[kineticParPrior$distribution == "unif"] <- 15 # default value to 2^-15:2^15
  
  for(exp_param in kineticPars$formulaName[!is.na(kineticPars$measured) & kineticPars$measured == TRUE]){
    kineticParPrior[kineticPars$formulaName == exp_param & !is.na(kineticPars$modelName), c(2:3)] <- median(met_abund[,colnames(met_abund) == kineticPars$modelName[kineticPars$formulaName == exp_param & !is.na(kineticPars$formulaName)]]) + c(-15,15)
  }#priors for measured metabolites (either in absolute or relative space) are drawn about the median
  
  # Prior for keq is centered around sum log(sub) - sum log(prod) - this deals with some species being in absolute space and others being absolute measurements
  kineticParPrior[kineticPars$SpeciesType == "keq", c(2:3)] <- median(rowSums(met_abund * c(rep(1,n_c)) %*% t(rxnSummary$rxnStoi[shmatch(colnames(met_abund), names(rxnSummary$rxnStoi))]))) + c(-20,20)
  
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
  
  
  #### Optimize l(par|X) using NLS ####
  
  rxn_data <- cbind(flux = flux$FVAmin + flux$FVAmax, met_abund, enzyme_abund)
  
  ##
  
  nls_parameters <- rbind(
    data.frame(formulaName = all_species$formulaName[all_species$SpeciesType == "Enzyme"], SpeciesType = "kcat", lower = 0, upper = Inf),
    data.frame(kineticPars[,colnames(kineticPars) %in% c("formulaName", "SpeciesType")], lower = 2^kineticParPrior$par_1, upper = 2^kineticParPrior$par_2)
  )
  nls_parameters$lower[nls_parameters$SpeciesType == "hillCoefficient"] <- 2^-3
  nls_parameters$upper[nls_parameters$SpeciesType == "hillCoefficient"] <- 2^3
  
  # generate a regression formula
  nls_formula <- gsub('^', 'flux ~ ', rxnEquations$full_kinetic_form_list)
  
  nls_parameters <- nls_parameters[sapply(nls_parameters$formulaName, function(x){grepl(x, nls_formula)}),]
  
  N_nls_runs <- 100
  nls_initializations <- matrix(NA, ncol = N_nls_runs, nrow = nrow(nls_parameters))
  
  for(i in 1:nrow(nls_parameters)){
    if(nls_parameters$SpeciesType[i] == "kcat"){
      nls_initializations[i,] <- 1
    }else{
      nls_initializations[i,] <- 2^runif(N_nls_runs, log2(nls_parameters[i,'lower']), log2(nls_parameters[i,'upper']))
    }
  }
  
  nls_runs <- list()
  for(j in 1:ncol(nls_initializations)){
    
    nls_runs[[j]] <- try(nls(formula = nls_formula,
                         data = rxn_data, algorithm = "port",
                         start = setNames(nls_initializations[,j], nls_parameters$formulaName),
                         lower = setNames(nls_parameters$lower, nls_parameters$formulaName),
                         upper = setNames(nls_parameters$upper, nls_parameters$formulaName)), silent = T)
  }
  
  valid_output <- sapply(nls_runs, function(x){class(x) == "nls"})
  if(sum(valid_output) == 0){next}else{
    
    NLS_logLik <- sapply(nls_runs[valid_output], function(x){
      
      fit_resid_error <- sqrt(mean((residuals(x) - mean(residuals(x)))^2))*sqrt(n_c/(n_c-1))
      #fit_resid_error <- sqrt(mean(residuals(x)^2) * (n_c/(n_c-1)))
      sum(dnorm(residuals(x), 0, fit_resid_error, log = T))
      
    }) %>% max()
    
    print("valid fit")
    NLS_fit_list[[rxN]] <- data.frame(rxN = rxN, logLik = NLS_logLik)
    
  }

}

NLS_fit_list <- rbind_all(NLS_fit_list)

NLS_logLik <- NLS_rxns %>% left_join(NLS_fit_list, by = c("row" = "rxN")) %>% filter(!is.na(logLik))

NLScompare <- NLS_logLik %>% left_join(all_reactionInfo %>% tbl_df() %>% select(rMech, ML), by = c("rxnID" = "rMech")) 

library(ggplot2)

ggplot(NLScompare, aes(x = logLik, y = ML)) + geom_point() + geom_abline(a = 0, b = 1)

load("flux_cache/modelComparison.Rdata")
all_reactionInfo %>% tbl_df()


