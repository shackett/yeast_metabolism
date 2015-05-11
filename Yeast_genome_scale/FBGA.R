#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2) #for visualization at the end
library(data.table)
library(nnls) #for non-negative regression used to fit kinetic parameters
library(ggplot2)
library(gplots) 
library(grid) # minor use to adjust layout of ggplot2 facets
library(coda) # for assessing mcmc convergence
library(dplyr)
source("FBA_lib.R")
library(scales) # for quick % conversion

options(stringsAsFactors = FALSE)

######## Import all of the species involved in a reaction and other reaction information ###########

load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)
n_c <- nrow(chemostatInfo)

##### Import list of metabolite abundances and rxn forms - **rxnf**

##### to construct a list of reaction information, including a mechanism.
# Match Boer metabolites to tIDs using chemical similarity (ChEBI + ChemmineR)
# Using BRENDA data, determine all species which could activate or inhibit our rxns and were ascertained via metabolomics
### KEGG RIDs -> EC -> modifiers, activators and ligands
### Determine the Ki and Km of inhibitors
# Use chemical similarity to match substrates and products for each reaction to determine where competitive inhibition of the product will act
# Generate a list of reaction form parameterizations
### Using reversible MM kinetics as the base form
### Additional reactions with BRENDA activators and inhibitors
### To search for novel allosteric modifiers in an unsupervised manner - allow for reaction with extensions an abundance profile that is not specified a priori

if(!file.exists("flux_cache/rxnf_formulametab.rdata")){
  source("reactionStructures.R")
}else{
  load("flux_cache/rxnf_formulametab.Rdata")
}

##### Import list of flux fluxes from FBA_run_full_reco.R - **flux_summary**
# IDs - reaction ID to reaction name correspondence for reactions with non-zero flux in the standard QP optimization
# cellularFluxes - intracellular fluxes (moles per hr per mL cellular volume) across all 25 nutrients conditions for reactions with non-zero flux in some condition with the standard QP optimization.
# fva_summary - rxns ~ conditions ~ min and maximum bound at QP solution - reactions are those that are moderately constrained under a majority of conditions
# total_flux_cast - rxn ~ conditions ~ min, max and QP solution -> same as fva_summary but with intersection of cellular fluxes boiled in

load("flux_cache/fluxSummaryQP.Rdata")


##### Associate enzymes with pathways ######

### Determine which pathways are associated with a protein ###
# This isn't used until FBGA.R, but is convenient to generate at this point

kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv") # > generated from keggIDparser.py KEGG IDs relative to yeast gene name (1gene -> 1 KEGG id, multiple  mapping between KEGG and genes)
if(is.null(kegg_enzyme_dict$PATHWAY)){
  gene_pathways(kegg_enzyme_dict); kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv")
}#generate a per-gene pathway annotation if one is not already generated
####
rxn_pathways <- read.delim("./flux_cache/reactionPathways.tsv")
rxn_enzyme_measured <- read.delim("./flux_cache/prot_matches.tsv")


####### Narrow reactions txwo those that are well-constrained & non-zero in a majority of conditions #####
### When considering reactions which carry zero flux under a subset of conditions, analyze reactions both when v = 0, and by removing conditions where v = 0 ###
rmCondList <- data.frame() # reactions which will be duplicated with some conditions removed


# valid reactions have well-constrained, non-zero fluxes and measured enzymes (a reaction possessing a minimal complement of metabolites is enforced in reactionStructures.R -> rxnf"
valid_rxns <- grep('r_[0-9]{4}', rownames(flux_summary$total_flux_cast), value = T) # valid flux
valid_rxns <- valid_rxns[valid_rxns %in% rxn_enzyme_measured$reaction[rxn_enzyme_measured$measured]] # valid proteins
valid_rxns <- valid_rxns[valid_rxns %in% unique(substr(names(rxnf), 1, 6))] # valid metabolites

rxnList <- rxnf[substr(names(rxnf),1,6) %in% valid_rxns]

for(rxN in 1:length(valid_rxns)){
  
  idx <- names(rxnList)[grep(valid_rxns[rxN],names(rxnList))]
  if(length(idx) ==0){print("no matches"); next}
  
  if(valid_rxns[rxN] %in% rxn_pathways$reactionID){
    Rxpathway <- rxn_pathways$pathway[rxn_pathways$reactionID == valid_rxns[rxN]]
  }else{Rxpathway <- ''}
  
  rxn_genes <- strsplit(rxnList[[idx[1]]]$genes, split = '/')[[1]]
  if(any(rxn_genes %in% kegg_enzyme_dict$SYST)){
    geneInfo <- kegg_enzyme_dict[chmatch(rxn_genes[rxn_genes %in% kegg_enzyme_dict$SYST], kegg_enzyme_dict$SYST),]
  }else{geneInfo <- NA}
  
  rxFlux <- as.data.frame(flux_summary$total_flux_cast[rownames(flux_summary$total_flux_cast) == valid_rxns[rxN],,])
  
  #Add the rxns with lacking fluxes to the rmCondList
  if (any(rxFlux$standardQP == 0)){
    rmCondList <- rbind(rmCondList,
                        data.frame(rxn = valid_rxns[rxN],
                                   cond = paste(rownames(rxFlux)[rxFlux$standardQP == 0],collapse=';'),
                                   source = 'zeroFluxes', nZero = sum(rxFlux$standardQP == 0)))
  }
  # manually flagged orotate overflow reactions and test removal
  if (valid_rxns[rxN] %in% c("r_0214", "r_0250", "r_0349", "r_0453", "r_4043")){
    rmCondList <- rbind(rmCondList,
                        data.frame(rxn = valid_rxns[rxN],
                                   cond = paste(grep('^U', rownames(rxFlux), value = T), collapse = ';'),
                                   source = 'overflow', nZero = length(grep('^U', rownames(rxFlux)))))
  }
  # manually flagged leucine dysregulated
  if (valid_rxns[rxN] %in% c("r_0096", "r_0097", "r_0352", "r_669", "r_0816")){
    rmCondList <- rbind(rmCondList,
                        data.frame(rxn = valid_rxns[rxN],
                                   cond = paste(grep('^L', rownames(rxFlux), value = T), collapse = ';'),
                                   source = 'dysregulation', nZero = length(grep('^L', rownames(rxFlux)))))
  }
  
  
  for (entry in idx){
    rxnList[[entry]]$reaction = flux_summary$IDs$Name[flux_summary$IDs$reactionID == valid_rxns[rxN]]
    rxnList[[entry]]$pathway = Rxpathway
    rxnList[[entry]]$geneInfo = geneInfo
    rxnList[[entry]]$flux <- rxFlux
  }
}


# Apply the changes of rmCondList
# add a copied entry without the mentioned conditions for all rows in the rm

rmCondList <- rmCondList[rmCondList$nZero < 10 & rmCondList$nZero > 2,]

for (i in 1:nrow(rmCondList)){
  rxn <- rmCondList$rxn[i]
  for (entry in names(rxnList)[grep(rxn,names(rxnList))]){
    nEntry <- paste(entry,'_rmCond',sep='')
    rxnList[[nEntry]] <- rxnList[[entry]]
    conds <- strsplit(rmCondList$cond[i],';')[[1]]
    rxnList[[nEntry]]$flux <- rxnList[[nEntry]]$flux[!rownames(rxnList[[nEntry]]$flux) %in% conds,]
    rxnList[[nEntry]]$rxnMet <- rxnList[[nEntry]]$rxnMet[!rownames(rxnList[[nEntry]]$rxnMet) %in% conds,]
    rxnList[[nEntry]]$all_species_SD <- rxnList[[nEntry]]$all_species_SD[!rownames(rxnList[[nEntry]]$all_species_SD) %in% conds,]
    rxnList[[nEntry]]$enzymeComplexes <- rxnList[[nEntry]]$enzymeComplexes[,!colnames(rxnList[[nEntry]]$enzymeAbund) %in% conds]
  
    rxnList[[nEntry]]$listEntry <- paste(rxnList[[nEntry]]$listEntry, '_rmCond',sep='')
    
    }
}

# remove irreversible reaction forms that are polarized in the wrong direction

irreversible_rxns <- data.frame(rxForm = names(rxnList)[grep('-im-', names(rxnList))])
irreversible_rxns$rxn <- substr(irreversible_rxns$rxForm, 1, 6)
irreversible_rxns$dir <- ifelse(grepl('-forward', irreversible_rxns$rxForm), "F", "R")

reaction_direction <- data.frame(npos = apply(flux_summary$total_flux_cast[,,'standardQP'], 1, function(x){sum(x > 0)}),
                                 nneg = apply(flux_summary$total_flux_cast[,,'standardQP'], 1, function(x){sum(x < 0)}))

irreversible_rxns$invalid <- (irreversible_rxns$rxn %in% rownames(reaction_direction)[reaction_direction$npos == 0] & irreversible_rxns$dir == "F")|
  (irreversible_rxns$rxn %in% rownames(reaction_direction)[reaction_direction$nneg == 0] & irreversible_rxns$dir == "R")

rxnList <- rxnList[!(names(rxnList) %in% irreversible_rxns$rxForm[irreversible_rxns$invalid])]
                       
### ensure that the ordering of conditions is the same ###

cond_mapping <- data.frame(standard = chemostatInfo$ChemostatCond, e = colnames(rxnList[[1]]$enzymeComplexes), m = rownames(rxnList[[1]]$rxnMet), f = rownames(rxnList[[1]]$flux))
if (!all(apply(apply(cond_mapping, c(1,2), toupper), 1, function(x){length(unique(x))}) == 1)){
  warning('There is a problem with the order of the conditions. (check cond_mapping)')
}

#### Determine which reaction have valid reaction mechanisms - all of them as of now ####

reactionForms <- sapply(rxnList, function(x){ifelse(!is.null(x$rxnForm), x$listEntry, NA)})  
rxnList_form <- rxnList[names(rxnList) %in% reactionForms]


rxnList_form <- rxnList_form[order(names(rxnList_form))] # order alpha-numerically so that the workload is spread out more evenly during optimization

#### save rxnList_form so that this self-sufficient list can be thrown at the cluster ###

# chunks to break rxnList into
chunk_size <- 100
chunk_assignment <- data.frame(set = names(rxnList_form), chunk = c(rep(1:floor(length(rxnList_form)/chunk_size), each = chunk_size), rep(ceiling(length(rxnList_form)/chunk_size), length(rxnList_form) %% chunk_size)))
print(paste("The number of parameter chunks is", ceiling(length(rxnList_form)/chunk_size), "submit this parameter when batch submitting processes in FluxOptim.sh"))

save(rxnList_form, cond_mapping, chunk_assignment, file = "paramOptim.Rdata")


# How well can flux be explained through regression of metabolites and enzyme on flux
flux_mass_action_regression()





##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@
#######* Import cluster parameter results #########
##@##@ Start here if loading parameter sets ##@##@#
##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@

load("paramOptim.Rdata")

param_sets <- list.files('FBGA_files/paramSets/')

param_run_info <- data.frame(file = param_sets, index = 1:length(param_sets), chunkNum =sapply(param_sets, function(x){as.numeric(sub('C', '', regmatches(x, regexec('C[0-9]+', x))[[1]]))}),
                             runNum = sapply(param_sets, function(x){as.numeric(sub('R', '', regmatches(x, regexec('R[0-9]+', x))[[1]]))}), n_samples = NA, sample_freq = NA, burn_in = NA)

### combine parameter runs together into a list indexed by reaction name ####

param_set_list <- list()

k <- 1
for(a_chunk in sort(unique(param_run_info$chunkNum))){
  
  for(an_index in param_run_info$index[param_run_info$chunkNum == a_chunk]){ #temporarily only loading run number 1 (1 tenth of data)
    load(paste(c("FBGA_files/paramSets/", param_run_info$file[an_index]), collapse = ""))
    param_run_info$n_samples[an_index] <- run_summary$markov_pars$n_samples
    param_run_info$sample_freq[an_index] <- run_summary$markov_pars$sample_freq
    param_run_info$burn_in[an_index] <- run_summary$markov_pars$burn_in
    
    for(rx in names(run_summary)[-1]){
      param_set_list[[as.character(k)]]$name <- data.frame(rx = rx, index = an_index, chunk = a_chunk)
      param_set_list[[as.character(k)]]$lik <- run_summary[[rx]]$likelihood
      param_set_list[[as.character(k)]]$MC <- run_summary[[rx]]$markovChain
      k <- k + 1
    }
  }
}

parSetInfo <- as.data.frame(t(sapply(param_set_list, function(x){x$name})))
parSetInfo$ML <- sapply(param_set_list, function(x){max(x$lik)})

#### Generate technical summaries of MCMC performance ####

if(length(unique(param_run_info$n_samples)) != 1 | length(unique(param_run_info$sample_freq)) != 1 | length(unique(param_run_info$burn_in)) != 1){
  stop("Run parameters differ between chains")
  }

all_rx <- unique(unlist(parSetInfo$rx))
chain_convergence <- sapply(all_rx, function(arxn){
  parSubset <- param_set_list[parSetInfo$rx == arxn]
  gelman.diag(as.mcmc.list(lapply(parSubset, function(x){as.mcmc(x$MC)})), autoburnin = F)$mpsrf # determine convergence of markov chains
  # gelman.plot(as.mcmc.list(lapply(parSubset, function(x){as.mcmc(x$MC)})))
})


#### Pull out relevent reaction information - such as reaction ID, number of parameters and maximum likelihood ####

load('flux_cache/metaboliteTables.RData')
npc <- ncol(metSVD$v)

# npar = number of enzymes (kcat) + number of metabolites (Km) + hill coefficients (which are not 1) 
# + hypothetical metabolite (trend governed by npc significant PCs) + Keq

reactionInfo <- data.frame(rMech = names(rxnList_form), reaction = sapply(names(rxnList_form), function(x){substr(x, 1, 6)}),
                           form = sapply(names(rxnList_form), function(x){substr(x, 8, 9)}), modification = sapply(names(rxnList_form), function(x){sub(substr(x, 1, 10), '', x)}),
                           ncond = sapply(rxnList_form, function(x){nrow(x$flux)}),
                           npar = sapply(rxnList_form, function(x){nrow(x$enzymeAbund) + nrow(x$rxnFormData) + sum(x$rxnFormData$Hill != 1) + ifelse(any(x$rxnFormData$SubstrateID == "t_metX"), npc, 0) + 1}),
                           MPSRF = chain_convergence)

reactionInfo$ML <- sapply(reactionInfo$rMech, function(x){max(parSetInfo$ML[parSetInfo$rx == x])}) # the maximum likelihood over all corresponding parameter sets

### Model comparison between more complex models with additional regulators or alterative reaction forms relative to reversible-MM using LRT and AICc ###

reaction_signif <- modelComparison(reactionInfo, rxnList_form)

reactionInfo <- reactionInfo %>% dplyr::select(-form) %>%
  left_join(reaction_signif %>% dplyr::select(-signifCode), by = c("reaction", "rMech"))

### Base reaction specific information ###

rxToPW <- NULL

for(rxN in grep('rm$', names(rxnList_form))){
  
  RXannot <- rxnList_form[[rxN]]$pathway
  RXannot <- sub('[_]{0,2}Metabolic pathways', '', RXannot)
  GENEannot <- rxnList_form[[rxN]]$geneInfo
  if(!(all(is.na(GENEannot)))){
    GENEannot$PATHWAY <- sub('[_]{0,2}Metabolic pathways', '', GENEannot$PATHWAY)
    if(all(GENEannot$PATHWAY == "")){
      GENEannot <- NA
      }
    }
  
  if(RXannot == "" | is.na(RXannot)){
    if(all(is.na(GENEannot))){
      pathways <- "Not annotated"
      RELannot <- NA
    }else{
      RELannot <- GENEannot$PATHWAY
      if(all(RELannot == "") |all(is.na(RELannot))){
      stop(paste("error with pathway assignment for reaction:", names(rxnList_form)[rxN]))  
      }
    }
  }else{
    RELannot <- RXannot
  }
  
  if(all(!is.na(RELannot))){
    pathways <- strsplit(RELannot, split = '__')[[1]]
    pathways <- pathways[pathways != ""]
  }
  
  rxToPW <- rbind(rxToPW, data.frame(rxN = rxN, rID = rxnList_form[[rxN]]$rxnID, reactionName = rxnList_form[[rxN]]$reaction, pathway = pathways))
    
}  

# Remove some minor pathway annotations
minorPW <- names(table(rxToPW$pathway))[table(rxToPW$pathway) <= 3]
otherPW <- rxToPW$rID[!(rxToPW$rID %in% unique(rxToPW$rID[!(rxToPW$pathway %in% minorPW)]))]

otherPWreactions <- rxToPW[rxToPW$rID %in% otherPW,]
otherPWreactions$pathway <- "Other"

rxToPW <- rbind(rxToPW[!(rxToPW$pathway %in% minorPW),], unique(otherPWreactions))

rxToPW <- rbind(rxToPW, data.frame(unique(rxToPW[,1:3]), pathway = "ALL REACTIONS"))

pathwaySet <- sort(table(rxToPW$pathway), decreasing = T)
pathwaySet <- data.frame(pathway = names(pathwaySet), members = unname(pathwaySet), display = paste(names(pathwaySet), ' (', unname(pathwaySet), ')', sep = ""))


#### Generate reaction plots and summaries ####

# For a subset of reactions

all_reactionInfo <- reactionInfo
reactionInfo <- all_reactionInfo %>% filter(modelType %in% c("rMM", "irreversible", "hypo met regulator") | Qvalue < 0.1)

load("companionFiles/PTcomparison_list.Rdata") # by-gene comparisons of protein and transcript abundance

shiny_flux_data <- list()
rxn_fits <- NULL
rxn_fit_params <- list()
fraction_flux_deviation <- NULL
MLdata <- NULL # Save summary of metabolic leverage
ELdata <- list() # Save full distribution of elasticities
TRdata <- NULL # Save summary transcriptional resposiveness
Hypo_met_candidates <- NULL

t_start = proc.time()[3]

# flag a few reactions where additional plots are created
custom_plotted_rxns <- c("r_1054-im-forward", "r_1054-rm","r_0962-rm", "r_0962-rm-t_0290-act-mm", "r_0816-rm", "r_0816-rm-t_0461-inh-uncomp", "r_0816-rm_rmCond", "r_0816-rm-t_0461-inh-uncomp_rmCond",
                         "r_0514-im-forward", "r_0514-rm", "r_0916-im-forward", "r_0916-rm", "r_0208-im-forward", "r_0208-rm")

#custom_plotted_rxns <- c(custom_plotted_rxns, reactionInfo$rMech[reactionInfo$modification %in% c("", "rmCond")])
custom_plotted_rxns <- unique(custom_plotted_rxns)

if(!(all(custom_plotted_rxns %in% reactionInfo$rMech))){stop("Some reaction mechanisms that you want to plot were not located")}

#arxn <- c("r_0214-rm")
#for(arxn in rxn_subset){
#for(arxn in custom_plotted_rxns){
for(arxn in reactionInfo$rMech){
  
  par_likelihood <- NULL
  par_markov_chain <- NULL
  
  parSubset <- param_set_list[parSetInfo$rx == arxn]
  
  for(i in 1:length(parSubset)){
    
    par_likelihood <- rbind(par_likelihood, data.frame(sample = 1:param_run_info$n_samples[parSubset[[i]]$name$index], likelihood = parSubset[[i]]$lik, index = parSubset[[i]]$name$index))
    par_markov_chain <- rbind(par_markov_chain, parSubset[[i]]$MC)
  }
  
  load(paste(c("FBGA_files/paramSets/", param_run_info$file[param_run_info$index == par_likelihood$index[1]]), collapse = ""))
  run_rxn <- run_summary[[arxn]]
  
  ### A couple of catches for poorly defined reactions or if optimization grossly fails ### 
  
  if(sum(is.finite(par_likelihood$likelihood)) == 0){# all parameter sets have a zero likelihood (e.g. SD = 0)
    warning(paste(arxn, "has zero valid parameter sets")); next
  }
  
  
  if(var(par_likelihood$likelihood) < 10^-10 | all(run_rxn$metabolites == 1)){
    warning(paste(arxn, "does not vary in likelihood, or has too many unmeasured species")); next
  } #skip underparameterized reactions - those with essentially no variation
  
  
  ### Determine how closely fluxes predicted from a parametric form match fluxes determined via constraint based modeling ###
  
  flux_fit <- flux_fitting(run_rxn, par_markov_chain, par_likelihood) #compare flux fitted using the empirical MLE of parameters
  rxn_fits <- rbind(rxn_fits, data.frame(rxn = arxn, flux_fit$fit_summary))
  rxn_fit_params[[arxn]] <- flux_fit$param_interval
  
  ### Take the fractional flux departure and dot-product between FBA and parametric vector
  
  vector_match <- data.frame(rxn = arxn, FFD = 1 - sum(abs(flux_fit$fitted_flux$fitted - (run_rxn$flux$FVAmin + run_rxn$flux$FVAmax)/2))/sum(abs(run_rxn$flux)),
                             dotProduct = sum(flux_fit$fitted_flux$fitted/sqrt(sum((flux_fit$fitted_flux$fitted)^2)) * (run_rxn$flux$FVAmin + run_rxn$flux$FVAmax)/2/sqrt(sum(((run_rxn$flux$FVAmin + run_rxn$flux$FVAmax)/2)^2))))
  vector_match$angle <- acos(vector_match$dotProduct) * 180/pi
  
  # Performance based on fraction of FVA intervals captured by parameteric 95% CI #
  fluxIntervals <- data.frame(VLB = run_rxn$flux$FVAmin, VUB = run_rxn$flux$FVAmax, PLB = flux_fit$fitted_flux$fitted - 2*flux_fit$fitted_flux$SD, PUB = flux_fit$fitted_flux$fitted + 2*flux_fit$fitted_flux$SD)
  fluxOverlap <- (mapply(function(VUB, PUB){min(VUB, PUB)}, VUB = fluxIntervals$VUB, PUB = fluxIntervals$PUB) - mapply(function(VLB, PLB){max(VLB, PLB)}, VLB = fluxIntervals$VLB, PLB = fluxIntervals$PLB))/
    mapply(function(VI, PI){min(VI, PI)}, VI = fluxIntervals$VUB - fluxIntervals$VLB, PI = fluxIntervals$PUB - fluxIntervals$PLB)
  fluxOverlap[fluxOverlap < 0] <- 0
  vector_match$"Interval Overlap" <- mean(fluxOverlap)
  
  fraction_flux_deviation <- rbind(fraction_flux_deviation, vector_match)
  
  ### Generate plots which show reaction information, species variation, flux fitting ... ###
  
  shiny_flux_data[[arxn]]$reactionInfo <- reaction_info_FBGA(rxnName) # Reaction information
  
  shiny_flux_data[[arxn]]$plotChoices$Likelihood <- likViolin(par_likelihood, run_summary$markov_pars) # Log-likelihoods of each markov chain
  
  species_plots <- species_plot(run_rxn, flux_fit, chemostatInfo)
  shiny_flux_data[[arxn]]$plotChoices <- append(species_plots, shiny_flux_data[[arxn]]$plotChoices)
  
  if("t_metX" %in% run_rxn$kineticPars$modelName){
    hypo_met_info <- hypoMetTrend(run_rxn, metSVD, tab_boer)
    
    Hypo_met_candidates <- rbind(Hypo_met_candidates, hypo_met_info[["specCorrQuantiles_all"]])
    
    shiny_flux_data[[arxn]]$plotChoices <- append(shiny_flux_data[[arxn]]$plotChoices, hypo_met_info[names(hypo_met_info) != "specCorrQuantiles_all"])
  }else if(any(run_rxn$kineticPars$SpeciesType == "hillCoefficient")){
    shiny_flux_data[[arxn]]$plotChoices$Hill <- hillPlot(run_rxn)
    }
  
  reaction_properties <- reactionProperties() # Evaluate elasticities for each markov sample.  Join metabolite variation with elasticites
  
  MLdata <- rbind(MLdata, reaction_properties$ML_summary)
  ELdata[[arxn]] <- reaction_properties$EL_summary
  
  reaction_plots <- reaction_properties$plots
  shiny_flux_data[[arxn]]$plotChoices <- append(shiny_flux_data[[arxn]]$plotChoices, reactionPropertiesPlots(reaction_plots))
  
  trans_res <- transcriptional_responsiveness()
  if("Plots" %in% names(trans_res)){ shiny_flux_data[[arxn]]$plotChoices <- append(shiny_flux_data[[arxn]]$plotChoices, trans_res$Plots) }
  if("TR" %in% names(trans_res)){ TRdata <- rbind(TRdata, trans_res$TR) }
  
  if(arxn %in% custom_plotted_rxns){
    # Generate custom plots for a subset of reactions
    customPlots(run_rxn, flux_fit, chemostatInfo)
    param_dist <- param_compare()
    ggsave(param_dist$bivariateHist, file = paste0("Figures/rxnPlots/", arxn, "-BiparamHist", ".pdf"),  width = min(ncol(par_markov_chain)*4.5, 49), height = min(ncol(par_markov_chain)*4.5, 49))
    ggsave(param_dist$univariateHist, file = paste0("Figures/rxnPlots/", arxn, "-UniparamHist", ".pdf"),  width = min(ncol(par_markov_chain)*5.5, 49), height = 6)
    }
  
  #if(grepl('0816', arxn)){
    # save parameter distributions for a subset of reactions
    #param_dist <- param_compare()
    #ggsave(param_dist, file = paste0("tmp/", arxn, "_paramHist.pdf"), width = 20, height = 20)
    #ggsave(shiny_flux_data[[arxn]]$plotChoices$Likelihood, file = paste0("tmp/", arxn, "_likViolin.pdf"), width = 16, height = 10)
  #}
  #shiny_flux_data[[arxn]]$plotChoices$"Parameter Comparison" <-  param_compare() # these plots makes the resulting list huge, but they are awesome ...
  
  # when all of the reaction forms for a given reaction have been plotted, then export them as a group
  rID = reactionInfo$reaction[reactionInfo$rMech == arxn]
  if(arxn == last(reactionInfo$rMech[reactionInfo$reaction == rID])){
    save(shiny_flux_data, file = paste0("shinyapp/reaction_data/", rID, "plots.Rdata"))
    shiny_flux_data <- list()  
  }
  
  if(which(reactionInfo$rMech == arxn) %% 10 == 0){
    cat(paste('\n',round((which(reactionInfo$rMech == arxn) / length(reactionInfo$rMech))*100, 2), "% complete - ", round((proc.time()[3] - t_start)/60, 0), " minutes elapsed", sep = ""))
  }
  
}; cat("\nDone!")

#for(a_name in names(shiny_flux_data)){
#  a_rxn_file <- shiny_flux_data[[a_name]][[2]]
#  for(a_plot in names(a_rxn_file)){
#    a_plot_name <- gsub('[^A-Za-z]', '', a_plot)
#    ggsave(file = paste0("tmp/", a_name, "==", a_plot_name, ".pdf"), plot = a_rxn_file[names(a_rxn_file) == a_plot][[1]])
#    }
#  }

# Summarize consistency of reaction forms for a single reaction
rxn_plot_list <- list()
for(rxn in unique(reactionInfo$reaction)){
  rxn_plot_list[[rxn]] <- modelComparisonPlots(rxn, reactionInfo, all_reactionInfo)
}

# Summarize consistency of reaction forms for all reactions in a pathway
pathway_plot_list <- list()
for(pw in pathwaySet$display){
  # iterate through pathways and plot pathway-level figures
  pathway_plot_list[[pw]] <- pathwayPlots(pw)
  }


#### Save lists which will be processed by Shiny app ####

save(pathwaySet, rxToPW, pathway_plot_list, rxn_plot_list, reactionInfo, file = "shinyapp/shinyData.Rdata")

# generate a minute version of shinyData that will load quickly when the App is being modified
#reactionInfo <- reactionInfo[1:20,]
#shiny_flux_data <- shiny_flux_data[names(shiny_flux_data) %in% reactionInfo$rMech]
#system.time(save(pathwaySet, rxToPW, reactionInfo, pathway_plot_list, shiny_flux_data, file = "shinyapp/shinySubData.Rdata"))


#### Save parameter estimates for further global analyses ####

save(rxn_fit_params, rxn_fits, reactionInfo, all_reactionInfo, MLdata, TRdata, fraction_flux_deviation, file = "flux_cache/paramCI.Rdata")
save(ELdata, file = "flux_cache/elasticityData.Rdata")

##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@###@###@###@###@###@###@###@
###################### * Summarize Reaction Behavior ##########################
##@##@ Start here if loading parameter estimates and kinetic summaries ##@##@##
##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@###@###@###@###@###@###@###@

##### Systems level comparison of optimized and external parameter values #####

load("flux_cache/paramCI.Rdata")
load("flux_cache/reconstructionWithCustom.Rdata")
reversibleRx <- read.table("companionFiles/reversibleRx.tsv", header = T)

# Test if there were any reactions that were discarded on the fly

if(sum(!(reactionInfo$rMech %in% rxn_fits$rxn)) != 0){
  reactionInfo <- reactionInfo %>% dplyr::filter(rMech %in% rxn_fits$rxn)
  }



### Determine which reactions to follow-up on based upon either
## A) Contain all substrates
## B) RM or regulation is well-fit despite missing substrates

load("paramOptim.Rdata")

RMMrxns <- reactionInfo$rMech[reactionInfo$modification == ""]
measure_exception <- c("H+", "H2O")

validRxnA <- sapply(RMMrxns, function(i){
  x <- rxnList_form[[i]]
  if(all(x$flux$standardQP >= 0)){
    substrates <- x$originMet[chmatch(names(x$rxnStoi)[x$rxnStoi < 0], names(x$originMet))]
  }else if(all(x$flux$standardQP <= 0)){
    substrates <- x$originMet[chmatch(names(x$rxnStoi)[x$rxnStoi > 0], names(x$originMet))]
  }else{
    substrates <- x$originMet
  }
  
  substrates <- substrates[names(substrates) %in% names(x$metNames)[!(x$metNames %in% measure_exception)]]
  if(!any(substrates == "nm")){
   x$rxnID
  }
})
validRxnA <- unlist(validRxnA) %>% unname()

###

validRxnB <- rxn_fits %>% tbl_df() %>% dplyr::select(rxn, parSpearman) %>% filter(parSpearman > 0.3) %>% left_join(reactionInfo %>% tbl_df() %>% dplyr::select(rxn = rMech, reaction, modification, Qvalue)) %>%
  filter(!grepl('t_metX', rxn)) %>% filter(is.na(Qvalue) | Qvalue < 0.1) %>% dplyr::select(reaction) %>% unlist() %>% unname() %>% unique()

valid_rxns <- union(validRxnA, validRxnB)
rmCond_rxns <- unique(reactionInfo$reaction[grep('rmCond', reactionInfo$modification)]) # reactions which carry zero flux under some conditions - consider only non-zero reactions

optimal_rxn_form <- sapply(valid_rxns, function(x){
  
  rx_forms <- reactionInfo[reactionInfo$reaction == x,]
  if(x %in% rmCond_rxns){
    rx_forms <- rx_forms[grep('rmCond', rx_forms$modification),]
    }
  rx_forms <- rx_forms %>% dplyr::filter(modelType %in% c("rMM", "regulator", "coopertivity", "2+ regulators"))
  
  # determine whether any regulator is significantly better than rMM
  bestReg <- rx_forms %>% filter(modelType == "regulator") %>% filter(ML == max(ML)) %>%
    mutate(AICc = 2*npar - 2*ML + 2*npar*(npar + 1)/(ncond - npar - 1))
  bestMulti <- rx_forms %>% filter(modelType %in% c("2+ regulators", "coopertivity")) %>% filter(Qvalue < 0.05) %>% filter(ML == max(ML)) %>%
    mutate(AICc = 2*npar - 2*ML + 2*npar*(npar + 1)/(ncond - npar - 1))
  
  # is the multiple regulation model substantially better than the single regulation model
  if(nrow(bestReg) == 1 & nrow(bestMulti) == 1){
    if(1 - 1/(exp((bestMulti$AICc - bestReg$AICc)/2) + 1) < 0.1){
      return(bestMulti$rMech)
    }
  }
  if(nrow(bestReg) == 1){
   if(bestReg$Qvalue < 0.05){
     return(bestReg$rMech)
   }
  }
  return(rx_forms$rMech[rx_forms$modelType == "rMM"])
  
})
  
# when multiple regulatory candidates have similar fits but one is a far better candidate based upon the literature, choose the literature-supported one
reactionInfo %>% filter(reaction == "r_0962")
optimal_rxn_form[grep('r_0962-rm', optimal_rxn_form)] <- "r_0962-rm-t_0290-act-mm" # pyruvate kinase activation by F16bisP over citrate inhibition
reactionInfo %>% filter(reaction == "r_0215") 
optimal_rxn_form[grep('r_0215-rm', optimal_rxn_form)] <- "r_0215-rm-t_0499-inh-uncomp_ultra" # aspartate kinase regulation by ultrasensitive inhbition by threonine

# remove reactions where major substrates are missing and a high-influence regulator is suggested likely mimicking this missing substrate
flawed_rxns <- c("r_0990", "r_0888", "r_0534", "r_0988", "r_0032", "r_0510")
# add several reactions that fail to capture some large aspect of the flux variability
flawed_rxns <- c(flawed_rxns, "r_0451", "r_0548", "r_0961")

# remove flawed reactions
optimal_rxn_form <- reactionInfo %>% filter(rMech %in% optimal_rxn_form) %>% filter(!(reaction %in% flawed_rxns)) %>% dplyr::select(rMech) %>% unlist() %>% unname()


### load manually annotated abbreviation of reaction name and pathway

fitReactionNames <- read.delim('companionFiles/fitReactionNames.txt')

if(!all(valid_rxns %in% fitReactionNames$reaction)){
 stop("some reactions that were tested need to be added to fitReactionNames.txt") 
}


#### Summary based on spearman correlation for MM and most significant regulator (if applicable) #####

bestModel <- reactionInfo[reactionInfo$rMech %in% optimal_rxn_form,] %>% dplyr::select(reaction, rMech, modelType)

complexReg <- bestModel %>% filter(!(modelType %in% c("rMM", "regulator"))) %>% dplyr::select(reaction) %>% unlist() %>% unname()
addedReg <- reactionInfo %>% filter(modelType == "regulator" & reaction %in% complexReg) %>% filter(!(reaction %in% rmCond_rxns & !grepl('rmCond$', rMech))) %>% group_by(reaction) %>%
  filter(ML == max(ML)) %>% dplyr::select(reaction, rMech, modelType)

allReg <- bestModel %>% filter(modelType != "rMM")%>% dplyr::select(reaction) %>% unlist() %>% unname()
addRMM <- reactionInfo %>% filter(modelType == "rMM" & reaction %in% allReg) %>% filter(!(reaction %in% rmCond_rxns & !grepl('rmCond$', rMech))) %>%
  dplyr::select(reaction, rMech, modelType)

# For complex reaction, also include rMM and single regulation
# For single regulation, include rMM
all_rxn_fits <- rbind(bestModel, addedReg, addRMM)

all_rxn_fits <- all_rxn_fits %>% left_join(rxn_fits %>% dplyr::select(rMech = rxn, spearman = parSpearman), by = "rMech") %>%
  mutate(modelType = ifelse(modelType %in% c("2+ regulators", "coopertivity"), "complex", modelType))

select_spearman_MMandReg <- all_rxn_fits %>% mutate(modelType = factor(modelType, levels = c("rMM", "regulator", "complex"))) %>% 
  group_by(reaction) %>% mutate(max_spearman = max(spearman)) %>% ungroup() %>% arrange(max_spearman, modelType) %>%
  mutate(spearman = ifelse(spearman < 0, 0, spearman))

select_spearman_MMandReg <- select_spearman_MMandReg %>% mutate(rMech = factor(rMech, levels = select_spearman_MMandReg$rMech))

barplot_theme_nox <- theme(text = element_text(size = 25), title = element_text(size = 25), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black"), axis.text.x = element_blank(),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )

ggplot(select_spearman_MMandReg, aes(x = rMech, y = spearman, fill = modelType)) + geom_bar(stat = "identity", color = "BLACK", width = 0.85) +
  barplot_theme_nox + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_fill_manual("", values = c("skyblue2", "orangered1", "green"), breaks = c("rMM", "regulator", "complex"),
                                                                            label = c("Reversible Michaelis-Menten", "Metabolite Regulator", "Complex Regulation")) +
  ggtitle('Correlation between measured and predicted flux') + expand_limits(y = c(0,1))
ggsave("Figures/MM_with_reg_spearman.pdf", height = 8, width = 12)

# Replot using a stacked barplot

regulation_cum_improvement <- all_rxn_fits %>% dplyr::select(-rMech) %>% spread(modelType, spearman) %>% 
  mutate(rMM = ifelse(rMM < 0, 0, rMM)) %>%
  mutate(regulator = ifelse(regulator < rMM, rMM, regulator)) %>%
  mutate(complex = ifelse(complex < regulator, regulator, complex)) %>%
  mutate(cum_regulator = regulator - rMM,
         cum_complex = complex - regulator)

regulation_cum_improvement <- regulation_cum_improvement %>% gather(modelType, spearman, -reaction) %>% filter(!is.na(spearman)) %>%
  filter(modelType %in% c("rMM", "cum_regulator", "cum_complex")) %>% mutate(modelType = factor(modelType, levels = c("rMM", "cum_regulator", "cum_complex"))) %>%
  left_join(fitReactionNames %>% dplyr::select(reaction, abbrev), by = "reaction")

# filter reactions that are never positively correlated
no_pos <- regulation_cum_improvement %>% group_by(reaction) %>% dplyr::summarize(no_pos = n() == sum(spearman <= 0))
regulation_cum_improvement <- regulation_cum_improvement %>% filter(!(reaction %in% no_pos$reaction[no_pos$no_pos]))

rMM_order <- regulation_cum_improvement %>% filter(modelType == "rMM") %>% arrange(spearman)

regulation_cum_improvement <- regulation_cum_improvement %>% mutate(reaction = factor(reaction, levels = rMM_order$reaction))

barplot_theme_withx <- barplot_theme_nox + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16))

ggplot(regulation_cum_improvement, aes(x = reaction, y = spearman, fill = modelType)) + geom_bar(stat = "identity", color = "black", width = 0.8) +
  ggtitle('Correlation between measured and predicted flux') +
  barplot_theme_withx + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.2,0.4, 0.6, 0.8,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0), breaks = regulation_cum_improvement$reaction, labels = regulation_cum_improvement$abbrev) +
  scale_fill_manual("", values = c("skyblue2", "gold1", "orangered1"), breaks = c("rMM", "cum_regulator", "cum_complex"),
                    label = c("Reversible Michaelis-Menten", "Metabolite Regulator", "Complex Regulation"))
ggsave("Figures/spearman_stack.pdf", height = 9, width = 14)

#### Comparison of reversible michaelis-menten and irreversible kinetics #####

spearman_reversibility <- data.frame(reactionInfo[,c('reaction', 'modification', 'Qvalue')], spearman = rxn_fits[,'parSpearman']) %>% tbl_df()  # all reactions
spearman_reversibility <- spearman_reversibility %>% dplyr::filter(is.na(Qvalue) | grepl('^(forward|reverse)', modification)) %>%
  dplyr::filter(reaction %in% valid_rxns)
spearman_reversibility <- spearman_reversibility %>% left_join(spearman_reversibility %>% group_by(reaction) %>% dplyr::summarize(rmCond_reaction = ifelse(sum(modification == 'rmCond') , T, F)))

rmCond_sets <- spearman_reversibility %>% ungroup() %>% filter(grepl('rmCond', modification)) %>% dplyr::select(reaction, modification) # Take the rmCond reactions when they exist and the normal reactions otherwise

spearman_reversibility <- rbind(spearman_reversibility %>% filter(rmCond_reaction == F),
spearman_reversibility %>% dplyr::filter(rmCond_reaction == T) %>% inner_join(rmCond_sets))

# spearman_reversibility %>% dplyr::filter(Qvalue < 0.1)

mm_form_counts <- spearman_reversibility %>% group_by(reaction) %>% summarize(count = n()) %>% dplyr::select(count) %>% unlist() %>% unname()
if(any(mm_form_counts > 2)){
 stop("choose either the forward or reverse irreversible rxn") 
}

spearman_reversibility <- spearman_reversibility %>% dplyr::mutate(Type = ifelse(grepl('^(forward|reverse)', modification), "Irreversible", "Reversible")) %>%
  group_by(reaction) %>% dplyr::mutate(maxSpear = max(spearman)) %>% ungroup() %>% dplyr::arrange(maxSpear) %>% dplyr::mutate(reaction = factor(reaction, levels = unique(reaction))) %>%
  dplyr::filter(maxSpear >= 0)

spearman_reversibility <- spearman_reversibility %>% group_by(reaction) %>% dplyr::mutate(revMax = if("Irreversible" %in% Type){
 if(spearman[Type == "Reversible"] >= spearman[Type == "Irreversible"]){T}else{F}
}else{T})
spearman_reversibility <- spearman_reversibility %>% dplyr::filter(!(Type == "Irreversible" & revMax == F))

# Determine improvement due to reversibility
rev_improvement <- spearman_reversibility %>% dplyr::filter(Type == "Irreversible") %>% rowwise() %>% dplyr::mutate(spearDiff = maxSpear - max(c(0, spearman))) %>% 
  dplyr::filter(spearDiff != 0) %>% dplyr::select(reaction, Type, spearman, spearDiff)
rev_fits <- spearman_reversibility %>% dplyr::filter(Type == "Reversible") %>% dplyr::select(reaction, Type, spearman)

stacked_spearman_rev <- rbind(
  rev_improvement %>% dplyr::select(-spearDiff), # lower correlation due to irreversibility
  rev_fits %>% dplyr::filter(reaction %in% rev_improvement$reaction) %>% dplyr::mutate(spearman = rev_improvement$spearDiff[rev_improvement$reaction == reaction]), # improved correlation due to reversibility
  rev_fits %>% dplyr::filter(!(reaction %in% rev_improvement$reaction)) %>% dplyr::mutate(Type = "Irreversible") # rev == irrev
)

stacked_spearman_rev <- stacked_spearman_rev %>% dplyr::filter(spearman >= 0)%>% dplyr::group_by(reaction) %>% dplyr::mutate(maxSpear = sum(spearman)) %>% 
  ungroup() %>% dplyr::arrange(maxSpear) %>%
  mutate(reaction = as.character(reaction)) %>% mutate(reaction = factor(reaction, levels = unique(reaction)))
stacked_spearman_rev$reaction = factor(stacked_spearman_rev$reaction, levels = unique(stacked_spearman_rev$reaction))

ggplot(stacked_spearman_rev, aes(x = reaction, y = spearman, fill = Type)) + geom_bar(stat = "identity", color = "black", width = 0.8) +
  ggtitle('Correlation between measured and predicted flux') +
  barplot_theme_withx + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0), breaks = fitReactionNames$reaction, labels = fitReactionNames$abbrev) + scale_fill_manual("", values = c("Reversible" = "skyblue2", "Irreversible" = "moccasin"), label = c("Irreversible Michaelis-Menten", "+ Reversibility"))

ggsave("Figures/spearman_stack_rev.pdf", height = 9, width = 14)

##### Replotting a few reactions for figures #####

if("param_set_list" %in% ls()){
allostery_affinity()  
}

##### Generate figure summarizing metabolic leverage for a condition #####

MLdata$Type[grepl("r_0302", MLdata$reaction) & MLdata$specie == "isocitrate"] <- "product" # aconitase is split into two reactions without a meaasured cis-aconitate so isocitrate acts like a product in the first rxn
MLdata <- data.table(MLdata %>% filter(conditions == "NATURAL"))
#metabolic_leverage_summary_plots("P0.05")

# Only looking at reactions that are well-fit
# Look at the best significant reaction form for reaction where CI overlap with flux carried is substantial (>50%)

#check_rxns <- setdiff(intersect(optimal_rxn_form, fraction_flux_deviation$rxn[fraction_flux_deviation$"Interval Overlap" > 0.5]), intersect(optimal_rxn_form, rxn_fits$rxn[rxn_fits$parSpearman > 0.6]))
#rxn_fits[rxn_fits$rxn %in% check_rxns,]
#fraction_flux_deviation[fraction_flux_deviation$rxn %in% check_rxns,]

#adequate_fit_optimal_rxn_form <- union(intersect(optimal_rxn_form, fraction_flux_deviation$rxn[fraction_flux_deviation$"Interval Overlap" > 0.5]), intersect(optimal_rxn_form, rxn_fits$rxn[rxn_fits$parSpearman > 0.6]))

# Save reactions where the optimal reaction form has a spearman correlation of > 0.6

adequate_fit_optimal_rxn_form <- intersect(optimal_rxn_form, rxn_fits$rxn[rxn_fits$parSpearman > 0.6])

adequate_rxn_form_data <- data.frame(rxnForm = adequate_fit_optimal_rxn_form) %>% mutate(reaction = substr(rxnForm, 1, 6)) %>% left_join(fitReactionNames, by = "reaction")

# check for cases where variable hill coefficient regulation may have been missing
# check arginosuccinate lyase |- Arginine (r_0207)
# glutamate dehydrogenase |- quinolinate (r_0471)

# significant_hill <- reactionInfo %>% filter(grepl('ultra', modification) & Qvalue < 0.1)

# summarize metabolic leverage

ML_inducibility_summary <- enzyme_control_source()

######## Split ML of well-fit reactions into enzyme and allosteric control ###########

# if any reactions flow in reverse, there substrate and product leverage needs to be flipped
reverse_rxns <- sapply(adequate_fit_optimal_rxn_form, function(x){all(rxnList_form[[x]]$flux$standardQP < 0)})
reverse_rxns <- names(reverse_rxns)[reverse_rxns]

ML_rxn_summary <- tbl_df(MLdata) %>% dplyr::filter(reaction %in% adequate_fit_optimal_rxn_form)

ML_rxn_summary <- rbind(
  ML_rxn_summary %>% filter(reaction %in% reverse_rxns, Type == "substrate") %>% mutate(Type = "product"),
  ML_rxn_summary %>% filter(reaction %in% reverse_rxns, Type == "product") %>% mutate(Type = "substrate"),
  ML_rxn_summary %>% filter(reaction %in% reverse_rxns, !(Type %in% c("substrate", "product"))),
  ML_rxn_summary %>% filter(!(reaction %in% reverse_rxns))
)

# Use metabolic leverage associated with parameter MLE
ML_rxn_summary <- ML_rxn_summary %>% mutate(Type = ifelse(Type %in% c("substrate", "product", "enzyme"), Type, "regulator")) %>%
  dplyr::select(reaction, specie, Type, condition, ML = MLE) %>% group_by(reaction, Type, condition) %>%
  dplyr::summarize(ML = sum(ML)) %>% group_by(reaction, Type) %>% dplyr::summarize(ML = median(ML)) %>%
  group_by(reaction) %>% dplyr::mutate(ML = ML / sum(ML))

ML_rxn_tall <- ML_rxn_summary %>% left_join(ML_inducibility_summary %>% dplyr::select(reaction, genes, reversibility)) %>%
  mutate(Type = factor(Type, levels = c("substrate", "product", "enzyme", "regulator")))

ML_rxn_summary <- tbl_df(dcast(ML_rxn_summary, reaction ~ Type, value.var = "ML", fill = 0)) # each rxn with assocaited enzyme, regulator and metabolic control fraction

ML_rxn_summary <- ML_rxn_summary %>% dplyr::mutate(rID = substr(reaction, 1,6)) %>% left_join(ML_inducibility_summary %>% dplyr::select(reaction, genes, reversibility)) %>%
  arrange(substrate + product)

ML_rxn_tall <- ML_rxn_tall %>% ungroup() %>% mutate(rxnForm = reaction, reaction = substr(reaction, 1, 6)) %>% 
  left_join(fitReactionNames %>% dplyr::select(reaction, abbrev)) %>%
  mutate(rxnForm = factor(rxnForm, levels = rev(ML_rxn_summary$reaction))) %>% arrange(rxnForm)

ML_rxn_tall <- ML_rxn_tall %>% mutate(reversibility = ifelse(reversibility == "T", "Kinetically Reversible", "Kinetically Irreversible"),
                                      reversibility = factor(reversibility, levels = c("Kinetically Reversible", "Kinetically Irreversible")))

TypeColors <- c("chartreuse3", "chartreuse", "deepskyblue", "orangered1")
  names(TypeColors) <- levels(ML_rxn_tall$Type)

barplot_facet_theme <- theme(text = element_text(size = 20, face = "plain"), title = element_text(size = 25), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black", size = 20),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank(),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)

ggplot(ML_rxn_tall, aes(x = rxnForm, y = ML, fill = Type, order = Type)) + geom_bar(stat = "identity", width = 0.85, color = "BLACK") +
  barplot_facet_theme + scale_y_continuous(expression('Metabolic Leverage: ' ~ frac("|"~epsilon[i]~"|"~sigma[i], sum("|"~epsilon[j]~"|"~sigma[j], "j = 1" , n))), expand = c(0,0)) +
  scale_fill_manual(values = TypeColors) + facet_grid(~ reversibility, scales = "free_x", space = "free_x") +
  scale_x_discrete("Reactions", breaks = ML_rxn_tall$rxnForm, labels =ML_rxn_tall$abbrev)
ggsave("Figures/metabolicLeverageBar.pdf", height = 10, width = 16)



color_key <- color_ternary(100)
#color_key <- color_simplex(200, 0, T, 6)

color_index <- mapply(function(x,y){
  which.min(abs(color_key$Table$enzyme - x) + abs(color_key$Table$allostery - y))
}, x = ML_rxn_summary$enzyme, y = ML_rxn_summary$regulator)

ML_rxn_summary <- ML_rxn_summary %>% mutate(rxnForm = reaction) %>% dplyr::select(-rID, -reaction, -genes) %>% left_join(adequate_rxn_form_data, by = "rxnForm")

# Look at two alternative ways of describing the 3 components forming the ternary space
# substrates & products, enzymes, regulators
# substrates, enzymes, regulators & products

# substrates & products, enzymes, regulators

ML_rxn_summary_1 <- ML_rxn_summary %>% cbind(color_key$Table %>% dplyr::slice(color_index) %>% dplyr::select(color)) %>%
  mutate(rxn_metabolite = substrate + product)

ML_rxn_ternaryPoints_1 <- ML_rxn_summary_1 %>%  mutate(x = (1/2)*(2*regulator + enzyme) / (regulator + enzyme + rxn_metabolite),
                                              y = sqrt(3)/2 * enzyme*(regulator + enzyme + rxn_metabolite))

summary(lm(ML_rxn_summary_1, formula = rxn_metabolite ~ reversibility))

color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(reversibility == "T"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(reversibility == "T"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
ggsave("Figures/MLcolorKey_REV.pdf", height = 9.1, width = 10.7)

# Irreversible
color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(reversibility == "F"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(reversibility == "F"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
  ggsave("Figures/MLcolorKey_FOR.pdf", height = 9.1, width = 10.7)

ML_rxn_ternaryPoints_labels <- ML_rxn_ternaryPoints_1 %>% filter(regulator != 0)

color_key$Figure_Color + 
  geom_point(data = ML_rxn_ternaryPoints_1, aes(x = x, y = y), size = 7, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1, aes(x = x, y = y), size = 6, shape = 21, fill = "WHITE") +
  geom_text(data = ML_rxn_ternaryPoints_labels, aes(x = x + 0.05, y = y + 0.05, label = abbrev), color = "WHITE")
ggsave("Figures/MLcolorKey.pdf", height = 9.1, width = 10.7)

control_layout <- ML_rxn_summary_1 %>% left_join(reactionInfo %>% dplyr::select(rxnForm = rMech, Name)) %>%
  dplyr::select(reaction, reaction.name, abbrev, rxnForm, form = Name, rxn_metabolite, enzyme, regulator, color)

write.table(control_layout, "flux_cache/control_layout.tsv", col.names = T, row.names = F, quote = F, sep = "\t")




######### Compare Km values found through optimization to those from literature ################

all_affinities <- read.delim("flux_cache/metaboliteAffinities.tsv") # all BRENDA substrates and regulators
all_affinities <- all_affinities[!is.na(all_affinities$log10mean),] # substrates and regulators with a ki or km

all_affinities <- lapply(1:nrow(all_affinities), function(i){
  data.frame(expand.grid(modelName = strsplit(all_affinities$tID[i], split = '/') %>% unlist(),
  rxn = strsplit(all_affinities$reactions[i], split = '/') %>% unlist()),
  all_affinities[i, !(colnames(all_affinities) %in% c('tID', 'reactions'))], row.names = NULL
  )
})
all_affinities <- do.call("rbind", all_affinities) %>% tbl_df() %>%
  mutate(modelName = as.character(modelName), rxn = as.character(rxn))

all_affinities <- all_affinities %>% filter(rxn %in% unique(reactionInfo$reaction))

substrate_affinities <- all_affinities %>% filter(speciesType == "substrate")
regulator_affinities <- all_affinities %>% filter(speciesType == "regulator")

### Using 95% confidence intervals for reaction parameters
yeast_only <- F
affinity_comparisons <- NULL

for(rxn in adequate_fit_optimal_rxn_form){
  
  rxData <- data.frame(rxn = reactionInfo$reaction[reactionInfo$rMech == rxn], rxnForm = rxn, 
                       cbind(rxn_fit_params[[rxn]]$kineticParPrior[rxn_fit_params[[rxn]]$kineticParPrior$SpeciesType == "Metabolite",] %>% dplyr::select(modelName, commonName, formulaName, measured),
                             rxn_fit_params[[rxn]]$param_interval[rxn_fit_params[[rxn]]$kineticParPrior$SpeciesType == "Metabolite",] %>% dplyr::select(parLB = get('X2.5.'), parUB = get('X97.5.'), parMLE = MLE, absoluteQuant),
                             rxn_fit_params[[rxn]]$param_species %>% dplyr::select(Subtype, medianAbund))
  )
  rxData <- rxData %>% mutate(absoluteQuant = as.logical(absoluteQuant))
  
  rxnData <- rbind(rxData %>% filter(Subtype %in% c("substrate", "product")) %>% left_join(substrate_affinities, by = c("rxn", "modelName")),
                   rxData %>% filter(!(Subtype %in% c("substrate", "product"))) %>% left_join(regulator_affinities, by = c("rxn", "modelName")))
  
  if(yeast_only){
    rxnData <- rxnData %>% filter(is.na(isYeast) | isYeast) 
  }else{
    rxnData <- rxnData %>% filter(is.na(isYeast) | !isYeast | (isYeast & formulaName %in% names(table(rxnData$formulaName))[unname(table(rxnData$formulaName)) == 1]))
  }
  
  if(length(unique(rxnData$EC[!is.na(rxnData$EC)])) != 1){
    # multiple EC numbers inform this reaction
    rxnData <- rxnData %>% group_by(formulaName) %>% filter(nQuant == max(nQuant) | is.na(nQuant))
  }
  
  affinity_comparisons <- rbind(affinity_comparisons, rxnData)
}

affinity_comparisons <- affinity_comparisons %>% tbl_df()

### compare substrate affinities for metabolites with absolute concentrations ###

absolute_comparison <- affinity_comparisons %>% filter(Subtype %in% c("substrate", "product"), absoluteQuant, !is.na(log10mean)) %>% mutate(log10mean = log10mean - 3) # convert from mM to M

boxplot_theme <- theme(text = element_text(size = 20), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black", size = 20),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )

absolute_comparison <- absolute_comparison %>% mutate(sdOflog10 = ifelse(is.na(sdOflog10), 0, sdOflog10), parCapture = log10(parLB) < log10mean & log10(parUB) > log10mean, parCaptureInterval = log10(parLB) < log10mean + 2*sdOflog10 & log10(parUB) > log10mean - 2*sdOflog10)

absolute_comparison$paramConsistency <- factor(mapply(x = absolute_comparison$parCapture, y = absolute_comparison$parCaptureInterval, function(x,y){
 if(x & y){
   "parameter interval contains BRENDA average"
 }else if(!x & y){
   "parameter interval consistent given BRENDA uncertainty"
 }else{
   "parameter estimate is inconsistent with BRENDA"
 }
}), levels = c("parameter interval contains BRENDA average",
               "parameter interval consistent given BRENDA uncertainty",
               "parameter estimate is inconsistent with BRENDA"))

ggplot(absolute_comparison, aes(x = log10mean, y = log10(parMLE), fill = paramConsistency, group = measured)) + geom_point(size = 6, shape = 21) + 
  geom_smooth(method = "lm", size = 2, color = "black", se = F) +
  boxplot_theme + scale_fill_brewer(palette = "Dark2") + scale_size_identity() +
  scale_x_continuous(expression("Literature consensus" ~ log[10] ~ "Affinity (M)")) +
  scale_y_continuous(expression("Inferred" ~ italic("in vivo") ~ log[10] ~ 'Affinity (M) 95% credicibility interval')) +
  guides(fill = guide_legend(nrow = 3))
ggsave("Figures/brendaConsistency.pdf", width = 8, height = 10)

summary(lm(log10(absolute_comparison$parMLE) ~ absolute_comparison$log10mean))
  

ggplot(absolute_comparison, aes(x = log10mean, y = log10(parMLE), ymin = log10(parLB), ymax = log10(parUB))) + geom_errorbar(aes(color = parCapture), size = 0.9, alpha = 0.7, width = 0.15) +
  geom_abline(a = 0, b = 1, size = 2) + boxplot_theme +
  scale_x_continuous(expression("Literature consensus" ~ log[10] ~ "Affinity (M)")) +
  scale_y_continuous(expression("Inferred" ~ italic("in vivo") ~ log[10] ~ "Affinity (M) 95% credicibility interval")) + coord_cartesian(ylim = c(-8,0)) +
  scale_color_manual(guide = "none", values = c("TRUE" = "chartreuse4", "FALSE" = "firebrick2"))
ggsave("Figures/brendaAffinity.pdf", width = 10, height = 10)
ggsave("Figures/brendaAffinity.eps", width = 10, height = 10, device=cairo_ps)

brendaAgree <- absolute_comparison %>% mutate(sdOflog10 = ifelse(is.na(sdOflog10), 0, sdOflog10), 
                                              parCapture = log10(parLB) < log10mean & log10(parUB) > log10mean, 
                                              parCaptureInterval = log10(parLB) < log10mean + 2*sdOflog10 & log10(parUB) > log10mean - 2*sdOflog10) %>%
  dplyr::select(rxn, rxnForm, modelName, nQuant, parCapture, parCaptureInterval)

brendaAgree$parCapture %>% table() %>% as.data.frame() %>% summarize(Freq[. == T] / sum(Freq))
brendaAgree$parCaptureInterval %>% table() %>% as.data.frame() %>% summarize(Freq[. == T] / sum(Freq))

### Look at metabolism-wide occupancy ###
# need to fully load up to *import cluster parameters* in this script

occupancy_comparison <- affinity_comparisons %>% filter(measured) %>% dplyr::select(-formulaName, -measured, -absoluteQuant, -(EC:isYeast), -(nQuant:speciesType))

# use concentrations of species across all conditions rather than just median

rxn_met_pairs <- occupancy_comparison %>% dplyr::select(rxnForm, modelName)

rxn_met_abundances <- lapply(1:nrow(rxn_met_pairs), function(i){
  metConc <- rxnList[[rxn_met_pairs[i,'rxnForm'] %>% unlist()]]$rxnMet[,rxn_met_pairs[i,'modelName'] %>% unlist()]
 cbind(rxn_met_pairs[i,'rxnForm'], rxn_met_pairs[i,'modelName'], logConc = metConc)
})
rxn_met_abundances <- do.call("rbind", rxn_met_abundances)

occupancy_comparison <- occupancy_comparison %>% left_join(rxn_met_abundances)

occupancy_comparison$Subtype[!(occupancy_comparison$Subtype %in% c("substrate", "product"))] <- "regulator"
occupancy_comparison$Specie <- NA

occupancy_comparison <- occupancy_comparison %>% mutate(log10S_km = log10(2^logConc / parMLE), log10occ = (2^logConc/(2^logConc + parMLE)))

occupancy_comparison <- occupancy_comparison %>% mutate(Subtype = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", Subtype, perl=TRUE))
occupancy_comparison <- occupancy_comparison %>% mutate(Subtype = factor(Subtype, levels = c("Substrate", "Product", "Regulator")))

# check whether distribution of affinities departs from null expectation

occupancy_dist_test <- occupancy_comparison %>% dplyr::select(Subtype, rxn, modelName, medianAbund, parMLE) %>% unique()
occupancy_dist_test <- occupancy_dist_test %>% mutate(QKeqDiff = medianAbund - log2(parMLE)) %>%
  group_by(Subtype) %>% summarize(KSp = ks.test(QKeqDiff, "punif", -15, 15)$p.value, counts = n())
occupancy_dist_test <- occupancy_dist_test %>% mutate(Sig = ifelse(KSp < 0.05, "*", ""), Sig = ifelse(KSp < 0.01, "**", Sig), Sig = ifelse(KSp < 0.001, "***", Sig)) %>%
  mutate(Sig = ifelse(Sig != "", paste0('\n',Sig), Sig), label = paste0(counts, Sig))

boxplot_theme_label_rotate <- boxplot_theme + theme(axis.title.y = element_text(angle = 0, size = 30), axis.title.x = element_blank(),
                                                    axis.text = element_text(size = 25))

ggplot() + geom_violin(data = occupancy_comparison, aes(x = Subtype, y = log10S_km), fill = "firebrick2", scale = "width") +
  geom_text(data = occupancy_dist_test, aes(x = Subtype, y = 0.5, label = label), color = "black", size = 8) +
  scale_y_continuous(expression(frac('[S]', K[M])), expand = c(0,0)) +
  boxplot_theme_label_rotate

ggplot() + geom_violin(data = occupancy_comparison, aes(x = Subtype, y = log10occ), fill = "darkgoldenrod1", scale = "width") +
  geom_text(data = occupancy_dist_test, aes(x = Subtype, y = 0.5, label = label), color = "black", size = 8) +
  scale_y_continuous(expression(frac('[S]', '[S]' + K[M])), breaks = seq(0, 1, by = 0.2), expand = c(0,0), limits = c(0,1)) +
  boxplot_theme_label_rotate
ggsave("Figures/speciesOccupancy.pdf", width = 8, height = 8)

#Iexp <- 2^(rnorm(100, 0, 1))
#Ivals <- 2^seq(-5,5,by = 0.01)

#Icomp <- expand.grid(conc = Iexp, affinity = Ivals) %>% tbl_df()
#Isummary <- Icomp %>% mutate(rate = affinity / (affinity + conc)) %>% group_by(affinity) %>%
#  summarize(average = mean(rate), rate_sd = sd(rate))

#ggplot2::qplot(x = Iexp) + scale_x_log10()
#ggplot(Isummary, aes(x = log2(affinity), y = average, ymin = average - 2*rate_sd, ymax = average + 2*rate_sd)) + geom_pointrange()
#ggplot(Isummary, aes(x = log2(affinity), y = rate_sd)) + geom_line()




######## Compare Keq to Gibbs free energy ########

# rxForm defined keq: prod[S] - prod[P]/keq
# calculate median absolute concentration or bounds based upon reasonable concentrations of small molecules
# relate optimized keq to these values
# Q = prod[P]/prod[S]
# Q/keq

free_energy_table <- adequate_rxn_form_data %>% left_join(reversibleRx %>% dplyr::select(reaction = rx, cc_dGr = CCdG, cc_dGrSD = CCdGsd, reversible = modelBound) %>%
                                                            mutate(reversible = ifelse(reversible == "reversible", "Kinetically reversible", "Kinetically irreversible")))
free_energy_table <- free_energy_table %>% mutate(pathway = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", pathway, perl=TRUE))

free_energy_table[,c("log_keq_LB")] <- log2(sapply(free_energy_table$rxnForm, function(x){rxn_fit_params[[x]]$param_interval$'X2.5.'[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("log_keq_UB")] <- log2(sapply(free_energy_table$rxnForm, function(x){rxn_fit_params[[x]]$param_interval$'X97.5.'[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("log_keq_MLE")] <- log2(sapply(free_energy_table$rxnForm, function(x){rxn_fit_params[[x]]$param_interval$MLE[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("logQmedian")] <- sapply(free_energy_table$rxnForm, function(x){as.numeric(rxn_fit_params[[x]]$param_interval$absoluteQuant[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"])})
free_energy_table <- free_energy_table %>% filter(!is.na(cc_dGr))

freebee_species <- c(water = "t_0399", proton = "t_0398", ammonium = "t_0233", diphosphate = "t_0332", co2 = "t_0249")
missing_species <- lapply(free_energy_table$rxnForm, function(x){
  unmeasuredSpec = rxn_fit_params[[x]]$kineticParPrior %>% filter(!is.na(measured) & !measured) %>% dplyr::select(SubstrateID = rel_spec, commonName)
  unmeasuredSpec = unmeasuredSpec %>% left_join(rxn_fit_params[[x]]$param_species %>% dplyr::select(SubstrateID, ReactionID, Subtype), by = "SubstrateID")
  unmeasuredSpec %>% filter(!(SubstrateID %in% freebee_species))
  })
missing_species <- do.call("rbind", missing_species)

# calculate Q for all conditions to compare relative to Keq

library(tidyr)

all_condition_Q <- lapply(1:nrow(free_energy_table), function(i){
  metStoi <- rxnList[[free_energy_table$rxnForm[i]]]$rxnStoi
  metConc <- rxnList[[free_energy_table$rxnForm[i]]]$rxnMet
  metConc[is.na(metConc)] <- 0
  metConc$condition <- rownames(metConc)
  
  metInfo <- metConc %>% gather(specie, conc, -condition) %>% mutate(specie = as.character(specie)) %>% left_join(data.frame(specie = names(metStoi), coef = unname(metStoi)), by = "specie")
  metInfo %>% tbl_df() %>% group_by(condition) %>% summarize(Q = sum(conc * coef)) %>% mutate(reaction = free_energy_table$reaction[i])
})
all_condition_Q <- do.call("rbind", all_condition_Q)

# keep track of Q relative to median Q
all_condition_Q <- free_energy_table %>% left_join(all_condition_Q, join = "reaction") %>% tbl_df() %>% dplyr::select(rxnForm, reaction, condition, Q, logQmedian)

# back to one estimate-per-reaction
free_energy_table <- free_energy_table %>% mutate(QkeqDiff_MLE = logQmedian - log_keq_MLE, QkeqDiff_UB = logQmedian - log_keq_LB, QkeqDiff_LB = logQmedian - log_keq_UB) %>% 
  dplyr::select(-c(log_keq_LB, log_keq_UB, log_keq_MLE, logQmedian))

# for reactions with always negative flux, flip Q-Keq
reverse_rxns <- free_energy_table %>% filter(QkeqDiff_LB > 0) %>% mutate(QkeqDiff_MLE = -1*QkeqDiff_MLE, QkeqDiff_UBtmp = QkeqDiff_UB, QkeqDiff_LBtmp = QkeqDiff_LB, QkeqDiff_LB = -1*QkeqDiff_UBtmp, QkeqDiff_UB = -1*QkeqDiff_LBtmp) %>%
  dplyr::select(-QkeqDiff_LBtmp, -QkeqDiff_UBtmp)
free_energy_table <- free_energy_table %>% filter(QkeqDiff_LB <= 0) %>% rbind(reverse_rxns)

free_energy_table$measuredProducts <- sapply(free_energy_table$reaction, function(x){
  # check missing substrates (because reaction flows backwards) or products otherwise
  missing_species %>% filter(ReactionID == x & Subtype == ifelse(x %in% reverse_rxns$reaction, "substrate", "product")) %>% nrow() != 0
})

# remove zero variance reactions - these will only exist if the substrate and product use the same measurement - e.g. G1P and G6P
zeroVarRxns <- all_condition_Q %>% group_by(reaction) %>% summarize(var = var(Q)) %>% filter(var == 0) %>% dplyr::select(reaction) %>% unlist() %>% as.character()
free_energy_table <- free_energy_table %>% filter(!(reaction %in% zeroVarRxns))
all_condition_Q <- all_condition_Q %>% filter(!(reaction %in% zeroVarRxns))

# order reactions according to free energy
free_energy_table <- free_energy_table %>% arrange(QkeqDiff_UB) %>% mutate(reaction = factor(reaction, levels = reaction))

# also flipping Q and logQmedian Q when keeping track of all conditions
all_conditions_rev <- all_condition_Q %>% filter(rxnForm %in% reverse_rxns$rxnForm) %>% mutate(Q = -1*Q, logQmedian = -1*logQmedian)
all_condition_Q <- all_condition_Q %>% filter(!(rxnForm %in% reverse_rxns$rxnForm)) %>% rbind(all_conditions_rev)

# comparing variable Q across conditions to the median condition where Q / keq was calculated
all_condition_Q <- all_condition_Q %>% mutate(reaction = factor(reaction, levels = free_energy_table$reaction)) %>% 
  left_join(free_energy_table %>% dplyr::select(rxnForm, reaction, reversible, genes, abbrev, pathway, measuredProducts, QkeqDiff_UB), join = "reaction") %>%
  mutate(Qkeq = (Q - logQmedian) + QkeqDiff_UB)

fast_C_lim <- all_condition_Q %>% filter(condition == "C0.30")

free_energy_theme <- boxplot_theme_label_rotate + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), strip.background = element_rect(fill = "coral"))

# Visualizing disequilibrium ratio

ggplot() + facet_grid(~ pathway, scale = "free_x", space = "free_x") +
  geom_violin(data = all_condition_Q, aes(x = reaction, y = 2^Qkeq, fill = reversible), scale = "width") +
  geom_errorbar(data = free_energy_table, aes(x = reaction, ymin = 2^QkeqDiff_LB, ymax = 2^QkeqDiff_UB), size = 1) +
  geom_point(data = fast_C_lim, aes(x = reaction, y = 2^Qkeq), fill = "chartreuse3", size = 4, shape = 21) + 
  geom_point(data = free_energy_table, aes(x = reaction, y = 2^QkeqDiff_MLE), fill = "darkblue", size = 4, shape = 21) +
  scale_y_log10(expression(frac(Q, K[eq]) ~ "=" ~ frac(v[r], v[f])), expand = c(0,0)) + coord_cartesian(ylim = c(10^-3, 1)) + free_energy_theme +
  scale_fill_discrete() +
  scale_x_discrete(breaks = free_energy_table$reaction, labels = free_energy_table$abbrev)
ggsave("Figures/reactionDisequilibrium.eps", width = 15, height = 8)
ggsave("Figures/reactionDisequilibrium.pdf", width = 15, height = 8)




