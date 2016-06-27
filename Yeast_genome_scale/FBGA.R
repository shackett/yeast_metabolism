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
library(tidyr)

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

add_pairwise_regulation <- T
add_all_single_metabolites <- T

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
  if (valid_rxns[rxN] %in% c("r_0096", "r_0097", "r_0352", "r_669",
                             "r_0115", "r_0759", "r_0118", "r_0818", "r_0816", "r_0208", "r_0207")){
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
    rxnList[[nEntry]]$rxnMet <- rxnList[[nEntry]]$rxnMet[!rownames(rxnList[[nEntry]]$rxnMet) %in% conds, colnames(rxnList[[nEntry]]$rxnMet), drop = F]
    rxnList[[nEntry]]$all_species_SD <- rxnList[[nEntry]]$all_species_SD[!rownames(rxnList[[nEntry]]$all_species_SD) %in% conds, colnames(rxnList[[nEntry]]$all_species_SD), drop = F]
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

all_reactionInfo <- reactionInfo %>%
  dplyr::select(-form) %>%
  left_join(reaction_signif %>%
              dplyr::select(-signifCode), by = c("reaction", "rMech"))

# remove some reaction forms that we don't want in the final analysis
all_reactionInfo <- all_reactionInfo %>%
  filter(!grepl('inh-noncomp|inh-comp', modification)) %>%
  filter(!modelType %in% c("hypo met regulator", "hypo met cooperativity", "irreversible"))

# add the top non-literature activator and inhibitor of each reaction
top_non_lit <- all_reactionInfo %>%
  filter(modelType == "non-literature supported met regulator") %>%
  mutate(type = ifelse(grepl('act', rMech), "act", "inh")) %>%
  group_by(reaction, ncond, type) %>%
  filter(Qvalue == min(Qvalue)) %>%
  filter(Qvalue < 0.1)

reactionInfo <- all_reactionInfo %>%
  filter(modelType %in% c("rMM", "regulator") |
           (modelType != "non-literature supported met regulator" & Qvalue < 0.1) |
           rMech %in% top_non_lit$rMech)

#save(list = c("reactionInfo", "all_reactionInfo", "param_set_list", "parSetInfo", "rxnList_form", "param_run_info", "metSVD", "tab_boer"), file = "flux_cache/modelComparison.Rdata")

load("flux_cache/modelComparison.Rdata")

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

# load("companionFiles/PTcomparison_list.Rdata") # by-gene comparisons of protein and transcript abundance

shiny_flux_data <- list()
rxn_fits <- NULL
rxn_fit_params <- list()
fraction_flux_deviation <- NULL
MLdata <- list() # Save summary of metabolic leverage
ELdata <- list() # Save full distribution of elasticities
Hypo_met_candidates <- NULL

t_start = proc.time()[3]

# flag a few reactions where additional plots are created
custom_plotted_rxns <- c("r_1054-rm","r_0962-rm", "r_0962-rm-t_0290-act-mm",
                         "r_0514-rm", "r_0916-rm", "r_0208-rm",
                         "r_0915-rm", "r_0915-rm-t_0234-inh-uncomp", # PRPP Amidotransferase
                         "r_0468-rm", "r_0468-rm-t_0495-inh-uncomp", # G5K
                         "r_0309-rm", "r_0309-rm-t_0652-act-mm", #CBS
                         "r_0962-rm", "r_0962-rm-t_0276-inh-uncomp", # Citrate -| Pyk
                         "r_0962-rm-pairwise-t_0276-inh-uncomp+t_0290-act-mm", # citrate + FBP regulating Pyk
                         "r_0959-rm", "r_0959-rm-t_0457-inh-uncomp", # PhePyr -| PDC
                         "r_0816-rm_rmCond", "r_0816-rm-t_0461-inh-uncomp_rmCond","r_0816-rm", "r_0816-rm-t_0461-inh-uncomp"# Ala -| OTCase
                         ) 

custom_plotted_rxns <- unique(custom_plotted_rxns)
if(!(all(custom_plotted_rxns %in% reactionInfo$rMech))){stop("Some reaction mechanisms that you want to plot were not located")}

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
  #run_rxn$specSD[,'t_0234'] <- 0.4567561
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
  
  MLdata[[arxn]] <- reaction_properties$ML_summary
  ELdata[[arxn]] <- reaction_properties$EL_summary
  
  reaction_plots <- reaction_properties$plots
  shiny_flux_data[[arxn]]$plotChoices <- append(shiny_flux_data[[arxn]]$plotChoices,
                                                reactionPropertiesPlots(reaction_plots))
  
  if(arxn %in% custom_plotted_rxns){
    
    if(arxn == "r_0962-rm-pairwise-t_0276-inh-noncomp+t_0290-act-mm"){
      # overwrite relative measurement with absolute value for N0.11. Only affects visualization.
      run_rxn$rxnSummary$rxnMet$"t_0276" <- 
        run_rxn$rxnSummary$rxnMet$"t_0276" + (log2(0.00812) - run_rxn$rxnSummary$rxnMet$"t_0276"[rownames(run_rxn$rxnSummary$rxnMet) == "N0.11"])
      write.table(run_rxn$rxnSummary$rxnMet, file = "Figures/FBPregulators.tsv", sep = "\t", col.names = T, row.names = T, quote = F)
    
      run_rxn$rxnSummary$originMet[names(run_rxn$rxnSummary$originMet) == "t_0276"] <- "abs"
      }
    
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
  
  gc()
  
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

save(rxn_fit_params, rxn_fits, fraction_flux_deviation, Hypo_met_candidates, file = "flux_cache/paramCI.Rdata")

ELdata <- lapply(1:length(ELdata), function(x){
  ELdata[[x]] %>% mutate(rMech = names(ELdata)[x])
}) %>% bind_rows %>% tbl_df
saveRDS(ELdata, file = "flux_cache/elasticityData.Rds")

MLdata <- MLdata %>% bind_rows
saveRDS(MLdata, file = "flux_cache/MLdata.Rds")

##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@###@###@###@###@###@###@###@
###################### * Summarize Reaction Behavior ##########################
##@##@ Start here if loading parameter estimates and kinetic summaries ##@##@##
##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@###@###@###@###@###@###@###@

##### Systems level comparison of optimized and external parameter values #####

load("flux_cache/modelComparison.Rdata")
load("flux_cache/paramCI.Rdata")
load("flux_cache/reconstructionWithCustom.Rdata")
reversibleRx <- read.table("companionFiles/reversibleRx.tsv", header = T)

# Test if there were any reactions that were discarded on the fly

if(sum(!(reactionInfo$rMech %in% rxn_fits$rxn)) != 0){
  reactionInfo <- reactionInfo %>% dplyr::filter(rMech %in% rxn_fits$rxn)
  }

# exclude regulation that does not improve fit
reactionInfo <- reactionInfo %>%
  filter(modelType == "rMM" | Qvalue < 0.1)

### Determine which reactions to follow-up on based upon either
### Either all substrates are measured, or only relatively unimportant species are missing

load("paramOptim.Rdata")

RMMrxns <- reactionInfo %>% filter(modelType == "rMM", ncond == max(ncond)) %>% dplyr::select(rMech) %>% unlist() %>% unname()

reaction_validity <- filter_reactions(reactionInfo, rxn_fits, rxnList_form)
valid_rxns <- reaction_validity$reaction[reaction_validity$include]

# load manually annotated abbreviation of reaction name and pathway (for display purposes)

fitReactionNames <- read.delim('companionFiles/fitReactionNames.txt')

if(!all(valid_rxns %in% fitReactionNames$reaction)){
 stop("some reactions that were tested need to be added to fitReactionNames.txt") 
}

#### How literature support affect significance ####

rxn_regulation <- regulation_lit_support(valid_rxns, all_reactionInfo)
literature_support <- rxn_regulation[['prob_reg']]

#### Misc summaries of literature support (can skip this section) ####

# total citations
rxn_regulation[["prob_reg"]] %>% ungroup() %>% dplyr::summarize(sce = sum(sce), other = sum(other)) # all citations
table(literature_support$reaction[literature_support$sce != 0]) %>% length() # # reaction with yeast annotations
table(literature_support$reaction) %>% length() # # reactions with some BRENDA support

# increased p(improve fit) for gold-standard and BRENDA regulation
set.seed(1234)
chisq.test(rxn_regulation[['GS_contingency']], simulate.p.value = T, B = 1e7)

# What fraction of gold-standard, BRENDA and all metabolites improve fit

GS_reactions <- unique(literature_support$reaction[literature_support$is_GS])

met_type_fit_improvement <- all_reactionInfo %>% tbl_df() %>%
  filter(modelType == "non-literature supported met regulator" | (rMech %in% literature_support$rMech)) %>%
  mutate(modelType = ifelse(rMech %in% literature_support$rMech[literature_support$is_GS], "gold-standard", modelType)) %>%
  group_by(reaction) %>%
  filter(ncond == min(ncond)) %>% ungroup() %>%
  mutate(tID = substr(modification, 1, 6), reg_type = substr(modification, 8, 10)) %>%
  filter(reaction %in% GS_reactions)

library(qvalue)

met_type_fit_improvement$Qvalue_combined <- qvalue(met_type_fit_improvement$Pvalue)$q
  
FDR_category_improve <- met_type_fit_improvement %>%
  mutate(is_sig = Qvalue < 0.1) %>%
  count(modelType, is_sig) %>%
  spread(key = is_sig, value = n) %>%
  mutate(frac = `TRUE` / (`TRUE` + `FALSE`))
  
matrix_FDR_category_improve <- as.matrix(FDR_category_improve %>% ungroup %>% select(`FALSE`, `TRUE`))

comparision_matrix_FDR_category_improve <- rbind(matrix_FDR_category_improve[1,], colSums(matrix_FDR_category_improve[-1,]))
chisq.test(comparision_matrix_FDR_category_improve, simulate.p.value = T, B = 1e7)

# numbers "sig" per reaciton

sum(FDR_category_improve$`TRUE`)/length(GS_reactions)

# p(model | literature)
summary(rxn_regulation[['fitted_model']])

# look at the overall effect of this prior

# all tested regulation
all_reactionInfo %>%
  filter(reaction %in% valid_rxns) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup %>%
  mutate(modelType = ifelse(rMech %in% literature_support$rMech[literature_support$is_GS], "gold-standard", modelType)) %>%
  count(modelType)

# improves fit
all_reactionInfo %>%
  filter(reaction %in% valid_rxns) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup %>%
  filter(modelType %in% c("regulator", "non-literature supported met regulator")) %>%
  mutate(modelType = ifelse(rMech %in% literature_support$rMech[literature_support$is_GS], "gold-standard", modelType)) %>%
  mutate(is_sig = Qvalue < 0.1) %>%
  count(is_sig, modelType)

# top overall
all_reactionInfo %>%
  filter(reaction %in% valid_rxns) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup %>%
  filter(is.na(Qvalue) | Qvalue < 0.1) %>% 
  filter(modelType %in% c("rMM", "regulator", "non-literature supported met regulator")) %>%
  mutate(modelType = ifelse(rMech %in% literature_support$rMech[literature_support$is_GS], "gold-standard", modelType)) %>%
  group_by(reaction) %>% filter(ML == max(ML)) %>% slice(1) %>%
  count(modelType)

# excluding non-literature  
all_reactionInfo %>%
  filter(reaction %in% valid_rxns) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup %>%
  filter(is.na(Qvalue) | Qvalue < 0.1) %>% 
  filter(modelType %in% c("rMM", "regulator")) %>%
  mutate(modelType = ifelse(rMech %in% literature_support$rMech[literature_support$is_GS], "gold-standard", modelType)) %>%
  group_by(reaction) %>% filter(ML == max(ML)) %>% slice(1) %>%
  count(modelType)
  
#### Use priors to reevaluate tested regulation ####
# We first determined which mechanisms meaningfully improve simple mechanisms
# Now, for each reaction use empirical priors to prioritize models and filter those models which are very unlikely
# Group models which are similar together (correlated metabolites & similar interactions with other clusters)

rMech_support <- filter_rMech_by_prior(valid_rxns, reactionInfo %>% filter(modelType != "non-literature supported met regulator"),
                                       literature_support, rxnList_form)

# For each rMech containing a regulator, specify whether the regulation is feed-back, feed-forward or cross-pathway

rMech_mode <- mode_of_regulation(rMech_support, rxnList_form)

rMech_support <- rMech_support %>% left_join(rMech_mode, by = "rMech") %>%
  left_join(rxn_fits %>% dplyr::select(rMech = rxn, pearson = parPearson, spearman = parSpearman), by = "rMech")

# compare best predicted regulation with AIC

top_candidates <- all_reactionInfo %>%
  filter(reaction %in% valid_rxns) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>%
  filter(modelType %in% c("rMM", "regulator")) %>%
  mutate(nreg = sum(modelType != "rMM")) %>%
  filter(is.na(Qvalue) | Qvalue < 0.1) %>% 
  mutate(modelType = ifelse(rMech %in% literature_support$rMech[literature_support$is_GS], "gold-standard", modelType)) %>%
  left_join(rMech_support %>% ungroup %>% select(rMech, AICc, pearson), by = "rMech") %>%
  mutate(less_than_10 = ifelse(nreg <= 10, T, F))

top_candidate_ML <- top_candidates %>% group_by(reaction) %>% filter(ML == max(ML))
top_candidate_AIC <- top_candidates %>% group_by(reaction) %>% filter(AICc == min(AICc)) %>% count(modelType)


#rMech_support %>% group_by(type) %>% dplyr::summarize(AIC_prob = sum(AIC_prob))

#### Look at the best supported reaction mechanisms for each reaction ####

optimal_reaction_form <- rMech_support %>%
  group_by(reaction) %>%
  filter(ncond == min(ncond)) %>%
  filter(if(spearman[type == "unregulated"] > 0.9){
    type == "unregulated"
    }else{
      AIC_prob == max(AIC_prob)
      })

# use cached optimal reaction forms if desired (and inputs have not changed)
# to deal with a couple of reactions that fluctuate due to non-determinism
#overwrite_with_cache <- T

#if(overwrite_with_cache){
  
#  current_rxn_form <- read.delim("flux_cache/all_reaction_table.tsv") %>% select(rMech) %>% mutate(reaction = substr(rMech, 1, 6))
#  inconsistent_rxns <- current_rxn_form %>% filter(!(rMech %in% optimal_reaction_form$rMech))
  
#  }
#}


#### Compare gold-standard regulation to significant regulation ####

GS_regulation <- read.delim("companionFiles/gold_standard_regulation.txt") %>%
    tbl_df() %>% mutate(modtype = ifelse(type == "inhibitor", "inh", type),
                        modtype = ifelse(type == "activator", "act", modtype)) %>%
  filter(reaction %in% valid_rxns)
GS_regulation <- GS_regulation %>% left_join(literature_support %>% dplyr::select(reaction, tID, modtype, is_sig), by = c("reaction", "tID", "modtype"))

GS_regulation <- GS_regulation %>% left_join(optimal_reaction_form %>% dplyr::select(reaction, rMech, reg_type = type, ncond, spearman), by = "reaction")

GS_regulation_support <- GS_regulation %>% rowwise() %>% mutate(optimal = grepl(paste(tID, modtype, sep = "-"), rMech)) %>%
  mutate(support = ifelse(optimal, "Strongest support", ""),
         support = ifelse(!is.na(is_sig) & is_sig & !optimal, "Supported", support),
         support = ifelse(spearman < 0.6, "No tested model is adequate", support),
         support = ifelse(support == "" & reg_type == "unregulated", "No regulation identified", support),
         support = ifelse(is.na(is_sig), "Not measured", support),
         support = ifelse(support == "", "Alternative regulation", support))

# Summary of gold-standard regulation significance

library(xtable)
print(xtable(GS_regulation_support %>% dplyr::select(reaction_name, regulator, type, support) %>% ungroup() %>% arrange(tolower(reaction_name), regulator)),
             include.rownames = F)

### Include (in)validated regulation ###

optimal_reaction_form <- rbind(optimal_reaction_form %>% filter(!(reaction %in% c("r_0466", "r_0816"))),
                               rMech_support %>% filter(rMech %in% c("r_0816-rm-t_0461-inh-uncomp_rmCond", "r_0466-rm-pairwise-t_0617-inh-uncomp+t_0234-inh-uncomp"))) %>%
  ungroup() %>% arrange(reaction)

append_GS <- data.frame(reaction_name = 'ornithine transcarbamylase', reaction = 'r_0816', regulator = 'alanine',
        tID = 't_0461', type = 'inhibitor', reference = NA, modtype = 'inh', is_sig = T, rMech = 'r_0816-rm-t_0461-inh-comp_rmCond',
        reg_type = NA, ncond = 20, spearman = NA, optimal = NA, support = NA)
GS_regulation_support <- rbind(GS_regulation_support, append_GS)

append_GS <- data.frame(reaction_name = 'glucose 6-phosphate dehydrogenase', reaction = 'r_0466', regulator = c("phosphoenolpyruvate", "AMP"),
        tID = c("t_0617", "t_0234"), type = 'inhibitor', reference = NA, modtype = 'inh', is_sig = F, rMech = "r_0466-rm-pairwise-t_0617-inh-uncomp+t_0234-inh-uncomp",
        reg_type = NA, ncond = 25, spearman = NA, optimal = NA, support = NA)
GS_regulation_support <- rbind(GS_regulation_support, append_GS)


#### Group reaction mechanisms into similarly behaving groups and filter low-confidence mechanisms ####
# save this for later once well-fit reactions are identified

reg_groups <- group_regulation(rMech_support, rxnList_form, GS_regulation_support)

#### Look at well fit reactions ####

adequate_fit_optimal_rxn_form <- optimal_reaction_form$rMech[optimal_reaction_form$pearson > 0.6]
  
adequate_rxn_form_data <- optimal_reaction_form %>%
  filter(rMech %in% adequate_fit_optimal_rxn_form) %>%
  left_join(fitReactionNames, by = "reaction")

total_tested_regulators <- all_reactionInfo %>% tbl_df() %>% filter(modelType == "regulator") %>% dplyr::select(reaction, modification) %>%
  tidyr::separate(modification, into = c("tID", "class", "subclass"), sep = "-") %>%
  dplyr::select(-subclass) %>% unique() %>% group_by(reaction) %>%
  dplyr::summarize(N_reg_tested = n())

regulation_summary_table <- adequate_rxn_form_data %>%
  left_join(reg_groups$regulation, by = c("reaction", "ncond")) %>%
  left_join(reactionInfo %>% dplyr::select(rMech, modelType, Name), by = "rMech") %>%
  left_join(total_tested_regulators, by = "reaction") %>%
  dplyr::select(`Reaction Name` = reaction.name, Abbreviation = abbrev, `Best Supported Regulation` = Name,
                Correlation = pearson, `Other Supported Regulation` = all_regulators,
                `\\# Regulators Tested` = N_reg_tested) %>%
  mutate(Correlation = round(Correlation, 2)) %>%
  arrange(`Reaction Name`)

write.table(regulation_summary_table, file = "Figures/all_reaction_table.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

#### Summarize all "best regulators" ####

sig_regulators <- reactionInfo %>% filter(rMech %in% adequate_fit_optimal_rxn_form) %>% filter(modelType != "rMM")

sig_regulators <- do.call("rbind", lapply(1:nrow(sig_regulators), function(i){data.frame(reaction = sig_regulators$reaction[i], rMech = sig_regulators$rMech[i], regulators = strsplit(sig_regulators$modification[i], split = "\\+")[[1]])})) %>%
  mutate(regulators = gsub('(pairwise-)|(_ultra)|(_rmCond)', '', regulators))

ranked_regulators <- literature_support %>% group_by(reaction) %>% arrange(desc(prob_sig)) %>% mutate(prior_rank = paste(1:n(), "of", n())) %>%
  ungroup() %>% dplyr::select(reaction, tID, modtype, is_GS, sce, other, prior_rank, FullName)

sig_regulators <- sig_regulators %>% separate(regulators, into = c("tID", "modtype", "subtype"), sep = "-") %>%
  left_join(ranked_regulators, by = c("reaction", "tID", "modtype"))

write.table(sig_regulators %>% arrange(desc(is_GS), desc(sce), desc(other)), file = "Figures/significant_regulation.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

#### Table with all regulators ####

all_tested_regulators <- literature_support %>%
  ungroup %>%
  select(rMech, reaction, tID, modtype, sce, other) %>%
  left_join(fitReactionNames %>% select(reaction, reaction.name), by = "reaction") %>%
  left_join(tab_boer %>% select(SpeciesType, SpeciesName), by = c("tID" = "SpeciesType")) %>%
  left_join(rMech_support %>% ungroup %>% select(rMech, AICc, AIC_prob), by = "rMech") %>%
  select(Reaction = reaction.name, Regulator = SpeciesName, RegType = modtype, sce, other, AICc, AICprob = AIC_prob) %>%
  arrange(Reaction, Regulator, RegType)

write.table(all_tested_regulators, file = "Figures/all_tested_regulators.tsv", row.names = F, col.names = T, quote = F, sep = "\t")



#### Summary based on spearman correlation for MM and most significant regulator (if applicable) #####

bestModel <- reactionInfo %>% filter(rMech %in% optimal_reaction_form$rMech) %>% dplyr::select(reaction, rMech, modelType)

complexReg <- bestModel %>% filter(!(modelType %in% c("rMM", "regulator"))) %>% dplyr::select(reaction) %>% unlist() %>% unname()
addedReg <- reactionInfo %>% filter(modelType == "regulator" & reaction %in% complexReg) %>% group_by(reaction) %>% filter(ncond == min(ncond)) %>%
  filter(ML == max(ML)) %>% dplyr::select(reaction, rMech, modelType)

allReg <- bestModel %>% filter(modelType != "rMM")%>% dplyr::select(reaction) %>% unlist() %>% unname()
addRMM <- reactionInfo %>% filter(modelType == "rMM" & reaction %in% allReg) %>% group_by(reaction) %>% filter(ncond == min(ncond)) %>%
  dplyr::select(reaction, rMech, modelType)

# For complex reaction, also include rMM and single regulation
# For single regulation, include rMM
all_rxn_fits <- rbind(bestModel, addedReg, addRMM)
all_rxn_fits <- all_rxn_fits %>% left_join(rxn_fits %>% dplyr::select(rMech = rxn, pearson = parPearson), by = "rMech")

# Replot using a stacked barplot

regulation_cum_improvement <- all_rxn_fits %>% dplyr::select(-rMech) %>% spread(modelType, pearson) %>% 
  mutate(rMM = ifelse(rMM < 0, 0, rMM)) %>%
  mutate(regulator = ifelse(regulator < rMM, rMM, regulator)) %>%
  mutate(cooperativity = ifelse(cooperativity < regulator, regulator, cooperativity)) %>%
  mutate(`2+ regulators` = ifelse(`2+ regulators` < regulator, regulator, `2+ regulators`)) %>%
  mutate(cum_regulator = regulator - rMM,
         cum_cooperativity = cooperativity - regulator,
         `cum_2+ regulators` = `2+ regulators` - regulator)

regulation_cum_improvement <- regulation_cum_improvement %>% gather(modelType, pearson, -reaction) %>% filter(!is.na(pearson)) %>%
  filter(modelType %in% c("rMM", "cum_regulator", "cum_cooperativity", "cum_2+ regulators")) %>% mutate(modelType = factor(modelType, levels = c("rMM", "cum_regulator", "cum_cooperativity", "cum_2+ regulators"))) %>%
  left_join(fitReactionNames %>% dplyr::select(reaction, abbrev), by = "reaction")

# filter reactions that are never positively correlated
no_pos <- regulation_cum_improvement %>% group_by(reaction) %>% dplyr::summarize(no_pos = n() == sum(pearson <= 0))
regulation_cum_improvement <- regulation_cum_improvement %>% filter(!(reaction %in% no_pos$reaction[no_pos$no_pos]))

rMM_order <- regulation_cum_improvement %>% filter(modelType == "rMM") %>% arrange(pearson)

regulation_cum_improvement <- regulation_cum_improvement %>% mutate(reaction = factor(reaction, levels = rMM_order$reaction))

barplot_theme_nox <- theme(text = element_text(size = 25), title = element_text(size = 25), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black"), axis.text.x = element_blank(),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )
barplot_theme_withx <- barplot_theme_nox + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16))

ggplot(regulation_cum_improvement, aes(x = reaction, y = pearson, fill = modelType)) + geom_bar(stat = "identity", color = "black", width = 0.8) +
  ggtitle('Correlation between measured and predicted flux') +
  barplot_theme_withx + scale_y_continuous(name = "Pearson correlation", expand = c(0,0), breaks = c(0,0.2,0.4, 0.6, 0.8,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0), breaks = regulation_cum_improvement$reaction, labels = regulation_cum_improvement$abbrev) +
  scale_fill_manual("", values = c("skyblue2", "gold1", "mediumorchid2", "orangered1"), breaks = c("rMM", "cum_regulator", "cum_cooperativity", "cum_2+ regulators"),
                    label = c("generalized Michaelis-Menten", "+ metabolite regulator", "with cooperativity", "+ second regulator"))
ggsave("Figures/Pearson_stack.pdf", height = 9, width = 14)


##### Replotting a few reactions for figures #####

#if("param_set_list" %in% ls()){
#allostery_affinity()  
#}


#### Split ML of well-fit reactions into enzyme and allosteric control ####

MLdata <- readRDS("flux_cache/MLdata.Rds")

ML_rxn_summary <- MLdata %>%
  filter(measure == "weighted_sensitivities_arith",
         conditions == "Natural") %>%
  left_join(reactionInfo %>% select(rMech, reaction), by = "rMech")
ML_rxn_summary$Type[grepl("r_0302", ML_rxn_summary$reaction) & ML_rxn_summary$specie == "isocitrate"] <- "product" # aconitase is split into two reactions without a meaasured cis-aconitate so isocitrate acts like a product in the first rxn

# all reactions
# either "all" or "supported"
rxn_subset <- "supported"

if(rxn_subset == "all"){
  
  adequate_fit_optimal_rxn_form[substr(adequate_fit_optimal_rxn_form, 1, 6) == "r_0962"] <- "r_0962-rm-pairwise-t_0276-inh-uncomp+t_0290-act-mm"
  ML_rxn_summary <- ML_rxn_summary %>% filter(rMech %in% adequate_fit_optimal_rxn_form)
  adequate_rxn_form_data <- rMech_support %>% filter(rMech %in% adequate_fit_optimal_rxn_form) %>% left_join(fitReactionNames, by = "reaction")
  
}else if(rxn_subset == "supported"){
  
  # verified regulation + significant rMM
  adequate_fit_optimal_rxn_form <- read.table("companionFiles/supported_rMechs.txt", header = T)$rMech
  if(!all(adequate_fit_optimal_rxn_form %in% ML_rxn_summary$rMech)){
    stop("some rMechs in supported_rMechs.txt are not defined")
  }
  ML_rxn_summary <- ML_rxn_summary %>%
    filter(rMech %in% adequate_fit_optimal_rxn_form)
  adequate_rxn_form_data <- rMech_support %>% filter(rMech %in% adequate_fit_optimal_rxn_form) %>% left_join(fitReactionNames, by = "reaction")
  
}else{
 stop("not a valid input for ''rxn_subset''") 
}


# if any reactions flow in reverse, their substrate and product leverage needs to be flipped
reverse_rxns <- sapply(adequate_fit_optimal_rxn_form, function(x){all(rxnList_form[[x]]$flux$standardQP < 0)})
reverse_rxns <- names(reverse_rxns)[reverse_rxns]

ML_rxn_summary <- rbind(
  ML_rxn_summary %>% filter(reaction %in% reverse_rxns, Type == "substrate") %>% mutate(Type = "product"),
  ML_rxn_summary %>% filter(reaction %in% reverse_rxns, Type == "product") %>% mutate(Type = "substrate"),
  ML_rxn_summary %>% filter(reaction %in% reverse_rxns, !(Type %in% c("substrate", "product"))),
  ML_rxn_summary %>% filter(!(reaction %in% reverse_rxns))
)

# Use metabolic leverage associated with parameter MAP estimate
ML_rxn_summary <- ML_rxn_summary %>%
  mutate(Type = ifelse(Type %in% c("substrate", "product", "enzyme"), Type, "regulator")) %>%
  dplyr::select(reaction, rMech, specie, name, Type, ML = metabolicLeverage, markovSample, logLik)
  
# summary table of individual species

round3 = function(x){round(x, 3)}

ML_uncertainty <- ML_rxn_summary %>% group_by(reaction, rMech, name, Type) %>%
  summarize(MAP = ML[which.max(logLik)],
            LB = quantile(ML, 0.025),
            UB = quantile(ML, 0.975),
            IQR = quantile(ML, 0.75) - quantile(ML, 0.25)) %>%
  mutate_each(funs(round3), c(MAP, LB, UB, IQR)) %>%
  left_join(fitReactionNames %>% dplyr::select(reaction, reaction.name, abbrev), by = "reaction") %>%
  ungroup %>%
  dplyr::select(reaction = reaction.name, abbrev, specie = name, type = Type, MAP:IQR) %>%
  mutate(type = factor(type, levels = c("substrate", "product", "enzyme", "regulator"))) %>%
  arrange(reaction, type)
write.table(ML_uncertainty, file = "Figures/MLsummary.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

# summarize the reversibility of individual reactions

reaction_info_with_reversibility <- adequate_rxn_form_data %>% dplyr::select(reaction, rMech, reaction.name, genes) %>%
  left_join(fitReactionNames %>% dplyr::select(reaction, abbrev), by = "reaction") %>%
  left_join(reversibleRx %>%
              mutate(reversibility = ifelse(modelBound == "reversible", T, F),
              thermo_irr = ifelse(CCdG + 1.96*CCdGsd < -5, T, F),
              reversibility = ifelse(reversibility == TRUE, "Kinetically Reversible", "Kinetically Irreversible"),
              reversibility = ifelse(reversibility == "Kinetically Irreversible" & thermo_irr == T, "Thermodynamically Irreversible", reversibility)) %>%
              dplyr::select(reaction = rx, reversibility), by = "reaction")

# summarize each class in each markov sample
ML_rxn_summary <- ML_rxn_summary %>%
  group_by(rMech, Type, markovSample, logLik) %>%
  summarize(ML = sum(ML))

# looking at MAP estimate
  
MAP_ML_rxn_summary <- ML_rxn_summary %>% ungroup %>% group_by(rMech) %>% filter(logLik == max(logLik)) %>%
  dplyr::select(rMech, Type, ML)

ML_rxn_tall <- MAP_ML_rxn_summary %>% left_join(reaction_info_with_reversibility, by = "rMech") %>%
  ungroup %>%
  mutate(Type = factor(Type, levels = c("substrate", "product", "enzyme", "regulator")))

MAP_ML_rxn_summary <- tbl_df(dcast(MAP_ML_rxn_summary, rMech ~ Type, value.var = "ML", fill = 0)) # each rxn with associated enzyme, regulator and metabolic control fraction

MAP_ML_rxn_summary <- MAP_ML_rxn_summary %>% left_join(reaction_info_with_reversibility, by = "rMech") %>%
  arrange(substrate + product)

ML_rxn_tall <- ML_rxn_tall %>%
  ungroup %>%
  mutate(rMech = factor(rMech, levels = rev(MAP_ML_rxn_summary$rMech))) %>%
  mutate(reversibility = factor(reversibility, levels = c("Kinetically Reversible", "Kinetically Irreversible", "Thermodynamically Irreversible"))) %>%
  arrange(Type)

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

ggplot(ML_rxn_tall %>% arrange(Type), aes(x = rMech, y = ML, fill = Type)) + geom_bar(stat = "identity", width = 0.85, color = "BLACK") +
  barplot_facet_theme + scale_y_continuous(expression('Metabolic Leverage: ' ~ frac("|"~epsilon[i]~"|"~sigma[i], sum("|"~epsilon[j]~"|"~sigma[j], "j = 1" , n))), expand = c(0,0)) +
  scale_fill_manual(values = TypeColors) + facet_grid(~ reversibility, scales = "free_x", space = "free_x") +
  scale_x_discrete("Reactions", breaks = ML_rxn_tall$rMech, labels = ML_rxn_tall$abbrev)
ggsave(paste0("Figures/", rxn_subset, "-metabolicLeverageBar.pdf"), height = 10, width = 16)

# Pie chart summarizing metabolic leverage component contributions

pie_theme <- theme(text = element_text(size = 20, face = "plain"), title = element_text(size = 25), 
                       panel.background = element_rect(fill = "gray100"), legend.position = "top", 
                       axis.ticks = element_blank(), axis.text = element_blank(),
                       panel.grid = element_blank(), axis.line = element_blank(), legend.title=element_blank())

ML_component_summary <- ML_rxn_tall %>% group_by(reversibility, Type) %>% dplyr::summarize(ML = sum(ML)) %>% group_by(reversibility) %>% mutate(ML = ML/sum(ML))

ggplot(ML_component_summary %>%
         ungroup %>%
         mutate(reversibility = factor(reversibility, levels = c("Kinetically Reversible", "Kinetically Irreversible", "Thermodynamically Irreversible"))) %>%
         arrange(reversibility), 
                                       aes(x = factor(1), y = ML, fill = Type)) + geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) + coord_polar(theta = "y") +
  facet_wrap(~ reversibility) + scale_fill_manual(values = TypeColors) + pie_theme +
  scale_x_discrete("") + scale_y_continuous("")
ggsave(paste0("Figures/", rxn_subset, "-metabolicLeveragePie.pdf"), height = 5, width = 7)

# Generate a ternary plot to show the contributions of metabolites and enzymes in controlling flux

color_key <- color_ternary(100)
#color_key <- color_simplex(200, 0, T, 6)

color_index <- mapply(function(x,y){
  which.min(abs(color_key$Table$enzyme - x) + abs(color_key$Table$allostery - y))
}, x = MAP_ML_rxn_summary$enzyme, y = MAP_ML_rxn_summary$regulator)

MAP_ML_rxn_summary <- MAP_ML_rxn_summary %>% mutate(rMech = reaction) %>% dplyr::select(-reaction, -genes) %>% left_join(adequate_rxn_form_data, by = "rMech")


# Look at two alternative ways of describing the 3 components forming the ternary space
# substrates & products, enzymes, regulators
# substrates, enzymes, regulators & products

# substrates & products, enzymes, regulators

MAP_ML_rxn_summary_1 <- MAP_ML_rxn_summary %>% cbind(color_key$Table %>% dplyr::slice(color_index) %>% dplyr::select(color)) %>%
  mutate(rxn_metabolite = substrate + product)

ML_rxn_ternaryPoints_1 <- MAP_ML_rxn_summary_1 %>%  mutate(x = (1/2)*(2*regulator + enzyme) / (regulator + enzyme + rxn_metabolite),
                                              y = sqrt(3)/2 * enzyme*(regulator + enzyme + rxn_metabolite))

ML_rxn_ternaryPoints_1 <- ML_rxn_ternaryPoints_1 %>% left_join(ML_rxn_tall %>% select(reaction, Dir = reversibility) %>% unique(), by = c("rMech" = "reaction"))

anova(lm(MAP_ML_rxn_summary_1, formula = regulator ~ reversibility))

color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(Dir == "Kinetically Reversible"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(Dir == "Kinetically Reversible"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
ggsave(paste0("Figures/", rxn_subset, "-MLcolorKey_kinRev.pdf"), height = 9.1, width = 10.7)

color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(Dir == "Kinetically Irreversible"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(Dir == "Kinetically Irreversible"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
ggsave(paste0("Figures/", rxn_subset, "-MLcolorKey_kinIrrev.pdf"), height = 9.1, width = 10.7)

color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(Dir == "Thermodynamically Irreversible"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1 %>% dplyr::filter(Dir == "Thermodynamically Irreversible"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
ggsave(paste0("Figures/", rxn_subset, "-MLcolorKey_thermIrrev.pdf"), height = 9.1, width = 10.7)


ML_rxn_ternaryPoints_labels <- ML_rxn_ternaryPoints_1 %>% filter(regulator != 0 | enzyme > 0.33)

color_key$Figure_Color + 
  geom_point(data = ML_rxn_ternaryPoints_1, aes(x = x, y = y), size = 7, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints_1, aes(x = x, y = y), size = 6, shape = 21, fill = "WHITE") +
  ggrepel::geom_label_repel(data = ML_rxn_ternaryPoints_labels, aes(x = x, y = y, label = abbrev.x), color = "BLACK", segment.size = 0, label.padding = unit(0.15, "lines"))
ggsave(paste0("Figures/", rxn_subset, "-MLcolorKey.pdf"), height = 9.1, width = 10.7)

### Visualize uncertainty by plotting a subset of random markov samples 

n_random_points <- 100
set.seed(1234)
random_markov_samples <- sample(1:max(ML_rxn_summary$markovSample), n_random_points)

Posterior_ML_rxn_summary <- ML_rxn_summary %>% filter(markovSample %in% random_markov_samples)

Posterior_ML_rxn_summary <- tbl_df(dcast(Posterior_ML_rxn_summary, rMech + markovSample ~ Type, value.var = "ML", fill = 0)) # each rxn with associated enzyme, regulator and metabolic control fraction

Posterior_ML_rxn_summary <- Posterior_ML_rxn_summary %>% left_join(reaction_info_with_reversibility, by = "rMech")

color_index <- mapply(function(x,y){
  which.min(abs(color_key$Table$enzyme - x) + abs(color_key$Table$allostery - y))
}, x = Posterior_ML_rxn_summary$enzyme, y = Posterior_ML_rxn_summary$regulator)

Posterior_ML_rxn_summary <- Posterior_ML_rxn_summary %>% mutate(rMech = reaction) %>% dplyr::select(-reaction, -genes)

# Look at two alternative ways of describing the 3 components forming the ternary space
# substrates & products, enzymes, regulators
# substrates, enzymes, regulators & products

# substrates & products, enzymes, regulators

Posterior_ML_rxn_summary <- Posterior_ML_rxn_summary %>%
  cbind(color_key$Table %>% dplyr::slice(color_index) %>% dplyr::select(color)) %>%
  mutate(rxn_metabolite = substrate + product) %>%
  mutate(x = (1/2)*(2*regulator + enzyme) / (regulator + enzyme + rxn_metabolite),
         y = sqrt(3)/2 * enzyme*(regulator + enzyme + rxn_metabolite))

reversibility_order = c("Kinetically Reversible", "Kinetically Irreversible", "Thermodynamically Irreversible")

ML_uncertainty_plots <- list()
for(a_reversibility in reversibility_order){

color_key$Figure_BW + 
  geom_point(data = Posterior_ML_rxn_summary %>% dplyr::filter(reversibility == a_reversibility), aes(x = x, y = y, fill = color), size = 4, shape = 21)

ggsave(paste0("Figures/", rxn_subset, "-", a_reversibility, "-uncertainty.pdf"), height = 9.1, width = 10.7)
}








#control_layout <- ML_rxn_summary_1 %>% dplyr::select(reaction, reaction.name, abbrev, rxn_metabolite, enzyme, regulator, color)
#write.table(control_layout, "flux_cache/control_layout.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

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

all_affinities <- all_affinities %>% filter(rxn %in% adequate_rxn_form_data$reaction)

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
  scale_y_continuous(expression("Inferred" ~ italic("in vivo") ~ log[10] ~ 'Affinity (M) MLE')) +
  guides(fill = guide_legend(nrow = 3))
ggsave("Figures/brendaConsistency.pdf", width = 8, height = 10)

summary(lm(log10(absolute_comparison$parMLE) ~ absolute_comparison$log10mean))
  
brendaAgree <- absolute_comparison %>% mutate(sdOflog10 = ifelse(is.na(sdOflog10), 0, sdOflog10), 
                                              parCapture = log10(parLB) < log10mean & log10(parUB) > log10mean, 
                                              parCaptureInterval = log10(parLB) < log10mean + 2*sdOflog10 & log10(parUB) > log10mean - 2*sdOflog10) %>%
  dplyr::select(rxn, rxnForm, modelName, nQuant, parCapture, parCaptureInterval)

brendaAgree$parCapture %>% table() %>% as.data.frame() %>% summarize(Freq[. == T] / sum(Freq))
brendaAgree$parCaptureInterval %>% table() %>% as.data.frame() %>% summarize(Freq[. == T] / sum(Freq))

### Look at metabolism-wide occupancy ###
# need to fully load up to *import cluster parameters* in this script

occupancy_comparison <- affinity_comparisons %>% filter(measured) %>% dplyr::select(-formulaName, -measured, -absoluteQuant, -(EC:isYeast), -(nQuant:speciesType))

#

#constrained_conc <- occupancy_comparison %>%
#  ungroup %>%
#  select(rxn, commonName, Subtype, parLB = parLB, parUB = parUB, medianAbund) %>%
#  filter(Subtype == "substrate") %>%
#  mutate(occLB = 2^medianAbund / (2^medianAbund + parUB),
#         occUB = 2^medianAbund / (2^medianAbund + parLB),
#         parLB = log2(parLB),
 #        parUB = log2(parUB),
#         posterior_range = parUB - parLB,
#         occ_range = occUB - occLB) %>% arrange(desc(posterior_range))



# use concentrations of species across all conditions rather than just median

rxn_met_pairs <- occupancy_comparison %>% dplyr::select(rxnForm, modelName)

rxn_met_abundances <- lapply(1:nrow(rxn_met_pairs), function(i){
  metConc <- rxnList_form[[rxn_met_pairs[i,'rxnForm'] %>% unlist()]]$rxnMet[,rxn_met_pairs[i,'modelName'] %>% unlist()]
 cbind(rxn_met_pairs[i,'rxnForm'], rxn_met_pairs[i,'modelName'], logConc = metConc)
})
rxn_met_abundances <- do.call("rbind", rxn_met_abundances)

occupancy_comparison <- occupancy_comparison %>% left_join(rxn_met_abundances, by = c("modelName", "rxnForm"))

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

####

load("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/flux_cache/elasticityData.Rdata")

ELdata <- lapply(affinity_comparisons$rxnForm, function(x){ELdata[[x]] %>% mutate(rMech = x)}) %>%
  bind_rows %>% tbl_df

rMech_elasticities <- ELdata %>%
  group_by(rMech, specie, markovSample) %>%
  dplyr::summarize(medElasticity = median(Elasticity))

rMech_elasticities <- rMech_elasticities %>%
  group_by(rMech, specie) %>%
  dplyr::summarize(LB = quantile(medElasticity, 0.025),
                   UB = quantile(medElasticity, 0.975)) %>%
  inner_join(affinity_comparisons %>% filter(measured) %>% select(rMech = rxnForm, specie = modelName, commonName, Subtype), by = c("rMech", "specie")) %>%
  mutate(elasticity_range = UB - LB)

    
  
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
                                                            mutate(reversible = ifelse(reversible == "reversible", "Kinetically reversible", "Kinetically irreversible"),
                                                                   reversible = ifelse(cc_dGr + 1.96*cc_dGrSD < -5 & reversible == "Kinetically irreversible", "Thermodynamically irreversible", reversible)), by = "reaction")
free_energy_table <- free_energy_table %>% mutate(pathway = gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", pathway, perl=TRUE))

free_energy_table[,c("log_keq_LB")] <- log2(sapply(free_energy_table$rMech, function(x){rxn_fit_params[[x]]$param_interval$'X2.5.'[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("log_keq_UB")] <- log2(sapply(free_energy_table$rMech, function(x){rxn_fit_params[[x]]$param_interval$'X97.5.'[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("log_keq_MLE")] <- log2(sapply(free_energy_table$rMech, function(x){rxn_fit_params[[x]]$param_interval$MLE[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("logQmedian")] <- sapply(free_energy_table$rMech, function(x){as.numeric(rxn_fit_params[[x]]$param_interval$absoluteQuant[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"])})
free_energy_table <- free_energy_table %>% filter(!is.na(cc_dGr))

freebee_species <- c(water = "t_0399", proton = "t_0398", ammonium = "t_0233", diphosphate = "t_0332", co2 = "t_0249")
missing_species <- lapply(free_energy_table$rMech, function(x){
  unmeasuredSpec = rxn_fit_params[[x]]$kineticParPrior %>% filter(!is.na(measured) & !measured) %>% dplyr::select(SubstrateID = rel_spec, commonName)
  unmeasuredSpec = unmeasuredSpec %>% left_join(rxn_fit_params[[x]]$param_species %>% dplyr::select(SubstrateID, ReactionID, Subtype), by = "SubstrateID")
  unmeasuredSpec %>% filter(!(SubstrateID %in% freebee_species))
  })
missing_species <- do.call("rbind", missing_species)

# calculate Q for all conditions to compare relative to Keq

library(tidyr)

all_condition_Q <- lapply(1:nrow(free_energy_table), function(i){
  metStoi <- rxnList_form[[free_energy_table$rMech[i]]]$rxnStoi
  metConc <- rxnList_form[[free_energy_table$rMech[i]]]$rxnMet
  metConc[is.na(metConc)] <- 0
  metConc$condition <- rownames(metConc)
  
  metInfo <- metConc %>% gather(specie, conc, -condition) %>% mutate(specie = as.character(specie)) %>% left_join(data.frame(specie = names(metStoi), coef = unname(metStoi)), by = "specie")
  metInfo %>% tbl_df() %>% group_by(condition) %>% summarize(Q = sum(conc * coef)) %>% mutate(reaction = free_energy_table$reaction[i])
})
all_condition_Q <- do.call("rbind", all_condition_Q)

# keep track of Q relative to median Q
all_condition_Q <- free_energy_table %>% left_join(all_condition_Q, join = "reaction") %>% tbl_df() %>% dplyr::select(rMech, reaction, condition, Q, logQmedian)

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
free_energy_table <- free_energy_table %>% ungroup() %>% arrange(QkeqDiff_UB) %>% mutate(reaction = factor(reaction, levels = reaction))

# also flipping Q and logQmedian Q when keeping track of all conditions
all_conditions_rev <- all_condition_Q %>% filter(rMech %in% reverse_rxns$rMech) %>% mutate(Q = -1*Q, logQmedian = -1*logQmedian)
all_condition_Q <- all_condition_Q %>% filter(!(rMech %in% reverse_rxns$rMech)) %>% rbind(all_conditions_rev)

# comparing variable Q across conditions to the median condition where Q / keq was calculated
all_condition_Q <- all_condition_Q %>% mutate(reaction = factor(reaction, levels = free_energy_table$reaction)) %>% 
  left_join(free_energy_table %>% dplyr::select(rMech, reaction, reversible, genes, abbrev, pathway, measuredProducts, QkeqDiff_UB), join = "reaction") %>%
  mutate(Qkeq = (Q - logQmedian) + QkeqDiff_UB)

fast_C_lim <- all_condition_Q %>% filter(condition == "C0.30")

free_energy_theme <- boxplot_theme_label_rotate + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), strip.background = element_rect(fill = "coral"))

# Visualizing disequilibrium ratio

ggplot() + facet_grid(~ pathway, scale = "free_x", space = "free_x") +
  geom_violin(data = all_condition_Q, aes(x = reaction, y = 2^Qkeq, fill = reversible), scale = "width") +
  #geom_errorbar(data = free_energy_table, aes(x = reaction, ymin = 2^QkeqDiff_LB, ymax = 2^QkeqDiff_UB), size = 1) +
  geom_point(data = fast_C_lim, aes(x = reaction, y = 2^Qkeq), fill = "chartreuse3", size = 4, shape = 21) + 
  #geom_point(data = free_energy_table, aes(x = reaction, y = 2^QkeqDiff_MLE), fill = "darkblue", size = 4, shape = 21) +
  scale_y_log10(expression(frac(Q, K[eq]) ~ "=" ~ frac(v[r], v[f])), expand = c(0,0)) + coord_cartesian(ylim = c(10^-3, 1)) + free_energy_theme +
  scale_fill_manual(values = c("Kinetically reversible" = "cornflowerblue", "Kinetically irreversible" = "goldenrod1", "Thermodynamically irreversible" = "firebrick1")) +
  scale_x_discrete(breaks = free_energy_table$reaction, labels = free_energy_table$abbrev)
ggsave("Figures/reactionDisequilibrium.eps", width = 15, height = 8)
ggsave("Figures/reactionDisequilibrium.pdf", width = 15, height = 8)

##### Comparing unsupervised search for hypothetical metabolites to literature-guided search #####

# Compare ranking of hypothetical regulators to brute-force method of testing all metabolites


all_met_fits <- all_reactionInfo %>% tbl_df() %>%
  filter(modelType %in% c("regulator", "non-literature supported met regulator"),
         reaction %in% valid_rxns) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>%
  dplyr::select(reaction, rMech, Pvalue, ML) %>% mutate(reg_type = substr(rMech, 11, 20)) %>% tidyr::separate(reg_type, into = c("tID", "reg_type"), sep = "-") %>%
  dplyr::select(reaction, tID, reg_type, Pvalue, ML)

hypo_met_fits <- Hypo_met_candidates %>% tbl_df() %>% mutate(reg_type = substr(rMech, 18, 20), reaction = substr(rMech, 1, 6)) %>%
  left_join(all_reactionInfo, by = c("rMech", "reaction")) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>%
  dplyr::select(reaction, tID = SpeciesType, reg_type, ED = `2.5%`)

all_hypo_compare <- all_met_fits %>% left_join(hypo_met_fits, by = c("reaction", "tID", "reg_type"))

# compute ranks for maximum likelihood of individual metabolites and for hypo_met approach
all_hypo_rankings <- all_hypo_compare %>% group_by(reaction) %>%
  arrange(desc(ML)) %>% mutate(MLrank = 1:n()) %>%
  arrange(ED) %>% mutate(EDrank = 1:n())

library(venneuler)

rank_co <- 40
all_hypo_rankings <- all_hypo_rankings %>% mutate(ML_low = MLrank < rank_co, ED_low = EDrank < rank_co)

all_hypo_rank_spearman <- all_hypo_rankings %>% group_by(reaction) %>%
  dplyr::summarize(spearman = cor.test(MLrank, EDrank, method = "spearman")$estimate, Pvalue = cor.test(MLrank, EDrank, method = "spearman")$p.value)

cor.test(all_hypo_rankings$MLrank, all_hypo_rankings$EDrank, method = "spearman")

as.data.frame(all_hypo_rankings)[colnames(as.data.frame(all_hypo_rankings)) %in% c("ML_low", "ED_low")]

venneuler(as.data.frame(all_hypo_rankings)[colnames(as.data.frame(all_hypo_rankings)) %in% c("ML_low", "ED_low")])

all_hypo_rankings %>% dplyr::select(ML_low, ED_low)
  #require(MASS)
  #require(dplyr)
  #require(ggplot2)
  
  #hex_theme <- theme(text = element_text(size = 23), title = element_text(size = 25), panel.background = element_rect(fill = "gray40"), legend.position = "top", 
  #                   panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line"), axis.text = element_text(color = "black"),
  #                   axis.ticks = element_line(size = 1, color = "black")) 

  #hex_bin_max <- 125
  #n_hex_breaks <- 7
  #hex_breaks <- c(1,2,4,8,16,32,64,125)
  
  #ggplot() + 
  #  scale_x_continuous("woot", expand = c(0.02,0.02)) + 
  #  scale_y_continuous("woot", expand = c(0.02,0.02)) + 
  #  scale_fill_gradientn(name = "Counts", colours = c("lightyellow", "orange", "firebrick1", "indianred4"), trans = "log", breaks = hex_breaks, labels = hex_breaks) +
  #  hex_theme + geom_hex(data = all_hypo_rankings %>% mutate(z = 1), aes(x = MLrank, y = EDrank, z = z), bins = 50)
  




ggplot(all_hypo_rankings, aes(x = MLrank, y = EDrank)) + geom_hex()


### Fit of hypothetical regulation ###

load("flux_cache/paramCI.Rdata")

# Look at a subset of reactions

# remove reactions where the search for hypothetical regualtion is nearly underdetermined
relevant_rxns <- setdiff(valid_rxns, reactionInfo %>% filter(ncond - npar < 5, modelType == "hypo met regulator") %>%
                           dplyr::select(reaction) %>% unlist() %>% unname() %>% unique())
relevant_rxns <- rbind(
  # select understudied reactions
  data.frame(reaction = relevant_rxns[relevant_rxns %in% c("r_0195", "r_0225", "r_0468", "r_0514", "r_0528", "r_1049")], class = "understudied"),
  # reactions not fit by regulation or rMM
  data.frame(reaction = relevant_rxns[relevant_rxns %in% (optimal_reaction_form %>% filter(!(rMech %in% adequate_fit_optimal_rxn_form)))$reaction], class = "poorly-fit")
)
unregulated_rxns <- (optimal_reaction_form %>% filter(reaction %in% relevant_rxns$reaction) %>% filter(type == "unregulated"))$reaction
relevant_rxns <- relevant_rxns %>% group_by(reaction) %>% dplyr::slice(1)







## Compare hypothetical regulator to significant regulators to rMM ###

Hypo_met_compare <- reactionInfo %>% tbl_df() %>% left_join(rxn_fits %>% dplyr::select(rMech = rxn, spearman = parSpearman), by = "rMech") %>%
  filter(reaction %in% relevant_rxns$reaction, modelType %in% c("rMM", "regulator", "hypo met regulator")) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup() %>%
  filter(!(reaction %in% unregulated_rxns & modelType == "regulator"))

# Is the hypothetical regulator much more likely than any tested regulatory metabolite

hypo_met_AICc <- Hypo_met_compare %>% filter(modelType %in% c("regulator", "hypo met regulator")) %>%
  filter(Qvalue < 0.1) %>%
  dplyr::select(rMech, reaction, ncond, npar, ML, modelType) %>%
  mutate(AICc = 2*npar - 2*ML + 2*npar*(npar + 1)/(ncond - npar - 1)) %>%
  group_by(reaction, modelType) %>% dplyr::summarize(AICc = min(AICc)) %>%
  spread(modelType, AICc)

meaningful_hypo_met <- hypo_met_AICc %>% mutate(rel_lik = exp((`hypo met regulator` - regulator)/2)) %>%
  mutate(hypo_met_role = ifelse(!is.na(`hypo met regulator`) & (rel_lik < 0.5 | is.na(regulator)), "improves regulator", "n.s.")) %>%
  mutate(hypo_met_role = ifelse(hypo_met_role == "n.s." & !is.na(`hypo met regulator`), "improves rMM", hypo_met_role))
  

# sort reactions by the significance of hypothetical regulators
Hypo_met_ranking <- Hypo_met_compare %>% filter(modelType == "hypo met regulator") %>% left_join(
  Hypo_met_compare %>% filter(modelType == "rMM") %>% dplyr::select(reaction, rmm_ML = ML), by = "reaction"
)

# look at the hypothetical regulator (activation or inhibition) that is strongest
Hypo_met_ranking <- Hypo_met_ranking %>% group_by(reaction) %>% filter(Qvalue == min(Qvalue), (ML - rmm_ML) == max(ML - rmm_ML)) %>%
  ungroup() %>% mutate(Qvalue = round(Qvalue, 4), ML_diff = ML - rmm_ML) %>% dplyr::arrange(desc(Qvalue), ML_diff) %>%
  left_join(meaningful_hypo_met %>% dplyr::select(reaction, hypo_met_role), by = "reaction") %>% rowwise() %>%
  mutate(hypo_met_role = ifelse(is.na(hypo_met_role), "n.s.", hypo_met_role))


Hypo_y_rMM <- Hypo_met_compare %>% filter(modelType == "rMM" | rMech %in% Hypo_met_ranking$rMech) %>% rowwise() %>%
  mutate(printType = ifelse(modelType == "hypo met regulator",
                            Hypo_met_ranking$hypo_met_role[Hypo_met_ranking$rMech %in% rMech],
                            modelType)) %>%
  mutate(reaction = factor(reaction, levels = Hypo_met_ranking$reaction),
         printType = factor(printType, levels = c("rMM", "n.s.", "improves rMM", "improves regulator"))) %>% ungroup()

# pull out regulation by ascertained metabolites as well

measured_spear <- Hypo_met_compare %>% filter(modelType == "regulator", Qvalue < 0.1) %>% dplyr::select(rMech, reaction, spearman) %>%
  group_by(reaction) %>% filter(spearman == max(spearman)) %>% dplyr::slice(1)

stacked_reg_hypo_comp <- rbind(Hypo_y_rMM %>% dplyr::select(reaction, modelType, spearman, printType),
      measured_spear %>% mutate(modelType = "regulator", printType = "regulator") %>% dplyr::select(reaction, modelType, spearman, printType))

stacked_reg_hypo_comp <- stacked_reg_hypo_comp %>% group_by(reaction) %>% mutate(regRef = spearman[modelType == "rMM"]) %>%
  mutate(spearman = ifelse(modelType == "regulator" & spearman < regRef, regRef, spearman)) %>%
  mutate(hypoRef = max(spearman[modelType %in% c("rMM", "regulator")])) %>%
  mutate(spearman = ifelse(modelType == "hypo met regulator" & spearman < hypoRef, hypoRef, spearman))

stacked_reg_hypo_comp <- stacked_reg_hypo_comp %>% mutate(spearman = ifelse(modelType == "regulator", spearman - regRef, spearman),
         spearman = ifelse(modelType == "hypo met regulator", spearman - hypoRef, spearman))

stacked_reg_hypo_comp <- stacked_reg_hypo_comp %>% ungroup() %>% arrange(regRef) %>%
  left_join(relevant_rxns, by = "reaction") %>%
  mutate(reaction = factor(reaction, levels = unique(reaction)),
         modelType = factor(modelType, levels = c("rMM", "regulator", "hypo met regulator")) )
  


# This method is underpowered particularly for reactions when the 
# literature can help us discriminate between quantitatively similar metabolites

barplot_theme <- theme(text = element_text(size = 22), title = element_text(size = 25), 
                       panel.background = element_rect(fill = "gray96"), legend.position = "top", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black"),
                       axis.text.x = element_text(angle = 60, hjust = 1),
                       panel.grid = element_blank(),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
)

ggplot(stacked_reg_hypo_comp, aes(x = reaction, y = spearman, order = modelType)) + geom_hline(size = 1, yintercept = c(0, 0.25, 0.5, 0.75, 1)) +
  geom_bar(aes(fill = printType), color = "black", size = 0.25, width = 0.9, stat = "identity", position = "stack") +
  scale_y_continuous("Spearman Correlation", expand = c(0,0)) +
  scale_fill_manual(values = c("rMM" = "gray70", "regulator" = "gray50", "n.s." = "darkgoldenrod2","improves rMM" = "coral1", "improves regulator" = "red"),
                    breaks = c("rMM", "regulator", "n.s.", "improves rMM", "improves regulator")) +
  scale_x_discrete(name = "Reactions", breaks = fitReactionNames$reaction, labels = fitReactionNames$abbrev) +
  barplot_theme + facet_grid(~class, scale = "free_x", space = "free_x", drop = T)
ggsave("Figures/hypoMetImprove.pdf", height = 7.6, width = 9)


### Look at hypothetical 

Hypo_met_candidates <- Hypo_met_candidates %>% tbl_df() %>% mutate(modelType = ifelse(grepl('act', rMech), 'activator', 'inhibitor')) %>%
  filter(!grepl('cc$', rMech)) %>% left_join(reactionInfo %>% dplyr::select(rMech, reaction, ncond, npar, Qvalue), by = "rMech") %>%
  filter(reaction %in% Hypo_y_rMM$reaction[Hypo_y_rMM$printType == "improves regulator"]) %>%
  group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup() %>%
  mutate(reaction = factor(reaction, levels = levels(Hypo_y_rMM$reaction)[levels(Hypo_y_rMM$reaction) %in% reaction]))

n_label = 5
x_cutoff = 0.8

all_hypo_met_candidates <- Hypo_met_candidates %>% dplyr::select(reaction, modelType, SpeciesName, `97.5%`) %>%
  group_by(reaction, modelType, `97.5%`) %>% dplyr::slice(n()) %>% # filter species w/ identical correlation (same measurement)
  group_by(reaction, modelType) %>% arrange(desc(`97.5%`)) %>% # arrange by correlation
  mutate(rank = 1:n(),
         label = ifelse(rank <= n_label, SpeciesName, '')) %>%
  mutate(ypos = which(levels(reaction) %in% reaction[1])+runif(n(), -0.2, 0.2)) %>%
  filter(`97.5%` > x_cutoff)

y_axis_label = data.frame(reaction = levels(all_hypo_met_candidates$reaction),
                          ypos = 1:length(levels(all_hypo_met_candidates$reaction))
                          ) %>% left_join(fitReactionNames, by = "reaction")

barplot_theme_facet <- barplot_theme + theme(panel.margin = unit(2, "lines"),
                                       panel.background = element_rect(fill = "gray96"),
                                       strip.background = element_rect(fill = "darkslategray2"))

ggplot(all_hypo_met_candidates, aes(x = `97.5%`, y = ypos, label = label)) + facet_grid(~ modelType) +
  geom_point(aes(color = ifelse(label != "", "RED", "BLACK"))) +
  geom_text(size = 3, angle = 20, hjust = -0.1, vjust = 0, color = "black") +
  coord_cartesian(xlim = c(0.8, 1), ylim = c(0.5, max(y_axis_label$ypos)+1.5)) +
  scale_y_discrete(name = "Reactions", breaks = y_axis_label$ypos, labels = y_axis_label$abbrev) +
  scale_x_continuous("Correlation of measured metabolites to optimal regulator") +
  barplot_theme_facet +
  scale_color_identity()
ggsave("Figures/hypoMetCandidates.pdf", height = 6, width = 10)


