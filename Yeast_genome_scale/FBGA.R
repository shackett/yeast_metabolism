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
  source("reactionStructures.r")
}else{
  load("flux_cache/rxnf_formulametab.rdata")
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
rmCondList <- data.frame() # reactions which will be duplicated with some conditions rmoved


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


### Pull out relevent reaction information - such as reaction ID, number of parameters and maximum likelihood ###

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

reactionInfo$changeP <- NA
reactionInfo$changeAIC <- NA

for(rx in c(1:nrow(reactionInfo))[reactionInfo$form != "rm" | !(reactionInfo$modification %in% c("", "rmCond"))]){ 
  # Both regulation and irreversible michaelis-menten kinetics are compared to reversible michaelis-menten kinetics
  rxn_eval <- reactionInfo[rx,]
  if(length(grep('rmCond', rxn_eval$modification)) == 0){
    rxn_ref <- reactionInfo[reactionInfo$reaction == rxn_eval$reaction & reactionInfo$form == "rm" & reactionInfo$modification == "",]
  }else{
    rxn_ref <- reactionInfo[reactionInfo$reaction == rxn_eval$reaction & reactionInfo$form == "rm" & reactionInfo$modification == "rmCond",]
  }
  
  likDiff <- rxn_eval$ML - rxn_ref$ML
  
  # calculate corrected AIC
  rxn_eval$AICc <- 2*rxn_eval$npar - 2*rxn_eval$ML + (2*rxn_eval$npar*(rxn_eval$npar + 1))/(rxn_eval$ncond - rxn_eval$npar - 1)
  rxn_ref$AICc <- 2*rxn_ref$npar - 2*rxn_ref$ML + (2*rxn_ref$npar*(rxn_ref$npar + 1))/(rxn_ref$ncond - rxn_ref$npar - 1)
  
  # determine the relative probability of a null or full model based on AIC
  if(rxn_eval$AICc < rxn_ref$AICc){
    ref_p <- exp((rxn_eval$AICc - rxn_ref$AICc)/2)
    reactionInfo$changeAIC[rx] <- 1 - 1/(ref_p + 1)
    }else{
      alt_p <- exp((rxn_ref$AICc - rxn_eval$AICc)/2)
      reactionInfo$changeAIC[rx] <- 1 - alt_p/(alt_p + 1)
      }
  
  if(rxn_eval$npar == rxn_ref$npar){
    # Alternative parameterization with same degrees of freedom
    reactionInfo$changeP[rx] <- 1/(exp(likDiff) + 1)
  }else if(rxn_eval$npar > rxn_ref$npar){
    # Alternative model is more complex
    reactionInfo$changeP[rx] <- 1 - pchisq(2*likDiff, rxn_eval$npar - rxn_ref$npar)
  }else{
    # Alternative model is less complex
    reactionInfo$changeP[rx] <- 1 - pchisq(-2*likDiff, rxn_ref$npar - rxn_eval$npar)
     cat(paste0("\nrx ", rx, ": ",-2*likDiff))
  }
}

### Compare each reaction form including a non-1 hill coefficient to a reaction form where this is fixed ###

reactionInfo$hillP <- NA

for(rx in grep('ultra', reactionInfo$modification)){
  rxn_eval <- reactionInfo[rx,]
  rxn_ref <- reactionInfo[reactionInfo$rMech %in% sub('_ultra', '', rxn_eval$rMech),] # find equivalent without possible ultra-sensitivity
  if(nrow(rxn_ref) == 0){
    print(paste(rxn_eval$rMech, "is not paired with a non-ultrasensitive form")); next
    }
  likDiff <- rxn_eval$ML - rxn_ref$ML
  
  reactionInfo$hillP[rx] <- 1 - pchisq(2*likDiff, rxn_eval$npar - rxn_ref$npar)
  }



### Identify reaction form modification which significantly improve the likelihood function ###
### comparisons are relative to reversible michaelis-menten kinetics ###

library(qvalue)
qStore <- list()

reactionInfo$Qvalue <- NA

# comparisons between CC and RM are pathological because they are so similar, resulting in high p-values - remove from consideration
qStore[["CC_RM"]] <- qvalue(reactionInfo$changeP[reactionInfo$form == "cc"], pi0.method = "bootstrap")
reactionInfo <- reactionInfo[!(reactionInfo$form == "cc"),]

# Comparison between forms of interest

# comparison between irreversible MM-kinetics and reversible MM-kinetics
# The p-value of an irreversible reaction reflect how much worse it fits relative to reversible MM-kinetics (i.e. a p-value of 1 if rMM = iMM)
test_subset <- grepl('^(forward|reverse)', reactionInfo$modification)
qStore[["MM_reversibility"]] <- qvalue(reactionInfo$changeP[test_subset])
reactionInfo$Qvalue[test_subset] <- qStore[["MM_reversibility"]]$q

### hypothetical metabolites with constrained hill coefficient vs. reversible-MM (AIC)
test_subset <- !is.na(reactionInfo$changeP) & c(1:nrow(reactionInfo)) %in% intersect(grep('metX', reactionInfo$modification), grep('ultra', reactionInfo$modification, invert = T))
qStore[["hREG_RM"]] <- qvalue(reactionInfo$changeAIC[test_subset])
reactionInfo$Qvalue[test_subset] <- qStore[["hREG_RM"]]$q

### hypothetical metabolites with an unconstrained hill coefficient (spike & slab prior) vs. both R-MM (AIC) and constrained hill (LRT)
test_subset <- !is.na(reactionInfo$changeP) & c(1:nrow(reactionInfo)) %in% intersect(grep('metX', reactionInfo$modification), grep('ultra', reactionInfo$modification, invert = F))
qStore[["hREGhill_RM"]] <- qvalue(reactionInfo$changeAIC[test_subset])
qStore[["hREGhill_hREG"]] <- qvalue(reactionInfo$hillP[test_subset])
# Q-values to note are those that have a significant effect relative to both MM and to a constrained hill coefficient
reactionInfo$Qvalue[test_subset] <- mapply(function(x,y){max(x,y)}, x = qStore[["hREGhill_RM"]]$q, y = qStore[["hREGhill_hREG"]]$q)

### literature metabolites with constrained hill coefficient vs. reversible-MM
test_subset <- !is.na(reactionInfo$changeP) & c(1:nrow(reactionInfo)) %in% intersect(grep('t_[0-9]{4}', reactionInfo$modification), grep('ultra', reactionInfo$modification, invert = T))
qStore[["lREG_RM"]] <- qvalue(reactionInfo$changeP[test_subset])
reactionInfo$Qvalue[test_subset] <- qStore[["lREG_RM"]]$q

### literature metabolites with an unconstrained hill coefficient (spike & slab prior) vs. both R-MM (AIC) and constrained hill (LRT)
test_subset <- !is.na(reactionInfo$changeP) & c(1:nrow(reactionInfo)) %in% intersect(grep('t_[0-9]{4}', reactionInfo$modification), grep('ultra', reactionInfo$modification, invert = F))
qStore[["lREGhill_RM"]] <- qvalue(reactionInfo$changeP[test_subset])
qStore[["lREGhill_lREG"]] <- qvalue(reactionInfo$hillP[test_subset])
reactionInfo$Qvalue[test_subset] <- mapply(function(x,y){max(x,y)}, x = qStore[["lREGhill_RM"]]$q, y = qStore[["lREGhill_lREG"]]$q)

### Test proposed complex regulation (combine with literature tests to generate a more stable distribution of p-values
manualRegulators <- read.delim('./companionFiles/manual_ComplexRegulation.txt')
test_subset <- reactionInfo$rMech %in% manualRegulators$TechnicalName

manualComplex_q <- qvalue(c(reactionInfo$changeP[test_subset], qStore[["lREGhill_RM"]]$pvalues))
manualComplex_q$pvalues <- manualComplex_q$pvalues[1:sum(test_subset)]; manualComplex_q$qvalues <- manualComplex_q$qvalues[1:sum(test_subset)]
qStore[["mREG"]] <- manualComplex_q
reactionInfo$Qvalue[test_subset] <- qStore[["mREG"]]$q

### Assign a symbolic indicator of significance based on q-value

signifCO <- data.frame(q_cutoff = c(0.1, 0.001, 0.00001), code = c("*", "**", "***"))

reactionInfo$signifCode <- sapply(reactionInfo$Qvalue, function(x){
  if(is.na(x) | x > signifCO$q_cutoff[1]){
    ""
  }else{
    rev(signifCO$code[x < signifCO$q_cutoff])[1]
  }
})



#### Setup indexing so that fits can be interactively analysed using a Shiny application #####
  
reactionInfo$Name = NA
reactionInfo$Name[reactionInfo$modification %in% c('', 'rmCond') & reactionInfo$form == "rm"] <- "Reversible michaelis-menten (default)"
reactionInfo$Name[reactionInfo$modification %in% c('', 'rmCond') & reactionInfo$form == "cc"] <- "Convenience kinetics"
reactionInfo$Name[grepl('^(forward|reverse)', reactionInfo$modification)] <- "Irreversible michaelis-menten"

for(rxN in grep('act|inh', reactionInfo$modification)){
  
  regName <- substr(reactionInfo$modification[rxN], 1, 6)
  rxRegs <- rxnList_form[[reactionInfo$rMech[rxN]]]$rxnFormData[rxnList_form[[reactionInfo$rMech[rxN]]]$rxnFormData$Type != "rct",]
  
  reactionInfo$Name[rxN] <- paste(sub('(^)([a-z])', '\\U\\2', rxRegs$Subtype[rxRegs$SubstrateID == regName], perl = T), 
                                  ifelse(rxRegs$Type[rxRegs$SubstrateID == regName] == "act", "activation", "inhibition"), 
                                  ifelse(length(grep('ultra', reactionInfo$modification[rxN])) == 0, "", "(variable hill)"),
                                  "by", unname(rxnList_form[[reactionInfo$rMech[rxN]]]$metNames[names(rxnList_form[[reactionInfo$rMech[rxN]]]$metNames) == regName]))
}

reactionInfo$Name[grep('rmCond', reactionInfo$modification)] <- sapply(reactionInfo$Name[grep('rmCond', reactionInfo$modification)], function(x){paste(x, "(zero flux reactions removed)")})
reactionInfo$Name[reactionInfo$rMech %in% manualRegulators$TechnicalName] <- manualRegulators$DisplayName[chmatch(reactionInfo$rMech[reactionInfo$rMech %in% manualRegulators$TechnicalName], manualRegulators$TechnicalName)] # rescue manually named reactions

reactionInfo$Name <- mapply(function(x,y){paste(x,y)}, x = reactionInfo$signifCode, y = reactionInfo$Name)

reactionInfo$Name <- gsub('^ ', '', reactionInfo$Name)
reactionInfo$Name <- gsub('  ', ' ', reactionInfo$Name)

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

reactionInfo$FullName <- mapply(function(x,y){paste(x, y, sep = " - ")}, x = rxToPW$reactionName[chmatch(reactionInfo$reaction, rxToPW$rID)], y = reactionInfo$Name)
  

pathwaySet <- sort(table(rxToPW$pathway), decreasing = T)
pathwaySet <- data.frame(pathway = names(pathwaySet), members = unname(pathwaySet), display = paste(names(pathwaySet), ' (', unname(pathwaySet), ')', sep = ""))


#### Generate reaction plots and summaries ####

load("companionFiles/PTcomparison_list.Rdata") # by-gene comparisons of protein and transcript abundance

shiny_flux_data <- list()
rxn_fits <- NULL
rxn_fit_params <- list()
fraction_flux_deviation <- NULL
MLdata <- NULL # Save summary of metabolic leverage
ELdata <- list() # Save full distribution of elasticities
TRdata <- NULL # Save summary transcriptional resposiveness

t_start = proc.time()[3]

# flag a few reactions where additional plots are created
custom_plotted_rxns <- c("r_1054-im-forward", "r_1054-rm","r_0962-rm", "r_0962-rm-t_0292-act-mm", "r_0962-rm-t_0290-act-mm", "r_0816-rm", "r_0816-rm-t_0461-inh-uncomp", "r_0816-rm_rmCond", "r_0816-rm-t_0461-inh-uncomp_rmCond",
                         "r_0514-im-forward", "r_0514-rm", "r_0916-im-forward", "r_0916-rm", "r_0208-im-forward", "r_0208-rm")
custom_plotted_rxns <- c("r_0816-rm", "r_0816-rm-t_0461-inh-uncomp", "r_0816-rm_rmCond", "r_0816-rm-t_0461-inh-uncomp_rmCond")

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
    shiny_flux_data[[arxn]]$plotChoices <- append(shiny_flux_data[[arxn]]$plotChoices, hypoMetTrend(run_rxn, metSVD, tab_boer))
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
# significant or default reaction forms

pathway_plot_list <- list()
for(pw in pathwaySet$display){
  # iterate through pathways and plot pathway-level figures
  pathway_plot_list[[pw]] <- pathwayPlots(pw)
  }

#### Save lists which will be processed by Shiny app ####

save(pathwaySet, rxToPW, reactionInfo, pathway_plot_list, file = "shinyapp/shinyData.Rdata")

# generate a minute version of shinyData that will load quickly when the App is being modified
#reactionInfo <- reactionInfo[1:20,]
#shiny_flux_data <- shiny_flux_data[names(shiny_flux_data) %in% reactionInfo$rMech]
#system.time(save(pathwaySet, rxToPW, reactionInfo, pathway_plot_list, shiny_flux_data, file = "shinyapp/shinySubData.Rdata"))


#### Save parameter estimates for further global analyses ####

save(rxn_fit_params, rxn_fits, reactionInfo, MLdata, TRdata, fraction_flux_deviation, file = "flux_cache/paramCI.Rdata")
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
  if(any(substrates == "nm")){
   x$rxnID
  }
})
validRxnA <- unlist(validRxnA) %>% unname()

###

validRxnB <- rxn_fits %>% tbl_df() %>% dplyr::select(rxn, parSpearman) %>% filter(parSpearman > 0.6) %>% left_join(reactionInfo %>% tbl_df() %>% dplyr::select(rxn = rMech, reaction, modification, Qvalue)) %>%
  filter(!grepl('t_metX', rxn)) %>% filter(is.na(Qvalue) | Qvalue < 0.05) %>% dplyr::select(reaction) %>% unlist() %>% unname() %>% unique()

valid_rxns <- union(validRxnA, validRxnB)
rmCond_rxns <- unique(reactionInfo$reaction[grep('rmCond', reactionInfo$modification)]) # reactions which carry zero flux under some conditions - consider only non-zero reactions

optimal_rxn_form <- sapply(valid_rxns, function(x){
  
  rx_forms <- reactionInfo[reactionInfo$reaction == x,]
  if(x %in% rmCond_rxns){
    rx_forms <- rx_forms[grep('rmCond', rx_forms$modification),]
    }
  rx_forms <- rx_forms %>% dplyr::filter(form != "im")
  
  rx_forms <- rx_forms[grep('t_metX', rx_forms$modification, invert = T),]
  rx_best_mod <- rx_forms$rMech[!is.na(rx_forms$Qvalue)][which.min(rx_forms$Qvalue[!is.na(rx_forms$Qvalue)])]
  if(length(rx_best_mod) == 0){
    rx_forms$rMech[rx_forms$modification == ""]
  }else if(rx_forms$Qvalue[rx_forms$rMech == rx_best_mod] < 0.05){
    rx_best_mod
  }else{
    rx_forms$rMech[rx_forms$modification %in% c("rmCond", "")]
  }
})


#reactionInfo[chmatch(optimal_rxn_form, reactionInfo$rMech),]

##### Summary based on interval overlap #####

interval_overlap_summary <- data.table(rxnForm = fraction_flux_deviation$rxn, intervalOverlap = fraction_flux_deviation$'Interval Overlap')
interval_overlap_summary$reaction = substr(interval_overlap_summary$rxnForm, 1, 6)
interval_overlap_summary <- interval_overlap_summary[interval_overlap_summary$reaction %in% valid_rxns,]

interval_overlap_summary$Type = NA
interval_overlap_summary$Type[interval_overlap_summary$rxnForm %in% reactionInfo$rMech[reactionInfo$modification %in% c("", "rmCond")]] <- "Substrates and Enzymes"
interval_overlap_summary$Type[grep('metX', interval_overlap_summary$rxnForm)] <- "+ hypothetical activator or inhibitor"
interval_overlap_summary <- interval_overlap_summary %>% dplyr::filter(!grepl('(forward|reverse)', rxnForm))
interval_overlap_summary$Type[is.na(interval_overlap_summary$Type)] <- "+ literature activator or inhibitor"

# Reduce to best overlapping basic reaction or regulations
interval_overlap_summary <- interval_overlap_summary[,list(predictionOverlap = max(intervalOverlap)), by = c("reaction", "Type")]
# add reactions with no literature supported regulation back with their RM values (and note with color)
interval_overlap_match <- interval_overlap_summary[!(interval_overlap_summary$reaction %in% interval_overlap_summary$reaction[interval_overlap_summary$Type %in% "+ literature activator or inhibitor"]) & interval_overlap_summary$Type == "Substrates and Enzymes",]
interval_overlap_match$Type <- "+ literature activator or inhibitor"

interval_overlap_summary <- rbind(data.frame(interval_overlap_summary, color = "firebrick1"), data.frame(interval_overlap_match, color = "darkgoldenrod1"))

# Order each type seperately
interval_overlap_summary$x <- NA
for(a_type in unique(interval_overlap_summary$Type)){
  interval_overlap_summary$x[interval_overlap_summary$Type == a_type][order(interval_overlap_summary$predictionOverlap[interval_overlap_summary$Type == a_type])] <- 1:sum(interval_overlap_summary$Type == a_type)
  }

interval_overlap_summary$Type <- factor(interval_overlap_summary$Type, levels = c("Substrates and Enzymes", "+ literature activator or inhibitor", "+ hypothetical activator or inhibitor"))
interval_overlap_summary$x <- factor(interval_overlap_summary$x)

interval_exp_summary <- data.table(interval_overlap_summary)
interval_exp_summary <- interval_exp_summary[,list(label = paste(signif(mean(predictionOverlap) * 100, 3), '% of fluxes\nare consistent', sep = '')) , by = "Type"]


barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(),
                         strip.background = element_rect(fill = "cornflowerblue"), panel.margin = unit(1.5, "lines"), axis.text.y = element_text(size = 20, color = "black"))


ggplot() + facet_grid(Type~.) + geom_bar(data = interval_overlap_summary, aes(x = x, y = predictionOverlap, fill = color), stat = "identity", position = "dodge", width = 0.85) + barplot_theme +
 scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_y_continuous(name = "Flux prediction ", expand = c(0,0)) + scale_fill_identity("Prediction Method") +
 geom_text(data = interval_exp_summary, aes(label = label), x = 5, y = 0.75, color = "black", hjust = 0) +
 ggtitle('Fraction of confidence intervals \n capturing experimental value')

ggsave("Figures/intervalOverlapSummary.pdf", height = 14, width = 10)



### Just looking at MM and signficant regulation ###

interval_overlap_summary <- data.table(rxnForm = fraction_flux_deviation$rxn, intervalOverlap = fraction_flux_deviation$'Interval Overlap')
interval_overlap_summary <- interval_overlap_summary[interval_overlap_summary$rxnForm %in% optimal_rxn_form,]
interval_overlap_summary$Type <- ifelse(reactionInfo$modification[chmatch(interval_overlap_summary$rxnForm, reactionInfo$rMech)] == "", "rMM", "regulator")
setkey(interval_overlap_summary, "intervalOverlap")
interval_overlap_summary$rxnForm <- factor(interval_overlap_summary$rxnForm, levels = interval_overlap_summary$rxnForm)

ggplot() + geom_bar(data = interval_overlap_summary , aes(x = rxnForm, y = intervalOverlap, fill = Type), stat = "identity", position = "dodge", width = 0.85) + barplot_theme + theme(axis.title.y = element_blank()) +
 scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_y_continuous(name = "", expand = c(0,0)) + scale_fill_brewer("", palette = "Set1", label = c("Allosteric Regulator", "Reversible Michaelis-Menten")) +
 ggtitle('Fraction of confidence intervals \n capturing experimental value')
ggsave("Figures/intervalOverlapOptimal.pdf", height = 7, width = 10)



##### Summary based on spearman correlation #####

spearman_MM <- rxn_fits$parSpearman[reactionInfo$modification == "" & substr(rxn_fits$rxn, 1, 6) %in% valid_rxns]
spearman_MM <- sort(spearman_MM[spearman_MM > 0])
spearman_MM_df <- data.frame(x = 1:length(spearman_MM), corr = spearman_MM)

ggplot() + geom_bar(data = spearman_MM_df, aes(x = x, y = spearman_MM), stat = "identity", position = "dodge", width = 0.85, fill = "coral") + barplot_theme +
 scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_y_continuous(name = "Spearman correlation: prediction ~ FBA", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1)) + scale_fill_identity() + expand_limits(y = c(0,1))
ggsave("Figures/MMspearmanCorr.pdf", height = 12, width = 13)

#### Summary based on spearman correlation for MM and most significant regulator (if applicable) #####

spearman_MMandReg <- data.frame(reactionInfo[,c('reaction', 'modification', 'Qvalue')], spearman = rxn_fits[,'parSpearman']) %>% tbl_df()  # all reactions
spearman_MMandReg <- spearman_MMandReg %>% filter(!grepl('t_metX', modification)) %>% filter(is.na(Qvalue) | Qvalue < 0.1) # filter hypothetical regulators and take MM or Qvalue < 0.1
spearman_MMandReg <- spearman_MMandReg %>% dplyr::filter(!(grepl('^(forward|reverse)', modification))) # filter irreversible MM-kinetics
spearman_MMandReg <- spearman_MMandReg %>% filter(reaction %in% valid_rxns)
spearman_MMandReg <- spearman_MMandReg %>% left_join(spearman_MMandReg %>% group_by(reaction) %>% dplyr::summarize(rmCond_reaction = ifelse(sum(modification == 'rmCond') , T, F)))

rmCond_sets <- spearman_MMandReg %>% ungroup() %>% filter(grepl('rmCond', modification)) %>% dplyr::select(reaction, modification) # Take the rmCond reactions when they exist and the normal reactions otherwise

spearman_MMandReg <- rbind(spearman_MMandReg %>% filter(rmCond_reaction == F),
spearman_MMandReg %>% filter(rmCond_reaction == T) %>% inner_join(rmCond_sets))

best_reg_form <- spearman_MMandReg %>% group_by(reaction) %>% dplyr::summarize(modification = modification[which.min(Qvalue)][1]) %>%
  filter(!is.na(modification))

select_spearman_MMandReg <- rbind(spearman_MMandReg %>% filter(modification %in% c("", "rmCond")),
                                  spearman_MMandReg %>% inner_join(best_reg_form))

select_spearman_MMandReg <- select_spearman_MMandReg %>% dplyr::mutate(regulator = ifelse(!(modification %in% c("", "rmCond")), T, F)) %>%
  dplyr::select(reaction, spearman, regulator)

select_spearman_MMandReg <- select_spearman_MMandReg %>% mutate(rxName = paste0(reaction, ifelse(regulator, "-reg", "-rm"))) %>%
  group_by(reaction) %>% dplyr::mutate(maxSpear = max(spearman)) %>% ungroup() %>% dplyr::arrange(maxSpear) %>%
  mutate(rxName = factor(rxName, levels = rxName))

all_neg_reactions <- select_spearman_MMandReg %>% group_by(reaction) %>% dplyr::summarize(npos = length(spearman[spearman > 0])) %>% filter(npos == 0)

select_spearman_MMandReg <- select_spearman_MMandReg %>% filter(!(reaction %in% all_neg_reactions$reaction))

barplot_theme_nox <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_blank(),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )

ggplot(select_spearman_MMandReg, aes(x = rxName, y = spearman, fill = regulator)) + geom_bar(stat = "identity", color = "BLACK", width = 0.85) +
  barplot_theme_nox + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_fill_manual("", values = c("skyblue2", "orangered1"), breaks = c(T, F), label = c("Metabolite Regulator", "Reversible Michaelis-Menten")) +
  ggtitle('Correlation between measured and predicted flux') + expand_limits(y = c(0,1))

ggsave("Figures/MM_with_reg_spearman.pdf", height = 8, width = 12)


# two plots

mm_spearman <- select_spearman_MMandReg %>% dplyr::filter(!regulator)

reg_spearman <- select_spearman_MMandReg %>% dplyr::filter(regulator) %>% dplyr::filter(spearman == maxSpear) %>% group_by(reaction) %>%
  dplyr::mutate(spearDiff = spearman - mm_spearman$spearman[mm_spearman$reaction == reaction])
reg_spearman$spearDiff[reg_spearman$spearDiff > reg_spearman$spearman] <- reg_spearman$spearman[reg_spearman$spearDiff > reg_spearman$spearman]

two_pane_spearman_summary <- rbind(
mm_spearman %>% dplyr::mutate(pane = "MM"),
rbind(mm_spearman %>% dplyr::filter(!(reaction %in% reg_spearman$reaction)), reg_spearman %>% dplyr::select(-spearDiff)) %>%
  dplyr::mutate(pane = "REG")
)

mm_order <- mm_spearman %>% dplyr::arrange(spearman) %>% dplyr::select(reaction) %>% unlist() %>% unname()
two_pane_spearman_summary <- two_pane_spearman_summary %>% dplyr::mutate(reaction = factor(reaction, levels = mm_order))

# increase facets spacing - color facets
ggplot(two_pane_spearman_summary, aes(x = reaction, y = spearman, fill = regulator)) + geom_bar(stat = "identity", color = "black", width = 0.8) +
  facet_grid(pane ~ .) + ggtitle('Correlation between measured and predicted flux') +
  barplot_theme_nox + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_fill_manual("", values = c("skyblue2", "orangered1"), breaks = c(T, F), label = c("Metabolite Regulator", "Reversible Michaelis-Menten"))

ggsave("Figures/spearman_two_pane.pdf", height = 9, width = 12)

# stack regulator

stacked_spearman <- rbind(mm_spearman %>% dplyr::select(reaction, spearman, regulator),
                          reg_spearman %>% dplyr::select(reaction, spearman = spearDiff, regulator))
stacked_spearman <- stacked_spearman %>% dplyr::mutate(reaction = factor(reaction, levels = mm_order))

ggplot(stacked_spearman %>% filter(spearman >= 0), aes(x = reaction, y = spearman, fill = regulator)) + geom_bar(stat = "identity", color = "black", width = 0.8) +
  ggtitle('Correlation between measured and predicted flux') +
  barplot_theme_nox + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_fill_manual("", values = c("skyblue2", "orangered1"), breaks = c(T, F), label = c("Metabolite Regulator", "Reversible Michaelis-Menten"))

ggsave("Figures/spearman_stack.pdf", height = 7, width = 12)

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
  barplot_theme_nox + scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_fill_manual("", values = c("Reversible" = "skyblue2", "Irreversible" = "moccasin"), label = c("Irreversible Michaelis-Menten", "+ Reversibility"))

ggsave("Figures/spearman_stack_rev.pdf", height = 7, width = 12)

##### Replotting a few reactions for figures #####

allostery_affinity()

##### Generate figure summarizing metabolic leverage for a condition #####

MLdata$Type[grepl("r_0302", MLdata$reaction) & MLdata$specie == "isocitrate"] <- "product" # aconitase is split into two reactions without a meaasured cis-aconitate so isocitrate acts like a product in the first rxn
MLdata <- data.table(MLdata %>% filter(conditions == "NATURAL"))
metabolic_leverage_summary_plots("P0.05")

# Only looking at reactions that are well-fit
# Look at the best significant reaction form for reaction where CI overlap with flux carried is substantial (>50%)

adequate_fit_optimal_rxn_form <- union(intersect(optimal_rxn_form, fraction_flux_deviation$rxn[fraction_flux_deviation$"Interval Overlap" > 0.5]), intersect(optimal_rxn_form, rxn_fits$rxn[rxn_fits$parSpearman > 0.6]))

# check for cases where variable hill coefficient regulation may have been missing
# check arginosuccinate lyase |- Arginine (r_0207)
# glutamate dehydrogenase |- quinolinate (r_0471)

significant_hill <- reactionInfo %>% filter(grepl('ultra', modification) & Qvalue < 0.1)

# when multiple regulatory candidates have similar fits but one is a far better candidate based upon the literature, choose the literature-supported one
reactionInfo %>% filter(reaction == "r_0962")
adequate_fit_optimal_rxn_form[grep('r_0962-rm', adequate_fit_optimal_rxn_form)] <- "r_0962-rm-t_0290-act-mm" # pyruvate kinase activation by F16bisP over citrate inhibition
reactionInfo %>% filter(reaction == "r_0215") 
adequate_fit_optimal_rxn_form[grep('r_0215-rm', adequate_fit_optimal_rxn_form)] <- "r_0215-rm-t_0499-inh-uncomp_ultra" # aspartate kinase regulation by ultrasensitive inhbition by threonine

# summarize metabolic leverage

ML_inducibility_summary <- enzyme_control_source()



######## Split ML of well-fit reactions into enzyme and allosteric control ###########

ML_rxn_summary <- tbl_df(MLdata) %>% dplyr::filter(reaction %in% adequate_fit_optimal_rxn_form) %>% 
  dplyr::mutate(Type = ifelse(Type %in% c("substrate", "product"), "rxn_metabolite", Type)) %>%
  dplyr::mutate(Type = ifelse(!(Type %in% c("rxn_metabolite", "enzyme")), "regulator", Type)) %>%
  dplyr::select(reaction, specie, Type, condition, ML = get("0.5")) %>% group_by(reaction, Type, condition) %>%
  dplyr::summarize(ML = sum(ML)) %>% group_by(reaction, Type) %>% dplyr::summarize(ML = median(ML)) %>%
  group_by(reaction) %>% dplyr::mutate(ML = ML / sum(ML))

ML_rxn_summary <- tbl_df(dcast(ML_rxn_summary, reaction ~ Type, value.var = "ML", fill = 0)) # each rxn with assocaited enzyme, regulator and metabolic control fraction

ML_rxn_summary <- ML_rxn_summary %>% dplyr::mutate(rID = substr(reaction, 1,6)) %>% left_join(ML_inducibility_summary %>% dplyr::select(reaction, genes, reversibility)) %>%
  left_join(reversibleRx %>% dplyr::select(rID = rx, CCdG, CCdGsd)) %>% arrange(rxn_metabolite)


color_key <- color_ternary(100)
#color_key <- color_simplex(200, 0, T, 6)

color_index <- mapply(function(x,y){
  which.min(abs(color_key$Table$enzyme - x) + abs(color_key$Table$allostery - y))
}, x = ML_rxn_summary$enzyme, y = ML_rxn_summary$regulator)

ML_rxn_summary <- ML_rxn_summary %>% cbind(color_key$Table %>% dplyr::slice(color_index) %>% dplyr::select(color))

ML_rxn_ternaryPoints <- ML_rxn_summary %>%  mutate(x = (1/2)*(2*regulator + enzyme) / (regulator + enzyme + rxn_metabolite),
                                              y = sqrt(3)/2 * enzyme*(regulator + enzyme + rxn_metabolite))

summary(lm(ML_rxn_summary, formula = rxn_metabolite ~ reversibility))

# Reversible
color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints %>% dplyr::filter(reversibility == "T"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints %>% dplyr::filter(reversibility == "T"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
ggsave("Figures/MLcolorKey_REV.pdf", height = 9.1, width = 10.7)

# Irreversible
color_key$Figure_BW + 
  geom_point(data = ML_rxn_ternaryPoints %>% dplyr::filter(reversibility == "F"), aes(x = x, y = y), size = 9, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints %>% dplyr::filter(reversibility == "F"), aes(x = x, y = y, fill = color), size = 8, shape = 21)
  ggsave("Figures/MLcolorKey_FOR.pdf", height = 9.1, width = 10.7)


# Name point according to the metabolic layout

rxn_names <- c('r_1838' = 'HCS', 'r_0988' = 'SDH', 'r_0915' = 'PPAT', 'r_0042' = 'DAHP synthase', 'r_0886' = 'PFK', 'r_0718' = 'ME',
  'r_0468' = 'G5K', 'r_0962' = 'PyK', 'r_4040' = 'PPK', 'r_0225' = 'ATP-PRTase', 'r_0887' = 'S7P PFK', 'r_0888' = 'PGM',
  'r_0195' = 'TPS', 'r_0214' = 'ATCase', 'r_0310' = 'CBL', 'r_0816' = 'OTCase', 'r_0450' = 'ALD', 'r_0215' = 'AspK',
  'r_0491' = 'G3PDH', 'r_0250' = 'CPS')

ML_rxn_ternaryPoints_labels <- ML_rxn_ternaryPoints %>% filter(regulator != 0) %>% rowwise() %>% mutate(label = rxn_names[names(rxn_names) == rID])

color_key$Figure <- color_key$Figure + 
  geom_point(data = ML_rxn_ternaryPoints, aes(x = x, y = y), size = 7, shape = 21, fill = "BLACK") +
  geom_point(data = ML_rxn_ternaryPoints, aes(x = x, y = y), size = 6, shape = 21, fill = "WHITE") +
  geom_text(data = ML_rxn_ternaryPoints_labels, aes(x = x + 0.05, y = y + 0.05, label = label), color = "WHITE")

#color_key$Figure + geom_point(data = ML_rxn_summary, aes(x = enzyme, y = regulator), color = "yellow")

#plot(1:45, 1:45, col = ML_rxn_summary$color, pch = 16, cex = 2)
ggsave("Figures/MLcolorKey.pdf", color_key$Figure, height = 9.1, width = 10.7)

ML_rxn_summary$reaction_info <- reactionInfo$FullName[chmatch(ML_rxn_summary$reaction, reactionInfo$rMech)]

write.table(ML_rxn_summary, "flux_cache/control_layout.tsv", col.names = T, row.names = F, quote = F, sep = "\t")




######### Compare Km values found through optimization to those from literature ################

all_affinities <- read.delim("flux_cache/metaboliteAffinities.tsv") # all BRENA substrates and regulators
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
    rxnData <- rxnData %>% filter(is.na(isYeast) | !isYeast)
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

boxplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black", size = 20),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
                       )


ggplot(absolute_comparison, aes(x = log10mean, y = log10(parMLE), ymin = log10(parLB), ymax = log10(parUB))) + geom_errorbar(size = 0.6, alpha = 0.5, width = 0.1) + geom_smooth(method = "lm", se = F, size = 2) +
  geom_abline(a = 0, b = 1, size = 2) + boxplot_theme +
  scale_x_continuous(expression("Literature" ~ log[10] ~ "Affinity (M)")) +
  scale_y_continuous(expression("Inferred" ~ italic("in vivo") ~ log[10] ~ "Affinity (M) 95% credicibility interval")) + coord_cartesian(ylim = c(-8,0))
ggsave("Figures/brendaAffinity.pdf", width = 8, height = 10)
ggsave("Figures/brendaAffinity.eps", width = 8, height = 10, device=cairo_ps)

brendaAgree <- absolute_comparison %>% mutate(sdOflog10 = ifelse(is.na(sdOflog10), 0, sdOflog10), parCapture = log10(parLB) < log10mean & log10(parUB) > log10mean, parCaptureInterval = log10(parLB) < log10mean + 2*sdOflog10 & log10(parUB) > log10mean - 2&sdOflog10) %>%
  dplyr::select(rxn, rxnForm, parCapture, parCaptureInterval)

brendaAgree$parCapture %>% table()
brendaAgree$parCaptureInterval %>% table()

### Look at metabolism-wide occupancy ###

occupancy_comparison <- affinity_comparisons %>% filter(measured) %>% dplyr::select(-formulaName, -measured, -absoluteQuant, -(EC:speciesType))

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

species_of_interest <- c("ATP", "ADP")
#species_of_interest <- c("ATP", "ADP", "AMP", "NADH", "NAD")
occupancy_comparison$Specie[occupancy_comparison$commonName %in% species_of_interest] <- occupancy_comparison$commonName[occupancy_comparison$commonName %in% species_of_interest]
occupancy_comparison$Specie[is.na(occupancy_comparison$Specie)] <- "The Rest"

occupancy_comparison <- occupancy_comparison %>% mutate(log10S_km = log10(2^logConc / parMLE), log10occ = (2^logConc/(2^logConc + parMLE)))

occupancy_comparison <- occupancy_comparison %>% mutate(Subtype = factor(Subtype, levels = c("substrate", "product", "regulator")))
occupancy_comparison$Specie <- factor(occupancy_comparison$Specie, levels = c(species_of_interest, "The Rest"))

spec_counts <- occupancy_comparison %>% group_by(Specie, Subtype, rxn) %>% summarize(counts = n()) %>% group_by(Specie, Subtype) %>% summarize(counts = n())

boxplot_theme_label_rotate <- boxplot_theme + theme(axis.title.y = element_text(angle = 0))

# S / Km
ggplot() + geom_violin(data = occupancy_comparison, aes(x = Specie, y = log10S_km), fill = "firebrick2", scale = "width") +
  facet_grid(Subtype ~ .) + geom_text(data = spec_counts, aes(x = Specie, y = 0, label = counts), color = "BLUE", size = 6) +
  scale_y_continuous(expression(frac(S, K[M])), breaks = seq(-5, 5, by = 1), labels = 10^seq(-5, 5, by = 1), expand = c(0,0)) +
  boxplot_theme_label_rotate
ggsave("Figures/speciesSKm.pdf", width = 10, height = 10)
ggsave("Figures/speciesSKm.eps", width = 10, height = 10, device=cairo_ps)

# S / S + Km
ggplot() + geom_violin(data = occupancy_comparison, aes(x = Specie, y = log10occ), fill = "firebrick2", scale = "width") +
  facet_grid(Subtype ~ .) + geom_text(data = spec_counts, aes(x = Specie, y = 0.5, label = counts), color = "BLUE", size = 6) +
  scale_y_continuous(expression(frac(S, S + K[M])), breaks = seq(0, 1, by = 0.2), label = sapply(seq(0, 1, by = 0.2) * 100, function(x){paste(x, "%", sep = "")})) +
  boxplot_theme_label_rotate
ggsave("Figures/speciesOccupancy.pdf", width = 10, height = 10)
ggsave("Figures/speciesOccupancy.eps", width = 10, height = 10, device=cairo_ps)


######## Compare Keq to Gibbs free energy ########

# rxForm defined keq: prod[S] - prod[P]/keq
# calculate median absolute concentration or bounds based upon reasonable concentrations of small molecules
# relate optimized keq to these values
# Q = prod[P]/prod[S]
# Q/keq

reaction_names <- read.delim('companionFiles/control_layout.tsv')
reaction_free_energy <- read.delim("companionFiles/cc_dG_python.tsv")

free_energy_table <- data.frame(RXform = adequate_fit_optimal_rxn_form, reaction = substr(adequate_fit_optimal_rxn_form, 1, 6))

free_energy_table <- free_energy_table %>% left_join( reversibleRx %>% dplyr::select(reaction = rx, cc_dGr = CCdG, cc_dGrSD = CCdGsd, reversible = modelBound) %>% mutate(reversible = ifelse(reversible == "reversible", "Reversible", "Irreversible")))

free_energy_table[,c("log_keq_LB")] <- log2(sapply(free_energy_table$RXform, function(x){rxn_fit_params[[x]]$param_interval$'X2.5.'[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("log_keq_UB")] <- log2(sapply(free_energy_table$RXform, function(x){rxn_fit_params[[x]]$param_interval$'X97.5.'[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("log_keq_MLE")] <- log2(sapply(free_energy_table$RXform, function(x){rxn_fit_params[[x]]$param_interval$MLE[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"]}))
free_energy_table[,c("logQmedian")] <- sapply(free_energy_table$RXform, function(x){as.numeric(rxn_fit_params[[x]]$param_interval$absoluteQuant[rxn_fit_params[[x]]$kineticParPrior$rel_spec == "keq"])})
free_energy_table <- free_energy_table %>% filter(!is.na(cc_dGr))

free_energy_table <- free_energy_table %>% left_join(reaction_names %>% dplyr::select(RXform = reaction, genes, abbrev, pathway))
free_energy_table$genes[is.na(free_energy_table$genes)] <- free_energy_table$abbrev[is.na(free_energy_table$abbrev)] <- free_energy_table$reaction[is.na(free_energy_table$abbrev)]
free_energy_table$pathway[is.na(free_energy_table$pathway)] <- "other"

freebee_species <- c(water = "t_0399", proton = "t_0398", ammonium = "t_0233", diphosphate = "t_0332", co2 = "t_0249")
missing_species <- lapply(free_energy_table$RXform, function(x){
  unmeasuredSpec = rxn_fit_params[[x]]$kineticParPrior %>% filter(!is.na(measured) & !measured) %>% dplyr::select(SubstrateID = rel_spec, commonName)
  unmeasuredSpec = unmeasuredSpec %>% left_join(rxn_fit_params[[x]]$param_species %>% dplyr::select(SubstrateID, ReactionID, Subtype), by = "SubstrateID")
  unmeasuredSpec %>% filter(!(SubstrateID %in% freebee_species))
  })
missing_species <- do.call("rbind", missing_species)

# calculate Q for all conditions to compare relative to Keq

library(tidyr)

all_condition_Q <- lapply(1:nrow(free_energy_table), function(i){
  metStoi <- rxnList[[free_energy_table$RXform[i]]]$rxnStoi
  metConc <- rxnList[[free_energy_table$RXform[i]]]$rxnMet
  metConc[is.na(metConc)] <- 0
  metConc$condition <- rownames(metConc)
  
  metInfo <- metConc %>% gather(specie, conc, -condition) %>% mutate(specie = as.character(specie)) %>% left_join(data.frame(specie = names(metStoi), coef = unname(metStoi)), by = "specie")
  metInfo %>% tbl_df() %>% group_by(condition) %>% summarize(Q = sum(conc * coef)) %>% mutate(reaction = free_energy_table$reaction[i])
})
all_condition_Q <- do.call("rbind", all_condition_Q)

# keep track of Q relative to median Q
all_condition_Q <- free_energy_table %>% left_join(all_condition_Q, join = "reaction") %>% tbl_df() %>% dplyr::select(RXform, reaction, condition, Q, logQmedian)

# back to one estimate-per-reaction
free_energy_table <- free_energy_table %>% mutate(QkeqDiff_MLE = logQmedian - log_keq_MLE, QkeqDiff_UB = logQmedian - log_keq_LB, QkeqDiff_LB = logQmedian - log_keq_UB) %>% dplyr::select(-c(log_keq_LB, log_keq_UB, log_keq_MLE, logQmedian))

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
all_conditions_rev <- all_condition_Q %>% filter(RXform %in% reverse_rxns$RXform) %>% mutate(Q = -1*Q, logQmedian = -1*logQmedian)
all_condition_Q <- all_condition_Q %>% filter(!(RXform %in% reverse_rxns$RXform)) %>% rbind(all_conditions_rev)

# comparing variable Q across conditions to the median condition where Q / keq was calculated
all_condition_Q <- all_condition_Q %>% mutate(reaction = factor(reaction, levels = free_energy_table$reaction)) %>% left_join(free_energy_table %>% dplyr::select(RXform, reaction, reversible, genes, abbrev, pathway, measuredProducts, QkeqDiff_UB), join = "reaction") %>%
  mutate(Qkeq = Q - logQmedian + QkeqDiff_UB)

fast_C_lim <- all_condition_Q %>% filter(condition == "C0.30")

free_energy_theme <- boxplot_theme_label_rotate + theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background = element_rect(fill = "coral"))

# Visualizing disequilibrium ratio

ggplot() + facet_grid(~ pathway, scale = "free_x", space = "free_x") +
  geom_violin(data = all_condition_Q, aes(x = reaction, y = 2^Qkeq, fill = reversible), scale = "width") +
  geom_errorbar(data = free_energy_table, aes(x = reaction, ymin = 2^QkeqDiff_LB, ymax = 2^QkeqDiff_UB), size = 1) +
  geom_point(data = fast_C_lim, aes(x = reaction, y = 2^Qkeq), fill = "chartreuse3", size = 3, shape = 21) + 
  geom_point(data = free_energy_table, aes(x = reaction, y = 2^QkeqDiff_MLE), fill = "darkblue", size = 3, shape = 21) +
  scale_y_log10(expression(frac(Q, K[eq]) ~ "=" ~ frac(v[r], v[f])), expand = c(0,0)) + coord_cartesian(ylim = c(10^-3, 1)) + free_energy_theme +
  scale_fill_discrete() +
  scale_x_discrete(breaks = free_energy_table$reaction, labels = free_energy_table$abbrev)

free_energy_table %>% mutate(Gstandard_LB = cc_dGr - 2*cc_dGrSD, Gstandard_UB = cc_dGr + 2*cc_dGrSD) %>% ungroup() # delta G / RT
(8.315 * (273 + 25) / 1000

##### Determine which metabolite or protein is most correlated with pathway flux for each pathway #####
# topology: metabolism_stoichiometry - set of all directed reactions

library(igraph)

load("flux_cache/yeast_stoi_directed.Rdata") # S: metabolism stoichiometry
load("flux_cache/reconstructionWithCustom.Rdata") # metabolic reconstruction files

carried_flux <- read.table("Flux_analysis/fluxCarriedSimple.tsv", header = T, sep = "\t")
S_carried <- S; colnames(S_carried) <- sub('_[FR]$', '', colnames(S_carried))
S_carried <- S_carried[,colnames(S_carried) %in% rownames(carried_flux)] # reactions which carry flux
S_carried <- S_carried[rowSums(S_carried != 0) != 0,] # metabolites which are utilized
S_carried <- S_carried[!(rownames(S_carried) %in% corrFile$SpeciesID[grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', corrFile$SpeciesName)]),] # remove cofactors and common species
#grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', corrFile$SpeciesName, value = T)


# convert stoichiometric matrix to a directed bipartite graph
S_graph <- melt(t(S_carried), varnames = c("source", "sink"))
S_graph <- S_graph[S_graph$value != 0,] # reactions - metabolite links
S_graph$source <- as.character(S_graph$source); S_graph$sink <- as.character(S_graph$sink) # class coercion

S_graph[S_graph$value < 0,] <- S_graph[S_graph$value < 0,][,c(2,1,3)] # for consumed metabolites, invert direction
S_graph <- S_graph[,-3]

S_igraph <- graph.data.frame(S_graph)
V(S_igraph)$type <- V(S_igraph) %in% colnames(S_carried)

sort(betweenness(S_igraph), decreasing = T)[1:4] # remove water, H+, ATP, ADP

# define a pathway operationally as a set of reactions satisfying:
# connected through other pathway members
# proportional flux through all members (pearson correlation = 1)

# define clusters of highly correlated reactions

fluxCorrelation <- abs(cor(t(carried_flux[rownames(carried_flux) %in% colnames(S_carried),])))
reaction_adjacency <- graph.adjacency(fluxCorrelation > 0.999)
reaction_clusters <- igraph::clusters(reaction_adjacency, mode = "weak") # find clusters of correlaated reacrtions

# identify connected subgraphs using a subset of reactions defined by these correlated reaction clusters (without restricting metabolities)

pathway_set <- list()
for(clust_n in 1:reaction_clusters$no){
  
  if(length(V(reaction_adjacency)$name[reaction_clusters$membership == clust_n]) < 3){next} # catch minute clusters
  
  # reduce stoichiometric matrix to the relevent set of reactions
  pw_subset <- S_carried[,colnames(S_carried) %in% V(reaction_adjacency)$name[reaction_clusters$membership == clust_n]]
  pw_subset <- pw_subset[,!duplicated(colnames(pw_subset))]
  pw_subset <- pw_subset[rowSums(pw_subset != 0) != 0,]
  
  if(nrow(pw_subset) < 4){next}
  
  # determine the weakly connected subgraphs within this set of reactions
  S_graph_subset <- melt(t(pw_subset), varnames = c("source", "sink"))
  S_graph_subset <- S_graph_subset[S_graph_subset$value != 0,] # reactions - metabolite links
  S_graph_subset$source <- as.character(S_graph_subset$source); S_graph_subset$sink <- as.character(S_graph_subset$sink) # class coercion

  S_graph_subset[S_graph_subset$value < 0,] <- S_graph_subset[S_graph_subset$value < 0,][,c(2,1,3)] # for consumed metabolites, invert direction
  S_graph_subset <- S_graph_subset[,-3]
  
  S_igraph_subset <- graph.data.frame(S_graph_subset)
  V(S_igraph_subset)$type <- V(S_igraph_subset) %in% colnames(pw_subset)

  reaction_clusters_subset <- igraph::clusters(S_igraph_subset, mode = "weak")
  # save clusters which involve more than 3 reactions
  cluster_rxns <- lapply(1:reaction_clusters_subset$no, function(x){
    clust_members <- V(S_igraph_subset)$name[reaction_clusters_subset$membership == x]
    clust_members[clust_members %in% colnames(pw_subset)]
    })
  
  for(a_clust in 1:reaction_clusters_subset$no){
    rxns <- cluster_rxns[[a_clust]]
    
    if(length(rxns) >= 3){
      pw_output <- pw_subset[,colnames(pw_subset) %in% rxns]
      pw_output <- pw_output[rowSums(pw_output != 0) != 0,]
      
      pathway_set <- append(pathway_set, list(stoi = pw_output))
      
      }
    }
  }

# load metabolomics and proteomics

metabolite_abundance <- read.delim('flux_cache/tab_boer_log2rel.txt')

enzyme_abund <- read.delim("./companionFiles/proteinAbundance.tsv")
rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
rxn_enzyme_groups <- read.delim("./flux_cache/rxn_enzyme_groups.tsv")
  
nbs <- 10000

 
rxn_pathways_carryFlux <- rxn_pathways[rxn_pathways$reactionID %in% rownames(carried_flux),]
rxn_pathways_class_list <- sapply(rxn_pathways_carryFlux$pathway, function(x){strsplit(x, split = '__')[[1]]})
names(rxn_pathways_class_list) <- rxn_pathways_carryFlux$reactionID
pw_counts <- table(unlist(rxn_pathways_class_list))

for(a_pathway in 1:length(pathway_set)){
  stoi_subset <- pathway_set[[a_pathway]]
  
  pw_flux <- apply(carried_flux[rownames(carried_flux) %in% colnames(stoi_subset),], 2, median)
  pw_flux <- pw_flux * ifelse(median(pw_flux) < 0, -1, 1) # set positive direction as direction carried in the majority of conditions
  
  # determine the name of this pathway using hypergeometric
  
  pw_annotations <- unname(unlist(rxn_pathways_class_list[colnames(stoi_subset)]))
  if(is.null(pw_annotations)){
    pw <- "Unspecified"
  }else{
    pw_annot <- table(pw_annotations)
    pw_p <- rep(NA, length(pw_annot))
    for(i in 1:length(pw_annot)){
      pw_p[i] <- 1 - phyper(q = pw_annot[i], m = pw_counts[names(pw_counts) == names(pw_annot[i])], n = sum(pw_counts) - pw_annot[i], k = sum(pw_annot))
    }
    pw <- names(pw_annot[pw_p == min(pw_p)][which.max(pw_annot[pw_p == min(pw_p)])])
  }
  
  rx_names <- unique(rxnFile$Reaction[rxnFile$ReactionID %in% colnames(stoi_subset)])  
  
# metabolites
  matching_met_index <- sapply(unique(corrFile$SpeciesType[corrFile$SpeciesID %in% rownames(stoi_subset)]), function(a_met){grep(a_met, metabolite_abundance$tID)})
  matching_met_df <- metabolite_abundance[unique(unlist(matching_met_index)),, drop = F]
  
  if(nrow(matching_met_df) > 0){
    matching_met_matrix <- 2^matching_met_df[,grep('[PCNLU][0-9.]{4}', colnames(matching_met_df))]
    if(!all(names(pw_flux) %in% colnames(matching_met_matrix))){stop("check ordering")}
    matching_met_matrix <- matching_met_matrix[,chmatch(colnames(matching_met_matrix), names(pw_flux)), drop = F]
    
    pathway_corr <- t(sapply(1:nrow(matching_met_matrix), function(a_row){
      # regress metabolite/enzyme on flux 
      # determine a null distribution using a pivotal t-statistic-based bootstrap
      
      pw_reg <- lm(pw_flux ~ matching_met_matrix[a_row,])
      reg_resids <- pw_reg$resid * sqrt(length(pw_flux)/pw_reg$df.residual)
      
      pivotal_t <- t(sapply(1:nbs, function(z){ 
        pw_null <- I(pw_reg$fitted + sample(reg_resids, replace = T))
        summary(lm(pw_null ~ matching_met_matrix[a_row,]))$coef[2,1:2]
        }))
      pivotal_dist <- pivotal_t[,1]/pivotal_t[,2]
      
      data.frame(pwnum = a_pathway, class = "metabolite", specie = matching_met_df$Metabolite[a_row], corr = cor(pw_flux, matching_met_matrix[a_row,]), pval = (1  - abs(0.5 - sum(pivotal_dist > 0)/nbs)*2) + (1/nbs))
      
      }))
      
    #pathway_corr <- t(sapply(1:nrow(matching_met_matrix), function(a_row){
    #  null = sapply(1:nperm, function(z){cor(pw_flux, sample(matching_met_matrix[a_row,]))})
    #  data.frame(pwnum = a_pathway, class = "metabolite", specie = matching_met_df$Metabolite[a_row], corr = cor(pw_flux, matching_met_matrix[a_row,]), pval = 1 - sum(cor(pw_flux, matching_met_matrix[a_row,]) > null)/(nperm + 1))
    #  }))
    
    }else{
      pathway_corr <- NULL
    }
  
  # enzymes
  enzyme_subset <- 2^enzyme_abund[rownames(enzyme_abund) %in% unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% colnames(stoi_subset)]),,drop = F]
  colnames(enzyme_subset) <- toupper(colnames(enzyme_subset))
  
  if(nrow(enzyme_subset) > 0){
    if(!all(names(pw_flux) %in% colnames(enzyme_subset))){stop("check ordering")}
    matching_enzyme_matrix <- enzyme_subset[,chmatch(colnames(enzyme_subset), names(pw_flux)), drop = F]
    
    enzyme_corr <- t(sapply(1:nrow(matching_enzyme_matrix), function(a_row){
      # regress metabolite/enzyme on flux 
      # determine a null distribution using a pivotal t-statistic-based bootstrap
      
      pw_reg <- lm(pw_flux ~ matching_enzyme_matrix[a_row,])
      reg_resids <- pw_reg$resid * sqrt(length(pw_flux)/pw_reg$df.residual)
      
      pivotal_t <- t(sapply(1:nbs, function(z){ 
        pw_null <- I(pw_reg$fitted + sample(reg_resids, replace = T))
        summary(lm(pw_null ~ matching_enzyme_matrix[a_row,]))$coef[2,1:2]
        }))
      pivotal_dist <- pivotal_t[,1]/pivotal_t[,2]
      
      data.frame(pwnum = a_pathway, class = "enzyme", specie = rownames(matching_enzyme_matrix)[a_row], corr = cor(pw_flux, matching_enzyme_matrix[a_row,]), pval = (1  - abs(0.5 - sum(pivotal_dist > 0)/nbs)*2) + (1/nbs))
      
      }))
    
    #enzyme_corr <- t(sapply(1:nrow(matching_enzyme_matrix), function(a_row){
    #  null = sapply(1:nperm, function(z){cor(pw_flux, sample(matching_enzyme_matrix[a_row,]))})
    #  data.frame(pwnum = a_pathway, class = "enzyme", specie = rownames(matching_enzyme_matrix)[a_row], corr = cor(pw_flux, matching_enzyme_matrix[a_row,]), pval = 1 - sum(cor(pw_flux, matching_enzyme_matrix[a_row,]) > null)/(nperm + 1))
    #  }))
    pathway_corr <- rbind(pathway_corr, enzyme_corr)
    }
  pathway_set[[a_pathway]] <- list(pw_name = pw, rx_names = rx_names, stoi = pathway_set[[a_pathway]], corr = pathway_corr, flux = pw_flux)
}

library(qvalue)

pw_associations <- sapply(1:length(pathway_set), function(x){pathway_set[[x]][[2]]})
pw_associations <- do.call(rbind, pw_associations)
pw_associations <- as.data.frame(apply(pw_associations, c(1,2), as.character))
pw_associations$qval <- qvalue(as.numeric(pw_associations$pval))$qvalues

# some metabolites track flux through pathways but many enzymes and metabolites are anti-correlated with flux-carried

sig_pw_associations <- pw_associations[pw_associations$qval < 0.05,]

sig_met_pw_associations <- sig_pw_associations[sig_pw_associations$class == "metabolite",]

# plot significant examples

pw_flux <- apply(carried_flux[rownames(carried_flux) %in% colnames(stoi_subset),], 2, median)
pw_flux <- pw_flux * ifelse(median(pw_flux) < 0, -1, 1) # set positive direction as direction carried in the majority of conditions










library(colorRamps)
heatmap.2(abs(cor(t(carried_flux))), trace = "none", col = blue2yellow(100))
heatmap.2((fluxCorrelation > 0.999)*1, trace = "none", col = green2red(2))
table(abs(cor(t(carried_flux))) > 0.999)



# igraph tcl-tk interface ?
tkplot(S_igraph)
tkplot(S_igraph)
layout.norm # rescale = F
# load
# add species
# save
# subset for relevant species