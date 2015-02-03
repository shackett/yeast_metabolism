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
    # cat(paste0("\nrx ", rx, ": ",-2*likDiff))
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
custom_plotted_rxns <- c("r_1054-im-forward", "r_1054-rm","r_0962-rm", "r_0962-rm-t_0292-act-mm", "r_0962-rm-t_0290-act-mm", "r_0816-rm_rmCond", "r_0816-rm-t_0461-inh-uncomp_rmCond")

#rxn_subset <- grep('1054|0962|0816', reactionInfo$rMech, value = T)
#for(arxn in rxn_subset){
for(arxn in custom_plotted_rxns){
#for(arxn in reactionInfo$rMech){
  
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
    cat(paste(round((which(reactionInfo$rMech == arxn) / length(reactionInfo$rMech))*100, 2), "% complete - ", round((proc.time()[3] - t_start)/60, 0), " minutes elapsed", sep = ""))
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

### Determine which reactions to follow-up on based upon either
## A) Contain all substrates
## B) RM is well-fit despite missing substrates

RMMrxns <- reactionInfo$rMech[reactionInfo$modification %in% c("rmCond", "")]

load("paramOptim.Rdata")
measure_exception <- c("H+", "H2O")
validRxnA <- sapply(rxnList_form[names(rxnList_form) %in% RMMrxns], function(x){
  if(all(x$flux$standardQP >= 0)){
    substrates <- x$originMet[chmatch(names(x$rxnStoi)[x$rxnStoi < 0], names(x$originMet))]
  }else if(all(x$flux$standardQP <= 0)){
    substrates <- x$originMet[chmatch(names(x$rxnStoi)[x$rxnStoi > 0], names(x$originMet))]
  }else{
    substrates <- x$originMet
  }
  
  substrates <- substrates[names(substrates) %in% names(x$metNames)[!(x$metNames %in% measure_exception)]]
  ifelse(any(substrates == "nm"), F, T)
})

validRxnB <- rxn_fits$parSpearman[chmatch(RMMrxns, rxn_fits$rxn)] > 0.6

valid_rxns <- substr(RMMrxns[validRxnA | validRxnB], 1, 6) %>% unique()
rmCond_rxns <- unique(reactionInfo$reaction[grep('rmCond', reactionInfo$modification)]) # reactions which carry zero flux under some conditions - consider only non-zero reactions

optimal_rxn_form <- sapply(valid_rxns, function(x){
  
  rx_forms <- reactionInfo[reactionInfo$reaction == x,]
  if(x %in% rmCond_rxns){
    rx_forms <- rx_forms[grep('rmCond', rx_forms$modification),]
    }
  
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

##### Summary based on spearman correlation for MM or significant regulator #####

spearman_MMandReg <- data.table(data.frame(reactionInfo[,c('reaction', 'modification', 'Qvalue')], spearman = rxn_fits[,'parSpearman'])) # all reactions
spearman_MMandReg <- spearman_MMandReg[grep('t_metX', spearman_MMandReg$modification, invert = T),] # filter hypothetical regulators

spearman_MMandReg <- spearman_MMandReg[spearman_MMandReg$Qvalue < 0.1 | is.na(spearman_MMandReg$Qvalue),] # filter for MM or regulation with Qvalue < 0.1
spearman_MMandReg$row <- 1:nrow(spearman_MMandReg)

spearman_MMandReg <- spearman_MMandReg[spearman_MMandReg[, row[which.max(spearman)], by = "reaction"]$V1,]
spearman_MMandReg <- spearman_MMandReg[spearman > 0,]
setkey(spearman_MMandReg, "spearman")
spearman_MMandReg$Type <- ifelse(spearman_MMandReg$modification %in% c("", "rmCond"), "rMM", "regulator")

spearman_MMandReg$reaction <- factor(spearman_MMandReg$reaction, levels = spearman_MMandReg$reaction)

ggplot() + geom_bar(data = spearman_MMandReg , aes(x = reaction, y = spearman, fill = Type), stat = "identity", position = "dodge", width = 0.85) + barplot_theme +
  scale_y_continuous(name = "Spearman correlation", expand = c(0,0), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_fill_brewer("", palette = "Set1", label = c("Allosteric Regulator", "Reversible Michaelis-Menten")) +
  ggtitle('Correlation between measured and predicted flux') + expand_limits(y = c(0,1))
ggsave("Figures/MM_reg_spearmanCorr.pdf", height = 12, width = 13)

##### Summary based on spearman correlation for MM and most significant regulator (if applicable) #####

spearman_MMandReg <- data.frame(reactionInfo[,c('reaction', 'modification', 'Qvalue')], spearman = rxn_fits[,'parSpearman']) %>% tbl_df()
spearman_MMandReg <- spearman_MMandReg %>% filter(!grepl('t_metX', modification)) %>% filter(is.na(Qvalue) | Qvalue < 0.1)
spearman_MMandReg <- spearman_MMandReg %>% filter(reaction %in% valid_rxns)
spearman_MMandReg <- spearman_MMandReg %>% left_join(spearman_MMandReg %>% group_by(reaction) %>% dplyr::summarize(rmCond_reaction = ifelse(sum(modification == 'rmCond') , T, F)))

rmCond_sets <- spearman_MMandReg %>% ungroup() %>% filter(grepl('rmCond', modification)) %>% dplyr::select(reaction, modification)
# Take the rmCond reactions when they exist and the normal reactions otherwise

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


##### Replotting a few reactions for figures #####

# TPI - michaelis-menten-kinetics
load("shinyapp/reaction_data/r_1054plots.Rdata")

reactionInfo %>% filter(reaction == "r_1054")

saved_plot <- ggplotGrob(shiny_flux_data[["r_1054-rm"]][['plotChoices']][['Flux and species']])
panes <- which(saved_plot[["layout"]][,"name"] == "panel")
grid.draw(saved_plot)



ggsave(plot = shiny_flux_data[["r_1054-rm"]][['plotChoices']][['Flux and species']], "Figures/TPIrmm.pdf", height = 14, width = 10)


# PyK - comparison of michaelis-menten kinetics with non-significant and significant regulation
load("shinyapp/reaction_data/r_0962plots.Rdata")

ggsave(plot = shiny_flux_data[["r_0962-rm"]][['plotChoices']][['Flux and species']], "Figures/PYKrmm.pdf", height = 14, width = 10) # RMM
ggsave(plot = shiny_flux_data[["r_0962-rm-t_0292-act-mm"]][['plotChoices']][['Flux and species']], "Figures/PYKf6p.pdf", height = 14, width = 10) # F6P
ggsave(plot = shiny_flux_data[["r_0962-rm-t_0290-act-mm"]][['plotChoices']][['Flux and species']], "Figures/PYKfbp.pdf", height = 14, width = 10) # FBP

reactionInfo[reactionInfo$rMech %in% c("r_0962-rm", "r_0962-rm-t_0292-act-mm", "r_0962-rm-t_0290-act-mm"),]

rxn_fits[rxn_fits$rxn %in% c("r_0962-rm", "r_0962-rm-t_0292-act-mm", "r_0962-rm-t_0290-act-mm"),]


# OTCase - comparsion of michaelis-menten kinetics and versus alanine
load("shinyapp/reaction_data/r_0816plots.Rdata")

ggsave(plot = shiny_flux_data[["r_0816-rm_rmCond"]][['plotChoices']][["Flux comparison"]], "Figures/OTCasermm.pdf", height = 6, width = 10) # RMM
ggsave(plot = shiny_flux_data[["r_0816-rm-t_0461-inh-uncomp_rmCond"]][['plotChoices']][["Flux comparison"]], "Figures/OTCaseAla.pdf", height = 6, width = 10) # ALA
shiny_flux_data[["r_0816-rm_rmCond"]][['plotChoices']][[2]]
ala_rxns <- reactionInfo %>% filter(reaction == "r_0816") %>% filter(grepl('L-alanine', FullName))

rxn_fits[rxn_fits$rxn %in% c("r_0816-rm_rmCond", "r_0816-rm-t_0461-inh-uncomp_rmCond"),]
run_rxn[["r_0816-rm-t_0461-inh-uncomp_rmCond"]]


par_likelihood <- NULL
par_markov_chain <- NULL

parSubset <- param_set_list[parSetInfo$rx == "r_0816-rm-t_0461-inh-uncomp_rmCond"]

for(i in 1:length(parSubset)){
  
  par_likelihood <- rbind(par_likelihood, data.frame(sample = 1:param_run_info$n_samples[parSubset[[i]]$name$index], likelihood = parSubset[[i]]$lik, index = parSubset[[i]]$name$index))
  par_markov_chain <- rbind(par_markov_chain, parSubset[[i]]$MC)
}


#load(paste(c("FBGA_files/paramSets/", param_run_info$file[param_run_info$index == par_likelihood$index[1]]), collapse = ""))
#run_rxn <- run_summary[["r_0816-rm-t_0461-inh-uncomp_rmCond"]]

otcase_alanine_ki <- data.frame(ki = 2^par_markov_chain[,'t_0461'])
otcase_alanine_ki_mle <- otcase_alanine_ki[which.max(par_likelihood$likelihood),]
otcase_alanine_ki_exp <- 14.8e-3

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 20, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                       axis.ticks.x = element_line(color = "black", size = 1), axis.ticks.y = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1)
                       )

ggplot(otcase_alanine_ki, aes(x = ki)) + geom_bar(binwidth = 0.05) + 
  geom_vline(xintercept = otcase_alanine_ki_mle, color = "RED", size = 2) + geom_vline(xintercept = otcase_alanine_ki_exp, color = "BLUE", size = 2) +
  scale_x_log10("Affinity", breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1), labels = c("100uM", "1mM", "10mM", "100mM", "1M", "10M"), expand = c(0,0)) + barplot_theme + 
  scale_y_continuous("Counts", expand = c(0,0)) + geom_blank(aes(y=1.1*..count..), binwidth = 0.05, stat="bin") +
  geom_text(data = data.frame(label = c("MLE", "Measured"), x = 8, y = c(87.5, 80), color = c("RED", "BLUE")), aes(x = x, y = y, label = label, color = color), size = 10) +
  scale_color_identity()
ggsave("Figures/OTCaseAla_dist.pdf", width = 7, height = 7)



##### Generate figure summarizing metabolic leverage for a condition #####

#library(stringr)

#ML_natural_vs_all <- MLdata %>% dplyr::select(reaction, specie:conditions, median = contains('0.5')) %>% ungroup() %>% filter(grepl('^[PCN]', condition)) %>%
#  group_by(specie, reaction, condition) %>%
#  summarize(ALL = median[conditions == "ALL"], NAT = median[conditions == "NATURAL"])

#ML_natural_vs_all <- ML_natural_vs_all %>% mutate(Diff = abs(ALL - NAT)) %>% group_by(specie, reaction) %>% summarize(diffSum = sum(Diff)) %>%
#  mutate(rxID = substr(reaction, 1, 6))

#ML_natural_vs_all_rxSummary <- ML_natural_vs_all %>% group_by(rxID, specie) %>% summarize(avg_diffSum = mean(diffSum), nrxns = length(diffSum)) %>% ungroup() %>% arrange(-avg_diffSum) %>% 
#  group_by(rxID) %>% filter(nrxns == max(nrxns))

#ML_natural_vs_all %>% ungroup() %>% arrange(-diffSum)

MLdata$Type[grepl("r_0302", MLdata$reaction) & MLdata$specie == "isocitrate"] <- "product" # aconitase is split into two reactions without a meaasured cis-aconitate so isocitrate acts like a product in the first rxn
MLdata <- data.table(MLdata %>% filter(conditions == "NATURAL"))
metabolic_leverage_summary_plots("P0.05")

# Only looking at reactions that are well-fit
# Look at the best significant reaction form for reaction where CI overlap with flux carried is substantial (>50%)

adequate_fit_optimal_rxn_form <- union(intersect(optimal_rxn_form, fraction_flux_deviation$rxn[fraction_flux_deviation$"Interval Overlap" > 0.5]), intersect(optimal_rxn_form, rxn_fits$rxn[rxn_fits$parSpearman > 0.6]))

enzyme_control_source <- function(){
  
  ##### Associating enzyme metabolic leverage with transcription factors, thermodynamics ... #####
  # Look at the ML of enzymes to identify cases where:
  # A) Potential inducibility is high verus low
  # B) Potential inducibility varies based upon nutrient condition
  
  ML_inducibility <- MLdata[reaction %in% adequate_fit_optimal_rxn_form,,]
  ML_inducibility <- ML_inducibility[Type == "enzyme",,]
  ML_inducibility <- ML_inducibility[,list(ML = sum(get("0.5")), nenzyme = length(get("0.5"))), by = c("reaction", "condition")]
  
  # collapse across conditions to mean(ML) and SD(ML)
  ML_inducibility_summary <- ML_inducibility[,list(ML_mean = mean(ML), ML_sd = sd(ML), ML_min = min(ML), ML_max = max(ML), ML_range = max(ML)-min(ML), nenzyme = nenzyme[1]), by = "reaction"]
  ML_inducibility_summary[,CV := ML_sd/ML_mean] 
  ML_inducibility_summary$rxn <- substr(ML_inducibility_summary$reaction, 1, 6)
  ML_inducibility_summary[,ML_logit := log(ML_mean/(1-ML_mean))]
  
  rxn_meta_info <- read.delim("flux_cache/rxnParYeast.tsv")
  aligned_ML_meta <- rxn_meta_info[chmatch(ML_inducibility_summary$rxn, rxn_meta_info$ReactionID),]
  
  # Determine pathway E.C. numbers
  # If multiple E.C. identifiers exist, take the first one (this will be the identifier associated with the kegg R ID if there is one)
  aligned_ML_meta$EC <- sapply(aligned_ML_meta$EC, function(x){
    strsplit(x, split = ",")[[1]][1]
  })
  
  aligned_ML_meta$EC1 <- sapply(aligned_ML_meta$EC, function(x){strsplit(x, split = "\\.")[[1]][1]})
  aligned_ML_meta$EC2 <- sapply(aligned_ML_meta$EC, function(x){paste(strsplit(x, split = "\\.")[[1]][1:2], collapse = ".")})
  
  # Determine pathway-by-reaction
  
  aligned_ML_meta$pathname <- sapply(aligned_ML_meta$pathname, function(x){strsplit(x, split = "__")[[1]][1]})
  
  rxn_meta_path_prune <- melt(lapply(aligned_ML_meta$pathname, function(x){
    strsplit(x, split = "__")
  }))
  
  rxn_meta_path_pruned <- rxn_meta_path_prune[rxn_meta_path_prune$value %in% names(table(rxn_meta_path_prune$value))[table(rxn_meta_path_prune$value) >= 4 & table(rxn_meta_path_prune$value) <= 20],]
  rxn_meta_path_pruned$rxn <- aligned_ML_meta$ReactionID[rxn_meta_path_pruned$L1]
  rxn_meta_path_pruned <- rbind(rxn_meta_path_pruned, data.frame(value = "Misc", L2 = 1, L1 = NA, rxn = aligned_ML_meta$ReactionID[!(aligned_ML_meta$ReactionID %in% rxn_meta_path_pruned$rxn)]))
  rxn_pathways_cast <- acast(rxn_meta_path_pruned, formula = rxn ~ value, value.var = "L2", fill = 0)
  
  # Determine reaction reversibility
  aligned_ML_meta$reversibility <- ifelse(reversibleRx$modelBound[chmatch(aligned_ML_meta$ReactionID, reversibleRx$rx)] == "greaterEqual", "F", "T")
  
  # Regression of EC and pathway on metabolic leverage
  cat("\nEnzyme control vs. EC\n")
  print(anova(lm(ML_inducibility_summary$ML_logit ~ aligned_ML_meta$EC1))) # associate with E.C. number
  cat("\nEnzyme control vs. pathway\n")
  print(anova(lm(ML_inducibility_summary$ML_logit ~ rxn_pathways_cast))) # assocaite with pathway
  cat("\nEnzyme control vs. rxn reversibility\n")
  print(anova(lm(ML_inducibility_summary$ML_logit ~ aligned_ML_meta$reversibility))) # associate with reversibility
  
  regulatory_contingency <- table(regulated = c(1:nrow(ML_inducibility_summary)) %in% grep('act|inh', ML_inducibility_summary$reaction), reversible = aligned_ML_meta$reversibility)
  cat("\nRegulatory reaction vs. rxn reversibility\n")
  print(chisq.test(regulatory_contingency, simulate.p.value = T))
  #regulatory_contingency[2,2]/sum(regulatory_contingency[,2])
  #regulatory_contingency[2,1]/sum(regulatory_contingency[,1])
  
  ML_inducibility_summary$reversibility <- aligned_ML_meta$reversibility
  setkey(ML_inducibility_summary, 'ML_mean')
  ML_inducibility_summary$rxn <- factor(ML_inducibility_summary$rxn, levels = ML_inducibility_summary$rxn)
  
  # Association with trancription factor targets
  
  # convert reaction enzymes to common names
  rxn_enzyme_groups <- read.delim("./flux_cache/rxn_enzyme_groups.tsv")
  library("org.Sc.sgd.db")
  c2o <- toTable(org.Sc.sgdCOMMON2ORF)
  
  # convert from systematic names to common names and then generate a compact summary of genes with consecutive numbers
  ML_reaction_enzymes <- sapply(ML_inducibility_summary$rxn, function(x){
    commonSubset <- sort(c2o$gene_name[chmatch(unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction == x]), c2o$systematic_name)])
    commonSubsetDF <- data.frame(a = regmatches(commonSubset, regexpr('^[A-Z]{3}', commonSubset)), n = regmatches(commonSubset, regexpr('[0-9]+', commonSubset)))
    
    gene_name_compact = NULL
    for(an_a in unique(commonSubsetDF$a)){
      if(length(commonSubsetDF$n[commonSubsetDF$a == an_a]) == 1){
        gene_name_compact <- c(gene_name_compact, paste(an_a, commonSubsetDF$n[commonSubsetDF$a == an_a], sep = ""))
      }else{
        q_seq <- as.numeric(commonSubsetDF$n[commonSubsetDF$a == an_a])
        q_group <- rep(1:length(q_seq))
        for(q_el in 1:(length(q_seq)-1)){
          if(q_seq[q_el] + 1 == q_seq[q_el + 1]){
            q_group[q_el + 1] <- q_group[q_el]
          }
        }
        group_track <- NULL
        for(a_group in unique(q_group)){
          if(length(q_seq[q_group == a_group]) == 1){
            group_track <- c(group_track, q_seq[q_group == a_group])
          }else{
            group_track <- c(group_track, paste(q_seq[q_group == a_group][1], q_seq[q_group == a_group][length(q_seq[q_group == a_group])], sep = "-"))
          }
        }
        gene_name_compact <- c(gene_name_compact, paste(an_a, paste(group_track, collapse = ", ") , sep = ""))
      }
    }
    return(c(expanded = paste(commonSubset, collapse = ", "), collapsed = paste(gene_name_compact, collapse = ", ")))
  })
  ML_reaction_enzymes <- as.data.frame(t(ML_reaction_enzymes))
  
  ML_inducibility_summary$genes <- ML_reaction_enzymes$collapsed
  
  # Determine transcription factors regulating reaction subset
  
  ML_gene_summaries <- data.frame(systematic = unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% ML_inducibility_summary$rxn]),
                                  common = c2o$gene_name[chmatch(unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% ML_inducibility_summary$rxn]), c2o$systematic_name)])
  
  # Import matrix relating transcription factors to their targets (http://www.yeastract.com/generateregulationmatrix.php)
  # Interaction based on "DNA binding and expression evidence", all genes considered #
  # expression or affinity : 10% non-zero
  # affinity : 2% non-zero
  
  # connect every gene involved in a reaction with whether it is a TF target
  TF_indirect <- as.matrix(read.delim("./companionFiles/yeast_TF_regulation.csv", sep = ";", row.names = 1)) # direct and indirect
  TF_indirect <- TF_indirect[,colnames(TF_indirect) %in% ML_gene_summaries$common]
  
  TF_direct <- as.matrix(read.delim("./companionFiles/Yeast_TF_affinity.csv", sep = ";", row.names = 1)) # direct targets
  TF_direct <- TF_direct[,colnames(TF_direct) %in% ML_gene_summaries$common]
  
  if(all(colnames(TF_indirect) == colnames(TF_direct))){
    
    rxn2gene <- matrix(0, nrow = nrow(ML_inducibility_summary), ncol = ncol(TF_direct))
    rownames(rxn2gene) <- ML_inducibility_summary$rxn; colnames(rxn2gene) <- colnames(TF_direct)
    for(i in 1:nrow(ML_reaction_enzymes)){
      rxn2gene[i,colnames(rxn2gene) %in% strsplit(ML_reaction_enzymes$expanded[i], split = ", ")[[1]]] <- 1
    }  
    
  }else{
    stop("TF_direct and TF_indirect gene complements differ -> rxn2gene transformation needs to be modified")
  }
  
  # convert matrix from TF ~ Gene to TF ~ Rxn
  
  TF_indirect_byrxn <- TF_indirect %*% t(rxn2gene); TF_indirect_byrxn[TF_indirect_byrxn != 0] <- 1
  TF_direct_byrxn <- TF_direct %*% t(rxn2gene); TF_direct_byrxn[TF_direct_byrxn != 0] <- 1
  
  # Reduce the number of transcription factors to those with a role in steady-state metabolism #
  # from FIRE, generate a subset of TFs shaping transcription across these conditions (based upon Brauer data)
  FIRE_TFs <- c("Msn2p", "Msn4p", "Gcn4p", "Bas1p", "Cbf1p", "Mbp1p", "Swi4p")
  
  TF_indirect_byrxn <- TF_indirect_byrxn[rownames(TF_indirect_byrxn) %in% FIRE_TFs,]
  TF_direct_byrxn <- TF_direct_byrxn[rownames(TF_direct_byrxn) %in% FIRE_TFs,]
  
  TF_ML_assoc <- data.frame(TF = c(rownames(TF_indirect_byrxn), rownames(TF_direct_byrxn)), Effect = c(rep("indirect", times = nrow(TF_indirect_byrxn)),rep("direct", times = nrow(TF_direct_byrxn))) , p = NA)
  
  for(i in 1:nrow(TF_ML_assoc)){
    if(TF_ML_assoc$Effect[i] == "indirect"){
      refVec <- TF_indirect_byrxn[rownames(TF_indirect_byrxn) == TF_ML_assoc$TF[i],]
    }else{
      refVec <- TF_direct_byrxn[rownames(TF_direct_byrxn) == TF_ML_assoc$TF[i],]
    }
    if(length(unique(refVec)) == 1){next}
    
    TF_ML_assoc$p[i] <- wilcox.test(ML_inducibility_summary$ML_mean[refVec == 1], ML_inducibility_summary$ML_mean[refVec == 0], alternative = "two.sided")$p.value
  }
  
  library(qvalue)
  
  TF_ML_assoc$q <- qvalue(TF_ML_assoc$p)$q
  TF_ML_assoc <- TF_ML_assoc[TF_ML_assoc$q < 0.1,] # look at TFs with an FDR of less than 0.1
  TF_ML_assoc$label <- paste(sub('p$', '', TF_ML_assoc$TF), TF_ML_assoc$Effect, sep = "-")
  
  TFsigSubset <- rbind(TF_indirect_byrxn[chmatch(TF_ML_assoc$TF[TF_ML_assoc$Effect == "indirect"], rownames(TF_indirect_byrxn)),],
                       TF_direct_byrxn[chmatch(TF_ML_assoc$TF[TF_ML_assoc$Effect == "direct"], rownames(TF_direct_byrxn)),])
  rownames(TFsigSubset) <- TF_ML_assoc$label
  
  TF_effect_melt <- data.table(melt(TFsigSubset))
  setnames(TF_effect_melt, colnames(TF_effect_melt), c("TF", "rxn", "target"))
  TF_effect_melt$rxn <- factor(TF_effect_melt$rxn, levels = levels(TF_effect_melt$rxn))
  TF_effect_melt$TF <- factor(TF_effect_melt$TF, levels = rev(sort(unique(as.character(TF_effect_melt$TF)))))
  TF_effect_melt$ypos <- max(ML_inducibility_summary$ML_max) + as.numeric(TF_effect_melt$TF)/25
  TF_effect_melt$fillCol <- ifelse(TF_effect_melt$target == 1, "chocolate1", "aliceblue")  
  
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                         panel.background = element_rect(fill = "gray90"), legend.position = "top", 
                         axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                         axis.text = element_text(color = "black"), axis.text.x = element_text(size = 18, angle = 75, hjust = 1, vjust = 1),
                         panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                         axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
  )
  
  write.table(ML_inducibility_summary, file = "flux_cache/ML_inducibility_summary.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
  
  ggplot() + geom_pointrange(data = ML_inducibility_summary, aes(x = rxn, y = ML_mean, ymin = ML_min, ymax = ML_max), size = 2, color = "blue3") +
    scale_x_discrete("Reactions", breaks = ML_inducibility_summary$rxn, labels = ML_inducibility_summary$genes) + scale_y_continuous("Enyzme Metabolic Leverage", breaks = seq(0,0.8, by = 0.2), expand = c(0,0)) +
    geom_raster(data = TF_effect_melt, aes(x = rxn, y = ypos, fill = fillCol)) + 
    geom_text(data = TF_effect_melt[TF_effect_melt$rxn == ML_inducibility_summary$rxn[1],], aes(x = rxn, y = ypos, label = TF), hjust = 0, size = 7) +
    barplot_theme + scale_fill_identity() + scale_color_identity() + expand_limits(y = 0)
  
  ggsave("Figures/MLstrength.pdf", height = 10, width = 14)
  
  return(ML_inducibility_summary)
  
}


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

test <- ML_rxn_ternaryPoints %>% dplyr::filter(reversibility == "T")

summary(lm(ML_rxn_summary, formula = rxn_metabolite ~ reversibility))


color_key$Figure + 
  geom_point(data = test, aes(x = x, y = y), size = 7, shape = 21, fill = "BLACK") +
  geom_point(data = test, aes(x = x, y = y), size = 6, shape = 21, fill = "WHITE")



# Name point according to the metabolic layout

rxn_names <- c('r_1838' = 'HCS', 'r_0988' = 'SDH', 'r_0915' = 'PPAT', 'r_0042' = 'DAHP synthase', 'r_0886' = 'PFK', 'r_0718' = 'ME',
  'r_0468' = 'G5K', 'r_0962' = 'PyK', 'r_4040' = 'PPK', 'r_0225' = 'ATP-PRTase', 'r_0887' = 'S7P PFK', 'r_0888' = 'PGM',
  'r_0195' = 'TPS', 'r_0214' = 'ATCase', 'r_0310' = 'CBL', 'r_0816' = 'OTCase', 'r_0450' = 'ALD', 'r_0215' = 'AspK',
  'r_0491' = 'G3PDH', 'r_0250' = 'CPS')

ML_rxn_ternaryPoints_labels <- ML_rxn_ternaryPoints %>% filter(regulator != 0) %>% mutate(label = rxn_names[names(rxn_names) == rID])

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

all_affinities <- read.delim("flux_cache/metaboliteAffinities.tsv")
all_affinities <- all_affinities[!is.na(all_affinities$log10mean),]

all_affinities <- all_affinities[unname(unique(unlist(sapply(unique(reactionInfo$reaction), function(x){grep(x, all_affinities$reactions)})))),] #only looking at reactions of interest


### Using 95% confidence intervals for reaction parameters
yeast_only <- F

affinity_comparisons <- NULL

for(rxn in names(rxn_fit_params)){
  
  rxData <- data.frame(rxn = reactionInfo$reaction[reactionInfo$rMech == rxn], rxnForm = rxn, 
                       specie = rownames(rxn_fit_params[[rxn]]$param_interval[rownames(rxn_fit_params[[rxn]]$param_interval) != "keq",]), rxn_fit_params[[rxn]]$param_interval[rownames(rxn_fit_params[[rxn]]$param_interval) != "keq",], 
                       lit_mean = NA, lit_SD = NA)
  rxData <- rxData[grep('^t_', rownames(rxData)),]
  rxData$absoluteQuant <- as.logical(rxData$absoluteQuant)
  
  rxData <- rxData[chmatch(rxn_fit_params[[rxn]]$param_species$SubstrateID, rxData$specie),]
  rxData$speciesType <- rxn_fit_params[[rxn]]$param_species$Type
  rxData$speciesSubtype <- rxn_fit_params[[rxn]]$param_species$Subtype
  rxData$medianAbund <- rxn_fit_params[[rxn]]$param_species$medianAbund
  
  # Pass BRENDA data to reaction-from-specific parameter estimates
  
  rxAffinities <- all_affinities[grep(reactionInfo$reaction[reactionInfo$rMech == rxn], all_affinities$reactions),]
  
  for(a_spec_n in c(1:nrow(rxData))[rxData$absoluteQuant]){
    modSubset <- rxAffinities[rxAffinities$tID == rxData$specie[a_spec_n],]
    
    if(rxData$speciesType[a_spec_n] == "rct"){
      output <- modSubset[modSubset$speciesType == "substrate" & modSubset$isYeast == yeast_only,]
    }
    
    if(rxData$speciesType[a_spec_n] == "inh"){
      output <- modSubset[modSubset$speciesType == "regulator" & modSubset$modtype == "inh" & modSubset$isYeast == yeast_only,]
    }
    
    if(rxData$speciesType[a_spec_n] == "act"){
      output <- modSubset[modSubset$speciesType == "regulator" & modSubset$modtype == "act" & modSubset$isYeast == yeast_only,]
    }
    
    if(nrow(output) > 1){
      # if multiple E.C. IDs collide when matched to reactions with the same modifiers, choose the better supported
      output <- output[which.max(output$nQuant),]
    }
    
    if(nrow(output) == 1){
      rxData[a_spec_n, c("lit_mean", "lit_SD")] <- output[,c("log10mean", "sdOflog10")]
    } 
  }
  affinity_comparisons <- rbind(affinity_comparisons, rxData)
}


### Use literature mean & SD of Km as a comparison

### Looking at species with absolue quantification ###

absolute_quant <- affinity_comparisons[!is.na(affinity_comparisons$lit_mean),]
absolute_quant$lit_mean <- absolute_quant$lit_mean - 3 # convert from log mM to log M

boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), 
                       legend.position = "right", panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), 
                       axis.line = element_blank(), axis.text = element_text(color = "black"))

### only look at rm kinetics and significant modifications ###

rxn_of_interest <- reactionInfo$rMech[reactionInfo$form == "rm" & reactionInfo$modification == "" | !is.na(reactionInfo$Qvalue) & reactionInfo$Qvalue < 0.1]

ggplot(absolute_quant[absolute_quant$rxnForm %in% rxn_of_interest,], aes(x = lit_mean, y = log10(MLE), color = speciesType)) + geom_point(size = 1.3) + stat_abline(size = 2, alpha = 0.5) +
  boxplot_theme + stat_smooth(method = "lm", fill = "black") +
  scale_color_brewer("Metabolite Class", palette = "Set1", labels = c("Inhibitors", 'Substrates/\nProducts')) + 
  scale_x_continuous(expression("Literature" ~ log[10] ~ "Affinity (M)")) +
  scale_y_continuous(expression("Inferred" ~ italic("in vivo") ~ log[10] ~ "Affinity (M)"))
ggsave("Figures/affinity_match.pdf", height = 10, width = 12)

summary(lm(data = absolute_quant[absolute_quant$rxnForm %in% rxn_of_interest,], formula = "log10(MLE) ~ lit_mean"))


##

### Look at the differences between predicted and literature values versus the quality of fit ###

par_lb <- absolute_quant$X2.5.
par_ub <- absolute_quant$X97.5.
lit_mean <- absolute_quant$lit_mean
lit_SD <- absolute_quant$lit_SD

dist_overlap <- function(par_lb, par_ub, lit_mean, lit_SD){
  sapply(1:length(par_lb), function(i){
    diff(pnorm(c(log10(par_lb[i]), log10(par_ub[i])), mean = lit_mean[i], sd = lit_SD[i]))/(log10(par_ub[i]) - log10(par_lb[i]))
  })
}
absolute_quant$dist_overlap <- dist_overlap(absolute_quant$X2.5., absolute_quant$X97.5., absolute_quant$lit_mean, absolute_quant$lit_SD)
absolute_quant$param_range <- log2(absolute_quant$X97.5.) - log2(absolute_quant$X2.5.)


absolute_quant_deviation <- absolute_quant[,c('rxn', 'rxnForm', 'MLE', 'lit_mean', 'lit_SD', 'speciesType', 'speciesSubtype', 'dist_overlap', 'param_range')]
absolute_quant_deviation$MLE <- log10(absolute_quant_deviation$MLE)
absolute_quant_deviation$Interval_Overlap <- fraction_flux_deviation$"Interval Overlap"[chmatch(absolute_quant_deviation$rxnForm, fraction_flux_deviation$rxn)]
absolute_quant_deviation$Deviation <- absolute_quant_deviation$MLE - absolute_quant_deviation$lit_mean
absolute_quant_deviation <- absolute_quant_deviation[absolute_quant$param_range < 10 & !is.na(absolute_quant_deviation$dist_overlap),]

scatter_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), 
      legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"),
      strip.background = element_rect(fill = "cyan"), axis.text = element_text(color = "black")) 


ggplot(absolute_quant_deviation[absolute_quant_deviation$rxnForm %in% rxn_of_interest,], aes(x = Interval_Overlap, y = dist_overlap, color = factor(speciesType))) + 
  geom_point(size = 4) + stat_smooth(method = "lm", size = 2) + scale_color_brewer("Metabolite Class", palette = "Set1", labels = c("Inhibitors", 'Substrates/\nProducts')) +
  scatter_theme + scale_x_continuous("Fit Quality") + scale_y_continuous('Overlap of affinity distribution\nand literature values') + 
  ggtitle('Accuracy of metabolite affinity predictions\n increases with fit quality')
ggsave("Figures/KmConsistencyVersusFit.pdf", height = 14, width = 14)



summary(lm(data = absolute_quant_deviation[absolute_quant_deviation$rxnForm %in% rxn_of_interest,], formula = "dist_overlap ~ Interval_Overlap * factor(speciesType)"))

### Look at metabolism-wide occupancy ###

mw_occupancy <- data.table(affinity_comparisons)[,list(rxnForm = rxnForm, specie = specie, log10S_km = log10(2^medianAbund/MLE), log10occ = (2^medianAbund/(2^medianAbund + MLE)), type = speciesSubtype)]

common_mets <- sort(table(mw_occupancy$specie), decreasing = T)[1:6]

mw_occupancy$chosenSpecies <- "The Rest"
mw_occupancy$chosenSpecies[mw_occupancy$specie %in% names(common_mets)] <- mw_occupancy$specie[mw_occupancy$specie %in% names(common_mets)]

mw_occupancy$chosenSpecies[mw_occupancy$chosenSpecies != "The Rest"] <- 
  sapply(corrFile$SpeciesName[chmatch(mw_occupancy$chosenSpecies[mw_occupancy$chosenSpecies != "The Rest"], corrFile$SpeciesType)], function(x){
  strsplit(x, split = ' \\[')[[1]][1]
})

mw_occupancy$type[mw_occupancy$type %in% c("competitive", "noncompetitive")] <- "non(competitive)"
mw_occupancy$type <- factor(mw_occupancy$type)

mw_occupancy <- mw_occupancy[!mw_occupancy$chosenSpecies %in% "H+",]

mw_occupancy$chosenSpecies <- factor(mw_occupancy$chosenSpecies, levels = names(sort(table(mw_occupancy$chosenSpecies), decreasing = T)))

# subset mw_occupany if desired #
mw_occupancy <- mw_occupancy[mw_occupancy$rxnForm %in% rxn_of_interest,]


table(mw_occupancy$type, mw_occupancy$chosenSpecies)
spec_counts <- melt(table(mw_occupancy$type, mw_occupancy$chosenSpecies))
colnames(spec_counts) <- c("type", "chosenSpecies", "Counts")

# S / Km
ggplot() + geom_violin(data = mw_occupancy, aes(x = chosenSpecies, y = log10S_km), fill = "firebrick2", scale = "width") +
  facet_grid(type ~ .) + geom_text(data = spec_counts, aes(x = chosenSpecies, y = 0, label = Counts), color = "BLUE", size = 6) +
  scale_y_continuous(expression(frac(S, K[M])), breaks = seq(-3, 3, by = 1), labels = 10^seq(-3, 3, by = 1)) +
  boxplot_theme
ggsave("Figures/speciesSKm.pdf", width = 12, height = 14)

# S / S + Km
ggplot() + geom_violin(data = mw_occupancy, aes(x = chosenSpecies, y = log10occ), fill = "firebrick2", scale = "width") +
  facet_grid(type ~ .) + geom_text(data = spec_counts, aes(x = chosenSpecies, y = 0.5, label = Counts), color = "BLUE", size = 6) +
  scale_y_continuous(expression(frac(S, S + K[M])), breaks = seq(0, 1, by = 0.2), label = sapply(seq(0, 1, by = 0.2) * 100, function(x){paste(x, "%", sep = "")})) +
  boxplot_theme
ggsave("Figures/speciesOccupancy.pdf", width = 12, height = 14)



######## Compare Keq to Gibbs free energy ########

# rxForm defined keq: prod[S] - prod[P]/keq
# calculate median absolute concentration or bounds based upon reasonable concentrations of small molecules
# relate optimized keq to these values
# Q = prod[P]/prod[S]
# Q/keq



reaction_free_energy <- read.delim("companionFiles/cc_dG_python.tsv")

free_energy_table <- data.frame(RXform = names(rxn_fit_params), reaction = substr(names(rxn_fit_params), 1, 6), cc_dGr = NA, cc_dGrSD = NA, log_keq_LB = NA, 
                                log_keq_UB = NA, log_keq_MLE = NA, logQmedian = NA, intervalOverlap = NA)

free_energy_table <- free_energy_table[free_energy_table$reaction %in% reaction_free_energy$reaction,]

free_energy_table$formType = ifelse(free_energy_table$RXform %in% reactionInfo$rMech[reactionInfo$form == "rm" & reactionInfo$modification == ""], "MM", "modification")

free_energy_table$cc_dGr = reaction_free_energy$dGr[chmatch(free_energy_table$reaction, reaction_free_energy$reaction)]
free_energy_table$cc_dGrSD = reaction_free_energy$dGrSD[chmatch(free_energy_table$reaction, reaction_free_energy$reaction)]

free_energy_table[,c("log_keq_LB")] <- log(sapply(free_energy_table$RXform, function(x){rxn_fit_params[[x]]$param_interval[rownames(rxn_fit_params[[x]]$param_interval) == "keq", 1]}))
free_energy_table[,c("log_keq_UB")] <- log(sapply(free_energy_table$RXform, function(x){rxn_fit_params[[x]]$param_interval[rownames(rxn_fit_params[[x]]$param_interval) == "keq", 2]}))
free_energy_table[,c("log_keq_MLE")] <- log(sapply(free_energy_table$RXform, function(x){rxn_fit_params[[x]]$param_interval$MLE[rownames(rxn_fit_params[[x]]$param_interval) == "keq"]}))
free_energy_table[,c("logQmedian")] <- sapply(free_energy_table$RXform, function(x){as.numeric(rxn_fit_params[[x]]$param_interval$absoluteQuant[rownames(rxn_fit_params[[x]]$param_interval) == "keq"])})
free_energy_table[,c("intervalOverlap")] <- sapply(free_energy_table$RXform, function(x){fraction_flux_deviation$"Interval Overlap"[fraction_flux_deviation$rxn == x]})

free_energy_table <- free_energy_table[free_energy_table$cc_dGrSD < 10,] #remove reactions with highly uncertain CC predictions

# Find the best fitting form (which doesn't have a hypothetical regulator) #
core_form_bool <- c(1:nrow(reactionInfo)) %in% grep('t_metX|rmCond', reactionInfo$modification, invert = T)
for(rxn in unique(free_energy_table$reaction)){
  free_energy_table$formType[free_energy_table$RXform == reactionInfo$rMech[reactionInfo$reaction == rxn & core_form_bool][which.max(reactionInfo$ML[reactionInfo$reaction == rxn & core_form_bool])]] <- "Best Fit"
  }




boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), 
                       legend.position = "right", panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), 
                       axis.line = element_blank(), axis.text = element_text(color = "black"))


energy_fit <- free_energy_table#[free_energy_table$formType == "Best Fit",]
energy_weight <- energy_fit$intervalOverlap; energy_weight[energy_weight < 0] <- 0

freeEfit <- lm((energy_fit$logQmedian - energy_fit$log_keq_MLE)*(8.315 * (273 + 25) / 1000) ~ energy_fit$cc_dGr, weights = energy_weight)

weight_fit <- coef(freeEfit)
summary(freeEfit)$coef[2,4]

ggplot(free_energy_table[free_energy_table$formType == "Best Fit",], aes(x = cc_dGr, xmin = cc_dGr - 2*cc_dGrSD, xmax = cc_dGr + 2*cc_dGrSD,
                              y = (logQmedian - log_keq_MLE)*(8.315 * (273 + 25) / 1000),
                              ymin = (logQmedian - log_keq_LB)*(8.315 * (273 + 25) / 1000),
                              ymax = (logQmedian - log_keq_UB)*(8.315 * (273 + 25) / 1000))) +
  geom_errorbar(aes(color = intervalOverlap), size = 1) + geom_errorbarh(aes(color = intervalOverlap), size = 1) + boxplot_theme +
  geom_abline(intercept = weight_fit[1], slope = weight_fit[2], size = 2, color = "BLUE") +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(expression("Component Contribution" ~ Delta[G]^o)) + 
  scale_y_continuous(expression("Parameter Inference" ~ Delta[G])) +
  scale_color_gradient2("Par ~ FBA interval overlap", low = "green", mid = "green", high = "red", limits = c(0,1))


### Remap axes ###
free_energy_direct_comp_table <- data.frame(predicted = (free_energy_table$logQmedian - free_energy_table$log_keq_MLE)*(8.315 * (273 + 25) / 1000),
                                            free_energy_table[,c('cc_dGr', 'intervalOverlap', 'RXform', 'formType')])
free_energy_direct_comp_table$Deviation <- free_energy_direct_comp_table$predicted - free_energy_direct_comp_table$cc_dGr

free_energy_subset <- free_energy_direct_comp_table#[free_energy_direct_comp_table$formType == "Best Fit",]
reversible_rxn_match <- reversibleRx[chmatch(substr(free_energy_subset$RXform, 1, 6), reversibleRx$rx),]

free_energy_subset$is_reversible <- ifelse(reversible_rxn_match$modelBound == "reversible", T, F)

ggplot(free_energy_subset, aes(x = cc_dGr, y = predicted, color = intervalOverlap)) + geom_point() + facet_wrap(~ is_reversible, ncol = 1) +
  scale_x_continuous(expression("Component Contribution" ~ Delta[G]^o)) + 
  scale_y_continuous(expression("Parameter Inference" ~ Delta[G])) +
  geom_abline(intercept = 0, slope = 1) + geom_point(x = 0, y = 0, size = 5, col = "BLACK") +
  scale_color_gradient("Par ~ FBA interval overlap", low = "green", high = "red", limits = c(0,1)) + boxplot_theme

chisq.test(table(free_energy_subset$predicted < 0, !free_energy_subset$is_reversible))
chisq.test(table(free_energy_subset$predicted < 0, free_energy_subset$cc_dGr < 0 ))

deviation_design_matrix <- as.matrix(data.frame(int = 1, reversible_effect = ifelse(free_energy_subset$is_reversible, 1, 0), 
                                                irreversible_slope = free_energy_subset$intervalOverlap * !free_energy_subset$is_reversible,
                                                reversible_slope = free_energy_subset$intervalOverlap * free_energy_subset$is_reversible))

summary(lm(abs(free_energy_subset$Deviation) ~ deviation_design_matrix + 0))

ggplot(free_energy_subset, aes(x = intervalOverlap, y = abs(Deviation), col = factor(is_reversible))) + geom_point() + stat_smooth(method = "lm", size = 2) +
  scale_y_continuous(expression("abs( Parameter Inference" ~ Delta[G] - ~ 'Component Contribution' ~ Delta[G]^o ~ ")")) + 
  scale_x_continuous("Par ~ FBA interval overlap") + scale_color_brewer("Reversible\nReaction?", palette = "Set1") + boxplot_theme
ggsave("Figures/gibbsDeviation.pdf", width = 12, height = 14)




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