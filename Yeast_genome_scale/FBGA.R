#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2) #for visualization at the end
library(data.table)
library(nnls) #for non-negative regression used to fit kinetic parameters
library(ggplot2)
library(gplots)
source("FBA_lib.R")

options(stringsAsFactors = FALSE)

######## Import all of the species involved in a reaction and other reaction information ###########

load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)
chemostatInfo <- chemostatInfo[!(chemostatInfo$condition %in% c("p0.05H1", "p0.05H2")),] # the 25 chemostat conditions of interest and their actual growth rates
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


##### Import enzyme abundances
#enzyme_abund <- read.delim("../ChemicalSpeciesQuant/Proteomics/proteinAbundance.tsv")
#rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
#enzyme_abund <- enzyme_abund[,c(16:20, 1:5, 11:15, 6:10, 21:25)] #reorder proteins so that they are the same order as fluxes and metabolites (otherwise a warning will be issued later)


##### Associate enzymes with pathways ######

kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv") # KEGG IDs relative to yeast gene name (1gene -> 1 KEGG id, multiple  mapping between KEGG and genes)
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
  
  if(rxnList[[idx[1]]]$genes %in% kegg_enzyme_dict$SYST){
    geneInfo <- kegg_enzyme_dict[chmatch(rxnList[[idx[1]]]$genes[rxnList[[idx[1]]]$genes %in% kegg_enzyme_dict$SYST], kegg_enzyme_dict$SYST),]
  }else{geneInfo <- NA}
  
  rxFlux <- as.data.frame(flux_summary$total_flux_cast[rownames(flux_summary$total_flux_cast) == valid_rxns[rxN],,])
  
  #Add the rxns with lacking fluxes to the rmCondList
  if (any(rxFlux$standardQP == 0)){
    rmCondList <- rbind(rmCondList ,data.frame(rxn = valid_rxns[rxN],
                       cond = paste(rownames(rxFlux)[rxFlux$standardQP == 0],collapse=';'),
                       source = 'zeroFluxes', nZero = sum(rxFlux$standardQP == 0)))
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

rmCondList <- rmCondList[rmCondList$nZero < 10,]

for (i in 1:nrow(rmCondList)){
  rxn <- rmCondList$rxn[i]
  for (entry in names(rxnList)[grep(rxn,names(rxnList))]){
    nEntry <- paste(entry,'_rmCond',sep='')
    rxnList[[nEntry]] <- rxnList[[entry]]
    conds <- strsplit(rmCondList$cond[i],';')[[1]]
    conds <- c(conds,toupper(conds))
    rxnList[[nEntry]]$flux <- rxnList[[nEntry]]$flux[!rownames(rxnList[[nEntry]]$flux) %in% conds,]
    rxnList[[nEntry]]$rxnMet <- rxnList[[nEntry]]$rxnMet[!rownames(rxnList[[nEntry]]$rxnMet) %in% conds,]
    rxnList[[nEntry]]$all_species_SD <- rxnList[[nEntry]]$all_species_SD[!rownames(rxnList[[nEntry]]$all_species_SD) %in% conds,]
    rxnList[[nEntry]]$enzymeComplexes <- rxnList[[nEntry]]$enzymeComplexes[,!colnames(rxnList[[nEntry]]$enzymeAbund) %in% conds]
  
    rxnList[[nEntry]]$listEntry <- paste(rxnList[[nEntry]]$listEntry, '_rmCond',sep='')
    
    }
}





### ensure that the ordering of conditions is the same ###

cond_mapping <- data.frame(standard = chemostatInfo$condition, e = colnames(rxnList[[1]]$enzymeComplexes), m = rownames(rxnList[[1]]$rxnMet), f = rownames(rxnList[[1]]$flux))

if (!all(toupper(cond_mapping$flux_cond) == cond_mapping$enzyme_cond & cond_mapping$enzyme_cond == cond_mapping$met_cond)){
  warning('There is a problem with the order of the conditions. (check cond_mapping)')
}

#### Determine which reaction have valiad reaction mechanisms - all of them as of now ####

reactionForms <- sapply(rxnList, function(x){ifelse(!is.null(x$rxnForm), x$listEntry, NA)})  
rxnList_form <- rxnList[names(rxnList) %in% reactionForms]


rxnList_form <- rxnList_form[order(names(rxnList_form))] # order alpha-numerically so that the workload is spread out more evenly during optimization

#### save rxnList_form so that this self-sufficient list can be thrown at the cluster ###

# chunks to break rxnList into
chunk_size <- 100
chunk_assignment <- data.frame(set = names(rxnList_form), chunk = c(rep(1:floor(length(rxnList_form)/chunk_size), each = chunk_size), rep(ceiling(length(rxnList_form)/chunk_size), length(rxnList_form) %% chunk_size)))
print(paste("The number of parameter chunks is", ceiling(length(rxnList_form)/chunk_size), "submit this parameter when batch submitting processes in FluxOptim.sh"))

save(rxnList_form, cond_mapping, chunk_assignment, file = "paramOptim.Rdata")




##### predict flux using linear regression using either metabolites and enzymes or their log - partition variance explained into that explained by metabolite and by enzymes ####

reaction_pred_log <- data.frame(nmetab = rep(NA, length(grep('rm$', names(rxnList)))), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, varExplainedJointly = NA, TSS = NA)
reaction_pred_linear <- data.frame(nmetab = rep(NA, length(grep('rm$', names(rxnList)))), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, varExplainedJointly = NA, TSS = NA)

for(rxN in grep('rm$', names(rxnList))){
  reaction_pred_linear$nenz[rxN] <- reaction_pred_log$nenz[rxN] <- length(rxnList[[rxN]]$enzymeComplexes[,1])
  reaction_pred_linear$nmetab[rxN] <- reaction_pred_log$nmetab[rxN] <- sum(!is.na(rxnList[[rxN]]$rxnMet))/25
  
  if(all(rxnList[[rxN]]$flux$standardQP >= 0)){
    rxFlux <- log2(rxnList[[rxN]]$flux$standardQP)
  } else if(all(rxnList[[rxN]]$flux$standardQP <= 0)){
    rxFlux <- log2(-1*rxnList[[rxN]]$flux$standardQP)
  } else{
    next
  } #if only forward flux consider its log, if only backwards flux consider the log of -1*flux, if the directionality changes then skip this rxn.
  rxFlux[!is.finite(rxFlux)] <- NA
  
  reaction_pred_log$nCond[rxN] <- reaction_pred_linear$nCond[rxN] <- sum(!is.na(rxFlux))
  
  if(reaction_pred_linear$nenz[rxN] != 0){
    enzyme_df = rxnList[[rxN]]$enzymeComplexes
    
    lenzymes <- as.matrix(data.frame(t(enzyme_df)))
    enzymes <- 2^lenzymes
  }else{lenzymes <- enzymes <- NULL}
  if(reaction_pred_linear$nmetab[rxN] != 0){
    metab_df = rxnList[[rxN]]$rxnMet
    
    lmetabs <- as.matrix(data.frame(metab_df[,colSums(!is.na(metab_df) != 0) != 0])); colnames(lmetabs) <- colnames(rxnList[[rxN]]$rxnMet)[colSums(!is.na(rxnList[[rxN]]$rxnMet) != 0) != 0]
    metabs <- 2^lmetabs
  }else{lmetabs <- metabs <- NULL}
  
  if(reaction_pred_linear$nenz[rxN] == 0 & reaction_pred_linear$nmetab[rxN] == 0){
    next
  }
  
  if(reaction_pred_linear$nenz[rxN] != 0){
    ### only enzymes ###
    ### prediction using log measures ###
    reaction_pred_log$Fenz[rxN] <- anova(lm(rxFlux ~ lenzymes))$F[1]
    reaction_pred_log$varExplainedEnzy[rxN] <- anova(lm(rxFlux ~ lenzymes))$Sum[1]
    reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lenzymes))$Sum)
    
    ### prediction using linear measures ###
    reaction_pred_linear$Fenz[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes))$F[1]
    reaction_pred_linear$varExplainedEnzy[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes))$Sum[1]
    reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes))$Sum)  
  }  
  
  if(reaction_pred_linear$nmetab[rxN] != 0){
    ### only metabolites ###
    ### prediction using log measures ###
    reaction_pred_log$Fmetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$F[1]
    reaction_pred_log$varExplainedMetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$Sum[1]
    reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lmetabs))$Sum)
    
    ### prediction using linear measures ###
    reaction_pred_linear$Fmetab[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ metabs))$F[1]
    reaction_pred_linear$varExplainedMetab[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ metabs))$Sum[1]
    reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux$standardQP ~ metabs))$Sum)
  }  
  
  if(reaction_pred_linear$nmetab[rxN] != 0 & reaction_pred_linear$nenz[rxN] != 0){
    ### both metabolites and enzymes ###
    reaction_pred_linear$varExplainedTotal[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes + metabs))$Sum[1:2])
    if(reaction_pred_linear$varExplainedTotal[rxN] < max(reaction_pred_linear$varExplainedMetab[rxN], reaction_pred_linear$varExplainedEnzy[rxN])){
      reaction_pred_linear$varExplainedTotal[rxN] <- max(reaction_pred_linear$varExplainedMetab[rxN], reaction_pred_linear$varExplainedEnzy[rxN])  
    }
    
    ### Variance explained in a full model versus reduced ones
    
    if(reaction_pred_linear$varExplainedTotal[rxN] > sum(reaction_pred_linear$varExplainedMetab[rxN] + reaction_pred_linear$varExplainedEnzy[rxN])){
      ### add variance explained in complete model to new class - jointly described
      reaction_pred_linear$varExplainedJointly[rxN] <- reaction_pred_linear$varExplainedTotal[rxN] - sum(reaction_pred_linear$varExplainedMetab[rxN] + reaction_pred_linear$varExplainedEnzy[rxN])
    }else{
      ### some variance is equally accounted for by metabolites or enzymes and should be pulled out and removed from both enzymes and mets
      reaction_pred_linear$varExplainedEither[rxN] <- sum(reaction_pred_linear$varExplainedMetab[rxN] + reaction_pred_linear$varExplainedEnzy[rxN]) - reaction_pred_linear$varExplainedTotal[rxN]
      reaction_pred_linear$varExplainedMetab[rxN] <- reaction_pred_linear$varExplainedMetab[rxN] - reaction_pred_linear$varExplainedEither[rxN]
      reaction_pred_linear$varExplainedEnzy[rxN] <- reaction_pred_linear$varExplainedEnzy[rxN] - reaction_pred_linear$varExplainedEither[rxN]
    }
    
    reaction_pred_log$varExplainedTotal[rxN] <- sum(anova(lm(rxFlux ~ lenzymes + lmetabs))$Sum[1:2])
    if(is.na(reaction_pred_log$varExplainedTotal[rxN])){next}
    
    if(reaction_pred_log$varExplainedTotal[rxN] < max(reaction_pred_log$varExplainedMetab[rxN], reaction_pred_log$varExplainedEnzy[rxN])){
      reaction_pred_log$varExplainedTotal[rxN] <- max(reaction_pred_log$varExplainedMetab[rxN], reaction_pred_log$varExplainedEnzy[rxN])  
    }
    
    if(reaction_pred_log$varExplainedTotal[rxN] > sum(reaction_pred_log$varExplainedMetab[rxN] + reaction_pred_log$varExplainedEnzy[rxN])){
      ### add variance explained in complete model to new class - jointly described
      reaction_pred_log$varExplainedJointly[rxN] <- reaction_pred_log$varExplainedTotal[rxN] - sum(reaction_pred_log$varExplainedMetab[rxN] + reaction_pred_log$varExplainedEnzy[rxN])
    }else{
      ### some variance is equally accounted for by metabolites or enzymes and should be pulled out and removed from both enzymes and mets
      reaction_pred_log$varExplainedEither[rxN] <- sum(reaction_pred_log$varExplainedMetab[rxN] + reaction_pred_log$varExplainedEnzy[rxN]) - reaction_pred_log$varExplainedTotal[rxN]
      reaction_pred_log$varExplainedMetab[rxN] <- reaction_pred_log$varExplainedMetab[rxN] - reaction_pred_log$varExplainedEither[rxN]
      reaction_pred_log$varExplainedEnzy[rxN] <- reaction_pred_log$varExplainedEnzy[rxN] - reaction_pred_log$varExplainedEither[rxN]
    }
    
  }
}


#hist(reaction_pred_linear$Fmetab)
#hist(reaction_pred_linear$Fenz)

reaction_pred_summary_log <- data.frame(N = reaction_pred_log$nCond, metaboliteVarianceExplained = reaction_pred_log$varExplainedMetab/reaction_pred_log$TSS, enzymeVarianceExplained = reaction_pred_log$varExplainedEnzy/reaction_pred_log$TSS, 
                                        varianceAmbiguouslyExplained = reaction_pred_log$varExplainedEither/reaction_pred_log$TSS, varianceJointlyExplained = reaction_pred_log$varExplainedJointly/reaction_pred_log$TSS)

reaction_pred_summary_linear <- data.frame(N = reaction_pred_linear$nCond, metaboliteVarianceExplained = reaction_pred_linear$varExplainedMetab/reaction_pred_linear$TSS, enzymeVarianceExplained = reaction_pred_linear$varExplainedEnzy/reaction_pred_linear$TSS,
                                           varianceAmbiguouslyExplained = reaction_pred_linear$varExplainedEither/reaction_pred_linear$TSS, varianceJointlyExplained = reaction_pred_linear$varExplainedJointly/reaction_pred_linear$TSS)

rownames(reaction_pred_summary_log) <- rownames(reaction_pred_summary_linear) <- names(rxnList)[grep('rm$', names(rxnList))]


residualDF <- reaction_pred_log$nCond - reaction_pred_log$nmetab - reaction_pred_log$nenz
residualDF <- residualDF[reaction_pred_summary_log$N >= 15]

reaction_pred_summary_log <- reaction_pred_summary_log[reaction_pred_summary_log$N >= 15,-1]
reaction_pred_summary_linear <- reaction_pred_summary_linear[reaction_pred_summary_linear$N >= 15,-1]

reaction_pred_summary_log <- reaction_pred_summary_log[order(apply(reaction_pred_summary_log, 1, sum, na.rm = TRUE)),]
reaction_pred_summary_linear <- reaction_pred_summary_linear[order(apply(reaction_pred_summary_linear, 1, sum, na.rm = TRUE)),]




reaction_pred_summary_plotter <- rbind(data.frame(modelType = "logFlux ~ logMetab + logEnzyme", melt(data.frame(index = c(1:length(reaction_pred_summary_log[,1])), reaction_pred_summary_log), id.vars = "index")),
                                       data.frame(modelType = "Flux ~ Metab + Enzyme", melt(data.frame(index = c(1:length(reaction_pred_summary_linear[,1])), reaction_pred_summary_linear), id.vars = "index")))

reaction_pred_summary_plotter <- reaction_pred_summary_plotter[!is.na(reaction_pred_summary_plotter$value),]

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                       panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), panel.margin = unit(1.5, "lines"), axis.text = element_text(color = "BLACK"),
                       strip.background = element_rect(fill = "chocolate1"), strip.text = element_text(vjust = 0.6)) 

rxnPredictionPlot <- ggplot(reaction_pred_summary_plotter, aes(x = factor(index), y = value, fill = as.factor(variable), color = "black")) + facet_grid(modelType ~ .)
rxnPredictionPlot + geom_bar(stat ="identity", width=0.75) + barplot_theme + geom_vline(aes(xintercept = 0), size = 0.5) + geom_hline(aes(yintercept = 0), size = 0.5) + 
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_y_continuous(name = "Fraction of variance explained", expand = c(0,0), limits = c(0,1)) +
  scale_fill_brewer("Prediction Method", palette = "Set2") + scale_color_identity()

#scale_fill_manual(values = c("enzymeVarianceExplained" = "sienna1", "metaboliteVarianceExplained" = "steelblue1", varianceAmbiguouslyExplained = "olivedrab3", varianceJointlyExplained = "red")) 

ggsave("varianceExplained.pdf", width = 20, height = 12)



##@##@##@###@###@##@##@##@###@###@##@##@##@###@###@

######## Import cluster parameter results #########
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
  
  for(an_index in param_run_info$index[param_run_info$chunkNum == a_chunk]){ #temporarely only loading run number 1 (1 tenth of data)
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

### Pull out relevent reaction information - such as reaction ID, number of parameters and maximum likelihood ###

load('flux_cache/metaboliteTables.RData')
npc <- ncol(metSVD$v)

parSetInfo <- as.data.frame(t(sapply(param_set_list, function(x){x$name})))
parSetInfo$ML <- sapply(param_set_list, function(x){max(x$lik)})

# npar = number of enzymes (kcat) + number of metabolites (Km) + hill coefficients (which are not 1) 
# + hypothetical metabolite (trend governed by npc significant PCs) + Keq

reactionInfo <- data.frame(rMech = names(rxnList_form), reaction = sapply(names(rxnList_form), function(x){substr(x, 1, 6)}),
                           form = sapply(names(rxnList_form), function(x){substr(x, 8, 9)}), modification = sapply(names(rxnList_form), function(x){sub(substr(x, 1, 10), '', x)}),
                           npar = sapply(rxnList_form, function(x){nrow(x$enzymeAbund) + nrow(x$rxnFormData) + sum(x$rxnFormData$Hill != 1) + ifelse(any(x$rxnFormData$SubstrateID == "t_metX"), npc, 0) + 1}))

reactionInfo$ML <- sapply(reactionInfo$rMech, function(x){max(parSetInfo$ML[parSetInfo$rx == x])})

reactionInfo$changeP <- NA

for(rx in c(1:nrow(reactionInfo))[reactionInfo$form != "rm" | !(reactionInfo$modification %in% c("", "rmCond"))]){
  rxn_eval <- reactionInfo[rx,]
  if(length(grep('rmCond', rxn_eval$modification)) == 0){
    rxn_ref <- reactionInfo[reactionInfo$reaction == rxn_eval$reaction & reactionInfo$form == "rm" & reactionInfo$modification == "",]
  }else{
    rxn_ref <- reactionInfo[reactionInfo$reaction == rxn_eval$reaction & reactionInfo$form == "rm" & reactionInfo$modification == "rmCond",]
  }
  
  likDiff <- rxn_eval$ML - rxn_ref$ML
  
  if(rxn_eval$npar == rxn_ref$npar){
    
    reactionInfo$changeP[rx] <- 1/(exp(likDiff) + 1)
    
  }else{
    
    reactionInfo$changeP[rx] <- 1 - pchisq(2*likDiff, rxn_eval$npar - rxn_ref$npar)
    
  }
}

### Identify reaction form modification which significantly improve the likelihood function ###
### comparisons are relative to reversible michaelis-menten kinetics ###

library(qvalue)
Qthresh <- 0.1

reactionInfo$Qvalue <- NA
reactionInfo$Qvalue[!is.na(reactionInfo$changeP)] <- qvalue(reactionInfo$changeP[!is.na(reactionInfo$changeP)])$q

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

for(rxN in grep('act|inh', reactionInfo$modification)){
  
  regName <- substr(reactionInfo$modification[rxN], 1, 6)
  rxRegs <- rxnList_form[[rxN]]$rxnFormData[rxnList_form[[rxN]]$rxnFormData$Type != "rct",]
  
  reactionInfo$Name[rxN] <- paste(sub('(^)([a-z])', '\\U\\2', rxRegs$Subtype[rxRegs$SubstrateID == regName], perl = T), 
                                  ifelse(rxRegs$Type[rxRegs$SubstrateID == regName] == "act", "activation", "inhibition"),
                                  "by", unname(rxnList_form[[rxN]]$metNames[names(rxnList_form[[rxN]]$metNames) == regName]))
  
}

reactionInfo$Name[grep('rmCond', reactionInfo$modification)] <- sapply(reactionInfo$Name[grep('rmCond', reactionInfo$modification)], function(x){paste(x, "(zero flux reactions removed)")})

reactionInfo$Name <- mapply(function(x,y){paste(x,y)}, x = reactionInfo$signifCode, y = reactionInfo$Name)



### Base reaction specific information ###

rxToPW <- NULL

for(rxN in grep('rm$', names(rxnList_form))){
  
  RXannot <- rxnList_form[[rxN]]$pathway
  GENEannot <- rxnList_form[[rxN]]$geneInfo
  
  if(RXannot == "" | is.na(RXannot)){
    if(all(is.na(GENEannot))){
      pathways <- "Not annotated"
      RELannot <- NA
    }else{
      RELannot <- GENEannot$PATHWAY
      if(RELannot == "" |is.na(RELannot)){
      print(rxN)  
      }
    }
  }else{
    RELannot <- RXannot
  }
  
  if(!is.na(RELannot)){
    pathways <- strsplit(RELannot, split = '__')[[1]]
  }
  
  rxToPW <- rbind(rxToPW, data.frame(rxN = rxN, rID = rxnList_form[[rxN]]$rxnID, reactionName = rxnList_form[[rxN]]$reaction, pathway = pathways))
    
}  

rxToPW <- rbind(rxToPW, data.frame(unique(rxToPW[,1:3]), pathway = "ALL REACTIONS"))

reactionInfo$FullName <- mapply(function(x,y){paste(x, y, sep = " - ")}, x = rxToPW$reactionName[chmatch(reactionInfo$reaction, rxToPW$rID)], y = reactionInfo$Name)
  



pathwaySet <- sort(table(rxToPW$pathway), decreasing = T)
pathwaySet <- data.frame(pathway = names(pathwaySet), members = unname(pathwaySet), display = paste(names(pathwaySet), ' (', unname(pathwaySet), ')', sep = ""))


#### Generate reaction plots and summaries ####

shiny_flux_data <- list()
rxn_fits <- NULL
rxn_fit_params <- list()
fraction_flux_deviation <- NULL
PLdata <- NULL


for(arxn in reactionInfo$rMech){
  
  par_likelihood <- NULL
  par_markov_chain <- NULL
  
  parSubset <- param_set_list[parSetInfo$rx == arxn]
  
  for(i in 1:length(parSubset)){
    
    par_likelihood <- rbind(par_likelihood, data.frame(sample = 1:param_run_info$n_samples[parSubset[[i]]$name$index], likelihood = parSubset[[i]]$lik, index = parSubset[[i]]$name$index))
    par_markov_chain <- rbind(par_markov_chain, parSubset[[1]]$MC)
  }
  
  load(paste(c("FBGA_files/paramSets/", param_run_info$file[param_run_info$index == par_likelihood$index[1]]), collapse = ""))
  run_rxn <- run_summary[[arxn]]
  
  ### A couple of catches for poorly defined reactions or if optimization grossly fails ### 
  
  if(sum(is.finite(par_likelihood$likelihood)) == 0){
    print(paste(arxn, "has zero valid parameter sets")); next
  } # all parameter sets have a zero likelihood (e.g. SD = 0)
  
  if(var(par_likelihood$likelihood) < 10^-10 | all(run_rxn$metabolites == 1)){
    print(paste(arxn, "does not vary in likelihood, or has too many unmeasured species")); next
  } #skip underparameterized reactions - those with essentially no variation
  
  
  
  
  ### Determine how closely fluxes predicted from a parametric form match fluxes determined via constraint based modeling ###
  
  flux_fit <- flux_fitting(run_rxn, par_markov_chain, par_likelihood) #compare flux fitted using the empirical MLE of parameters
  rxn_fits <- rbind(rxn_fits, data.frame(rxn = arxn, flux_fit$fit_summary))
  rxn_fit_params[[arxn]] <- flux_fit$param_interval
  rxn_fit_params[[arxn]]$absoluteQuant <- ifelse(rownames(rxn_fit_params[[arxn]]) %in% names(run_rxn$rxnSummary$originMet)[run_rxn$rxnSummary$originMet == "abs"], T, F)
    
  
  ### Take the fractional flux departure and dot-product between FBA and parametric vector
  
  vector_match <- data.frame(rxn = arxn, FFD = 1 - sum(abs(flux_fit$fitted_flux$fitted - (run_rxn$flux$FVAmin + run_rxn$flux$FVAmax)/2))/sum(abs(run_rxn$flux)),
                             dotProduct = sum(flux_fit$fitted_flux$fitted/sqrt(sum((flux_fit$fitted_flux$fitted)^2)) * (run_rxn$flux$FVAmin + run_rxn$flux$FVAmax)/2/sqrt(sum(((run_rxn$flux$FVAmin + run_rxn$flux$FVAmax)/2)^2))))
  vector_match$angle <- acos(vector_match$dotProduct) * 180/pi
  
  fraction_flux_deviation <- rbind(fraction_flux_deviation, vector_match)
  
  
  ### Generate plots which show reaction information, species variation, flux fitting ... ###
  
  shiny_flux_data[[arxn]]$reactionInfo <- reaction_info_FBGA(rxnName) # Reaction information
  
  shiny_flux_data[[arxn]]$plotChoices$Likelihood <- likViolin(par_likelihood, run_summary$markov_pars) # Log-likelihoods of each markov chain
  
  species_plots <- species_plot(run_rxn, flux_fit, chemostatInfo)
  
  shiny_flux_data[[arxn]]$plotChoices <- append(species_plots, shiny_flux_data[[arxn]]$plotChoices)
  
  reaction_plots <- reactionProperties()
  PLdata <- rbind(PLdata, reaction_plots$PL_summary)
  shiny_flux_data[[arxn]]$plotChoices <- append(species_plots, reaction_plots[!(names(reaction_plots) %in% c("PL_summary"))])
  
  
  if("t_metX" %in% run_rxn$kineticPars$modelName){
    shiny_flux_data[[arxn]]$plotChoices <- append(shiny_flux_data[[arxn]]$plotChoices, hypoMetTrend(run_rxn, metSVD, tab_boer))
  }
  
  
}












# significant or default reaction forms

pathway_plot_list <- list()
for(pw in pathwaySet$display){
  # iterate through pathways and plot pathway-level figures
  pathway_plot_list[[pw]] <- pathwayPlots(pw)
  }

# order reactions according to the fit of the best reaction form (currently using DotProduct)




#### Save lists which will be processed by Shiny app ####

print(object.size((shiny_flux_data)), units = "auto")

save(pathwaySet, rxToPW, reactionInfo, pathway_plot_list, shiny_flux_data, file = "shinyapp/shinyData.Rdata")





reactionInfo
shiny_flux_data
rxToPW
pathwaySet


rxToPW
fraction_flux_deviation
rxn_fits




######### Compare Km values found through optimization to those from literature ################

all_affinities <- read.delim("flux_cache/metaboliteAffinities.tsv")

### Using 95% confidence intervals for reaction parameters
rxn_fit_params



### Use literature mean & SD of Km as a comparison

