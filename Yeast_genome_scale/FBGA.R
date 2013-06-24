#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Dropbox/Rabino/git/FBA_SRH/Yeast_genome_scale")

library(reshape2) #for visualization at the end
library(nnls) #for non-negative regression used to fit kinetic parameters
library(ggplot2)
library(gplots)
source("FBA_lib.R")

options(stringsAsFactors = FALSE)

######## Import all of the species involved in a reaction and other reaction information ###########

load("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)

##### import list of metabolite abundances and rxn forms

chemostatInfo <- chemostatInfo[!(chemostatInfo$condition %in% c("p0.05H1", "p0.05H2")),] # the 25 chemostat conditions of interest and their actual growth rates


load("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/rxnf_formulametab.rdata") #### run Vito script which using the stoichiometric matrix and FBA_run files to construct a list of reaction information, including a mechanism.

##### Import fluxes
load("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/fluxSummaryQP.Rdata") #load flux through each reaction

##### Import enzyme abundances
enzyme_abund <- read.delim("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/relAbundMatrix.tsv")
#### load protLMfile.Rdata when rerun in order to get all of the genes that a non-unique peptide-set could correspond to
#measured_genes <- unlist(sapply(colnames(prot_abund_final), function(x){strsplit(x, "/")}))

##### Associate enzymes with pathways ######

kegg_enzyme_dict <- read.delim("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/yeastNameDict.tsv") # KEGG IDs relative to yeast gene name (1gene -> 1 KEGG id, multiple  mapping between KEGG and genes)
if(is.null(kegg_enzyme_dict$PATHWAY)){
  genes_to_pathways = read.delim("http://rest.kegg.jp/link/pathway/sce", header = FALSE); colnames(genes_to_pathways) <- c("gene", "pathwayCode")
  pathway_names = read.delim("http://rest.kegg.jp/list/pathway/sce", header = FALSE); colnames(pathway_names) <- c("pathwayCode", "pathway")
  genes_to_pathways$gene <- gsub('sce:', '', genes_to_pathways$gene)
  pathway_names$pathway <- sapply(pathway_names$pathway, function(x){strsplit(x, split = " - Sacc")[[1]][1]})
  pathway_names$pathway <- sapply(pathway_names$pathway, function(x){strsplit(x, split = " - yeast")[[1]][1]})
  pathway_match <- merge(genes_to_pathways, pathway_names)
  kegg_enzyme_dict$PATHWAY <- sapply(kegg_enzyme_dict$SYST, function(x){
    paste(pathway_match$pathway[pathway_match$gene == x], collapse = "__")
  })
  write.table(kegg_enzyme_dict, "../KEGGrxns/yeastNameDict.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}#generate a per-gene pathway annotation if one is not already generated


##### For each reaction in the consensus reconstruction, determine which pathways are associated with its enzymes ######

rxnFile <- read.delim('/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/rxn_yeast.tsv', stringsAsFactors = FALSE)
met_genes <- data.frame(reaction = unique(rxnFile$ReactionID), genes = NA, pathway = NA)

for(rxN in 1:length(met_genes[,1])){
  rxSubset <- rxnFile[rxnFile$ReactionID == met_genes$reaction[rxN],]
  met_genes$genes[rxN] <- paste(rxSubset$MetName[is.na(rxSubset$StoiCoef)], collapse = "/")
  gene_subset <- strsplit(paste(rxSubset$MetName[is.na(rxSubset$StoiCoef)], collapse = ":"), split = ":")[[1]]
  met_genes$pathway[rxN] <- paste(unique(strsplit(paste(kegg_enzyme_dict[kegg_enzyme_dict$SYST %in% gene_subset,]$PATHWAY, collapse = "__"), "__")[[1]]), collapse = "__")
}

###### Create a list containing model and experimental information for all reactions in the consensus reconstruction #####

rxnList_all <- rxnf

for(rxN in c(1:length(rxnList_all))){
  kegg_subset <- met_genes[met_genes$reaction == substr(names(rxnList_all),1,6)[rxN],]
  if(length(kegg_subset[,1]) == 0){next}
  rxnList_all[[rxN]]$pathway = kegg_subset$pathway
  rxnList_all[[rxN]]$genes = kegg_subset$genes
  rxnList_all[[rxN]]$enzymeAbund = enzyme_abund[enzyme_abund$Gene %in% strsplit(kegg_subset$genes, split = '[/:]')[[1]],]
}



####### Narrow the previous list to only rxns which carry flux #####

rxnList <- rxnf[substr(names(rxnf),1,6) %in% flux_summary$IDs$reactionID]

for(rxN in c(1:length(flux_summary$IDs[,1]))[grep('r_', flux_summary$IDs$reactionID)]){
  kegg_subset <- met_genes[met_genes$reaction == flux_summary$IDs$reactionID[rxN],]
  idx <- names(rxnList)[grep(flux_summary$IDs$reactionID[rxN],names(rxnList))]
  if(length(kegg_subset[,1]) == 0 | length(idx) ==0){next}
  
  for (entry in names(rxnList)[grep(flux_summary$IDs$reactionID[rxN],names(rxnList))]){
    rxnList[[entry]]$reaction = flux_summary$IDs$Name[rxN]
    rxnList[[entry]]$pathway = kegg_subset$pathway
    rxnList[[entry]]$genes = kegg_subset$genes
    rxnList[[entry]]$enzymeAbund = enzyme_abund[enzyme_abund$Gene %in% strsplit(kegg_subset$genes, split = '[/:]')[[1]],]
    rxnList_all[[entry]]$flux <- rxnList[[entry]]$flux <- flux_summary$ceulluarFluxes[rxN,]
  }
  
}
#save(rxnList_all, file = "all_rxnList.Rdata")


##### predict flux using linear regression using either metabolites and enzymes or their log - partition variance explained into that explained by metabolite and by enzymes ####

reaction_pred_log <- data.frame(nmetab = rep(NA, length(rxnList)), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, TSS = NA)
reaction_pred_linear <- data.frame(nmetab = rep(NA, length(rxnList)), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, TSS = NA)


cond_mapping <- data.frame(flux_cond = names(rxnList[[1]]$flux), 
                           enzyme_reordering = sapply(toupper(rownames(rxnList[[1]]$rxnMet)), 
                                                      function(x){
                                                        c(1:length(colnames(rxnList[[1]]$enzymeAbund)))[colnames(rxnList[[1]]$enzymeAbund) == x]})) 
cond_mapping$enzyme_cond = colnames(rxnList[[1]]$enzymeAbund)[cond_mapping$enzyme_reordering]
cond_mapping$met_cond = rownames(rxnList[[1]]$rxnMet)

# check if everything is correctly ordered
if (!all(toupper(cond_mapping$flux_cond) == cond_mapping$enzyme_cond & cond_mapping$enzyme_cond == cond_mapping$met_cond)){
  warning('There is a problem with the order of the conditions. (check cond_mapping)')
}

for(rxN in 1:length(rxnList)){
  reaction_pred_linear$nenz[rxN] <- reaction_pred_log$nenz[rxN] <- length(rxnList[[rxN]]$enzymeAbund[,1])
  reaction_pred_linear$nmetab[rxN] <- reaction_pred_log$nmetab[rxN] <- sum(!is.na(rxnList[[rxN]]$rxnMet))/25
  
  if(all(rxnList[[rxN]]$flux >= 0)){
    rxFlux <- log2(rxnList[[rxN]]$flux)
  } else if(all(rxnList[[rxN]]$flux <= 0)){
    rxFlux <- log2(-1*rxnList[[rxN]]$flux)
  } else{
    next
  } #if only forward flux consider its log, if only backwards flux consider the log of -1*flux, if the directionality changes then skip this rxn.
  rxFlux[!is.finite(rxFlux)] <- NA
  
  reaction_pred_log$nCond[rxN] <- reaction_pred_linear$nCond[rxN] <- sum(!is.na(rxFlux))
  
  if(reaction_pred_linear$nenz[rxN] != 0){
    enzyme_df = rxnList[[rxN]]$enzymeAbund[,cond_mapping$enzyme_reordering]
    
    lenzymes <- as.matrix(data.frame(t(enzyme_df))); colnames(lenzymes) <- rxnList[[rxN]]$enzymeAbund$Gene
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
    reaction_pred_linear$Fenz[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes))$F[1]
    reaction_pred_linear$varExplainedEnzy[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes))$Sum[1]
    reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux ~ enzymes))$Sum)  
  }  
  
  if(reaction_pred_linear$nmetab[rxN] != 0){
    ### only metabolites ###
    ### prediction using log measures ###
    reaction_pred_log$Fmetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$F[1]
    reaction_pred_log$varExplainedMetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$Sum[1]
    reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lmetabs))$Sum)
    
    ### prediction using linear measures ###
    reaction_pred_linear$Fmetab[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ metabs))$F[1]
    reaction_pred_linear$varExplainedMetab[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ metabs))$Sum[1]
    reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux ~ metabs))$Sum)
  }  
  
  if(reaction_pred_linear$nmetab[rxN] != 0 & reaction_pred_linear$nenz[rxN] != 0){
    ### both metabolites and enzymes ###
    reaction_pred_linear$varExplainedTotal[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux ~ enzymes + metabs))$Sum[1:2])
    reaction_pred_log$varExplainedTotal[rxN] <- sum(anova(lm(rxFlux ~ lenzymes + lmetabs))$Sum[1:2])
    
    reaction_pred_linear$varExplainedEither[rxN] <- sum(reaction_pred_linear$varExplainedMetab[rxN], reaction_pred_linear$varExplainedEnzy[rxN]) - reaction_pred_linear$varExplainedTotal[rxN]
    reaction_pred_linear$varExplainedEnzy[rxN] <- reaction_pred_linear$varExplainedEnzy[rxN] - reaction_pred_linear$varExplainedEither[rxN]
    reaction_pred_linear$varExplainedMetab[rxN] <- reaction_pred_linear$varExplainedMetab[rxN] - reaction_pred_linear$varExplainedEither[rxN]
    
    reaction_pred_log$varExplainedEither[rxN] <- sum(reaction_pred_log$varExplainedMetab[rxN], reaction_pred_log$varExplainedEnzy[rxN]) - reaction_pred_log$varExplainedTotal[rxN]
    reaction_pred_log$varExplainedEnzy[rxN] <- reaction_pred_log$varExplainedEnzy[rxN] - reaction_pred_log$varExplainedEither[rxN]
    reaction_pred_log$varExplainedMetab[rxN] <- reaction_pred_log$varExplainedMetab[rxN] - reaction_pred_log$varExplainedEither[rxN]
  }
}


qplot(reaction_pred_linear$Fmetab)
qplot(reaction_pred_linear$Fenz)

reaction_pred_log$varExplainedJointly <- NA; reaction_pred_log$varExplainedJointly[!is.na(reaction_pred_log$varExplainedEither) & reaction_pred_log$varExplainedEither < 0] <- -1*reaction_pred_log$varExplainedEither[!is.na(reaction_pred_log$varExplainedEither) & reaction_pred_log$varExplainedEither < 0]
reaction_pred_linear$varExplainedJointly <- NA; reaction_pred_linear$varExplainedJointly[!is.na(reaction_pred_linear$varExplainedEither) & reaction_pred_linear$varExplainedEither < 0] <- -1*reaction_pred_linear$varExplainedEither[!is.na(reaction_pred_linear$varExplainedEither) & reaction_pred_linear$varExplainedEither < 0]
reaction_pred_log$varExplainedEither[!is.na(reaction_pred_log$varExplainedEither) & reaction_pred_log$varExplainedEither < 0] <- NA
reaction_pred_linear$varExplainedEither[!is.na(reaction_pred_linear$varExplainedEither) & reaction_pred_linear$varExplainedEither < 0] <- NA

reaction_pred_summary_log <- data.frame(N = reaction_pred_log$nCond, metaboliteVarianceExplained = reaction_pred_log$varExplainedMetab/reaction_pred_log$TSS, enzymeVarianceExplained = reaction_pred_log$varExplainedEnzy/reaction_pred_log$TSS, 
                                        varianceAmbiguouslyExplained = reaction_pred_log$varExplainedEither/reaction_pred_log$TSS, varianceJointlyExplained = reaction_pred_log$varExplainedJointly/reaction_pred_log$TSS)
reaction_pred_summary_linear <- data.frame(N = reaction_pred_linear$nCond, metaboliteVarianceExplained = reaction_pred_linear$varExplainedMetab/reaction_pred_linear$TSS, enzymeVarianceExplained = reaction_pred_linear$varExplainedEnzy/reaction_pred_linear$TSS,
                                           varianceAmbiguouslyExplained = reaction_pred_linear$varExplainedEither/reaction_pred_linear$TSS, varianceJointlyExplained = reaction_pred_linear$varExplainedJointly/reaction_pred_linear$TSS)
rownames(reaction_pred_summary_log) <- rownames(reaction_pred_summary_linear) <- names(rxnList)


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
                       panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line = element_blank()) 

rxnPredictionPlot <- ggplot(reaction_pred_summary_plotter, aes(x = factor(index), y = value, fill = as.factor(variable), color = "black")) + facet_grid(modelType ~ .)
rxnPredictionPlot + geom_bar(binwidth = 1) + barplot_theme + geom_vline(aes(xintercept = 0), size = 0.5) + geom_hline(aes(yintercept = 0), size = 0.5) + 
  scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_y_continuous(name = "Fraction of variance explained", expand = c(0,0), limits = c(0,1)) +
  scale_fill_brewer("Prediction Method", palette = "Set2") + scale_color_identity()

#scale_fill_manual(values = c("enzymeVarianceExplained" = "sienna1", "metaboliteVarianceExplained" = "steelblue1", varianceAmbiguouslyExplained = "olivedrab3", varianceJointlyExplained = "red")) 

ggsave("varianceExplained.pdf", width = 20, height = 12)


### compare variance explained using NNLS vs LS vs rxn form
### flux vs predicted flux colored by condition
### metabolite abundance ~ DR colored by condition, faceted over metabolites
### enzyme abundance ~ DR colored by condition, faceted over enzymes
### import conditions list
### analyze metabolomics data - convert to the same DR in proteomics/flux





#### Determine which reaction have valid reaction mechanisms ####

reactionForms <- sapply(rxnList, function(x){ifelse(!is.null(x$rxnForm), x$rxnID, NA)})  
rxnList_form <- rxnList[substr(names(rxnList),1,6) %in% reactionForms]
n_c <- length(rxnList_form[[1]]$flux)

#### save rxnList_form so that this self-sufficient list can be thrown at the cluster ###


##### looking at subset of reactions where a kinetic form was produced because most/all substrates were ascertained ####

load("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/paramOptim.Rdata")

run_summary <- list() #### MCMC run output and formatted inputs

markov_pars <- list()
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 10000 #how many total markov samples are desired
markov_pars$burn_in <- 500 #how many initial samples should be skipped

markov_pars <- list()
markov_pars$sample_freq <- 5 #what fraction of markov samples are reported (this thinning of samples decreases sample autocorrelation)
markov_pars$n_samples <- 1000 #how many total markov samples are desired
markov_pars$burn_in <- 500 #how many initial samples should be skipped


run_summary$markov_pars <- markov_pars

#save(rxnList_form, cond_mapping, file = "paramOptim.Rdata")





###### import cluster parameter results #####

param_sets <- list.files('FBGA_files/paramSets/')

param_run_info <- data.frame(file = param_sets, index = 1:length(param_sets), runNum = sapply(param_sets, function(x){as.numeric(regmatches(x, regexec('[0-9]+', x))[[1]])}), n_samples = NA, sample_freq = NA, burn_in = NA)

param_set_list <- list()
for(an_index in param_run_info$index){
  load(paste(c("FBGA_files/paramSets/", param_run_info$file[an_index]), collapse = ""))
  param_run_info$n_samples[an_index] <- run_summary$markov_pars$n_samples
  param_run_info$sample_freq[an_index] <- run_summary$markov_pars$sample_freq
  param_run_info$burn_in[an_index] <- run_summary$markov_pars$burn_in
  param_set_list[[an_index]] <- run_summary
}

names(param_set_list[[1]])[-1]


## reconstruct environment of FBGA
# load("/home/vitoz/Dropbox/Rabino/RabinowitzData/Data_files/param_set_list.Rdata")
# names(param_set_list[[1]])[-1]
# 
# index = 1:length(param_set_list)
# param_run_info = data.frame(index)
# 
# for(an_index in param_run_info$index){
#   run_summary <- param_set_list[[an_index]]
#   param_run_info$n_samples[an_index] <- run_summary$markov_pars$n_samples
#   param_run_info$sample_freq[an_index] <- run_summary$markov_pars$sample_freq
#   param_run_info$burn_in[an_index] <- run_summary$markov_pars$burn_in
# }

all_rxns <- names(rxnList_form)
rxn_fits <- NULL
rxn_fit_params <- list()

for(arxn in all_rxns){
  #which runs contain the desired reaction
  rxn_found <- unlist(lapply(param_set_list, function(x){arxn %in% names(x)}))
  if(!all(rxn_found)){print(paste(c(arxn, "missing in", paste(param_run_info$file[!rxn_found], collapse = " & ")), collapse = " "))}
  if(all(!rxn_found)){print("skip"); next}
  
  par_likelihood <- NULL
  par_markov_chain <- NULL
  
  for(i in param_run_info$index[rxn_found]){
    run_rxn <- param_set_list[[i]][names(param_set_list[[i]]) == arxn][[1]]
    
    par_likelihood <- rbind(par_likelihood, data.frame(sample = 1:param_run_info$n_samples[i], likelihood = run_rxn$likelihood, index = param_run_info$index[i]))
    par_markov_chain <- rbind(par_markov_chain, run_rxn$markovChain)
  }
  
  if(var(par_likelihood$likelihood) < 10^-10 | all(run_rxn$metabolites == 1)){next} #skip underparameterized reactions - those with essentially no variation
  
  
  
  if(length(run_rxn$enzymes[1,]) != 0){
    flux_fit <- flux_fitting(run_rxn, par_markov_chain, par_likelihood) #compare flux fitted using the empirical MLE of parameters
    rxn_fits <- rbind(rxn_fits, data.frame(rxn = arxn, flux_fit$fit_summary))
    rxn_fit_params[[arxn]] <- flux_fit$param_interval  
  }
  
  
  pdf(file = paste(c("/home//vitoz/Dropbox/Rabino/git//vzcode//fbga_plots/", arxn, ".pdf"), collapse = ""), width = 20, height = 20)
  
  print(ggplot() + geom_violin(data = par_likelihood, aes(x = factor(index), y = likelihood), fill = "RED"))
  print(param_compare(run_rxn, par_markov_chain, par_likelihood))
  species_plots <- species_plot(run_rxn, flux_fit, chemostatInfo)
  print(species_plots[[1]])
  print(species_plots[[2]])
  print(species_plots[[3]])
  print(species_plots[[4]])
  dev.off()
  
}

save(param_set_list, file = "param_set_list.Rdata")


under_determined_rxnFits = rxn_fits[rxn_fits$residDF > 10,]
rownames(under_determined_rxnFits) <- under_determined_rxnFits$rxn
under_determined_rxnFits <- under_determined_rxnFits[,-c(1:2)]
under_determined_rxnFits_frac <- under_determined_rxnFits/under_determined_rxnFits$TSS
under_determined_rxnFits_frac <- under_determined_rxnFits_frac[,!(colnames(under_determined_rxnFits_frac) %in% c("TSS", "LS_met", "LS_enzyme"))]
under_determined_rxnFits_frac <- under_determined_rxnFits_frac[rowSums(is.na(under_determined_rxnFits_frac[,1:3])) == 0,]
under_determined_rxnFits_frac$bestFit = ifelse(under_determined_rxnFits_frac[,1] > under_determined_rxnFits_frac[,2], TRUE, "LS/NNLS max")
under_determined_rxnFits_frac$bestFit[under_determined_rxnFits_frac$bestFit == "TRUE"] <- ifelse(apply(under_determined_rxnFits_frac[under_determined_rxnFits_frac$bestFit == "TRUE",1:3], 1, which.max) == 1, "parametric fit", "parameteric fit > NNLS")

under_determined_rxnFits_frac <- under_determined_rxnFits_frac[order(under_determined_rxnFits_frac$bestFit),]
barplot_ordering <- 1:length(under_determined_rxnFits_frac[,1])
for(afit in unique(under_determined_rxnFits_frac$bestFit)){
  orderfit <- if(afit %in% c("parametric fit", "parameteric fit > NNLS")){"parametricFit"}else{"LS"}
  barplot_ordering[under_determined_rxnFits_frac$bestFit == afit] <- barplot_ordering[under_determined_rxnFits_frac$bestFit == afit][order(under_determined_rxnFits_frac[,colnames(under_determined_rxnFits_frac) == orderfit][under_determined_rxnFits_frac$bestFit == afit])]
}
under_determined_rxnFits_frac <- under_determined_rxnFits_frac[barplot_ordering,]

under_determined_rxnFits_frac$reaction <- factor(rownames(under_determined_rxnFits_frac), levels = rownames(under_determined_rxnFits_frac))
rxnFit_frac_melt <- melt(under_determined_rxnFits_frac, id.vars = c("reaction", "bestFit"))

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                       panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 90), axis.line = element_blank()) 


rxnFit_frac_plot <- ggplot(data = rxnFit_frac_melt, aes(x = reaction, y = value, fill = variable)) + facet_wrap(~bestFit, ncol = 1, scales = "free_x")
rxnFit_frac_plot + geom_bar(stat = "identity", position = "dodge") + barplot_theme + scale_x_discrete(name = "Reactions", expand = c(0,0)) +
  scale_y_continuous(name = "Fraction of variance explained", expand = c(0,0), limits = c(0,1)) + scale_fill_brewer("Prediction Method", palette = "Set2")
ggsave("parFitQuality.pdf", height = 20, width = 20)


# make another histogram that shows the explained variance with parametric vs nnls
rxn_fits$shareNNLS <- rxn_fits$NNLS / rxn_fits$TSS
rxn_fits$shareNNLS[rxn_fits$shareNNLS <0] <- 0
rxn_fits$shareParam <- rxn_fits$parametricFit/ rxn_fits$TSS
rxn_fits$shareParam [rxn_fits$shareParam <0] <- 0
hist(rxn_fits$shareNNLS,100)
hist(rxn_fits$shareParam,100)

rxn_fits$rxn[rxn_fits$shareParam >0.6]
rxn_fits <- rxn_fits[order(rxn_fits$shareParam),]
fil <- rxn_fits$shareParam > rxn_fits$shareNNLS
rxn_fits$rxn[fil] 
hist(rxn_fits$shareParam[fil])

plotRxn <- melt(rxn_fits[fil,c('rxn','shareParam','shareNNLS')])
colnames(plotRxn) <- c('reaction','fit','ESS')
plotRxn$fit <- as.character(plotRxn$fit)
plotRxn$reaction <- factor(plotRxn$reaction ,levels=unique(plotRxn$reaction))
plotRxn$fit[plotRxn$fit == 'shareParam'] <- 'parametric'
plotRxn$fit[plotRxn$fit == 'shareNNLS'] <- 'NNLS'

ggplot(data=plotRxn,aes(x=reaction,y=ESS))+
  geom_bar(stat="identity",alpha=1,aes(fill=fit), position = 'dodge')+
  theme(axis.text.x = element_text(angle=90,vjust=.50))


ggsave('ESShist.png',width=9,height=5)













