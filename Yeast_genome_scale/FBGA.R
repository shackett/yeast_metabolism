#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2) #for visualization at the end
library(nnls) #for non-negative regression used to fit kinetic parameters
library(ggplot2)
library(gplots)
source("FBA_lib.R")

options(stringsAsFactors = FALSE)

######## Import all of the species involved in a reaction and other reaction information ###########

load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)
chemostatInfo <- chemostatInfo[!(chemostatInfo$condition %in% c("p0.05H1", "p0.05H2")),] # the 25 chemostat conditions of interest and their actual growth rates

##### Import list of metabolite abundances and rxn forms

load("rxnf_formulametab.rdata") #### run Vito script which using the stoichiometric matrix and FBA_run files to construct a list of reaction information, including a mechanism.


##### Import fluxes
load("fluxSummaryQP.Rdata") #load flux through each reaction


##### Import enzyme abundances
enzyme_abund <- read.delim("../ChemicalSpeciesQuant/Proteomics/proteinAbundance.tsv")
rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
enzyme_abund <- enzyme_abund[,c(16:20, 1:5, 11:15, 6:10, 21:25)] #reorder proteins so that they are the same order as fluxes and metabolites (otherwise a warning will be issued later)


##### Associate enzymes with pathways ######

kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv") # KEGG IDs relative to yeast gene name (1gene -> 1 KEGG id, multiple  mapping between KEGG and genes)
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


##### For each reaction in the consensus reconstruction, determine which pathways and genes are associated ######

rxnFile <- read.delim('rxn_yeast.tsv', stringsAsFactors = FALSE)
met_genes <- data.frame(reaction = unique(rxnFile$ReactionID), genes = NA, pathway = NA)

for(rxN in 1:length(met_genes[,1])){
  rxSubset <- rxnFile[rxnFile$ReactionID == met_genes$reaction[rxN],]
  gene_subset <- strsplit(paste(rxSubset$MetName[is.na(rxSubset$StoiCoef)], collapse = ":"), split = ":")[[1]]
  if(length(gene_subset) == 0){
    met_genes$genes[rxN] <- ""
  }else{
    met_genes$genes[rxN] <- paste(unique(sort(gene_subset)), collapse = "/")
  }
  
  met_genes$pathway[rxN] <- paste(unique(strsplit(paste(kegg_enzyme_dict[kegg_enzyme_dict$SYST %in% gene_subset,]$PATHWAY, collapse = "__"), "__")[[1]]), collapse = "__")
}

###### Create a list containing model and experimental information for all reactions in the consensus reconstruction #####

rxnList_all <- rxnf

for(rxN in c(1:length(rxnList_all))){
  kegg_subset <- met_genes[met_genes$reaction == substr(names(rxnList_all),1,6)[rxN],]
  if(length(kegg_subset[,1]) == 0){next}
  
  rxnList_all[[names(rxnList_all)[rxN]]]$pathway = kegg_subset$pathway
  rxnList_all[[names(rxnList_all)[rxN]]]$genes = kegg_subset$genes
  
  rxnList_all[[names(rxnList_all)[rxN]]]$enzymeAbund = enzyme_abund[unlist(sapply(strsplit(kegg_subset$genes, split = '[/:]')[[1]], function(x){grep(x, rownames(enzyme_abund))})),]
  
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
    rxnList[[entry]]$enzymeAbund = enzyme_abund[unlist(sapply(strsplit(kegg_subset$genes, split = '[/:]')[[1]], function(x){grep(x, rownames(enzyme_abund))})),]
    rxnList_all[[entry]]$flux <- rxnList[[entry]]$flux <- flux_summary$cellularFluxes[rxN,]
  }
  
}
save(rxnList_all, file = "all_rxnList.Rdata")


### ensure that the ordering of conditions is the same ###

cond_mapping <- data.frame(standard = chemostatInfo$condition, e = colnames(enzyme_abund), m = rownames(rxnList[[1]]$rxnMet), f = colnames(flux_summary$cellularFluxes))

if (!all(toupper(cond_mapping$flux_cond) == cond_mapping$enzyme_cond & cond_mapping$enzyme_cond == cond_mapping$met_cond)){
  warning('There is a problem with the order of the conditions. (check cond_mapping)')
}



##### predict flux using linear regression using either metabolites and enzymes or their log - partition variance explained into that explained by metabolite and by enzymes ####

reaction_pred_log <- data.frame(nmetab = rep(NA, length(grep('rm$', names(rxnList)))), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, varExplainedJointly = NA, TSS = NA)
reaction_pred_linear <- data.frame(nmetab = rep(NA, length(grep('rm$', names(rxnList)))), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, varExplainedJointly = NA, TSS = NA)

for(rxN in grep('rm$', names(rxnList))){
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
    enzyme_df = rxnList[[rxN]]$enzymeAbund
    
    lenzymes <- as.matrix(data.frame(t(enzyme_df))); rownames(lenzymes) <- colnames(enzyme_abund)
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


hist(reaction_pred_linear$Fmetab)
hist(reaction_pred_linear$Fenz)


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


### remove the fraction of variance accounted for by either metabolites or enzymes from each of them ####

#reaction_pred_summary_log$metaboliteVarianceExplained[!is.na(reaction_pred_summary_log$varianceJointlyExplained)] <- reaction_pred_summary_log$metaboliteVarianceExplained[!is.na(reaction_pred_summary_log$varianceJointlyExplained)] - reaction_pred_summary_log$varianceJointlyExplained[!is.na(reaction_pred_summary_log$varianceJointlyExplained)]
#reaction_pred_summary_log$enzymeVarianceExplained[!is.na(reaction_pred_summary_log$varianceJointlyExplained)] <- reaction_pred_summary_log$enzymeVarianceExplained[!is.na(reaction_pred_summary_log$varianceJointlyExplained)] - reaction_pred_summary_log$varianceJointlyExplained[!is.na(reaction_pred_summary_log$varianceJointlyExplained)]
                                                   
#reaction_pred_summary_linear$metaboliteVarianceExplained[!is.na(reaction_pred_summary_linear$varianceJointlyExplained)] <- reaction_pred_summary_linear$metaboliteVarianceExplained[!is.na(reaction_pred_summary_linear$varianceJointlyExplained)] - reaction_pred_summary_linear$varianceJointlyExplained[!is.na(reaction_pred_summary_linear$varianceJointlyExplained)]
#reaction_pred_summary_linear$enzymeVarianceExplained[!is.na(reaction_pred_summary_linear$varianceJointlyExplained)] <- reaction_pred_summary_linear$enzymeVarianceExplained[!is.na(reaction_pred_summary_linear$varianceJointlyExplained)] - reaction_pred_summary_linear$varianceJointlyExplained[!is.na(reaction_pred_summary_linear$varianceJointlyExplained)]


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


### compare variance explained using NNLS vs LS vs rxn form
### flux vs predicted flux colored by condition
### metabolite abundance ~ DR colored by condition, faceted over metabolites
### enzyme abundance ~ DR colored by condition, faceted over enzymes
### import conditions list
### analyze metabolomics data - convert to the same DR in proteomics/flux





#### Determine which reaction have valid reaction mechanisms ####

reactionForms <- sapply(rxnList, function(x){ifelse(!is.null(x$rxnForm), x$listEntry, NA)})  
rxnList_form <- rxnList[names(rxnList) %in% reactionForms]
n_c <- nrow(chemostatInfo)

#### save rxnList_form so that this self-sufficient list can be thrown at the cluster ###

# chunks to break rxnList into
chunk_size <- 100
chunk_assignment <- data.frame(set = names(rxnList_form), chunk = c(rep(1:floor(length(rxnList_form)/chunk_size), each = chunk_size), rep(ceiling(length(rxnList_form)/chunk_size), length(rxnList_form) %% chunk_size)))
print(paste("The number of parameter chunks is", ceiling(length(rxnList_form)/chunk_size), "submit this parameter when batch submitting processes in FluxOptim.sh"))

save(rxnList_form, cond_mapping, chunk_assignment, file = "paramOptim.Rdata")


##### looking at subset of reactions where a kinetic form was produced because most/all substrates were ascertained ####






###### import cluster parameter results #####

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
    
parSetInfo <- as.data.frame(t(sapply(param_set_list, function(x){x$name})))
parSetInfo$ML <- sapply(param_set_list, function(x){max(x$lik)})

reactionInfo <- data.frame(rMech = names(rxnList_form), reaction = sapply(names(rxnList_form), function(x){substr(x, 1, 6)}),
  form = sapply(names(rxnList_form), function(x){substr(x, 8, 9)}), modification = sapply(names(rxnList_form), function(x){sub(substr(x, 1, 10), '', x)}),
  npar = sapply(rxnList_form, function(x){nrow(x$enzymeAbund) + nrow(x$rxnFormData) + 1}))

reactionInfo$ML <- sapply(reactionInfo$rMech, function(x){max(parSetInfo$ML[parSetInfo$rx == x])})

reactionInfo$changeP <- NA

for(rx in c(1:nrow(reactionInfo))[reactionInfo$form != "rm" | reactionInfo$modification != ""]){
  rxn_eval <- reactionInfo[rx,]
  rxn_ref <- reactionInfo[reactionInfo$reaction == rxn_eval$reaction & reactionInfo$form == "rm" & reactionInfo$modification == "",]
  likDiff <- rxn_eval$ML - rxn_ref$ML
  
  if(rxn_eval$npar == rxn_ref$npar){
    
    reactionInfo$changeP[rx] <- 1/(exp(likDiff) + 1)
    
  }else{
    
    reactionInfo$changeP[rx] <- 1 - pchisq(2*likDiff, rxn_eval$npar - rxn_ref$npar)
    
  }
}

### Identify reaction form modification which significantly improve the likelihood function ###

library(qvalue)
Qthresh <- 0.1

reactionInfo$Qvalue <- NA
reactionInfo$Qvalue[!is.na(reactionInfo$changeP)] <- qvalue(reactionInfo$changeP[!is.na(reactionInfo$changeP)])$q

### reduce reactions of interest to primary forms and reparameterizations ###

reactionInfo <- reactionInfo[is.na(reactionInfo$Qvalue) | (!is.na(reactionInfo$Qvalue) & reactionInfo$Qvalue < Qthresh),]

rxn_fits <- NULL
rxn_fit_params <- list()

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
  
  
  if(var(par_likelihood$likelihood) < 10^-10 | all(run_rxn$metabolites == 1)){next} #skip underparameterized reactions - those with essentially no variation
  
  
  
  if(length(run_rxn$enzymes[1,]) != 0){
    flux_fit <- flux_fitting(run_rxn, par_markov_chain, par_likelihood) #compare flux fitted using the empirical MLE of parameters
    rxn_fits <- rbind(rxn_fits, data.frame(rxn = arxn, flux_fit$fit_summary))
    rxn_fit_params[[arxn]] <- flux_fit$param_interval  
  }
  
  
  pdf(file = paste(c("FBGA_files/outputPlots/", arxn, ".pdf"), collapse = ""), width = 20, height = 20)
  
  reaction_info_FBGA(rxnName)
  
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
save(rxn_fit_params, file = "paramDist.Rdata")

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

 
  










