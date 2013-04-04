#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2) #for visualization at the end
library(nnls) #for non-negative regression used to fit kinetic parameters
library(ggplot2)
library(gplots)

options(stringsAsFactors = FALSE)

######## Import all of the species involved in a reaction and other reaction information ###########

load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)

##### import list of metabolite abundances and rxn forms

chemostatInfo <- chemostatInfo[!(chemostatInfo$condition %in% c("p0.05H1", "p0.05H2")),] # the 25 chemostat conditions of interest and their actual growth rates


load("rxnf_formulametab.rdata") #### run Vito script which using the stoichiometric matrix and FBA_run files to construct a list of reaction information, including a mechanism.

##### Import fluxes
load("fluxSummaryQP.Rdata") #load flux through each reaction

##### Import enzyme abundances
enzyme_abund <- read.delim("../ChemicalSpeciesQuant/Proteomics/run_Deg/relAbundMatrix.tsv")
#### load protLMfile.Rdata when rerun in order to get all of the genes that a non-unique peptide-set could correspond to
#measured_genes <- unlist(sapply(colnames(prot_abund_final), function(x){strsplit(x, "/")}))

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


##### For each reaction in the consensus reconstruction, determine which pathways are associated with its enzymes ######

rxnFile <- read.delim('rxn_yeast.tsv', stringsAsFactors = FALSE)
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
  kegg_subset <- met_genes[met_genes$reaction == names(rxnList_all)[rxN],]
  if(length(kegg_subset[,1]) == 0){next}
  
  rxnList_all[[names(rxnList_all)[rxN]]]$pathway = kegg_subset$pathway
  rxnList_all[[names(rxnList_all)[rxN]]]$genes = kegg_subset$genes
  rxnList_all[[names(rxnList_all)[rxN]]]$enzymeAbund = enzyme_abund[enzyme_abund$Gene %in% strsplit(kegg_subset$genes, split = '[/:]')[[1]],]
  }



####### Narrow the previous list to only rxns which carry flux #####

rxnList <- rxnf[names(rxnf) %in% flux_summary$IDs$reactionID]

for(rxN in c(1:length(flux_summary$IDs[,1]))[grep('r_', flux_summary$IDs$reactionID)]){
  kegg_subset <- met_genes[met_genes$reaction == flux_summary$IDs$reactionID[rxN],]
  if(length(kegg_subset[,1]) == 0){next}
  
  rxnList[[flux_summary$IDs$reactionID[rxN]]]$reaction = flux_summary$IDs$Name[rxN]
  rxnList[[flux_summary$IDs$reactionID[rxN]]]$pathway = kegg_subset$pathway
  rxnList[[flux_summary$IDs$reactionID[rxN]]]$genes = kegg_subset$genes
  rxnList[[flux_summary$IDs$reactionID[rxN]]]$enzymeAbund = enzyme_abund[enzyme_abund$Gene %in% strsplit(kegg_subset$genes, split = '[/:]')[[1]],]
  rxnList_all[[flux_summary$IDs$reactionID[rxN]]]$flux <- rxnList[[flux_summary$IDs$reactionID[rxN]]]$flux <- flux_summary$ceulluarFluxes[rxN,]
}
#save(rxnList_all, file = "all_rxnList.Rdata")


##### predict flux using linear regression using either metabolites and enzymes or their log - partition variance explained into that explained by metabolite and by enzymes ####

reaction_pred_log <- data.frame(nmetab = rep(NA, length(rxnList)), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, TSS = NA)
reaction_pred_linear <- data.frame(nmetab = rep(NA, length(rxnList)), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, TSS = NA)


cond_mapping <- data.frame(flux_cond = names(rxnList[[1]]$flux), met_reordering = c(11:15, 1:5, 6:10, 16:20, 21:25), enzyme_reordering = sapply(toupper(chemostatInfo$condition), function(x){c(1:length(colnames(rxnList[[1]]$enzymeAbund)))[colnames(rxnList[[1]]$enzymeAbund) == x]})) #remove this when metabolite conditions are remapped onto flux/protein ones prior to this script
cond_mapping$met_cond = rownames(rxnList[[1]]$rxnMet)[cond_mapping$met_reordering] #remove this when new metabolite abundance are introduced
cond_mapping$enzyme_cond = colnames(rxnList[[1]]$enzymeAbund)[cond_mapping$enzyme_reordering]

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
    metab_df = rxnList[[rxN]]$rxnMet[cond_mapping$met_reordering,]
    
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
  scale_fill_manual(values = c("enzymeVarianceExplained" = "sienna1", "metaboliteVarianceExplained" = "steelblue1", varianceAmbiguouslyExplained = "olivedrab3", varianceJointlyExplained = "red")) + scale_color_identity()
  

ggsave("varianceExplained.pdf", width = 20, height = 12)


### compare variance explained using NNLS vs LS vs rxn form
### flux vs predicted flux colored by condition
### metabolite abundance ~ DR colored by condition, faceted over metabolites
### enzyme abundance ~ DR colored by condition, faceted over enzymes
### import conditions list
### analyze metabolomics data - convert to the same DR in proteomics/flux





#### Determine which reaction have valid reaction mechanisms ####

reactionForms <- sapply(rxnList, function(x){ifelse(!is.null(x$rxnForm), x$rxnID, NA)})  
rxnList_form <- rxnList[names(rxnList) %in% reactionForms]
n_c <- length(rxnList_form[[1]]$flux)

#### save rxnList_form so that this self-sufficient list can be thrown at the cluster ###


##### looking at subset of reactions where a kinetic form was produced because most/all substrates were ascertained ####

load("paramOptim.Rdata")

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



########### Functions ##########

par_draw <- function(updates){
  #### update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ####
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
      } else if(kineticParPrior$distribution[par_n] == "unif"){
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
  
  predOcc <- model.matrix(occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*enzyme_abund #occupany of enzymes * relative abundance of enzymes
  
  flux_fit <- nnls(enzyme_activity, flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  fit_resid_error <- sqrt(mean((flux_fit$resid - mean(flux_fit$resid))^2))
  
  sum(dnorm(flux, flux_fit$fitted, fit_resid_error, log = TRUE))
  
  }
current_pars -> proposed_params

############# Body ###########

for(rxN in 1:length(rxnList_form)){
  
  rxnSummary <- rxnList_form[[rxN]]
  
  occupancyEq <- as.formula(paste("~", sub(paste(paste("E", rxnSummary$rxnID, sep = "_"), " \\* ", paste("V", rxnSummary$rxnID, sep = "_"), sep = ""), "1", rxnSummary$rxnForm)[2], sep = " "))
  
  ### Create a data.frame describing the relevent parameters for the model ###
  kineticPars <- data.frame(rel_spec = c(rxnSummary$enzymeAbund[,1], colnames(rxnSummary$rxnMet)), 
  SpeciesType = c(rep("Enzyme", times = length(rxnSummary$enzymeAbund[,1])), rep("Metabolite", times = length(colnames(rxnSummary$rxnMet)))), modelName = NA, commonName = NA, formulaName = NA, measured = NA)
  kineticPars$formulaName[kineticPars$SpeciesType == "Enzyme"] <- paste("E", rxnSummary$rxnID, sep = "_")
  kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){rxnSummary$metsID2tID[names( rxnSummary$metsID2tID) == x]}))
  kineticPars$commonName[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){rxnSummary$metNames[names(rxnSummary$metNames) == x]}))
  kineticPars$commonName[kineticPars$SpeciesType == "Enzyme"] <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  kineticPars$formulaName[kineticPars$SpeciesType == "Metabolite"] <- paste("K", rxnSummary$rxnID, kineticPars$modelName[kineticPars$SpeciesType == "Metabolite"], sep = "_")
  
  all_species <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, occupancyEq)) != 0}) | kineticPars$SpeciesType == "Enzyme",]
  
  kineticPars <- kineticPars[sapply(kineticPars$formulaName, function(ele_used){length(grep(ele_used, occupancyEq)) != 0}),] #remove species which dont appear in the reaction equation
  kineticPars <- rbind(kineticPars, c("keq", "keq", NA, NA, paste("Keq", rxnSummary$rxnID, sep = ""), NA))
  
  ### Create a matrix containing the metabolites and enzymes 
  enzyme_abund <- t(rxnSummary$enzymeAbund[,cond_mapping$enzyme_reordering]); colnames(enzyme_abund) <- kineticPars$rel_spec[kineticPars$SpeciesType == "Enzyme"]
  met_abund <- rxnSummary$rxnMet[cond_mapping$met_reordering,]
  met_abund <- met_abund[,colnames(met_abund) %in% kineticPars$rel_spec]
  
  if(length(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]) <= 1){
    met_abund <- data.frame(met_abund)
    colnames(met_abund) <- kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"]
    kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- !all(is.na(met_abund))
    }else{
      kineticPars$measured[kineticPars$SpeciesType == "Metabolite"] <- unname(sapply(kineticPars$rel_spec[kineticPars$SpeciesType == "Metabolite"], function(x){(apply(is.na(met_abund), 2, sum) == 0)[names((apply(is.na(met_abund), 2, sum) == 0)) == x]}))
      }
  
  
  ### set missing data to invariant across conditions
  met_abund[,!as.logical(kineticPars$measured[kineticPars$rel_spec %in% colnames(met_abund)])] <- 1
  met_abund <- met_abund^2
  colnames(met_abund) <- unname(sapply(colnames(met_abund), function(x){kineticPars$modelName[kineticPars$rel_spec == x]}))
  
  enzyme_abund <- enzyme_abund^2
  
  flux <- rxnSummary$flux/median(abs(rxnSummary$flux[rxnSummary$flux != 0])) #flux, scaled to a prettier range
  
  #### write flux parameters to a list ####
  run_summary[[names(rxnList_form)[rxN]]]$metabolites <- met_abund
  run_summary[[names(rxnList_form)[rxN]]]$enzymes <- enzyme_abund
  run_summary[[names(rxnList_form)[rxN]]]$flux <- flux
  run_summary[[names(rxnList_form)[rxN]]]$occupancyEq <- occupancyEq
  run_summary[[names(rxnList_form)[rxN]]]$all_species <- all_species
  run_summary[[names(rxnList_form)[rxN]]]$kineticPars <- kineticPars
  run_summary[[names(rxnList_form)[rxN]]]$rxnSummary <- rxnSummary
  
  
  
  kineticParPrior <- data.frame(distribution = rep(NA, times = length(kineticPars[,1])), par_1 = NA, par_2 = NA) #par1/2 of a uniform are the lower bound and upper bound; par1/2 of a normal are the mean and variance
  kineticParPrior$distribution <- "unif"
  kineticParPrior$par_1 <- -10; kineticParPrior$par_2 <- 10
  
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
  
  colnames(markov_par_vals) <- ifelse(kineticPars$SpeciesType == "keq", "keq", kineticPars$commonName)
  
  colnames(markov_par_vals) <- unname(sapply(colnames(markov_par_vals), function(name_int){
    if(length(strsplit(name_int, split = "")[[1]]) >= 25){
      split_name <- strsplit(name_int, split = "")[[1]]
      split_pois <- c(1:length(split_name))[split_name %in% c(" ", "-")][which.min(abs(20 - c(1:length(split_name)))[split_name %in% c(" ", "-")])]
      split_name[split_pois] <- "\n"
      paste(split_name, collapse = "")
      }else{name_int}
    }))
  
  run_summary[[names(rxnList_form)[rxN]]]$markovChain <- markov_par_vals
  run_summary[[names(rxnList_form)[rxN]]]$likelihood <- lik_track
  
  }
current_pars <- markov_par_vals[which.max(lik_track),]

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

all_rxns <- names(rxnList_form)
rxn_fits <- NULL

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
  }
  
  
  pdf(file = paste(c("FBGA_files/outputPlots/", arxn, ".pdf"), collapse = ""), width = 20, height = 20)
  
  print(ggplot() + geom_violin(data = par_likelihood, aes(x = factor(index), y = likelihood), fill = "RED"))
  print(param_compare(run_rxn, par_markov_chain))
  
  dev.off()
  
  }

under_determined_rxnFits = rxn_fits[rxn_fits$residDF > 10,]
rownames(under_determined_rxnFits) <- under_determined_rxnFits$rxn
under_determined_rxnFits <- under_determined_rxnFits[,-c(1:2)]
under_determined_rxnFits_frac <- under_determined_rxnFits/under_determined_rxnFits$TSS
under_determined_rxnFits_frac <- under_determined_rxnFits_frac[,!(colnames(under_determined_rxnFits_frac) %in% c("TSS", "LS_met", "LS_enzyme"))]

under_determined_rxnFits_frac <- under_determined_rxnFits_frac[order(apply(under_determined_rxnFits_frac, 1, max, na.rm = TRUE)),]
under_determined_rxnFits_frac$reaction <- factor(rownames(under_determined_rxnFits_frac), levels = rownames(under_determined_rxnFits_frac))
rxnFit_frac_melt <- melt(under_determined_rxnFits_frac, id.vars = "reaction")

rxnFit_frac_plot <- ggplot(data = rxnFit_frac_melt, aes(x = reaction, y = value, fill = variable))
rxnFit_frac_plot + geom_bar(position = "dodge")



flux_fitting <- function(run_rxn, par_markov_chain, par_likelihood){
  # predict flux based upon parameter sets to determine how much variance in flux can be accounted for using the prediction
  param_interval <- exp(apply(par_markov_chain, 2, function(x){quantile(x, probs = c(0.025, 0.975))}))
  param_interval <- data.frame(cbind(t(param_interval), median = exp(apply(par_markov_chain, 2, median))))
  
  par_stack <- rep(1, n_c) %*% t(exp(par_markov_chain[which.max(par_likelihood$likelihood),])); colnames(par_stack) <- run_rxn$kineticPars$formulaName
  occupancy_vals <- data.frame(run_rxn$metabolites, par_stack)
  predOcc <- model.matrix(run_rxn$occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
 
  enzyme_activity <- (predOcc %*% t(rep(1, sum(run_rxn$all_species$SpeciesType == "Enzyme"))))*run_rxn$enzymes #occupany of enzymes * relative abundance of enzymes
  flux_fit <- nnls(enzyme_activity, run_rxn$flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  
  fit_summary <- data.frame(residDF = sum(run_rxn$flux != 0) - length(par_stack[1,]), parametricFit = NA, NNLS = NA, LS = NA, LS_met = NA, LS_enzyme = NA, TSS = NA)
  
  ### using flux fitted from the median parameter set, how much variance is explained
  fit_summary$parametricFit = anova(lm(run_rxn$flux ~ flux_fit$fitted))$S[1]
  
  run_rxn$rxnSummary$metsID2tID
  run_rxn$rxnSummary$rxnStoi
  
  
  ### using flux fitted using non-negative least squares regression, how much variance is explained ### metabolite abundances are corrected for whether the metabolite is a product (*-1) or reactant (*1)
  NNLSmetab <- run_rxn$metabolites * -1*(rep(1, n_c) %*% t(sapply(colnames(run_rxn$metabolites), function(x){run_rxn$rxnSummary$rxnStoi[names(run_rxn$rxnSummary$rxnStoi) == names(run_rxn$rxnSummary$metsID2tID)[run_rxn$rxnSummary$metsID2tID == x]]})))
    
  if(all(run_rxn$flux < 0)){
    NNLSanova <- anova(lm(-1*run_rxn$flux ~ nnls(as.matrix(data.frame(NNLSmetab, run_rxn$enzymes)), -1*run_rxn$flux)$fitted))
    fit_summary$NNLS <- ifelse(length(NNLSanova$S) != 1, NNLSanova$S[1], NA)
    }else{
      NNLSanova <- anova(lm(run_rxn$flux ~ nnls(as.matrix(data.frame(NNLSmetab, run_rxn$enzymes)), run_rxn$flux)$fitted))
      fit_summary$NNLS = ifelse(length(NNLSanova$S) != 1, NNLSanova$S[1], NA)
    }
  
  
  ### using LS regression, how much variance is explained 
  fit_summary$LS_met = anova(lm(run_rxn$flux ~ run_rxn$metabolites))$S[1]
  fit_summary$LS_enzyme = anova(lm(run_rxn$flux ~ run_rxn$enzymes))$S[1]
  fit_summary$LS = sum(anova(lm(run_rxn$flux ~ run_rxn$metabolites + run_rxn$enzymes))$S[1:2])
  fit_summary$TSS = sum(anova(lm(run_rxn$flux ~ run_rxn$metabolites + run_rxn$enzymes))$S)
  
  output <- list()
    output$fit_summary <- fit_summary
    output$param_interval <- param_interval
    output$fitted_flux <- flux_fit$fitted
  output
  }


param_compare <- function(run_rxn, par_markov_chain, par_likelihood){
  # visualize the joint and marginal distribution of parameter values from the markov chain
  
  par_combinations <- expand.grid(1:length(run_rxn$kineticPars[,1]), 1:length(run_rxn$kineticPars[,1]))
  like_comparison <- ifelse(par_combinations[,1] == par_combinations[,2], TRUE, FALSE)
  
  max_likelihood <- par_markov_chain[which.max(par_likelihood$likelihood),]
  
  par_comp_like <- NULL
  for(i in 1:sum(like_comparison)){
    par_comp_like <- rbind(par_comp_like, data.frame(xval = par_markov_chain[,par_combinations[like_comparison,][i,1]], parameter_1 = colnames(par_markov_chain)[par_combinations[like_comparison,][i,1]],
         parameter_2 = colnames(par_markov_chain)[par_combinations[like_comparison,][i,1]]))
      }
  
  par_comp_dissimilar <- NULL
  for(i in 1:sum(!like_comparison)){
    par_comp_dissimilar <- rbind(par_comp_dissimilar, data.frame(xval = par_markov_chain[,par_combinations[!like_comparison,][i,1]], yval = par_markov_chain[,par_combinations[!like_comparison,][i,2]], 
          parameter_1 = colnames(par_markov_chain)[par_combinations[!like_comparison,][i,1]], parameter_2 = colnames(par_markov_chain)[par_combinations[!like_comparison,][i,2]]))
      }
  
  MLEbarplot <- data.frame(xval = max_likelihood[par_combinations[like_comparison,1]], parameter_1 = colnames(par_markov_chain)[par_combinations[like_comparison,1]],
      parameter_2 = colnames(par_markov_chain)[par_combinations[like_comparison,1]])
  MLEpoints <- data.frame(xval = max_likelihood[par_combinations[!like_comparison,1]], yval = max_likelihood[par_combinations[!like_comparison,2]],
      parameter_1 = colnames(par_markov_chain)[par_combinations[!like_comparison,1]], parameter_2 = colnames(par_markov_chain)[par_combinations[!like_comparison,2]])
  
  
  
  #### determine the maximum bin from the histogram so that values can be scaled to the bivariate histogram values ###

  par_hist_binwidth = 0.2
  
  max_density <- max(apply(par_markov_chain, 2, function(x){max(table(round(x/par_hist_binwidth)))}))
  
  par_comp_dissimilar$yval <- par_comp_dissimilar$yval*(max_density/20) + max_density/2
  density_trans_inv <- function(x){x*(max_density/20) + max_density/2}
  density_trans <- function(x){(x - max_density/2)/(max_density/20)}
  
  
  
  hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), 
      legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line"), axis.title = element_blank()) 

  #labels = c(0:floor(max_density/50))*50, breaks = c(0:floor(max_density/50))*50
  
  ggplot() + geom_hex(data = par_comp_dissimilar, aes(x = xval, y = yval)) + geom_bar(data = par_comp_like, aes(x = xval), binwidth = par_hist_binwidth, col = "black") + facet_grid(parameter_2 ~ parameter_1, scales = "fixed") + hex_theme +
    scale_fill_gradientn(name = "Counts", colours = c("white", "darkgoldenrod1", "chocolate1", "firebrick1", "black")) +
    scale_x_continuous(NULL, expand = c(0.02,0.02)) + scale_y_continuous(NULL, expand = c(0.01,0.01), labels = density_trans, breaks = density_trans_inv(seq(-10, 10, by = 5))) +
    geom_vline(data = MLEbarplot, aes(xintercept = xval), col = "cornflowerblue", size = 2) + geom_point(data = MLEpoints, aes(x = xval, y = yval), size = 2, col = "cornflowerblue")

  }











