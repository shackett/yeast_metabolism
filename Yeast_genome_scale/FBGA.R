#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

library(reshape2) #for visualization at the end
library(nnls) #for non-negative regression used to fit kinetic parameters
library(ggplot2)

options(stringsAsFactors = FALSE)

n_c <- 25

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
    markov_pars$n_samples <- 100 #how many total markov samples are desired
    markov_pars$burn_in <- 0 #how many initial samples should be skipped


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

n_c <- length(rxnList_form[[1]]$flux)

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









flux_fitting <- function(x){
  # predict flux based upon parameter sets to determine how much variance in flux can be accounted for using the prediction
  param_interval <- exp(apply(markov_par_vals, 2, function(x){quantile(x, probs = c(0.025, 0.975))}))
  param_interval <- data.frame(cbind(t(param_interval), median = exp(apply(markov_par_vals, 2, median))))
  
  par_stack <- rep(1, n_c) %*% t(param_interval$median); colnames(par_stack) <- kineticPars$formulaName
  occupancy_vals <- data.frame(met_abund, par_stack)
  predOcc <- model.matrix(occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*enzyme_abund #occupany of enzymes * relative abundance of enzymes
  flux_fit <- nnls(enzyme_activity, flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  
  fit_summary <- data.frame(parametricFit = NA, NNLS = NA, LS = NA, LS_met = NA, LS_enzyme = NA, TSS = NA)
  
  ### using flux fitted from the median parameter set, how much variance is explained
  fit_summary$parametricFit = anova(lm(flux ~ flux_fit$fitted))$S[1]
  
  ### using flux fitted using non-negative least squares regression, how much variance is explained
  fit_summary$NNLS = anova(lm(flux ~ nnls(as.matrix(data.frame(met_abund, enzyme_abund)), flux)$fitted))$S[1]
  
  ### using LS regression, how much variance is explained 
  fit_summary$LS_met = anova(lm(flux ~ met_abund))$S[1]
  fit_summary$LS_enzyme = anova(lm(flux ~ enzyme_abund))$S[1]
  fit_summary$LS = sum(anova(lm(flux ~ met_abund + enzyme_abund))$S[1:2])
  fit_summary$TSS = sum(anova(lm(flux ~ met_abund + enzyme_abund))$S)
  }


param_compare <- function(kineticPars){
  # visualize the joint and marginal distribution of parameter values from the markov chain
  
  par_combinations <- expand.grid(1:length(kineticPars[,1]), 1:length(kineticPars[,1]))
  like_comparison <- ifelse(par_combinations[,1] == par_combinations[,2], TRUE, FALSE)
  
  par_comp_like <- NULL
  for(i in 1:sum(like_comparison)){
    par_comp_like <- rbind(par_comp_like, data.frame(xval = markov_par_vals[,par_combinations[like_comparison,][i,1]], parameter_1 = colnames(markov_par_vals)[par_combinations[like_comparison,][i,1]],
         parameter_2 = colnames(markov_par_vals)[par_combinations[like_comparison,][i,1]]))
      }
  
  par_comp_dissimilar <- NULL
  for(i in 1:sum(!like_comparison)){
    par_comp_dissimilar <- rbind(par_comp_dissimilar, data.frame(xval = markov_par_vals[,par_combinations[!like_comparison,][i,1]], yval = markov_par_vals[,par_combinations[!like_comparison,][i,2]], 
          parameter_1 = colnames(markov_par_vals)[par_combinations[!like_comparison,][i,1]], parameter_2 = colnames(markov_par_vals)[par_combinations[!like_comparison,][i,2]]))
      }
  
  #### determine the maximum bin from the histogram so that values can be scaled to the bivariate histogram values ###

  par_hist_binwidth = 0.2
  
  max_density <- max(apply(markov_par_vals, 2, function(x){max(table(round(x/par_hist_binwidth)))}))
  
  
  par_comp_dissimilar$yval <- par_comp_dissimilar$yval*(max_density/20) + max_density/2
  density_trans_inv <- function(x){x*(max_density/20) + max_density/2}
  density_trans <- function(x){(x - max_density/2)/(max_density/20)}
  
  density_trans <- function(){function(x) (x - max_density/2)/(max_density/20)}
  
  
  hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), 
      legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line"), axis.title = element_blank()) 

  ggplot() + geom_hex(data = par_comp_dissimilar, aes(x = xval, y = yval)) + geom_bar(data = par_comp_like, aes(x = xval), binwidth = par_hist_binwidth, col = "black") + facet_grid(parameter_2 ~ parameter_1, scales = "fixed") + hex_theme +
    scale_fill_gradientn(name = "Counts", colours = c("white", "darkgoldenrod1", "chocolate1", "firebrick1", "black"), labels = c(0:floor(max_density/50))*50, breaks = c(0:floor(max_density/50))*50) +
    scale_x_continuous(NULL, expand = c(0.02,0.02)) + scale_y_continuous(NULL, expand = c(0.01,0.01), labels = density_trans, breaks = density_trans_inv(seq(-10, 10, by = 5)))

  
  ggsave("mcmc_output.pdf", width = 20, height = 20)
  
  
     
  }






########## TOY ##############

# species involved in reaction
n_c <- 100
metabs <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(metabs) <- c("fu", "bar")
enzyme <- matrix(rgamma(50, 1, 2), nrow = n_c, ncol = 2); colnames(enzyme) <- c("yfg1", "yfg2")
# Model parameters, and correspondence between parameters and species - i.e. kcat - species 1 matches the activity of the first enzyme catalyzing the reaction. km - species 1 means this the affinity of the first metabolite (column 1) for  the enzyme
assocConst <- data.frame(name = c("yfg1_kcat", "fu_km", "C"), type = c("kcat", "km", "ki"), specie = c(1, 1, 2), priorMean = NA, priorSD = NA) # all of these parameters are lognormally distributed, the parameterization that I am using for this is the parameter MLE (prior mean) and the log standard deviation (priorSD)
assocConst$priorSD <- assocConst$priorMean <- c(1, 0.25, 0.8)
assocConstTRUE <- assocConst; assocConstTRUE$priorMean <- c(1, 0.5, 0.6)

trueFlux <- (enzyme[,1] * assocConstTRUE$priorMean[1] * metabs[,1] / (metabs[,1] + assocConstTRUE$priorMean[2])) #simulate measured fluxes from the rxn equation form assuming that metabolites and enzymes are measured accurately
measuredFlux <- trueFlux*rlnorm(n_c, 0, 0.25) #simulate measured fluxes from the rxn equation form with added lognormal noise

parNum <- length(assocConst[,1]) #how many parameters are there in the model
plot(measuredFlux ~ metabs[,1])
plot(measuredFlux ~ enzyme[,1])


# function predicting flux from provided paramters
#exampleformula
#enzyme activity
#substrate occupancy

rxnEqn <- as.formula(" ~ I(yfg1_kcat * yfg1 * fu / (fu + fu_km)) + 0")



 ## Genetic algorithm parameters ##
N = 5000 # the number of parameter sets competing for ultimate fit
mu = 0.1*length(assocConst[,1]) # the lambda mutation rate across all parameters 
generations <- 100
genMeanlogFit <- rep(NA, generations)

#initialize with parameters drawn from the parameters prior
gaConstants <- sapply(1:parNum, function(initConstN){
  #rlnorm(N, assocConst[initConstN,colnames(assocConst) == "priorMean"], assocConst[initConstN,colnames(assocConst) == "priorSD"])
  exp(rnorm(N, log(assocConst[initConstN,colnames(assocConst) == "priorMean"]), assocConst[initConstN,colnames(assocConst) == "priorSD"]))
  })
colnames(gaConstants) <- assocConst$name

for(genN in 1:generations){
  ### generations of selection (implicitely includes drift - because sampling is proportional to fitness but N is finite), and mutation
  
  ### Mutation ###
  
  mutations <- rpois(N, mu)
  mutations <- sapply(mutations, function(x){min(c(x, length(assocConst[,1])))}) #floor to the maximum number of parameters
  mutInds <- c(1:N)[mutations != 0]
  
  newPar <- sapply(1:length(mutInds), function(mutInd){
    replpar <- sample(1:parNum, mutations[mutInd], replace = TRUE)
    indPars <- gaConstants[mutInd,]
    for(mut in replpar){
      indPars[mut] <- exp(rnorm(1, log(assocConst[mut,colnames(assocConst) == "priorMean"]), assocConst[mut,colnames(assocConst) == "priorSD"]))
        #rlnorm(1, assocConst[mut,colnames(assocConst) == "priorMean"], assocConst[mut,colnames(assocConst) == "priorSD"])
      }
    indPars
    })
  
  gaConstants[mutInds,] <- t(newPar)
  
  ######
  
  indPrior <- apply(sapply(1:parNum, function(parN){
    dnorm(log(gaConstants[,parN]), log(assocConst$priorMean[parN]),  assocConst$priorSD[parN])
    #dlnorm(gaConstants[,parN], assocConst$priorMean[parN], assocConst$priorSD[parN])
  }), 1, prod) #presolve the prior probability of a model (of parameter values)

  ### Selection ###
  
  indLogFit <- sapply(1:N, indFitnessFxn)
  
  rFit <- exp(indLogFit - max(indLogFit)) #fitness relative to most fit individual
  progeny <- rowSums(rmultinom(N, 1, rFit))
  progeny <- unlist(sapply(1:N, function(x){rep(x, progeny[x])}))
  
  genMeanlogFit[genN] <- mean(indLogFit[progeny]) #save the average log fitness of individuals

  gaConstants <- gaConstants[progeny,] #replace current population with selected progeny
  
  if(genN/10 == floor(genN/10)){
    print(paste("Generation", genN, "complete.  The mean log fitness is", signif(genMeanlogFit[genN], 5), sep = " "))
    }
  
  }

plot(genMeanlogFit, pch = 16, col = "RED", cex = 0.3)

finalFit <- t(sapply(1:N, indFitnessFxn, splitFit = TRUE))
colnames(finalFit) <- c("fit", "prior")

# visualize the parameter correlation matrix
library(colorRamps)
library(gplots)
parCorrMat <- cor(log(gaConstants))
#parSpearCorrMat <- cor(log(gaConstants), method = "spearman")
heatmap.2(parCorrMat, symm = TRUE, col = blue2yellow(1000), dendrogram = "row", symkey = TRUE)

# visualize parameter distributions colored by score

par_fit <- data.frame(gaConstants, finalFit)
par_fit_stack <- melt(par_fit, measure.vars = assocConst$name)
par_fit_stack$logVal = log(par_fit_stack$value)
#par_fit_stack2 <- melt(par_fit_stack, measure.vars = c("fit", "prior"))
#colnames(par_fit_stack2) <- c("Parameter", "ParValue", "BayesFactorComponent", "BFVal"); par_fit_stack2$BFVal <- exp(par_fit_stack2$BFVal)
#dist_plotter <- ggplot(par_fit_stack, aes(x =  logVal)) + facet_grid(variable ~ .) + geom_abline
#dist_plotter + geom_histogram(par_fit_stack, aes(x = logVal), binwidth = 0.01)

colnames(par_fit_stack) <- c("Fit", "Prior", "name", "Value", "lnValue")
maxCatVal <- sapply(assocConst$name, function(x){max(table(floor(par_fit_stack$lnValue[par_fit_stack$name == x] * 100)))})
assocConst$CImin <- log(assocConst$priorMean) - 2 * assocConst$priorSD; assocConst$CImax <- log(assocConst$priorMean) + 2 * assocConst$priorSD
assocConst$maxBinCount = maxCatVal


hist_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "none", 
  panel.grid.minor = element_blank(), legend.key.width = unit(6, "line"), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan")) 

dist_plotter <- ggplot() + facet_grid(name ~ ., scales = "free_y") + hist_theme
dist_plotter + geom_histogram(data = par_fit_stack, aes(x = lnValue), binwidth = 0.1) + geom_vline(data = assocConst, aes(xintercept = log(priorMean), colour = "RED", size = 2)) +
  geom_errorbarh(data = assocConst, aes(x = log(priorMean), xmin = CImin, xmax = CImax, y = maxCatVal/2, height = maxBinCount/10, size = 2, colour = "RED"))

ggsave("gatoy.pdf", height = 10, width = 10)







#dist_plotter <- ggplot(par_fit_stack2, aes(x =  ParValue, fill = BFVal)) + facet_grid(BayesFactorComponent ~ Parameter)
#dist_plotter + geom_histogram() + scale_fill_gradient(name = "woot", low = "black", high = "firebrick1")

#scale_fill_gradient(name = "Counts", low = "black", high = "firebrick1", trans = "log", breaks = hex_breaks, labels = hex_breaks)

# for pairs of parameters exceeding a correlation cutoff plot a bivariate histogram
# find a good example first
corCO <- .2

hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "gray90"), legend.position = "top", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line")) 

parPoint <- sapply(rownames(parCorrMat), function(x){paste(x, colnames(parCorrMat), sep = "/")})
parAssoc <- parPoint[upper.tri(parCorrMat)][abs(parCorrMat[upper.tri(parCorrMat)]) > corCO]
if(length(parAssoc) != 0){
  
  for(parSet in parAssoc){
    bivariateDF <- data.frame(log(gaConstants[,c(1:parNum)[assocConst$name %in% strsplit(parSet, split = '/')[[1]]]]))
    colnames(bivariateDF) <- c("P1", "P2")
    print(ggplot(bivariateDF, aes(x = P1, y = P2)) + geom_hex() + scale_x_continuous(assocConst$name[assocConst$name %in% strsplit(parSet, split = '/')[[1]]][1], expand = c(0.02,0.02)) +
      scale_y_continuous(assocConst$name[assocConst$name %in% strsplit(parSet, split = '/')[[1]]][2], expand = c(0.02,0.02)) + hex_theme +
      scale_fill_gradient(name = "Counts", low = "black", high = "firebrick1"))
    }
   
     
  }






indFitnessFxn <- function(i, splitFit = FALSE){
  indParams <- t(t(rep(1, n_c))) %*% gaConstants[i,]; colnames(indParams) <- assocConst$name
  indDat <- data.frame(metabs, enzyme, indParams)
  
  predFlux <- model.matrix(rxnEqn, data = indDat)[,1] #evaluates the combination of data and parameters in indDat according to the reaction equation
  
  # setting the mean of measured flux to the mean of the predicted flux is equivalent to scaling flux magnitude to account for kcat (when weights are not used this is true, otherwise the weighted mean needs to be determined)
  # Determine how to propogate weights associated with measurement through to get a weight on predicted flux
  
  lpredFluxcent <- log(predFlux) + (mean(log(measuredFlux)) - mean(log(predFlux))) #centered log flux
  
  SSE <- sum((lpredFluxcent - log(measuredFlux))^2) #determine the sum of squared error between predicted log flux and observed log flux
  var_fit <- SSE / (n_c - 1) #correct estimate of variance according to number of fitted parameters + 1
  
  if(splitFit == FALSE){
  sum(dnorm(lpredFluxcent, log(measuredFlux), sqrt(var_fit), log = TRUE)) + log(indPrior[i]) # bayes factor for this parameter set is returned
    }else{
      c(sum(dnorm(lpredFluxcent, log(measuredFlux), sqrt(var_fit), log = TRUE)), log(indPrior[i]))
      }
  }








