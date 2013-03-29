#### Determining the extent to which flux predicted based on metabolomics and proteomics can be 
## aligned to experimentally measured flux ###

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

options(stringsAsFactors = FALSE)

n_c <- 25

#import list of metabolite abundances and rxn forms

load("rxnf_formulametab.rdata")

#import fluxes
load("fluxSummaryQP.Rdata") #load flux through each reaction

#import enzyme abundances
enzyme_abund <- read.delim("../ChemicalSpeciesQuant/Proteomics/run_Deg/relAbundMatrix.tsv")
#### load protLMfile.Rdata when rerun in order to get all of the genes that a non-unique peptide-set could correspond to
#measured_genes <- unlist(sapply(colnames(prot_abund_final), function(x){strsplit(x, "/")}))

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



rxnFile <- read.delim('rxn_yeast.tsv', stringsAsFactors = FALSE)
met_genes <- data.frame(reaction = unique(rxnFile$ReactionID), genes = NA, pathway = NA)

for(rxN in 1:length(met_genes[,1])){
  rxSubset <- rxnFile[rxnFile$ReactionID == met_genes$reaction[rxN],]
  met_genes$genes[rxN] <- paste(rxSubset$MetName[is.na(rxSubset$StoiCoef)], collapse = "/")
  gene_subset <- strsplit(paste(rxSubset$MetName[is.na(rxSubset$StoiCoef)], collapse = ":"), split = ":")[[1]]
  met_genes$pathway[rxN] <- paste(unique(strsplit(paste(kegg_enzyme_dict[kegg_enzyme_dict$SYST %in% gene_subset,]$PATHWAY, collapse = "__"), "__")[[1]]), collapse = "__")
  }

###### all reactions #####

rxnList_all <- rxnf

for(rxN in c(1:length(rxnList_all))){
  kegg_subset <- met_genes[met_genes$reaction == names(rxnList_all)[rxN],]
  if(length(kegg_subset[,1]) == 0){next}
  
  rxnList_all[[names(rxnList_all)[rxN]]]$pathway = kegg_subset$pathway
  rxnList_all[[names(rxnList_all)[rxN]]]$genes = kegg_subset$genes
  rxnList_all[[names(rxnList_all)[rxN]]]$enzymeAbund = enzyme_abund[enzyme_abund$Gene %in% strsplit(kegg_subset$genes, split = '[/:]')[[1]],]
  }



####### rxns which carry flux #####

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


##### predict flux ####

reaction_pred_log <- data.frame(nmetab = rep(NA, length(rxnList)), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedMetab = NA, varExplainedEnzy = NA, TSS = NA)
reaction_pred_linear <- data.frame(nmetab = rep(NA, length(rxnList)), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedMetab = NA, varExplainedEnzy = NA, TSS = NA)

cond_mapping <- data.frame(flux_cond = names(rxnList[[1]]$flux), met_reordering = c(11:15, 1:5, 6:10, 16:20, 21:25), enzyme_reordering = c(17:21, 2:6, 12:16, 7:11, 23:27)) #remove this when metabolite conditions are remapped onto flux/protein ones prior to this script
cond_mapping$met_cond = rownames(rxnList[[1]]$rxnMet)[cond_mapping$met_reordering]
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
    } else if(reaction_pred_linear$nenz[rxN] != 0 & reaction_pred_linear$nmetab[rxN] == 0){
    ### only enzymes ###
      ### prediction using log measures ###
      reaction_pred_log$Fenz[rxN] <- anova(lm(rxFlux ~ lenzymes))$F[1]
      reaction_pred_log$varExplainedEnzy[rxN] <- anova(lm(rxFlux ~ lenzymes))$Sum[1]
      reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lenzymes))$Sum)
      
      ### prediction using linear measures ###
      reaction_pred_linear$Fenz[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes))$F[1]
      reaction_pred_linear$varExplainedEnzy[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes))$Sum[1]
      reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux ~ enzymes))$Sum)  
      
    } else if(reaction_pred_linear$nenz[rxN] == 0 & reaction_pred_linear$nmetab[rxN] != 0){
    ### only metabolites ###
      ### prediction using log measures ###
      reaction_pred_log$Fmetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$F[1]
      reaction_pred_log$varExplainedMetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$Sum[1]
      reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lmetabs))$Sum)
      
      ### prediction using linear measures ###
      reaction_pred_linear$Fmetab[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ metabs))$F[1]
      reaction_pred_linear$varExplainedMetab[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ metabs))$Sum[1]
      reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux ~ metabs))$Sum)
      
    } else{
    ### both metabolites and enzymes ###
      reaction_pred_linear$Fmetab[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes + metabs))$F[2]
      reaction_pred_linear$Fenz[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes + metabs))$F[1]
      reaction_pred_linear$varExplainedMetab[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes + metabs))$Sum[2]
      reaction_pred_linear$varExplainedEnzy[rxN] <- anova(lm(rxnList[[rxN]]$flux ~ enzymes + metabs))$Sum[1]
      reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux ~ enzymes + metabs))$Sum)
      
      reaction_pred_log$Fmetab[rxN] <- anova(lm(rxFlux ~ lenzymes + lmetabs))$F[2]
      reaction_pred_log$Fenz[rxN] <- anova(lm(rxFlux ~ lenzymes + lmetabs))$F[1]
      reaction_pred_log$varExplainedMetab[rxN] <- anova(lm(rxFlux ~ lenzymes + lmetabs))$Sum[2]
      reaction_pred_log$varExplainedEnzy[rxN] <- anova(lm(rxFlux ~ lenzymes + lmetabs))$Sum[1]
      reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lenzymes + lmetabs))$Sum)
      
    }
  }

qplot(reaction_pred_linear$Fmetab)
qplot(reaction_pred_linear$Fenz)

reaction_pred_summary_log <- data.frame(N = reaction_pred_log$nCond, metaboliteVarianceExplained = reaction_pred_log$varExplainedMetab/reaction_pred_log$TSS, enzymeVarianceExplained = reaction_pred_log$varExplainedEnzy/reaction_pred_log$TSS)
reaction_pred_summary_linear <- data.frame(N = reaction_pred_linear$nCond, metaboliteVarianceExplained = reaction_pred_linear$varExplainedMetab/reaction_pred_linear$TSS, enzymeVarianceExplained = reaction_pred_linear$varExplainedEnzy/reaction_pred_linear$TSS)

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
  scale_fill_manual(values = c("enzymeVarianceExplained" = "sienna1", "metaboliteVarianceExplained" = "steelblue1")) + scale_color_identity()
  

ggsave("varianceExplained.pdf", width = 20, height = 12)
        
library(reshape2) #for visualization at the end


# form a Km hyperparameter using the distribution of [met]/km estimates




# look at subset of reactions which carry flux
# using 






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
dist_plotter + geom_histogram(data = par_fit_stack, aes(x = lnValue), binwidth = 0.01) + geom_vline(data = assocConst, aes(xintercept = log(priorMean), colour = "RED", size = 2)) +
  geom_errorbarh(data = assocConst, aes(x = log(priorMean), xmin = CImin, xmax = CImax, y = maxCatVal/2, height = maxCatVal/10, size = 2, colour = "RED"))







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








