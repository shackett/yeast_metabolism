#### Libraries ####
library(gplots)
library(ggplot2)
library(data.table)
library(reshape2)
library(RColorBrewer)

#### Options ####

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")
source("FBA_lib.R")

options(stringsAsFactors = FALSE)


# Specify whether growth optimization should be through linear programming (LP) - maximizing biomass given nutrients or 
# through quadratic programming (QP) - optimally matching experimental boundary fluxes to optimized ones.
QPorLP <- "QP" # LP and PhPP no longer supported
if(QPorLP == "QP"){
  #ln -s /Library/gurobi510/mac64/lib/libgurobi51.so libgurobi51.so
  library(gurobi) # QP solver interface
  }else{print("try the QP, all the cool kids are doing it")}


###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing valid reactions, species (their composition) both from the core SBML model and supplemented manual annotations #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

### Load SBML files describing metabolites, rxn stoichiometry ###

inputFilebase = "yeast"
rxnFile = read.delim(paste("rxn_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
rxnparFile = read.delim(paste("rxn_par_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
corrFile = read.delim(paste("spec_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
specparFile = read.delim(paste("species_par_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
compFile <- read.delim(paste("comp_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
fluxDirFile <- read.delim(paste("flux_dir_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)

### Add additional reactions ###

customRx <- "customRxns.txt"
customList <- parse_custom("customRxns.txt")

rxnFile <- rbind(rxnFile, customList$rxnFile)
rxnparFile <- rbind(rxnparFile, customList$rxnparFile)
corrFile <- rbind(corrFile, customList$corrFile)
specparFile <- rbind(specparFile, customList$specparFile)
fluxDirFile <- rbind(fluxDirFile, customList$fluxDirFile)

### Determine unique metabolites and reactions ###

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)

#### Load or write stoichiometry matrix of reactions and their altered metabolites ####

if(!file.exists("flux_cache/yeast_stoi.Rdata")){write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}
load("flux_cache/yeast_stoi.Rdata")

named_stoi <- stoiMat # create a version of the stoichiometric matrix with common names for rxns and mets to allow for easier searching
met_dict <- metIDtoSpec(rownames(named_stoi)); met_dict <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
rxn_dict <- rxnIDtoEnz(colnames(named_stoi)); rxn_dict <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})
rownames(named_stoi) <- met_dict; colnames(named_stoi) <- rxn_dict

### Elemental composition of metabolites ###

if(!file.exists("flux_cache/stoiMetsComp.tsv")){elemental_composition(metabolites)}
modelMetComp <- read.delim("flux_cache/stoiMetsComp.tsv", header = TRUE)





#rxnFileAppend <- read.delim("customRxns.txt", sep = "\t")
#rxnFile <- rbind(rxnFile, rxnFileAppend)

#customMets <- read.delim("customMets.txt", sep = "\t")
#corrFileAppend <- customMets[,colnames(customMets) %in% c("SpeciesID", "SpeciesName", "SpeciesType", "Compartment")]
#corrFile <- rbind(corrFile, corrFileAppend)

#rxnparFileAppend <- customMets[,colnames(customMets) %in% c("SpeciesID", "SpeciesType", "Evidence")]
#colnames(rxnparFile) <- colnames(rxnparFileAppend)
#rxnparFile <- rbind(rxnparFile, rxnparFileAppend)

#MetCompAppend <- modelMetComp[1:nrow(customMets),]
#MetCompAppend$ID <- customMets$SpeciesID
#MetCompAppend$name <- customMets$SpeciesName
#MetCompAppend[,!(colnames(MetCompAppend) %in% c("ID", "name"))] <- NA

#modelMetComp <- rbind(modelMetComp, MetCompAppend)

###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing boundary conditions and reaction reversibility #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

#compositionFile <- read.delim("../Yeast_comp_energy.txt") #energy required to assimilate biomass components
nutrientFile <- read.delim("Boer_nutrients.txt")[1:6,1:6]; nutrientCode <- data.frame(nutrient = colnames(nutrientFile)[-1], shorthand = c("n", "p", "c", "L", "u"))
rownames(nutrientFile) <- nutrientFile[,1]; nutrientFile <- nutrientFile[,-1]
nutModelNames <- data.frame(commonName = rownames(nutrientFile), modelName = sapply(rownames(nutrientFile), function(x){paste(x, '[extracellular]')}))
rownames(nutrientFile) <- nutModelNames$modelName

load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)


### add model annotations of flux direction and reaction reversibility as well as Vito's implementation of Elad's component contribution free energy prediction###

reversibleRx <- data.frame(rx = reactions, reversible = 0, CCdG = NA, CCdGsd = NA, CCdGdir = NA, modelRev = NA, modelBound = NA, manual = NA, rxFlip = NA, annotComment = NA)

reversibleRx$modelRev = fluxDirFile$Reversible[chmatch(reversibleRx$rx, fluxDirFile$ReactionID)]
reversibleRx$modelBound = fluxDirFile$FluxBound[chmatch(reversibleRx$rx, fluxDirFile$ReactionID)]

#table(fluxDirFile$Reversible, fluxDirFile$FluxBound)

reversibleRx$reversible = ifelse(reversibleRx$modelBound == "greaterEqual", 1, 0) # directionality is coded as -1: irreversibly backward, 0: reversible, 1: irreversibly forward



#ccPred <- read.delim("cc_dG_matlab.tsv")
#ccPred$dir <- 0
#ccPred$dir[ccPred$dGr - 1.96*ccPred$dGrSD > 30] <- -1
#ccPred$dir[ccPred$dGr + 1.96*ccPred$dGrSD < -30] <- 1

#reversibleRx[chmatch(ccPred$reaction, reversibleRx$rx), colnames(reversibleRx) %in% c("CCdG", "CCdGsd", "CCdGdir")] <- ccPred[,-1]

#reversibleRx$reversible[!is.na(reversibleRx$CCdGdir)] <- reversibleRx$CCdGdir[!is.na(reversibleRx$CCdGdir)]

#read in manually flipped and directed reactions

#thermAnnotate = read.delim("thermoAnnotate.txt", header = TRUE, sep = "\t")
#for(rxN in 1:nrow(thermAnnotate)){
  #flip reaction direction (and free energy) if stated directionality is unconventional  
#  if(!is.na(thermAnnotate$flip[rxN]) & thermAnnotate$flip[rxN]){
#    stoiMat[,colnames(stoiMat) == thermAnnotate$reaction[rxN]] <- stoiMat[,colnames(stoiMat) == thermAnnotate$reaction[rxN]]*-1
#    reversibleRx$reversible[reversibleRx$rx == thermAnnotate$reaction[rxN]] <- reversibleRx$reversible[reversibleRx$rx == thermAnnotate$reaction[rxN]]*-1
#    }
  
  #manually define reaction direction
#  reversibleRx$rxFlip[reversibleRx$rx == thermAnnotate$reaction[rxN]] <- thermAnnotate$flip[rxN]
#  reversibleRx$manual[reversibleRx$rx == thermAnnotate$reaction[rxN]] <- thermAnnotate$direction[rxN]
#  }
#reversibleRx$reversible[!is.na(reversibleRx$manual)] <- reversibleRx$manual[!is.na(reversibleRx$manual)]








#### Associate rxns with enzyme ascertainment in proteomics s.t. flux can favor measured pathways ####
### Associate proteins with reactions and format complexes

rxn_enzyme_groups = NULL
for(rxN in 1:nrow(rxnparFile)){
  
  enzymeCombos <- rxnparFile$Enzymes[rxN]
  
  if(enzymeCombos == ""){
    next
  } # no enzymes
  
  enzymeCombos <- strsplit(enzymeCombos, ' OR ')[[1]]
  
  for(combo in 1:length(enzymeCombos)){
    rxn_enzyme_groups <- rbind(rxn_enzyme_groups, data.frame(reaction = rxnparFile$ReactionID[rxN], group = combo, enzyme = regmatches(enzymeCombos[combo], gregexpr('[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}', enzymeCombos[combo]))[[1]]))
  }
}



enzyme_abund <- read.delim("../ChemicalSpeciesQuant/Proteomics/proteinAbundance.tsv")
rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]

prot_matches <- sapply(reactions, function(x){
  rxMatches <- rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction == x]
  length(rownames(enzyme_abund)[rownames(enzyme_abund) %in% rxMatches]) != 0
  })

### cache files and pass to reaction equation formulation script ###

write.table(rxn_enzyme_groups, file = "flux_cache/rxn_enzyme_groups.tsv", sep = "\t", col.names = T, row.names = F, quote = F) # a data.frame indicating how proteins form catalytic units (monimers, dimers...)
write.table(prot_matches, file = "flux_cache/prot_matches.tsv", sep = "\t", col.names = T, row.names = F, quote = F) # boolean vector indicating whether a reaction's proteins were ascertained via proteomics




### Determine which metabolic pathways are associated with a reaction ###

reactionMatches <- data.frame(reactionID = rxnparFile$ReactionID[grep('reaction/', rxnparFile$Annotation)], KEGGrxnID = sapply(grep('reaction/', rxnparFile$Annotation, value = T), function(x){regmatches(x, regexpr('R[0-9]+', x))}))

reactions_to_pathways = read.delim("http://rest.kegg.jp/link/reaction/pathway", header = FALSE); colnames(reactions_to_pathways) <- c("pathwayCode", "RID")
reactions_to_pathways$RID = sapply(reactions_to_pathways$RID, function(x){sub('rn:', '', x)})

rxPathways = read.delim("http://rest.kegg.jp/list/pathway", header = FALSE); colnames(rxPathways) = c("pathwayCode", "pathway")
reactions_to_pathways$pathway = rxPathways$pathway[chmatch(reactions_to_pathways$pathwayCode, rxPathways$pathwayCode)]

reactionMatches$pathway = sapply(reactionMatches$KEGGrxnID, function(x){
  PS = reactions_to_pathways$pathway[reactions_to_pathways$RID == x]
  PS <- PS[!is.na(PS)]
  paste(PS, collapse = "__")
})


### Determine which pathways are associated with a protein ###

kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv") # > generated from keggIDparser.py KEGG IDs relative to yeast gene name (1gene -> 1 KEGG id, multiple  mapping between KEGG and genes)
if(is.null(kegg_enzyme_dict$PATHWAY)){
  gene_pathways(kegg_enzyme_dict); kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv")
}#generate a per-gene pathway annotation if one is not already generated




# favor flux through chosen central carbon metabolism pathways by reducing the L1 penalization of |v|

pathways <- sapply(reactionMatches$pathway, function(x){strsplit(x, '__')[[1]]})
unq_pathways <- unique(unlist(pathways))
centralC <- c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)")

centralCmatch <- unlist(lapply(pathways, function(x){
  sum(x %in% centralC) != 0
  }))

centralCrxnMatch <- rep(FALSE, times = length(reactions))
centralCrxnMatch[reactions %in% reactionMatches$reactionID[centralCmatch]] <- TRUE

rxn_search('exchange', is_rxn = F)

centralCrxnMatch[names(centralCrxnMatch) %in% 'r_0449'] <- FALSE ### remove FBPase's annotation as Glyc/Gluconeogenic in order to direct excess ATP hydrolysis through a designated reaction (r_4042) 

# favor flux through cytosolic ATPase, + transport of ATP, ADP + Pi to restrict the wasting of excess energy to a single reaction

ATPbreakdownRxns <- c("r_4042")#, "r_0249", "r_1149", "r_1150", "r_1167", "r_1168", "r_1459", "r_1460", "r_1461", "r_1462", "r_1463", "r_1868")


prot_penalty <- (1 - prot_matches)/2 + (1 - centralCrxnMatch)/2 # penalization by fraction of non-measured enzymes and favor central C metabolism
prot_penalty[(reactions %in% ATPbreakdownRxns)] <- 0.1

####
#save(list = ls(), file = "FBAinputFiles.Rdata")




#### Define the treatment in terms of nutrient availability and auxotrophies ####

treatment_par <- list()
n_c <- 25
for(i in 1:n_c){
  #define nutrient uptake and excretion rate - soft matches on (using maximal available for now)
  #measured_bounds <- data.frame(nutrient = rownames(nutrientFile), conc_per_t = nutrientFile[,colnames(nutrientFile) == nutrientCode$nutrient[nutrientCode$shorthand == chemostatInfo$limitation[i]]]*chemostatInfo$actualDR[i])
  
  measured_bounds <- mediaSummary[mediaSummary$condition == chemostatInfo$condition[i],]
  measured_bounds <- rbind(measured_bounds, data.frame(condition = chemostatInfo$condition[i], specie = rownames(nutrientFile)[!(rownames(nutrientFile) %in% measured_bounds$specie)], change = NA, sd = NA, lb = 0, 
     ub = nutrientFile[,colnames(nutrientFile) == nutrientCode$nutrient[nutrientCode$shorthand == chemostatInfo$limitation[i]]][!(rownames(nutrientFile) %in% measured_bounds$specie)], type = "uptake", density = NA))
  
  # multiply steady-state concentrations by DR to get the uptake/excretion rates
  measured_bounds$change <- measured_bounds$change*chemostatInfo$actualDR[i]
  measured_bounds$sd <- measured_bounds$sd*chemostatInfo$actualDR[i]
  measured_bounds$lb <- measured_bounds$lb*chemostatInfo$actualDR[i]
  measured_bounds$ub <- measured_bounds$ub*chemostatInfo$actualDR[i]
  
  #remove phosphate because empirical uptake rates far exceed capacity of biomass assimilation
  measured_bounds <- data.frame(measured_bounds)
  measured_bounds[measured_bounds$specie == "phosphate", colnames(measured_bounds) %in% c("change", "sd")] <- NA
  measured_bounds <- data.table(measured_bounds)
  
  
  treatment_par[[chemostatInfo$condition[i]]][["nutrients"]] <- measured_bounds
  
  #define ura3 and leu2 auxotrophies
  if(chemostatInfo$limitation[i] == "L"){treatment_par[[chemostatInfo$condition[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID))}
  if(chemostatInfo$limitation[i] == "u"){treatment_par[[chemostatInfo$condition[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID))}
  if(chemostatInfo$limitation[i] %in% c("c", "p", "n")){treatment_par[[chemostatInfo$condition[i]]][["auxotrophies"]] <- NA}
  
  #define observed fluxes per culture volume #scaled to intracellular volume later
  biomass_match <- data.frame(specie = comp_by_cond$compositionFile$MetName, AltName = comp_by_cond$compositionFile$AltName,change = unname(-1*comp_by_cond$cultureMolarity[,colnames(comp_by_cond$cultureMolarity) == chemostatInfo$condition[i]]*chemostatInfo$actualDR[i]))
  biomass_list <- list()
  
  for(component in unique(comp_by_cond$compositionFile$varCategory)){
    principal_costs <- biomass_match[comp_by_cond$compositionFile$varCategory %in% component,]
    
    if(component == "Maintenance ATP hydrolysis"){
      total_costs <- principal_costs
      }else{
        #costs of monomer assimilation incorporated into biomass flux
        energetic_costs <- as.matrix(comp_by_cond$biomassExtensionE[comp_by_cond$biomassExtensionE$name %in% principal_costs$AltName,-1])
        energetic_costs_aggregate <- t(principal_costs$change) %*% energetic_costs; colnames(energetic_costs_aggregate) <- colnames(comp_by_cond$biomassExtensionE)[-1]
        total_costs <- rbind(principal_costs, data.frame(specie = colnames(energetic_costs_aggregate), AltName = colnames(energetic_costs_aggregate), change = t(unname(energetic_costs_aggregate)))[energetic_costs_aggregate != 0,])
      }
    biomass_list[[component]]$exchange = total_costs
    
    # define the accuracy of a constraint in terms of the coefficient of variation - sd over mean
    if(component %in% colnames(comp_by_cond$CV_table)){
      biomass_list[[component]]$SD = as.numeric(subset(comp_by_cond$CV_table, comp_by_cond$CV_table$condition == chemostatInfo$condition[i], component))
      }else{
        biomass_list[[component]]$SD = 1/10
      }
    }
  
  
  treatment_par[[chemostatInfo$condition[i]]][["boundaryFlux"]] = biomass_list
  }
possibleAuxotrophies = c(as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID)), as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID)))



#### Determine the compartmentation of each reaction ####

compartment <- sapply(reactions, function(x){rxnFile$Compartment[rxnFile$ReactionID == x][1]})



#### Define species involved in boundary-conditions ####

## extract the metabolite ID corresponding to the extracellular introduction of nutrients ##

sources <- c("D-glucose", "ammonium", '^phosphate \\[extracellular\\]', "sulphate", "uracil", "L-leucine")
		
resource_matches <- lapply(sources, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

boundary_met <- NULL
for(x in 1:length(sources)){
boundary_met <- rbind(boundary_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


## extract the IDs of excreted metabolites ##



excreted <- c('ethanol \\[extracellular\\]', 'acetate \\[extracellular\\]', 'glycerol \\[extracellular\\]')

resource_matches <- lapply(excreted, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

excreted_met <- NULL
for(x in 1:length(excreted)){
excreted_met <- rbind(excreted_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


## extract the metabolite ID corresponding to cytosolic metabolites being assimilated into biomass ##

sinks <- unique(c(comp_by_cond$compositionFile$AltName, colnames(comp_by_cond$biomassExtensionE)[-1]))


resource_matches <- lapply(sinks, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile, reduceByLength = T)

comp_met <- NULL
for(x in 1:length(sinks)){
  comp_met <- rbind(comp_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "cytoplasm"],])
}



## freely exchanging metabolites through extracellular compartment ##

free_flux <- c("carbon dioxide", "oxygen", "H2O", "H+")

resource_matches <- lapply(free_flux, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

freeExchange_met <- NULL
for(x in 1:length(free_flux)){
  freeExchange_met <- rbind(freeExchange_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}




#### Search for rxns with only products or reactants ####

#is.unbalanced <- rep(NA, times = length(stoiMat[1,]))

#for(i in 1:length(stoiMat[1,])){
#	is.unbalanced[i] <- ifelse((length(stoiMat[,i][stoiMat[,i] > 0]) != 0) & (length(stoiMat[,i][stoiMat[,i] < 0]) != 0), FALSE, TRUE)
#	}

#rem.unbalanced <- colnames(stoiMat)[is.unbalanced]



#### Remove generic reactions - those which are not mass balanced or are generalizations of a class of species ####



#labelz <- c("isa", "protein production", "biomass production", "growth", "lipid production", "IPC synthase")
#aggregate_rxns <- NULL

#	aggregate_rxns <- union(aggregate_rxns, rxn_search(labelz[l], named_stoi, is_rxn = TRUE, index = TRUE))
#	}

#rem.aggregate <- colnames(stoiMat)[aggregate_rxns]
#rxn_search(named_stoi, labelz[l], is_rxn = TRUE, index = TRUE)

#look for rxns that produce CO2 and make them irreversible

#carb_match <- rxn_search(named_stoi, "carbon dioxide", is_rxn = FALSE, index = TRUE)

#co_two_producing_rx <- apply(stoiMat[met_dict == "carbon dioxide",carb_match] < 0, 2, sum) == 0
#co_two_producing_rx <- names(co_two_producing_rx)[co_two_producing_rx]

#reversibleRx$reversible[reversibleRx$rx %in% co_two_producing_rx] <- 1




#### output files ####



if(QPorLP == "QP"){
  growth_rate <- data.frame(cond = chemostatInfo$condition[1:n_c], limit = chemostatInfo$limitation[1:n_c], dr = chemostatInfo$actualDR[1:n_c], L1 = NA, L2 = NA)
  }

flux_vectors <- list()
#save(stoiMat, rxnFile, rxnparFile, corrFile, compFile, metComp, reversibleRx, comp_by_cond, nutrientFile, chemostatInfo, file = "condition_model_setup.Rdata") #save a .Rdata file to generate reaction formulae

######################## Set up the equality and inequality constriants for FBA ################

#remove reactions which are defective 
S_rxns = stoiMat#[,!(colnames(stoiMat) %in% c(rem.unbalanced, rem.aggregate))]


#split reactions which can carry forward and reverse flux into two identical reactions
validrxns <- reversibleRx[reversibleRx$rx %in% colnames(S_rxns),]
validrxnsFluxdir <- sapply(colnames(S_rxns), function(rxdir){
  validrxns$reversible[validrxns$rx == rxdir]
})

stoiRxSplit <- NULL
for(rxsplit in 1:length(validrxnsFluxdir)){
  if(unname(validrxnsFluxdir[rxsplit]) == 0){out <- rbind(c(paste(c(names(validrxnsFluxdir)[rxsplit], "F"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "F"), c(paste(c(names(validrxnsFluxdir)[rxsplit], "R"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "R"))}
  if(unname(validrxnsFluxdir[rxsplit]) == 1){out <- c(paste(c(names(validrxnsFluxdir)[rxsplit], "F"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "F")}  
  if(unname(validrxnsFluxdir[rxsplit]) == -1){out <- c(paste(c(names(validrxnsFluxdir)[rxsplit], "R"), collapse = '_'), names(validrxnsFluxdir)[rxsplit], "R")}  
  stoiRxSplit <- rbind(stoiRxSplit, out)
}
stoiRxSplit <- as.data.frame(cbind(stoiRxSplit)); colnames(stoiRxSplit) <- c("rxDesignation", "reaction", "direction")

S_rxns_split <- sapply(stoiRxSplit$reaction, function(rx){
  S_rxns[,colnames(S_rxns) == rx]
})
colnames(S_rxns_split) <- stoiRxSplit$rxDesignation

##### Stoichiometry and bounds of boundary fluxes #######

## Unconstrained Chemical Influx - e.g. water, gases ##

free_rxns <- sapply(freeExchange_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})

freeS <- matrix(0, ncol = length(free_rxns), nrow = length(metabolites))
for(i in 1:length(free_rxns)){
  freeS[free_rxns[i],i] <- 1
}

freeRxSplit <- data.frame(rxDesignation = c(paste(paste(c(freeExchange_met$SpeciesName), "boundary"), "F", sep = '_'), paste(paste(c(freeExchange_met$SpeciesName), "boundary"), "R", sep = '_')), 
      reaction = paste(c(freeExchange_met$SpeciesName, freeExchange_met$SpeciesName), "boundary"), direction = rep(c("F", "R"), each = length(freeExchange_met$SpeciesName))
      )

freeS_split <- cbind(freeS, freeS)
colnames(freeS_split) <- freeRxSplit$rxDesignation


## Nutrient Influx ##

nutrient_rxns <- sapply(boundary_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})
nutrientS <- matrix(0, ncol = length(nutrient_rxns), nrow = length(metabolites))
for(i in 1:length(nutrient_rxns)){
  nutrientS[nutrient_rxns[i],i] <- 1
}

nutrientRxSplit <- data.frame(rxDesignation = c(paste(c(boundary_met$SpeciesName), "boundary_offset"), paste(c(boundary_met$SpeciesName), "boundary_match_F"), paste(c(boundary_met$SpeciesName), "boundary_match_R")), 
      reaction = rep(paste(c(boundary_met$SpeciesName), "boundary"), times = 3), direction = c(rep("F", length(nutrient_rxns)*2), rep("R", length(nutrient_rxns))))

nutrientS_split <- cbind(nutrientS, nutrientS, nutrientS)
colnames(nutrientS_split) <- nutrientRxSplit$rxDesignation



## Excreted Metabolite Efflux - non-negative ##

efflux_rxns <- sapply(excreted_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})
effluxS <- matrix(0, ncol = length(efflux_rxns), nrow = length(metabolites))
for(i in 1:length(efflux_rxns)){
  effluxS[efflux_rxns[i],i] <- -1
}

effluxRxSplit <- data.frame(rxDesignation = c(paste(c(excreted_met$SpeciesName), "boundary_offset"), paste(c(excreted_met$SpeciesName), "boundary_match_F"), paste(c(excreted_met$SpeciesName), "boundary_match_R")), 
      reaction = rep(paste((excreted_met$SpeciesName), "boundary"), times = 3), direction = c(rep("F", length(efflux_rxns)*2), rep("R", length(efflux_rxns))))

effluxS_split <- cbind(effluxS, effluxS, effluxS)
colnames(effluxS_split) <- effluxRxSplit$rxDesignation


## Composition fxn ##

## generate a conversion matrix between each variance categories components and their indecies ##

biomassS <- matrix(0, ncol = length(unique(comp_by_cond$compositionFile$varCategory)), nrow = length(metabolites))

biomassRxSplit <- data.frame(rxDesignation = c(paste(c(unique(comp_by_cond$compositionFile$varCategory)), "comp_offset"), paste(c(unique(comp_by_cond$compositionFile$varCategory)), "comp_match_F"), paste(c(unique(comp_by_cond$compositionFile$varCategory)), "comp_match_R")), 
      reaction = rep(paste((unique(comp_by_cond$compositionFile$varCategory)), "composition"), times = 3), direction = c(rep("F", length(biomassS[1,])*2), rep("R", length(biomassS[1,]))))

biomassConv <- list()
for(a_rxn in biomassRxSplit$rxDesignation){
  biomassConv[[a_rxn]]$varCategory <- strsplit(a_rxn, split = ' comp')[[1]][1]
  category_species <- treatment_par[[1]][["boundaryFlux"]][[biomassConv[[a_rxn]]$varCategory]]$exchange$AltName
  category_species <- sapply(category_species, function(x){comp_met[comp_met$SpeciesName == x,]$SpeciesID[1]}) 
  met_row <- sapply(unname(category_species), function(x){c(1:length(metabolites))[metabolites == x] })
  biomassConv[[a_rxn]]$conversion <- data.frame(name = names(category_species), ID = unname(category_species), index = unname(met_row))
  
  }

biomassS_split <- cbind(biomassS, biomassS, biomassS)
colnames(biomassS_split) <- biomassRxSplit$rxDesignation





S <- cbind(S_rxns_split, freeS_split, nutrientS_split, effluxS_split, biomassS_split)

Sinfo <- rbind(stoiRxSplit, freeRxSplit, nutrientRxSplit, effluxRxSplit, biomassRxSplit)

S <- S * t(t(rep(1, times = nrow(S)))) %*% t(ifelse(Sinfo$direction == "R", -1, 1)) #invert stoichiometry for backwards flux

############ Gv >= h - bounds ########

## previous splitting of reversible reactions, allows restriction of each reaction's flux to be either non-negative or non-positive (depending on constrained direction) ##

#directG <- diag(ifelse(Sinfo$direction == "F", 1, -1))

## bounding maximum nutrient uptake rates ##

#influxG <- sapply(unique(nutrientRxSplit$reaction), function(rxchoose){
#  foo <- rep(0, length(Sinfo[,1]))
#  index_match <- c(1:length(Sinfo[,1]))[Sinfo$reaction == rxchoose]
#  foo[index_match] <- ifelse(Sinfo$direction[index_match] == "F", 1, -1)
#  foo
#  }) ### the rhs of these bounds are the maximal uptake rates, and the sense is <=

## since for the gurobi setup nutrient fluxes were split into 3 fluxes in order to accomidate the QP without an offset setup, the linear combination of uptake fluxes should instead be bounded
## this is difficult to do with this solver so instead for each molecule of a nutrient that is taken up an additional 'bookkeeping nutrient' will also be taken up.  The efflux of this bookkeeping
## nutrient through a single reaction can then be bounded by the empirical maximum nutrient influx rate, maintenance of global flux balance will enforce that the sum of actual nutrient influxes can't be higher
## than that of the directly bounded bookkeeping nutrient efflux

trackingS <- matrix(0, ncol =  length(S[1,]), nrow = length(unique(nutrientRxSplit$reaction)))
rownames(trackingS) <- paste(unique(nutrientRxSplit$reaction), "_bookkeeping", sep = "")

for(arxn in unique(nutrientRxSplit$reaction)){
  trackingS[rownames(trackingS) == paste(arxn, "_bookkeeping", sep = ""), Sinfo$reaction == arxn] <- ifelse(Sinfo$direction[Sinfo$reaction == arxn] == "F", 1, -1)
  }

## add in a reaction which siphons off each bookkeeping nutrient

S <- rbind(S, trackingS)

bookkeepingRx <- data.frame(rxDesignation = paste(unique(nutrientRxSplit$reaction), "_bookkeeping", sep = ""),
    reaction = paste(unique(nutrientRxSplit$reaction), "_bookkeeping", sep = ""), direction = "F")

bookkeepingS <- matrix(0, ncol = length(unique(nutrientRxSplit$reaction)), nrow = length(S[,1]))
colnames(bookkeepingS) <- bookkeepingRx$rxDesignation

for(arxn in unique(nutrientRxSplit$reaction)){
  bookkeepingS[rownames(S) == paste(arxn, "_bookkeeping", sep = ""),colnames(bookkeepingS) == paste(arxn, "_bookkeeping", sep = "")] <- -1
  }

Sinfo <- rbind(Sinfo, bookkeepingRx)
S <- cbind(S, bookkeepingS)

################ F - flux balance ############

Fzero <- rep(0, times = length(S[,1]))


########### Quadratic programming to match nutrient uptake/excretion rates and produce biomass ####

if(QPorLP == "QP"){

  library(Matrix) #for sparse matrix class
  library(gurobi) #interface for gurobi solver
  #ln -s /Library/gurobi510/mac64/lib/libgurobi51.so libgurobi51.so #a symbolic link was necessary to get gurobi to find its C++ code

  qpModel <- list()
  #qpparams <- list(Presolve=2, OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)
  qpparams <- list(OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)

  
  flux_elevation_factor <- 1000
  flux_penalty <- 1000/(flux_elevation_factor)
  
  qpModel$A <- Matrix(S)
  qpModel$rhs <- Fzero #flux balance
  qpModel$sense <- rep("=", times = length(S[,1])) #global flux balance
  qpModel$lb <- rep(0, times = length(S[1,])) #all fluxes are greater than zero
  qpModel$ub <- rep(Inf , times = length(S[1,])) #overwrite bookkeeping fluxes with empirical maximum rates
  
  qpModel$Q <- diag(rep(0, length(S[1,]))) #min t(v)Qv
  
  qpModel$obj <- rep(flux_penalty, length(S[1,])) #min c * v where v >= 0 to 
  qpModel$obj[grep('^r_', Sinfo$rxDesignation, invert = TRUE)] <- 0
  # lessen penalization for rxns with measured proteins and central C rxns
  qpModel$obj[qpModel$obj != 0] <- qpModel$obj[qpModel$obj != 0] * prot_penalty[chmatch(Sinfo$reaction[qpModel$obj != 0], names(prot_penalty))]
  
  ### QP-specific output_files ###
  
  residual_flux_stack <- NULL
  composition_balance <- NULL
  
  ### iterate through conditions and optimize the fit of fluxes to boundary and biomass conditions ###
  
  for(treatment in 1:n_c){
    
    cond_bound <- rep(Inf, times = length(S[1,]))
    cond_bound[Sinfo$reaction %in% treatment_par[[treatment]]$auxotrophies] <- 0 #auxotrophies have a maximum flux of zero
    
    cond_nutrients <- treatment_par[[treatment]]$nutrients[treatment_par[[treatment]]$nutrients$type == "uptake",]
    cond_nutrients$index <- sapply(paste(cond_nutrients$specie, "boundary_bookkeeping"), function(x){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == x]})
    
    cond_bound[cond_nutrients$index] <- cond_nutrients$ub #maximal nutrient fluxes set as [nutrient]*DR
    
    qpModel$ub <- cond_bound #hard bound the maximal flux through each reaction - Inf except for nutrient absorption
    
    ## constrain the offset fluxes to exactly equal the expected flux ##
    #lb = ub = 1 for composition, exchange rates for nutrient/excreted mets
    
    cond_boundary_rxns <- data.frame(Sinfo[grep('offset', Sinfo[,1]),], index = grep('offset', Sinfo[,1]))
    cond_boundary_rxns$rate <- NA
    for(nutrient in treatment_par[[treatment]]$nutrients$specie){
      cond_boundary_rxns$rate[cond_boundary_rxns$reaction == paste(nutrient, "boundary")] <- treatment_par[[treatment]]$nutrients$change[treatment_par[[treatment]]$nutrients$specie == nutrient]
      }
    cond_boundary_rxns$rate[cond_boundary_rxns$reaction %in% paste(unique(comp_by_cond$compositionFile$varCategory), "composition")] <- 1
    
    qpModel$lb[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    qpModel$lb[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- 0
    
    qpModel$ub[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    qpModel$ub[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- Inf
    
    ## quadratic matching of exchange fluxes and production of biomass components ##
    
    matchedSpecies <- data.frame(Sinfo[grep('match', Sinfo[,1]),], index = grep('match', Sinfo[,1]))
    matchedSpecies$Precision <- NA
    
    ## input the precision of each media specie - uptake and excretion ##
    
    for(nutrientSpec in treatment_par[[treatment]]$nutrients$specie){
      matchedSpecies$Precision[matchedSpecies$reaction == paste(nutrientSpec, "boundary")] <- (1/treatment_par[[treatment]]$nutrients$sd[treatment_par[[treatment]]$nutrients$specie == nutrientSpec])^2
      }
    
    stacked_comp_offset <- rbind(matchedSpecies[,1:4], cond_boundary_rxns[,1:4])
    
    for(biomassSpec in names(treatment_par[[treatment]]$boundaryFlux)){
      ## overwrite diagonal elements of Q with precision (inverse variance) ##
      matchedSpecies$Precision[matchedSpecies$reaction == paste(biomassSpec, "composition")] <- (1/treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$SD)^2
      
      ### overwrite stoichiometry for each biomass reaction to reflect the actual portions consumed such that the expected flux is 1.
      for(biomassSpecRx in names(biomassConv)[grep(paste('^', biomassSpec, sep = ""), names(biomassConv))]){
        
        comp_replacement <- data.frame(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange, biomassConv[[biomassSpecRx]])
        
        qpModel$A[comp_replacement$conversion.index,stacked_comp_offset$index[stacked_comp_offset$rxDesignation == biomassSpecRx]] <- comp_replacement$change*ifelse(sum(grep('_R$', biomassSpecRx)) != 0, -1, 1)
        
        }
    }
    #qpModel$A[rowSums(qpModel$A[,stacked_comp_offset$index] != 0) != 0,stacked_comp_offset$index]
    
    ## correct precisions for flux elevation factor
    matchedSpecies$Precision <- matchedSpecies$Precision / flux_elevation_factor^2
    
    ## for some species, no precision estimate is provided because they weren't directly measured.
    ##  By enforcing some arbitrary penalty on these reactions the unbounded stoichiometrically equivalent 'offset fluxes' should carry flux
    matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- 1
    
    diag(qpModel$Q)[matchedSpecies$index] <- matchedSpecies$Precision
    
    qpModel$lb <- qpModel$lb * flux_elevation_factor
    qpModel$ub <- qpModel$ub * flux_elevation_factor
    
    
   
    solvedModel <- gurobi(qpModel, qpparams) #solve with barrier algorithm
    
    if(solvedModel$status == "NUMERIC"){
      alt_parms <- qpparams; alt_parms$method <- 1
      solvedModel <- gurobi(qpModel, alt_parms) #solve with dual simplex
      }
    
    
    ### Optional flux checks ###
    
    #z <- maxFlux() # which reactions can carry flux (slow)
    #z <- loosenFlux(balanceStoi) # allow for free flux through a defined reaction of choice (with provided stoichiometry)
    #z$x[grep('r_0249', Sinfo$reaction)]
    #FF <- data.frame(rx = 'r_0745_F', flux = 1e1) #complex I ETC
    #FF <- data.frame(rx = 'r_0249_F', flux = 0.001) 
    #z <- forcedFlux(FF) #Force flux through a reaction to evaluate where it might be plugged up
    #ffAnalyze <- data.frame(rx = Sinfo$rxDesignation, force = z$x, reference = solvedModel$x)
    #ffAnalyze[,-1] <- apply(ffAnalyze[,-1], c(1,2), function(x){max(10^-10, x)})
    #ffAnalyze$log2Fluxchange <- log2(ffAnalyze$force) - log2(ffAnalyze$reference)
    #plot(log2(ffAnalyze$force) ~ log2(ffAnalyze$reference))
    #ffAnalyze[abs(ffAnalyze$log2Fluxchange) > 5,][order(ffAnalyze[abs(ffAnalyze$log2Fluxchange) > 5,]$log2Fluxchange),]
    
    
    
    #metFeed <- data.frame(mets = metabolites, obj = NA)
    #for(metN in 1:nrow(metFeed)){
    #  metFeedbalanceStoi <- data.frame(specie = metFeed$mets[metN], stoi = 1)
    #  metFeed$obj[metN] <- loosenFlux(balanceStoi, justObj = T)
    #  }
    #metIDtoSpec(metFeed$mets[metFeed$obj < 1])
    
    
    
    ### outputs ###
    
    collapsedFlux <- sapply(unique(Sinfo$reaction), function(frcombo){
      frindices <- c(1:length(Sinfo[,1]))[Sinfo$reaction == frcombo]
      sum(ifelse(Sinfo$direction[frindices] == "F", 1, -1) * solvedModel$x[frindices]/flux_elevation_factor)
    })
    
    # set fluxes below 10^-10 to zero, hist(log10(collapsedFlux)) indicates that there are many fluxes which are non-zero because of failure to round to zero
    collapsedFlux[abs(collapsedFlux) < 10^-10] <- 0
    
    # display contributions to L2 penalty 
    constrainedFlux <- data.frame(Sinfo[grep('match|book', Sinfo$rxDesignation),], flux = solvedModel$x[grep('match|book', Sinfo$rxDesignation)], Prec = diag(qpModel$Q)[grep('match|book', Sinfo$rxDesignation)],
    lb = qpModel$lb[grep('match|book', Sinfo$rxDesignation)], ub = qpModel$ub[grep('match|book', Sinfo$rxDesignation)])
    constrainedFlux$penalty <- constrainedFlux$flux^2 * constrainedFlux$Prec
    
    # deviations between allowable fluxes and empirical fluxes
    residualFlux <- data.table(reactions = names(collapsedFlux[grep('[boundary|composition]$', names(collapsedFlux))]), net_flux = unname(collapsedFlux[grep('[boundary|composition]$', names(collapsedFlux))]))
    
    experimental_flux <- rep(0, nrow(residualFlux))
    experimental_flux[chmatch(cond_boundary_rxns$reaction[!is.na(cond_boundary_rxns$rate)], residualFlux$reaction)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    
    residualFlux$experimental <- experimental_flux
    
    experimental_precision <- sapply(residualFlux$reaction, function(x){
      matchedSpecies$Precision[matchedSpecies$reaction == x][1]
    })
    experimental_precision[experimental_precision == 1] <- NA
    residualFlux$sd <- 1/sqrt(experimental_precision * flux_elevation_factor^2)
    residualFlux[,resid_st := (experimental - net_flux)/sd,]
    
    
    growth_rate$L1[treatment] <- sum(solvedModel$x*qpModel$obj)
    growth_rate$L2[treatment] <- sum(t(solvedModel$x) %*% qpModel$Q %*% t(t(solvedModel$x)))
    
    
    flux_vectors[[names(treatment_par)[treatment]]]$"flux" <- collapsedFlux
    flux_vectors[[names(treatment_par)[treatment]]]$"constraints" <- constrainedFlux
    flux_vectors[[names(treatment_par)[treatment]]]$"residual" <- residualFlux
    
    residual_flux_stack <- rbind(residual_flux_stack, cbind(condition = chemostatInfo$condition[treatment], residualFlux))
    
    ### Compare the elemental composition of all boundary constraints ### 
    
    boundary_label <- data.frame(reaction = residualFlux$reactions, color = NA)
    
    boundary_label$color[boundary_label$reaction %in% paste(free_flux, "boundary")] <- brewer.pal(sum(boundary_label$reaction %in% paste(free_flux, "boundary")), "Dark2")
    boundary_label$color[is.na(boundary_label$color)][grep('boundary', boundary_label$reaction[is.na(boundary_label$color)])] <- brewer.pal(length(grep('boundary', boundary_label$reaction[is.na(boundary_label$color)])), "Set3")
    boundary_label$color[grep('composition', boundary_label$reaction)] <- c(brewer.pal(length(grep('composition', boundary_label$reaction))-1, "Pastel2"), "brown1")
    
    boundary_stoichiometry <- qpModel$A[,sapply(residualFlux$reactions, function(x){(1:nrow(Sinfo))[Sinfo$reaction == x & Sinfo$direction == "F"][1]})]
    boundary_stoichiometry <- boundary_stoichiometry[rowSums(boundary_stoichiometry != 0) != 0,]
    boundary_stoichiometry <- boundary_stoichiometry[grep('bookkeeping', rownames(boundary_stoichiometry), invert = T),] 
    boundary_stoichiometry <- as.matrix(boundary_stoichiometry)
    
    boundary_elements <- convert_to_elemental(boundary_stoichiometry, modelMetComp)
    
    rownames(boundary_stoichiometry) <- boundary_elements$name
    ele_df <- NULL
    for(an_ele in c("C", "H", "N", "O", "P", "S")){
      an_ele_df <- data.frame(reactions = residualFlux$reactions, element = an_ele, netFlux = residualFlux$net_flux * t(boundary_stoichiometry) %*% boundary_elements[,colnames(boundary_elements) == an_ele],
      experimentalFlux = residualFlux$experimental * t(boundary_stoichiometry) %*% boundary_elements[,colnames(boundary_elements) == an_ele])
      
      an_ele_df$exchange <- ifelse(an_ele_df$netFlux <= 0, "Product", "Reactant") 
      an_ele_df$netFlux <- an_ele_df$netFlux * ifelse(an_ele_df$netFlux < 0, -1, 1) 
      an_ele_df$experimentalFlux <- an_ele_df$experimentalFlux * ifelse(an_ele_df$experimentalFlux < 0, -1, 1) 
      
      ele_df <- rbind(ele_df, an_ele_df)
      }
    ele_df$exchange <- factor(ele_df$exchange, levels = c("Reactant", "Product"))
    ele_df$color <- boundary_label$color[chmatch(ele_df$reactions, boundary_label$reaction)]
    ele_df$condition <- chemostatInfo$condition[treatment]
    
    composition_balance <- rbind(composition_balance, ele_df)
    
    }  
  } 

#reaction_info('r_0249')
#rxnFile[rxnFile$ReactionID == "r_0005",]
#trackedMet = '^NAD\\(\\+\\)'
#trackMetConversion('^2-oxoglutarate', T)


########### Mass balance of individual conditions ####### 

ele_df_melt <- data.table(melt(composition_balance, id.vars = c("reactions", "condition", "element", "exchange", "color")))
ele_df_melt[,limitation := chemostatInfo$limitation[chmatch(ele_df_melt$condition, chemostatInfo$condition)],]
ele_df_melt[,GR := chemostatInfo$DRgoal[chmatch(ele_df_melt$condition, chemostatInfo$condition)],]
ele_df_melt[,shortvar := factor(ifelse(variable == "netFlux", "OPT", "EXP"), levels = c("EXP", "OPT")),]


barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "bottom", 
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.key.width = unit(3, "line"),
    strip.background = element_rect(fill = "coral1"), panel.margin = unit(1, "lines"))

pdf(file = "massBalance.pdf", height = 15, width = 15)
for(an_element in unique(ele_df_melt$element)){
  elemental_mat <- subset(ele_df_melt,element == an_element,)
  
  involvedRxns <- elemental_mat[,sum(value != 0), by = reactions]
  elemental_mat <- elemental_mat[reactions %in% involvedRxns[V1 != 0,]$reactions,]
  
  elemental_barplot <- ggplot(elemental_mat, aes(x = exchange, y = value, fill = color)) + barplot_theme + facet_grid(limitation + shortvar ~ GR) + scale_fill_identity(name = "Class", guide = guide_legend(nrow = ifelse(length(unique(elemental_mat$reactions)) < 5, 2, 5)), labels = boundary_label$reaction[boundary_label$reaction %in% involvedRxns[V1 != 0,]$reactions], breaks = boundary_label$color[boundary_label$reaction %in% involvedRxns[V1 != 0,]$reactions]) +
    scale_y_continuous("Flux per element (moles / hour * mL cellular volume)") + scale_x_discrete("", expand = c(0,0)) + ggtitle(paste("Experimental and Optimized -- ", an_element, " -- mass balance"))
  print(elemental_barplot + geom_bar(stat = "identity", position = "stack"))  
  }
dev.off()


    


########### Generate a flux matrix of reactions carrying flux in some condition, as which deviations account for the L2 penalty #####

rxNames <- data.frame(reactionID = unique(Sinfo$reaction), Name = unique(Sinfo$reaction))
rxNames$Name[grep('r_[0-9]+', rxNames$Name)] = unname(rxnIDtoEnz(rxNames$Name[grep('r_[0-9]+', rxNames$Name)]))

#rxNames <- unique(Sinfo$reaction); rxNames[grep('r_[0-9]+', rxNames)] <- unname(rxnIDtoEnz(rxNames[grep('r_[0-9]+', rxNames)]))
fluxMat <- matrix(NA, ncol = n_c, nrow = length(flux_vectors[[1]]$flux)); colnames(fluxMat) <- names(flux_vectors); rownames(fluxMat) <- rxNames$Name

for(i in 1:n_c){
  fluxMat[,i] <- flux_vectors[[i]]$flux
  }
rxNames <- rxNames[rowSums(fluxMat) != 0,]
fluxMat <- fluxMat[rowSums(fluxMat) != 0,]

fluxMat_per_cellVol <- fluxMat / t(t(rep(1, length(fluxMat[,1])))) %*% chemostatInfo$VolFrac_mean[1:n_c] # moles per h*mL cell volume

####

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "top", 
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 20, angle = 70, vjust = 0.5, colour = "BLACK"), legend.key.width = unit(3, "line")) 


ggplot(residual_flux_stack[!is.na(residual_flux_stack$sd),], aes(x = factor(condition), y = resid_st*-1)) + facet_wrap(~ reactions, ncol = 3) + geom_bar(width = 0.75, stat = "identity") + barplot_theme + scale_y_continuous("predicted - experimental flux / sd(experimental)") + ggtitle("Fit of experimental and optimized flux")
ggsave("flux_residuals.pdf", width = 18, height = 20)
  

### for now only look at reactions which carry flux in >= 80% of conditions
#rxNames <- rxNames[rowSums(fluxMat_per_cellVol == 0) <= 5,]
#fluxMat_per_cellVol <- fluxMat_per_cellVol[rowSums(fluxMat_per_cellVol == 0) <= 5,]




std_flux <- fluxMat_per_cellVol / (t(t(apply(fluxMat_per_cellVol, 1, sd)*ifelse((apply(fluxMat_per_cellVol, 1, function(x){median(x[x != 0])}) >= 0), 1, -1))) %*% rep(1, n_c))

std_flux_rxnName <- std_flux
rownames(std_flux_rxnName) <- apply(rxNames, 1, function(x){ifelse(x[1] == x[2], x[1], paste(x, collapse = '_'))})
std_flux_rxnName <- std_flux_rxnName[grep('boundary|composition', rownames(std_flux_rxnName), invert = T),]



write.output(std_flux_rxnName, "Flux_analysis/fluxCarriedHM.tsv")

rawFlux <- fluxMat
rownames(rawFlux) <- rxNames$reactionID
rawFlux <- rawFlux[grep('bookkeeping', rownames(rawFlux), invert = T),]

write.table(rawFlux, file = "Flux_analysis/fluxCarriedSimple.tsv", quote = F, sep = "\t", col.names = TRUE, row.names = T)


npc <- 5
impSVD <- svd(std_flux_rxnName, nu = npc, nv = npc)
impSVD_pcs <- impSVD$v
colnames(impSVD_pcs) <- paste("PC", c(1:npc))
rownames(impSVD_pcs) <- colnames(std_flux_rxnName)

heatmap.2(t(impSVD_pcs), Colv = FALSE, Rowv = FALSE, trace = "none", col = greenred(100), dendrogram = "none", colsep = c(5,10,15,21), denscol = "white")

heatmap.2(std_flux, trace = "none", Colv = F)





flux_summary <- list()
flux_summary$IDs = rxNames
flux_summary$cellularFluxes = fluxMat_per_cellVol


save(flux_summary, file = "fluxSummaryQP.Rdata")




### visulaize flux through individual rxns ###

set.seed(1337)
rxnsChosen <- sample(nrow(std_flux), 10)

rxnsMelt <- melt(fluxMat_per_cellVol[rxnsChosen,])
colnames(rxnsMelt) <- c("Rxn", "Condition", "Flux")

scatter_theme <- barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
  panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 90), axis.line = element_blank(), strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "blue1")) 


ggplot(rxnsMelt, aes(x = factor(Condition), y = Flux)) + geom_point(size = 3, col = "coral1")+  facet_wrap( ~ Rxn, ncol = 2, scale = "free_y") + scatter_theme + scale_x_discrete("Condition") + scale_y_continuous("Flux Carried") + ggtitle("Flux through random reactions")







  






###### output fluxes so that they can be visualzied using S. cerevisae cellular overview #####

choice_conditions <- names(flux_vectors)#[grep(0.05, names(flux_vectors))]
cond_rownames <- names(flux_vectors[[1]]$flux)
cond_flux <- matrix(NA, ncol = length(choice_conditions), nrow = length(cond_rownames)); rownames(cond_flux) <- cond_rownames; colnames(cond_flux) <- choice_conditions

for(cond in choice_conditions){
  cond_flux[,colnames(cond_flux) == cond] <- flux_vectors[[c(1:length(flux_vectors))[names(flux_vectors) == cond]]]$flux
  }
cond_flux <- cond_flux[rowSums(cond_flux != 0) != 0,]

#which proteins are associated with each rxn
rxnEnzymes <- as.data.frame(rxnIDtoSGD(rownames(cond_flux)))
colnames(rxnEnzymes) <- c("compartment", "genes")

#visualize each compartment seperately
for(a_compartment in unique(rxnEnzymes$compartment)[!(unique(rxnEnzymes$compartment) %in% c("exchange", NA))]){
  comp_name <- compFile$compName[compFile$compID == a_compartment]
  comp_fluxes <- cond_flux[!is.na(rxnEnzymes$compartment) & rxnEnzymes$compartment == a_compartment,]
  comp_enzymes <- rxnEnzymes$genes[!is.na(rxnEnzymes$compartment) & rxnEnzymes$compartment == a_compartment]
  if(length(comp_enzymes) <= 1){next}
  
  #write each rxns flux, attributing flux to each associated enzyme
  
  #visualizing fluxes
  comp_outputDF <- NULL
  
  for(rxnN in c(1:length(comp_enzymes))){
    rxGenes = strsplit(comp_enzymes[rxnN], ':')[[1]]
    if(length(rxGenes) == 0){next}else{
      tmpMat <- matrix(comp_fluxes[rxnN,], ncol = length(choice_conditions), nrow = length(rxGenes), byrow = T)
      colnames(tmpMat) <- choice_conditions
      comp_outputDF <- rbind(comp_outputDF, data.frame(Enzyme = rxGenes, tmpMat))
      }
    }
  #flux relative to biomass objective
  #comp_outputDF[,-1] <- comp_outputDF[,-1]/(t(t(rep(1, length(comp_outputDF[,1])))) %*% sapply(choice_conditions, function(cmatch){growth_rate$growth[growth_rate$cond == cmatch]}))
  #max(abs(range(comp_outputDF[,-1])))
  colnames(comp_outputDF)[1] <- paste("$", colnames(comp_outputDF)[1], sep = "")
  
  #visualizing which reactions carry flux in a given condition
  ternary_outputDF <- comp_outputDF; ternary_outputDF[,-1][ternary_outputDF[,-1] < 0] <- -1; ternary_outputDF[,-1][ternary_outputDF[,-1] > 0] <- 1
  
  write.table(comp_outputDF, file = paste(c("SGDprojectionFiles/", comp_name, "fluxes.tsv"), collapse = ""), sep = "\t", row.names = F, col.names = T,  quote = F)
  write.table(ternary_outputDF, file = paste(c("SGDprojectionFiles/", comp_name, "TernaryFlux.tsv"), collapse = ""), sep = "\t", row.names = F, col.names = T,  quote = F)  
}







#colorz <- rep(c(1:5), each = 6)
#plot(growth_rate$growth, col = colorz, , xlab = "condition", ylab = "growth rate", pch = 16)
#legend("topleft", unique(growth_rate$limit), text.col = c(1:5))


#uncooperative rxns

#rxnSet <- reversibleRx[,1][reversibleRx[,1] %in% colnames(S)[apply((thermoG[c(1:length(thermoG[,1]))[!(c(1:length(thermoG[,1])) %in% c(1:28,30:32,34:35,37:44,47:58,60:64,66:71,73:74,76:77,79:82))],]) != 0, 2, sum) != 0]]

#rxnum <- 7
#rxnstoi <- S[S[,colnames(S) == rxnSet[rxnum]] != 0 ,colnames(S) == rxnSet[rxnum]]
#names(rxnstoi) <- metIDtoSpec(names(rxnstoi))
#rxnstoi

save(flux_vectors, growth_rate, treatment_par, file = "Flux_analysis/knitrFluxFilez.Rdata")
save(rxnFile, rxnparFile, corrFile, compFile, metComp, chemostatInfo, nutrientFile, reversibleRx, file = "Flux_analysis/knitrNetFilez.Rdata")











#heatmap of fluxes-per-unit growth

flux_per_gr <- reduced_flux_mat/matrix(growth_rate$growth, ncol = length(reduced_flux_mat[1,]), nrow = length(reduced_flux_mat[,1]), byrow = TRUE)
flux_per_gr <- flux_per_gr[apply(flux_per_gr, 1, sd) != 0,]

heatmap.2(t(scale(t(flux_per_gr), TRUE, TRUE)), trace = "n", cexRow = 0.05, col = green2red(1000))

####### BRIDGE TO NETWORK LAYOUT ########

  

# generate the stoichiometry matrix from rxns carrying flux
#load("totalStoi.Rdata")
update_layout = FALSE

Scollapse <- sapply(unique(Sinfo$reaction), function(x){
    first_match <- c(1:nrow(Sinfo))[Sinfo$reaction == x][1]
    S[,first_match]*ifelse(Sinfo$direction[first_match] == "F", 1, -1)
    })
Scollapse <- Scollapse[grep('bookkeeping', rownames(Scollapse), invert = T),grep('bookkeeping', colnames(Scollapse), invert = T)]



relevant_species <- list()
relevant_species$reactions = rxNames[grep('bookkeeping', rxNames$reactionID, invert = T),]
relevant_species$metabolites = rownames(Scollapse)[rowSums(Scollapse[,colnames(Scollapse) %in% rxNames$reactionID] != 0) != 0]

  

if(update_layout == TRUE){
  
  #updating an existing position file
  metSty.old <- read.delim("metSty.tsv", sep = "\t", header = TRUE)
  rxnSty.old <- read.delim("rxnSty.tsv", sep = "\t", header = TRUE)
  
  
  
  Stotal <- Scollapse[,colnames(Scollapse)[colnames(Scollapse) %in% union(relevant_species$reactions$reactionID, rxnSty.old$ReactionID)]]
  Stotal <- Stotal[apply(Stotal != 0, 1, sum) != 0,]
  
  #new rxns
  new_rxns <- relevant_species$reactions$reactionID[!(relevant_species$reactions$reactionID %in% rxnSty.old$ReactionID)]
  new_mets <- rownames(Stotal)[!(rownames(Stotal) %in% metSty.old$SpeciesID)]
  
  if(length(new_mets) != 0){
    
    metSty_bind <- as.data.frame(matrix(NA, ncol = length(metSty.old[1,]), nrow = length(new_mets)))
    colnames(metSty_bind) <- colnames(metSty.old)
    for(i in 1:length(new_mets)){
      metSty_bind[i,c(1:3)] <- corrFile[corrFile$SpeciesID %in% new_mets[i],][,c(1,2,4)]
    }
    metSty = rbind(metSty.old, metSty_bind)
    
  }else{
    metSty = rbind(metSty.old) 
  }
  
  if(length(new_rxns) != 0){
    
    rxnSty_bind <- as.data.frame(matrix(NA, ncol = length(rxnSty.old[1,]), nrow = length(new_rxns)))
    colnames(rxnSty_bind) <- colnames(rxnSty.old)
    for(i in 1:length(new_rxns)){
      if(new_rxns[i] %in% rxnFile$ReactionID){
        rxnSty_bind[i,c(1:3)] <- rxnFile[rxnFile$ReactionID %in% new_rxns[i],][1,][c(2,1,3)]
      }else{
        rxnSty_bind$ReactionID[i] <- new_rxns[i]
      }}
    rxnSty_bind$Reaction[is.na(rxnSty_bind$Reaction)] <- rxnSty_bind$ReactionID[is.na(rxnSty_bind$Reaction)]
    
    rxnSty = rbind(rxnSty.old, rxnSty_bind)
    
  }else{
    rxnSty = rxnSty.old
  }
	
  
  save(Stotal, metSty, rxnSty, reversibleRx, corrFile, relevant_species, file = "totalStoiAux.Rdata")

  write.table(metSty, file = "metSty.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(rxnSty, file = "rxnSty.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

  }







######### compare boundary fluxes with their limit #########

limiting_fluxes <- matrix(nrow = length(names(treatment_par)), ncol = length(treatment_par[[1]]$nutrients[,1])); rownames(limiting_fluxes) <- names(treatment_par); colnames(limiting_fluxes) <- treatment_par[[1]]$nutrients[,1]

for(treatment in 1:length(names(treatment_par))){
limiting_fluxes[treatment,]	<-(reduced_flux_mat[sapply(sapply(treatment_par[[treatment]]$nutrients$nutrient, function(x){paste(x, "boundary")}), function(x){c(1:length(reduced_flux_mat[,1]))[rownames(reduced_flux_mat) %in% x]}), treatment]*-1)/treatment_par[[treatment]]$nutrients$conc_per_t
	}
image(t(limiting_fluxes))
#library(xtable)
#xtable(limiting_fluxes)



	