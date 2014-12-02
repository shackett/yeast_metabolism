setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

#### Libraries ####
library(gplots)
library(ggplot2)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(grid)
library(dplyr)

source("FBA_lib.R")

options(stringsAsFactors = FALSE)


#### Options ####

# Specify whether growth optimization should be through linear programming (LP) - maximizing biomass given nutrients or 
# through quadratic programming (QP) - optimally matching experimental boundary fluxes to optimized ones.
QPorLP <- "QP" # LP and PhPP no longer supported
if(QPorLP == "QP"){
  print("This is the standard run")
}else if (QPorLP == 'checkfeas'){
  print("This run identifies reactions and metabolites that are infeasible in the final S matrix")
}else{stop("try QP or checkfeas, all the cool kids are doing it")}

# if the gurobi solver is used via python, what kind of model should be solved
# which modus should be used for the problem
modeGurobi = 'python'
pythonMode = 'simple' # simple,thdyn, dir or ll (loopless)
FVA = 'T' # Should flux variblility analysis be performed
FVAreduced = 'T' # If FVA is performed should it only be performed on reactions which carried non-zero flux in some condition during QP
# this requires that the script be run twice - once to generate "rawFlux" and again to use this set of reactions to defined the reduced
# reaction stoichiometry
useCluster ='write' # can have 'F' for false, 'write' for write the cluster input or 'load' load cluster output


###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing valid reactions, species (their composition) both from the core SBML model and supplemented manual annotations #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

### Load SBML files describing metabolites, rxn stoichiometry ###

rxnFile = read.delim("companionFiles/rxn_yeast.tsv")
rxnparFile = read.delim("companionFiles/rxn_par_yeast.tsv")
corrFile = read.delim("companionFiles/spec_yeast.tsv")
specparFile = read.delim("companionFiles/species_par_yeast.tsv")
compFile = read.delim("companionFiles/comp_yeast.tsv")
fluxDirFile = read.delim("companionFiles/flux_dir_yeast.tsv")

### Add additional reactions ###

customList <- parse_custom("companionFiles/customRxns.txt")

rxnFile <- rbind(rxnFile, customList$rxnFile)
rxnparFile <- rbind(rxnparFile, customList$rxnparFile)
corrFile <- rbind(corrFile, customList$corrFile)
specparFile <- rbind(specparFile, customList$specparFile)
fluxDirFile <- rbind(fluxDirFile, customList$fluxDirFile)

save(rxnFile, rxnparFile, corrFile, fluxDirFile, fluxDirFile, file = "flux_cache/reconstructionWithCustom.Rdata")

### Determine unique metabolites and reactions ###

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)

#### Load or write stoichiometry matrix of reactions and their altered metabolites ####

if(!file.exists("flux_cache/yeast_stoi.Rdata")){write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}
load("flux_cache/yeast_stoi.Rdata")
if(nrow(stoiMat) != length(metabolites) | ncol(stoiMat) != length(reactions)){
  write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)
}

named_stoi <- stoiMat # create a version of the stoichiometric matrix with common names for rxns and mets to allow for easier searching
met_dict <- metIDtoSpec(rownames(named_stoi)); met_dict <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
rxn_dict <- rxnIDtoEnz(colnames(named_stoi)); rxn_dict <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})
rownames(named_stoi) <- met_dict; colnames(named_stoi) <- rxn_dict

### Elemental composition of metabolites ###

if(!file.exists("flux_cache/stoiMetsComp.tsv")){elemental_composition(metabolites)}
modelMetComp <- read.delim("flux_cache/stoiMetsComp.tsv", header = TRUE)
# re-run elemental_composition(metabolites) if the model changes and new metabolites are introduced



###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
##### Load files describing boundary conditions and reaction reversibility #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

#compositionFile <- read.delim("../Yeast_comp_energy.txt") #energy required to assimilate biomass components
nutrientFile <- read.delim("companionFiles/Boer_nutrients.txt")[1:6,1:6]; nutrientCode <- data.frame(nutrient = colnames(nutrientFile)[-1], shorthand = c("N", "P", "C", "L", "U"))
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

### append directionality with manual annotation in several cases ###
manualDirectionality <- read.delim("companionFiles/thermoAnnotate.txt")
reversibleRx$manual[chmatch(manualDirectionality$Reaction, reversibleRx$rx)] <- manualDirectionality$Direction
reversibleRx$reversible[!is.na(reversibleRx$manual)] <- reversibleRx$manual[!is.na(reversibleRx$manual)]


### add the component contribution dG predictions # used in thermodynamic flux analysis - not for polarizing reaction direction
ccPred <- read.delim("companionFiles/cc_dG_python.tsv") # generated from votti/vzcode/gen_input_comp_cont.R on github
reversibleRx[chmatch(ccPred$reaction, reversibleRx$rx), colnames(reversibleRx) %in% c("CCdG", "CCdGsd")] = ccPred[,-1]

### cache reaction directionality and thermodynamics
write.table(reversibleRx, file = "companionFiles/reversibleRx.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

#### Associate rxns with enzyme ascertainment in proteomics s.t. flux can favor measured pathways ####
### Associate proteins with reactions and format complexes

rxn_enzyme_groups = NULL
for(rxN in 1:nrow(rxnparFile)){
  
  enzymeCombos <- rxnparFile$Enzymes[rxN]
  
  if(enzymeCombos == ""){
    next
  } # no enzymes
  
  enzymeCombos <- strsplit(enzymeCombos, ' OR ')[[1]] # split complexes demarcated by OR
  
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

### cache files so that they are available when reaction equations are created ###

write.table(rxn_enzyme_groups, file = "flux_cache/rxn_enzyme_groups.tsv", sep = "\t", col.names = T, row.names = F, quote = F) # a data.frame indicating how proteins form catalytic units (monimers, dimers...)
write.table(data.frame(reaction = names(prot_matches), measured = unname(prot_matches)), file = "flux_cache/prot_matches.tsv", sep = "\t", col.names = T, row.names = F, quote = F) # boolean vector indicating whether a reaction's proteins were ascertained via proteomics




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

write.table(reactionMatches, "./flux_cache/reactionPathways.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

### Determine which pathways are associated with a protein ###
# Is this still used ?
#kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv") # > generated from keggIDparser.py KEGG IDs relative to yeast gene name (1gene -> 1 KEGG id, multiple  mapping between KEGG and genes)
#if(is.null(kegg_enzyme_dict$PATHWAY)){
#  gene_pathways(kegg_enzyme_dict); kegg_enzyme_dict <- read.delim("../KEGGrxns/yeastNameDict.tsv")
#}#generate a per-gene pathway annotation if one is not already generated
####



# favor flux through chosen central carbon metabolism pathways by reducing the L1 penalization of |v|

pathways <- sapply(reactionMatches$pathway, function(x){strsplit(x, '__')[[1]]})
unq_pathways <- unique(unlist(pathways))
centralC <- c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)")

centralCmatch <- unlist(lapply(pathways, function(x){
  sum(x %in% centralC) != 0
}))

centralCmatch[reactionMatches$reactionID %in% c("r_0713", "r_0962")] <- FALSE # remove MDH and PyK from this list to remove futile cycling

centralCrxnMatch <- rep(FALSE, times = length(reactions))
centralCrxnMatch[reactions %in% reactionMatches$reactionID[centralCmatch]] <- TRUE

ETCrxns <- c("r_0226", "r_0438", "r_0439", "r_0773", "r_1021", "r_0831", "r_0832") # electron transport reactions

centralCrxnMatch[reactions %in% ETCrxns] <- TRUE

# favor flux through cytosolic ATPase, + transport of ATP, ADP + Pi to restrict the wasting of excess energy to a single reaction

ATPbreakdownRxns <- c("r_4042", "r_1110", "r_1111", "r_1661", "r_3543", "r_3585", "r_3601", "r_3651", "r_3666",
                      "r_1244", "r_1245", "r_2005", "r_2008", "r_3537", "r_3605", "r_3649", "r_3663", "r_3940", "r_3961")

prot_penalty <- (1 - prot_matches)/2 + (1 - centralCrxnMatch)/2 # penalization by fraction of non-measured enzymes and favor central C metabolism
prot_penalty[reactions %in% ATPbreakdownRxns] <- 0.1

freeTransportRxns = c("r_1277", "r_2096", "r_1978", "r_1979", "r_1696", "r_1697") # transport of water, carbon dioxide and oxygen
prot_penalty[reactions %in% freeTransportRxns] <- 0



#### Define the treatment in terms of nutrient availability and auxotrophies ####

treatment_par <- list()
n_c <- nrow(chemostatInfo)
for(i in 1:n_c){
  #define nutrient uptake and excretion rate - soft matches on (using maximal available for now)
  #measured_bounds <- data.frame(nutrient = rownames(nutrientFile), conc_per_t = nutrientFile[,colnames(nutrientFile) == nutrientCode$nutrient[nutrientCode$shorthand == chemostatInfo$limitation[i]]]*chemostatInfo$actualDR[i])
  
  measured_bounds <- mediaSummary[mediaSummary$condition == chemostatInfo$ChemostatCond[i],]
  measured_bounds <- rbind(measured_bounds, data.frame(condition = chemostatInfo$ChemostatCond[i], specie = rownames(nutrientFile)[!(rownames(nutrientFile) %in% measured_bounds$specie)], change = NA, sd = NA, lb = 0, 
                                                       ub = nutrientFile[,colnames(nutrientFile) == nutrientCode$nutrient[nutrientCode$shorthand == chemostatInfo$Limitation[i]]][!(rownames(nutrientFile) %in% measured_bounds$specie)], type = "uptake", density = NA))
  
  # multiply steady-state concentrations by DR to get the uptake/excretion rates
  measured_bounds$change <- measured_bounds$change*chemostatInfo$actualDR[i]
  measured_bounds$sd <- measured_bounds$sd*chemostatInfo$actualDR[i]
  measured_bounds$lb <- measured_bounds$lb*chemostatInfo$actualDR[i]
  measured_bounds$ub <- measured_bounds$ub*chemostatInfo$actualDR[i]
  
  # remove phosphate because empirical uptake rates far exceed capacity of biomass assimilation
  measured_bounds <- data.frame(measured_bounds)
  measured_bounds[measured_bounds$specie == "phosphate [extracellular]", colnames(measured_bounds) %in% c("change", "sd")] <- NA
  measured_bounds <- data.table(measured_bounds)
  
  # define approximate oxygen uptake using a RQ of 5 for non-carbon-limited culture and < 5 for carbon-limited culture
  # employed by approximating vCO2 as 5/4 [ethanol + actetate]
  
  if(chemostatInfo$Limitation[i] == "C"){
    oxygen_uptake = (measured_bounds$change[measured_bounds$specie %in% c("D-glucose [extracellular]")]/5 + 
                       sum(measured_bounds$change[measured_bounds$specie %in% c("ethanol [extracellular]", "acetate [extracellular]")])/2)/2    
    
    oxygen_bounds <- data.frame(condition = chemostatInfo$ChemostatCond[i], specie = 'oxygen [extracellular]', change = 0, sd = Inf, 
                                lb = oxygen_uptake,
                                ub = Inf, type = "uptake", density = NA)
  } else {
    oxygen_uptake = (measured_bounds$change[measured_bounds$specie %in% c("D-glucose [extracellular]")]/5 + 
                       sum(measured_bounds$change[measured_bounds$specie %in% c("ethanol [extracellular]", "acetate [extracellular]")])/4)/2    
    
    oxygen_bounds <- data.frame(condition = chemostatInfo$ChemostatCond[i], specie = 'oxygen [extracellular]', 
                                change = oxygen_uptake,
                                sd = NA, lb = 0, ub = Inf, type = "uptake", density = NA)
    oxygen_bounds$sd <- oxygen_bounds$change/5
  }
  
  measured_bounds <- rbind(measured_bounds, oxygen_bounds)
  
  treatment_par[[chemostatInfo$ChemostatCond[i]]][["nutrients"]] <- measured_bounds
  
  #define ura3 and leu2 auxotrophies
  if(chemostatInfo$Limitation[i] == "L"){treatment_par[[chemostatInfo$ChemostatCond[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID))}
  if(chemostatInfo$Limitation[i] == "U"){treatment_par[[chemostatInfo$ChemostatCond[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID))}
  if(chemostatInfo$Limitation[i] %in% c("C", "P", "N")){treatment_par[[chemostatInfo$ChemostatCond[i]]][["auxotrophies"]] <- NA}
  
  #define observed fluxes per culture volume #scaled to intracellular volume later
  biomass_match <- data.frame(specie = comp_by_cond$compositionFile$MetName, AltName = comp_by_cond$compositionFile$AltName,change = unname(-1*comp_by_cond$cultureMolarity[,colnames(comp_by_cond$cultureMolarity) == chemostatInfo$ChemostatCond[i]]*chemostatInfo$actualDR[i]), index = comp_by_cond$compositionFile$index)
  biomass_list <- list()
  
  for(component in unique(comp_by_cond$compositionFile$varCategory)){
    principal_costs <- biomass_match[comp_by_cond$compositionFile$varCategory %in% component,]
    
    if(component == "Maintenance ATP hydrolysis"){
      total_costs <- principal_costs[,colnames(principal_costs) != "index"]
    }else{
      #costs of monomer assimilation incorporated into biomass flux
      energetic_costs <- as.matrix(comp_by_cond$biomassExtensionE[principal_costs$index,-1])
      energetic_costs_aggregate <- t(principal_costs$change) %*% energetic_costs; colnames(energetic_costs_aggregate) <- colnames(comp_by_cond$biomassExtensionE)[-1]
      total_costs <- rbind(principal_costs[,colnames(principal_costs) != "index"], data.frame(specie = colnames(energetic_costs_aggregate), AltName = colnames(energetic_costs_aggregate), change = t(unname(energetic_costs_aggregate)))[energetic_costs_aggregate != 0,])
    }
    biomass_list[[component]]$exchange = total_costs
    
    # define the accuracy of a constraint in terms of the coefficient of variation - sd over mean
    if(component %in% colnames(comp_by_cond$CV_table)){
      biomass_list[[component]]$SD = as.numeric(subset(comp_by_cond$CV_table, comp_by_cond$CV_table$condition == chemostatInfo$ChemostatCond[i], component))
    }else{
      biomass_list[[component]]$SD = 1/10
    }
  }
  
  treatment_par[[chemostatInfo$ChemostatCond[i]]][["boundaryFlux"]] = biomass_list
}

possibleAuxotrophies = c(as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID)), as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID)))



#### Determine the compartmentation of each reaction ####

compartment <- sapply(reactions, function(x){rxnFile$Compartment[rxnFile$ReactionID == x][1]})


#### Define species involved in boundary-conditions ####

## extract the metabolite ID corresponding to the extracellular introduction of nutrients ##
boundary_met <- treatment_par[[1]]$nutrients %>% filter(type == "uptake") %>% 
  select(SpeciesName = specie) %>% tbl_df() %>% left_join(corrFile)
if(nrow(boundary_met) != nrow(treatment_par[[1]]$nutrients %>% filter(type == "uptake"))){stop("A valid match was not found for all absorbed metabolites")}


## extract the IDs of excreted metabolites ##
excreted_met <- treatment_par[[1]]$nutrients %>% filter(type == "excretion") %>% 
  select(SpeciesName = specie) %>% tbl_df() %>% left_join(corrFile)
if(nrow(excreted_met) != nrow(treatment_par[[1]]$nutrients %>% filter(type == "excretion"))){stop("A valid match was not found for all excreted metabolites")}


## extract the metabolite ID corresponding to cytosolic metabolites being assimilated into biomass ##
sinks <- unique(c(comp_by_cond$compositionFile$AltName[comp_by_cond$compositionFile$FluxType == "Boundary"], colnames(comp_by_cond$biomassExtensionE)[-1]))
comp_met <- data.frame(SpeciesName = sinks) %>% tbl_df() %>% left_join(corrFile)
if(nrow(comp_met) != length(sinks)){stop("A valid match was not found for all sinks")}


## freely exchanging metabolites through extracellular compartment ##
free_flux <- c("carbon dioxide [extracellular]", "H2O [extracellular]", "H+ [extracellular]")
freeExchange_met <- data.frame(SpeciesName = free_flux) %>% tbl_df() %>% left_join(corrFile)


##  a list of reactions corresponding to a measured fluxes is formed ## 
## this method is not fully implemented but could be used to include MFA, KFP ... fluxes

interalFluxes <- read.delim("companionFiles/internalFluxMatches.csv", sep = ",")
interal_flux_list <- list()

for(i in c(1:nrow(interalFluxes))[!is.na(interalFluxes$Consumed_species)]){
  
  specSubset <- stoiMat[rownames(stoiMat) %in% corrFile$SpeciesID[grep(paste('^', interalFluxes$Consumed_species[i], sep = ""), corrFile$SpeciesName)],]
  specSubset <- specSubset[,colSums(specSubset < 0) != 0] # all reactions consuming the specie
  specSubset <- specSubset[,colSums(specSubset > 0) == 0] # remove transport reactions
  
  interal_flux_list[[interalFluxes$Measured_specie[i]]] <- data.frame(reactions = colnames(specSubset), coef = colSums(abs(specSubset)), positive = interalFluxes$is_positive[i])
}
for(i in c(1:nrow(interalFluxes))[!is.na(interalFluxes$Measured_reactions)]){
  interal_flux_list[[interalFluxes$Measured_specie[i]]] <- data.frame(reactions = strsplit(interalFluxes$Measured_reactions[i], split = ",")[[1]], coef = 1, positive = interalFluxes$is_positive[i])
}


#### Remove reactions which are incompatable with our method - e.g. an invariant biomass function ####

labelz <- c("biomass")
aggregate_rxns <- NULL
for(l in length(labelz)){
  aggregate_rxns <- union(aggregate_rxns, rxn_search(labelz[l], named_stoi, is_rxn = TRUE, index = TRUE))
}

## Remove lipid breakdown reactions ##
# 1) reactions are annotated as a hydrolase or lipase.  2) reactions produce a fatty acid annotated with the isa fatty acid tag

fatty_acids <- data.frame(name = c("myristate", "palmitate", "palmitoleate", "stearate", "oleate"), tID = NA)
fatty_acids$tID <- sapply(fatty_acids$name, function(x){unique(corrFile$SpeciesType[grep(paste('^', x, sep = ""), corrFile$SpeciesName)])})
all_FA <- corrFile$SpeciesID[corrFile$SpeciesType %in% fatty_acids$tID]

HLsubset <- stoiMat[metabolites %in% all_FA, grep('hydrolase|lipase', colnames(named_stoi))]
lipase_reactions <- colnames(HLsubset)[colSums(HLsubset) != 0]





###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
####### Setup matrices defining the stoichiometry of each reaction and how reactions will be constrained #######
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@



######################## Set up the equality and inequality constriants for FBA ################

# remove reactions which are defective 
S_rxns = stoiMat[,!(colnames(stoiMat) %in% c(aggregate_rxns, lipase_reactions))]

# remove reactions which cannot carry flux - see "Identify reactions that are infeasible"
if(file.exists("flux_cache/infeasibleRxnMet.txt") & QPorLP != 'checkfeas'){
  infRxMet <- read.table("flux_cache/infeasibleRxnMet.txt"); colnames(infRxMet) <- "ID"
  infRxMet <- rbind(infRxMet, 'r_0163') # remove several reactions which wont be used and cycle (ADH reverse)
  
  infRxMet$type <- substr(infRxMet$ID, 1, 1)
  
  reversibleRx <- reversibleRx[!(reversibleRx$rx %in% infRxMet$ID[infRxMet$type == "r"]),]
  metabolites <- metabolites[!(rownames(S_rxns) %in% infRxMet$ID[infRxMet$type == "s"])]
  S_rxns <- S_rxns[!(rownames(S_rxns) %in% infRxMet$ID[infRxMet$type == "s"]), !(colnames(S_rxns) %in% infRxMet$ID[infRxMet$type == "r"])]
  
}


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

## oxygen is directly coupled with ETC, to avoid un-physiological cycling ##
oxygen_stoi <- S_rxns[,grep('r_0438', colnames(S_rxns))]
oxygen_stoi[names(oxygen_stoi) == "s_1278"] <- 0
nutrientS[,boundary_met$SpeciesName == "oxygen [extracellular]"] <- oxygen_stoi # oxygen uptake == complex IV flux
S_rxns_split[grep('s_2817', rownames(S_rxns_split)), grep('r_218[23]', colnames(S_rxns_split))] <- 0 # make oxygen use for FA-desaturation a freebee.


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
## the stoichiometery will be populated on a by-condition basis ##

biomassS <- matrix(0, ncol = length(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), nrow = length(metabolites))

biomassRxSplit <- data.frame(rxDesignation = c(paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), "comp_offset"), 
                                               paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), "comp_match_F"),
                                               paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])), "comp_match_R")), 
                             reaction = rep(paste((unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])),"composition"), times = 3), 
                             direction = c(rep("F", length(biomassS[1,])*2), rep("R", length(biomassS[1,]))))

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

############ Gv >= h : bounds ########

## previous splitting of reversible reactions, allows restriction of each reaction's flux to be either non-negative or non-positive (depending on constrained direction) ##

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


## Add pseudoreactions which track some combination of fluxes ##
## rIDs are either identified by being pre-specified (Measured_reactions) or involve consumption of a specific specie (Consumed_species)

if(length(interal_flux_list) != 0){
  
  internal_trackingS <- matrix(0, nrow = length(interal_flux_list), ncol = ncol(S))
  colnames(internal_trackingS) <- colnames(S)
  rownames(internal_trackingS) <- paste(names(interal_flux_list), "usage_bookkeeping")
  
  for(i in 1:length(interal_flux_list)){
    Frxns <- interal_flux_list[[i]]$reactions[interal_flux_list[[i]]$positive]
    Revrxns <- interal_flux_list[[i]]$reactions[!interal_flux_list[[i]]$positive]
    matched_reactions <- c(1:nrow(Sinfo))[(Sinfo$reaction %in% Frxns & Sinfo$direction == "F") | Sinfo$reaction %in% Revrxns]
    
    internal_trackingS[i,matched_reactions] <- interal_flux_list[[i]]$coef[chmatch(Sinfo$reaction[matched_reactions], interal_flux_list[[i]]$reactions)]
    
  }
  
  S <- rbind(S, internal_trackingS)
  
  
  internalRxSplit <- data.frame(rxDesignation = c(paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])), "comp_offset"), 
                                                  paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])), "comp_match_F"),
                                                  paste(c(unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])), "comp_match_R")), 
                                reaction = rep(paste((unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])),"composition"), times = 3), 
                                direction = c(rep("F", length(interal_flux_list)*2), rep("R", length(interal_flux_list))))
  
  
  internal_bookkeeping <- matrix(0, nrow = nrow(S), ncol = nrow(internalRxSplit))
  colnames(internal_bookkeeping) <- internalRxSplit$rxDesignation
  
  for(j in 1:ncol(internal_bookkeeping)){
    internal_bookkeeping[rownames(S) %in% paste(strsplit(internalRxSplit$reaction[j], ' composition')[[1]], "usage_bookkeeping"), j] <- -1
  }
  
  Sinfo <- rbind(Sinfo, internalRxSplit)
  S <- cbind(S, internal_bookkeeping)
}

################ F : flux balance ############

Fzero <- rep(0, times = length(S[,1]))


############### Identify reactions that are infeasible by the structure of the S matrix ##########
if(!file.exists("flux_cache/infeasibleRxnMet.txt") | QPorLP == 'checkfeas'){
  ## There are compounds that are only formed and never consumed or only consumed but never produced
  
  S_tmp <- S # S_tmp needs to have temporary consumption rates filled in
  S_tmp_blank <- S_tmp[,colnames(S_tmp) %in% names(biomassConv)]
  for(j in 1:ncol(S_tmp_blank)){
    biomass_fill <- biomassConv[[colnames(S_tmp_blank)[j]]]$conversion
    biomass_fill$coef <- treatment_par[[1]]$boundaryFlux[[strsplit(colnames(S_tmp_blank)[j], split = ' comp')[[1]][1]]]$exchange$change
    S_tmp_blank[biomass_fill$index,j] <- biomass_fill$coef  
  }
  
  S_tmp_blank <- S_tmp_blank *  rep(1, nrow(S_tmp)) %*% t(ifelse(Sinfo$direction[colnames(S_tmp) %in% names(biomassConv)] == "F", 1, -1))
  
  S_tmp[,colnames(S_tmp) %in% names(biomassConv)] <- S_tmp_blank
  
  cpdfil_old = rep(T,nrow(S_tmp))
  cpdfil = rep(F,nrow(S_tmp))
  while (any(cpdfil != cpdfil_old)){
    cpdfil_old = cpdfil
    cpdfil = apply(S_tmp,1,function(x){all(x>=0)||all(x<=0)}) #metabolites which are only consumed or produced
    rxnfil = apply(S_tmp[cpdfil,],2,function(x){any(x!=0)}) #reactions involving these metabolites are flagged
    S_tmp[,rxnfil]=0 #flagged reactions are zeroed
  }
  rxnfil = apply(S_tmp,2,function(x){all(x==0)})
  infRxnMet = c(rownames(S_tmp)[cpdfil],unique(substr(colnames(S_tmp)[rxnfil & grepl('^r_',colnames(S_tmp))],1,6)))
  
  write(infRxnMet,file='flux_cache/infeasibleRxnMet.txt')
}

save(S, file = 'flux_cache/yeast_stoi_directed.Rdata')




###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@
####### Quadratic programming to match nutrient uptake/excretion rates and produce biomass #####
###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@###@

if(QPorLP == "QP"){
  
  library(Matrix) #for sparse matrix class
  library(gurobi) #interface for gurobi solver
  #ln -s /Library/gurobi510/mac64/lib/libgurobi51.so libgurobi51.so #a symbolic link was necessary to get gurobi to find its C++ code
  
  qpModel <- list()
  #qpparams <- list(Presolve=2, OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)
  qpparams <- list(OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)
  
  
  flux_elevation_factor <- 10^6
  flux_penalty <- 500/(flux_elevation_factor)
  #flux_penalty <- 1000/(flux_elevation_factor)
  
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
  qpModel$obj[qpModel$obj == 0] <- 10^-8 # reduce zero L1 penalties to a minute value, to aid in finding feasible solutions
  
  ### QP-specific output_files ###
  
  residual_flux_stack <- NULL
  composition_balance <- NULL
  fva_summary <- NULL
  growth_rate <- data.frame(cond = chemostatInfo$ChemostatCond, limit = chemostatInfo$Limitation, dr = chemostatInfo$actualDR, L1 = NA, L2 = NA)
  flux_vectors <- list()
  
  ### iterate through conditions and optimize the fit of fluxes to boundary and biomass conditions ###
  
  for(treatment in 1:n_c){
    
    # Define the maximal flux through reactions which are intrinsically constrained by the concentration of nutrients in the media, or are auxotrophies
    
    cond_bound <- rep(Inf, times = length(S[1,]))
    cond_bound[Sinfo$reaction %in% treatment_par[[treatment]]$auxotrophies] <- 0 #auxotrophies have a maximum flux of zero
    
    cond_nutrients <- treatment_par[[treatment]]$nutrients[treatment_par[[treatment]]$nutrients$type == "uptake",]
    cond_nutrients$index <- sapply(paste(cond_nutrients$specie, "boundary_bookkeeping"), function(x){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == x]})
    
    cond_bound[cond_nutrients$index] <- cond_nutrients$ub #maximal nutrient fluxes set as [nutrient]*DR
    
    qpModel$ub <- cond_bound #hard bound the maximal flux through each reaction - Inf except for nutrient absorption
    
    # Define lower bounds on flux through reactions - zero unless otherwise informed
    
    cond_min <- rep(0, times = length(S[1,]))
    cond_min[cond_nutrients$index] <- cond_nutrients$lb
    
    qpModel$lb <- cond_min
    
    # Setup quadratic penalty
    
    ## constrain the offset fluxes to exactly equal the expected flux ##
    # exchange rates for nutrient/excreted mets, lb = ub = Glucose uptake = fluxComp for composition
    
    cond_boundary_rxns <- data.frame(Sinfo[grep('offset', Sinfo[,1]),], index = grep('offset', Sinfo[,1]))
    cond_boundary_rxns$rate <- NA
    for(nutrient in treatment_par[[treatment]]$nutrients$specie){
      cond_boundary_rxns$rate[cond_boundary_rxns$reaction == paste(nutrient, "boundary")] <- treatment_par[[treatment]]$nutrients$change[treatment_par[[treatment]]$nutrients$specie == nutrient]
    }
    fluxComp = cond_boundary_rxns$rate[ cond_boundary_rxns$rxDesignation == 'D-glucose [extracellular] boundary_offset']
    cond_boundary_rxns$rate[cond_boundary_rxns$reaction %in% paste(unique(comp_by_cond$compositionFile$varCategory), "composition")] <- fluxComp
    
    qpModel$lb[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    qpModel$ub[cond_boundary_rxns$index][!is.na(cond_boundary_rxns$rate)] <- cond_boundary_rxns$rate[!is.na(cond_boundary_rxns$rate)]
    
    qpModel$lb[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- 0
    qpModel$ub[cond_boundary_rxns$index][is.na(cond_boundary_rxns$rate)] <- Inf
    
    ## quadratic matching of exchange fluxes and production of biomass components ##
    
    matchedSpecies <- data.frame(Sinfo[grep('match', Sinfo[,1]),], index = grep('match', Sinfo[,1]))
    matchedSpecies$Precision <- NA
    
    ## input the precision of each media specie - uptake and excretion ##
    
    for(nutrientSpec in treatment_par[[treatment]]$nutrients$specie){
      matchedSpecies$Precision[matchedSpecies$reaction == paste(nutrientSpec, "boundary")] <-
        (1/treatment_par[[treatment]]$nutrients$sd[treatment_par[[treatment]]$nutrients$specie == nutrientSpec])^2
    }
    
    stacked_comp_offset <- rbind(matchedSpecies[,1:4], cond_boundary_rxns[,1:4])
    
    # normalize all mass reactions to fluxComp (glucose uptake)
    for(biomassSpec in names(treatment_par[[treatment]]$boundaryFlux)){
      ## overwrite diagonal elements of Q with precision (inverse variance) ##
      matchedSpecies$Precision[matchedSpecies$reaction == paste(biomassSpec, "composition")] <- (1/(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$SD*fluxComp))^2
      biomassSpec_regex <- gsub('\\[', '\\\\[', biomassSpec)
      biomassSpec_regex <- gsub('\\]', '\\\\]', biomassSpec_regex)
      
      if(biomassSpec %in% unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Boundary"])){
        ### overwrite stoichiometry for each biomass reaction to reflect the actual portions consumed such that the expected flux is equal to fluxComp ##
        
        for(biomassSpecRx in names(biomassConv)[grep(paste0('^', biomassSpec_regex), names(biomassConv))]){
          
          comp_replacement <- data.frame(treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange, biomassConv[[biomassSpecRx]])
          
          qpModel$A[comp_replacement$conversion.index,stacked_comp_offset$index[stacked_comp_offset$rxDesignation == biomassSpecRx]] <- comp_replacement$change*ifelse(sum(grep('_R$', biomassSpecRx)) != 0, -1, 1) / fluxComp
        }
      }
      
      if(biomassSpec %in% unique(comp_by_cond$compositionFile$varCategory[comp_by_cond$compositionFile$FluxType == "Internal"])){
        ### overwrite internal fluxes ##
        for(biomassSpecRx in names(biomassConv)[grep(paste('^', biomassSpec_regex, sep = ""), names(biomassConv))]){
          
          comp_replacement <- internalRxSplit[internalRxSplit$reaction %in% paste(biomassSpec, "composition"),]
          comp_replacement$change <- treatment_par[[treatment]]$boundaryFlux[[biomassSpec]]$exchange$change * ifelse(comp_replacement$direction == "F", 1, -1) / fluxComp
          
          qpModel$A[rownames(S) == paste(biomassSpec, "usage_bookkeeping"), chmatch(internalRxSplit$rxDesignation[internalRxSplit$reaction %in% paste(biomassSpec, "composition")], colnames(S))] <- comp_replacement$change
          
        }
      }
      
    }
    
    #qpModel$A[rowSums(qpModel$A[,stacked_comp_offset$index] != 0) != 0,stacked_comp_offset$index]
    
    ## correct precisions for flux elevation factor
    matchedSpecies$Precision <- matchedSpecies$Precision / flux_elevation_factor^2
    
    ## for some species, no precision estimate is provided because they weren't directly measured.
    ## By enforcing some arbitrary penalty on these reactions the unbounded stoichiometrically equivalent 'offset fluxes' should carry flux
    matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- 1
    #matchedSpecies$Precision[is.na(matchedSpecies$Precision)] <- median(matchedSpecies$Precision,na.rm=T)
    
    diag(qpModel$Q)[matchedSpecies$index] <- matchedSpecies$Precision
    
    qpModel$lb <- qpModel$lb * flux_elevation_factor
    
    # Set the upper bounds to maximally maxMult times the glucose uptake rate
    maxMult = 100
    qpModel$ub <- sapply(qpModel$ub,function(x){min(maxMult*fluxComp, x)})
    qpModel$ub <- qpModel$ub * flux_elevation_factor
    
    
    ### Run standard quadratic programming ###
    
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
    
    # adjust composition fluxes to remove glucose uptake rate and restore units
    
    collapsedFlux[grep('composition', names(collapsedFlux))] <- sapply(grep('composition', names(collapsedFlux)), function(x){
      x <- collapsedFlux[x]
      x/fluxComp * abs(treatment_par[[names(treatment_par)[treatment]]]$"boundaryFlux"[[sub(' composition', '', names(x))]]$exchange$change[1])
    })
    
    
    # save summaries of flux carried and deviations between allowable and experimentally-measured fluxes
    flux_vectors[[names(treatment_par)[treatment]]]$"flux" <- collapsedFlux
    flux_vectors[[names(treatment_par)[treatment]]]$"constraints" <- constrainedFlux
    flux_vectors[[names(treatment_par)[treatment]]]$"residual" <- residualFlux
    
    residual_flux_stack <- rbind(residual_flux_stack, cbind(condition = chemostatInfo$ChemostatCond[treatment], residualFlux))
    
    ### Compare the elemental composition of all boundary constraints ### 
    
    boundary_label <- data.frame(reaction = residualFlux$reactions, color = NA)
    
    boundary_label$color[boundary_label$reaction %in% paste(free_flux, "[extracellular] boundary")] <- brewer.pal(sum(boundary_label$reaction %in% paste(free_flux, "[extracellular] boundary")), "Dark2")
    boundary_label$color[is.na(boundary_label$color)][grep('boundary', boundary_label$reaction[is.na(boundary_label$color)])] <- brewer.pal(length(grep('boundary', boundary_label$reaction[is.na(boundary_label$color)])), "Set3")
    boundary_label$color[grep('composition', boundary_label$reaction)] <- c(brewer.pal(8, "Pastel2"), brewer.pal(length(grep('composition', boundary_label$reaction)) - 8, "YlOrBr"))
    
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
    ele_df$condition <- chemostatInfo$ChemostatCond[treatment]
    
    composition_balance <- rbind(composition_balance, ele_df)
    
    
    ### Determine a range of feasible fluxes using python-gurobi based Flux variability analysis ###
    if (modeGurobi == 'python'){
      
      ### output will depend on useCluster:
      ### write: create files in Gurobi_python/test_files to run using qp_fba_clust.py
      ### load: read pythout.txt files containing solved FVA models
      ### otherwise... run locally: prohibitively slow for FVA
      
      # to submit all treatments to the que use "cetus_script.sh XXfolderXX" - where XXfolderXX is relative to the login directory
      # e.g. cetus_script.sh FBA/FBA_python/pythonVA
      
      # reduce FVA to the union of reactions carrying flux via standard QP - read previous reactions
      
      if(FVAreduced){
        if(file.exists('flux_cache/fluxCarriedSimple.tsv')){
          rx_with_flux <- rownames(read.delim('flux_cache/fluxCarriedSimple.tsv'))
          
          if( all(names(collapsedFlux)[collapsedFlux != 0] %in% rx_with_flux) ){
            
            validRx <- Sinfo$reaction %in% rx_with_flux
            validMets <- rowSums(qpModel$A[,validRx] != 0) != 0
            
            reduced_qpModel <- list()
            reduced_qpModel$A <- qpModel$A[validMets,validRx]
            reduced_qpModel$rhs <- qpModel$rhs[validMets]
            reduced_qpModel$sense <- qpModel$sense[validMets]
            reduced_qpModel$lb <- qpModel$lb[validRx]
            reduced_qpModel$ub <- qpModel$ub[validRx]
            reduced_qpModel$Q <- qpModel$Q[validRx,validRx]
            reduced_qpModel$obj <- qpModel$obj[validRx]
            
            solved_reduced_Model <- gurobi(reduced_qpModel, qpparams) #solve with barrier algorithm
            
            pythout = FVA_setup(reduced_qpModel, useCluster)
            
          }else{
            warning("\'fluxCarriedSimple.tsv\' is out of date, rerun after it is generated (finishing running script and it probably will be)")
          }
        }else{
          warning("FVA files will not be generated until \'fluxCarriedSimple.tsv\' has been generated (finishing running script and it probably will be)")
        }
        
      }else{
        pythout = FVA_setup(qpModel, useCluster)
      }
      
      if(is.null(pythout)){print(paste(treatment, "Gurobi files prepared for cluster use", collapse = " "))}else{
        
        modelOut = FVA_read(pythout)
        
        fva_sum <- modelOut[,grep('FVA[min|max]|asID|offset|runStatus|violation', colnames(modelOut))]
        
        ## add offset back in
        fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))] <- fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))] + t(t(ifelse(!is.na(fva_sum$offset), fva_sum$offset, 0))) %*% rep(1, length(grep('^FVA[min|max]', colnames(fva_sum))))
        ## remove flux elevation factor
        fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))] <- fva_sum[,grep('^FVA[min|max]', colnames(fva_sum))]/flux_elevation_factor
        fva_sum <- fva_sum[,colnames(fva_sum) != "offset"]
        fva_sum <- data.table(melt(fva_sum, id.vars = "asID"))
        
        fva_sum[,dataType := 
                  if(length(grep('^FVA[min|max]', variable)) != 0){"optim"}else if(length(grep('^runStatus', variable)) != 0){"status"}else if(length(grep('^violation', variable)) != 0){"violation"}
                , by = variable]
        
        levels(fva_sum$variable) <- sub('[A-Za-z.]+FVA', 'FVA', levels(fva_sum$variable))
        
        fva_sum <- data.table(dcast(fva_sum, "asID + variable ~ dataType", value.var = "value"))
        fva_sum$optim <- as.numeric(fva_sum$optim)
        
        fva_sum[, dataType:= strsplit(as.character(variable), split = "_")[[1]][1], by = variable]
        fva_sum[, boundType:= strsplit(as.character(variable), split = "_")[[1]][1], by = variable]
        fva_sum[, logLikelihood:= strsplit(as.character(variable), split = "_")[[1]][2], by = variable]
        fva_sum[, relLikelihood:= round(as.numeric(logLikelihood) - max(as.numeric(logLikelihood)), 5), by = logLikelihood]
        fva_sum$treatment = treatment
        fva_sum = fva_sum[,list(treatment, asID, boundType, relLikelihood, optim, status, violation),]
        
        fva_sum$optim[grep('comp', fva_sum$asID)] <- sapply(grep('comp', fva_sum$asID), function(x){
          x = data.frame(ID = sub(' comp', '', fva_sum$asID[x]), flux = fva_sum$optim[x])
          x[1,2]/fluxComp * abs(treatment_par[[names(treatment_par)[treatment]]]$"boundaryFlux"[[x[1,1]]]$exchange$change[1])
        })
        
        fva_summary <- rbind(fva_summary, fva_sum)
        
      }
    }
  }
} 





#reaction_info('r_0249')
#rxnFile[rxnFile$ReactionID == "r_0005",]
#trackedMet = '^NAD\\(\\+\\)'
#trackMetConversion('^trehalose', F)
#trackMetConversion('Gly\\-tRNA\\(Gly\\)', T)

########### Mass balance of individual conditions ####### 

ele_df_melt <- data.table(melt(composition_balance, id.vars = c("reactions", "condition", "element", "exchange", "color")))
ele_df_melt[,limitation := chemostatInfo$Limitation[chmatch(ele_df_melt$condition, chemostatInfo$ChemostatCond)],]
ele_df_melt[,GR := chemostatInfo$DRgoal[chmatch(ele_df_melt$condition, chemostatInfo$ChemostatCond)],]
ele_df_melt[,shortvar := factor(ifelse(variable == "netFlux", "OPT", "EXP"), levels = c("EXP", "OPT")),]


barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "bottom", 
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.key.width = unit(3, "line"),
    strip.background = element_rect(fill = "coral1"), panel.margin = unit(1, "lines"))

pdf(file = "massBalance.pdf", height = 15, width = 15)
for(an_element in unique(ele_df_melt$element)){
  elemental_mat <- subset(ele_df_melt,element == an_element,)
  
  involvedRxns <- elemental_mat[,sum(value != 0), by = reactions]
  elemental_mat <- elemental_mat[reactions %in% involvedRxns[V1 != 0,]$reactions,]
  
  elemental_barplot <- ggplot(elemental_mat, aes(x = exchange, y = value, fill = color)) + barplot_theme + facet_grid(limitation + shortvar ~ GR) + scale_fill_identity(name = "Class", guide = guide_legend(nrow = ifelse(length(unique(elemental_mat$reactions)) < 5, 2, 6)), labels = boundary_label$reaction[boundary_label$reaction %in% involvedRxns[V1 != 0,]$reactions], breaks = boundary_label$color[boundary_label$reaction %in% involvedRxns[V1 != 0,]$reactions]) +
    scale_y_continuous("Flux per atom (moles / hour * mL cellular volume)") + scale_x_discrete("", expand = c(0,0)) + ggtitle(paste("Experimental and Optimized -- ", an_element, " -- mass balance"))
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

#write.output(fluxMat_per_cellVol, "Flux_analysis/rawFluxHM.tsv")

#sort(apply(fluxMat_per_cellVol, 1, function(x){median(abs(x)[x != 0])}), decreasing = T)[1:200]

### observed RQ ###

#-1*fluxMat_per_cellVol[rownames(fluxMat_per_cellVol) == "carbon dioxide [extracellular] boundary"] / fluxMat_per_cellVol[rownames(fluxMat_per_cellVol) == "oxygen [extracellular] boundary"]

#carbon_dioxide_ex <- fluxMat_per_cellVol[rownames(fluxMat_per_cellVol) == "carbon dioxide [extracellular] boundary"]
#small_mol_excretion <- apply(fluxMat_per_cellVol[rownames(fluxMat_per_cellVol) %in% c("ethanol [extracellular] boundary", "acetate [extracellular] boundary"),], 2, sum)
#oxyen_used <- fluxMat_per_cellVol[rownames(fluxMat_per_cellVol) == "oxygen [extracellular] boundary_bookkeeping",]

#qplot(y = -1*carbon_dioxide_ex, x = small_mol_excretion, color = chemostatInfo$Limitation) +
#  geom_abline(intercept = 0, slope = 1)

#qplot(y = oxyen_used, x = -1*carbon_dioxide_ex, color = chemostatInfo$Limitation) +
#  geom_abline(intercept = 0, slope = 1)




###
std_flux <- fluxMat_per_cellVol / (t(t(apply(fluxMat_per_cellVol, 1, sd)*ifelse((apply(fluxMat_per_cellVol, 1, function(x){median(x[x != 0])}) >= 0), 1, -1))) %*% rep(1, n_c))

std_flux_rxnName <- std_flux
rownames(std_flux_rxnName) <- apply(rxNames, 1, function(x){ifelse(x[1] == x[2], x[1], paste(x, collapse = '_'))})
std_flux_rxnName <- std_flux_rxnName[grep('boundary|composition', rownames(std_flux_rxnName), invert = T),]

####
# add a column for the L1 penalty and write to a data.frame
std_flux_rxnName_withL1 = data.frame(L1penalty = -3*prot_penalty[chmatch(substr(rownames(std_flux_rxnName), 1, 6), names(prot_penalty))], std_flux_rxnName)
rownames(std_flux_rxnName_withL1) = rownames(std_flux_rxnName)
write.output(std_flux_rxnName_withL1, "flux_cache/fluxCarriedHM.tsv")

rawFlux <- fluxMat
rownames(rawFlux) <- rxNames$reactionID
#rawFlux <- rawFlux[grep('bookkeeping', rownames(rawFlux), invert = T),]

write.table(rawFlux, file = "flux_cache/fluxCarriedSimple.tsv", quote = F, sep = "\t", col.names = TRUE, row.names = T) # all reactions with non-zero flux in some condition


npc <- 5
impSVD <- svd(std_flux_rxnName, nu = npc, nv = npc)
impSVD_pcs <- impSVD$v
colnames(impSVD_pcs) <- paste("PC", c(1:npc))
rownames(impSVD_pcs) <- colnames(std_flux_rxnName)

heatmap.2(t(impSVD_pcs), Colv = FALSE, Rowv = FALSE, trace = "none", col = greenred(100), dendrogram = "none", colsep = c(5,10,15,20), denscol = "white")

heatmap.2(std_flux, trace = "none", Colv = F, col = greenred(100))




#### Determine which reactions have smooth flux patterns within each condition ####
# use F-test

smooth_flux_Fstat <- data.frame(rxn = rownames(fluxMat_per_cellVol), Fstat = NA, valid_flux = NA)
cond_model_matrix <- model.matrix(data = chemostatInfo, ~ factor(Limitation) + factor(Limitation)*actualDR)

for(i in 1:nrow(fluxMat_per_cellVol)){
  smooth_flux_Fstat$Fstat[i] <- anova(lm(fluxMat_per_cellVol[i,] ~ cond_model_matrix + 0))[1,5]
  }

smooth_flux_Fstat$valid_flux <- smooth_flux_Fstat$Fstat < 0.0001

std_flux_all <- fluxMat_per_cellVol / (t(t(apply(fluxMat_per_cellVol, 1, sd)*ifelse((apply(fluxMat_per_cellVol, 1, function(x){median(x[x != 0])}) >= 0), 1, -1))) %*% rep(1, n_c))

heatmap.2(std_flux_all[smooth_flux_Fstat$valid_flux,], Colv = FALSE, trace = "none", col = greenred(100), dendrogram = "none", colsep = c(5,10,15,20), denscol = "white")
heatmap.2(std_flux_all[!smooth_flux_Fstat$valid_flux,], Colv = FALSE, trace = "none", col = greenred(100), dendrogram = "none", colsep = c(5,10,15,20), denscol = "white")

write.output(std_flux_rxnName_withL1[smooth_flux_Fstat$valid_flux[grep('boundary|composition', rownames(std_flux_rxnName), invert = T)],], "flux_cache/fluxCarriedHMsmooth.tsv")



#### Summarize data for further analysis ####

flux_summary <- list()
flux_summary$IDs = rxNames
flux_summary$cellularFluxes = fluxMat_per_cellVol[,cond_order_check(colnames(fluxMat_per_cellVol))] # save cellular fluxes and verify correct ordering
flux_summary$QPresiduals <- residual_flux_stack
flux_summary$boundaryFluxes <- treatment_par

#### summarize flux variability analysis fluxes
#fva_summary_bu -> fva_summary
table(fva_summary$status) # 2 = optimal, 13 = suboptimal
fva_summary$optim[!(fva_summary$status %in% c(2,13))] <- NA # FVA results from optimization which did not find an optimal or suboptimal solution 

fva_dt <- data.table(fva_summary)
fva_dt_summary <- fva_dt[,list(diff = optim[boundType == "FVAmax"] - optim[boundType == "FVAmin"], avgFlux = (abs(optim[boundType == "FVAmax"]) + abs(optim[boundType == "FVAmin"]))/2, viol = max(violation)), by = c("treatment", "asID")]
fva_dt_summary <- fva_dt_summary[!is.na(diff) & asID %in% flux_summary$IDs$reactionID,,]
fva_dt_summary[,fractionDeparture := diff/avgFlux,]
byrxSummary <- fva_dt_summary[,list(departure = median(fractionDeparture), weighted_departure = sum((fractionDeparture * avgFlux)/sum(avgFlux))), by = "asID"]


fva_summary_df <- acast(fva_summary, formula = "asID ~ treatment ~ boundType", value.var = "optim")
cell_dens_scaling <- array(rep(chemostatInfo$VolFrac_mean[1:n_c], each = nrow(fva_summary_df)), dim = dim(fva_summary_df)) # scale flux-per-L culture/h, to flux-per-mL intracellular volume/h
fva_summary_df <- fva_summary_df/cell_dens_scaling

#### Determine which reactions' fluxes are sufficiently constrained to merit further analysis

constrained_rxns <- byrxSummary[departure < 0.3 | weighted_departure < 0.2, asID]
byrxSummary$plot_color <- ifelse(byrxSummary$asID %in% constrained_rxns, "RED", "BLUE")
ggplot(byrxSummary, aes(x = departure, y = weighted_departure, color = plot_color)) + geom_point() + 
  scale_color_identity() + scale_x_continuous("median[maximum flux - minimum flux / sup|max, min|]") + scale_y_continuous("median[maximum flux - minimum flux / sup|max, min|] - weighted by flux magnitude") +
  ggtitle("Flux variability at solution, scaled by magnitude")


fva_summary_df <- fva_summary_df[rownames(fva_summary_df) %in% intersect(constrained_rxns, flux_summary$IDs$reactionID[smooth_flux_Fstat$valid_flux]),,]
flux_summary$fva_summary_df = fva_summary_df
colnames(fva_summary_df) <- colnames(flux_summary$fva_summary_df) <- colnames(fluxMat_per_cellVol)[cond_order_check(colnames(fluxMat_per_cellVol))]


# verify comparable fluxes between standard QP and FVA approaches

sharedRx <- grep('^r_[0-9]{4}', intersect(flux_summary$IDs$reactionID[smooth_flux_Fstat$valid_flux], rownames(flux_summary$fva_summary_df)), value = T)

standard_QP <- flux_summary$cellularFluxes[flux_summary$IDs$reactionID %in% sharedRx,]; rownames(standard_QP) <- sharedRx
fva_QP <- fva_summary_df[rownames(fva_summary_df) %in% sharedRx,,]

SQPmelt <- melt(standard_QP); colnames(SQPmelt) <- c("rID", "condition", "flux"); SQPmelt$model = "standardQP"
FVAQPmelt <- melt(fva_QP); colnames(FVAQPmelt) <- c("rID", "condition", "model", "flux")
allFluxMelt <- rbind(SQPmelt[,chmatch(colnames(SQPmelt), colnames(FVAQPmelt))], FVAQPmelt)
compaFluxes <- dcast(allFluxMelt, formula = "rID + condition ~ model", value.var = "flux")

# is the flux calculated using the standard QP formulation captured between the FVA upper and lower bound
compaFluxes$fluxCap <- between(x = compaFluxes$standardQP, lower = compaFluxes$FVAmin, upper = compaFluxes$FVAmax, incbounds = T)
# fractional departure between QP solution and closest bound relative to ub-lb
compaFluxes$rangefluxDeparture <- apply(abs(compaFluxes$standardQP - data.frame(compaFluxes$FVAmax, compaFluxes$FVAmin))/(compaFluxes$FVAmax - compaFluxes$FVAmin), 1, min)
# same thing relative to maximum flux
compaFluxes$magfluxDeparture <- apply(abs(compaFluxes$standardQP - data.frame(compaFluxes$FVAmax, compaFluxes$FVAmin))/apply(data.frame(abs(compaFluxes$FVAmax), abs(compaFluxes$FVAmin)), 1, max, na.rm = T), 1, min)

total_flux_cast <- acast(allFluxMelt, formula = "rID ~ condition ~ model", value.var = "flux")
total_flux_cast <- total_flux_cast[,cond_order_check(colnames(total_flux_cast)),]
colnames(total_flux_cast) <- toupper(colnames(total_flux_cast))

# Dealing with some infinite fluxes (very few) introduced through FVA #

invalid_flux <- melt(apply(total_flux_cast[,,names(total_flux_cast[1,1,]) != "standardQP"], c(1,2), function(x){any(is.na(x))}))
invalid_flux <- invalid_flux[invalid_flux$value,]

for(i in 1:nrow(invalid_flux)){
  
  fva_rep <- total_flux_cast[rownames(total_flux_cast) == invalid_flux[i,1], colnames(total_flux_cast) == invalid_flux[i,2],]
  fva_rep[is.na(fva_rep)] <- fva_rep[names(fva_rep) == "standardQP"]
  
  total_flux_cast[rownames(total_flux_cast) == invalid_flux[i,1], colnames(total_flux_cast) == invalid_flux[i,2],] <- fva_rep
  
  }



#plot(sort(compaFluxes$magfluxDeparture[!compaFluxes$fluxCap])[1:3500])

# save flux summaries and pipe into FBGA.R #
flux_summary$total_flux_cast <- total_flux_cast

save(flux_summary, file = "flux_cache/fluxSummaryQP.Rdata")

# save some objects for use by the manuscript


save(stoiMat, rxnparFile, corrFile, compFile, metComp, chemostatInfo, nutrientFile, reversibleRx, file = "Flux_analysis/knitrNetFilez.Rdata")




### visulaize flux through individual rxns ###

set.seed(1337)
rxnsChosen <- sample(nrow(std_flux), 10)

rxnsMelt <- melt(fluxMat_per_cellVol[rxnsChosen,])
colnames(rxnsMelt) <- c("Rxn", "Condition", "Flux")

scatter_theme <- barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
  panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 90), axis.line = element_blank(), strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "blue1")) 


ggplot(rxnsMelt, aes(x = factor(Condition), y = Flux)) + geom_point(size = 3, col = "coral1") +  facet_wrap( ~ Rxn, ncol = 2, scale = "free_y") + scatter_theme + scale_x_discrete("Condition") + scale_y_continuous("Flux Carried") + ggtitle("Flux through random reactions")







  






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



	