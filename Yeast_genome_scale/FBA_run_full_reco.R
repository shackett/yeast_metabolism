#### Libraries ####
library(lpSolve)
library(limSolve)
library(gplots)
library(ggplot2)
source("FBA_lib.R")



#### Options ####

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

options(stringsAsFactors = FALSE)
inputFilebase = "yeast"

# Specify whether growth optimization should be through linear programming (LP) - maximizing biomass given nutrients or 
# through quadratic programming (QP) - optimally matching experimental boundary fluxes to optimized ones.
QPorLP <- "LP"

#calculate shadow prices to perform phenotypic phase plane analysis
shadow_prices = FALSE
if(QPorLP == "QP"){shadow_prices = FALSE} #don't calculate shadow prices when QP problem is being solved 
if(shadow_prices){
  treatmentPartials <- list()
  library(lpSolveAPI)
  generatePhPP <- TRUE; PhPPgenerated <- FALSE
}else{generatePhPP = TRUE}
# When optimization function changes, PhPP analysis becomes more complicated as composition becomes a function of nutrient conditions.  If experimental composition measurements over nutrient conditions were dense enough a
# surface could be fitted and an interesting PhPP analysis could be performed.  Analysis of shadow prices is still useful when partials are evaluated only at observed experimental conditions.
 


#### Load SBML files describing metabolites, rxn stoichiometry ####
rxnFile = read.delim(paste("rxn_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
rxnparFile = read.delim(paste("species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
corrFile = read.delim(paste("spec_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
compFile <- read.delim(paste("comp_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)

#### Load files describing boundary conditions and reaction reversibility from ecoli ####
metComp <- read.delim("METeleComp.tsv", stringsAsFactors = FALSE)
compositionFile <- read.csv2("../Yeast_comp.csv", sep = ",", stringsAsFactors = FALSE)
nutrientFile <- read.delim("Boer_nutrients.txt")[1:6,1:6]; nutrientCode <- data.frame(nutrient = colnames(nutrientFile)[-1], shorthand = c("n", "p", "c", "L", "u"))
rownames(nutrientFile) <- nutrientFile[,1]; nutrientFile <- nutrientFile[,-1]
reversibleRx <- read.delim("../EcoliYeastMatch/revRxns.tsv", sep = "\t", header = TRUE)
load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata") #load condition specific boundary fluxes and chemostat info (actual culture DR)


#### Load weizman free energy files and mapping designations ####
rxDesignations <- read.delim("../KEGGrxns/yeastNameDict.tsv", sep = "\t", header = T)
rxFreeEnergy <- read.delim("../KEGGrxns/kegg_reactions_PGC_ph5.0.csv", sep = ",", header = T); colnames(rxFreeEnergy) <- c("KEGGID", "freeEnergykJ_mol", "pH", "ionicStrength", "Note")
rxFreeEnergy$KEGGIDreformat <- sapply(rxFreeEnergy$KEGGID, function(KID){
    paste(c("K", rep(0, 5 - length(strsplit(as.character(KID), "")[[1]])), KID), collapse = "")
    })
reversibleRx$WEIZfreeE = NA; reversibleRx$WEIZdir = NA; reversibleRx$rxFlip = NA; reversibleRx$manual = NA; annotComment = NA

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)
enzymes <- rxnIDtoGene(reactions)

#### Load or write stoichiometry matrix of reactions and their altered metabolites ####

if(file.exists("yeast_stoi.R")){
  load("yeast_stoi.R")
} else {write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}



#### For reactions where proteins match multiple KEGG reactions, manually choose which is the proper match ####
if(!file.exists("manualKEGGrxns.Rdata")){
manualRxKEGGmatch <- NULL
for(rx in reactions){
  if(enzymes[names(enzymes) == rx] == ""){next}
  rxEnzymes <- strsplit(enzymes[names(enzymes) == rx], split = '/')[[1]]
  rxMatches <- rxFreeEnergy[rxFreeEnergy$KEGGIDreformat %in% unique(rxDesignations[rxDesignations$SYST %in% rxEnzymes,]$KEGG),]
  
  if(length(rxMatches$KEGGIDreformat) > 1){
    print(c("Multiple KEGG IDs match the following reaction"))
    print(reaction_info(rx))
    print(c("Which is the correct match?"))
    for(i in 1:length(rxMatches[,1])){
      print(paste(c(i, rxMatches$KEGGIDreformat[i], rxDesignations$NAME[rxDesignations$KEGG == rxMatches$KEGGIDreformat[i]][1], rxMatches$freeEnergykJ_mol[i]), collapse = " : "))
      }
    response <- as.numeric(readline(promp = "Well..."))
    
    manualRxKEGGmatch <- rbind(manualRxKEGGmatch, data.frame(reaction = rx, KEGG = rxMatches$KEGGIDreformat[response]))
    print("----------------------------------------")
    
    }
  #save(manualRxKEGGmatch, file = "manualKEGGrxns.Rdata")
}
}else{
  load("manualKEGGrxns.Rdata")
}

# Match unambiguous reactions to KEGG-associated free energy - multiple enzymes involved in a reaction or enzymes associated with multiple reactions degenerate this relationship - encorporate manual matching which was done on line 37
for(rx in reactions){
  if(enzymes[names(enzymes) == rx] == ""){next}
  rxEnzymes <- strsplit(enzymes[names(enzymes) == rx], split = '/')[[1]]
  rxMatches <- rxFreeEnergy[rxFreeEnergy$KEGGIDreformat %in% unique(rxDesignations[rxDesignations$SYST %in% rxEnzymes,]$KEGG),]
  if(length(rxMatches[,1]) == 1){
    reversibleRx$WEIZfreeE[reversibleRx$rx == rx] <- rxMatches$freeEnergykJ_mol
    
    }else{
    if(rx %in% manualRxKEGGmatch$reaction[!is.na(manualRxKEGGmatch$KEGG)]){
      reversibleRx$WEIZfreeE[reversibleRx$rx == rx] <-  rxMatches$freeEnergykJ_mol[rxMatches$KEGGIDreformat == manualRxKEGGmatch$KEGG[manualRxKEGGmatch$reaction == rx]]
    }}
}

#read in manually flipped and directed reactions (for flipped reations, the SM direction will be flipped and EC and Weiz direction will remain the same.  This is because SM (species matching) already flipped the free
#sign of the free energy of a reaction if substrates and products were reversed.  

thermAnnotate = read.delim("thermoAnnotate.txt", header = TRUE, sep = "\t")
for(rxN in 1:length(thermAnnotate[,1])){
  #flip reaction direction (and free energy) if stated directionality is unconventional  
  if(!is.na(thermAnnotate$flip[rxN]) & thermAnnotate$flip[rxN]){
    stoiMat[,colnames(stoiMat) == thermAnnotate$reaction[rxN]] <- stoiMat[,colnames(stoiMat) == thermAnnotate$reaction[rxN]]*-1
    reversibleRx[reversibleRx$rx == thermAnnotate$reaction[rxN],colnames(reversibleRx) %in% c("reversible", "SMfreeE")] <- reversibleRx[reversibleRx$rx == thermAnnotate$reaction[rxN],colnames(reversibleRx) %in% c("reversible", "SMfreeE")]*-1
    }
  
  #manually define reaction direction
  reversibleRx$rxFlip[reversibleRx$rx == thermAnnotate$reaction[rxN]] <- thermAnnotate$flip[rxN]
  reversibleRx$manual[reversibleRx$rx == thermAnnotate$reaction[rxN]] <- thermAnnotate$direction[rxN]
  }
reversibleRx$reversible[!is.na(reversibleRx$manual)] <- reversibleRx$manual[!is.na(reversibleRx$manual)]






#### Define the treatment in terms of nutrient availability and auxotrophies ####

dilution_rates <- seq(0.05, 0.3, 0.05)

treatment_par <- list()

for(i in 1:length(nutrientFile[1,])){
  for(j in 1:length(dilution_rates)){
	treatment_par[[paste(colnames(nutrientFile)[i], dilution_rates[j], collapse = "")]][["nutrients"]] <- data.frame(nutrient = rownames(nutrientFile), conc_per_t = nutrientFile[,i]*dilution_rates[j], stringsAsFactors = FALSE)
	
	#leu2
	if(colnames(nutrientFile)[i] == "Leucine"){
		treatment_par[[paste(colnames(nutrientFile)[i], dilution_rates[j], collapse = "")]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID))
		}
	#ura3
	if(colnames(nutrientFile)[i] == "Uracil"){
		treatment_par[[paste(colnames(nutrientFile)[i], dilution_rates[j], collapse = "")]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID))		}
	if(is.null(treatment_par[[paste(colnames(nutrientFile)[i], dilution_rates[j], collapse = "")]][["auxotrophies"]])){
		treatment_par[[paste(colnames(nutrientFile)[i], dilution_rates[j], collapse = "")]][["auxotrophies"]] <- NA
		}
	
	}}

treatment_par <- list()
n_c <- 25
for(i in 1:n_c){
  #define nutrient uptake rates (using maximal available for now)
  treatment_par[[chemostatInfo$condition[i]]][["nutrients"]] <- data.frame(nutrient = rownames(nutrientFile), conc_per_t = nutrientFile[,colnames(nutrientFile) == nutrientCode$nutrient[nutrientCode$shorthand == chemostatInfo$limitation[i]]]*chemostatInfo$actualDR[i])
  
  #define ura3 and leu2 auxotrophies
  if(chemostatInfo$limitation[i] == "u"){treatment_par[[chemostatInfo$condition[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID))}
  if(chemostatInfo$limitation[i] == "L"){treatment_par[[chemostatInfo$condition[i]]][["auxotrophies"]] <- as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID))}
  if(chemostatInfo$limitation[i] %in% c("c", "p", "n")){treatment_par[[chemostatInfo$condition[i]]][["auxotrophies"]] <- NA}
  
  #define observed fluxes per culture volume #eventually scale to the intracellular volume where these fluxes occur
  treatment_par[[chemostatInfo$condition[i]]][["boundaryFlux"]] = comp_by_cond$cultureMolarity[,colnames(comp_by_cond$cultureMolarity) == chemostatInfo$condition[i]]*chemostatInfo$actualDR[i]
  }

if(shadow_prices){
  #### If Phenotypic phase plane analysis (PhPP) is being performed, establish which condititions will be compared ####
  # which nutrients will be compared in a pairwise manner
  nutrientComp <- c("ammonium", "phosphate", "D-glucose")
  maxNutrientConc <- sapply(c(1:length(nutrientFile[,1])), function(x){nutrientFile[x, apply(nutrientFile, 1, which.max)[x]]})
  names(maxNutrientConc) <- rownames(nutrientFile); maxNutrientConc <- maxNutrientConc[rowSums(nutrientFile == 0) == 0]
  
  PhPPcond <- list()
  
  #generate all pairs of conditions
  nutrientPairs <- NULL
  for(i in 1:length(nutrientComp)){
    if(i == length(nutrientComp)){next} 
    for(j in (i+1):length(nutrientComp)){
      nutrientPairs <- rbind(nutrientPairs, c(nutrientComp[i], nutrientComp[j]))
    }}
  
  # generate concentration of each varying nutrient on an exponentially sampled grid
  conc_samples <- 20
  orderedNutrient <- sapply(nutrientComp, function(x){maxNutrientConc[names(maxNutrientConc) == x]}); names(orderedNutrient) <- sapply(names(orderedNutrient), function(y){strsplit(y, split = '\\.')[[1]][1]})
  nutrientGrad <- sapply(orderedNutrient, function(nconc){
    nconc * 1/exp((c(0:20)*(log(1000)/conc_samples)))
  })
  
  for(a_pairN in 1:length(nutrientPairs[,1])){
    nutLattice <- expand.grid(nutrientGrad[,colnames(nutrientGrad) == nutrientPairs[a_pairN,1]], nutrientGrad[,colnames(nutrientGrad) == nutrientPairs[a_pairN,2]])
    colnames(nutLattice) <- nutrientPairs[a_pairN,]
    
    invariantNut <- t(t(rep(1, length(nutLattice[,1])))) %*% maxNutrientConc[!(names(maxNutrientConc) %in% nutrientPairs[a_pairN,])]
    colnames(invariantNut) <- names( maxNutrientConc[!(names(maxNutrientConc) %in% nutrientPairs[a_pairN,])])
    
    PhPPcond$nutrients[[a_pairN]] <- nutrientPairs[a_pairN,]
    PhPPcond$gradient[[a_pairN]] <- cbind(nutLattice, invariantNut)
  }
}




#### Determine the compartmentation of each reaction ####

compartment <- sapply(reactions, function(x){rxnFile$Compartment[rxnFile$ReactionID == x][1]})

#### Define species involved in boundary-conditions ####

## extract the metabolite ID corresponding to the extracellular introduction of nutrients ##

sources <- c("D-glucose", "ammonium", "phosphate", "sulphate", "uracil", "L-leucine")
		
resource_matches <- lapply(sources, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

boundary_met <- NULL
for(x in 1:length(sources)){
boundary_met <- rbind(boundary_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


## extract the IDs of excreted metabolites ##

excreted <- c("acetate", "ethanol", "succinate(2-)", "(R)-lactate", "L-alanine", "L-glutamate")

resource_matches <- lapply(excreted, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

excreted_met <- NULL
for(x in 1:length(excreted)){
excreted_met <- rbind(excreted_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


## extract the metabolite ID corresponding to cytosolic metabolites being assimilated into biomass ##

sinks <- compositionFile$AltName

resource_matches <- lapply(sinks, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

comp_met <- NULL
for(x in 1:length(sinks)){
  comp_met <- rbind(comp_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "cytoplasm"],])
}

## freely exchanging metabolites through extracellular compartment ##

free_flux <- c("carbon dioxide", "oxygen", "water")

resource_matches <- lapply(free_flux, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

freeExchange_met <- NULL
for(x in 1:length(free_flux)){
  freeExchange_met <- rbind(freeExchange_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


skip_me = TRUE
######### Confirm mass balance of reactions ############
if(skip_me == FALSE){

met_chebi <- unlist(sapply(metabolites, metToCHEBI))
met_chebi_comp <- metComp[metComp$ID %in% met_chebi,]
met_chebi_comp <- met_chebi_comp[,c(TRUE, TRUE, apply(met_chebi_comp[,-c(1,2)], 2, sum) != 0)]

ele_comp <- lapply(met_chebi, function(x){
	if(length(met_chebi_comp[met_chebi_comp$ID %in% x,]) != 0){
	met_chebi_comp[met_chebi_comp$ID %in% x,]}
	})

#matrix of elemental abundance of species corresponding to rows of the stoichiometric matrix
ele_comp_mat <- matrix(NA, nrow = length(stoiMat[,1]), ncol = length(ele_comp[[1]])-2); rownames(ele_comp_mat) <- rownames(stoiMat); colnames(ele_comp_mat) <- names(ele_comp[[1]])[-c(1:2)]
for(el in 1:length(stoiMat[1,])){
	if(length(unlist(ele_comp[[el]][-c(1:2)])) != 0){
		ele_comp_mat[el,] <- unlist(ele_comp[[el]][-c(1:2)])
		}
	}

missed_compounds <- met_chebi[!is.na(met_chebi)][names(met_chebi[!is.na(met_chebi)]) %in% rownames(ele_comp_mat)[apply(is.na(ele_comp_mat), 1, sum) != 0]]

#create a data.frame of species with a chebi ID, but without a chemical formula match in my composition file
#missed_df <- data.frame(metIDtoSpec(names(missed_compounds)), names(missed_compounds), missed_compounds)
#unique_missed_df <- as.data.frame(matrix(NA, nrow = length(unique(missed_compounds)), ncol = 3))
#for(i in 1:length(unique(missed_compounds))){
#	unique_missed_df[i,] <- missed_df[missed_df[,3] == sort(unique(missed_compounds))[i],][1,]
#	}

#unique_missed_df[,3] %in% chem_form$COMPOUND_ID


#out <- as.data.frame(matrix(0, ncol = length(colnames(ele_comp_mat)), nrow = length(unique_missed_df[,1])))
#colnames(out) <- c(colnames(ele_comp_mat))

#write.table(cbind(unique_missed_df, out), row.names = FALSE, col.names = TRUE, file = "woot", sep = "\t")
add.chebi.comp <- read.delim("chebi_srh_curated.tsv", stringsAsFactors = FALSE)
add.chebi.comp <- add.chebi.comp[is.na(add.chebi.comp$generic),]
add.chebi.comp$internal_ID <- sapply(add.chebi.comp$internal_ID, function(x){corrFile$SpeciesType[corrFile$SpeciesID == x]})

for(i in 1:length(ele_comp_mat[,1])){
	
	if(corrFile$SpeciesType[corrFile$SpeciesID == rownames(ele_comp_mat)[i]] %in% add.chebi.comp$internal_ID){
		ele_comp_mat[i,] <- unlist(add.chebi.comp[add.chebi.comp$internal_ID %in% corrFile$SpeciesType[corrFile$SpeciesID == rownames(ele_comp_mat)[i]], -c(1:3, length(add.chebi.comp[1,]))])
		}}
		
		

mass_balanced <- matrix(NA, nrow = length(stoiMat[1,]), ncol = length(ele_comp_mat[1,]) + 1); rownames(mass_balanced) <- colnames(stoiMat); colnames(mass_balanced) <- c("missingIDs", colnames(ele_comp_mat))
mass_balanced <- as.data.frame(mass_balanced, stringsAsFactors = FALSE)

for(colnum in c(1:length(stoiMat[1,]))){
	
	if(length(stoiMat[,colnum][stoiMat[,colnum] != 0]) > 1){
	if(sum(is.na(ele_comp_mat[stoiMat[,colnum] != 0,][,1])) != 0){
		#some species are missing IDs
		#ignore those species and determine for which species the reaction is balanced
		mass_balanced$missingIDs[colnum] <- TRUE
	
		defined_spec <- !is.na(ele_comp_mat[stoiMat[,colnum] != 0,][,1])
		if(sum(defined_spec) > 1){
		mass_balanced[colnum, -1] <- t(ele_comp_mat[stoiMat[,colnum] != 0,][defined_spec,])%*% (stoiMat[,colnum][stoiMat[,colnum] != 0][defined_spec]) == 0
		}
		}else{
			
			mass_balanced$missingIDs[colnum] <- FALSE
			
			if(is.vector(ele_comp_mat[stoiMat[,colnum] != 0,]) == FALSE){
			mass_balanced[colnum, -1] <- t(ele_comp_mat[stoiMat[,colnum] != 0,]) %*% (stoiMat[,colnum][stoiMat[,colnum] != 0]) == 0
			}
			}}}
	
#reactions to use: those that have transformed metabolites
good_rxns <- colnames(stoiMat)[!is.na(mass_balanced[,2])]	
	
add_rxns <- mass_balanced[colnames(stoiMat) %in% good_rxns,][mass_balanced[colnames(stoiMat) %in% good_rxns,]$missingIDs == TRUE,]	
	
#reactions carrying flux with elemental composition information	
non_na_MB <- mass_balanced[rownames(mass_balanced) %in% rownames(reduced_flux_mat),][!is.na(mass_balanced[rownames(mass_balanced) %in% rownames(reduced_flux_mat),]$P),]

#reactions carrying flux that are not mass balanced for an element
non_na_MB_stoi <- stoiMat[apply(stoiMat[,colnames(stoiMat) %in% rownames(non_na_MB[non_na_MB$P == FALSE,])] != 0, 1, sum) != 0,colnames(stoiMat) %in% rownames(non_na_MB[non_na_MB$P == FALSE,])]

#for met
rxnparFile[rxnparFile[,1] == corrFile$SpeciesType[corrFile$SpeciesID == "s_0504"],]

mb.ids <- sapply(rownames(non_na_MB_stoi)[is.na(ele_comp_mat[rownames(ele_comp_mat) %in% rownames(non_na_MB_stoi),][,1])], function(x){corrFile$SpeciesType[corrFile$SpeciesID == x]})
#for rxn 
rxnFile[rxnFile$ReactionID %in% "r_1279",]
#destroy glycine-cleavage complex (lipoylprotein)

stoiMat <- stoiMat[,!is.na(mass_balanced[,2])]

rxnIDtoEnz(colnames(stoiMat)[is.na(mass_balanced[,2])])




rxnparFile[rxnparFile[,1] %in% unique(mb.ids),]

metIDtoSpec(rownames(non_na_MB_stoi)[is.na(ele_comp_mat[rownames(ele_comp_mat) %in% rownames(non_na_MB_stoi),][,1])])
metIDtoSpec("s_0504")

	
met_dict <- metIDtoSpec(rownames(non_na_MB_stoi))
rxn_dict <- rxnIDtoEnz(colnames(non_na_MB_stoi))

rownames(non_na_MB_stoi) <- sapply(c(1:length(non_na_MB_stoi[,1])), function(x){met_dict[x][[1]]})
colnames(non_na_MB_stoi) <- sapply(c(1:length(non_na_MB_stoi[1,])), function(x){rxn_dict[x][[1]]})

rxnparFile[rxnparFile[,1] %in% unique(mb.ids),]


}

#metIDtoSpec("s_0334")

#rxnparFile[,3][rxnparFile[,1] %in% corrFile$SpeciesType[corrFile$SpeciesID %in% "r_0267"]]

#taken from chebi
#chem_form <- read.delim("../Yeast_reconstruction/Sequences/chemical_data.tsv", stringsAsFactors = FALSE)
#chem_form <- chem_form[chem_form$SOURCE == "ChEBI",]


#chem_form[chem_form$ID %in% met_chebi,]
#rxnparFile[grep("7814", rxnparFile[,3]),]



#### Search for un-balanced rxns ######

is.unbalanced <- rep(NA, times = length(stoiMat[1,]))

for(i in 1:length(stoiMat[1,])){
	is.unbalanced[i] <- ifelse((length(stoiMat[,i][stoiMat[,i] > 0]) != 0) & (length(stoiMat[,i][stoiMat[,i] < 0]) != 0), FALSE, TRUE)
	}

rem.unbalanced <- colnames(stoiMat)[is.unbalanced]



#### Remove generic reactions - those which are not mass balanced or are generalizations of a class of species ####

named_stoi <- stoiMat
met_dict <- metIDtoSpec(rownames(named_stoi))
met_dict <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
rxn_dict <- rxnIDtoEnz(colnames(named_stoi))
rxn_dict <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})

rownames(named_stoi) <- met_dict
colnames(named_stoi) <- rxn_dict

labelz <- c("isa", "protein production", "biomass production", "growth", "lipid production", "IPC synthase")
aggregate_rxns <- NULL

for(l in 1:length(labelz)){
	aggregate_rxns <- union(aggregate_rxns, rxn_search(named_stoi, labelz[l], is_rxn = TRUE, index = TRUE))
	}

rem.aggregate <- colnames(stoiMat)[aggregate_rxns]




#look for rxns that produce CO2 and make them irreversible

carb_match <- rxn_search(named_stoi, "carbon dioxide", is_rxn = FALSE, index = TRUE)

co_two_producing_rx <- apply(stoiMat[met_dict == "carbon dioxide",carb_match] < 0, 2, sum) == 0
co_two_producing_rx <- names(co_two_producing_rx)[co_two_producing_rx]

reversibleRx$reversible[reversibleRx$rx %in% co_two_producing_rx] <- 1






growth_rate <- data.frame(limit = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[1]}), dr = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[2]}), growth = NA)

flux_vectors <- list()
######################## Set up the linear equations for FBA #######################

for(treatment in 1:length(names(treatment_par))){

#remove reactions which are defective 
S_rxns = stoiMat[,!(colnames(stoiMat) %in% c(treatment_par[[treatment]]$auxotrophies, rem.unbalanced, rem.aggregate))]

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
stoiRxSplit <- as.data.frame(cbind(stoiRxSplit, 0)); colnames(stoiRxSplit) <- c("rxDesignation", "reaction", "direction", "bound")

S_rxns_split <- sapply(stoiRxSplit$reaction, function(rx){
  S_rxns[,colnames(S_rxns) == rx]
  })
colnames(S_rxns_split) <- stoiRxSplit$rxDesignation

##### Stoichiometry and bounds of boundary fluxes #######

## Nutrient Influx & Unconstrained Chemical Influx ##

influx_rxns <- sapply(c(boundary_met$SpeciesID, freeExchange_met$SpeciesID), function(id){c(1:length(metabolites))[metabolites == id]})

influxS <- matrix(0, ncol = length(influx_rxns), nrow = length(metabolites))
for(i in 1:length(influx_rxns)){
	influxS[influx_rxns[i],i] <- 1
	}

influxSplit <- data.frame(rxDesignation = c(paste(paste(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName), "boundary"), "F", sep = '_'), paste(paste(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName), "boundary"), "R", sep = '_')), 
   reaction = paste(c(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName), c(boundary_met$SpeciesName, freeExchange_met$SpeciesName)), "boundary"), direction = rep(c("F", "R"), each = length(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName))), 
  bound = c(sapply(boundary_met$SpeciesName, function(bmet){treatment_par[[treatment]]$nutrients$conc_per_t[treatment_par[[treatment]]$nutrients$nutrient %in% bmet]}), rep(0, times = length(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName)) + length(freeExchange_met$SpeciesName)))
  )
                                                     
influxS_split <- cbind(influxS, influxS)
colnames(influxS_split) <- influxSplit$rxDesignation

## Excreted Metabolite Efflux - non-negative ##


efflux_rxns <- sapply(excreted_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})

effluxS <- matrix(0, ncol = length(efflux_rxns), nrow = length(metabolites))
for(i in 1:length(efflux_rxns)){
	effluxS[efflux_rxns[i],i] <- -1
	}

## Composition fxn ##

compVec <- rep(0, times = length(metabolites))
for(i in 1:length(comp_met$SpeciesID)){
  compVec[rownames(stoiMat) == comp_met$SpeciesID[i]] <- as.numeric(compositionFile$StoiCoef)[i]
}

sinkStoi <- cbind(effluxS, compVec); colnames(sinkStoi) <- c(paste(excreted_met$SpeciesName, "boundary"), "composition")
sinkSplit <- data.frame(rxDesignation = c(paste(excreted_met$SpeciesName, "boundary"), "composition"), reaction = c(paste(excreted_met$SpeciesName, "boundary"), "composition"), direction = "F", bound = 0)

S <- cbind(S_rxns_split, influxS_split, sinkStoi)
Sinfo <- rbind(stoiRxSplit, influxSplit, sinkSplit)


################ F - flux balance ############

Fzero <- rep(0, times = length(S[,1]))

############ Gv >= h - bounds ########

## previous splitting of reversible reactions, allows restriction of each reaction's flux to be either non-negative or non-positive (depending on constrained direction) ##

directG <- diag(ifelse(Sinfo$direction == "F", 1, -1))
directh <- rep(0, length(Sinfo[,1]))

## bounding maximum nutrient uptake rates ##


influxG <- t(sapply(c(1:length(Sinfo[,1]))[Sinfo$rxDesignation %in% paste(paste(boundary_met$SpeciesName, "boundary"), "F", sep = '_')], function(rxchoose){
  foo <- rep(0, length(Sinfo[,1])); foo[rxchoose] <- -1; foo
  }))
influxh <- -1 * as.numeric(Sinfo[Sinfo$rxDesignation %in% paste(paste(boundary_met$SpeciesName, "boundary"), "F", sep = '_'),]$bound)

Gtot <- rbind(directG, influxG)
htot <- c(directh, influxh) 	


############### costFxn - indicates the final rxn in S ######
#lp loss fxn
flux_penalty = 0.001
costFxn = c(rep(0, times = length(S[1,]) -1), -1) + ifelse(Sinfo$direction == "F", 1, -1) * flux_penalty

#qp loss fxn
#costFxn = c(rep(0, times = length(S[1,]) -1), 1)


######## use linear programming to maximize biomass #######

linp_solution <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = costFxn, ispos = FALSE)

#reconstruct reaction-specific flux from the F or R reaction which actually carried flux

collapsedFlux <- sapply(unique(Sinfo$reaction), function(frcombo){
  sum(linp_solution$X[names(linp_solution$X) %in% Sinfo$rxDesignation[Sinfo$reaction == frcombo]])
})

flux_vectors[[names(treatment_par)[treatment]]] <- collapsedFlux

growth_rate$growth[treatment] <- linp_solution$solutionNorm*-1



if(shadow_prices){
  #if the dual solution (shadow prices) are desired use an alternative solver which reports these values
  #these shadow prices indicate how much the addition of a metabolite change growth rate.
  
  costFxn = c(rep(0, times = length(S[1,]) -1), -1) #+ rep(1, times = length(S[1,])) * flux_penalty #removed flux penalty so that partials growth / 
  #partial metabolite change is the dual solution rather than this combined with flux resistance.
  
  #solve with lpSolveAPI
  pS <- S; pS <- pS * t(t(rep(1, times = length(S[,1])))) %*% t(ifelse(Sinfo$direction == "R", -1, 1)) #all v must be greater than zero
  lpObj <- make.lp(nrow = 0, ncol = length(S[1,]))
  for(i in 1:length(S[,1])){
    add.constraint(lpObj, pS[i,], "=", 0)
    }
  for(i in 1:length(influxh)){
    set.bounds(lpObj, upper = influxh*-1, columns = apply(influxG, 1, function(x){c(1:length(x))[x == -1]}))
    }
  set.objfn(lpObj, costFxn)
  
  #solve LP problem
  solve(lpObj)
  
  dual_solution <- get.dual.solution(lpObj)[1:length(S[,1])] 
  names(dual_solution) <- rownames(S)

  treatmentPartials[[names(treatment_par)[treatment]]] <- dual_solution
  
  #generate a phenotypic phase plane if one hasn't previously been generated
  if(generatePhPP & (PhPPgenerated == FALSE) & is.na(treatment_par[[treatment]]$auxotrophies)){
    
    for(nutComp in 1:length(PhPPcond$nutrients)){
      nutrientConditions <- PhPPcond$gradient[[nutComp]]
      
      nutIndices <- sapply(colnames(nutrientConditions), function(rxmatch){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == paste(paste(rxmatch, "boundary"), "F", sep = "_")]})
      
      
      for(nutCompConc in 1:lenght(nutrientConditions[,1])){
        set.bounds(lpObj, upper = nutrientConditions[nutCompConc,], columns = nutIndices)
        solve(lpObj)  
        get.dual.solution(lpObj)[1:length(S[,1])]
        get.objective(lpObj)
      }
      
    }
    
    
    
    
    
    }
  
  }


}





###### output fluxes so that they can be visualzied using S. cerevisae cellular overview #####

choice_conditions <- names(flux_vectors)[grep(0.05, names(flux_vectors))]
cond_rownames <- unique(unlist(sapply(choice_conditions, function(treatment){
names(flux_vectors[[c(1:length(flux_vectors))[names(flux_vectors) == treatment]]]) 
})))
cond_flux <- matrix(NA, ncol = length(choice_conditions), nrow = length(cond_rownames)); rownames(cond_flux) <- cond_rownames; colnames(cond_flux) <- choice_conditions


for(cond in choice_conditions){
  one_flux <- flux_vectors[[c(1:length(flux_vectors))[names(flux_vectors) == cond]]]
  cond_flux[sapply(names(one_flux), function(name_match){c(1:length(cond_rownames))[cond_rownames == name_match]}), choice_conditions == cond] <- unname(one_flux)
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
  comp_outputDF[,-1] <- comp_outputDF[,-1]/(t(t(rep(1, length(comp_outputDF[,1])))) %*% sapply(choice_conditions, function(cmatch){growth_rate$growth[rownames(growth_rate) == cmatch]}))
    #max(abs(range(comp_outputDF[,-1])))
  colnames(comp_outputDF) <- paste("$", colnames(comp_outputDF), sep = "")
  
  #visualizing which reactions carry flux in a given condition
  ternary_outputDF <- comp_outputDF; ternary_outputDF[,-1][ternary_outputDF[,-1] < 0] <- -1; ternary_outputDF[,-1][ternary_outputDF[,-1] > 0] <- 1
  
  write.table(comp_outputDF, file = paste(c("SGDprojectionFiles/", comp_name, "fluxes.tsv"), collapse = ""), sep = "\t", row.names = F, col.names = T,  quote = F)
  write.table(ternary_outputDF, file = paste(c("SGDprojectionFiles/", comp_name, "TernaryFlux.tsv"), collapse = ""), sep = "\t", row.names = F, col.names = T,  quote = F)  
}

suspects <- comp_outputDF[,1][abs(comp_outputDF[,3]) > 0.005]
suspectID <- rxnEnzymes[sapply(rxnEnzymes$genes, function(matcher){
  sum(strsplit(matcher, split = ':')[[1]] %in% suspects) != 0
  }),]
                          

acondFlux <- flux_vectors[[c(1:length(flux_vectors))[names(flux_vectors) == "Glucose 0.05"]]]


reaction_info("r_0708")
rxn_search(named_stoi, "glutamate deh", T)
qplot(comp_outputDF[,4])
heatmap.2(as.matrix(comp_outputDF[,-1]), symbreaks = TRUE)
comp_outputDF[,1][comp_outputDF[,4] < -60]


#### Analysing shadow prices ###

shadowCond <- sort(names(treatmentPartials))
shadowMets <- NULL; for(met in 1:length(shadowCond)){shadowMets <- union(shadowMets, names(treatmentPartials[[met]]))}; shadowMets <- sort(shadowMets)

shadowMetAbund <- matrix(NA, nrow = length(shadowMets), ncol = length(shadowCond)); colnames(shadowMetAbund) <- shadowCond; rownames(shadowMetAbund) <- shadowMets
for(j in shadowCond){
  condPrice <- treatmentPartials[[j]]
  shadowMetAbund[,colnames(shadowMetAbund) == j] <- condPrice[sapply(shadowMets, function(x){c(1:length(condPrice))[names(condPrice) == x]})]
  }


nonNULLshadow <- shadowMetAbund[rowSums(shadowMetAbund != 0) != 0,]
nonNULLshadow[nonNULLshadow < -3] <- -3; nonNULLshadow[nonNULLshadow > 3] <- 3
heatmap.2(nonNULLshadow, Colv = FALSE, trace = "n", dendrogram = "row", cexRow = 0.5, col = green2red(100))











#determine whether reactions can carry flux under the current boundary conditons

rxnGR <- data.frame(rxn = colnames(S), fluxF = NA, boundedF = NA, fluxR = NA, boundedR = NA)
x <- t(sapply(rxnGR[1:5,1], optimize_rx_flux))
rxnGR[,2:5] <- t(sapply(rxnGR[,1], optimize_rx_flux))

optimize_rx_flux <- function(rx){
	
	cost1 <- ifelse(colnames(S) %in% rx, -1, 0)
	linp_solution1 <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = cost1, ispos = FALSE)
	
	cost2 <- ifelse(colnames(S) %in% rx, 1, 0)
	linp_solution2 <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = cost2, ispos = FALSE)
	
	c(linp_solution1$solutionNorm*-1, linp_solution1$IsError, linp_solution2$solutionNorm*-1, linp_solution2$IsError)
	
	}


no_flux <- rxnGR$rxn[!((rxnGR$fluxF != 0 | rxnGR$boundedF == 1) | (rxnGR$fluxR != 0 | rxnGR$boundedR == 1))]
#rxns that can't carry flux 
rxnIDtoEnz(no_flux)


rxnGR$rxn[rxnGR$growth == 0 & rxnGR$bounded == 0]




#colorz <- rep(c(1:5), each = 6)
#plot(growth_rate$growth, col = colorz, , xlab = "condition", ylab = "growth rate", pch = 16)
#legend("topleft", unique(growth_rate$limit), text.col = c(1:5))


#uncooperative rxns

#rxnSet <- reversibleRx[,1][reversibleRx[,1] %in% colnames(S)[apply((thermoG[c(1:length(thermoG[,1]))[!(c(1:length(thermoG[,1])) %in% c(1:28,30:32,34:35,37:44,47:58,60:64,66:71,73:74,76:77,79:82))],]) != 0, 2, sum) != 0]]

#rxnum <- 7
#rxnstoi <- S[S[,colnames(S) == rxnSet[rxnum]] != 0 ,colnames(S) == rxnSet[rxnum]]
#names(rxnstoi) <- metIDtoSpec(names(rxnstoi))
#rxnstoi

save(flux_vectors, growth_rate, treatment_par, file = "Flux_analysis/SweaveFluxFilez.Rdata")
save(rxnFile, rxnparFile, corrFile, compFile, metComp, compositionFile, nutrientFile, reversibleRx, file = "Flux_analysis/SweaveNetFilez.Rdata")





all_rxns <- NULL
for(treatment in 1:length(names(treatment_par))){
	all_rxns <- union(all_rxns, names(flux_vectors[[names(flux_vectors)[treatment]]]))
	}

specified_rxns <- c(grep("r_", all_rxns), grep("boundary", all_rxns))

ordered_rxns <- all_rxns[c(specified_rxns, c(1:length(all_rxns))[!(all_rxns %in% all_rxns[specified_rxns])])]

cond_fluxes <- matrix(NA, ncol = length(names(treatment_par)), nrow = length(all_rxns))
rownames(cond_fluxes) <- ordered_rxns
colnames(cond_fluxes) <- names(treatment_par)

for(treatment in 1:length(names(treatment_par))){
	fluxes <- flux_vectors[[names(treatment_par)[treatment]]]

	cond_fluxes[,treatment] <- unlist(sapply(rownames(cond_fluxes), function(x){ifelse(is.null(fluxes[names(fluxes) %in% x]),NA,fluxes[names(fluxes) %in% x])}))

	}

reduced_flux_mat <- cond_fluxes[apply(cond_fluxes != 0, 1, sum) != 0,]
reduced_flux_mat <- reduced_flux_mat[!is.na(apply(reduced_flux_mat, 1, sd, na.rm = TRUE)),]
renamed_reduced_flux <- reduced_flux_mat; rownames(renamed_reduced_flux) <- c(rxnIDtoEnz(rownames(reduced_flux_mat)[grep("r_", rownames(reduced_flux_mat))]), rownames(reduced_flux_mat)[grep("r_", rownames(reduced_flux_mat), invert = TRUE)])

std.reduced_flux_mat <- (renamed_reduced_flux - apply(renamed_reduced_flux, 1, mean, na.rm = TRUE))/apply(renamed_reduced_flux, 1, sd, na.rm = TRUE)

#write.table(reduced_flux_mat, file = "carriedFlux.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

heatmap.2(std.reduced_flux_mat, Colv = FALSE, trace = "n", dendrogram = "row", cexRow = 0.05, col = green2red(1000))

heatmap.2(std.reduced_flux_mat[1:50,], Colv = FALSE, trace = "n", dendrogram = "row", cexRow = 0.5, col = green2red(1000))



#look at the subset of rxns carrying flux in ura and leu

renamed_reduced_flux[apply(renamed_reduced_flux[,!(growth_rate$limit %in% c("Leucine", "Uracil"))] != 0, 1, sum) == 0,]


#heatmap of fluxes-per-unit growth

flux_per_gr <- reduced_flux_mat/matrix(growth_rate$growth, ncol = length(reduced_flux_mat[1,]), nrow = length(reduced_flux_mat[,1]), byrow = TRUE)
flux_per_gr <- flux_per_gr[apply(flux_per_gr, 1, sd) != 0,]

heatmap.2(t(scale(t(flux_per_gr), TRUE, TRUE)), trace = "n", cexRow = 0.05, col = green2red(1000))

####### BRIDGE TO NETWORK LAYOUT ########



# generate the stoichiometry matrix from rxns carrying flux
load("totalStoi.Rdata")
first_write = FALSE

if(first_write == TRUE){
#writing a position file from scratch
Stotal <- S[,colnames(S) %in% rownames(reduced_flux_mat)]
Stotal <- Stotal[apply(Stotal != 0, 1, sum) != 0,]

metSty <- data.frame(SpeciesID = rep(NA, times = length(Stotal[,1])), SpeciesName = rep(NA, times = length(Stotal[,1])), Compartment = rep(NA, times = length(Stotal[,1])), x = rep(NA, times = length(Stotal[,1])), y = rep(NA, times = length(Stotal[,1]))) 
rxnSty <- data.frame(ReactionID = rep(NA, times = length(Stotal[1,])), Reaction = rep(NA, times = length(Stotal[1,])), Compartment = rep(NA, times = length(Stotal[1,])), x1 = rep(NA, times = length(Stotal[1,])), y1 = rep(NA, times = length(Stotal[1,])), x2 = rep(NA, times = length(Stotal[1,])), y2 = rep(NA, times = length(Stotal[1,]))) 

for(i in 1:length(Stotal[,1])){
	metSty[i,c(1:3)] <- corrFile[corrFile$SpeciesID %in% rownames(Stotal)[i],][,c(1,2,4)]
	}
for(i in 1:length(Stotal[1,])){
	if(colnames(Stotal)[i] %in% rxnFile$ReactionID){
		rxnSty[i,c(1:3)] <- rxnFile[rxnFile$ReactionID %in% colnames(Stotal)[i],][1,][c(2,1,3)]
		}else{
			rxnSty$ReactionID[i] <- colnames(Stotal)[i]
			}}
}else{
	#updating an existing position file
	metSty.old <- read.delim("metSty.tsv", sep = "\t", header = TRUE)
	rxnSty.old <- read.delim("rxnSty.tsv", sep = "\t", header = TRUE)
	
	Stotal <- S[,colnames(S) %in% union(rownames(reduced_flux_mat), rxnSty.old$ReactionID)]
	Stotal <- Stotal[apply(Stotal != 0, 1, sum) != 0,]
	
	#new rxns
	new_rxns <- rownames(reduced_flux_mat)[!(rownames(reduced_flux_mat) %in% rxnSty.old$ReactionID)]
	new_mets <- rownames(Stotal)[!(rownames(Stotal) %in% metSty.old$SpeciesID)]
	
	metSty_bind <- as.data.frame(matrix(NA, ncol = length(metSty.old[1,]), nrow = length(new_mets)))
	colnames(metSty_bind) <- colnames(metSty.old)
	for(i in 1:length(new_mets)){
		metSty_bind[i,c(1:3)] <- corrFile[corrFile$SpeciesID %in% new_mets[i],][,c(1,2,4)]
		}
	
	rxnSty_bind <- as.data.frame(matrix(NA, ncol = length(rxnSty.old[1,]), nrow = length(new_rxns)))
	colnames(rxnSty_bind) <- colnames(rxnSty.old)
	for(i in 1:length(new_rxns)){
	if(new_rxns[i] %in% rxnFile$ReactionID){
		rxnSty_bind[i,c(1:3)] <- rxnFile[rxnFile$ReactionID %in% new_rxns[i],][1,][c(2,1,3)]
		}else{
			rxnSty_bind$ReactionID[i] <- new_rxns[i]
			}}
	
	metSty = rbind(metSty.old, metSty_bind)
	rxnSty = rbind(rxnSty.old, rxnSty_bind)
	
	}

save(Stotal, metSty, rxnSty, reversibleRx, file = "totalStoiAux.Rdata")

write.table(metSty, file = "metSty.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(rxnSty, file = "rxnSty.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)






######### compare boundary fluxes with their limit #########

limiting_fluxes <- matrix(nrow = length(names(treatment_par)), ncol = length(treatment_par[[1]]$nutrients[,1])); rownames(limiting_fluxes) <- names(treatment_par); colnames(limiting_fluxes) <- treatment_par[[1]]$nutrients[,1]

for(treatment in 1:length(names(treatment_par))){
limiting_fluxes[treatment,]	<-(reduced_flux_mat[sapply(sapply(treatment_par[[treatment]]$nutrients$nutrient, function(x){paste(x, "boundary")}), function(x){c(1:length(reduced_flux_mat[,1]))[rownames(reduced_flux_mat) %in% x]}), treatment]*-1)/treatment_par[[treatment]]$nutrients$conc_per_t
	}
image(t(limiting_fluxes))
#library(xtable)
#xtable(limiting_fluxes)


########## evaluating rxns ########

#query_rxns <- c("protein production")

#eval_mat <- stoiMat[apply(stoiMat[,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns] != 0, 1, sum) != 0,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))

######### evaluate metabolite #########

#eval_mat <- stoiMat[,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))

reduced_flux_mat[rownames(reduced_flux_mat) %in% c("r_0393", "r_0394"),]


query_met <- c("stearoyl")

eval_mets(query_met, TRUE)

metToCHEBI("s_1002")
#
########






named_stoi <- stoiMat
met_dict <- metIDtoSpec(rownames(named_stoi))
rxn_dict <- rxnIDtoEnz(colnames(named_stoi))

rownames(named_stoi) <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
colnames(named_stoi) <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})

aggregate_rxns <- c(colnames(rxn_search(named_stoi, "isa", is_rxn = TRUE)), colnames(rxn_search(named_stoi, "protein production", is_rxn = TRUE)))

##########


#flux through regexed rxn
renamed_reduced_flux[rownames(renamed_reduced_flux) %in% names(rxn_search(named_stoi, "g", is_rxn = FALSE)[1,]),]



growth_rate$growth

"uracil boundary"

growth_rate$growth[30]/(treatment_par$'Uracil 0.3'$nutrients$conc_per_t[treatment_par$'Uracil 0.3'$nutrients$nutrient == "uracil"]/as.numeric(compositionFile[compositionFile$MetName %in% "UMP",]$StoiCoef))



colnames(stoiMat)



#### looking at extracellular rxns

extra <- stoiMat[,compartment == "c_05"]
extra_stoi <- extra[apply(extra != 0, 1, sum) != 0,]

met_dict <- metIDtoSpec(rownames(extra_stoi))
rxn_dict <- rxnIDtoEnz(colnames(extra_stoi))

rownames(extra_stoi) <- sapply(c(1:length(extra_stoi[,1])), function(x){met_dict[x][[1]]})
colnames(extra_stoi) <- sapply(c(1:length(extra_stoi[1,])), function(x){rxn_dict[x][[1]]})

#look at all reactions involving a metabolite from a particular compartment, that actually carry flux

extra <- stoiMat[rownames(stoiMat) %in% corrFile[corrFile$Compartment == "c_02",]$SpeciesID,apply(stoiMat[rownames(stoiMat) %in% corrFile[corrFile$Compartment == "c_02",]$SpeciesID,] != 0, 2, sum) != 0]

extra_red <- stoiMat[,colnames(stoiMat) %in% colnames(extra)[colnames(extra) %in% rownames(reduced_flux_mat)]]
extra_red <- extra_red[apply(extra_red != 0, 1, sum) != 0,]

met_dict <- metIDtoSpec(rownames(extra_red))
rxn_dict <- rxnIDtoEnz(colnames(extra_red))

rownames(extra_red) <- sapply(c(1:length(extra_red[,1])), function(x){met_dict[x][[1]]})
colnames(extra_red) <- sapply(c(1:length(extra_red[1,])), function(x){rxn_dict[x][[1]]})

 
 

#fxns


	
	