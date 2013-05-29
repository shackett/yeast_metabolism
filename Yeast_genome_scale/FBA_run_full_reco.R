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
inputFilebase = "yeast"

# Specify whether growth optimization should be through linear programming (LP) - maximizing biomass given nutrients or 
# through quadratic programming (QP) - optimally matching experimental boundary fluxes to optimized ones.
QPorLP <- "QP"

#calculate shadow prices to perform phenotypic phase plane analysis
shadow_prices = FALSE
if(QPorLP == "QP"){shadow_prices = FALSE
  #ln -s /Library/gurobi510/mac64/lib/libgurobi51.so libgurobi51.so
  library(gurobi) # QP solver interface
  } #don't calculate shadow prices when QP problem is being solved 
if(shadow_prices){
  library(rgl)
  library(lpSolveAPI) # LP solver interface
  treatmentPartials <- list()
  generatePhPP <- TRUE; PhPPgenerated <- FALSE
}else{generatePhPP = FALSE}

# When optimization function changes, PhPP analysis becomes more complicated as composition becomes a function of nutrient conditions.  If experimental composition measurements over nutrient conditions were dense enough a
# surface could be fitted and an interesting PhPP analysis could be performed.  Analysis of shadow prices is still useful when partials are evaluated only at observed experimental conditions.
 


#### Load SBML files describing metabolites, rxn stoichiometry ####
rxnFile = read.delim(paste("rxn_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
rxnparFile = read.delim(paste("species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
corrFile = read.delim(paste("spec_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
compFile <- read.delim(paste("comp_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)

#### Load files describing boundary conditions and reaction reversibility from ecoli ####
metComp <- read.delim("METeleComp.tsv", stringsAsFactors = FALSE)
#compositionFile <- read.delim("../Yeast_comp_energy.txt") #energy required to assimilate biomass components
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

#named_stoi <- stoiMat
#rownames(named_stoi) <- metIDtoSpec(rownames(named_stoi))
#rxn_search(named_stoi, 'glycerol', is_rxn = FALSE)

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
  
  #define observed fluxes per culture volume #perhaps eventually scale to the intracellular volume where these fluxes occur
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
        biomass_list[[component]]$SD = 1/50
      }
    }
  
  
  treatment_par[[chemostatInfo$condition[i]]][["boundaryFlux"]] = biomass_list
  }
possibleAuxotrophies = c(as.character(unique(rxnFile[grep("isopropylmalate dehydrogenase", rxnFile$Reaction),]$ReactionID)), as.character(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID)))

#### During LP should the similarity of shadow prices across the nutrient landscape be investigated ####

if(generatePhPP){
  ## If Phenotypic phase plane analysis (PhPP) is being performed, establish which condititions will be compared ##
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
  conc_samples <- 100
  orderedNutrient <- sapply(nutrientComp, function(x){maxNutrientConc[names(maxNutrientConc) == x]}); names(orderedNutrient) <- sapply(names(orderedNutrient), function(y){strsplit(y, split = '\\.')[[1]][1]})
  nutrientGrad <- sapply(orderedNutrient, function(nconc){
    nconc * 1/exp((c(0:conc_samples)*(log(1000)/conc_samples)))
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



excreted <- unique(mediaSummary$specie[mediaSummary$type == "excretion"])

resource_matches <- lapply(excreted, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

excreted_met <- NULL
for(x in 1:length(excreted)){
excreted_met <- rbind(excreted_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


## extract the metabolite ID corresponding to cytosolic metabolites being assimilated into biomass ##

sinks <- unique(c(comp_by_cond$compositionFile$AltName, colnames(comp_by_cond$biomassExtensionE)))


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

#### Confirm mass balance of reactions ############
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

#write.table(data.frame(ID = metabolites, name = unname(metIDtoSpec(metabolites)),ele_comp_mat), file = "../ChemicalSpeciesQuant/stoiMetsComp.tsv", sep = "\t", col.names = TRUE, row.names = FALSE) #dump for boundary condition determination in boundaryDataBlender.R

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

#### Search for rxns with only products or reactants ####

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

stoiMat[,rxn_search(named_stoi, "growth", is_rxn = T, index = T)]


#look for rxns that produce CO2 and make them irreversible

carb_match <- rxn_search(named_stoi, "carbon dioxide", is_rxn = FALSE, index = TRUE)

co_two_producing_rx <- apply(stoiMat[met_dict == "carbon dioxide",carb_match] < 0, 2, sum) == 0
co_two_producing_rx <- names(co_two_producing_rx)[co_two_producing_rx]

reversibleRx$reversible[reversibleRx$rx %in% co_two_producing_rx] <- 1


#save(stoiMat, rxnFile, rxnparFile, corrFile, compFile, metComp, reversibleRx, comp_by_cond, nutrientFile, chemostatInfo, file = "condition_model_setup.Rdata") #save a .Rdata file to generate reaction formulae

#### output files ####


if(QPorLP == "LP"){
  growth_rate <- data.frame(cond = chemostatInfo$condition[1:n_c], limit = chemostatInfo$limitation[1:n_c], dr = chemostatInfo$actualDR[1:n_c], growth = NA)
  }
if(QPorLP == "QP"){
  growth_rate <- data.frame(cond = chemostatInfo$condition[1:n_c], limit = chemostatInfo$limitation[1:n_c], dr = chemostatInfo$actualDR[1:n_c], L1 = NA, L2 = NA)
  }

flux_vectors <- list()
  

######################## Set up the equality and inequality constriants for FBA ################

#remove reactions which are defective 
S_rxns = stoiMat[,!(colnames(stoiMat) %in% c(rem.unbalanced, rem.aggregate))]

colnames(stoiMat)

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
  
  #conversionMat <- matrix(0, ncol = length(met_row), nrow = length(metabolites))
  #for(j in 1:length(met_row)){
  #  conversionMat[met_row[j],j] <- 1
  #  }
  
  #biomassConv[[a_rxn]]$conversionMatrix <- conversionMat

  }

biomassS_split <- cbind(biomassS, biomassS, biomassS)
colnames(biomassS_split) <- biomassRxSplit$rxDesignation





S <- cbind(S_rxns_split, freeS_split, nutrientS_split, effluxS_split, biomassS_split)

Sinfo <- rbind(stoiRxSplit, freeRxSplit, nutrientRxSplit, effluxRxSplit, biomassRxSplit)

S <- S * t(t(rep(1, times = length(S[,1])))) %*% t(ifelse(Sinfo$direction == "R", -1, 1)) #invert stoichiometry for backwards flux

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


  
########### Linear programming to maximize growth given nutrient availability and calculate dual solution / shadow prices #######################

if(QPorLP == "LP"){
  
  library(lpSolveAPI) # interface to lpSolve
  
  lpObj <- make.lp(nrow = 0, ncol = length(S[1,]))
  for(i in 1:length(S[,1])){
    add.constraint(lpObj, S[i,], "=", 0)
  }
  
  for(treatment in 1:n_c){
  
    delete.column(lpObj, length(S[1,])) #delete the current composition function so that it can be overwritten
    
    compVec <- rep(0, times = length(metabolites))
    for(i in 1:length(comp_met$SpeciesID)){
      compVec[rownames(stoiMat) == comp_met$SpeciesID[i]] <- treatment_par[[treatment]]$"boundaryFlux"[i] 
    }  
    
    add.column(lpObj, compVec) #add back the condition-specific molar anabolic rates
  
    cond_nut_bound <- data.frame(index = sapply(paste(paste(boundary_met$SpeciesName, "boundary"), "F", sep = '_'), function(x){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == x]}), treatment_par[[treatment]]$nutrients)
    
    set.bounds(lpObj, upper = cond_nut_bound$conc_per_t, columns = cond_nut_bound$index) #set maximize flux into nutrient import as the DR*concentration
      
    
    # set condition-specific auxotrophies - overwrite possible auxotrophies with Inf max flux
    auxoIndecies <- c(1:length(Sinfo[,1]))[Sinfo$reaction %in% possibleAuxotrophies]
    set.bounds(lpObj, upper = ifelse(Sinfo[auxoIndecies,]$reaction %in% treatment_par[[treatment]]$auxotrophies, 0, Inf), columns = auxoIndecies)
    
    ## Determine flux distributions and growth rate - jointly maximizing growth and minimizing ##
    
    flux_penalty = 0.001 # resistance on fluxes to minimize feutality
    
    costFxn = c(rep(0, times = length(S[1,]) -1), -1) + c(rep(1, times = length(S[1,])-1), 0)*flux_penalty
      
    set.objfn(lpObj, costFxn)
    solve(lpObj)
    growth_rate$growth[treatment] <- get.objective(lpObj)*-1
    
    objective_flux <- get.variables(lpObj)
    
    collapsedFlux <- sapply(unique(Sinfo$reaction), function(frcombo){
      frindices <- c(1:length(Sinfo[,1]))[Sinfo$reaction == frcombo]
      sum(ifelse(Sinfo$direction[frindices] == "F", 1, -1) * objective_flux[frindices])
    })
    
    flux_vectors[[names(treatment_par)[treatment]]]$"flux" <- collapsedFlux
    
    ## Determine shadow prices ##
    
    costFxn = c(rep(0, times = length(S[1,]) -1), -1)
    set.objfn(lpObj, costFxn)
    solve(lpObj)
    
    flux_vectors[[names(treatment_par)[treatment]]]$"shadowPrices" <- get.dual.solution(lpObj)[2:(length(S[,1])+1)] 
    }  
  
  if(generatePhPP & (PhPPgenerated == FALSE) & is.na(treatment_par[[treatment]]$auxotrophies)){ #this code needs to be fixed - should be straightforward to modify lpObj
    
    for(nutComp in 1:length(PhPPcond$nutrients)){
      nutrientConditions <- PhPPcond$gradient[[nutComp]]
      
      nutIndices <- sapply(colnames(nutrientConditions), function(rxmatch){c(1:length(Sinfo[,1]))[Sinfo$rxDesignation == paste(paste(rxmatch, "boundary"), "F", sep = "_")]})
      
      nutShadow <- matrix(NA, nrow = length(nutrientConditions[,1]), ncol = length(S[,1]))
      nutGR <- rep(NA, times = length(nutrientConditions[,1]))
            
      for(nutCompConc in 1:length(nutrientConditions[,1])){
        set.bounds(lpObj, upper = nutrientConditions[nutCompConc,], columns = nutIndices)
        solve(lpObj)  
        nutShadow[nutCompConc,] <- get.dual.solution(lpObj)[1:length(S[,1])]
        nutGR[nutCompConc] <- get.objective(lpObj)
        }
      
      xlattice <- sort(unique(nutrientConditions[,1]))
      ylattice <- sort(unique(nutrientConditions[,2]))
      zlattice <- matrix(nutGR, conc_samples + 1, conc_samples + 1, byrow = TRUE)[(conc_samples + 1):1,(conc_samples + 1):1]
        
        
      PhPPcond$nutShadow[[nutComp]] <- nutShadow
      PhPPcond$nutGR[[nutComp]] <- nutGR
      
      #determine how to color different regions by shadow price clustering
      #heatmap.2(t(nutShadow), trace = "none")
      
      gradColor <- colorRampPalette(c("aliceblue", "firebrick1"))(1000)
      nutDF <- cbind(nutrientConditions, GR = abs(nutGR), color = gradColor[ceiling(nutDF$GR/max(nutDF$GR)*1000)])
      
      PhPPcond$plotDF[[nutComp]] <- nutDF
      
      #plot3d(x = nutDF[,colnames(nutDF) == PhPPcond$nutrient[[nutComp]][1]], y = nutDF[,colnames(nutDF) == PhPPcond$nutrient[[nutComp]][2]], z = nutDF$GR, col = nutDF$color, xlab = paste(PhPPcond$nutrient[[nutComp]][1], "M/hr"), ylab = paste(PhPPcond$nutrient[[nutComp]][2], "M/hr"), zlab = "Biomass")
      #nutDFtmp <- nutDF; colnames(nutDFtmp)[1:2] <- c("x", "y")
      #wireframe(GR ~ x*y, data = nutDFtmp, xlab = paste(PhPPcond$nutrient[[nutComp]][1], "M/hr"), ylab = paste(PhPPcond$nutrient[[nutComp]][2], "M/hr"), zlab = "Biomass", drape = TRUE)
      }
    PhPPgenerated <- TRUE #only do this part once
    }  
  }    

rxNames <- unique(Sinfo$reaction); rxNames[grep('r_[0-9]+', rxNames)] <- unname(rxnIDtoEnz(rxNames[grep('r_[0-9]+', rxNames)]))
fluxMat <- matrix(NA, ncol = n_c, nrow = length(flux_vectors[[1]]$flux)); colnames(fluxMat) <- names(flux_vectors); rownames(fluxMat) <- rxNames

shadowMat <- matrix(NA, ncol = n_c, nrow = length(flux_vectors[[1]]$shadowPrices)); colnames(shadowMat) <- names(flux_vectors); rownames(shadowMat) <- unname(metIDtoSpec(rownames(S)))

for(i in 1:n_c){
  shadowMat[,i] <- flux_vectors[[i]]$shadowPrices
  fluxMat[,i] <- flux_vectors[[i]]$flux
  }

########### Quadratic programming to match nutrient uptake/excretion rates and produce biomass ####

if(QPorLP == "QP"){

  library(gurobi) #interface for gurobi solver
  #ln -s /Library/gurobi510/mac64/lib/libgurobi51.so libgurobi51.so #a symbolic link was necessary to get gurobi to find its C++ code

  qpModel <- list()
  #qpparams <- list(Presolve=2, OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)
  qpparams <- list(OptimalityTol = 10^-9, FeasibilityTol = 10^-9, BarConvTol = 10^-16)

  
  #flux_elevation_factor <- 10^9 # multiply boundary fluxes by a factor so that minute numbers don't effect tolerance
  #flux_penalty <- 10^-6
  
  #flux_elevation_factor <- 10^9 # multiply boundary fluxes by a factor so that minute numbers don't effect tolerance
  #flux_penalty <- 1/flux_elevation_factor
  
  flux_elevation_factor <- 1000
  flux_penalty <- 1/flux_elevation_factor
 # flux_penalty <- 0.001
  
  qpModel$A <- S
  qpModel$rhs <- Fzero #flux balance
  qpModel$sense <- rep("=", times = length(S[,1])) #global flux balance
  qpModel$lb <- rep(0, times = length(S[1,])) #all fluxes are greater than zero
  qpModel$ub <- rep(Inf , times = length(S[1,])) #overwrite bookkeeping fluxes with empirical maximum rates
  
  qpModel$Q <- diag(rep(0, length(S[1,]))) #min t(v)Qv
  
  qpModel$obj <- rep(flux_penalty, length(S[1,])) #min c * v where v >= 0 to 
  qpModel$obj[grep('^r_', Sinfo$rxDesignation, invert = TRUE)] <- 0
  
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
        
        #if(sum(grep('_R$', biomassSpecRx))){
        #  qpModel$ub[stacked_comp_offset$index[stacked_comp_offset$rxDesignation == biomassSpecRx]] <- 0.2
        #  }else if(sum(grep('_F$', biomassSpecRx))){
        #    qpModel$ub[stacked_comp_offset$index[stacked_comp_offset$rxDesignation == biomassSpecRx]] <- 4
        #    }
        
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
    
    ##### add an additional rxn which introduces a metabolite of interest ####
    
    #rxnFile[grep('NADP\\(\\+\\)$', rxnFile$MetName),][1:10,]
    #ATP: 446
    #ADP: 400
    #pi: 1207
    #NAD+: s_1082
    #NADH: s_1087
    #H+: s_0764_b
    #NADPH: s_1096
    #NADP+: s_1091
    
    
    #balanceStoi <- data.frame(specie = c("s_0446", "s_0400", "s_1207"), stoi = c(1, -1, -1))
    #balanceStoi <- data.frame(specie = c("s_1087", "s_1082", "s_0764_b"), stoi = c(1, -1, -1))
    balanceStoi <- data.frame(specie = c("s_1096", "s_1091", "s_0764_b"), stoi = c(1, -1, -1))
    
    
    
    balanceVec <- rep(0, times = nrow(loose_model$A))
    balanceVec[chmatch(balanceStoi$specie, rownames(loose_model$A))] <- balanceStoi$stoi
    
    loose_model <- qpModel
    loose_model$A <- cbind(loose_model$A, balanceVec)
    loose_model$lb <- c(loose_model$lb, -Inf)
    loose_model$ub <- c(loose_model$ub, Inf)
    loose_model$Q <- diag(c(diag(loose_model$Q), 0))
    loose_model$obj <- c(loose_model$obj, 0)
    
    solvedModel <- gurobi(loose_model, qpparams)
    
    if(solvedModel$status == "NUMERIC"){
      growth_rate$L1[treatment] <- NA
      growth_rate$L2[treatment] <- NA
      next
      }
    ##
    
    #solvedModel <- gurobi(qpModel, qpparams)
    
    #metIDtoSpec(rownames(qpModel$A)[order(solvedModel$pi, decreasing = T)][1:100])
    #solvedModel$pi[order(solvedModel$pi, decreasing = T)][1:10]
    
    
    
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
    
    
    #range(qpModel$A %*% solvedModel$x)
    
    growth_rate$L1[treatment] <- sum(solvedModel$x*loose_model$obj)
    growth_rate$L2[treatment] <- sum(t(solvedModel$x) %*% loose_model$Q %*% t(t(solvedModel$x)))
    
    
    #growth_rate$L1[treatment] <- sum(solvedModel$x*qpModel$obj)
    #growth_rate$L2[treatment] <- sum(t(solvedModel$x) %*% qpModel$Q %*% t(t(solvedModel$x)))
    
    
    flux_vectors[[names(treatment_par)[treatment]]]$"flux" <- collapsedFlux
    flux_vectors[[names(treatment_par)[treatment]]]$"constraints" <- constrainedFlux
    flux_vectors[[names(treatment_par)[treatment]]]$"residual" <- residualFlux
    
    residual_flux_stack <- rbind(residual_flux_stack, cbind(condition = chemostatInfo$condition[treatment], residualFlux))
    
    ### Compare the elemental composition of all boundary constraints ### 
    
    modelMetComp <- read.table("../ChemicalSpeciesQuant/stoiMetsComp.tsv", header = TRUE)
    
    boundary_label <- data.frame(reaction = residualFlux$reactions, color = NA)
    boundary_label$color[grep('boundary', boundary_label$reaction)] <- brewer.pal(length(grep('boundary', boundary_label$reaction)), "Set3")
    boundary_label$color[grep('composition', boundary_label$reaction)] <- brewer.pal(length(grep('composition', boundary_label$reaction)), "Pastel2")
    
    boundary_stoichiometry <- qpModel$A[,sapply(residualFlux$reactions, function(x){(1:nrow(Sinfo))[Sinfo$reaction == x & Sinfo$direction == "F"][1]})]
    boundary_stoichiometry <- boundary_stoichiometry[rowSums(boundary_stoichiometry != 0) != 0,]
    boundary_stoichiometry <- boundary_stoichiometry[grep('bookkeeping', rownames(boundary_stoichiometry), invert = T),] 
    
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
    
    
    #ele_df_melt <- melt(ele_df, id.vars = c("reactions", "condition", "element", "exchange", "color"))
    
    #barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "bottom", 
    #  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.key.width = unit(3, "line"),
    #  strip.background = element_rect(fill = "coral1"))
    
    #ggplot(ele_df, aes(x = factor(1), y = netFlux, fill = color)) + geom_bar(stat = "identity", position = "stack") + barplot_theme + facet_grid(element ~ exchange, scale = "free") + scale_x_discrete("", expand = c(0,0)) + scale_fill_identity(name = "Class", guide = guide_legend(nrow = 4), labels = boundary_label$reaction, breaks = boundary_label$color)
    #ggplot(ele_df, aes(x = exchange, y = experimentalFlux, fill = color)) + geom_bar(stat = "identity", position = "stack") + barplot_theme + facet_grid(element ~ exchange, scale = "free") + scale_x_discrete("", expand = c(0,0)) + scale_fill_identity(name = "Class", guide = guide_legend(nrow = 4), labels = boundary_label$reaction, breaks = boundary_label$color)
    
    #ggplot(ele_df_melt, aes(x = exchange, y = value, fill = color)) + geom_bar(stat = "identity", position = "stack") + barplot_theme + facet_grid(element ~ variable + exchange, scale = "free") + scale_x_discrete("", expand = c(0,0)) + scale_fill_identity(name = "Class", guide = guide_legend(nrow = 5), labels = boundary_label$reaction, breaks = boundary_label$color)
    
    
    
    }  
  }  
reaction_info('r_0005')
rxnFile[rxnFile$ReactionID == "r_0005",]
trackedMet = '^NAD\\(\\+\\)'

trackMetConversion('^acetate$')
trackMetConversion('ergosta-5,7,24\\(28\\)-trien-3beta-ol')
trackMetConversion('^NAD\\(\\+\\)', T)
trackMetConversion('ATP$', T)
trackMetConversion('ubiquinol-6', T)
trackMetConversion('cytochrome', T)

rxn_search(named_stoi, 'cytochrome')

trackedMet <- 'ferricytochrome'







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


ggplot(residual_flux_stack[!is.na(residual_flux_stack$sd),], aes(x = factor(condition), y = resid_st*-1)) + facet_wrap(~ reactions, ncol = 2, scale = "free_y") + geom_bar(stat = "identity") + barplot_theme + scale_y_continuous("predicted - experimental flux / sd(experimental") + ggtitle("Fit of experimental and optimized flux")
ggsave("flux_residuals.pdf", width = 8, height = 20)
  
residual_flux_melted <- dcast(residual_flux_stack, formula = reactions ~ condition, value.var = "resid_st")
  residual_flux_stack$resid_st

residual_flux_stack[residual_flux_stack[,2] == "AA flux composition",]

### for now only look at reactions which carry flux in >= 80% of conditions
#rxNames <- rxNames[rowSums(fluxMat_per_cellVol == 0) <= 5,]
#fluxMat_per_cellVol <- fluxMat_per_cellVol[rowSums(fluxMat_per_cellVol == 0) <= 5,]

heatmap.2(fluxMat_per_cellVol / (t(t(apply(fluxMat_per_cellVol, 1, sd))) %*% rep(1, n_c)), trace = "none", Colv = F)


flux_summary <- list()
flux_summary$IDs = rxNames
flux_summary$ceulluarFluxes = fluxMat_per_cellVol


save(flux_summary, file = "fluxSummaryQP.Rdata")


heatmap.2(fluxMat_per_cellVol / (t(t(apply(abs(fluxMat_per_cellVol[rowSums(fluxMat_per_cellVol) != 0,]), 1, mean))) %*% rep(1, n_c)), trace = "none")


high_flux <- order(rowSums(abs(fluxMat_per_cellVol)), decreasing = T)
rownames(fluxMat_per_cellVol)[high_flux][1:50]


row_check <- 100
for(row_check in 1:100){
if(sum(fluxMat_per_cellVol[row_check,] <= 0) == 0){
  plot_bounds <- c(0, max(fluxMat_per_cellVol[row_check,]))
  } else if(sum(fluxMat_per_cellVol[row_check,] <= 0) == n_c){
    plot_bounds <- c(min(fluxMat_per_cellVol[row_check,]), 0)
} else{plot_bounds <- c(min(fluxMat_per_cellVol[row_check,]), max(fluxMat_per_cellVol[row_check,]))}
print(plot(fluxMat_per_cellVol[row_check,] ~ chemostatInfo$actualDR[1:n_c], pch = 16, cex = 2, col = factor(chemostatInfo$limitation[1:n_c]), xlim = c(0, 0.3), ylim = plot_bounds, xlab = "DR", ylab = "flux per cell volume (moles/mL cell volume per hr)"))
}

  







  






###### output fluxes so that they can be visualzied using S. cerevisae cellular overview #####

choice_conditions <- names(flux_vectors)[grep(0.05, names(flux_vectors))]
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


	
	