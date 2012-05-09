library(lpSolve)
library(limSolve)

setwd("/Users/seanhackett/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")
source("FBA_lib.R")

#inputFilebase = "yeast_GS"
inputFilebase = "yeast"


rxnFile = read.delim(paste("rxn_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
rxnparFile = read.delim(paste("species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
corrFile = read.delim(paste("spec_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
compFile <- read.delim(paste("comp_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)

compositionFile <- read.csv2("../Yeast_comp.csv", sep = ",", stringsAsFactors = FALSE)
nutrientFile <- read.delim("Boer_nutrients.txt")[1:6,1:6]
rownames(nutrientFile) <- nutrientFile[,1]; nutrientFile <- nutrientFile[,-1]

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)

######### fxns to convert between IDs and species ######

metIDtoSpec <- function(meta){
	sapply(meta, function(x){
		corrFile$SpeciesName[corrFile$SpeciesID == x]
		})}

rxnIDtoEnz <- function(rxn){
	sapply(rxn, function(x){
		rxnFile$Reaction[rxnFile$ReactionID == x][1]
		})}












######### treatment ###########

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

#load or write stoichiometry matrix of reactions and their altered metabolites

if(file.exists("yeast_stoi.R")){
	load("yeast_stoi.R")
	} else {write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}

############ preserving compartmentation ############

#reactions = unique(rxnFile$ReactionID)

#reactions are compartment specific

compartment <- sapply(reactions, function(x){rxnFile$Compartment[rxnFile$ReactionID == x][1]})

#extract the metabolite ID corresponding to the extracellular introduction of nutrients

sources <- c("D-glucose", "ammonium", "phosphate", "sulphate", "uracil", "L-leucine")

perfect.match <- function(source, query, corrFile){
all_char <- "[[:graph:][:space:]]"
	
tmp <- corrFile[grep(source, query)[!(grep(source, query) %in% union(grep(paste(all_char, source, sep = ""), query), grep(paste(source, all_char, sep = ""), query)))],]
if(length(tmp[,1]) == 0){tmp <- corrFile[grep(source, query, fixed = TRUE),]}
tmp
	}
	
	

resource_matches <- lapply(sources, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

boundary_met <- NULL
for(x in 1:length(sources)){
boundary_met <- rbind(boundary_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}

#extract the metabolite ID corresponding to cytosolic metabolites being assimilated into biomass

sinks <- compositionFile$AltName

resource_matches <- lapply(sinks, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

comp_met <- NULL
for(x in 1:length(sinks)){
comp_met <- rbind(comp_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "cytoplasm"],])
}

#freely exchanging metabolites through extracellular compartment

free_flux <- c("carbon dioxide", "oxygen", "water")

resource_matches <- lapply(free_flux, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

freeExchange_met <- NULL
for(x in 1:length(free_flux)){
freeExchange_met <- rbind(freeExchange_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}


#### search for un-balanced rxns ######

is.unbalanced <- rep(NA, times = length(stoiMat[1,]))

for(i in 1:length(stoiMat[1,])){
	is.unbalanced[i] <- ifelse((length(stoiMat[,i][stoiMat[,i] > 0]) != 0) & (length(stoiMat[,i][stoiMat[,i] < 0]) != 0), FALSE, TRUE)
	}

rem.unbalanced <- colnames(stoiMat)[is.unbalanced]


boundary_put <- stoiMat[,is.unbalanced][apply(abs(stoiMat[,is.unbalanced]), 1, sum) != 0,]

rownames(boundary_put) <- metIDtoSpec(rownames(boundary_put))
colnames(boundary_put) <- rxnIDtoEnz(colnames(boundary_put))

all_rxns <- rxnIDtoEnz(colnames(stoiMat))
#all_rxns[grep("biomass", all_rxns)]

stoiMat[colnames(stoiMat) %in% "r_1812"][stoiMat[colnames(stoiMat) %in% "r_1812"] != 0]
metIDtoSpec(rownames(stoiMat)[stoiMat[colnames(stoiMat) %in% "r_1812"] != 0])












growth_rate <- data.frame(limit = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[1]}), dr = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[2]}), growth = NA)

flux_vectors <- list()
######################## Set up the linear equations for FBA #######################

for(treatment in 1:length(names(treatment_par))){

#remove reactions which are defective (such 
S_rxns = stoiMat[,!(colnames(stoiMat) %in% c(treatment_par[[treatment]]$auxotrophies, rem.unbalanced))]

#added reactions for boundary fluxes

##### Nutrient Influx #####
##### Unconstrained Chemical Influx #####

influx_rxns <- c(1:length(metabolites))[metabolites %in% c(boundary_met$SpeciesID, freeExchange_met$SpeciesID)]

influxS <- matrix(0, ncol = length(influx_rxns), nrow = length(metabolites))
for(i in 1:length(influx_rxns)){
	influxS[influx_rxns[i],i] <- -1
	}
	
##### Composition fxn #####

compVec <- rep(0, times = length(metabolites))
for(i in 1:length(comp_met$SpeciesID)){
	compVec[rownames(stoiMat) == comp_met$SpeciesID[i]] <- as.numeric(compositionFile$StoiCoef)[i]
	}
	
S <- cbind(S_rxns, influxS, compVec)		
colnames(S) <- c(colnames(S_rxns), sapply(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName), function(x){paste(x, "boundary")}), "composition")

################ F - flux balance ############

Fzero <- rep(0, times = length(S[,1]))


#dim nconst x nrxns
#the reactions corresponding to the boundary fluxes added in influx reactions for nutrients are given a value of -1.

influxG <- matrix(0, ncol = length(S[1,]), nrow = length(boundary_met$SpeciesID))
influxh <- NULL

for(i in 1:length(boundary_met$SpeciesName)){
influxG[i, c(1:length(S[1,]))[colnames(S) %in% paste(as.character(boundary_met$SpeciesName)[i], "boundary")]] <- -1
influxh	<- c(influxh, -1*treatment_par[[treatment]]$nutrients$conc_per_t[treatment_par[[treatment]]$nutrients$nutrient %in% boundary_met$SpeciesName[i]])
	
	}

############### costFxn - indicates the final rxn in S ######

costFxn = c(rep(0, times = length(S[1,]) -1), 1)

######## use linear programming to maximize biomass #######

linp_solution <- linp(E = S, F = Fzero, G = influxG, H = influxh, Cost = costFxn, ispos = FALSE)


flux_vectors[[names(treatment_par)[treatment]]] <- linp_solution$X

growth_rate$growth[treatment] <- linp_solution$solutionNorm

}


v <- linp(E = S, F = Fzero, G = influxG, H = influxh, Cost = costFxn, ispos = FALSE)$X
v[v!=0]


v <- rep(1, times = length(S[1,]))

influxG%*%v >= influxh
v%*%costFxn


influxH <- rep(0, times = length(boundary_met$SpeciesID))

c(1:length(metabolites))[metabolites %in% as.character(boundary_met$SpeciesID)]



#Gv >= h






}











######### approximate the concentrations of source nutrients as a function of the rate of compounds accumulating and their composition.

accumulants <- data.frame(ID = rownames(stoiMat)[stoiMat[,colnames(stoiMat) == "R_biomass"] != 0], S_name = NA, ChEBI_name = NA, common_name = NA, values = stoiMat[,colnames(stoiMat) == "R_biomass"][stoiMat[,colnames(stoiMat) == "R_biomass"] != 0])

for(i in 1:length(accumulants[,1])){
	accumulants$S_name[i] = as.character(corrFile$SpeciesType[corrFile$SpeciesID %in% accumulants$ID[i]])
	}

parFile_chebi = rxnparFile[grep("chebi", rxnparFile[,3]),]

for (i in 1:length(parFile_chebi[,1])){
	parFile_chebi[i,4] <- unlist(strsplit(as.character(parFile_chebi[i,3]), split = "CHEBI:"))[2]
	}

parFile_chebi <- parFile_chebi[,-3]
parFile_chebi[,1] <- as.character(parFile_chebi[,1])

for(i in 1:length(accumulants[,1])){
	accumulants$ChEBI_name[i] <- parFile_chebi[,3][parFile_chebi[,1] == accumulants$S_name[i]]
	accumulants$common_name[i] <- as.character(parFile_chebi[,2][parFile_chebi[,1] == accumulants$S_name[i]])
	}


chebiComp = read.delim("METeleComp.tsv")
chebiComp = chebiComp[chebiComp$ID %in% parFile_chebi[,3],]

chebiComp <- cbind(chebiComp[,c(1,2)], chebiComp[3:length(chebiComp[1,])][,apply(chebiComp[3:length(chebiComp[1,])], 2, sum) != 0])

elements = colnames(chebiComp)[3:length(chebiComp[1,])]
accumulant_comp = matrix(data = 0, ncol = length(elements), nrow = length(accumulants[,1]))
accumulant_comp = as.data.frame(accumulant_comp)
colnames(accumulant_comp) = elements
rownames(accumulant_comp) = accumulants$ChEBI_name


for (i in 1:length(accumulant_comp[,1])){
	
	if(accumulants$ChEBI_name[i] %in% chebiComp$ID){
		accumulant_comp[i,] = chebiComp[3:length(chebiComp[1,])][chebiComp$ID %in% accumulants$ChEBI_name[i],]
		} else {
		accumulant_comp[i,] = rep(NA, times = length(elements))		}
	}

#glucan, mannan, glycogen - C6H10O5
accumulant_comp[c("37671", "28808", "28087"),] = rep(c(6, 0, 10, 0, 0, 0, 5, 0, 0), each = 3)

elemental_accum = apply(accumulants$value * accumulant_comp, 2, sum)



source_comp = matrix(data = 0, ncol = length(elements), nrow = 11)
colnames(source_comp) = elements
rownames(source_comp) = c("C6H12O6", "(NH4)2SO4", "KH2PO4", "CaCl2", "NaCl", "MgSO4", "KCl", "water", "CO2", "O2", "H+")

source_comp[1,c(1,3,7)] = c(6,12,6)
source_comp[2,c(3,5,7,9)] = c(8, 2, 4, 1)
source_comp[3,c(3,4,7,8)] = c(2,1,4,1)	
source_comp[5,6] = 1
source_comp[6,c(7,9)] = c(4,1)
source_comp[7,4] = 1
source_comp[8,c(3,7)] = c(2,1)
source_comp[9,c(1,7)] = c(1,2)
source_comp[10,7] = 2	
source_comp[11,3] = 1

tsource_comp = t(source_comp)

#H = elemental_accum
H = c(-1*elemental_accum, rep(0, 11))

#G = source_comp
G = rbind(tsource_comp, diag(11)*c(rep(1, times = 8), -1, 1, 1))
G[,9] <- G[,9]*-1

Cost = apply(source_comp, 1, sum)
Cost[9] <- -3


solver = linp(G=G, H=H, Cost = Cost, ispos = FALSE)

ex = solver$X

G %*% ex >= 0




#solver$X %*% source_comp


