library(lpSolve)
library(limSolve)

setwd("/Users/seanhackett/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")
source("FBA_lib.R")

#inputFilebase = "yeast_GS"
inputFilebase = "yeast"


rxnFile = read.delim(paste("rxn_", inputFilebase, ".tsv", sep = ""))
#rxnparFile = read.delim(paste("species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE)
corrFile = read.delim(paste("spec_", inputFilebase, ".tsv", sep = ""))

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)

#load or write stoichiometry matrix of reactions and their altered metabolites

if(file.exists("yeast_stoi.R")){
	load("yeast_stoi.R")
	} else {write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}






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


