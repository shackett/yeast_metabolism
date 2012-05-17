library(lpSolve)

rxnFile = read.table("/Users/seanhackett/Desktop/RabinowitzLab/FBA_SRH/Yeast_reconstruction/rxn_yeast_nocomp.tsv", sep = "\t", colClasses = c("factor", "factor", "factor", "numeric"), header = TRUE, fill = TRUE)
parFile = read.delim("/Users/seanhackett/Desktop/RabinowitzLab/FBA_SRH/Yeast_reconstruction/species_par_yeast_nocomp.tsv")
corrFile = read.delim("/Users/seanhackett/Desktop/RabinowitzLab/FBA_SRH/Yeast_reconstruction/spec_yeast_nocomp.tsv")

reactions = unique(rxnFile$ReactionID)

rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)

#create the stoichiometry matrix, indicating the change of metabolites (rows) per chemical reaction (columns)

stoiMat <- matrix(data = 0, ncol = length(reactions), nrow = length(metabolites))
rownames(stoiMat) <- metabolites
colnames(stoiMat) <- reactions

for (i in 1:length(rxnFile[,1])){
	stoiMat[c(1:length(metabolites))[metabolites == rxnStoi$Metabolite[i]], c(1:length(reactions))[reactions == rxnStoi$ReactionID[i]]] <- rxnStoi$StoiCoef[i]	
	}

#specify the lower and upper bounds of each reaciton flux

boundaries <- matrix(data = NA, ncol = 2, nrow = length(reactions))
rownames(boundaries) <- reactions
colnames(boundaries) <- c("lowerBound", "upperBound")

for (i in 1:length(reactions)){
	rxnPars <- parFile[parFile$ReactionID %in% reactions[i],]
	if(length(rxnPars[,1]) != 4){
		print(i)
		}
		
	boundaries[i,1] <- rxnPars$Value[rxnPars$ParID == "LOWER_BOUND"][1]
	boundaries[i,2] <- rxnPars$Value[rxnPars$ParID == "UPPER_BOUND"][1]
	}
	
optimization_vec = stoiMat[,colnames(stoiMat) == "R_biomass"]
const_mat <- stoiMat[,!(colnames(stoiMat) %in% "R_biomass")]
	
	
lp(direction = "min", const.dir = "<="








