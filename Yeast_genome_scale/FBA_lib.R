
write_stoiMat = function(metabolites, reactions, corrFile, rxnFile, internal_names = FALSE){

#create the stoichiometry matrix, indicating the change of metabolites (rows) per chemical reaction (columns)

	stoiMat <- matrix(data = 0, ncol = length(reactions), nrow = length(metabolites))

	metName = rep(NA, times = length(metabolites))
	for (i in 1:length(metabolites)){
		metName[i] = as.character(corrFile$SpeciesName[corrFile$SpeciesID %in% metabolites[i]][1])
		}
	rxnName = rep(NA, times = length(reactions))
	for (i in 1:length(reactions)){
		rxnName[i] = as.character(rxnFile$Reaction[rxnFile$ReactionID %in% reactions[i]][1])
		}
	rownames(stoiMat) <- metName	
	colnames(stoiMat) <- rxnName

	for (i in 1:length(rxnFile[,1])){
		stoiMat[c(1:length(metabolites))[metabolites == rxnStoi$Metabolite[i]], c(1:length(reactions))[reactions == rxnStoi$ReactionID[i]]] <- rxnStoi$StoiCoef[i]		}

	if (internal_names == FALSE){
		save(stoiMat, file = "yeast_stoi.R")} else {
		rownames(stoiMat) <- metabolites	
		colnames(stoiMat) <- reactions
		save(stoiMat, file = "yeast_stoi.R")
			}
	}

