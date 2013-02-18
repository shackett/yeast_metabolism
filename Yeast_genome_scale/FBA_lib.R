
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

perfect.match <- function(source, query, corrFile){
all_char <- "[[:graph:][:space:]]"
	
tmp <- corrFile[grep(source, query)[!(grep(source, query) %in% union(grep(paste(all_char, source, sep = ""), query), grep(paste(source, all_char, sep = ""), query)))],]
if(length(tmp[,1]) == 0){tmp <- corrFile[grep(source, query, fixed = TRUE),]}
tmp
	}


rxn_search = function(stoiMat, search_string, is_rxn = TRUE, index = FALSE){
	#search by metabolite or reactant and return all reactions and nonzero metabolites.
  #stoiMat rows and columns must be named with metabolite and enzyme common names: named_stoi
	if (is_rxn == TRUE){
		colz = grep(search_string, colnames(stoiMat), fixed = TRUE)
		} else {
		met = grep(search_string, rownames(stoiMat), fixed = TRUE)
		if (length(met) == 1){
			colz = c(1:length(stoiMat[1,]))[stoiMat[met,] != 0]
			} else {
			colz = c(1:length(stoiMat[1,]))[apply(stoiMat[met,], 2, is.not.zero)]
		}}
	
	if(length(colz) == 0){
		print("no hits")
		} else {
			if(index == TRUE){
				colz
				} else {
			
			rxns = stoiMat[,colz]
			if(is.vector(rxns)){
				c(colz, rxns[rxns != 0])
				} else {
					output <- rbind(colz, rxns[apply(rxns, 1, is.not.zero),])
					colnames(output) = colnames(stoiMat)[colz]
					output
					}}
		}
	}

flip.rxn <- function(reactions, joint.stoi){
	stoi <- joint.stoi
	stoi[,colnames(joint.stoi) %in% reactions] <- stoi[,colnames(joint.stoi) %in% reactions]*-1
	stoi
	}
	
is.not.zero = function(vec){
	length(vec[vec!=0]) != 0
	}	
	
######### fxns to convert between IDs and species ######

metIDtoSpec <- function(meta){
	sapply(meta, function(x){
		corrFile$SpeciesName[corrFile$SpeciesID == x]
		})}

rxnIDtoEnz <- function(rxn){
	sapply(rxn, function(x){
		rxnFile$Reaction[rxnFile$ReactionID == x][1]
		})}
		
rxnIDtoGene <- function(rxns){
  sapply(rxns, function(rx){
    paste(unique(strsplit(paste(rxnFile[rxnFile$ReactionID == rx,]$MetName[is.na(rxnFile[rxnFile$ReactionID == rx,]$StoiCoef)], collapse = ':'), ':')[[1]]), collapse = '/')
    })}

metToCHEBI <- function(mets){
	#associate species IDs and CHEBI ids where available/applicable
	if(length(grep("chebi", rxnparFile[,3][rxnparFile[,1] == corrFile$SpeciesType[corrFile$SpeciesID %in% mets]])) == 0){
		NA
		}else{
	unlist(strsplit(rxnparFile[,3][rxnparFile[,1] == corrFile$SpeciesType[corrFile$SpeciesID %in% mets]], split = "%3A"))[2]	
	}}

#rxnIDtoSGD <- function(rxnIDs){
  #output the compartment where a reaction occurs followed by all of the genes involved in the rxn
  
 # output <- t(sapply(rxnIDs, function(rxnID){
  #  tmp <- rxnFile[rxnFile$ReactionID == rxnID,]
  #  c(tmp$Compartment[1], paste(unique(strsplit(paste(tmp$MetName[is.na(tmp$StoiCoef)], collapse = ':'), ':')[[1]]), collapse = ':'))
  #}))
  



eval_mets <- function(query_met, grep_it = FALSE){
	#find all of the reactions that a greped or exact matched metabolite participates in and then get all of the other metabolites also in those reactions
	if(grep_it == TRUE){
		met_matches <- grep(query_met, metIDtoSpec(rownames(stoiMat)))
		}else{
			met_matches <- c(1:length(stoiMat[,1]))[metIDtoSpec(rownames(stoiMat)) %in% query_met]
			}
	if(length(met_matches) == 0){print("miss")}
	if(length(met_matches) == 1){
		eval_mat <- stoiMat[apply(stoiMat[,stoiMat[met_matches,] != 0] != 0, 1, sum) != 0,stoiMat[met_matches,] != 0]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))
		}else{
		eval_mat <- stoiMat[apply(stoiMat[,apply(stoiMat[met_matches,] != 0, 2, sum) != 0] != 0, 1, sum) != 0, apply(stoiMat[met_matches,] != 0, 2, sum) != 0]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))
		}
	eval_mat
	}		
  
  
  
reaction_info <- function(rxnName){
    
  ##### write a function to list:
  # reactants -> products
  # Reaction name and designation
  # Thermodynamics
  # KEGG and EC reactions name
    
  rxnStoi <- stoiMat[,colnames(stoiMat) == rxnName][stoiMat[,colnames(stoiMat) == rxnName] != 0]
  speciesNames <- metIDtoSpec(names(rxnStoi))
  rxnDir <- reversibleRx$reversible[reversibleRx$rx == rxnName]
  if(rxnDir == 1){rxnDir <- " -> "}
  if(rxnDir == 0){rxnDir <- " <=> "}
  if(rxnDir == -1){rxnDir <- " <- "}
    
  substrate_prep <- paste(sapply(c(1:length(rxnStoi[rxnStoi < 0])), function(x){
    tmp <- (rxnStoi[rxnStoi < 0] * -1)[x]
    if(tmp == 1){tmp <- ''}
    paste(tmp, speciesNames[rxnStoi < 0][x])
  }), collapse = ' + ')
    
  product_prep <- paste(sapply(c(1:length(rxnStoi[rxnStoi > 0])), function(x){
    tmp <- (rxnStoi[rxnStoi > 0])[x]
    if(tmp == 1){tmp <- ''}
    paste(tmp, speciesNames[rxnStoi > 0][x])
  }), collapse = ' + ')
    
  rxList <- list()
  rxList$reaction = unname(rxnIDtoEnz(rxnName))
  rxList$enzymes = unname(rxnIDtoGene(rxnName))
  rxList$stoichiometry = paste(substrate_prep, rxnDir, product_prep)
  rxList$thermo = reversibleRx[reversibleRx$rx == rxnName,]
  rxList  
  }  
  