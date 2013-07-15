
########### Function used in flux prediction and rxn species conversions ###########

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

metIDtoSpec <- function(meta, includeComp = F){
	if(includeComp){
    sapply(meta, function(x){
       paste(corrFile[corrFile$SpeciesID == x, colnames(corrFile) %in% c("SpeciesName", "Compartment")] , collapse = "-")
       })
    }else{
     sapply(meta, function(x){
       corrFile$SpeciesName[corrFile$SpeciesID == x]
       })
    }
}

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

rxnIDtoSGD <- function(rxnIDs){
  #output the compartment where a reaction occurs followed by all of the genes involved in the rxn  
  
  output <- t(sapply(rxnIDs, function(rxnID){
    tmp <- rxnFile[rxnFile$ReactionID == rxnID,]
    c(tmp$Compartment[1], paste(unique(strsplit(paste(tmp$MetName[is.na(tmp$StoiCoef)], collapse = ':'), ':')[[1]]), collapse = ':'))
  }))
output
}

write.output <- function(tab, output){
  #write an output table in a way where all columns will have a column name
  #rownames go to 1st column, 1st column name is ``Gene''
  
  tab <- data.frame(Gene = rownames(tab), tab, stringsAsFactors = FALSE)
  write.table(tab, file = output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


convert_to_elemental <- function(redStoi, modelMetComp){
  # convert a stoichiometric matrix (with named columns) and "S" IDs as rownames to their elemental constituents
  require(data.table)
  
  redStoi <- boundary_stoichiometry
  specieEle <- modelMetComp[chmatch(rownames(redStoi), modelMetComp$ID),]
  
  specieEle[specieEle$name == "(1->3)-beta-D-glucan", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 10, 5) #one oxygen shared bc of condensation
  specieEle[specieEle$name == "glycogen", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 10, 5)
  specieEle[specieEle$name == "mannan", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 10, 5)
  specieEle[specieEle$name == "D-glucose", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 12, 6)
  specieEle[specieEle$name == "polyphosphate", colnames(specieEle) %in% c("O", "P")] <- c(4, 1)
  specieEle[is.na(specieEle)] <- 0
  specieEle[,-c(1:2)] <- apply(specieEle[,-c(1:2)], c(1,2), as.numeric)
  specieEle <- specieEle[,c(TRUE, TRUE, colSums(specieEle[,-c(1:2)]) != 0)]
  
  atomicMasses <- data.frame(element = c("C", "H", "N", "O", "P", "S"), mass = c(12.0107, 1.00794, 14.00674, 15.9994, 30.973761, 32.066))
  
  if(all(atomicMasses$element == colnames(specieEle[,-c(1,2)]))){
    specieEle$MW <- t(t(specieEle[,-c(1,2)])) %*% t(t(c(atomicMasses[,2])))
    specieEle
    }else{
      print("atomicMasses are mismatched, additional elemental masses need to be included")
      }
  
  }



trackMetConversion <- function(trackedMet, allRxns = FALSE){
  
  ## for a metabolite of interest, seperate reactions which carry flux consuming or producing it into each compartment and display a summary ##
  if(allRxns){
    flux <- collapsedFlux
    }else{
    flux <- collapsedFlux[collapsedFlux != 0]
    }
  
  reducedStoi <- matrix(0, ncol = length(flux), nrow = nrow(stoiMat))
  rownames(reducedStoi) <- rownames(stoiMat)
  colnames(reducedStoi) <- names(flux)
  
  # add model stoichiometry
  index_match <- chmatch(names(flux), colnames(stoiMat))
  reducedStoi[,!is.na(index_match)] <- stoiMat[,index_match[!is.na(index_match)]]
  
  boundaryStoi <- qpModel$A[,sapply(names(flux)[is.na(index_match)], function(x){(1:nrow(Sinfo))[Sinfo$reaction == x & Sinfo$direction == "F"][1]})]
  boundaryStoi <- boundaryStoi[(rownames(boundaryStoi) %in% rownames(stoiMat)),] #remove bookkeeping flux
  reducedStoi[,is.na(index_match)] <- boundaryStoi
  
  metNames <- metIDtoSpec(rownames(reducedStoi))
  
  print(paste("Species Matched :", paste(unique(unname(metNames[grep(trackedMet, metNames)])), collapse = "/")))
  
  if(length(grep(trackedMet, metNames)) == 1){
    involvedRxns <- reducedStoi[,reducedStoi[grep(trackedMet, metNames),] != 0]
    }else{
      involvedRxns <- reducedStoi[,colSums(reducedStoi[grep(trackedMet, metNames),] != 0) != 0]
      }
  
  
  reactionNames <- data.frame(ID = colnames(involvedRxns), reaction = rxnIDtoEnz(colnames(involvedRxns)), compartment = unname(rxnIDtoSGD(colnames(involvedRxns))[,1]))
  reactionNames$compartment[is.na(reactionNames$compartment)] <- "boundary"
  
  # balance within individual compartments
  compartment_balance <- as.data.frame(matrix(0, ncol = length(grep(trackedMet, metNames)), nrow = length(unique(reactionNames$compartment))))
  rownames(compartment_balance) <- unique(reactionNames$compartment)
  colnames(compartment_balance) <- mapply(function(x,y){paste(x, y)}, x = metNames[grep(trackedMet, metNames)], y = compFile$compName[chmatch(rxnFile$Compartment[rxnFile$Compartment != "exchange"][chmatch(rownames(reducedStoi)[grep(trackedMet, metNames)], rxnFile$Metabolite[rxnFile$Compartment != "exchange"])], compFile$compID)])
  
  #visualize each compartment seperately
  for(a_compartment in unique(reactionNames$compartment)){
    comp_name <- ifelse(sum(grep('^c_', a_compartment)) == 1, compFile$compName[compFile$compID == a_compartment], a_compartment)
    as.matrix(reducedStoi[grep(trackedMet, metNames),])
    if(length(grep(trackedMet, metNames)) == 1){
      comp_fluxes <- cbind(reactionNames[reactionNames$compartment == a_compartment,], flux = flux[reducedStoi[grep(trackedMet, metNames),] != 0][reactionNames$compartment == a_compartment], stoi = NA, effect = NA)
      }else{
        comp_fluxes <- cbind(reactionNames[reactionNames$compartment == a_compartment,], flux = flux[colSums(reducedStoi[grep(trackedMet, metNames),] != 0) != 0][reactionNames$compartment == a_compartment], stoi = NA, effect = NA)
        }
    
    comp_stoi <- involvedRxns[,reactionNames$compartment == a_compartment]
    
    if(a_compartment == "boundary" | nrow(comp_fluxes) == 1){
      print(toupper(a_compartment))
      print(comp_fluxes)
      if(nrow(comp_fluxes) == 1){
        compartment_balance[rownames(compartment_balance) == a_compartment,] <- (comp_stoi * comp_fluxes$flux)[grep(trackedMet, metNames)]
        }else{
          compartment_balance[rownames(compartment_balance) == a_compartment,] <- (comp_stoi %*% comp_fluxes$flux)[grep(trackedMet, metNames)]
          }
      next  
    }
    
    compartment_balance[rownames(compartment_balance) == a_compartment,] <- (comp_stoi %*% comp_fluxes$flux)[grep(trackedMet, metNames)]
    
    for(a_rxn in 1:nrow(comp_fluxes)){
      
      rxnStoi <- comp_stoi[comp_stoi[,a_rxn] != 0,a_rxn]*ifelse(comp_fluxes$flux[a_rxn] >= 0, 1, -1)
      if(a_compartment == "exchange"){
        speciesNames <- metIDtoSpec(names(rxnStoi), T)
        }else{
          speciesNames <- metIDtoSpec(names(rxnStoi))
          }
      
      rxnDir <- reversibleRx$reversible[reversibleRx$rx == comp_fluxes$ID[a_rxn]]
      
      metOfInterest <- unname(rxnStoi)[grep(trackedMet, speciesNames)]
      
      if(all(metOfInterest < 0)){
        comp_fluxes$effect[a_rxn] <- "consumed"
        }else if(all(metOfInterest > 0)){
          comp_fluxes$effect[a_rxn] <- "produced"
          }else{
            comp_fluxes$effect[a_rxn] <- "exchanged"
            }
      
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
      
      comp_fluxes$stoi[a_rxn] = paste(substrate_prep, rxnDir, product_prep)
    }
    
    comp_fluxes <- comp_fluxes[order(comp_fluxes$effect, abs(comp_fluxes$flux), decreasing = T),]
    comp_fluxes$flux <- abs(comp_fluxes$flux)
    
    print(toupper(comp_name))
    writeLines('\n')
    
    for(a_change in c("produced", "consumed", "exchanged")){
      print(comp_fluxes[comp_fluxes$effect == a_change,])
      }
    writeLines('\n\n\n')
  }
  print(compartment_balance)
}

######## Diagnostic flux prediction methods ##########

maxFlux <- function(){
  ## Using the current QP model, for each standard reaction (those in the initial model) determine the maximum flux that can be carried through each reaction.
  ## output the maximum flux for reactions with a finite maximum flux
  
  boundModel <- qpModel
  
  SinfoBounds <- Sinfo[boundModel$obj != 0,]; SinfoBounds$index = (1:length(boundModel$obj))[boundModel$obj != 0]; SinfoBounds$maxFlux = NA
  
  for(a_rxn in SinfoBounds$index){
    
    #for each constrained reaction switch from penalizing flux to giving a bonus
    #overwrite
    boundModel$obj[a_rxn] <- 10
    boundModel$A[,a_rxn] <- qpModel$A[,a_rxn]*-1 #reverse stoichiometry and allow for negative flux
    boundModel$lb[a_rxn] <- -Inf
    boundModel$ub[a_rxn] <- 0
    
    #evaluate
    solved_bound <- gurobi(boundModel, list(OutputFlag = 0))
    
    if(solved_bound$status == "NUMERIC"){
      #rerun with simplex
      solved_bound <- gurobi(boundModel, list(method = 1, OutputFlag = 0))
      SinfoBounds$maxFlux[a_rxn] <- ifelse(solved_bound$status == "UNBOUNDED", Inf, unname(solved_bound$x[a_rxn]*-1))
      }else{
        SinfoBounds$maxFlux[a_rxn] <- ifelse(solved_bound$status == "INF_OR_UNBD", Inf, unname(solved_bound$x[a_rxn]*-1))
        }
    #print(solved_bound$status)
    
    #restore
    boundModel$obj[a_rxn] <- qpModel$obj[a_rxn]
    boundModel$A[,a_rxn] <- qpModel$A[,a_rxn]
    boundModel$lb[a_rxn] <- qpModel$lb[a_rxn]
    boundModel$ub[a_rxn] <- qpModel$ub[a_rxn]
  }
  
  flux_bounds <- SinfoBounds[!is.na(SinfoBounds$maxFlux) & SinfoBounds$maxFlux != Inf,]
  flux_bounds
}


loosenFlux <- function(balanceStoi, justObj = FALSE){
  # allow for unbounded flux of an additional reaction to evaluate effects on flux
  
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
  #balanceStoi <- data.frame(specie = c("s_1096", "s_1091", "s_0764_b"), stoi = c(1, -1, -1))
  
  
  loose_model <- qpModel
  
  balanceVec <- rep(0, times = nrow(loose_model$A))
  balanceVec[chmatch(balanceStoi$specie, rownames(loose_model$A))] <- balanceStoi$stoi
  
  loose_model$A <- cbind(loose_model$A, balanceVec)
  loose_model$lb <- c(loose_model$lb, -Inf)
  loose_model$ub <- c(loose_model$ub, Inf)
  loose_model$Q <- diag(c(diag(loose_model$Q), 0))
  loose_model$obj <- c(loose_model$obj, 0)
  
  solvedModel <- gurobi(loose_model, list(OutputFlag = 0))
  
   if(solvedModel$status == "NUMERIC"){
      #rerun with simplex
      solvedModel <- gurobi(loose_model, list(method = 1, OutputFlag = 0))
      }
  if(justObj){
    solvedModel$objval
    }else{
      solvedModel
      }
  
}


forcedFlux <- function(forcedRx){
  # force the flux through reactions (given by reaction designation) to exceed a provided value
  require(data.table)
  
  forced_model <- qpModel
  
  # bound a flux by the forced value and if reversible bound the complementary flux to prevent compensatory reverse flux
  forced_model$lb[chmatch(FF$rx, Sinfo$rxDesignation)] <- forcedRx$flux
  forced_model$ub[chmatch(c(sub('_F', '_R', forcedRx$rx[grep('_F$', forcedRx$rx)]), sub('_R', '_F', forcedRx$rx[grep('_R$', forcedRx$rx)])), Sinfo$rxDesignation)]
  
  solvedModel <- gurobi(forced_model, qpparams)
  
  solvedModel
}





####### Reaction and metabolite information ##########


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


reaction_info_FBGA <- function(rxnName){
  
  require(gplots)
  
  #### similar to reaction_info but using a run_rxn file and dumping output as flat text
  rxnInfo <- run_rxn$rxnSummary
  
  rxnStoi <- rxnInfo$rxnStoi
  speciesNames <- unname(sapply(names(rxnStoi), function(x){rxnInfo$metNames[names(rxnInfo$metNames) == x]}))
  speciesName_lengths <- sapply(speciesNames, function(x){length(strsplit(x, '')[[1]])})
    
  rxnDir <- 0 #the rxn directionality was not transmitted
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
  
  if(length(strsplit(paste(substrate_prep, rxnDir, product_prep), "")[[1]]) > 60){
    
    if(length(strsplit(substrate_prep, "")[[1]]) > 60){
      substrates <- strsplit(substrate_prep, '\\+')[[1]]
      sub_long <- which.max(speciesName_lengths[rxnStoi < 0])
      
      substrates <- c(substrates[-sub_long], substrates[sub_long])
      substrate_re <- NULL
      for(i in 1:length(substrates)){
        if(i < length(substrates) - 1){
          substrate_re <- c(substrate_re, substrates[i], " +")
        }else if(i == length(substrates) - 1){
          substrate_re <- c(substrate_re, substrates[i], " +\n")
        }else{
          substrate_re <- c(substrate_re, substrates[i])
        }
      }
      substrate_prep <- paste(substrate_re, collapse = "")
      }
    
    if(length(strsplit(product_prep, "")[[1]]) > 60){
      
      products <- strsplit(product_prep, '\\+')[[1]]
      prod_long <- which.max(speciesName_lengths[rxnStoi > 0])
      
      products <- c(products[-prod_long], products[prod_long])
      products_re <- NULL
      for(i in 1:length(products)){
        if(i < length(products) - 1){
          products_re <- c(products_re, products[i], " +")
        }else if(i == length(products) - 1){
          products_re <- c(products_re, products[i], " +\n")
        }else{
          products_re <- c(products_re, products[i])
        }
      }
      product_prep <- paste(products_re, collapse = "")
      }
    
    rxn_stoi <- paste(substrate_prep, "\n\t", rxnDir, product_prep)
    }else{
    rxn_stoi <- paste(substrate_prep, rxnDir, product_prep)
    }
  rxn_stoi <- gsub(' +', ' ', rxn_stoi)
  rxn_stoi <- gsub('^ ', '', rxn_stoi)
  rxn_stoi <- gsub(' $', '', rxn_stoi)
  rxn_stoi <- gsub('\n ', '\n', rxn_stoi)
  
  textplot(paste(c(rxnInfo$reaction,
      rxnInfo$genes,
      rxn_stoi,             
      paste(strsplit(rxnInfo$pathway, split = "__")[[1]], collapse = "\n")), collapse = "\n\n")
      , cex = 3, valign = "top", halign = "left")
  }  



########## Functions used in optimization of fitted flux versus actual flux ############

species_plot <- function(run_rxn, flux_fit, chemostatInfo){
  
  require(data.table)
  
  #generate a list of four plots:
  #1: Flux predicted from FBA versus parametric form - plot1
  #2: Flux predicted from FBA versus parametric form - plot2
  #3: Flux predicted from FBA ~ GR
  #4: Relationship between metabolites/enzyme levels and condition
  #5: Relationship between metabolite/enzyme levels and flux carried
  
  output_plots <- list()
  
  scatter_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "right", 
      panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan"),
      legend.key.size = unit(3, "line"), legend.text = element_text(size = 40, face = "bold"))

  chemoInfoSubset <- chemostatInfo[chmatch(rownames(flux_fit$fitted_flux), chemostatInfo$condition),]
  DRordering <- data.frame(DRs = c(0.05, 0.11, 0.16, 0.22, 0.30), order = 1:5)
  chemoInfoSubset$DRorder <- DRordering$order[chmatch(as.character(chemoInfoSubset$DRgoal), as.character(DRordering$DRs))]
    
  
  flux_plot_alt <- data.frame(FBA = run_rxn$flux, Parametric = flux_fit$fitted_flux, condition = chemoInfoSubset$limitation, DR = chemoInfoSubset$DRorder)
  flux_range <- range(c(0, flux_plot_alt$FBA, flux_plot_alt$Parametric))
  
  output_plots$flux_plot1 <- ggplot() + scatter_theme +
    geom_text(data = flux_plot_alt, aes(x = FBA, y = Parametric, col = condition, label = DR), size = 10) + scale_size_identity() + 
    scale_color_brewer("Limitation", palette = "Set1") + ggtitle("Flux determined using FBA versus parametric form") + geom_abline(intercept = 0, slope = 1, size = 2) + 
    xlim(flux_range[1], flux_range[2]) + ylim(flux_range[1], flux_range[2])
  
  flux_plot <- data.frame(FBA = run_rxn$flux, Parametric = flux_fit$fitted_flux, condition = chemoInfoSubset$limitation, DR = chemoInfoSubset$actualDR)
  #output_plots$flux_plot2 <- ggplot() + geom_path(data = flux_plot, aes(x = FBA, y = Parametric, col = condition, size = 2)) + scatter_theme +
  #  geom_point(data = flux_plot, aes(x = FBA, y = Parametric, col = condition, size = DR*50)) + scale_size_identity() + 
  #  scale_color_brewer("Limitation", palette = "Set2") + ggtitle("Flux determined using FBA versus parametric form") + geom_abline(intercept = 0, slope = 1, size = 2) +
  #  xlim(flux_range[1], flux_range[2]) + ylim(flux_range[1], flux_range[2])
  
  FBA_flux_range <- range(c(0, flux_plot_alt$FBA))
  output_plots$FBA_flux <- ggplot() + geom_path(data = flux_plot, aes(x = DR, y = FBA, col = condition, size = 2)) + scatter_theme +
    geom_point(data = flux_plot, aes(x = DR, y = FBA, col = condition, size = DR*50)) + scale_size_identity() + 
    scale_color_brewer("Limitation", palette = "Set2") + ggtitle("Flux determined using FBA") +
    ylim(FBA_flux_range[1], FBA_flux_range[2])
    
  all_species <- data.frame(run_rxn$enzymes, run_rxn$metabolites)
  colnames(all_species) <- run_rxn$all_species$commonName
  
  all_species_tab <- run_rxn$all_species[colSums(all_species != 1) != 0,]
  all_species <- all_species[,colSums(all_species != 1) != 0]
  
  species_df <- melt(data.frame(all_species, condition = chemoInfoSubset$limitation, DR = chemoInfoSubset$actualDR, flux = run_rxn$flux), id.vars = c("condition", "DR", "flux"))
  
  output_plots$species <- ggplot() + geom_path(data = species_df, aes(x = DR, y = value, col = condition, size = 2)) + facet_wrap(~ variable, scale = "free_y") +
    scatter_theme + scale_size_identity() + scale_color_brewer("Limitation", palette = "Set2") + scale_y_continuous("Relative concentration") +
    ggtitle("Relationship between metabolites/enzyme levels and condition") + expand_limits(y = 0)
  
  output_plots$flux_species <- ggplot() + geom_path(data = species_df, aes(x = value, y = flux, col = condition, size = 2)) + facet_wrap(~ variable, scale = "free") +
    scatter_theme + scale_size_identity() + scale_color_brewer("Limitation", palette = "Set2") + scale_x_continuous("Relative concentration") +
    geom_point(data = species_df, aes(x = value, y = flux, col = condition, size =  DR*30)) + ggtitle("Relationship between metabolite/enzyme levels and flux carried") + expand_limits(y = 0)
  
  output_plots
  }






flux_fitting <- function(run_rxn, par_markov_chain, par_likelihood){
  
  require(data.table)
  require(nnls)
  
  # predict flux based upon parameter sets to determine how much variance in flux can be accounted for using the prediction
  param_interval <- exp(apply(par_markov_chain, 2, function(x){quantile(x, probs = c(0.025, 0.975))}))
  param_interval <- data.frame(cbind(t(param_interval), median = exp(apply(par_markov_chain, 2, median)), MLE = exp(par_markov_chain[which.max(par_likelihood$likelihood),])))
  
  n_c <- nrow(run_rxn$metabolites )
  
  # take the maximum likelihood estimate
  par_stack <- rep(1, n_c) %*% t(exp(par_markov_chain[which.max(par_likelihood$likelihood),]))
  colnames(par_stack) <- run_rxn$kineticPars$formulaName
  occupancy_vals <- data.frame(run_rxn$metabolites, par_stack)
  predOcc <- model.matrix(run_rxn$occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  
  enzyme_activity <- (predOcc %*% t(rep(1, sum(run_rxn$all_species$SpeciesType == "Enzyme"))))*run_rxn$enzymes #occupany of enzymes * relative abundance of enzymes
  flux_fit <- nnls(enzyme_activity, run_rxn$flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  
  fit_summary <- data.frame(residDF = sum(run_rxn$flux != 0) - length(par_stack[1,]), parametricFit = NA, NNLS = NA, LS = NA, LS_met = NA, LS_enzyme = NA, TSS = NA)
  
  ### using LS regression, how much variance is explained
  fit_summary$LS_met = anova(lm(run_rxn$flux ~ run_rxn$metabolites))$S[1]
  fit_summary$LS_enzyme = anova(lm(run_rxn$flux ~ run_rxn$enzymes))$S[1]
  fit_summary$LS = sum(anova(lm(run_rxn$flux ~ run_rxn$metabolites + run_rxn$enzymes))$S[1:2])
  fit_summary$TSS = sum(anova(lm(run_rxn$flux ~ run_rxn$metabolites + run_rxn$enzymes))$S)
  
  ### using flux fitted from the MLE parameter set, how much variance is explained
  #fit_summary$parametricFit = anova(lm(run_rxn$flux ~ flux_fit$fitted+0))$S[1]
  
  fit_summary$parametricFit <- fit_summary$TSS-sum((flux_fit$residuals)^2)
  
  ### using flux fitted using non-negative least squares regression, how much variance is explained ### metabolite abundances are corrected for whether the metabolite is a product (*-1) or reactant (*1)
  NNLSmetab <- run_rxn$metabolites * -1*(rep(1, n_c) %*% t(unname(run_rxn$rxnSummary$rxnStoi)[chmatch(colnames(run_rxn$metabolites), names(run_rxn$rxnSummary$rxnStoi))]))
  ### add potential activators and inhibitors
  if(length(run_rxn$rxnSummary$rxnFormData$Subtype[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))]) != 0){
    modifiers = data.frame(modifier = run_rxn$rxnSummary$rxnFormData$SubstrateID[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))], direction = ifelse(run_rxn$rxnSummary$rxnFormData$Type[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))] == "act", 1, -1))
    
    if(nrow(modifiers) > 1){print("fix code here to allow for multiple modifiers")}
    
    modifier_effect <- matrix(run_rxn$metabolites[,chmatch(modifiers$modifier, colnames(run_rxn$metabolites))])
    colnames(modifier_effect) <- paste(c(modifiers$modifier, "mod"), collapse = "")
    NNLSmetab <- cbind(NNLSmetab, modifier_effect)
  }
  NNLSmetab <- NNLSmetab[,apply(NNLSmetab, 2, function(x){var(x) != 0})]
  
  if(all(run_rxn$flux < 0)){
    
    nnls_fit <- nnls(as.matrix(data.frame(NNLSmetab, run_rxn$enzymes)), -1*run_rxn$flux)
    tpred <-nnls_fit$residuals
    fit_summary$NNLS <- fit_summary$TSS-sum((tpred)^2)
  
    fit_summary$nnlsPearson <- cor(-1*nnls_fit$fitted, run_rxn$flux, method = "pearson")
    fit_summary$nnlsSpearman <- cor(-1*nnls_fit$fitted, run_rxn$flux, method = "spearman")
  
    
  }else{
    
    nnls_fit <- nnls(as.matrix(data.frame(NNLSmetab, run_rxn$enzymes)), run_rxn$flux)
    tpred <- nnls_fit$residuals
    fit_summary$NNLS <- fit_summary$TSS-sum((tpred)^2)
    
    fit_summary$nnlsPearson <- cor(nnls_fit$fitted, run_rxn$flux, method = "pearson")
    fit_summary$nnlsSpearman <- cor(nnls_fit$fitted, run_rxn$flux, method = "spearman")
  
  }
   
  ### correlations
  fit_summary$parPearson <- cor(flux_fit$fitted, run_rxn$flux, method = "pearson")
  fit_summary$parSpearman <- cor(flux_fit$fitted, run_rxn$flux, method = "spearman")
  
  output <- list()
  output$fit_summary <- fit_summary
  output$param_interval <- param_interval
  output$fitted_flux <- flux_fit$fitted
  output
}


calcElast <- function(run_rxn, par_markov_chain, par_likelihood){
  
  require(reshape2)
  
  rxnData <- run_rxn$rxnSummary
  rxn <- rxnData$rxnID
  
  # Calculate the linear B enyzme coefficients by nnls
  met_abund <- run_rxn$metabolites
  
  enz_abund <- run_rxn$enzymes
  colnames(enz_abund) <- rownames(rxnData$enzymeAbund)
  
  # get create the parameternames
  tIDs <-  colnames(met_abund)
  
  # remove all tIDs not included in the reaction
  tIDs <- tIDs[tIDs %in% rxnData$rxnFormData$SubstrateID]  
  
  parNames <- unname(sapply(tIDs,function(x){
    paste(c('K_',rxnData$rxnID,'_',x),collapse='')
  }))
  
  parNames <- c(parNames,paste('Keq',rxnData$rxnID,sep=''))
  
  # find the MLE parameter values
  nEnz <- nrow(rxnData$enzymeAbund)
  n_c <- length(rxnData$flux)
  filML <- par_likelihood$likelihood == max(par_likelihood$likelihood)
  params <- exp(par_markov_chain[filML,])
  
  par_stack <- rep(1, n_c) %*% t(params); 
  colnames(par_stack) <- parNames
  
  
  # reshape to formulas in order to be suitable for calculating the sensitivities
  # remove the ~I(
  form <- paste(sub(paste(paste("E", rxnData$rxnID, sep = "_"), " \\* ", paste("V", rxnData$rxnID, sep = "_"), sep = ""), "1", rxnData$rxnForm)[2], sep = " ")
  form <- substr(form,3,nchar(form)-5)
  #form <- gsub('E_r_....\\ \\*\\ V_r_....','',form)
  
  ##Calculate the elasticities
  # write the input matrix
  
  elInpmat <- data.frame(met_abund,par_stack)
  
  predOcc <- with(elInpmat,eval(parse(text=form)))
  names(predOcc) <- rownames(elInpmat)
  
  colnames(elInpmat) <- gsub('-','',colnames(elInpmat))
  colnames(elInpmat) <- gsub('\\.','',colnames(elInpmat))
  
  elastMat <- elInpmat
  elastMat[,] <- NA
  
  eq <- eval(parse(text = paste('expression(',form,')')))
  
  #calculate elasticities
  for (fac in colnames(elInpmat)){
    #take the devirate
    dform <- D(eq,fac)
    elastMat[,fac] <- with(elInpmat,eval(dform))
  }
  
  # until now we have the sensitivities, for the elasticities divide by the value at the point,
  # and multiply by the parameter value at the point
  elastMat <- elastMat /  predOcc
  elastMat <- elastMat * elInpmat
  
  
  ## prepare things for plotting
  plotmat <- data.frame(t(elastMat))
  plotmat$ID <- row.names(plotmat)
  plotmat <- melt(plotmat, id = 'ID',variable.name='Condition',value.name='Elasticity')
  colnames(plotmat) <- c("ID", "Condition", "Elasticity")
  
  # Add the original parameter values
  
  parmat <- data.frame(t(elInpmat))
  parmat$ID <- row.names(parmat)
  parmat <- melt(parmat, id = 'ID',variable.name='Condition',value.name='Value')
  plotmat$Value <- parmat$Value
  colnames(parmat) <- c("ID", "Condition", "Value")
  
  plotmat$Type <- 'Parameter'
  plotmat[grep('^t',plotmat$ID),'Type']  <- 'Metabolite'
  
  plotmat <- plotmat[c(grep('^t',plotmat$ID),grep('^K',plotmat$ID)),]
  
  # which metabolite is the parameter associated with?
  plotmat$metID <- sapply(plotmat$ID,function(x){
    substr(x,nchar(x)-5,nchar(x))
  }) 
  
  #delete the meatbolites without a parameter
  plotmat <- plotmat[!(!(plotmat$ID %in% plotmat[plotmat$Type == 'Parameter',colnames(plotmat)=='metID']) &
                         plotmat$Type == 'Metabolite'), ]
  
  plotmat$keq <- sapply(plotmat$ID,function(x){substr(x,1,4) == 'Keqr'})
  plotmat <- plotmat[order(plotmat$metID,plotmat$Condition,plotmat$Type), ]
  
  # calculate occupancies
  plotmat$Occupancy[plotmat$Type == 'Metabolite'] <- log10(plotmat$Value[plotmat$Type == 'Metabolite']/plotmat$Value[plotmat$Type == 'Parameter' & !plotmat$keq ])
  
  # get the names of the parameters
  plotmat$Name <-sapply(plotmat$ID,function(x){
    x <-substr(x,nchar(x)-5,nchar(x))
    rxnData$metNames[x]})
  
  fil <- is.na(plotmat$Name)           
  plotmat$Name[fil] <- plotmat$ID[fil]
  plotmat$Name <- unname(unlist(plotmat$Name))
  
  # split long names to multiple line names
  plotmat$Name<- unname(sapply(plotmat$Name, function(name_int){
    if(length(strsplit(name_int, split = "")[[1]]) >= 25){
      split_name <- strsplit(name_int, split = "")[[1]]
      split_pois <- c(1:length(split_name))[split_name %in% c(" ", "-")][which.min(abs(20 - c(1:length(split_name)))[split_name %in% c(" ", "-")])]
      split_name[split_pois] <- "\n"
      paste(split_name, collapse = "")
    }else{name_int}
  }))
  
  # make a cap for the elasticities
  plotmat$Elasticity[plotmat$Elasticity > 10] <- 10
  plotmat$Elasticity[plotmat$Elasticity < -10] <- -10
  
  # make a background colouring matrix
  plotmat$Limitation <- unname(unlist(sapply(plotmat$Condition,function(x){
    substring(x,1,1)
  })))
  plotmat$Growthlab <- unname(unlist(sapply(plotmat$Condition,function(x){
    substring(x,2,nchar(as.character(x)))
  })))
  plotmat$Growth<- as.numeric(as.factor(plotmat$Growthlab))
  plotmat$Groups <- paste(plotmat$Limitation,plotmat$Type,sep='')
  
  backmat <- data.frame(xstart = seq(0.5,20.5,20/length(unique(plotmat$Limitation))), xend = seq(5.5,25.5,20/length(unique(plotmat$Limitation))), Limitation = unique(plotmat$Limitation))
  
  # reorder stuff, to get the factor order reasonable
  
  fil <- grep('Keqr_',plotmat$Name)
  plotmat <- rbind(plotmat[-fil,],plotmat[fil,])
  
  plotmat$Name <- factor(plotmat$Name,levels = unique(plotmat$Name))
  # Plot it
  ### Plot the elasticies ###
  p <- ggplot()+  
    geom_rect(data = backmat, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,fill = Limitation), alpha = 0.2) +
    scale_fill_brewer(palette='Pastel1',guide = 'none')+
    geom_line(data=plotmat,aes(x=Growth,y=Elasticity,colour=Name,linetype=Type))+
    facet_grid(~Limitation )+
    scale_colour_brewer(palette='Set1')+
    scale_x_continuous(breaks=1:5,labels=unique(plotmat$Growthlab),expand=c(0.05,0.05))+
    theme(axis.text.x = element_text(angle=90,vjust=.50),panel.background= element_rect(fill = '#EDEDED'))+
    labs(title = rxnData$reaction,x='dilution rate [1/h]')
  
  if(max(abs(plotmat$Elasticity))>1){
    p<- p+scale_y_continuous(breaks=floor(min(plotmat$Elasticity)):ceiling(max(plotmat$Elasticity)),expand=c(0.05,0.05))
  }
  
  elPlots <- list()
  elPlots$elast <- p
  
  ### Plot the occupancies ###
  p <- ggplot()+  
    geom_rect(data = backmat, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,fill = Limitation), alpha = 0.2) +
    scale_fill_brewer(palette='Pastel1',guide = 'none')+
    geom_line(data=plotmat[plotmat$Type == 'Metabolite',],aes(x=Growth,y=Occupancy,colour=Name,linetype=Type))+
    facet_grid(~Limitation )+
    scale_colour_brewer(palette='Set1')+
    scale_y_continuous(limits = c(-max(abs(plotmat$Occupancy[plotmat$Type == 'Metabolite'])),max(abs(plotmat$Occupancy[plotmat$Type == 'Metabolite']))),expand=c(0.05,0.05))+
    scale_x_continuous(breaks=1:5,labels=unique(plotmat$Growthlab),expand=c(0.05,0.05))+
    theme(axis.text.x = element_text(angle=90,vjust=.50),panel.background= element_rect(fill = '#EDEDED'))+
    labs(title = rxnData$reaction,y='log10( [met] / kM )',x='dilution rate [1/h]')
  
  elPlots$occ <- p
  return(elPlots)
} 



param_compare <- function(run_rxn, par_markov_chain, par_likelihood){
  
  ### write common names ###
  
  rename_table <- data.frame(tID = colnames(par_markov_chain), commonName = NA, commonPrint = NA)
  rename_table$commonName[rename_table$tID %in% names(run_rxn$rxnSummary$metNames)] <- unname(run_rxn$rxnSummary$metNames)[chmatch(rename_table$tID[rename_table$tID %in% names(run_rxn$rxnSummary$metNames)], names(run_rxn$rxnSummary$metNames))]
  rename_table$commonName[rename_table$tID == "keq"] <- "Keq"
  rename_table$commonPrint <- nameReformat(names = rename_table$commonName, totalChar = 120)
  
  named_par_markov_chain <- par_markov_chain
  colnames(named_par_markov_chain) <- rename_table$commonPrint
  
  # visualize the joint and marginal distribution of parameter values from the markov chain
  
  par_combinations <- expand.grid(1:length(run_rxn$kineticPars[,1]), 1:length(run_rxn$kineticPars[,1]))
  like_comparison <- ifelse(par_combinations[,1] == par_combinations[,2], TRUE, FALSE)
  
  max_likelihood <- named_par_markov_chain[which.max(par_likelihood$likelihood),]
  
  par_comp_like <- NULL
  for(i in 1:sum(like_comparison)){
    par_comp_like <- rbind(par_comp_like, data.frame(xval = named_par_markov_chain[,par_combinations[like_comparison,][i,1]], parameter_1 = colnames(named_par_markov_chain)[par_combinations[like_comparison,][i,1]],
         parameter_2 = colnames(named_par_markov_chain)[par_combinations[like_comparison,][i,1]]))
      }
  
  par_comp_dissimilar <- NULL
  for(i in 1:sum(!like_comparison)){
    par_comp_dissimilar <- rbind(par_comp_dissimilar, data.frame(xval = named_par_markov_chain[,par_combinations[!like_comparison,][i,1]], yval = named_par_markov_chain[,par_combinations[!like_comparison,][i,2]], 
          parameter_1 = colnames(named_par_markov_chain)[par_combinations[!like_comparison,][i,1]], parameter_2 = colnames(named_par_markov_chain)[par_combinations[!like_comparison,][i,2]]))
      }
  
  MLEbarplot <- data.frame(xval = max_likelihood[par_combinations[like_comparison,1]], parameter_1 = colnames(named_par_markov_chain)[par_combinations[like_comparison,1]],
      parameter_2 = colnames(named_par_markov_chain)[par_combinations[like_comparison,1]])
  MLEpoints <- data.frame(xval = max_likelihood[par_combinations[!like_comparison,1]], yval = max_likelihood[par_combinations[!like_comparison,2]],
      parameter_1 = colnames(named_par_markov_chain)[par_combinations[!like_comparison,1]], parameter_2 = colnames(named_par_markov_chain)[par_combinations[!like_comparison,2]])
  
  
  
  #### determine the maximum bin from the histogram so that values can be scaled to the bivariate histogram values ###

  par_hist_binwidth = 0.2
  
  max_density <- max(apply(named_par_markov_chain, 2, function(x){max(table(round(x/par_hist_binwidth)))}))
  
  density_trans_inv <- function(x){x*(max_density/20) + max_density/2}
  density_trans <- function(x){(x - max_density/2)/(max_density/20)}
  
  par_comp_dissimilar$yval <- density_trans_inv(par_comp_dissimilar$yval)
  MLEpoints$yval <- density_trans_inv(MLEpoints$yval)
  
  
  hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), 
      legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line"), axis.title = element_blank()) 

  ggplot() + geom_hex(data = par_comp_dissimilar, aes(x = xval, y = yval)) + geom_bar(data = par_comp_like, aes(x = xval), binwidth = par_hist_binwidth, col = "black") + facet_grid(parameter_2 ~ parameter_1, scales = "fixed") + hex_theme +
    scale_fill_gradientn(name = "Counts", colours = c("white", "darkgoldenrod1", "chocolate1", "firebrick1", "black")) +
    scale_x_continuous(NULL, expand = c(0.02,0.02)) + scale_y_continuous(NULL, expand = c(0.01,0.01), labels = density_trans, breaks = density_trans_inv(seq(-10, 10, by = 5))) +
    geom_vline(data = MLEbarplot, aes(xintercept = xval), col = "cornflowerblue", size = 2) + geom_point(data = MLEpoints, aes(x = xval, y = yval), size = 2, col = "cornflowerblue")

  }


calcElast <- function(run_rxn, par_markov_chain, par_likelihood){
  
  rxnData <- run_rxn$rxnSummary
  rxn <- rxnData$rxnID
  
  # Calculate the linear B enyzme coefficients by nnls
  met_abund <- run_rxn$metabolites
  n_c <- nrow(met_abund)
  
  enz_abund <- run_rxn$enzymes
  colnames(enz_abund) <- rownames(rxnData$enzymeAbund)
  
  # get create the parameternames
  tIDs <- colnames(met_abund)
  
  # remove all tIDs not included in the reaction
  tIDs <- tIDs[tIDs %in% rxnData$rxnFormData$SubstrateID]
  
  parNames <- unname(sapply(tIDs,function(x){
    paste(c('K_',rxnData$rxnID,'_',x),collapse='')
  }))
  
  parNames <- c(parNames,paste('Keq',rxnData$rxnID,sep=''))
  
  # find the MLE parameter values
  nEnz <- nrow(rxnData$enzymeAbund)
  n_c <- length(rxnData$flux)
  filML <- par_likelihood$likelihood == max(par_likelihood$likelihood)
  params <- exp(par_markov_chain[filML,])
  
  par_stack <- rep(1, n_c) %*% t(params);
  colnames(par_stack) <- parNames
  
  
  # reshape to formulas in order to be suitable for calculating the sensitivities
  # remove the ~I(
  form <- paste(sub(paste(paste("E", rxnData$rxnID, sep = "_"), " \\* ", paste("V", rxnData$rxnID, sep = "_"), sep = ""), "1", rxnData$rxnForm)[2], sep = " ")
  form <- substr(form,3,nchar(form)-5)
  #form <- gsub('E_r_....\\ \\*\\ V_r_....','',form)
  
  ##Calculate the elasticities
  # write the input matrix
  
  elInpmat <- data.frame(met_abund,par_stack)
  
  predOcc <- with(elInpmat,eval(parse(text=form)))
  names(predOcc) <- rownames(elInpmat)
  
  colnames(elInpmat) <- gsub('-','',colnames(elInpmat))
  colnames(elInpmat) <- gsub('\\.','',colnames(elInpmat))
  
  elastMat <- elInpmat
  elastMat[,] <- NA
  
  eq <- eval(parse(text = paste('expression(',form,')')))
  
  #calculate elasticities
  for (fac in colnames(elInpmat)){
    #take the devirate
    dform <- D(eq,fac)
    elastMat[,fac] <- with(elInpmat,eval(dform))
  }
  
  # until now we have the sensitivities, for the elasticities divide by the value at the point,
  # and multiply by the parameter value at the point
  elastMat <- elastMat / predOcc
  elastMat <- elastMat * elInpmat
  
  
  ## prepare things for plotting
  plotmat <- data.frame(t(elastMat))
  plotmat$ID <- row.names(plotmat)
  plotmat <- melt(plotmat, id = 'ID',variable.name='Condition',value.name='Elasticity')
  
  # Add the original parameter values
  
  parmat <- data.frame(t(elInpmat))
  parmat$ID <- row.names(parmat)
  parmat <- melt(parmat, id = 'ID',variable.name='Condition',value.name='Value')
  plotmat$Value <- parmat$Value
  
  plotmat$Type <- 'Parameter'
  plotmat[grep('^t',plotmat$ID),'Type'] <- 'Metabolite'
  
  plotmat <- plotmat[c(grep('^t',plotmat$ID),grep('^K',plotmat$ID)),]
  
  # which metabolite is the parameter associated with?
  plotmat$metID <- sapply(plotmat$ID,function(x){
    substr(x,nchar(x)-5,nchar(x))
  })
  
  #delete the meatbolites without a parameter
  plotmat <- plotmat[!(!(plotmat$ID %in% plotmat[plotmat$Type == 'Parameter',colnames(plotmat)=='metID']) &
                         plotmat$Type == 'Metabolite'), ]
  
  plotmat$keq <- sapply(plotmat$ID,function(x){substr(x,1,4) == 'Keqr'})
  plotmat <- plotmat[order(plotmat$metID,plotmat$Condition,plotmat$Type), ]
  
  # calculate occupancies
  plotmat$Occupancy[plotmat$Type == 'Metabolite'] <- log10(plotmat$Value[plotmat$Type == 'Metabolite']/plotmat$Value[plotmat$Type == 'Parameter' & !plotmat$keq ])
  
  # get the names of the parameters
  plotmat$Name <-sapply(plotmat$ID,function(x){
    x <-substr(x,nchar(x)-5,nchar(x))
    rxnData$metNames[x]})
  
  fil <- is.na(plotmat$Name)
  plotmat$Name[fil] <- plotmat$ID[fil]
  plotmat$Name <- unname(unlist(plotmat$Name))
  
  # split long names to multiple line names
  plotmat$Name<- unname(sapply(plotmat$Name, function(name_int){
    if(length(strsplit(name_int, split = "")[[1]]) >= 25){
      split_name <- strsplit(name_int, split = "")[[1]]
      split_pois <- c(1:length(split_name))[split_name %in% c(" ", "-")][which.min(abs(20 - c(1:length(split_name)))[split_name %in% c(" ", "-")])]
      split_name[split_pois] <- "\n"
      paste(split_name, collapse = "")
    }else{name_int}
  }))
  
  # make a cap for the elasticities
  plotmat$Elasticity[plotmat$Elasticity > 10] <- 10
  plotmat$Elasticity[plotmat$Elasticity < -10] <- -10
  
  # make a background colouring matrix
  plotmat$Limitation <- unname(unlist(sapply(plotmat$Condition,function(x){
    substring(x,1,1)
  })))
  plotmat$Growthlab <- unname(unlist(sapply(plotmat$Condition,function(x){
    substring(x,2,nchar(as.character(x)))
  })))
  plotmat$Growth<- as.numeric(as.factor(plotmat$Growthlab))
  plotmat$Groups <- paste(plotmat$Limitation,plotmat$Type,sep='')
  
  backmat <- data.frame(xstart = seq(0.5,20.5,5), xend = seq(5.5,25.5,5), Limitation = unique(plotmat$Limitation))
  
  # reorder stuff, to get the factor order reasonable
  
  fil <- grep('Keqr_',plotmat$Name)
  plotmat <- rbind(plotmat[-fil,],plotmat[fil,])
  
  plotmat$Name <- factor(plotmat$Name,levels = unique(plotmat$Name))
  # Plot it
  ### Plot the elasticies ###
  p <- ggplot()+
    geom_rect(data = backmat, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,fill = Limitation), alpha = 0.2) +
    scale_fill_brewer(palette='Pastel1',guide = 'none')+
    geom_line(data=plotmat,aes(x=Growth,y=Elasticity,colour=Name,linetype=Type))+
    facet_grid(~Limitation )+
    scale_colour_brewer(palette='Set1')+
    scale_x_continuous(breaks=1:5,labels=unique(plotmat$Growthlab),expand=c(0.05,0.05))+
    theme(axis.text.x = element_text(angle=90,vjust=.50),panel.background= element_rect(fill = '#EDEDED'))+
    labs(title = rxnData$reaction,x='dilution rate [1/h]')
  
  if(max(abs(plotmat$Elasticity))>1){
    p<- p+scale_y_continuous(breaks=floor(min(plotmat$Elasticity)):ceiling(max(plotmat$Elasticity)),expand=c(0.05,0.05))
  }
  
  elPlots <- list()
  elPlots$elast <- p
  
  ### Plot the occupancies ###
  p <- ggplot()+
    geom_rect(data = backmat, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,fill = Limitation), alpha = 0.2) +
    scale_fill_brewer(palette='Pastel1',guide = 'none')+
    geom_line(data=plotmat[plotmat$Type == 'Metabolite',],aes(x=Growth,y=Occupancy,colour=Name,linetype=Type))+
    facet_grid(~Limitation )+
    scale_colour_brewer(palette='Set1')+
    scale_y_continuous(limits = c(-max(abs(plotmat$Occupancy[plotmat$Type == 'Metabolite'])),max(abs(plotmat$Occupancy[plotmat$Type == 'Metabolite']))),expand=c(0.05,0.05))+
    scale_x_continuous(breaks=1:5,labels=unique(plotmat$Growthlab),expand=c(0.05,0.05))+
    theme(axis.text.x = element_text(angle=90,vjust=.50),panel.background= element_rect(fill = '#EDEDED'))+
    labs(title = rxnData$reaction,y='log10( [met] / kM )',x='dilution rate [1/h]')
  
  elPlots$occ <- p
  return(elPlots)
} 




par_draw <- function(updates){
  #### update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ####
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
      } else if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
      }
    }
  draw
  }


lik_calc <- function(proposed_params){
  #### determine the likelihood of predicted flux as a function of metabolite abundance and kinetics parameters relative to actual flux ####
  
  par_stack <- rep(1, n_c) %*% t(proposed_params); colnames(par_stack) <- kineticPars$formulaName
  par_stack <- exp(par_stack)
  occupancy_vals <- data.frame(met_abund, par_stack)
  
  predOcc <- model.matrix(occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*enzyme_abund #occupany of enzymes * relative abundance of enzymes
  
  flux_fit <- nnls(enzyme_activity, flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  fit_resid_error <- sqrt(mean((flux_fit$resid - mean(flux_fit$resid))^2))
  
  sum(dnorm(flux, flux_fit$fitted, fit_resid_error, log = TRUE))
  
  }


nameReformat <- function(names, totalChar){
  ### split names into multiple lines if they are too long
  
  nameLengths <- data.frame(name = rename_table$commonName, nchar = sapply(rename_table$commonName, function(x){length(strsplit(x, "")[[1]])}))
  
  if(!all(nameLengths$nchar < totalChar/nrow(nameLengths))){
    for(i in c(1:nrow(nameLengths))[nameLengths$nchar >= totalChar/nrow(nameLengths)]){
      longName <- nameLengths$name[i]
      #find location of space or - 
      split_pos <- rev(gregexpr('[ -]', longName)[[1]])[which.min(abs(rev(gregexpr('[ -]', longName)[[1]]) - nameLengths$nchar[i]/2))]
      substr(longName, split_pos, split_pos) <- "#"
      nameLengths$name[i] <- gsub("#", '\\\n', longName)
    }
  }
  nameLengths$name
}
  