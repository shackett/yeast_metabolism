
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

  rxStoi <- rxnFile[!is.na(rxnFile$StoiCoef),]
	for (i in 1:nrow(rxStoi)){
		stoiMat[c(1:length(metabolites))[metabolites == rxStoi$Metabolite[i]], c(1:length(reactions))[reactions == rxStoi$ReactionID[i]]] <- as.numeric(rxStoi$StoiCoef[i])
	  }

	if (internal_names == FALSE){
		save(stoiMat, file = "flux_cache/yeast_stoi.Rdata")} else {
		rownames(stoiMat) <- metabolites	
		colnames(stoiMat) <- reactions
		save(stoiMat, file = "flux_cache/yeast_stoi.Rdata")
			}
	}

parse_custom = function(customRx){
  
  ### read supplemental reaction files and append to input data.frames ###
  
  outputList <- list()
  
  inputFile = read.table(customRx, header = F, sep = "\t", fill = T, blank.lines.skip = F)
  
  ### Species-level annotation ###
  
  input_bound = c(1:nrow(inputFile))[inputFile[,1] == "!Species"] + c(1,-1)
  spec_input = inputFile[(input_bound[1]+1):input_bound[2],colSums(inputFile[(input_bound[1]+1):input_bound[2],] != "") != 0]
  colnames(spec_input) <- inputFile[input_bound[1],colSums(inputFile[(input_bound[1]+1):input_bound[2],] != "") != 0]
  
  corr = spec_input[,colnames(spec_input) %in% c("SpeciesID", "SpeciesName", "SpeciesType", "Compartment")]
  outputList$corrFile = corr[match(unique(corr$SpeciesType), corr$SpeciesType),]

  spec_par = spec_input[,colnames(spec_input) %in% c("SpeciesID", "SpeciesName", "Annotation")]
  outputList$specparFile = spec_par
  
  ### Reaction stoichiometry annotation ###
  
  input_bound = c(1:nrow(inputFile))[inputFile[,1] == "!Reactions"] + c(1,-1)
  rxn_input = inputFile[(input_bound[1]+1):input_bound[2],]
  colnames(rxn_input) <- inputFile[input_bound[1],]
  
  outputList$rxnFile = rxn_input
  
  ### Reaction flux and references annotation ###
  
  input_bound = c(1:nrow(inputFile))[inputFile[,1] == "!Reaction_Parameters"] + c(1,-1)
  rxPar_input = inputFile[(input_bound[1]+1):input_bound[2],apply(inputFile[(input_bound[1]+1):input_bound[2],], 2, function(x){sum(x[!is.na(x)] != "")}) != 0]
  colnames(rxPar_input) <- inputFile[input_bound[1],apply(inputFile[(input_bound[1]+1):input_bound[2],], 2, function(x){sum(x[!is.na(x)] != "")}) != 0]
  
  outputList$rxnparFile = rxPar_input[,colnames(rxPar_input) %in% c("ReactionID", "Enzymes", "Annotation")]
  outputList$fluxDirFile = rxPar_input[,colnames(rxPar_input) %in% c("ReactionID", "Reversible", "FluxBound")]
  
  return(outputList)
  
  }


perfect.match <- function(source, query, corrFile, reduceByLength = F){
  all_char <- "[[:graph:][:space:]]"
  
  tmp <- corrFile[grep(source, query)[!(grep(source, query) %in% union(grep(paste(all_char, source, sep = ""), query), grep(paste(source, all_char, sep = ""), query)))],]
  if(length(tmp[,1]) == 0){tmp <- corrFile[grep(source, query, fixed = TRUE),]}
  
  if(reduceByLength){
    tmp[which.min(abs(sapply(tmp$SpeciesName, function(x){length(unlist(strsplit(x, split = "")))}) - length(unlist(strsplit(source, split = ""))))),]
    
    }else{
      tmp
    }
  
}


rxn_search = function(search_string, stoiMat = named_stoi, is_rxn = TRUE, index = FALSE){
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
				reactions[colz]
				} else {
			
			rxns = stoiMat[,colz]
			if(is.vector(rxns)){
				c(reactions[colz], rxns[rxns != 0])
				} else {
					output <- rbind(reactions[colz], rxns[apply(rxns, 1, is.not.zero),])
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

refine_stoi <- function(stoiMat, reversibleRx, modelMetComp, measured_bounds){
  
  require(data.table)
  
  #### Two tasks using stoichiometry as a graph: 1) Look at stoichiometric connections and identify loops that only differ in cofactor usage. 2) Find and remove orphaned nodes ####
  
  validSinks = c(comp_by_cond$compositionFile$AltName[comp_by_cond$compositionFile$varCategory != "Maintenance ATP hydrolysis"], unique(measured_bounds$specie[measured_bounds$type == "excretion"]))
  validSinks = modelMetComp$ID[chmatch(validSinks, modelMetComp$name)]
  
  ### remove reaction which cannot carry flux because of stoichiometry structure ###
  
  #extraneous_mets <- modelMetComp$ID[sort(union(c(1:nrow(modelMetComp))[!is.na(modelMetComp$C) & modelMetComp$C ==0], grep('^NAD|^ATP|dioxide',  modelMetComp$name)))]
  #core_rxnFile <- rxnFile[!is.na(rxnFile$StoiCoef) & !(rxnFile$Metabolite %in% extraneous_mets),]
  core_rxnFile <- rxnFile[!is.na(rxnFile$StoiCoef),]
  core_rxnFile$reversible = reversibleRx$reversible[chmatch(core_rxnFile$ReactionID, reversibleRx$rx)]

  edgeDF <- rbind(
    data.frame(source = core_rxnFile$Metabolite[core_rxnFile$StoiCoef < 0 & core_rxnFile$reversible == 1], dest = core_rxnFile$ReactionID[core_rxnFile$StoiCoef < 0 & core_rxnFile$reversible == 1]),
    data.frame(source = core_rxnFile$ReactionID[core_rxnFile$StoiCoef > 0 & core_rxnFile$reversible == 1], dest = core_rxnFile$Metabolite[core_rxnFile$StoiCoef > 0 & core_rxnFile$reversible == 1]),
    data.frame(source = core_rxnFile$ReactionID[core_rxnFile$reversible == 0], dest = core_rxnFile$Metabolite[core_rxnFile$reversible == 0]),
    data.frame(source = core_rxnFile$Metabolite[core_rxnFile$reversible == 0], dest = core_rxnFile$ReactionID[core_rxnFile$reversible == 0])
    )
  
  g <- graph.empty() + vertices(union(edgeDF$source, edgeDF$dest))
  g[from=edgeDF$source, to=edgeDF$dest] <- TRUE
  
  #### find all enzymes which can reach a valid biomass sink ###
  
  validSinks <- validSinks[validSinks %in% V(g)$name] #some inorganic or cofactor species removed
  validSinks = data.frame(sID = validSinks, vertex = c(1:vcount(g))[chmatch(validSinks, V(g)$name)])
  
  
  feedable_vertices <- NULL
  for(a_sink in 1:nrow(validSinks)){
    feedable_vertices[[a_sink]] <- subcomponent(g, validSinks$vertex[a_sink], "in")
    }
  
  all_feedable_vertices <- sort(unique(unlist(feedable_vertices)))
  
  orphaned_rxns <- rxnIDtoEnz(grep('^r_', V(g)$name[-all_feedable_vertices], value = T)) #reactions which cannot carry flux into a biomass component
  orphaned_mets <- metIDtoSpec(grep('^s_', V(g)$name[-all_feedable_vertices], value = T)) #metabolites which have no valid paths to biomass
  
  ### Identify loops - strongly connected clusters 
  
  extraneous_mets <- modelMetComp$ID[sort(union(c(1:nrow(modelMetComp))[!is.na(modelMetComp$C) & modelMetComp$C == 0], grep('^NAD|^[CGA]TP \\[|^[CGA]DP \\[|^[CGA]MP \\[|dioxide|biomass|^lipid|^coenzyme A',  modelMetComp$name)))]
  core_rxnFile <- rxnFile[!is.na(rxnFile$StoiCoef) & !(rxnFile$Metabolite %in% extraneous_mets),]
  core_rxnFile$reversible = reversibleRx$reversible[chmatch(core_rxnFile$ReactionID, reversibleRx$rx)]
  core_rxnFile <- core_rxnFile[grep('biomass', core_rxnFile$Reaction, invert = T),] # remove boundary reactions
  
  edgeDF <- rbind(
    data.frame(source = core_rxnFile$Metabolite[core_rxnFile$StoiCoef < 0 & core_rxnFile$reversible == 1], dest = core_rxnFile$ReactionID[core_rxnFile$StoiCoef < 0 & core_rxnFile$reversible == 1]),
    data.frame(source = core_rxnFile$ReactionID[core_rxnFile$StoiCoef > 0 & core_rxnFile$reversible == 1], dest = core_rxnFile$Metabolite[core_rxnFile$StoiCoef > 0 & core_rxnFile$reversible == 1]),
    data.frame(source = core_rxnFile$ReactionID[core_rxnFile$reversible == 0], dest = core_rxnFile$Metabolite[core_rxnFile$reversible == 0]),
    data.frame(source = core_rxnFile$Metabolite[core_rxnFile$reversible == 0], dest = core_rxnFile$ReactionID[core_rxnFile$reversible == 0])
    )
  
  g <- graph.empty() + vertices(union(edgeDF$source, edgeDF$dest))
  g[from=edgeDF$source, to=edgeDF$dest] <- TRUE
  
  ### remove high betweenness pathological reactions ###
  g_between <- betweenness(g)
  high_btwness <- sort(g_between, decreasing = T)[1:100]
  plot(high_btwness)
  
  g_degree <- degree(g)
  high_degree <- sort(g_degree, decreasing = T)[1:100]
  metIDtoSpec(names(high_degree))
  
  ### remove nodes that are isolated through weak connectivity - not even connected ... ###
  
  graph_clusters_weak <- clusters(g, mode = "weak")
  isolated_rxn_ex <- c(1:vcount(g))[graph_clusters_weak$membership == which(graph_clusters_weak$csize == 5)[1]]
  isolated_rxn_edges <- edgeDF[unique(unlist(sapply(isolated_rxn_ex, function(x){incident(g, x, mode = "all")}))),] # get all adjacent edges
  edgeDF[edgeDF$source %in% unlist(isolated_rxn_edges) | edgeDF$dest %in% unlist(isolated_rxn_edges),] # confirm weak connectivity
  
  g <- delete.vertices(g, c(1:vcount(g))[graph_clusters_weak$membership != which.max(graph_clusters_weak$csize)])
  
  ### isolate strongly connected modules - identifying loops ###
  
  
  graph_clusters_strong <- clusters(g, mode = "strong")
  
  strong_clusters <- data.frame(set = 0, vertex = V(g)$name[graph_clusters_strong$membership %in% c(1:graph_clusters_strong$no)[graph_clusters_strong$csize != 1]], cluster = graph_clusters_strong$membership[graph_clusters_strong$membership %in% c(1:graph_clusters_strong$no)[graph_clusters_strong$csize != 1]])
  gred <- delete.vertices(g, c(1:vcount(g))[V(g)$name %in% names(high_btwness)])
  
  clusters(gred, mode = "strong")
  
}



######### Functions to convert between IDs and species ######

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
	if(length(grep("chebi", specparFile$Annotation[specparFile$SpeciesID == mets])) == 0){
		NA
		}else{
	  unlist(strsplit(specparFile$Annotation[specparFile$SpeciesID == mets], split = "CHEBI:"))[2]	
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

elemental_composition <- function(metabolites){
    
    ### Determine the elemental composition of all metabolites with valid ChEBI IDs ###
    ### & write to a table ###
  
    metComp <- read.delim('../Yeast_reconstruction/Sequences/METeleComp.tsv') # from seq_parser.py > elements_formula
    
    met_chebi <- unlist(sapply(metabolites, metToCHEBI))
    
    
    met_chebi_comp <- metComp[metComp$ID %in% met_chebi,]
    met_chebi_comp <- met_chebi_comp[,c(TRUE, TRUE, apply(met_chebi_comp[,-c(1,2)], 2, sum) != 0)]
    
    ele_comp <- lapply(met_chebi, function(x){
      if(length(met_chebi_comp[met_chebi_comp$ID %in% x,]) != 0){
        met_chebi_comp[met_chebi_comp$ID %in% x,]}
    })
    
    #matrix of elemental abundance of species corresponding to rows of the stoichiometric matrix
    ele_comp_mat <- matrix(NA, nrow = length(metabolites), ncol = length(ele_comp[[1]])-2); rownames(ele_comp_mat) <- metabolites; colnames(ele_comp_mat) <- names(ele_comp[[1]])[-c(1:2)]
    for(el in 1:length(metabolites)){
      if(length(unlist(ele_comp[[el]][-c(1:2)])) != 0){
        ele_comp_mat[el,] <- unlist(ele_comp[[el]][-c(1:2)])
      }
    }
    
    write.table(data.frame(ID = metabolites, name = unname(metIDtoSpec(metabolites)),ele_comp_mat), file = "flux_cache/stoiMetsComp.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = F) #dump for boundary condition determination in boundaryDataBlender.R
  }

convert_to_elemental <- function(redStoi, modelMetComp){
  # convert a stoichiometric matrix (with named columns) and "S" IDs as rownames to their elemental constituents
  require(data.table)
  
  redStoi <- boundary_stoichiometry
  specieEle <- modelMetComp[chmatch(rownames(redStoi), modelMetComp$ID),]
  
  specieEle[specieEle$name == "(1->3)-beta-D-glucan [cytoplasm]", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 10, 5) #one oxygen shared bc of condensation
  specieEle[specieEle$name == "glycogen [cytoplasm]", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 10, 5)
  specieEle[specieEle$name == "mannan [cytoplasm]", colnames(specieEle) %in% c("C", "H", "O")] <- c(6, 10, 5)
  specieEle[specieEle$name == "polyphosphate [cytoplasm]", colnames(specieEle) %in% c("O", "P")] <- c(4, 1)
  specieEle[specieEle$name == "complex sphingolipid [cytoplasm]", colnames(specieEle) %in% c("C", "H", "N", "O")] <- c(16, 35, 1, 2)  # this is primarily a begnign standin so the MW is not zero resulting in conversion problems

  specieEle[is.na(specieEle)] <- 0
  specieEle[,-c(1:2)] <- apply(specieEle[,-c(1:2)], c(1,2), as.numeric)
  specieEle <- specieEle[,c(TRUE, TRUE, colSums(specieEle[,-c(1:2)]) != 0)]
  
  atomicMasses <- data.frame(element = c("C", "H", "N", "O", "P", "R", "S"), mass = c(12.0107, 1.00794, 14.00674, 15.9994, 30.973761, 0, 32.066))
  atomicMasses <- atomicMasses[atomicMasses$element %in% colnames(specieEle[,-c(1,2)]),]
  
  if(all(atomicMasses$element == colnames(specieEle[,-c(1,2)]))){
    specieEle$MW <- t(t(specieEle[,-c(1,2)])) %*% t(t(c(atomicMasses[,2])))
    specieEle
    }else{
      print("atomicMasses are mismatched, additional elemental masses need to be included")
      }
  
  }


gene_pathways <- function(kegg_enzyme_dict){
  #generate a per-gene pathway annotation if one is not already generated
  
  genes_to_pathways = read.delim("http://rest.kegg.jp/link/pathway/sce", header = FALSE); colnames(genes_to_pathways) <- c("gene", "pathwayCode")
  pathway_names = read.delim("http://rest.kegg.jp/list/pathway/sce", header = FALSE); colnames(pathway_names) <- c("pathwayCode", "pathway")
  genes_to_pathways$gene <- gsub('sce:', '', genes_to_pathways$gene)
  pathway_names$pathway <- sapply(pathway_names$pathway, function(x){strsplit(x, split = " - Sacc")[[1]][1]})
  pathway_names$pathway <- sapply(pathway_names$pathway, function(x){strsplit(x, split = " - yeast")[[1]][1]})
  pathway_match <- merge(genes_to_pathways, pathway_names)
  kegg_enzyme_dict$PATHWAY <- sapply(kegg_enzyme_dict$SYST, function(x){
    paste(pathway_match$pathway[pathway_match$gene == x], collapse = "__")
  })
  write.table(kegg_enzyme_dict, "../KEGGrxns/yeastNameDict.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }




trackMetConversion <- function(trackedMet, allRxns = T){
  
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
  reducedStoi[chmatch(rownames(boundaryStoi), rownames(reducedStoi)),is.na(index_match)] <- as.matrix(boundaryStoi)
  
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

##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@
######## Diagnostic flux prediction methods ##########
##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@


FVA_setup <- function(aModel, useCluster){

  ### Setup a model that can be used to run Gurobi through python - specifically implementing flux variability analysis ###
  
  # This is a character vector with length:
  # tables are sent row by row + 1 row for colnames, 1 for rownames
  # vectors as 1 row
  # each row has the variable name as first entry
  # variables: A(tab), rhs(vec), sense(vec), lb(vec), ub(vec), Q(tab), obj(vec), parameters(vec)
  # namesrxndG(vec),CI95rxndG(vec),namesthdynMet(vec),lbthdynMet(vec),ubthdynMet(vec),pythonMode(vec)  
  
  if (useCluster != 'load'){
    
    # 10 vec + 4 row/col names = 14 single rows
    lenPyt <- 15 + nrow(aModel$A) + nrow(aModel$Q)
    
    pythonDat <- vector(mode="character",length=lenPyt)
    pos =1
    
    pythonDat[pos] = paste(c('mode',pythonMode),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('FVA',FVA),collapse='\t')
    pos <- pos+1
    
    for(i in 1:nrow(aModel$A)){
      pythonDat[pos] = paste(c('A',aModel$A[i,]),collapse='\t')
      pos <- pos+1
    }
    pythonDat[pos] = paste(c('rownamesA',rownames(aModel$A)),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('colnamesA',colnames(aModel$A)),collapse='\t')
    pos <- pos+1
    
    for(i in 1:nrow(aModel$Q)){
      pythonDat[pos] = paste(c('Q',aModel$Q[i,]),collapse='\t')
      pos <- pos+1
    }
    pythonDat[pos] = paste(c('rownamesQ',rownames(aModel$Q)),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('colnamesQ',colnames(aModel$Q)),collapse='\t')
    pos <- pos+1
    
    pythonDat[pos] = paste(c('rhs',aModel$rhs),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('sense',aModel$sense),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('lb',aModel$lb),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('ub',aModel$ub),collapse='\t')
    pos <- pos+1
    pythonDat[pos] = paste(c('obj',aModel$obj),collapse='\t')
  }  
  
  if(!file.exists("./Gurobi_python")){
    dir.create("./Gurobi_python")
  }
  if(!file.exists("./Gurobi_python/test_files")){
    dir.create("./Gurobi_python/test_files")
  }
    
  if (useCluster == 'write'){
    write(pythonDat,file=paste('./Gurobi_python/test_files/pythonDat_',treatment,'.txt',sep=''),sep='\t')
    pythout = NULL
  } else if(useCluster == 'load'){
    pythout = (read.table(paste('./Gurobi_python/test_files/pythout_',treatment,'.txt',sep=''),sep='$'))[,1]
  } else {
    pythout <- system('python ./Gurobi_python/qp_fba_clust.py',intern=T,input=pythonDat)
  }
  return(pythout)
}


FVA_read <- function(pythout){
  
  ### Parse flux variability model which was run in python ###
  
  idxS = which(pythout == 'output_start')+1
  idxE = which(pythout == 'output_end')-1
  if (length(idxE)==0){
    idxE = length(pythout)
  }
  
  solvedModel = list()
  modelOut <- data.frame(param = sapply(pythout[idxS:idxE],function(x){strsplit(x,'\t')[[1]][1]}))
  modelOut$value <- as.numeric(sapply(pythout[idxS:idxE],function(x){strsplit(x,'\t')[[1]][2]}))
  modelOut$runStatus <- sapply(pythout[idxS:idxE],function(x){sub('status_', '', strsplit(x,'\t')[[1]][3])}) # extract optimization status codes for FVA
  modelOut$violation <- as.numeric(sapply(pythout[idxS:idxE],function(x){sub('viol_', '', strsplit(x,'\t')[[1]][4])}))
  
  solvedModel$x <- modelOut$value[ modelOut$param %in% colnames(qpModel$A)]
  print(1)
  modelOut$Type <- sapply(modelOut$param,function(x){
    if (substr(x,1,5) =='RTlnc'){
      return('conc')
    } else if (substr(x,1,3) =='dGr'){
      if(substr(x,nchar(x),nchar(x))=='R'){
        return('dGrR')
      } else {
        return('dGrF')
      }
    } else if (substr(x,nchar(x)-1,nchar(x)) == 'sw'){
      return('sw')
    } else if (substr(x,nchar(x)-5,nchar(x)) == 'offset'){
      return('offset')
    } else if(substr(x,nchar(x),nchar(x))=='F'){
      return('rxnF')
    } else if(substr(x,nchar(x),nchar(x))=='R'){
      return('rxnR')
    } else if (substr(x,nchar(x)-2,nchar(x)) == 'min'){
      return(paste('FVAmin',gsub('(.*_FVA_)||(_min)','',x),sep='_'))
    } else if (substr(x,nchar(x)-2,nchar(x)) == 'max'){
      return(paste('FVAmax',gsub('(.*_FVA_)||(_max)','',x),sep='_'))
    } else return ('rxnF')
  })
  print(1.5)
  # what is the type described by the parameter (is it associated with a reaction or metabolites (eg RTln))
  modelOut$asType <- sapply(modelOut$param,function(x){
    if(substr(x,1,5) =='RTlnc'){
      return('met')
    } else {
      return('rxn')
    }
  })
  
  modelOut$asID <- sapply(1:length(modelOut$param),function(x){
    name =modelOut$param[x]
    name = gsub('_match','',name)
    if (modelOut$asType[x] == 'rxn'){
      idx = regexec('r_...._',name)[[1]][1]
      if (idx >0){
        return(substr(name,idx,idx+5))
      } else {
        if (modelOut$Type[x] == 'sw'){
          return(substr(name,1,nchar(name)-3))
        } else if (modelOut$Type[x] == 'offset'){
          return(substr(name,1,nchar(name)-7))
        } else if (modelOut$Type[x] %in% c('dGrR','dGrF')){
          return(substr(name,5,nchar(name)-2))
        } else if (grepl('FVA',modelOut$Type[x])){
          return(strsplit(name,'_FVA')[[1]][1])
        } else {
          return(substr(name,1,nchar(name)-2))
        }
      }
    }else { # metabolites
      idx = regexec('_s_',name)[[1]][1]
      return(substr(name,idx+1,nchar(name)))
    }
  })
  
  modelOut = reshape(modelOut,idvar = c('asID','asType'),timevar='Type',drop =c('param'),direction ='wide')
  colnames(modelOut) <- gsub('value\\.','',colnames(modelOut))
  modelOut <- modelOut[,colSums(!(is.na(modelOut)|modelOut=="") ) != 0] # remove optimization status codes for fluxes which were not posed as an individual optimizaiton problem
  
  modelOut$rxnNet <- (sapply(modelOut$rxnF,function(x){max(c(x,0),na.rm=T)})
                      -sapply(modelOut$rxnR,function(x){max(c(x,0),na.rm=T)}) +
                        sapply(modelOut$offset,function(x){max(c(x,0),na.rm=T)}))/flux_elevation_factor
  modelOut$rxnNet[abs(modelOut$rxnNet) < 1e-19] = 0
  modelOut$name <- sapply(modelOut$asID,function(x){
    if (substr(x,1,2) == 'r_'){
      return(rxnIDtoEnz(x))
    } else if (substr(x,1,2) == 's_'){
      return(metIDtoSpec(x))
    } else return(x)
  })
  modelOut$comparID <- sapply(modelOut$asID,function(x){
    if (substr(x,1,2) == 'r_'){
      return(rxnFile$Compartment[ rxnFile$ReactionID == x][1])
    } else (return(NA))
  })
  
  modelOut
  }


FVA_summary <- function(modelOut){
  
  ### Using python-gurobi output create a clean output showing flux bounds relative to objective value ###
  
  idx <- grepl('FVA',colnames(modelOut))
  modelOut[,idx] <- modelOut[,idx]/flux_elevation_factor
  modelOut[,idx][abs(modelOut[,idx])< 1e-19 & !is.na(modelOut[,idx])] = 0
  #modelOut$FVAvar <- (modelOut$FVAmax-modelOut$FVAmin)/abs(modelOut$rxnNet)
  
  # look at backwards reactions
  fil = modelOut$rxnNet < 0
  
  # write reactions out for fluxviz
  fil <- modelOut$asType =='rxn' & abs(modelOut$rxnNet) != 0
  fil <- modelOut$asType =='rxn' & abs(modelOut$rxnNet) >0.01
  #hist(modelOut$rxnNet[fil],xlim=c(-0.02,0.02),ylim=c(0,10),breaks=2000,freq=F)
  #write.table(modelOut[fil,c('asID','rxnNet')],file='fluxviz.val',sep='\t',
  # row.names=F,col.names=F,quote=F)
  modelOutL[[treatment]] <- modelOut
  # finite fva: a = modelOut[ !is.finite(modelOut$FVAmax) | !is.finite(modelOut$FVAmin) ,]
  # reactions by name:a = modelOut[ grep('glucose',modelOut$name) ,]
  #well constraint rxn:
  print(2)
  
  ## make the 3D FVA matrix
  FVAcuts <- sapply(colnames(modelOut)[grepl('FVA.*min',colnames(modelOut))],function(x){
    strsplit(x,'_')[[1]][2]
  })
  FVAcuts <- as.numeric(FVAcuts)
  names(FVAcuts) <- colnames(modelOut)[grepl('FVA.*min',colnames(modelOut))]
  FVAcuts <- FVAcuts-min(FVAcuts)
  FVAcuts <- sort(FVAcuts)
  FVAval <- array(0, dim=c(3,nrow(modelOut),length(FVAcuts)))
  dimnames(FVAval)[[1]] <- c('FVAmin','FVAmax','FVAvar')
  dimnames(FVAval)[[2]] <- rownames(modelOut)
  dimnames(FVAval)[[3]] <- as.character(FVAcuts)
  
  FVAval['FVAmin',,as.character(FVAcuts)] <- as.matrix(modelOut[,names(FVAcuts)])
  FVAval['FVAmax', ,as.character(FVAcuts)] <- as.matrix(modelOut[,gsub('min','max',names(FVAcuts))])
  FVAval['FVAvar', , ] <- sapply(as.character(FVAcuts),function(i){
    (FVAval['FVAmax', ,i ]-FVAval['FVAmin', , i])/apply(FVAval[c('FVAmax','FVAmin'), ,i ],2,mean)})
  FVAL[[treatment]] <- FVAval
  rm(FVAval) 
  
}














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
  
  ### Return reaction name and stoichiometry ###
  
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
  
  if(length(rxnStoi[rxnStoi < 0]) == 0){
    substrate_prep = ""
    rxnDir <- " <- "
  }else{
    substrate_prep <- paste(sapply(c(1:length(rxnStoi[rxnStoi < 0])), function(x){
      tmp <- (rxnStoi[rxnStoi < 0] * -1)[x]
      if(tmp == 1){tmp <- ''}
      paste(tmp, speciesNames[rxnStoi < 0][x])
    }), collapse = ' + ')
  }  
  
  if(length(rxnStoi[rxnStoi > 0]) == 0){
    product_prep = ""
    rxnDir <- " -> "
  }else{
    product_prep <- paste(sapply(c(1:length(rxnStoi[rxnStoi > 0])), function(x){
      tmp <- (rxnStoi[rxnStoi > 0])[x]
      if(tmp == 1){tmp <- ''}
      paste(tmp, speciesNames[rxnStoi > 0][x])
    }), collapse = ' + ')
  }
  
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
  
  return(paste(c(rxnInfo$reaction,
      rxnInfo$genes,
      rxn_stoi,             
      paste(strsplit(rxnInfo$pathway, split = "__")[[1]], collapse = "\n")), collapse = "\n\n"))
  
  #textplot(paste(c(rxnInfo$reaction,
  #    rxnInfo$genes,
  #    rxn_stoi,             
  #    paste(strsplit(rxnInfo$pathway, split = "__")[[1]], collapse = "\n")), collapse = "\n\n")
  #    , cex = 3, valign = "top", halign = "left")
  }


cond_order_check <- function(query){
  
  require(plyr)
  
  # determine how to align a vector query sequence to match the standard naming in the reference -> first capitalize all letters and match numbers
  
  reference <- data.frame(standard = c("P", "C", "N", "L", "U"), boer = c("PO4", "GLU", "NH4", "LEU", "URA"))
  reference$standard <- factor(reference$standard)

  condSummary <- data.frame(name = query, index = 1:length(query), limitation = regmatches(query, regexpr('^[A-Za-z]', query)), DR = regmatches(query, regexpr('[0-9.]{3,4}', query)))
  condSummary$limitation <- factor(toupper(condSummary$limitation), levels = reference$standard)
  condSummary$DR <- sprintf("%.2f", as.numeric(condSummary$DR))
  
  if(!all(reference$limitation %in% reference$standard)){
    return(print("limtation names are invalid"))
  }else{
    resort <- arrange(condSummary, limitation, DR)
    resort$newOrder <- 1:length(query)
    rank(resort$index)
    }
  
}



########## Functions used in optimization of fitted flux versus actual flux ############

modelComparison <- function(reactionInfo, rxnList_form){
  
  require(dplyr)
  require(qvalue)
  
  # Evaluate nested comparisons between simpler and more complicated models
  # alternative parametrization of simple kinetics (with the same # of parameters)
  # 1 regulators vs. rMM
  # 1 regulator w/ cooperativity vs. 1 regulator & rMM
  # 2 regulators vs. each single regulator
  
  # rMech | rID | modelType | tID | regulationType | ML | npar
  # rMech is a unique ID with multiple rows if multiple regulators exist
  
  regulator_info <- lapply(rxnList_form, function(x){
    x$rxnFormData %>% dplyr::select(SubstrateID, Hill, Subtype, form = EqType) %>% filter(!(Subtype %in% c("substrate", "product"))) %>%
      mutate(rMech = x$listEntry, reaction = x$rxnID, isoenzyme_specific = ifelse(all(is.na(x$rxnFormData$enzymeInvolved)), F, T))
  })
  regulator_info <- do.call("rbind", regulator_info)
  
  reaction_info_comparison <- reactionInfo %>% tbl_df() %>% dplyr::select(rMech, reaction, ML, ncond, npar)
  
  # for each reaction, the type of reaction form (e.g. rMM, irreversible kinetics) will determine how reaction forms are compared
  
  reaction_info_comparison$modelType <- sapply(rxnList_form, function(x){
    # At some point this should be converted to a directed graph representation, to allow more sophisticated models to be posed
    # this should be done when reaction forms are originally generated
    
    if(setequal(x$rxnFormData$Subtype, c("substrate", "product"))){
      # simple kinetics
      if(x$rxnFormData$EqType[1] == "rm" & all(is.na(x$rxnFormData$enzymeInvolved))){
        "rMM" 
      }else{
        "alternative_simple_kinetics" 
      }
    }else if(length(setdiff(c("substrate", "product"), x$rxnFormData$Subtype)) != 0){
      # irreversible kinetics
      "irreversible"
    }else{
      regulators <- x$rxnFormData[!(x$rxnFormData$Subtype %in% c("substrate", "product")),]
      if(any(regulators$Hill == 0)){
        "cooperativity" 
      }else if(nrow(regulators) == 1){
        "regulator"
      }else{
        "2+ regulators"
      }
    }
  })
  
  regulator_info <- regulator_info %>% left_join(reaction_info_comparison %>% dplyr::select(rMech, modelType, ncond), by = "rMech")
  
  # for each reaction, point to its reference kinetic form(s)
  
  reaction_info_comparison$parentReaction <- sapply(1:nrow(reaction_info_comparison), function(i){
    if(reaction_info_comparison$modelType[i] == "rMM"){
      NA
    }else{
      base_model <- reaction_info_comparison %>% filter(reaction == reaction_info_comparison$reaction[i],
                                                        modelType == "rMM",
                                                        ncond == reaction_info_comparison$ncond[i])
      rMM_form <- base_model$rMech
      
      if(reaction_info_comparison$modelType[i] %in% c("alternative_simple_kinetics", "irreversible", "regulator")){
        rMM_form
      }else{
        
        # nested regulation
        evaluated_reg <- regulator_info %>% filter(rMech == reaction_info_comparison$rMech[i]) %>% dplyr::select(-rMech, -reaction)
        # all regulation of this reaction
        rxn_reg <- regulator_info %>% filter(reaction ==  reaction_info_comparison$reaction[i], ncond == reaction_info_comparison$ncond[i],
                                             !isoenzyme_specific)
        
        if(reaction_info_comparison$modelType[i] == "cooperativity"){
          if(nrow(evaluated_reg) > 1){
            warning("multiple reaction allostery not currently supported") 
          }else{
            rxn_reg <- rxn_reg %>% filter(modelType == "regulator")
            reg_parent <- inner_join(evaluated_reg %>% mutate(Hill = 1), rxn_reg, by = c("SubstrateID", "Subtype", "form", "Hill"))$rMech
          }
        }
        if(reaction_info_comparison$modelType[i] == "2+ regulators"){
          if(nrow(evaluated_reg) > 2){
            warning("3+ reactions not supported") 
          }else{
            parent_reg <- inner_join(evaluated_reg, rxn_reg %>% filter(modelType == "regulator"), by = c("SubstrateID", "Hill", "Subtype", "form"))$rMech
            if(length(parent_reg) != 2){
              warning("non-standard formatting - parental reactions not found")
            }else{
              reg_parent <- paste(parent_reg, collapse = ",") 
            }
          }
        }
        paste(c(rMM_form, reg_parent), collapse = ",")
      }
    }
  })
  
  # add some additional modelTypes
  # treat a freely inferred metabolite (w and w/o hill seperately from real metabolites)
  reaction_info_comparison <- reaction_info_comparison %>% mutate(modelType = ifelse(grepl('t_metX', rMech), paste("hypo met", modelType), modelType))
  
  # Compare full and reduced models using AICc and LRT
  model_comparison_list <- list()
  
  for(i in 1:nrow(reaction_info_comparison)){
    
    if(is.na(reaction_info_comparison$parentReaction[i])){
      # RMM 
      
    }else{
      
      daughter_rxn <- reaction_info_comparison %>% dplyr::slice(i)
      
      parent_rxn_names <- strsplit(reaction_info_comparison$parentReaction[i], split = ",")[[1]]
      parent_rxns <-  reaction_info_comparison %>% filter(rMech %in% parent_rxn_names)
      
      # determine the relative probability of eaach null versus full model based on AIC
      
      daughter_rxn <- daughter_rxn %>% mutate(AICc = 2*npar - 2*ML + 2*npar*(npar + 1)/(ncond - npar - 1))
      parent_rxns <- parent_rxns %>% mutate(AICc = 2*npar - 2*ML + 2*npar*(npar + 1)/(ncond - npar - 1),
                                            changeAIC = 1 - 1/(exp((daughter_rxn$AICc - AICc)/2) + 1)
                                            )
      
      parent_rxns <- parent_rxns %>% mutate(likDiff = daughter_rxn$ML - ML)
      
      # comparing models using LRT
      
      parent_rxns <- rbind(
        parent_rxns %>% filter(npar == daughter_rxn$npar) %>% mutate(changeP = 1/(exp(likDiff) + 1)), # Alternative parameterization with same degrees of freedom
        parent_rxns %>% filter(npar < daughter_rxn$npar) %>% mutate(changeP = 1 - pchisq(2*likDiff, daughter_rxn$npar - npar)), # Alternative model is more complex
        parent_rxns %>% filter(npar > daughter_rxn$npar) %>% mutate(changeP = 1 - pchisq(-2*likDiff, npar - daughter_rxn$npar)) # Alternative model is less complex
      )
      
      model_comparison_list[[i]] <- parent_rxns %>% dplyr::mutate(daughter_rMech = as.character(daughter_rxn$rMech), daughterModel = as.character(daughter_rxn$modelType)) %>%
        dplyr::select(reaction, daughter_rMech, daughterModel, parent_rMech = rMech, parentModel = modelType, changeAIC, changeP)
    }
  }
  
  # Multiple hypothesis correction
  # p-values are analyzed seperately for each daughterModel X parentModel
  
  model_comparison_list <- do.call("rbind", model_comparison_list)
  
  #tmp -> model_comparison_list
  
  # when looking at the effect of 2 regulators, look at the smallest improvement relative to simpler models (this is conservative)
  model_comparison_list <- rbind(
    model_comparison_list %>% filter(!(daughterModel == "2+ regulators" & parentModel == "regulator")),
    model_comparison_list %>% filter(daughterModel == "2+ regulators", parentModel == "regulator") %>% group_by(reaction, daughter_rMech, daughterModel) %>%
      dplyr::summarize(changeAIC = max(changeAIC), changeP = max(changeP)) %>% mutate(parentModel = "regulator", parent_rMech = "combined_regulators")
  )
  
  # qvalue doesn't play nicely with dplyr, so breaking the whole list into pieces and later reforming it
  # when looking at "hypothetical metabolites" (governed by principle components), use AICc for model comparison because the large number of fitted parameters
  # results in an misposed LRT
  comparison_groups <- model_comparison_list %>% group_by(daughterModel, parentModel) %>% dplyr::summarize(N = n()) %>%
    dplyr::mutate(metricUsed = ifelse(daughterModel %in% c("hypo met regulator", "hypo met cooperativity"), "changeAIC", "changeP"))
  comparison_groups$pi0 <- NA; comparison_groups$Nsig <- NA
  comparison_list <- list()
  
  for(i in 1:nrow(comparison_groups)){
    
    if(comparison_groups$N[i] < 5){next}
    
    comparison_subset <- model_comparison_list %>% filter(daughterModel == comparison_groups$daughterModel[i], parentModel == comparison_groups$parentModel[i])
    
    qvalue_object <- qvalue(comparison_subset[,comparison_groups$metricUsed[i]] %>% unlist() %>% unname(), pi0.method = "bootstrap")
    comparison_groups$pi0[i] <- qvalue_object$pi0
    comparison_subset$Qvalue <- qvalue_object$q
    
    comparison_groups$Nsig[i] <- sum(qvalue_object$q < 0.1)
    
    comparison_list[[i]] <- comparison_subset
    
  }
  
  all_models <- do.call("rbind", comparison_list)
  
  # Each reaction is assigned a consensus measure of significance - based on the least significant comparison between reaction forms
  all_models <- all_models %>% group_by(reaction, daughter_rMech, daughterModel) %>% dplyr::summarize(Qvalue = max(Qvalue)) %>%
    dplyr::select(reaction, rMech = daughter_rMech, modelType = daughterModel, Qvalue)
  
  # add back the baseline rMM reaction forms
  all_models <- rbind(all_models, reaction_info_comparison %>% filter(modelType == "rMM") %>% dplyr::select(reaction, rMech, modelType) %>% dplyr::mutate(Qvalue = NA)) %>% ungroup()
  
  signifCO <- data.frame(q_cutoff = c(0.1, 0.001, 0.00001), code = c("*", "**", "***"))
  
  all_models$signifCode <- sapply(all_models$Qvalue, function(x){
    if(is.na(x) | x > signifCO$q_cutoff[1]){
      ""
    }else{
      rev(signifCO$code[x < signifCO$q_cutoff])[1]
    }
  })
  
  # clean-up reaction naming
  
  all_models <- all_models %>% mutate(Name = NA) %>% mutate(Name = ifelse(modelType == "rMM", "Reversible michaelis-menten (default)", Name),
                                                            Name = ifelse(modelType == "alternative_simple_kinetics", "Convenience kinetics", Name),
                                                            Name = ifelse(modelType == "irreversible", "Irreversible michaelis-menten", Name))
  
  for(i in c(1:nrow(all_models))[is.na(all_models$Name)]){
    
    rxn_reg <- regulator_info %>% filter(rMech == all_models$rMech[i]) %>%
      left_join(data.frame(SubstrateID = names(rxnList_form[[all_models$rMech[i]]]$metNames), Name = unname(rxnList_form[[all_models$rMech[i]]]$metNames)), by = "SubstrateID")
    
    subnames <- apply(rxn_reg, 1, function(x){
      paste(x['Subtype'], ifelse(x['Subtype'] %in% c("cc", "mm"), "activation", "inhibition"), ifelse(x['Hill'] == 0, "(variable hill)", ""), "by", x['Name'])
    })
    
    all_models$Name[i] <- paste(subnames, collapse = " + ")
  }
  
  all_models <- rbind( 
    all_models %>% filter(!grepl('rmCond', rMech)),
    all_models %>% filter(grepl('rmCond', rMech)) %>% mutate(Name = sub('$', ' (zero flux / overflow conditions removed)', Name))
  )
  
  manualRegulators <- read.delim('./companionFiles/manual_ComplexRegulation.txt')
  
  all_models <- rbind(
    all_models %>% filter(!(rMech %in% manualRegulators$TechnicalName)),
    all_models %>% filter(rMech %in% manualRegulators$TechnicalName) %>% dplyr::select(-Name) %>% left_join(manualRegulators %>% dplyr::select(rMech = TechnicalName, Name = DisplayName) %>% unique(), by = "rMech")
  )
  
  all_models <- all_models %>% dplyr::mutate(Name = paste(Name, signifCode),
                                             Name = gsub('^ ', '', Name),
                                             Name = gsub('  ', ' ', Name))
  
  all_models$FullName <- sapply(all_models$rMech, function(a_mech){
    rxnList_form[[a_mech]]$reaction
  })
  all_models <- all_models %>% mutate(FullName = paste(FullName, Name, sep = " - "))
  
  # all_models %>% filter(is.na(Qvalue) | Qvalue < 0.1) %>% arrange(reaction) %>% View()
  return(all_models)
  
}




species_plot <- function(run_rxn, flux_fit, chemostatInfo){
  
  require(data.table)
  require(ggplot2)
  require(RColorBrewer)
  require(reshape2)
  require(grid)
  
  #generate a list of four plots:
  #1: Flux predicted from FBA ~ GR
  #2: Flux predicted from FBA and parametric form by condition
  #3: Relationship between metabolites/enzyme levels and condition
  #4: Relationship between metabolites/enzyme levels showing dynamic range
  #5: Relationship between metabolite/enzyme levels and flux carried
  #6: Flux, substrates, products, enzymes and regulators each get their own pane and are visualized together
  
  output_plots <- list()
  
  scatter_theme <- theme(text = element_text(size = 50, face = "bold"), title = element_text(size = 40, face = "bold"), panel.background = element_rect(fill = "azure"), 
      panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
      legend.text = element_text(size = 40, face = "bold"), axis.text = element_text(color = "black"))

  flux <- run_rxn$flux
  n_c <- nrow(flux)
  
  met_abund <- run_rxn$rxnSummary$rxnMet
  met_abund <- met_abund[,colSums(is.na(met_abund)) == 0, drop = F]
  
  enzyme_abund <- run_rxn$enzymes
  
  Chemoconds <- data.frame(name = factor(rownames(enzyme_abund), levels = rownames(enzyme_abund)),
                           Limitation = factor(substr(rownames(enzyme_abund),1,1), levels = unique(substr(rownames(enzyme_abund),1,1))),
                           DR = gsub('^[A-Z]', '', rownames(enzyme_abund)))
  Chemoconds$actualDR <- chemostatInfo$actualDR[chmatch(as.character(Chemoconds$name), chemostatInfo$ChemostatCond)]
    
  
  DRordering <- data.frame(DRs = as.character(c("0.05", "0.11", "0.16", "0.22", "0.30")), order = 1:5)
  Chemoconds$DRorder <- DRordering$order[chmatch(as.character(Chemoconds$DR), DRordering$DRs)]
  
  
  #### Just looking at FBA fit ####
  
  flux_plot_FBA <- data.table(METHOD = "FBA", FLUX = (flux$FVAmin + flux$FVAmax)/2, LB = flux$FVAmin, UB = flux$FVAmax, condName = Chemoconds$name, condition = Chemoconds$Limitation, DR = Chemoconds$actualDR)
  
  flux_range <- range(c(0, flux_plot_FBA$LB, flux_plot_FBA$UB))
  
  flux_plot_condLabel <- flux_plot_FBA[,list(x = DR[which.max(DR)], y = FLUX[which.max(DR)]), by = condition]
  
  output_plots$"FBA flux" <- ggplot() + geom_hline(y = 0, size = 2) + scatter_theme +
    geom_path(data = flux_plot_FBA, aes(x = DR, y = FLUX, col = condition, group = condition), size = 2) +
    geom_linerange(data = flux_plot_FBA, aes(x = DR, y = FLUX, col = condition, ymin = LB, ymax = UB), size = 2) + 
    geom_point(data = flux_plot_FBA, aes(x = DR, y = FLUX, col = condition, size = sqrt(DR)*16)) +
    geom_text(data = flux_plot_condLabel, aes(x = x, y = y, label = condition, family = "mono", fontface = "bold"), color = "black") +
    scale_size_identity() + scale_color_brewer(guide = "none", palette = "Set1") + ggtitle("Flux determined using FBA") +
    scale_y_continuous("Flux Carried", limits = flux_range)
  
  
  #### Comparision of parametric and FBA fit ####
  
  ci_theme <- theme(text = element_text(size = 30, face = "bold"), title = element_text(size = 30, face = "bold"), panel.background = element_rect(fill = "azure"), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
                    legend.text = element_text(size = 40, face = "bold"), axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90))
  
  flux_plot_FBA <- data.frame(METHOD = "FBA", FLUX = (flux$FVAmin + flux$FVAmax)/2, LB = flux$FVAmin, UB = flux$FVAmax, condName = Chemoconds$name, condition = Chemoconds$Limitation, DR = Chemoconds$DR)
  flux_plot_PAR <- data.frame(METHOD = "PAR", FLUX = flux_fit$fitted_flux$fitted, LB = flux_fit$fitted_flux$fitted - 2*flux_fit$fitted_flux$SD, UB = flux_fit$fitted_flux$fitted + 2*flux_fit$fitted_flux$SD, condName = Chemoconds$name, condition = Chemoconds$Limitation, DR = Chemoconds$DR)
  flux_plot_comp <- rbind(flux_plot_FBA, flux_plot_PAR)
  
  flux_range <- range(c(0, flux_plot_comp$LB, flux_plot_comp$UB))
  
  flux_plot_header <- data.frame(label = c("FBA-determined flux at optimum", "Parametric fit (95% CI)"), METHOD = c("FBA", "PAR"), x = 1, y = max(flux_range) - diff(flux_range)/15 * c(0,1))
  
  output_plots$"Flux comparison" <- ggplot() + geom_hline(y = 0, size = 2) + geom_pointrange(data = flux_plot_comp, aes(x = condName, y = FLUX, ymin = LB, ymax = UB, col = factor(METHOD)), size = 1.5, alpha = 0.8) +
    geom_text(data = flux_plot_header, aes(x = x, y = y, label = label, col = factor(METHOD)), size = 10, hjust = 0) +
    ci_theme + scale_color_brewer("Method", guide = "none", palette = "Set1") +
    scale_y_continuous("Relative Flux", limits = flux_range) + scale_x_discrete("Conditions")
  
  
  #### Changes in species abundance across dilution rates ###
  
  all_species <- data.frame(met_abund, enzyme_abund)
  all_species_SD <- run_rxn$specSD
  
  if('t_metX' %in% run_rxn$kineticPars$modelName){
    
    ### hypothetical metabolite overwriten with value from PCs ###
    
    mle_pars <- par_markov_chain[which.max(par_likelihood$likelihood),]
    releventPCs <- metSVD$v[rownames(metSVD$v) %in% rownames(flux),]
    npc <- ncol(releventPCs)
    
    all_species$t_metX <- t(mle_pars[run_rxn$kineticParPrior$SpeciesType == "PCL"] %*% diag(metSVD$d[1:npc]) %*% t(releventPCs))
    all_species_SD <- cbind(all_species_SD, data.frame("t_metX" = 0))
    all_species_SD <- all_species_SD[,chmatch(colnames(all_species), colnames(all_species_SD))]
    
  }
  
  if(any(colnames(all_species) != colnames(all_species_SD))){
    stop(paste(arxn, "has misaligned metabolite point-estimates and standard deviations"))
    }
  
  colnames(all_species) <- colnames(all_species_SD) <- run_rxn$all_species$commonName[chmatch(colnames(all_species), run_rxn$all_species$rel_spec)]
  
  #all_species_tab <- run_rxn$all_species[colSums(all_species != 1) != 0,]
  #all_species <- all_species[,colSums(all_species != 1) != 0]
  
  scatter_theme <- theme(text = element_text(size = 23, face = "bold", color = "black"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "right", 
                         panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan"),
                         legend.text = element_text(size = 20, face = "bold"), axis.text = element_text(color = "black"), legend.key = element_rect(fill = "white"), panel.border = element_rect(colour = "black", fill = NA))
  
  species_df <- melt(cbind(all_species, data.frame(condition = Chemoconds$Limitation, DR = Chemoconds$actualDR, FLB = flux$FVAmin, FUB = flux$FVAmax)), id.vars = c("condition", "DR", "FLB", "FUB"), value.name = "RA")
  species_df$SD <- melt(all_species_SD, value.name = "SD", measure.vars = colnames(all_species_SD))$SD
  species_df$LB <- species_df$RA - 2*species_df$SD
  species_df$UB <- species_df$RA + 2*species_df$SD
  
  scatter_theme_facetRotate <- scatter_theme + theme(strip.text.y = element_text(angle = 0))
  
  output_plots$"Metabolite and enzyme abundance" <- ggplot(species_df, aes(x = DR, y = RA, col = condition)) + geom_path(size = 2, alpha = 0.7) + geom_errorbar(aes(ymin = LB, ymax = UB), size = 1, alpha = 0.7) +
    facet_wrap(~ variable, scale = "free_y", ncol = 2) +
    scatter_theme + scale_size_identity() + scale_color_brewer("", palette = "Set1") + scale_y_continuous(expression(log[2] ~ "relative or absolute concentration")) +
    ggtitle("Metabolite and enzyme relative abundance") + scale_x_continuous("Dilution Rate (1/h)")
  
  output_plots$"Metabolite and enzyme fold changes" <- ggplot(species_df, aes(x = DR, y = RA, col = condition)) + geom_path(size = 2, alpha = 0.7) + geom_errorbar(aes(ymin = LB, ymax = UB), size = 1, alpha = 0.7) +
    facet_grid(variable ~ ., scale = "free_y", space = "free_y") +
    scatter_theme_facetRotate + scale_size_identity() + scale_color_brewer("", palette = "Set1") + scale_y_continuous(expression(log[2] ~ "relative or absolute concentration")) +
    ggtitle("Metabolite and enzyme relative abundance") + scale_x_continuous("Dilution Rate (1/h)", expand = c(0.01, 0))
  
  output_plots$"Flux ~ species" <- ggplot(species_df, aes(x = (FLB + FUB)/2, y = 2^RA, xmin = FLB, xmax = FUB, ymin = 2^LB, ymax = 2^UB, col = condition)) + geom_path(size = 2, alpha = 0.7) +
    facet_wrap( ~ variable, scale = "free_y", ncol = 2) + geom_point(aes(size = sqrt(DR)*14), alpha = 0.7) +
    geom_errorbar(size = 1, alpha = 0.7) + geom_errorbarh(size = 1, alpha = 0.7) +
    expand_limits(x = 0, y = 0) + scatter_theme + theme(axis.text.x = element_text(angle = 90)) + scale_size_identity() + scale_color_brewer("", palette = "Set1") + scale_y_continuous("Relative or absolute concentration") + 
    ggtitle("Relationship between metabolite, enzyme levels and flux carried") + scale_x_continuous("Relative flux carried")
  
  #### Combining flux carried and species into a single faceted plot ####

  species_df <- melt(cbind(all_species, condition = Chemoconds$name), id.vars = c("condition"), value.name = "RA")
  species_df$SD <- melt(all_species_SD, value.name = "SD", measure.vars = colnames(all_species_SD))$SD
  species_df <- data.table(species_df)
  species_df$metOrigin <- sapply(species_df$variable, function(x){
    if(x %in% unname(run_rxn$rxnSummary$metNames)){
      run_rxn$rxnSummary$originMet[chmatch(names(run_rxn$rxnSummary$metNames)[run_rxn$rxnSummary$metNames == x], names(run_rxn$rxnSummary$originMet))]
    }else{
     "rel" 
    }
  })
  species_df[,LB :=  2^(RA - 2*SD),]
  species_df[,UB :=  2^(RA + 2*SD),]
  species_df[,RA :=  2^RA,]
  
  # first scale absolute units by powers of 10 such that mean is between 1&10
  species_df[,units := as.character(if(all(metOrigin == "rel")){NA}else{
    meanPower <- format(mean(abs(RA)), scientific = T)
    regmatches(meanPower, regexpr('[-0-9]+$', meanPower))}),by = "variable"]
  species_df$RA[!is.na(species_df$units)] <- species_df$RA[!is.na(species_df$units)]/(10^as.numeric(species_df$units[!is.na(species_df$units)]))
  species_df$UB[!is.na(species_df$units)] <- species_df$UB[!is.na(species_df$units)]/(10^as.numeric(species_df$units[!is.na(species_df$units)]))
  species_df$LB[!is.na(species_df$units)] <- species_df$LB[!is.na(species_df$units)]/(10^as.numeric(species_df$units[!is.na(species_df$units)]))
  
  # Assign species to substrates, products, regulators, enzymes
  species_df$Pane <- NA
  species_df$Pane[species_df$variable %in% run_rxn$kineticPars$commonName[chmatch(run_rxn$rxnSummary$rxnFormData$SubstrateID[run_rxn$rxnSummary$rxnFormData$Subtype == "substrate"], run_rxn$kineticPars$modelName)]] <- "Substrate"
  species_df$Pane[species_df$variable %in% run_rxn$kineticPars$commonName[chmatch(run_rxn$rxnSummary$rxnFormData$SubstrateID[run_rxn$rxnSummary$rxnFormData$Subtype == "product"], run_rxn$kineticPars$modelName)]] <- "Product"
  species_df$Pane[species_df$variable %in% run_rxn$kineticPars$commonName[chmatch(run_rxn$rxnSummary$rxnFormData$SubstrateID[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))], run_rxn$kineticPars$modelName)]] <- "Regulator"
  species_df$Pane[species_df$variable %in% colnames(enzyme_abund)] <- "Enzyme"
  species_df <- data.table(species_df)
  
  # scale relative units by arbitrary multiplier to align them as closely as possible to absolute data
  for(a_Pane in unique(species_df$Pane)){
    species_subset <- species_df[Pane == a_Pane,,]
    if(any(species_subset$metOrigin == "abs")){
      abs_mean <- mean(species_subset$RA[species_subset$metOrigin == "abs"])
      for(relSpecies in unique(species_subset$variable[species_subset$metOrigin == "rel"])){
        species_df$UB[species_df$variable == relSpecies] <- species_df$UB[species_df$variable == relSpecies] * abs_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$LB[species_df$variable == relSpecies] <- species_df$LB[species_df$variable == relSpecies] * abs_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$RA[species_df$variable == relSpecies] <- species_df$RA[species_df$variable == relSpecies] * abs_mean/mean(species_df$RA[species_df$variable == relSpecies])
      }
    }else{
      rel_mean <- mean(species_subset$RA)
      for(relSpecies in unique(species_subset$variable)){
        species_df$UB[species_df$variable == relSpecies] <- species_df$UB[species_df$variable == relSpecies] * rel_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$LB[species_df$variable == relSpecies] <- species_df$LB[species_df$variable == relSpecies] * rel_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$RA[species_df$variable == relSpecies] <- species_df$RA[species_df$variable == relSpecies] * rel_mean/mean(species_df$RA[species_df$variable == relSpecies])
      }
    }
  }
  
  # truncate upper bounds at max(RA) * 1.2
  species_df$Truncated <- F
  for(a_Pane in unique(species_df$Pane)){
    species_subset <- species_df[Pane == a_Pane,,]
    maxRA <- max(species_subset$RA)
    if(any(species_subset$UB > maxRA*1.3)){
      species_df$Truncated[species_df$UB > maxRA*1.3 & species_df$Pane == a_Pane] <- T
      species_df$UB[species_df$UB > maxRA*1.3 & species_df$Pane == a_Pane] <-  maxRA*1.3
      }
    }
  
  flux_plot <- data.table(melt(flux_plot_comp, id.vars = c("condName", "FLUX", "LB", "UB"), measure.vars = "METHOD"))
  flux_plot <- flux_plot[,-which(colnames(flux_plot) == "variable"), with = F]
  flux_plot$variable <- as.character("Flux")
  flux_plot$Truncated <- F
  setnames(flux_plot, old = c("condName", "FLUX", "value", "variable"), new = c("condition", "RA", "variable", "Pane"))
  
  all_changes <- rbind(species_df[,chmatch(colnames(flux_plot), colnames(species_df)), with = F], flux_plot)
  all_changes$Pane <- factor(all_changes$Pane, levels = c("Substrate", "Product", "Enzyme", "Regulator", "Flux"))
  
  all_changes$condition <- factor(all_changes$condition, levels = rownames(met_abund))
  all_changes$Limitation <- factor(substr(all_changes$condition, 1, 1))
  all_changes$pathGroup <- factor(paste(all_changes$Limitation, all_changes$variable, sep = "_"))
  
  # Define species color
  
  variable_colors <- data.frame(unique(as.data.frame(all_changes)[,c("variable", "Pane")]), color = NA)
  
  if(sum(variable_colors$Pane %in% c("Substrate", "Product")) > 8){
    variable_colors$color[variable_colors$Pane %in% c("Substrate", "Product")] <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel2")[c(1:(sum(variable_colors$Pane %in% c("Substrate", "Product")) - 8))])
  }else{
    variable_colors$color[variable_colors$Pane %in% c("Substrate", "Product")] <- brewer.pal(8, "Dark2")[1:sum(variable_colors$Pane %in% c("Substrate", "Product"))]
    }
  variable_colors$color[variable_colors$Pane == "Flux"] <- brewer.pal(8, "Set1")[1:2]
  variable_colors$color[variable_colors$Pane == "Regulator"] <- rev(brewer.pal(9, "Reds"))[c(3,5,7,4,6)][1:sum(variable_colors$Pane == "Regulator")]
  variable_colors$color[variable_colors$Pane == "Enzyme"] <- rev(brewer.pal(9, "YlGnBu"))[c(1,3,5,7,2,4,6,8,9)][1:sum(variable_colors$Pane == "Enzyme")]
  
  variable_colors <- variable_colors[chmatch(levels(all_changes$variable), as.character(variable_colors$variable)),]
  
  # Define species labels
  
  variable_labels <- data.frame(unique(as.data.frame(all_changes)[,c("variable", "Pane")]), units = NA, label = NA)
  variable_labels$units <- species_df$units[chmatch(as.character(variable_labels$variable), as.character(species_df$variable))]
  
  
  if(!all(is.na(variable_labels$units))){
    abs_units <- data.frame(units = as.numeric(variable_labels$units[!is.na(variable_labels$units)]), referenceSymbol = NA, referencePower = NA, scaled =  NA)
    
    abs_units[abs_units$units >= 0, 2] <- "M"; abs_units[abs_units$units >= 0, 3] <- 0
    abs_units[data.table::between(abs_units$units, -3, -1, incbounds = T), 2] <- "mM"; abs_units[data.table::between(abs_units$units, -3, -1, incbounds = T), 3] <- -3
    abs_units[data.table::between(abs_units$units, -6, -4, incbounds = T), 2] <- "uM"; abs_units[data.table::between(abs_units$units, -6, -4, incbounds = T), 3] <- -6
    abs_units[abs_units$units < -7, 2] <- "nM"; abs_units[abs_units$units < -7, 3] <- -9
    
    abs_units$scaled <- 10^(abs_units$units - as.numeric(abs_units$referencePower))
    abs_units$scaled[abs_units$scaled == 1] <- ""
    
    variable_labels$units[!is.na(variable_labels$units)] <- paste0(abs_units$scaled, abs_units$referenceSymbol)
    }
  variable_labels$units[is.na(variable_labels$units)] <- "a.u."
  
  variable_labels$label <- paste0(variable_labels$variable, ' (', variable_labels$units, ')')
  variable_labels$label[variable_labels$variable == "FBA"] <- "Observed flux"
  variable_labels$label[variable_labels$variable == "PAR"] <- "Michaelis-Menten flux (95% CI)"
  tmp <- variable_labels$y[variable_labels$variable == "FBA"]
  variable_labels$y[variable_labels$variable == "FBA"] <- variable_labels$y[variable_labels$variable == "PAR"]
  variable_labels$y[variable_labels$variable == "PAR"] <- tmp
  
  variable_labels$x = "P0.05"
  variable_labels$y = NA
  
  for(a_Pane in unique(variable_labels$Pane)){
    PaneMaxVal <- max(all_changes$UB[all_changes$Pane == a_Pane])
    variable_labels$y[variable_labels$Pane == a_Pane] <- PaneMaxVal - c(1:sum(variable_labels$Pane == a_Pane))*PaneMaxVal*0.12
    }
  
  ci_theme <- theme(text = element_text(size = 30, face = "bold"), title = element_text(size = 30, face = "bold"), panel.background = element_rect(fill = "azure"), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
                    axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90),
                    panel.margin = unit(1, "lines"))
  
  output_plots$"Flux and species" <- ggplot() + geom_path(data = all_changes, aes(x = condition, y = RA, col = variable, group = pathGroup), size = 2, alpha = 1) +
    geom_pointrange(data = all_changes, aes(x = condition, y = RA, ymin = LB, ymax = UB, col = factor(variable)), size = 1.2, alpha = 0.5) +
    geom_point(data = all_changes, aes(x = condition, y = RA, col = factor(variable)), size = 5) +
    geom_text(data = all_changes[all_changes$Truncated,], aes(x = condition, y = UB, col = variable), label = "+", size = 8, alpha = 0.5) +
    geom_text(data = variable_labels, aes(x = x, y = y, col = variable, label = label), hjust = 0) +
    geom_vline(x = 0) + geom_hline(y = 0) +
    facet_grid(Pane ~ ., scale = "free_y") +
    scale_color_manual(guide = "none", labels = variable_colors$variable, values = variable_colors$color) +
    ci_theme + scale_size_identity() +
    scale_y_continuous("Concentration and Flux", expand = c(0,0)) + scale_x_discrete("Conditions") + expand_limits(y = 0)
  
  
  return(output_plots)
  }


likViolin <- function(par_likelihood, markovPars){
  
  require(ggplot2)
  
  ### show the distribution of likelihood values for each markov chain - similarity of their shapes is a good indication of convergence to the posterior
  ### strong deviations should be corrected by decreasing autocorrelation (by decreasing how often a markov sample is saved) or increasing the lenght of burnin or number of samples
  
  boxplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 20, angle = 60, vjust = 0.6), axis.line = element_blank(), strip.background = element_rect(fill = "darkseagreen2"),
                         strip.text = element_text(size = 25, colour = "darkblue"))
  
  likRange <- range(par_likelihood$likelihood)
  parList <- data.frame(x = 2, y = likRange[1] + c(likRange[2] - likRange[1])/10 * c(0:2), parameter = unname(mapply(function(x,y){paste(x, y, sep = ": ")}  , x = names(unlist(markovPars)), y = unname(unlist(markovPars)))))
  
  ggplot() + geom_violin(data = par_likelihood, aes(x = factor(index), y = likelihood), fill = "darkblue") + geom_text(data = parList, aes(x = x, y = y, label = parameter)) +
    boxplot_theme + scale_x_discrete("Chain Number") + scale_y_continuous("Log-likelihood")
  
}


hypoMetTrend <- function(run_rxn, metSVD, tab_boer){
  
  ##### Looking at a hypothetical "omtimal" regulator #####
  # For each reaction evaluate a hypothetical activator and inhibitor which best explains flux
  # hypothetical regulators variation is governed by PCs
  # U is inferred based on draws from gaussian densities
  # Dt(V) is the principle components of the metabolomics matrix
  
  # Return:
  # Plots:
  # - Trace of hypothetical regulators relative abundance across conditions
  # - Metabolites which look most like this trend
  # - Distribution of hill coefficients (if this is a variable)
  # Summary:
  # - Correlation of each metabolite to median RA of hypothetical metabolite  
  
  #### Generate informative plots showing the metabolic trends which inhibitors or activators of a reaction may follow ####
  
  require(reshape2)
  require(ggplot2)
  
  boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 90, hjust = 1), axis.line = element_blank(),
                         axis.text = element_text(color = "black"))
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), axis.text.x = element_text(size = 12, angle = 90), axis.line = element_blank()) 
  
  scatter_theme <- theme(text = element_text(size = 23, face = "bold", color = "black"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "right", 
                         panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan"),
                         legend.text = element_text(size = 20, face = "bold"), axis.text = element_text(color = "black"), legend.key = element_rect(fill = "white"))
  
  
  hypoMetPlots <- list()
  
  ### Calculate distribution of interpretable measures ###
  
  releventPCs <- metSVD$v[rownames(metSVD$v) %in% rownames(run_rxn$rxnSummary$flux),]
  npc <- ncol(releventPCs)
  
  reconstructedDraws <- par_markov_chain[,run_rxn$kineticParPrior$SpeciesType == "PCL"] %*% diag(metSVD$d[1:npc]) %*% t(releventPCs)
  reconstructedDraws <- reconstructedDraws - apply(reconstructedDraws, 1, mean)
  
  ### Create a violin plot, with the MLE highlighted ###
  
  mod_type <- run_rxn$rxnSummary$rxnFormData$Type[run_rxn$rxnSummary$rxnFormData$SubstrateID == "t_metX"]
  mod_type <- ifelse(mod_type == "act", "allosteric activator", "allosteric inhibitor")
  
  PCconds <- data.frame(name = factor(colnames(reconstructedDraws), levels = colnames(reconstructedDraws)), Limitation = factor(substr(colnames(reconstructedDraws),1,1)), DR = gsub('^[A-Z]', '', colnames(reconstructedDraws)))
  
  PCreco_melt <- melt(cbind(PCconds, t(reconstructedDraws)), id.vars = colnames(PCconds))
  
  ### the MLE ###
  
  mle_pars <- par_markov_chain[which.max(par_likelihood$likelihood),]
  mle_RA <- mle_pars[run_rxn$kineticParPrior$SpeciesType == "PCL"] %*% diag(metSVD$d[1:npc]) %*% t(releventPCs)
  mle_RA <- mle_RA - mean(mle_RA)
  
  mle_info <- data.frame(PCconds, value = c(mle_RA))
  
  hypoMetPlots$"Hypothetical Regulator Trend" <- ggplot() + geom_violin(data = PCreco_melt, aes(x = name, y = value, fill = factor(Limitation))) + scale_y_continuous(expression(log[2]~"Relative Concentration")) +
    geom_point(data = mle_info, aes(x = name, y = value), size = 3, shape = 21, color = "BLACK", fill = "RED") + scale_fill_brewer(guide = "none", palette = "Pastel1") +
    ggtitle("Hypothetical Regulator Trend") + boxplot_theme + scale_x_discrete("Nutrient Condition")
  
  
  ### Show distribution of hill coefficients ###
  
  if(any(run_rxn$kineticPars$SpeciesType == "hillCoefficient")){
    
    hillDF <- data.frame(hill = par_markov_chain[,run_rxn$kineticParPrior$rel_spec == "t_metX" & run_rxn$kineticParPrior$SpeciesType == "hillCoefficient"])
    
    alloPrior <- run_rxn$kineticParPrior[run_rxn$kineticParPrior$rel_spec == "t_metX" & run_rxn$kineticParPrior$SpeciesType == "hillCoefficient",]
    allo_binwidth <- 0.05
    nrand <- 10000
    SpSldensity <- data.frame(hill = c(rnorm(nrand*(1-alloPrior$par_3), alloPrior$par_1, alloPrior$par_2), rep(0, nrand*alloPrior$par_3))) # generate the expected null from the spike and slab prior
    
    hypoMetPlots$Hill <- ggplot() + stat_bin(data = SpSldensity, aes(x = 2^hill, y = ..density.., geom="line", position="identity"), binwidth = 0.05, fill="blue", alpha = 0.7) +
      stat_bin(data = hillDF, aes(x = 2^hill, y = ..density.., geom="line", position="identity"), binwidth = 0.05, fill = "RED", alpha = 0.7) +
      geom_vline(x = 2^mle_pars[run_rxn$kineticParPrior$rel_spec == "t_metX" & run_rxn$kineticParPrior$SpeciesType == "hillCoefficient"], size = 2, alpha = 0.5) +
      geom_text(aes(x = 0, y = 9.8, label = 'Posterior distribution', color = "RED"), hjust = 0, size = 8) +
      geom_text(aes(x = 0, y = 9, label = 'Prior / expectation', color = "BLUE"), hjust = 0, size = 8) +
      scale_x_continuous("Hill coefficient", limits = c(0, max(2^hillDF$hill))) + scale_y_continuous("Density", expand = c(0,0)) +
      barplot_theme + ggtitle("Hill coefficient posterior samples") + scale_color_identity() + scale_fill_identity() + scale_size_identity()
  }
  
  ### What measured metabolites most resemble these trends ###
  
  specCorrQuantiles <- t(apply(cor(t(tab_boer[,colnames(tab_boer) %in% colnames(reconstructedDraws)]), t(reconstructedDraws)), 1, function(x){quantile(x, probs = c(0.025, 0.5, 0.975))}))
  
  # saving all correlations
  specCorrQuantiles_all <- cbind(rMech = run_rxn$rxnSummary$listEntry, tab_boer[,colnames(tab_boer) %in% c("SpeciesName", "SpeciesType")], specCorrQuantiles)
  rownames(specCorrQuantiles_all) <- NULL
  hypoMetPlots$specCorrQuantiles_all <- specCorrQuantiles_all
  
  # only looking at top hits
  specCorrQuantiles <- specCorrQuantiles[order(specCorrQuantiles[,2], decreasing = T),][1:4,]
  
  specCorrQuantiles_m <- melt(t(specCorrQuantiles)); colnames(specCorrQuantiles_m) <- c("Quantile", "variable", "Correlation")
  specCorrQuantiles_m$y <- rep(2:4, times = nrow(specCorrQuantiles))
  specCorrQuantiles_m$text <- mapply(function(x, y){paste(x,y, sep = ": ")}, x = specCorrQuantiles_m$Quantile, y = round(specCorrQuantiles_m$Correlation,2))
  
  specCorrQuantiles_m <- rbind(specCorrQuantiles_m, data.frame(Quantile = NA, variable = rownames(specCorrQuantiles), Correlation = NA, y = 1, text = "Correlation quantiles"))
  
  PCsimilaritySubset <- tab_boer[rownames(tab_boer) %in% rownames(specCorrQuantiles),colnames(tab_boer) %in% colnames(reconstructedDraws)]
  PCsimilaritySubset <- PCsimilaritySubset[chmatch(rownames(PCsimilaritySubset), rownames(specCorrQuantiles)),]
  subsetRange <- t(apply(PCsimilaritySubset, 1, range))
  
  specCorrQuantiles_m$yadj <- mapply(function(x,y){subsetRange[rownames(subsetRange) == x][2] - diff(subsetRange[rownames(subsetRange) == x])/10*y}, x = specCorrQuantiles_m$variable, y = specCorrQuantiles_m$y)
  specCorrQuantiles_m$variable <- factor(specCorrQuantiles_m$variable, levels = rownames(specCorrQuantiles))
  
  similarGeneMelt <- data.table(melt(cbind(PCconds, t(PCsimilaritySubset)), id.vars = colnames(PCconds)))
  similarGeneMelt$Limitation <- factor(similarGeneMelt$Limitation, levels = unique(chemostatInfo$Limitation))
  similarGeneMelt$variable <- factor(similarGeneMelt$variable, levels = rownames(specCorrQuantiles))
  
  hypoMetPlots$candidateMetabolites <- ggplot() + geom_path(data = similarGeneMelt, aes(x = DR, y = value, col = Limitation, size = 2, group = Limitation)) + facet_wrap(~ variable, scale = "free_y", ncol = 2) +
    geom_text(data = specCorrQuantiles_m, aes(x = 1, y = yadj, label = text), hjust = 0) +
    scatter_theme + scale_size_identity() + scale_color_brewer("", palette = "Set1") + scale_y_continuous("Relative concentration") +
    ggtitle("Most highly correlated metabolite to predicted regulator")
              
  return(hypoMetPlots)
  
}


hillPlot <- function(run_rxn){
  
  #### Generate a histogram of hill coefficients for an allosteric regulator ####
  
  require(reshape2)
  require(ggplot2)
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), axis.text.x = element_text(size = 12, angle = 90), axis.line = element_blank()) 
  
  # a hill coefficient of 0 is a flag for a parameter
  mod_type <- run_rxn$rxnSummary$rxnFormData$Type[run_rxn$rxnSummary$rxnFormData$Hill == 0]
  mod_type <- ifelse(mod_type == "act", "allosteric activator", "allosteric inhibitor")
  
  ### the MLE ###
  
  hillDF <- data.frame(hill = par_markov_chain[,run_rxn$kineticParPrior$SpeciesType == "hillCoefficient"])
  hill_MLE <- hillDF$hill[which.max(par_likelihood$likelihood)]
  
  alloPrior <- run_rxn$kineticParPrior[run_rxn$kineticParPrior$SpeciesType == "hillCoefficient",]
  allo_binwidth <- 0.05
  nrand <- 10000
  SpSldensity <- data.frame(hill = c(rnorm(nrand*(1-alloPrior$par_3), alloPrior$par_1, alloPrior$par_2), rep(0, nrand*alloPrior$par_3)))
  
  hillPlot <- ggplot() + stat_bin(data = SpSldensity, aes(x = 2^hill, y = ..density.., geom="line", position="identity"), binwidth = 0.05, fill="blue", alpha = 0.5) +
    stat_bin(data = hillDF, aes(x = 2^hill, y = ..density.., geom="line", position="identity"), binwidth = 0.05, fill = "RED", alpha = 0.5) +
    geom_vline(x = 2^hill_MLE, size = 2, alpha = 0.5) +
    geom_text(aes(x = 0, y = 9.8, label = 'Posterior distribution', color = "RED"), hjust = 0, size = 8) +
    geom_text(aes(x = 0, y = 9, label = 'Prior / expectation', color = "BLUE"), hjust = 0, size = 8) +
    scale_x_continuous("Hill coefficient", limits = c(0, max(2^hillDF$hill))) + scale_y_continuous("Density", expand = c(0,0)) +
    barplot_theme + ggtitle("Hill coefficient posterior samples") + scale_color_identity() + scale_fill_identity() + scale_size_identity()
  
  return(hillPlot)
  
}



flux_fitting <- function(run_rxn, par_markov_chain, par_likelihood){
  
  ### Look at the MLE parameter set - the posterior sample with the highest likelihood ###
  ######## Determine overall fit - how well flux(par) matches FBA-determined flux #########
  #### Optimization to determine flux and sd(flux) is carried out as in FluxOptimClus.R ##
  
  require(data.table)
  require(nnls)
  
  n_c <- nrow(run_rxn$metabolites )
  
  # predict flux based upon parameter sets to determine how much variance in flux can be accounted for using the prediction
  param_summary <- list()
  
  param_interval <- 2^(apply(par_markov_chain, 2, function(x){quantile(x, probs = c(0.025, 0.975))}))
  param_interval <- data.frame(t(param_interval), median = 2^(apply(par_markov_chain, 2, median)), MLE = 2^(par_markov_chain[which.max(par_likelihood$likelihood),]), row.names = NULL)
  
  param_interval$absoluteQuant <- ifelse(run_rxn$kineticParPrior$rel_spec %in% names(run_rxn$rxnSummary$originMet)[run_rxn$rxnSummary$originMet == "abs"] & run_rxn$kineticParPrior$SpeciesType == "Metabolite", T, F)
  # store logQ in the keq absolute quant position - meaningful in relation to keq
  param_interval$absoluteQuant[run_rxn$kineticParPrior$SpeciesType == "keq"] <- as.character((run_rxn$kineticParPrior$par_2[run_rxn$kineticParPrior$SpeciesType == "keq"] + run_rxn$kineticParPrior$par_1[run_rxn$kineticParPrior$SpeciesType == "keq"])/2)
  
  # Save posterior distribution summaries and meaning of parameters
  param_summary$param_interval <- param_interval
  param_summary$kineticParPrior <- run_rxn$kineticParPrior
  
  param_summary$param_species <- run_rxn$rxnSummary$rxnFormData
  median_species_abund <- apply(run_rxn$metabolites, 2, median)
  param_summary$param_species$medianAbund <- median_species_abund[chmatch(param_summary$param_species$SubstrateID, names(median_species_abund))]
  
  
  ### Reproduce the flux and sd(flux) of the most likely parameter set ###
  
  flux <- run_rxn$flux
  met_abund <- run_rxn$metabolites
  enzyme_abund <- run_rxn$enzymes
  
  mle_pars <- par_markov_chain[which.max(par_likelihood$likelihood),]
  
  if('t_metX' %in% run_rxn$kineticPars$modelName){
  
    ### hypothetical metabolite overwriten with value from PCs ###
    releventPCs <- metSVD$v[rownames(metSVD$v) %in% rownames(flux),]
    npc <- ncol(releventPCs)
  
    if(npc != sum(run_rxn$kineticParPrior$SpeciesType == "PCL")){
      print(paste("Inconsistent number of principal components for rx", arxn, "likely version problem -> rerun parameter sets"))
      return(NA)
      }
    
    met_abund$t_metX <- t(mle_pars[run_rxn$kineticParPrior$SpeciesType == "PCL"] %*% diag(metSVD$d[1:npc]) %*% t(releventPCs))
    
  }
  
  mle_pars <- mle_pars[run_rxn$kineticPars$SpeciesType != "PCL"]
  par_stack <- rep(1, n_c) %*% t(2^(mle_pars)); colnames(par_stack) <- run_rxn$kineticPars$formulaName[run_rxn$kineticPars$SpeciesType != "PCL"]
  
  occupancy_vals <- data.frame(met_abund, par_stack)
  
  #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  #occupany of enzymes * relative abundance of enzymes
  kinetically_differing_isoenzymes <- any(names(run_rxn$occupancyEq[["l_occupancyExpression"]]) %in% rownames(run_rxn$rxnSummary$enzymeComplexes))
  if(!(kinetically_differing_isoenzymes)){
    predOcc <- eval(run_rxn$occupancyEq[["l_occupancyExpression"]], occupancy_vals) #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
    enzyme_activity <- (predOcc %*% t(rep(1, ncol(enzyme_abund))))*2^enzyme_abund #occupany of enzymes * relative abundance of enzymes
  }else{
    occEqtn_complex_match <- data.frame(complex = colnames(enzyme_abund), occEqtn = NA)
    occEqtn_complex_match$occEqtn[occEqtn_complex_match$complex %in% names(run_rxn$occupancyEq[["l_occupancyExpression"]])] <- occEqtn_complex_match$complex[occEqtn_complex_match$complex %in% names(run_rxn$occupancyEq[["l_occupancyExpression"]])]
    occEqtn_complex_match$occEqtn[is.na(occEqtn_complex_match$occEqtn)] <- "other"
    
    enzyme_activity <- NULL
    for(isoenzyme in names(run_rxn$occupancyEq[["l_occupancyExpression"]])){
      predOcc <- eval(run_rxn$occupancyEq[["l_occupancyExpression"]][[isoenzyme]], occupancy_vals)
      enzyme_activity <- cbind(enzyme_activity, predOcc %*% t(rep(1, sum(occEqtn_complex_match$occEqtn == isoenzyme)))*2^enzyme_abund[,colnames(enzyme_abund) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]])
    }
  }
  
  # fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  # flux objective is set as the average of the minimal and maximal allowable flux flowing through the reaction at the optimal solution
  
  flux_fit <- nnls(enzyme_activity, (flux$FVAmax + flux$FVAmin)/2) 
  
  # determine the sd of the fitted measures
  
  nnlsCoef <- t(t(rep(1, n_c)))  %*% flux_fit$x; colnames(nnlsCoef) <- run_rxn$all_species$formulaName[run_rxn$all_species$SpeciesType == "Enzyme"]
  
  all_components <- data.frame(occupancy_vals, enzyme_abund, nnlsCoef)
  
  # partial derivatives of each measured specie in a condition
  
  measured_partials <- run_rxn$occupancyEq$kinetic_form_partials[names(run_rxn$occupancyEq$kinetic_form_partials) %in% c(run_rxn$kineticPars$modelName[run_rxn$kineticPars$measured & !is.na(run_rxn$kineticPars$measured)], colnames(run_rxn$enzymes))]
  
  comp_partials <- matrix(NA, nrow = n_c, ncol = length(measured_partials))
  colnames(comp_partials) <- names(measured_partials)
  
  for(j in 1:ncol(comp_partials)){
    comp_partials[,j] <- with(all_components, eval(run_rxn$occupancyEq$kinetic_form_partials[[j]]))
  }
  
  # calculate the fitted standard deviation after first finding the by-condition residual covariance matrix
  
  flux_SD <- rep(NA, n_c)
  for(i in 1:n_c){
    sampleCov <- run_rxn$specCorr * t(t(run_rxn$specSD[i,])) %*% run_rxn$specSD[i,]
    flux_SD[i] <- sqrt(t(comp_partials[i,]) %*% sampleCov %*% t(t(comp_partials[i,])))
  }
  
  
  ##### Comparison of flux determined through constraint-based modeling and parameteric fitting ######
  
  
  fit_summary <- data.frame(residDF = nrow(run_rxn$flux) - ncol(par_markov_chain), parametricFit = NA, NNLS = NA, LS = NA, LS_met = NA, LS_enzyme = NA, TSS = NA)
  
  v_metab <- colnames(met_abund)[colSums(met_abund == 0) != n_c] # coerce metabolite abundances to a matrix
  met_abund <- as.matrix(met_abund[,colSums(met_abund == 0) != n_c]) 
  colnames(met_abund) <- v_metab
  
  avg_flux <- (flux$FVAmax + flux$FVAmin)/2
  
  ### using LS regression, how much variance is explained.  Not really a fair comparison
  fit_summary$LS_met = anova(lm(avg_flux ~ met_abund))$S[1]
  fit_summary$LS_enzyme = anova(lm(avg_flux ~ enzyme_abund))$S[1]
  fit_summary$LS = sum(anova(lm(avg_flux ~ met_abund + enzyme_abund))$S[1:2])
  fit_summary$TSS = sum(anova(lm(avg_flux ~ met_abund + enzyme_abund))$S)
  
  ### using flux fitted from the MLE parameter set, how much variance is explained
  
  fit_summary$parametricFit <- fit_summary$TSS-sum((flux_fit$residuals)^2)
  
  ### using flux fitted using non-negative least squares regression, how much variance is explained ### metabolite abundances are corrected for whether the metabolite is a product (*-1) or reactant (*1)
  NNLSmetab <- 2^met_abund * -1*(rep(1, n_c) %*% t(unname(run_rxn$rxnSummary$rxnStoi)[chmatch(colnames(met_abund), names(run_rxn$rxnSummary$rxnStoi))]))
  
  ### add potential activators and inhibitors
  if(length(run_rxn$rxnSummary$rxnFormData$Subtype[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))]) != 0){
    modifiers = data.frame(modifier = run_rxn$rxnSummary$rxnFormData$SubstrateID[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))], direction = ifelse(run_rxn$rxnSummary$rxnFormData$Type[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))] == "act", 1, -1))
    modifier_effect <- (2^run_rxn$metabolites)[,chmatch(modifiers$modifier, colnames(run_rxn$metabolites)),drop = F] * t(t(rep(1,n_c))) %*% modifiers$direction
    colnames(modifier_effect) <- paste(modifiers$modifier, "mod", sep = "")
    NNLSmetab <- cbind(NNLSmetab, modifier_effect)
  }
  NNLSmetab <- NNLSmetab[,apply(NNLSmetab, 2, function(x){var(x) != 0})]
  
  if(all(avg_flux < 0)){
    
    nnls_fit <- nnls(as.matrix(data.frame(NNLSmetab, 2^enzyme_abund)), -1*avg_flux)
    tpred <- nnls_fit$residuals
    fit_summary$NNLS <- fit_summary$TSS-sum((tpred)^2)
  
    fit_summary$nnlsPearson <- c(cor(-1*nnls_fit$fitted, avg_flux, method = "pearson"))
    fit_summary$nnlsSpearman <- c(cor(-1*nnls_fit$fitted, avg_flux, method = "spearman"))
    
  }else{
    
    nnls_fit <- nnls(as.matrix(data.frame(NNLSmetab, 2^enzyme_abund)), avg_flux)
    tpred <- nnls_fit$residuals
    fit_summary$NNLS <- fit_summary$TSS-sum((tpred)^2)
    
    fit_summary$nnlsPearson <- c(cor(nnls_fit$fitted, avg_flux, method = "pearson"))
    fit_summary$nnlsSpearman <- c(cor(nnls_fit$fitted, avg_flux, method = "spearman"))
  
  }
   
  ### correlations
  fit_summary$parPearson <- c(cor(flux_fit$fitted, avg_flux, method = "pearson"))
  fit_summary$parSpearman <- c(cor(flux_fit$fitted, avg_flux, method = "spearman"))
  
  # correct pearson correlation for measurement noise propagated through reaction eqtn # assume accuracy of flux for now
  # for simplicity currently assume fixed noise - equal to the median SD - ideally this should be replaced with a weighted measure which
  # directly accounts for the defined differences in propagated experimental noise == flux_SD
  corrected_pearson <- fit_summary$parPearson / sqrt(var(avg_flux)/(var(avg_flux) + median(flux_SD)^2))
  fit_summary$parPearson_c <- ifelse(corrected_pearson > 0, min(corrected_pearson, 1), max(corrected_pearson, -1))
    
  output <- list()
  output$fit_summary <- fit_summary
  output$param_interval <- param_summary
  output$fitted_flux <- data.frame(fitted = flux_fit$fitted, SD = flux_SD)
  output
}


modelComparisonPlots <- function(rxn, reactionInfo, all_reactionInfo){
  
  require(ggplot2)
  require(dplyr)
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 20, color = "black", hjust = 1), axis.line = element_blank(),
                         strip.background = element_rect(fill = "cornflowerblue"), axis.ticks.y = element_line(size = 1, color = "black"))
  
  
  # Look at alternative models for a single reaction 
  rxInfo_subset <- reactionInfo %>% filter(reaction %in% rxn) %>% filter(is.na(Qvalue) | Qvalue < 0.1)
  
  # pull in fit information
  rxInfo_subset <- rxInfo_subset %>% left_join(rxn_fits %>% dplyr::select(rMech = rxn, parSpearman), by = "rMech") %>%
    left_join(fraction_flux_deviation %>% dplyr::select(rMech = rxn, dotProduct, angle, get("Interval Overlap")), by = "rMech") %>%
    mutate(modelType = ifelse(modelType %in% c("2+ regulators", "cooperativity"), "coop/2+reg", modelType))
  
  maskedReg <- all_reactionInfo %>% filter(reaction == rxn) %>% filter(!(rMech %in% rxInfo_subset$rMech)) %>%
    filter(!(modelType %in% c("alternative_simple_kinetics", "irreversible", "hypo met regulator", "hypo met cooperativity"))) %>%
    mutate(modelType = ifelse(modelType %in% c("2+ regulators", "cooperativity"), "coop/2+reg", modelType)) %>%
    group_by(reaction, modelType, ncond) %>% dplyr::summarize(N = n())
  
  if(nrow(maskedReg) == 0){
    neglectedNull <- expand.grid(modelType = c("regulator", "coop/2+reg"),  ncond = unique(rxInfo_subset$ncond), stringsAsFactors = F)  
  }else{
    neglectedNull <- expand.grid(modelType = c("regulator", "coop/2+reg"),  ncond = unique(rxInfo_subset$ncond), stringsAsFactors = F) %>% anti_join(maskedReg %>% dplyr::select(modelType, ncond), by = c("modelType", "ncond"))
  }
  
  if(nrow(neglectedNull) != 0){
    maskedReg <- rbind(maskedReg, 
                       data.frame(neglectedNull, reaction = rxn, N = 0))
  }
  
  barplot_theme <- barplot_theme + theme(axis.text.x = element_text(size = min(8 + 80/(nrow(rxInfo_subset)+nrow(maskedReg)), 24), color = "black", vjust = 0.5, hjust = 1, angle = 90))
  
  pathwayPlot_list <- list()
  
  for(plotType in c("SpCorr", "DotProduct")){
    
    #### Define plotting variables ####
    if(plotType == "SpCorr"){
      plotDat <- rxInfo_subset %>% dplyr::select(Name, modelType, metric = parSpearman, ncond)
      good_dir = "+"
    }
    if(plotType == "DotProduct"){
      plotDat <- rxInfo_subset %>% dplyr::select(Name, modelType, metric = dotProduct, ncond)
      good_dir = "+"
    }
    #if(plotType == "Angle"){
    #  plotDat <- rxInfo_subset %>% dplyr::select(Name, modelType, metric = angle, ncond)
    #  good_dir = "-"
    #}
    #if(plotType == "Capture"){
    #  plotDat <- rxInfo_subset %>% dplyr::select(Name, modelType, metric = get("Interval Overlap"), ncond)
    #  good_dir = "+"
    #}
    
    #### Order reactions first by lowest/highest value for a rxn and then within a reaction ####
    
    nullRegulation <- maskedReg %>% ungroup() %>% mutate(Name = paste("Not improved by", modelType)) %>% left_join(plotDat %>% filter(modelType == "rMM") %>% dplyr::select(ncond, metric), by = "ncond") %>%
      ungroup() %>% dplyr::mutate(modelType = paste("null", modelType))
    
    plotDat <- rbind(plotDat, nullRegulation %>% dplyr::select(-c(reaction, N)))
    
    classAesthetics <- data.frame(modelType = c("irreversible", "rMM", "alternative_simple_kinetics", "null regulator", "regulator", "null coop/2+reg",
                                                "coop/2+reg", "hypo met regulator", "hypo met cooperativity"),
                                  fill = c("cornflowerblue", "dodgerblue3", "blueviolet", "white", "firebrick1", "white", "limegreen", "darkgoldenrod1", "darkgoldenrod3"),
                                  color = c("black", "black", "black", "firebrick1", "black", "limegreen", "black", "black", "black"))
    
    plotDat <- plotDat %>% left_join(classAesthetics, by = "modelType") %>%
      mutate(modelType = factor(modelType, levels = classAesthetics$modelType))
    
    
    if(good_dir == "+"){
      
      plotDat <- plotDat %>% group_by(modelType) %>% arrange(metric)
      
    }else{
      
      plotDat <- plotDat %>% group_by(modelType) %>% arrange(desc(metric))
      
    }
    
    plotDat <- plotDat %>% ungroup() %>% mutate(Name = ifelse(modelType == "coop/2+reg", gsub('\\+', '+\n', Name), Name)) %>%
      mutate(Name = factor(Name, levels = unique(Name)))
    
    if(plotType == "SpCorr"){
      
      compPlot <- ggplot(data = plotDat, aes(x = Name, y = metric)) +
        geom_bar(aes(color = color, fill = fill), stat = "identity", position = "dodge", width = 0.75) + 
        barplot_theme + scale_x_discrete(name = "Reactions") +
        scale_y_continuous(expression("Spearman correlation between "~ V^FBA ~ "&" ~ V^PAR), expand = c(0,0)) + expand_limits(y = 1) +
        geom_text(data = nullRegulation, aes(x = Name, y = metric/2, label = N)) +
        scale_color_identity() + scale_fill_identity() +
        geom_hline(yintercept = -Inf, size = 2) +
        geom_vline(xintercept = -Inf, size = 2) + scale_size_identity()
        
    }
    
    if(plotType == "DotProduct"){
      
      compPlot <- ggplot(data = plotDat, aes(x = Name, y = metric)) +
        geom_bar(aes(color = color, fill = fill), stat = "identity", position = "dodge", width = 0.75) + 
        barplot_theme + scale_x_discrete(name = "Reactions") +
        scale_y_continuous(name = "Unit dot product", expand = c(0,0)) +
        coord_cartesian(ylim = c(round(min(plotDat$metric) - 0.1, 1),1)) +
        geom_text(data = nullRegulation, aes(x = Name, y = (metric + round(min(plotDat$metric) - 0.1, 1))/2, label = N)) +
        scale_color_identity() + scale_fill_identity() +
        geom_hline(yintercept = round(min(plotDat$metric) - 0.1, 1), size = 2) +
        geom_vline(xintercept = -Inf, size = 2) + scale_size_identity()
        
    }  
    
    if(length(unique(plotDat$ncond)) != 1){
     compPlot <- compPlot + facet_wrap(~ ncond, scale = "free_x")
    }
    
    pathwayPlot_list[[plotType]] <- compPlot
  }
  return(pathwayPlot_list)
}
  
  

pathwayPlots <- function(pathway){
  
  require(ggplot2)
  require(dplyr)
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = 20, color = "black", hjust = 1), axis.line = element_blank(),
                         strip.background = element_rect(fill = "cornflowerblue"), axis.ticks.y = element_line(size = 1, color = "black"))
  
  
  #### Generate plots which show how well the optimization is performing according to different metrics ####
  #### pathway gives a subset of reactions, of which a subset which are significant ########################
  #### after FDR correction (or are the default form) will be used #########################################
  
  reactions <- rxToPW$rID[rxToPW$pathway == pathwaySet$pathway[pathwaySet$display == pathway]]
  
  rxInfo_subset <- reactionInfo %>% filter(reaction %in% reactions) %>% filter(modelType %in% c("rMM", "irreversible", "regulator", "cooperativity", "2+ regulators"))
  
  noReg <- rxInfo_subset %>% filter(modelType == "rMM" | (modelType == "irreversible" & (!is.na(Qvalue) & Qvalue < 0.1)))
  oneReg <- rxInfo_subset %>% filter(modelType == "regulator") %>% group_by(reaction) %>% filter(ML == max(ML))
  
  multiReg <- rxInfo_subset %>% filter(modelType %in% c("cooperativity", "2+ regulators")) %>% left_join(oneReg %>% dplyr::select(reaction, ncond, one_reg_ML = ML), by = c("ncond", "reaction"))
  multiReg <- multiReg %>% filter(ML > one_reg_ML) %>% group_by(reaction) %>% filter(ML == max(ML)) %>% dplyr::select(-one_reg_ML)
  
  rxInfo_plotted <- rbind(noReg, oneReg, multiReg) 
  
  rxInfo_plotted <- rxInfo_plotted %>% left_join(rxn_fits %>% dplyr::select(rMech = rxn, parSpearman), by = "rMech") %>%
    left_join(fraction_flux_deviation %>% dplyr::select(rMech = rxn, dotProduct, angle, get("Interval Overlap")), by = "rMech") %>%
    mutate(modelType = ifelse(modelType %in% c("2+ regulators", "cooperativity"), "coop/2+reg", modelType))
  
  # Determine how many rxnForms that improve fit are masked (because there is a stronger fit)
  maskedReg <- rxInfo_subset %>% filter(!(rMech %in% rxInfo_plotted$rMech) & modelType != "irreversible") %>% mutate(modelType = ifelse(modelType %in% c("2+ regulators", "cooperativity"), "coop/2+reg", modelType)) %>%
    group_by(reaction, modelType, ncond) %>% dplyr::summarize(N = n())
  
  barplot_theme <- barplot_theme + theme(axis.text.x = element_text(size = min(8 + 80/(nrow(rxInfo_subset)+nrow(maskedReg)), 24), color = "black", vjust = 0.5, hjust = 1, angle = 90))
  
  pathwayPlot_list <- list()
  
  for(plotType in c("SpCorr", "DotProduct")){
    
    #### Define plotting variables ####
    if(plotType == "SpCorr"){
      plotDat <- rxInfo_plotted %>% dplyr::select(reaction, FullName, modelType, ncond, metric = parSpearman)
      good_dir = "+"
    }
    if(plotType == "DotProduct"){
      plotDat <- rxInfo_plotted %>% dplyr::select(reaction, FullName, modelType, ncond, metric = dotProduct)
      good_dir = "+"
    }
    
    plotDat <- plotDat %>% filter(!is.na(metric))
    #### Order reactions first by lowest/highest value for a rxn and then within a reaction ####
    
    nullRegulation <- maskedReg %>% ungroup() %>% mutate(FullName = paste(reaction, "also improved by", modelType)) %>% left_join(plotDat %>% filter(modelType == "rMM") %>% dplyr::select(reaction, ncond, metric), by = c("reaction", "ncond")) %>%
      ungroup() %>% dplyr::mutate(modelType = paste("null", modelType))
    
    plotDat <- rbind(plotDat, nullRegulation %>% dplyr::select(-c(N)))
    
    classAesthetics <- data.frame(modelType = c("irreversible", "rMM", "alternative_simple_kinetics", "null regulator", "regulator", "null coop/2+reg",
                                                "coop/2+reg", "hypo met regulator", "hypo met cooperativity"),
                                  fill = c("cornflowerblue", "dodgerblue3", "blueviolet", "white", "firebrick1", "white", "limegreen", "darkgoldenrod1", "darkgoldenrod3"),
                                  color = c("black", "black", "black", "firebrick1", "black", "limegreen", "black", "black", "black"))
    
    plotDat <- plotDat %>% left_join(classAesthetics, by = "modelType") %>%
      mutate(modelType = factor(modelType, levels = classAesthetics$modelType))
    
    if(good_dir == "+"){
      
      rMM_order <- plotDat %>% filter(modelType == "rMM") %>% group_by(reaction) %>% dplyr::summarize(max_metric = max(metric)) %>% arrange(max_metric)
      plotDat <- plotDat %>% mutate(reaction = factor(reaction, levels = rMM_order$reaction)) %>% ungroup() %>% arrange(reaction, ncond, modelType)
      
    }else{
      
      rMM_order <- plotDat %>% filter(modelType == "rMM") %>% group_by(reaction) %>% dplyr::summarize(min_metric = min(metric)) %>% arrange(desc(min_metric))
      
      plotDat <- plotDat %>% mutate(reaction = factor(reaction, levels = rMM_order$reaction)) %>% ungroup() %>% arrange(reaction, modelType)
      
    }
    
    plotDat <- plotDat %>% ungroup() %>% mutate(FullName = ifelse(modelType == "coop/2+reg", gsub('\\+', '+\n', FullName), FullName)) %>%
      mutate(FullName = factor(FullName, levels = unique(FullName)))
    
    nullRegulation <- nullRegulation %>% mutate(reaction = factor(reaction, levels = levels(plotDat$reaction)))
    
    if(plotType == "SpCorr"){
      
      compPlot <- ggplot(data = plotDat, aes(x = FullName, y = metric)) +
        geom_bar(aes(color = color, fill = fill), stat = "identity", position = "dodge", width = 0.75) + 
        barplot_theme + scale_x_discrete(name = "Reactions") +
        scale_y_continuous(expression("Spearman correlation between "~ V^FBA ~ "&" ~ V^PAR), expand = c(0,0)) + expand_limits(y = 1) +
        geom_text(data = nullRegulation, aes(x = FullName, y = metric/2, label = N)) +
        scale_color_identity() + scale_fill_identity() +
        geom_hline(yintercept = -Inf, size = 2) +
        geom_vline(xintercept = -Inf, size = 2) + scale_size_identity() +
        facet_grid(~ reaction + ncond, scale = "free_x", space = "free_x")
        
    }
    
    if(plotType == "DotProduct"){
      
      compPlot <- ggplot(data = plotDat, aes(x = FullName, y = metric,)) +
        geom_bar(aes(color = color, fill = fill), stat = "identity", position = "dodge", width = 0.75) + 
        barplot_theme + scale_x_discrete(name = "Reactions") +
        scale_y_continuous(name = "Unit dot product", expand = c(0,0)) +
        coord_cartesian(ylim = c(round(min(plotDat$metric) - 0.1, 1),1)) +
        scale_color_identity() + scale_fill_identity() +
        geom_hline(yintercept = round(min(plotDat$metric) - 0.1, 1), size = 2) +
        geom_vline(xintercept = -Inf, size = 2) + scale_size_identity() +
        facet_grid(~ reaction + ncond, scale = "free_x", space = "free_x") +
        geom_text(data = nullRegulation, aes(x = FullName, y = (metric + round(min(plotDat$metric) - 0.1, 1))/2, label = N))
        
    }  
    
    pathwayPlot_list[[plotType]] <- compPlot
  }
  return(pathwayPlot_list)
}


customPlots <- function(run_rxn, flux_fit, chemostatInfo){
  
  require(data.table)
  require(ggplot2)
  require(RColorBrewer)
  require(reshape2)
  require(grid)
  
  ### Generate plots that are relevent for a subset of reaction: ###
  #1: Flux predicted from FBA and parametric form by condition
  #2: Substrates, products, enzymes and regulators each get their own pane and are visualized together
  
  output_plots <- list()
  
  flux <- run_rxn$flux
  n_c <- nrow(flux)
  
  met_abund <- run_rxn$rxnSummary$rxnMet
  met_abund <- met_abund[,colSums(is.na(met_abund)) == 0, drop = F]
  
  enzyme_abund <- run_rxn$enzymes
  
  Chemoconds <- data.frame(name = factor(rownames(enzyme_abund), levels = rownames(enzyme_abund)),
                           Limitation = factor(substr(rownames(enzyme_abund),1,1), levels = unique(substr(rownames(enzyme_abund),1,1))),
                           DR = gsub('^[A-Z]', '', rownames(enzyme_abund)))
  Chemoconds$actualDR <- chemostatInfo$actualDR[chmatch(as.character(Chemoconds$name), chemostatInfo$ChemostatCond)]
    
  
  DRordering <- data.frame(DRs = as.character(c("0.05", "0.11", "0.16", "0.22", "0.30")), order = 1:5)
  Chemoconds$DRorder <- DRordering$order[chmatch(as.character(Chemoconds$DR), DRordering$DRs)]
  
  #### Comparision of parametric and FBA fit ####
  
  ci_theme <- theme(text = element_text(size = 30, face = "bold"), title = element_text(size = 30, face = "bold"), panel.background = element_rect(fill = "azure"), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
                    legend.text = element_text(size = 40, face = "bold"), axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90))
  
  flux_plot_FBA <- data.frame(METHOD = "FBA", FLUX = (flux$FVAmin + flux$FVAmax)/2, LB = flux$FVAmin, UB = flux$FVAmax, condName = Chemoconds$name, condition = Chemoconds$Limitation, DR = Chemoconds$DR)
  flux_plot_PAR <- data.frame(METHOD = "PAR", FLUX = flux_fit$fitted_flux$fitted, LB = flux_fit$fitted_flux$fitted - 2*flux_fit$fitted_flux$SD, UB = flux_fit$fitted_flux$fitted + 2*flux_fit$fitted_flux$SD, condName = Chemoconds$name, condition = Chemoconds$Limitation, DR = Chemoconds$DR)
  flux_plot_comp <- rbind(flux_plot_FBA, flux_plot_PAR)
  
  flux_range <- range(c(0, flux_plot_comp$FLUX))
  
  flux_plot_header <- data.frame(label = c("FBA-determined flux at optimum", "Parametric fit (95% CI)"), METHOD = c("FBA", "PAR"), x = 1, y = max(flux_range) - diff(flux_range)/15 * c(0,1))
  
  flux_plot_data <- list(flux_plot_comp = flux_plot_comp, flux_plot_header = flux_plot_header, flux_range = flux_range)
  
  #### Combining flux carried and species into a single faceted plot ####

  all_species <- data.frame(met_abund, enzyme_abund)
  all_species_SD <- run_rxn$specSD
  colnames(all_species) <- colnames(all_species_SD) <- run_rxn$all_species$commonName[chmatch(colnames(all_species), run_rxn$all_species$rel_spec)]
  
  species_df <- melt(cbind(all_species, data.frame(condition = Chemoconds$Limitation, DR = Chemoconds$actualDR, FLB = flux$FVAmin, FUB = flux$FVAmax)), id.vars = c("condition", "DR", "FLB", "FUB"), value.name = "RA")
  species_df$SD <- melt(all_species_SD, value.name = "SD", measure.vars = colnames(all_species_SD))$SD
  species_df$LB <- species_df$RA - 2*species_df$SD
  species_df$UB <- species_df$RA + 2*species_df$SD
  
  
  species_df <- melt(cbind(all_species, condition = Chemoconds$name), id.vars = c("condition"), value.name = "RA")
  species_df$SD <- melt(all_species_SD, value.name = "SD", measure.vars = colnames(all_species_SD))$SD
  species_df <- data.table(species_df)
  species_df$metOrigin <- sapply(species_df$variable, function(x){
    if(x %in% unname(run_rxn$rxnSummary$metNames)){
      run_rxn$rxnSummary$originMet[chmatch(names(run_rxn$rxnSummary$metNames)[run_rxn$rxnSummary$metNames == x], names(run_rxn$rxnSummary$originMet))]
    }else{
     "rel" 
    }
  })
  species_df[,LB :=  2^(RA - 2*SD),]
  species_df[,UB :=  2^(RA + 2*SD),]
  species_df[,RA :=  2^RA,]
  
  # first scale absolute units by powers of 10 such that mean is between 1&10
  species_df[,units := as.character(if(all(metOrigin == "rel")){NA}else{
    meanPower <- format(mean(abs(RA)), scientific = T)
    regmatches(meanPower, regexpr('[-0-9]+$', meanPower))}),by = "variable"]
  species_df$RA[!is.na(species_df$units)] <- species_df$RA[!is.na(species_df$units)]/(10^as.numeric(species_df$units[!is.na(species_df$units)]))
  species_df$UB[!is.na(species_df$units)] <- species_df$UB[!is.na(species_df$units)]/(10^as.numeric(species_df$units[!is.na(species_df$units)]))
  species_df$LB[!is.na(species_df$units)] <- species_df$LB[!is.na(species_df$units)]/(10^as.numeric(species_df$units[!is.na(species_df$units)]))
  
  # Assign species to substrates, products, regulators, enzymes
  species_df$Pane <- NA
  species_df$Pane[species_df$variable %in% run_rxn$kineticPars$commonName[chmatch(run_rxn$rxnSummary$rxnFormData$SubstrateID[run_rxn$rxnSummary$rxnFormData$Subtype == "substrate"], run_rxn$kineticPars$modelName)]] <- "Substrate"
  species_df$Pane[species_df$variable %in% run_rxn$kineticPars$commonName[chmatch(run_rxn$rxnSummary$rxnFormData$SubstrateID[run_rxn$rxnSummary$rxnFormData$Subtype == "product"], run_rxn$kineticPars$modelName)]] <- "Product"
  species_df$Pane[species_df$variable %in% run_rxn$kineticPars$commonName[chmatch(run_rxn$rxnSummary$rxnFormData$SubstrateID[!(run_rxn$rxnSummary$rxnFormData$Subtype %in% c("substrate", "product"))], run_rxn$kineticPars$modelName)]] <- "Regulator"
  species_df$Pane[species_df$variable %in% colnames(enzyme_abund)] <- "Enzyme"
  species_df <- data.table(species_df)
  
  # scale relative units by arbitrary multiplier to align them as closely as possible to absolute data
  for(a_Pane in unique(species_df$Pane)){
    species_subset <- species_df[Pane == a_Pane,,]
    if(any(species_subset$metOrigin == "abs")){
      abs_mean <- mean(species_subset$RA[species_subset$metOrigin == "abs"])
      for(relSpecies in unique(species_subset$variable[species_subset$metOrigin == "rel"])){
        species_df$UB[species_df$variable == relSpecies] <- species_df$UB[species_df$variable == relSpecies] * abs_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$LB[species_df$variable == relSpecies] <- species_df$LB[species_df$variable == relSpecies] * abs_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$RA[species_df$variable == relSpecies] <- species_df$RA[species_df$variable == relSpecies] * abs_mean/mean(species_df$RA[species_df$variable == relSpecies])
      }
    }else{
      rel_mean <- mean(species_subset$RA)
      for(relSpecies in unique(species_subset$variable)){
        species_df$UB[species_df$variable == relSpecies] <- species_df$UB[species_df$variable == relSpecies] * rel_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$LB[species_df$variable == relSpecies] <- species_df$LB[species_df$variable == relSpecies] * rel_mean/mean(species_df$RA[species_df$variable == relSpecies])
        species_df$RA[species_df$variable == relSpecies] <- species_df$RA[species_df$variable == relSpecies] * rel_mean/mean(species_df$RA[species_df$variable == relSpecies])
      }
    }
  }
  
  # truncate upper bounds at max(RA) * 1.2
  species_df$Truncated <- F
  for(a_Pane in unique(species_df$Pane)){
    species_subset <- species_df[Pane == a_Pane,,]
    maxRA <- max(species_subset$RA)
    if(any(species_subset$UB > maxRA*1.3)){
      species_df$Truncated[species_df$UB > maxRA*1.3 & species_df$Pane == a_Pane] <- T
      species_df$UB[species_df$UB > maxRA*1.3 & species_df$Pane == a_Pane] <-  maxRA*1.3
      }
    }
  
  flux_plot <- data.table(melt(flux_plot_comp, id.vars = c("condName", "FLUX", "LB", "UB"), measure.vars = "METHOD"))
  flux_plot <- flux_plot[,-which(colnames(flux_plot) == "variable"), with = F]
  flux_plot$variable <- as.character("Flux")
  flux_plot$Truncated <- F
  setnames(flux_plot, old = c("condName", "FLUX", "value", "variable"), new = c("condition", "RA", "variable", "Pane"))
  
  all_changes <- rbind(species_df[,chmatch(colnames(flux_plot), colnames(species_df)), with = F], flux_plot)
  all_changes$Pane <- factor(all_changes$Pane, levels = c("Substrate", "Product", "Enzyme", "Regulator", "Flux"))
  
  all_changes$condition <- factor(all_changes$condition, levels = rownames(met_abund))
  all_changes$Limitation <- factor(substr(all_changes$condition, 1, 1))
  all_changes$pathGroup <- factor(paste(all_changes$Limitation, all_changes$variable, sep = "_"))
  
  # Define species color
  
  variable_colors <- data.frame(unique(as.data.frame(all_changes)[,c("variable", "Pane")]), color = NA)
  
  if(sum(variable_colors$Pane %in% c("Substrate", "Product")) > 8){
    variable_colors$color[variable_colors$Pane %in% c("Substrate", "Product")] <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel2")[c(1:(sum(variable_colors$Pane %in% c("Substrate", "Product")) - 8))])
  }else{
    variable_colors$color[variable_colors$Pane %in% c("Substrate", "Product")] <- brewer.pal(8, "Dark2")[1:sum(variable_colors$Pane %in% c("Substrate", "Product"))]
    }
  variable_colors$color[variable_colors$Pane == "Flux"] <- brewer.pal(8, "Set1")[1:2]
  variable_colors$color[variable_colors$Pane == "Regulator"] <- rev(brewer.pal(9, "Reds"))[c(3,5,7,4,6)][1:sum(variable_colors$Pane == "Regulator")]
  variable_colors$color[variable_colors$Pane == "Enzyme"] <- rev(brewer.pal(9, "YlGnBu"))[c(1,3,5,7,2,4,6,8,9)][1:sum(variable_colors$Pane == "Enzyme")]
  
  variable_colors <- variable_colors[chmatch(levels(all_changes$variable), as.character(variable_colors$variable)),]
  
  # Define species labels
  
  variable_labels <- data.frame(unique(as.data.frame(all_changes)[,c("variable", "Pane")]), units = NA, label = NA)
  variable_labels$units <- species_df$units[chmatch(as.character(variable_labels$variable), as.character(species_df$variable))]
  
  
  if(!all(is.na(variable_labels$units))){
    abs_units <- data.frame(units = as.numeric(variable_labels$units[!is.na(variable_labels$units)]), referenceSymbol = NA, referencePower = NA, scaled =  NA)
    
    abs_units[abs_units$units >= 0, 2] <- "M"; abs_units[abs_units$units >= 0, 3] <- 0
    abs_units[data.table::between(abs_units$units, -3, -1, incbounds = T), 2] <- "mM"; abs_units[data.table::between(abs_units$units, -3, -1, incbounds = T), 3] <- -3
    abs_units[data.table::between(abs_units$units, -6, -4, incbounds = T), 2] <- "uM"; abs_units[data.table::between(abs_units$units, -6, -4, incbounds = T), 3] <- -6
    abs_units[abs_units$units < -7, 2] <- "nM"; abs_units[abs_units$units < -7, 3] <- -9
    
    abs_units$scaled <- 10^(abs_units$units - as.numeric(abs_units$referencePower))
    abs_units$scaled[abs_units$scaled == 1] <- ""
    
    variable_labels$units[!is.na(variable_labels$units)] <- paste0(abs_units$scaled, abs_units$referenceSymbol)
    }
  variable_labels$units[is.na(variable_labels$units)] <- "a.u."
  
  variable_labels$label <- paste0(variable_labels$variable, ' (', variable_labels$units, ')')
  variable_labels$label[variable_labels$variable == "FBA"] <- "Observed flux"
  variable_labels$label[variable_labels$variable == "PAR"] <- "Michaelis-Menten flux (95% CI)"
  tmp <- variable_labels$y[variable_labels$variable == "FBA"]
  variable_labels$y[variable_labels$variable == "FBA"] <- variable_labels$y[variable_labels$variable == "PAR"]
  variable_labels$y[variable_labels$variable == "PAR"] <- tmp
  
  variable_labels$x = "P0.05"
  variable_labels$y = NA
  
  for(a_Pane in unique(variable_labels$Pane)){
    PaneMaxVal <- max(all_changes$UB[all_changes$Pane == a_Pane])
    variable_labels$y[variable_labels$Pane == a_Pane] <- PaneMaxVal - c(1:sum(variable_labels$Pane == a_Pane))*PaneMaxVal*0.12
    }
  
  multiple_species_plot_data <- list(all_changes = all_changes, variable_labels = variable_labels)
  
  # Flux Summary
  
  ci_theme <- theme(text = element_text(size = 30, face = "bold"), title = element_text(size = 30, face = "bold"), panel.background = element_rect(fill = "azure"), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
                    legend.text = element_text(size = 40, face = "bold"), axis.text = element_text(color = "black", size = 30), axis.text.x = element_blank(),
                    axis.ticks.y = element_line(size = 3, color = "black"), axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, "lines"))
  
  ggplot() + geom_hline(y = 0, size = 3) + geom_vline(x = 0, size = 3) +
    geom_path(data = flux_plot_data$flux_plot_comp %>% filter(METHOD == "FBA"), aes(x = condName, y = FLUX, col = factor(METHOD), group = condition), size = 4, alpha = 1) +
    geom_point(data = flux_plot_data$flux_plot_comp %>% filter(METHOD == "FBA"), aes(x = condName, y = FLUX, col = factor(METHOD)), size = 10, alpha = 1) +
    ci_theme + scale_color_manual("Method", guide = "none", values= c("FBA" = "firebrick1", "PAR" = "darkblue")) +
    scale_y_continuous("Relative Flux", limits = c(flux_range[1],flux_range[2]*1.05), expand = c(0,0)) + scale_x_discrete("Conditions") +
    scale_size_identity()
  
  ggsave(paste0("Figures/rxnPlots/", arxn, "-flux", ".pdf"), height = 6, width = 10)
  
  corr_label <- data.frame(x = "U0.11", y = flux_range[which.max(abs(flux_range))]*0.95, label = paste("r^2", "~'=", round(flux_fit$fit_summary$parPearson^2, 2), "'"))
  
  flux_plot_data$flux_plot_comp$UB[flux_plot_data$flux_plot_comp$UB >= flux_range[2]*1.05] <- flux_range[2]*1.05
  flux_plot_data$flux_plot_comp$LB[flux_plot_data$flux_plot_comp$LB <= flux_range[1]] <- flux_range[1]
  
  ggplot() + geom_hline(y = 0, size = 3) + geom_vline(x = 0, size = 3) +
    geom_path(data = flux_plot_data$flux_plot_comp %>% filter(METHOD == "FBA"), aes(x = condName, y = FLUX, col = factor(METHOD), group = condition), size = 4, alpha = 1) +
    geom_point(data = flux_plot_data$flux_plot_comp %>% filter(METHOD == "FBA"), aes(x = condName, y = FLUX, col = factor(METHOD)), size = 10, alpha = 1) +
    geom_path(data = flux_plot_data$flux_plot_comp %>% filter(METHOD == "PAR"), aes(x = condName, y = FLUX, col = factor(METHOD), group = condition), size = 4, alpha = 1) +
    geom_text(data = corr_label, aes(x = x, y = y, label = label), size = 8, parse = T) +
    geom_pointrange(data = flux_plot_data$flux_plot_comp %>% filter(METHOD == "PAR"), aes(x = condName, y = FLUX, ymin = LB, ymax = UB, col = factor(METHOD)), size = 2.6, alpha = 1, width = 10) +
    ci_theme + scale_color_manual("Method", guide = "none", values= c("FBA" = "firebrick1", "PAR" = "darkblue")) +
    scale_y_continuous("Relative Flux", limits = c(flux_range[1],flux_range[2]*1.05), expand = c(0,0)) + scale_x_discrete("Conditions") +
    scale_size_identity()
  
  ggsave(paste0("Figures/rxnPlots/", arxn, "-fluxMatch", ".pdf"), height = 6, width = 10)

  # Species Summary
  
  ci_theme <- theme(text = element_text(size = 30, face = "bold"), title = element_text(size = 30, face = "bold"), panel.background = element_rect(fill = "azure"), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
                    legend.text = element_text(size = 40, face = "bold"), axis.text = element_text(color = "black", size = 30), axis.text.x = element_blank(),
                    axis.ticks.y = element_line(size = 3, color = "black"), axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, "lines"),
                     panel.margin = unit(1, "lines"))
  
  
  multiple_species_plot_data$all_changes <- multiple_species_plot_data$all_changes %>% filter(Pane != "Flux")
  multiple_species_plot_data$variable_labels <- multiple_species_plot_data$variable_labels %>% filter(Pane != "Flux")
  
  
  ggplot() + geom_path(data = multiple_species_plot_data$all_changes, aes(x = condition, y = RA, col = variable, group = pathGroup), size = 4, alpha = 1) +
    geom_point(data = multiple_species_plot_data$all_changes, aes(x = condition, y = RA, col = factor(variable)), size = 10) +
    geom_text(data = multiple_species_plot_data$variable_labels, aes(x = x, y = y, col = variable, label = label), hjust = 0) +
    geom_blank(data = multiple_species_plot_data$all_changes, aes(y = RA*1.05)) + 
    geom_vline(x = 0, size = 3) + geom_hline(y = 0, size = 3) +
    facet_grid(Pane ~ ., scale = "free_y") +
    scale_color_manual(guide = "none", labels = variable_colors$variable, values = variable_colors$color) +
    ci_theme + scale_size_identity() +
    scale_y_continuous("Concentration and Flux", expand = c(0,0)) + scale_x_discrete("Conditions") + expand_limits(y = 0)
 
  ggsave(paste0("Figures/rxnPlots/", arxn, "-metabolites", ".pdf"), height = 2 + 3*length(unique(multiple_species_plot_data$all_changes$Pane)), width = 10)
  
}





reactionProperties <-  function(){
  
  ### Generate plots which combine both MCMC parameter estimates and experimental data ###
  
  output <- list()
  output_plots <- list() # the output
  
  require(data.table)
  require(ggplot2)
  require(scales) # for converting fractions to percent
  require(reshape2)
  require(dplyr)
  
  ### Define and align experimental data ###
  
  n_c <- nrow(run_rxn$metabolites)
  
  flux <- run_rxn$flux
  met_abund <- run_rxn$rxnSummary$rxnMet
  met_abund[is.na(met_abund)] <- 0
  enzyme_abund <- run_rxn$enzymes
  
  mle_pars <- par_markov_chain[which.max(par_likelihood$likelihood),]
  
  if('t_metX' %in% run_rxn$kineticPars$modelName){
  
    ### hypothetical metabolite overwriten with value from PCs ###
    releventPCs <- metSVD$v[rownames(metSVD$v) %in% rownames(met_abund),]
    npc <- ncol(releventPCs)
  
    if(npc != sum(run_rxn$kineticParPrior$SpeciesType == "PCL")){
      print(paste("Inconsistent number of principal components for rx", arxn, "likely version problem -> rerun parameter sets"))
      return(NA)
      }
    
    met_abund$t_metX <- t(mle_pars[run_rxn$kineticParPrior$SpeciesType == "PCL"] %*% diag(metSVD$d[1:npc]) %*% t(releventPCs))
    
  }
  
  dist_pars <- par_markov_chain[,run_rxn$kineticPars$SpeciesType != "PCL"]
  colnames(dist_pars) <- run_rxn$kineticPars$formulaName[run_rxn$kineticPars$SpeciesType != "PCL"]
  
  ### Calculate occupancy ###
  
  measured_mets <- met_abund[,colSums(met_abund == 0) != n_c, drop = F]
  
  if("hillCoefficient" %in% run_rxn$kineticPars$SpeciesType){
    hill_species = matrix(1, ncol = ncol(measured_mets), nrow = nrow(dist_pars))
    allosteric_spec <- run_rxn$rxnSummary$rxnFormData$SubstrateID[run_rxn$rxnSummary$rxnFormData$Subtype == "allosteric"]
    hill_species[,colnames(measured_mets) == allosteric_spec] <- 2^dist_pars[,colnames(dist_pars) == "h_allo"]
    }else{hill_species = 1}
  
  measured_met_affinity <- dist_pars[,chmatch(colnames(measured_mets), run_rxn$kineticPars$rel_spec), drop = F]
  
  affinity_array <- array(data = NA, dim = c(nrow(dist_pars), n_c, ncol(measured_met_affinity)))
  dimnames(affinity_array) <- list(markovSample = c(1:nrow(dist_pars)), condition = rownames(measured_mets), metabolite = run_rxn$kineticPars$formatted[chmatch(colnames(measured_met_affinity), run_rxn$kineticPars$formulaName)])
  
  for(cond in 1:n_c){
    
    # 1 / ( (Km/[L])^h  + 1)
    
    affinity_array[,cond,] <- 1/((2^measured_met_affinity / (rep(1, nrow(measured_met_affinity)) %*% t(t(2^measured_mets[cond,]))))^hill_species + 1)
    
    }
  
  affinity_array_melt <- data.table(melt(affinity_array))
  
  affinity_summary <- affinity_array_melt[,list(LB = boxplot.stats(value)$stats[1], LH = boxplot.stats(value)$stats[2], median = boxplot.stats(value)$stats[3],
                                                UH = boxplot.stats(value)$stats[4], UB = boxplot.stats(value)$stats[5]), by = c("metabolite", "condition")]
  affinity_summary$condition <- factor(affinity_summary$condition, levels = chemostatInfo$ChemostatCond[chemostatInfo$ChemostatCond %in% affinity_summary$condition])
  
  output_plots$affinity_summary <- affinity_summary
  
  ### Calculate flux elasticity over all parameter sets, conditions and non-constant species ###
  
  # take partial derivatives of elasticity_calc
  elast_equation <- eval(parse(text = paste('expression(', run_rxn$occupancyEq$elasticity_calc, ')')))
  
  elast_partials <- list()
  for(spec in names(run_rxn$occupancyEq$kinetic_form_partials)){
    elast_partials[[spec]] <- D(elast_equation, spec)
  }
  
  all_species <- data.frame(2^met_abund, 2^enzyme_abund)
  all_species <- all_species[,chmatch(names(elast_partials), colnames(all_species))]
  
  # calculate sensitivities #  partial F / partial S
  kinetically_differing_isoenzymes <- any(names(run_rxn$occupancyEq[["l_occupancyExpression"]]) %in% rownames(run_rxn$rxnSummary$enzymeComplexes))
  mcmc_elasticity <- sapply(1:nrow(dist_pars), function(x){
    
    pars <- dist_pars[x,]
    dist_par_stack <- rep(1, n_c) %*% t(2^(pars))
    
    occupancy_vals <- data.frame(met_abund, dist_par_stack)
    
    #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
    #occupany of enzymes * relative abundance of enzymes
    
    if(!(kinetically_differing_isoenzymes)){
      predOcc <- eval(run_rxn$occupancyEq[["l_occupancyExpression"]], occupancy_vals) #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
      enzyme_activity <- (predOcc %*% t(rep(1, ncol(enzyme_abund))))*2^enzyme_abund #occupany of enzymes * relative abundance of enzymes
    }else{
      occEqtn_complex_match <- data.frame(complex = colnames(enzyme_abund), occEqtn = NA)
      occEqtn_complex_match$occEqtn[occEqtn_complex_match$complex %in% names(run_rxn$occupancyEq[["l_occupancyExpression"]])] <- occEqtn_complex_match$complex[occEqtn_complex_match$complex %in% names(run_rxn$occupancyEq[["l_occupancyExpression"]])]
      occEqtn_complex_match$occEqtn[is.na(occEqtn_complex_match$occEqtn)] <- "other"
      
      enzyme_activity <- NULL
      for(isoenzyme in names(run_rxn$occupancyEq[["l_occupancyExpression"]])){
        predOcc <- eval(run_rxn$occupancyEq[["l_occupancyExpression"]][[isoenzyme]], occupancy_vals)
        enzyme_activity <- cbind(enzyme_activity, predOcc %*% t(rep(1, sum(occEqtn_complex_match$occEqtn == isoenzyme)))*2^enzyme_abund[,colnames(enzyme_abund) %in% occEqtn_complex_match$complex[occEqtn_complex_match$occEqtn == isoenzyme]])
      }
    }
    
    # Calculate Vmax as during standard fitting
    
    flux_fit <- nnls(enzyme_activity, (flux$FVAmax + flux$FVAmin)/2)
    nnlsCoef <- t(t(rep(1, n_c)))  %*% flux_fit$x; colnames(nnlsCoef) <- run_rxn$all_species$formulaName[run_rxn$all_species$SpeciesType == "Enzyme"]
    
    # Setup all species in linear-space
    
    occupancy_vals <- data.frame(2^met_abund, dist_par_stack)
    all_components <- data.frame(occupancy_vals, 2^enzyme_abund, nnlsCoef)
    
    # partial derivatives of each measured specie in a condition
    
    comp_partials <- matrix(NA, nrow = n_c, ncol = length(elast_partials))
    colnames(comp_partials) <- names(elast_partials)
    
    for(j in 1:ncol(comp_partials)){
      comp_partials[,j] <- with(all_components ,eval(elast_partials[[j]]))
    }
    
    # partial F / partial S * [S]/F
    as.matrix(comp_partials * all_species/flux_fit$fitted)
      
  }, simplify = "array")
  
  flux_names <- names(elast_partials)
  flux_names[flux_names %in% run_rxn$kineticPars$rel_spec] <- run_rxn$kineticPars$formatted[chmatch(flux_names[flux_names %in% run_rxn$kineticPars$rel_spec], run_rxn$kineticPars$rel_spec)]
  
  dimnames(mcmc_elasticity) <- list(condition = rownames(measured_mets), specie = flux_names, markovSample = c(1:nrow(dist_pars)))
  
  output$EL_summary <- melt(mcmc_elasticity, value.name = "Elasticity") # Save full distributions of elasticitiies for MCA
  
  ### Plot elasticity of log-abundances ###
  
  flux_elasticity_melt <- data.table(melt(mcmc_elasticity))
  
  elasticity_summary <- flux_elasticity_melt[,list(LB = boxplot.stats(value)$stats[1], LH = boxplot.stats(value)$stats[2], median = boxplot.stats(value)$stats[3],
                                                   UH = boxplot.stats(value)$stats[4], UB = boxplot.stats(value)$stats[5]), by = c("specie", "condition")]
  elasticity_summary$condition <- factor(elasticity_summary$condition, levels = chemostatInfo$ChemostatCond[chemostatInfo$ChemostatCond %in% elasticity_summary$condition])
  
  output_plots$elasticity_summary <- elasticity_summary
  
  
  ### Physiological leverage: absolute partial correlation weighted by across-condition SD
  
  all_exp_species <- data.frame(enzyme_abund, measured_mets[, colnames(measured_mets) %in% run_rxn$kineticPars$rel_spec[run_rxn$kineticPars$formatted %in% flux_names], drop = F])
  colnames(all_exp_species)[colnames(all_exp_species) %in% run_rxn$kineticPars$rel_spec] <- run_rxn$kineticPars$formatted[chmatch(colnames(all_exp_species)[colnames(all_exp_species) %in% run_rxn$kineticPars$rel_spec], run_rxn$kineticPars$rel_spec)]
  
  all_exp_species <- all_exp_species[,chmatch(colnames(mcmc_elasticity), colnames(all_exp_species))]
  
  # Calculate the physiological SD from total Variance - within Condition Variance
  
  #condSD <- run_rxn$specSD
  #colnames(condSD)[colnames(condSD) %in% run_rxn$kineticPars$rel_spec] <- run_rxn$kineticPars$formatted[chmatch(colnames(condSD)[colnames(condSD) %in% run_rxn$kineticPars$rel_spec], run_rxn$kineticPars$rel_spec)]
  #WICSS <- apply(condSD^2, 2, sum)
  #OCSS <- apply((all_exp_species - rep(1,nrow(all_exp_species)) %*% t(apply(all_exp_species, 2, mean)))^2, 2, sum)
  
  # weighting by sd of log abundances is similar to weighting by IQR (with better stability)
  # Across all conditions
  aSDweightedElasticity <- abs(mcmc_elasticity) * array(rep(1, n_c) %*% t(apply(all_exp_species, 2, sd)), dim = dim(mcmc_elasticity)) 
  # Across natural conditions
  nSDweightedElasticity <- abs(mcmc_elasticity[grep('^[PCN]', rownames(mcmc_elasticity)),,]) *
    array(rep(1, length(grep('^[PCN]', rownames(mcmc_elasticity)))) %*% t(apply(all_exp_species[grep('^[PCN]', rownames(mcmc_elasticity)),], 2, sd)), dim = c(length(grep('^[PCN]', rownames(mcmc_elasticity))), dim(mcmc_elasticity)[2:3])) 
  
  we_melt <- rbind(tbl_df(melt(aSDweightedElasticity)) %>% dplyr::mutate(conditions = "ALL"),
                   tbl_df(melt(nSDweightedElasticity)) %>% dplyr::mutate(conditions = "NATURAL"))
  
  we_melt <- we_melt %>% group_by(condition, markovSample, conditions) %>% dplyr::mutate(total_we = sum(value)) %>% filter(total_we != 0) %>%
    dplyr::mutate(physiological_leverage = value/total_we)
  
  we_summary <- we_melt %>% group_by(specie, condition, conditions) %>% dplyr::summarize(LB = boxplot.stats(physiological_leverage)[['stats']][1], LH = boxplot.stats(physiological_leverage)[['stats']][2], median = boxplot.stats(physiological_leverage)[['stats']][3],
                                                                                  UH = boxplot.stats(physiological_leverage)[['stats']][4], UB = boxplot.stats(physiological_leverage)[['stats']][5])
  we_summary <- we_summary %>% ungroup() %>% mutate(Limitation = factor(substr(condition, 1, 1), levels = c("P", "C", "N", "L", "U")),
                                                    DR = factor(sub('[A-Z]', '', condition)),
                                                    condition = factor(condition, levels = chemostatInfo$ChemostatCond[chemostatInfo$ChemostatCond %in% condition]))
  
  output_plots$we_summary <- we_summary
  
  ### Summarize physiological leverage to look for general trends over reactions ###
  
  # leverage implied by the most likely parameter set
  ML_MLE <- we_melt %>% filter(markovSample == which.max(par_likelihood$likelihood)) %>% ungroup() %>% dplyr::select(condition, specie, conditions, MLE = physiological_leverage)
  # distribution of leverage over posterior credibility interval
  ML_summary <- we_melt %>% group_by(specie, condition, conditions) %>% dplyr::summarize("0.025" = quantile(physiological_leverage, probs = 0.025), "0.5" = quantile(physiological_leverage, probs = 0.5), 
                              "0.975" = quantile(physiological_leverage, probs = 0.975)) %>% ungroup()
  
  ML_summary <- ML_summary %>% left_join(ML_MLE, by = c("specie", "condition", "conditions"))
  
  ML_summary$Type <- NA
  ML_summary$Type[ML_summary$specie %in% colnames(enzyme_abund)] <- "enzyme"
  
  ML_summary$Type[is.na(ML_summary$Type)] <- run_rxn$rxnSummary$rxnFormData$Subtype[chmatch(run_rxn$kineticPars$rel_spec[chmatch(as.character(ML_summary$specie[is.na(ML_summary$Type)]), run_rxn$kineticPars$formatted)], run_rxn$rxnSummary$rxnFormData$SubstrateID)]
  ML_summary$reaction = arxn
  
  ML_summary$condition <- factor(ML_summary$condition, levels = chemostatInfo$ChemostatCond[chemostatInfo$ChemostatCond %in% ML_summary$condition])
  
  output$ML_summary = ML_summary
  
  output$plots <- output_plots
  return(output)
  
}



reactionPropertiesPlots <- function(reaction_plots){
  
  # Plot output from reactionProperties.  This was split into two functions, because when trying to return plots from reactionProperties, MASSIVE pdfs were generated without a clear cause #
  
  boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 90, hjust = 1), axis.line = element_blank(),
                         axis.text = element_text(color = "black"), strip.background = element_rect(fill = "cornsilk1"), strip.text = element_text(color = "darkblue"))
  
  
  output_plots <- NULL
  
  output_plots$"Metabolite Occupancy" <- ggplot(reaction_plots$affinity_summary, aes(x = condition, ymin = LB, lower = LH, middle = median, upper = UH, ymax = UB, fill = "cornsilk1", color = "brown1")) + facet_wrap(~ metabolite, scales = "free_y") + 
    geom_boxplot(stat = "identity") + boxplot_theme + scale_fill_identity() + scale_color_identity() + expand_limits(y=c(0,1)) +
    scale_x_discrete("Experimental Condition") + scale_y_continuous(name = expression("Metabolite Occupancy: "~ frac("[M]"^"h", "[M]"^"h" + K[M]^"h")), expand = c(0,0))
  
  output_plots$"Flux Elasticity" <- ggplot(reaction_plots$elasticity_summary, aes(x = condition, ymin = LB, lower = LH, middle = median, upper = UH, ymax = UB, fill = "cornsilk1", color = "brown1")) + facet_wrap(~ specie, scales = "free_y", ncol = 2) + 
    geom_boxplot(stat = "identity") + boxplot_theme + scale_fill_identity() + scale_color_identity() + expand_limits(y=0) +
    scale_x_discrete("Experimental Condition") + scale_y_continuous(name = expression("Elasticity: "~ frac(rho~"F", rho~S)~frac("[S]","F")))
  
  output_plots$"Metabolic Leverage (all)" <- ggplot(reaction_plots$we_summary %>% filter(conditions == "ALL"), aes(x = condition, ymin = LB, lower = LH, middle = median, upper = UH, ymax = UB, fill = factor(specie))) + 
    geom_boxplot(stat = "identity") + expand_limits(y=c(0,1)) + boxplot_theme +
    scale_y_continuous("Metabolic Leverage", labels = percent_format(), expand = c(0,0)) + scale_x_discrete("Experimental Condition") +
    scale_fill_discrete() + facet_wrap(~ specie, ncol = 2)
  
  output_plots$"Metabolic Leverage (natural)" <- ggplot(reaction_plots$we_summary  %>% filter(conditions == "NATURAL"), aes(x = condition, ymin = LB, lower = LH, middle = median, upper = UH, ymax = UB, fill = factor(specie))) + 
    geom_boxplot(stat = "identity") + expand_limits(y=c(0,1)) + boxplot_theme +
    scale_y_continuous("Metabolic Leverage", labels = percent_format(), expand = c(0,0)) + scale_x_discrete("Experimental Condition") +
    scale_fill_discrete() + facet_wrap(~ specie, ncol = 2)
  
  
  return(output_plots)
  
}


transcriptional_responsiveness <- function(){
  

  ### Generate a plot which compare enzyme levels to transcript levels across these conditons ###
  ### Also, relate expected change in flux to transcriptional control ###
  ### TR = R2fm * R2pt * ML
  ### R2fm - Fraction of variance in flux accounted for by reaction mechanism
  ### R2pt - Fraction of variance in protein abundance due to linear changes in transcripts
  ### ML - Metabolic leverage - a hypothesis regarding the anticipated instantaneous change in flux when a protein is altered
  
  require(ggplot2)
  require(data.table)
  
  ci_theme <- theme(text = element_text(size = 30, face = "bold"), title = element_text(size = 30, face = "bold"), panel.background = element_rect(fill = "azure"), 
                    panel.grid.minor = element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "cyan"),
                    legend.text = element_text(size = 40, face = "bold"), axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90))
  
  TRplots <- list()
  
  #### Determine TR ####
  
  # Determine whether the enzyme can be catalyzed by a single enzyme (as opposed to a complex)
  enzyme_groups <- run_rxn$rxnSummary$enzGroup
  enzyme_groups <- enzyme_groups[enzyme_groups %in% names(table(enzyme_groups))[table(enzyme_groups) == 1]]
  
  R2fm <- c(flux_fit$fit_summary$parPearson_c^2)
  
  PTsubset <- PTcomparison[names(PTcomparison) %in% unique(names(enzyme_groups))]
  if(length(PTsubset) != 0){
    
    PTcorr <- unlist(lapply(PTsubset, function(x){x$PTcorr}))
    
    ### Determine transcriptional effectiveness for reactions where the reaction mechanism is predictive and transcriptional regulation of proteins has been demonstrated ###
    
    if(flux_fit$fit_summary$parPearson_c > 0 & any(PTcorr > 0)){
      
      positive_enzymes <- PTcorr[PTcorr > 0]
      
      positive_enzyme_groups <- enzyme_groups[names(enzyme_groups) %in% names(positive_enzymes)]
      
      positive_summary <- data.table(enzyme = names(positive_enzyme_groups), group = positive_enzyme_groups)
      positive_summary[,PTcorr := positive_enzymes[names(positive_enzymes) == enzyme], by = "enzyme"]
      positive_summary[,specie := paste(group, enzyme, sep = "_")]
      
      MLdata_subset <- reaction_properties$ML_summary %>% filter(conditions == "ALL") %>% as.data.frame() %>% data.table()
      MLdata_subset$specie <- as.character(MLdata_subset$specie)
      
      TR <- MLdata_subset[,list(q0.025 = get("0.025")*positive_summary$PTcorr[positive_summary$specie == specie]^2 * R2fm,
                                q0.5 = get("0.5")*positive_summary$PTcorr[positive_summary$specie == specie]^2 * R2fm,
                                q0.975 = get("0.975")*positive_summary$PTcorr[positive_summary$specie == specie]^2 * R2fm), 
                          by = c("specie", "condition", "reaction")]
      TR[,TU := 1 - positive_summary$PTcorr[positive_summary$specie == specie]^2, by = "specie"] #transcriptional uncertainty
      TR[,MU := 1 - R2fm,] # mechanistic uncertainty
      TR$limitation <- substr(levels(TR$condition)[TR$condition], 1, 1)
      TR$limitation <- factor(TR$limitation, levels = c("P", "C", "N", "L", "U"))
      
      TRcompliment <- TR[,list(source = "transcriptional uncertainty", fraction = (1 - (1-TU)*(1-MU))*TU/(TU + MU), alpha = 0.6), by = c("specie", "condition")]
      TRcompliment <- rbind(TRcompliment, TR[,list(source = "mechanistic uncertainty", fraction = (1 - (1-TU)*(1-MU))*MU/(TU + MU), alpha = 0.6), by = c("specie", "condition")])
      TRcompliment <- rbind(TRcompliment, TR[,list(source = "total controllable", fraction = (1-TU)*(1-MU), alpha = 0), by = c("specie", "condition")])
      TRcompliment$source <- factor(TRcompliment$source, levels = c("total controllable", "mechanistic uncertainty", "transcriptional uncertainty"))
      TRcompliment$condition <- factor(TRcompliment$condition, levels = levels(TR$condition))
      
      TRcomp_label <- TRcompliment[condition == "N0.16", list(condition = "N0.16", source = "mechanistic uncertainty", y = fraction[source == "total controllable"]+0.04), by = "specie"]
      TRcomp_label <- rbind(TRcomp_label, TRcompliment[condition == "N0.16", list(condition = "N0.16", source = "transcriptional uncertainty", y = sum(fraction[source %in% c("mechanistic uncertainty", "total controllable")])+0.04), by = "specie"])
      
      TRcompliment <- TRcompliment %>% tbl_df() %>% arrange(source)
      
      TRplots$Plots$"Transcriptional Responsiveness" <- ggplot(data = TR, aes(x = condition)) + geom_bar(data = TRcompliment, aes(y = fraction, fill = source, alpha = alpha), stat = "identity") + facet_wrap(~specie) +
        scale_alpha_identity("") + geom_pointrange(data = TR, aes(ymin = q0.025, y = q0.5, ymax = q0.975, color = limitation), size = 1) +
        geom_text(data = TRcomp_label, aes(label = source, y = y), size = 8) + ci_theme + scale_color_brewer("", guide = "none", palette = "Set1") +
        scale_y_continuous("Predicted Transcriptional Responsiveness", expand = c(0,0), labels = percent_format()) + scale_fill_brewer("Source of Departure", guide = "none", palette = "Set2") +
        ggtitle("Flux responsiveness to transcriptional changes") + scale_x_discrete("Nutrient Condition")
      
      TRplots$TR <- TR
      
    }
  }

  #### Generate Plots comparing proteins and transcripts ####
  
  scatter_theme <- theme(text = element_text(size = 23, face = "bold", color = "black"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), 
                         panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan"),
                         legend.text = element_text(size = 20, face = "bold"), axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90))
  
  
  enzyme_groups <- run_rxn$rxnSummary$enzGroup
  PTsubset <- PTcomparison[names(PTcomparison) %in% unique(names(enzyme_groups))]
  
  if(length(PTsubset) != 0){
    
    PTdt <- data.table(melt(lapply(PTsubset, function(x){
      x$PTabund
    }), id.vars = "Condition"))
    setnames(PTdt, c("L1", "variable", "value"), c("Gene", "Specie", "Relative Abundance"))
    PTdt$Limitation <- factor(substr(PTdt$Condition, 1, 1), levels = c("P", "C", "N", "L", "U"))
    PTdt$Condition <- factor(PTdt$Condition, levels = levels(MLdata$condition))
    
    TRplots$Plots$"Enzyme and Transcript Levels" <- ggplot(PTdt, aes(x = Condition, y = get("Relative Abundance"), color = Limitation, group = Limitation)) + geom_path(size = 2) + geom_point(size = 5) +
      facet_grid(Specie ~ Gene) + scatter_theme +
      scale_color_brewer("", guide = "none", palette = "Set1") + scale_y_continuous(expression(log[2]~" Relative Abundance")) +
      scale_x_discrete("Nutrient Condition") + ggtitle('Comparison of transcripts, \n enzymes and their ratio')
    
  }
  
  return(TRplots)
  
}



param_compare <- function(){
  
  ### Generate bivariate histograms of mcmc parameter estimates and marginal histograms ###
  
  require(dplyr)
  require(tidyr)
  
  ### write common names ###
  
  rename_table <- data.frame(tID = colnames(par_markov_chain), commonName = NA, commonPrint = NA)
  rename_table$commonName[rename_table$tID %in% names(run_rxn$rxnSummary$metNames)] <- unname(run_rxn$rxnSummary$metNames)[chmatch(rename_table$tID[rename_table$tID %in% names(run_rxn$rxnSummary$metNames)], names(run_rxn$rxnSummary$metNames))]
  rename_table$commonName[rename_table$tID == "keq"] <- "Keq"
  rename_table$commonPrint[!is.na(rename_table$commonName)] <- nameReformat(names = rename_table$commonName[!is.na(rename_table$commonName)], totalChar = 120)
  #print(rename_table)
  #break
  named_par_markov_chain <- par_markov_chain
  colnames(named_par_markov_chain) <- rename_table$commonPrint
  
  rename_table <- rename_table %>% left_join(run_rxn$rxnSummary$rxnFormData %>% dplyr::select(tID = SubstrateID, Subtype), by = "tID") %>%
    mutate(Subtype = ifelse(tID == "keq", "keq", Subtype),
           Subtype = ifelse(Subtype %in% c("substrate", "product", "keq"), Subtype, "regulator"),
           Subtype = factor(Subtype, levels = c("substrate", "product", "regulator", "keq"))) %>%
    arrange(Subtype)
  
  facet_ordering_key <- rename_table %>% dplyr::select(commonPrint, Subtype) %>% arrange(Subtype) %>%
    mutate(commonPrint = factor(commonPrint, levels = commonPrint))
  
  # adjust parameters which are bounded by log-uniform distributions into equal bounds
  
  named_par_markov_chain[,run_rxn$kineticParPrior$distribution == "unif"] <- named_par_markov_chain[,run_rxn$kineticParPrior$distribution == "unif"] -
    t(t(rep(1, nrow(named_par_markov_chain)))) %*% (run_rxn$kineticParPrior[run_rxn$kineticParPrior$distribution == "unif",'par_1'] + run_rxn$kineticParPrior[run_rxn$kineticParPrior$distribution == "unif",'par_2'])/2
  
  named_par_markov_chain <- named_par_markov_chain[,!is.na(rename_table$commonName)] # remove parameters which aren't affinities/keq
  
  # visualize the joint and marginal distribution of parameter values from the markov chain
  
  par_combinations <- expand.grid(1:ncol(named_par_markov_chain), 1:ncol(named_par_markov_chain))
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
  
  par_hist_binwidth = 0.4
  
  max_density <- max(apply(named_par_markov_chain, 2, function(x){max(table(round(x/par_hist_binwidth)))}))
  
  density_trans_inv <- function(x){x*(max_density/30) + max_density/2}
  density_trans <- function(x){(x - max_density/2)/(max_density/30)}
  
  par_comp_dissimilar$yval <- density_trans_inv(par_comp_dissimilar$yval)
  MLEpoints$yval <- density_trans_inv(MLEpoints$yval)
  
  # extend axes to bounds of prior
  
  bound_expand <- rename_table %>% left_join(run_rxn$kineticParPrior %>% dplyr::select(tID = rel_spec, lb = par_1, ub = par_2), by = "tID") %>% mutate(center = (lb + ub)/2, lb = lb - center, ub = ub - center) %>%
    gather(bound, value, lb:ub) %>% mutate(parameter_1 = commonPrint, parameter_2 = commonPrint) %>%
    mutate(value_trans = density_trans_inv(value))
  
  # add factors to paramers so that they are displayed in a logical order
  par_comp_like <- par_comp_like %>% mutate(parameter_1 = factor(parameter_1, levels = facet_ordering_key$commonPrint), parameter_2 = factor(parameter_2, levels = facet_ordering_key$commonPrint))
  par_comp_dissimilar <- par_comp_dissimilar%>% mutate(parameter_1 = factor(parameter_1, levels = facet_ordering_key$commonPrint), parameter_2 = factor(parameter_2, levels = facet_ordering_key$commonPrint))
  MLEbarplot <- MLEbarplot %>% mutate(parameter_1 = factor(parameter_1, levels = facet_ordering_key$commonPrint), parameter_2 = factor(parameter_2, levels = facet_ordering_key$commonPrint))
  MLEpoints <- MLEpoints %>% mutate(parameter_1 = factor(parameter_1, levels = facet_ordering_key$commonPrint), parameter_2 = factor(parameter_2, levels = facet_ordering_key$commonPrint))
  bound_expand <- bound_expand %>% mutate(parameter_1 = factor(parameter_1, levels = facet_ordering_key$commonPrint), parameter_2 = factor(parameter_2, levels = facet_ordering_key$commonPrint))
  
  hex_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                     panel.background = element_rect(fill = "gray92"), legend.position = "left", 
                     axis.ticks = element_line(color = "black", size = 1),
                     axis.text = element_text(color = "black", size = 20),
                     panel.grid.minor = element_blank(), panel.grid.major = element_line(size = 1),
                     axis.line = element_line(color = "black", size = 1), legend.key.height = unit(4, "line"),
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = element_rect(fill = "coral"),
                     panel.margin = unit(1.5, "lines"), axis.title = element_blank()
  )
  
  output_plots <- list()
  
  output_plots[["bivariateHist"]] <- ggplot() + geom_hex(data = par_comp_dissimilar, aes(x = xval, y = yval)) + facet_grid(parameter_2 ~ parameter_1, scales = "free", space = "free") + hex_theme +
    scale_fill_gradientn(name = "Counts", colours = c("white", "darkgoldenrod1", "chocolate1", "firebrick1", "black"), trans = "log10", breaks = c(1,3,10,30,100)) +
    scale_x_continuous(expression(log[2]), expand = c(0, 0), breaks = seq(-20, 20, by = 5)) + scale_y_continuous(NULL, expand = c(0, 0), labels = density_trans, breaks = density_trans_inv(seq(-20, 20, by = 5))) +
    geom_point(data = MLEpoints, aes(x = xval, y = yval), size = 4, col = "cornflowerblue") +
    geom_vline(data = MLEbarplot, aes(xintercept = xval), col = "cornflowerblue", size = 2) +  geom_bar(data = par_comp_like, aes(x = xval), binwidth = par_hist_binwidth, col = "black") +
    geom_blank(data = bound_expand, aes(x = value, y = value_trans))
  
  quantiles <- par_comp_like %>% dplyr::select(parameter_1, xval) %>% group_by(parameter_1) %>% summarize(LH = quantile(xval, probs = c(0.025)), UH = quantile(xval, probs = c(0.975))) %>%
    gather(bound, xval, -parameter_1)
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray92"), legend.position = "top", 
                       axis.ticks = element_line(color = "black", size = 1),
                       axis.text = element_text(color = "black", size = 20),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank(), panel.margin = unit(1.5, "lines")
                       )
  
  output_plots[["univariateHist"]] <- ggplot() + geom_vline(data = quantiles, aes(xintercept = xval), color = "blue", size = 2) + geom_vline(data = MLEbarplot, aes(xintercept = xval), color = "red", size = 2) +
    geom_bar(data = par_comp_like, aes(x = xval), col = "black", binwidth = 0.5) + 
    geom_blank(data = bound_expand, aes(x = value)) +
    geom_blank(data = par_comp_like, aes(x = xval, y=1.05*..count..), stat="bin", binwidth = 0.5) +
    facet_grid(~ parameter_1, scale = "free_x", space = "free_x") + barplot_theme +
    scale_y_continuous(NULL, expand = c(0, 0), breaks = seq(0, 1000, by = 50)) +
    scale_x_continuous(expression(log[2]), expand = c(0, 0), breaks = seq(-20, 20, by = 5))
  
  return(output_plots)
  
}

metabolic_leverage_summary_plots <- function(nutrient_cond = "P0.05"){
  
  # Generate 3 plots for a nutrient condition and dump ggsave pdfs #
  # Summarize metabolic leverage - divide fraction of leverage by reaction into substrate, product and enzyme
  # Label the fraction of variance explained by the parameteric fit based on relative metabolic leverage and identify the unexplained fraction of variance
  # Label the CI overlap of the parameteric fit and FBA flux based on relative metabolic leverage and identify the unexplained fraction of variance  
  # For the latter two plots, include both the RM leverage and the addition of the most significant putatively known regulator (if one exists) 
  
  # if nutrient condition == median - use median ML, otherwise provide nutrient condition
  
  if(nutrient_cond == "median"){
    MLdata_reduced <- MLdata[, list(q0.5 = median(get("0.5"))), by = c("specie", "Type", "reaction")]
  }else if(nutrient_cond %in% unique(MLdata$condition)){
    MLdata_reduced <- MLdata[condition == nutrient_cond, list(q0.5 = get("0.5")), by = c("specie", "Type", "reaction")]
  }else{
    stop("invalid nutrient condition")
  }
  
  ##### Summarize weighted-elasticities / metabolic leverage #####
  MLsummary <- MLdata_reduced[reaction %in% reactionInfo$rMech[reactionInfo$form == "rm" & reactionInfo$modification == ""],]
  MLsummary <- MLsummary[substr(reaction, 1, 6) %in% valid_rxns,] # only look at well defined reactions
  
  MLsummary <- MLsummary[,list(leverage = sum(q0.5)), by = c("Type", "reaction")]
  MLsummary$Type <- factor(MLsummary$Type, levels = c("substrate", "product", "enzyme"))
  MLsummary[,totalLeverage := sum(leverage), by = "reaction"]
  MLsummary[,leverage := leverage/totalLeverage]
  MLsummary$VarExplained <- rxn_fits$parPearson[chmatch(MLsummary$reaction, rxn_fits$rxn)]
  
  enz_leverage <- data.frame(reaction = MLsummary[Type == "enzyme",reaction], rank = NA)
  enz_leverage$rank[order(MLsummary[Type == "enzyme",leverage])] <- 1:nrow(enz_leverage)
  
  MLsummary$rank <- factor(enz_leverage$rank[chmatch(MLsummary$reaction, enz_leverage$reaction)])
  setkeyv(MLsummary, c("Type", "rank"))
  setkey(MLsummary)
  
  ggplot(MLsummary, aes(x = rank, y = leverage, fill = Type)) + geom_bar(stat = "identity", width = 0.85) +
    barplot_theme + scale_y_continuous(expression('Metabolic Leverage: ' ~ frac("|"~epsilon[i]~"|"~sigma[i], sum("|"~epsilon[j]~"|"~sigma[j], "j = 1" , n))), expand = c(0,0)) +
    scale_fill_brewer(palette = "Set1") + scale_x_discrete("Reactions")
  ggsave("Figures/metabolicLeverage.pdf", height = 6, width = 10)
  
  
  ###### Summarize metabolic leverage ######
  # Reversible MM kinetics and RM + best significant regulator
  # Indicate fraction of missing explanatory power
  
  # consider reactions with a positive correlation between the rm form and flux carried
  
  
  MLsummary <- MLdata_reduced[reaction %in% optimal_rxn_form[rxn_fits$parPearson[chmatch(optimal_rxn_form, rxn_fits$rxn)] > 0],,] # optimal kinetic form and pearson corr > 0
  MLsummary <- MLsummary[,list(leverage = sum(q0.5)), by = c("Type", "reaction")]
  MLsummary$Type[MLsummary$Type %in% c("allosteric", "noncompetitive", "uncompetitive", "competitive", "mm", "cc")] <- "regulatory"
  
  # partitioning explained variance into ML contributions - keeping relative fractions proportional but adjusting the sum
  # MLfrac_species = ML * r2
  # MLfrac_residual = sum(ML) * (1-r2)
  # divide both by r2
  # MLfrac_species = ML
  # MLfrac_residual = sum(ML) * (1-r2)/r2
  
  residualLeverage <- MLsummary[,list(Type = "residual", leverage = sum(leverage) * (1-rxn_fits$parPearson[chmatch(reaction, rxn_fits$rxn)]^2)/rxn_fits$parPearson[chmatch(reaction, rxn_fits$rxn)]^2)
                                , by = "reaction"]
  
  MLsummary <- rbind(MLsummary, residualLeverage, use.names = T)
  
  MLsummary$Type <- factor(MLsummary$Type, levels = c("enzyme", "substrate", "product", "regulatory", "residual"))
  MLsummary[,totalLeverage := sum(leverage), by = "reaction"]
  MLsummary[,leverage := leverage/totalLeverage]
  
  MLsummary[,rxn := reactionInfo$reaction[reactionInfo$rMech == reaction], by = "reaction"]
  MLsummary[,rxnType := ifelse(any(Type %in% "regulatory"), "Regulator", "RM"), by = "reaction"]
  MLsummary$rxnType <- factor(MLsummary$rxnType, levels = c("RM", "Regulator"))
  
  
  
  rank_var_explained <- data.frame(reaction = unique(MLsummary$rxn), rank = NA)
  rank_var_explained$residLev <- sapply(rank_var_explained$reaction, function(x){min(MLsummary$leverage[MLsummary$Type == "residual" & MLsummary$rxn == x], na.rm = T)})
  rank_var_explained$rank[order(rank_var_explained$residLev, decreasing = T)] <- 1:nrow(rank_var_explained)
  
  
  MLsummary$rxn <- factor(MLsummary$rxn, levels = rank_var_explained$reaction[order(rank_var_explained$rank)])
  MLsummary <- MLsummary[order(MLsummary$rxn, MLsummary$rxnType, MLsummary$Type),]
  MLsummary$reaction <- factor(MLsummary$reaction, levels = unique(MLsummary$reaction))
  
  TypeColors <- c("firebrick1", "blue2", "chartreuse3", "darkorange", "gray90")
  names(TypeColors) <- levels(MLsummary$Type)
  
  ggplot(MLsummary, aes(x = reaction, y = leverage, fill = Type)) + geom_bar(stat = "identity", width = 0.85) +
    barplot_theme + scale_y_continuous(expression('Metabolic Leverage: ' ~ frac("|"~epsilon[i]~"|"~sigma[i], sum("|"~epsilon[j]~"|"~sigma[j], "j = 1" , n))), expand = c(0,0)) +
    scale_fill_manual(values = TypeColors)
  ggsave("Figures/metabolicLeverage_variance.pdf", height = 8, width = 12)
  
  
  
  ##### Same summary with capture probability #####
  
  MLsummary <- MLdata_reduced[reaction %in% optimal_rxn_form[fraction_flux_deviation$"Interval Overlap"[chmatch(optimal_rxn_form, fraction_flux_deviation$rxn)] > 0],,]
  MLsummary <- MLsummary[,list(leverage = sum(q0.5)), by = c("Type", "reaction")]
  MLsummary$Type[MLsummary$Type %in% c("allosteric", "noncompetitive", "uncompetitive", "competitive", "mm", "cc")] <- "regulatory"
  
  residualLeverage <- MLsummary[,list(Type = "residual", leverage = sum(leverage) * (1-fraction_flux_deviation$"Interval Overlap"[chmatch(reaction, fraction_flux_deviation$rxn)])/fraction_flux_deviation$"Interval Overlap"[chmatch(reaction, fraction_flux_deviation$rxn)])
                                , by = "reaction"]
  
  MLsummary <- rbind(MLsummary, residualLeverage, use.names = T)
  
  
  MLsummary$Type <- factor(MLsummary$Type, levels = c("enzyme", "substrate", "product", "regulatory", "residual"))
  MLsummary[,totalLeverage := sum(leverage), by = "reaction"]
  MLsummary[,leverage := leverage/totalLeverage]
  
  MLsummary[,rxn := reactionInfo$reaction[reactionInfo$rMech == reaction], by = "reaction"]
  MLsummary[,rxnType := ifelse(any(Type %in% "regulatory"), "Regulator", "RM"), by = "reaction"]
  MLsummary$rxnType <- factor(MLsummary$rxnType, levels = c("RM", "Regulator"))
  
  rank_var_explained <- data.frame(reaction = unique(MLsummary$rxn), rank = NA)
  rank_var_explained$residLev <- sapply(rank_var_explained$reaction, function(x){min(MLsummary$leverage[MLsummary$Type == "residual" & MLsummary$rxn == x])})
  rank_var_explained$rank[order(rank_var_explained$residLev, decreasing = T)] <- 1:nrow(rank_var_explained)
  
  
  MLsummary$rxn <- factor(MLsummary$rxn, levels = rank_var_explained$reaction[order(rank_var_explained$rank)])
  MLsummary <- MLsummary[order(MLsummary$rxn, MLsummary$rxnType, MLsummary$Type),]
  MLsummary$reaction <- factor(MLsummary$reaction, levels = unique(MLsummary$reaction))
  
  TypeColors <- c("firebrick1", "blue2", "chartreuse3", "darkorange", "gray90")
  names(TypeColors) <- levels(MLsummary$Type)
  
  
  ggplot(MLsummary, aes(x = reaction, y = leverage, fill = Type)) + geom_bar(stat = "identity", width = 0.85) +
    barplot_theme + scale_y_continuous(expression('Metabolic Leverage: ' ~ frac("|"~epsilon[i]~"|"~sigma[i], sum("|"~epsilon[j]~"|"~sigma[j], "j = 1" , n))), expand = c(0,0)) +
    scale_fill_manual(values = TypeColors)
  ggsave("Figures/metabolicLeverage_intervalOverlap.pdf", height = 8, width = 12)

  ### Using interval_overlap summary determine reactions with > 50% CI capture - look at the ML sources of these reactions ###
  
  MLsummary <- MLsummary[Type != "residual" & reaction %in% MLsummary$reaction[MLsummary$Type == "residual" & MLsummary$leverage < 0.5],]
  MLsummary[,totalLeverage := sum(leverage), by = "reaction"]
  MLsummary[,leverage := leverage/totalLeverage]
  
  rank_var_explained <- data.frame(reaction = unique(MLsummary$rxn), rank = NA)
  rank_var_explained$enzML <- sapply(rank_var_explained$reaction, function(x){min(MLsummary$leverage[MLsummary$Type == "enzyme" & MLsummary$rxn == x])})
  rank_var_explained$rank[order(rank_var_explained$enzML, decreasing = T)] <- 1:nrow(rank_var_explained)
  
  MLsummary$rxn <- factor(MLsummary$rxn, levels = rank_var_explained$reaction[order(rank_var_explained$rank)])
  MLsummary <- MLsummary[order(MLsummary$rxn, MLsummary$rxnType, MLsummary$Type),]
  MLsummary$reaction <- factor(MLsummary$reaction, levels = unique(MLsummary$reaction))
  
  ggplot(MLsummary, aes(x = reaction, y = leverage, fill = Type)) + geom_bar(stat = "identity", width = 0.85) +
    barplot_theme + scale_y_continuous(expression('Metabolic Leverage: ' ~ frac("|"~epsilon[i]~"|"~sigma[i], sum("|"~epsilon[j]~"|"~sigma[j], "j = 1" , n))), expand = c(0,0)) +
    scale_fill_brewer(palette = "Set1") + scale_x_discrete("Reactions")
  ggsave("Figures/metabolicLeverage_goodFit.pdf", height = 6, width = 10)
  
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
  
  nameLengths <- data.frame(name = names, nchar = sapply(names, function(x){length(strsplit(x, "")[[1]])}))
  
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

########## Functions used for summarizing behavior of sets of reactions ############

filter_reactions <- function(reactionInfo, rxn_fits, rxnList_form){
  
  # Some fits were tested despite missing substrates, for these reactions,
  # we don't want a regulator substituting for missing substrates
  
  reaction_quality <- NULL
  
  for(a_reaction  in unique(reactionInfo$reaction)){
    
    rMM_form <- reactionInfo$rMech[reactionInfo$reaction == a_reaction & reactionInfo$modelType == "rMM" & reactionInfo$ncond == max(reactionInfo$ncond)]
    
    rMM_data <- rxnList_form[[rMM_form]]
    # test to see whether any two metabolites use the same data
    
    metAbundances <- rMM_data$rxnMet[,colSums(is.na(rMM_data$rxnMet)) == 0, drop = F]
    metAbundances <- metAbundances - matrix(apply(metAbundances, 2, mean), ncol = ncol(metAbundances), nrow = nrow(metAbundances), byrow = T)
    
    if(any(c(dist(t(metAbundances))) < 10^-12)){
      reaction_quality <- rbind(reaction_quality, data.frame(reaction = a_reaction, missing = NA, include = F, reason = "shared measurement"))
      next
    }
    
    reversibility <- ifelse(reversibleRx$modelBound[reversibleRx$rx == a_reaction] == "reversible", "reversible", "irreversible")
    
    reaction_stoi <- data.frame(tID = names(rMM_data$rxnStoi), stoi = unname(rMM_data$rxnStoi))
    if(all(rMM_data$flux$standardQP >= 0)){
      reaction_stoi$role = ifelse(reaction_stoi$stoi < 0, "substrate", "product")
    }else if(all(rMM_data$flux$standardQP <= 0)){
      reaction_stoi$role = ifelse(reaction_stoi$stoi > 0, "substrate", "product")
    }else{
      reaction_stoi$role = "substrate"
    }
    
    reaction_stoi <- reaction_stoi %>% left_join(data.frame(tID = names(rMM_data$originMet), measured = unname(rMM_data$originMet)), by = "tID")
    reaction_stoi <- reaction_stoi %>% left_join(data.frame(tID = names(rMM_data$metNames), name = unname(rMM_data$metNames)), by = "tID")
    
    # freebie metabolites
    measure_exception <- c("H+", "H2O", "ammonium")
    reaction_stoi <- reaction_stoi %>% filter(!(name %in% measure_exception))
    
    # Major species measured
    if(reversibility == "irreversible" & all(reaction_stoi$measured[reaction_stoi$role == "substrate"] != "nm")){
      reaction_quality <- rbind(reaction_quality, data.frame(reaction = a_reaction, missing = sum(reaction_stoi$measured == "nm"), include = T, reason = "irreversible and measured substrates"))
      next
    }
    if(reversibility == "reversible" & all(reaction_stoi$measured != "nm")){
      reaction_quality <- rbind(reaction_quality, data.frame(reaction = a_reaction, missing = 0, include = T, reason = "reversible and all species measured"))
      next
    }
    # No substrates are measured
    if(all(reaction_stoi$measured[reaction_stoi$role == "substrate"] == "nm")){
      reaction_quality <- rbind(reaction_quality, data.frame(reaction = a_reaction, missing = sum(reaction_stoi$measured[reaction_stoi$role == "substrate"] == "nm"), include = F, reason = "all substrates missing"))
      next
    }
    
    # Some relevant species missing - either the fit is still good but needs to be checked or the reaction should be discarded
    
    fit_rMechs <- reactionInfo %>% filter(reaction == a_reaction) %>% filter(modelType %in% c("rMM", "regulator", "cooperativity", "2+ regulators")) %>%
      left_join(rxn_fits %>% dplyr::select(rMech = rxn, spearman = parSpearman), by = "rMech") %>% filter(is.na(Qvalue) | Qvalue < 0.05)
    
    if(any(fit_rMechs$spearman > 0.3)){
      reaction_quality <- rbind(reaction_quality, data.frame(reaction = a_reaction, missing = sum(reaction_stoi$measured == "nm"), include = NA, reason = "some species missing: manually inspect"))
    }else{
      reaction_quality <- rbind(reaction_quality, data.frame(reaction = a_reaction, missing = sum(reaction_stoi$measured == "nm"), include = F, reason = "missing species and all rMech poorly fit"))
    }
  }
  
  reaction_quality %>% filter(!is.na(include) & !include)
  
  suspect_reactions <- reaction_quality %>% filter(is.na(include))
  
  suspect_rMechs <- reactionInfo %>% filter(reaction %in% suspect_reactions$reaction) %>% filter(modelType %in% c("rMM", "regulator", "cooperativity", "2+ regulators")) %>%
    filter(is.na(Qvalue) | Qvalue < 0.05)
  
  # For each reaction where the best reaction form is Michaelis-Menten, either this form is sufficient
  # despite missing a substrate or the missing substrate is likely not particularly important
  
  #ML_suspect_rMechs <- suspect_rMechs %>% group_by(reaction) %>% dplyr::summarize(bestModel = modelType[which.max(ML)])
  #ML_suspect_rMechs <- ML_suspect_rMechs$reaction[ML_suspect_rMechs$bestModel == "rMM"]
  
  #suspect_reactions %>% filter(!(reaction %in% ML_suspect_rMechs))
  
  missing_major_substrates <- c("r_0534", # Hexokinase, no glucose
                                "r_0799", # no dTDP
                                "r_0800", # no GDP
                                "r_0818", # no N-acetyl ornithine substrate
                                "r_0988", # Saccharopine DH, no saccharopine
                                "r_1055", # missing major substrate
                                "r_0510") # fluxes are suspect
  
  reaction_quality$include[reaction_quality$reaction %in% missing_major_substrates] <- FALSE
  reaction_quality$include[is.na(reaction_quality$include)] <- TRUE
  
  return(reaction_quality)
  
}

regulation_lit_support <- function(relevant_rxns, all_reactionInfo){
  
  # Load BRENDA regulation and gold-standard regulation and compare to supported
  # reaction forms with a single regulator
  
  # Test
  # How often does gold-standard regulation improve fit relative to other regulation
  # Is p(signif) related to literature representation
  # - output model fitting p(signif)
  # - output tbl_df() with each regulator * reaction * inhibitor/activator
  # -- summarize literature support
  # -- fitted probability of significance given literature and GS
  
  require(dplyr)
  
  return_list <- list()
  
  all_affinities <- read.delim("flux_cache/metaboliteAffinities.tsv")
  all_regulation <- all_affinities %>% tbl_df() %>% filter(speciesType == "regulator")
  
  all_regulation_unfold <- apply(all_regulation, 1, function(x){
    x <- unlist(x)
    cbind(expand.grid(tID = strsplit(x['tID'], split = '/')[[1]], reaction = strsplit(x['reactions'], split = '/')[[1]], stringsAsFactors = F),
          t(x[c('EC', 'isYeast', 'modtype', 'nQual')]))
  })
  all_regulation <- do.call("rbind", all_regulation_unfold) %>% tbl_df() %>%
    mutate(nQual = as.numeric(nQual),
           isYeast = gsub(' ', '', isYeast), isYeast = as.logical(isYeast),
           isYeast = ifelse(isYeast, "sce", "other"))
  
  all_regulation <- all_regulation %>% group_by(EC, reaction, tID, isYeast, modtype) %>% dplyr::summarize(nQual = sum(nQual)) %>%
    spread(key = isYeast, value = nQual, fill = 0)
  
  # If records from multiple sources exist (due to a reaction with multiple E.C. #s), take the best annotated reaction
  all_regulation <- all_regulation %>% group_by(reaction, tID, modtype) %>% dplyr::summarize(other = max(other), sce = max(sce))
  
  # remove rmCond reactions, only look at non-pathological reactions (fit reasonably using some tested model)
  
  tested_regulation <- all_reactionInfo %>% filter(modelType == "regulator") %>% tbl_df() %>%
    group_by(reaction) %>% filter(ncond == min(ncond)) %>% ungroup() %>%
    filter(reaction %in% relevant_rxns) %>%
    mutate(modtype = ifelse(grepl('-inh-', modification), 'inh', 'act'),
           tID = regmatches(modification, regexpr('t_[0-9]{4}', modification)))
  
  # if multiple models of activaton or inhibition are tested, just look at the best
  tested_regulation <- tested_regulation %>% group_by(reaction, tID, modtype) %>% filter(ML == max(ML))
  
  tested_regulation <- tested_regulation %>% dplyr::select(rMech, reaction, tID, modtype, Qvalue, FullName) %>%
    mutate(is_sig = ifelse(Qvalue < 0.05, T, F))
  
  tested_regulation <- tested_regulation %>% left_join(all_regulation, by = c("reaction", "tID", "modtype")) %>%
    mutate(other = ifelse(is.na(other), 0, other),
           sce = ifelse(is.na(sce), 1, sce))
  
  # look at gold standard validated regulation #
  GS_regulation <- read.delim("companionFiles/gold_standard_regulation.txt") %>%
    tbl_df() %>% mutate(modtype = ifelse(type == "inhibitor", "inh", type),
                        modtype = ifelse(type == "activator", "act", modtype)) %>%
    mutate(is_GS = T) %>% dplyr::select(tID, reaction, modtype, is_GS)
  
  tested_regulation <- tested_regulation %>% left_join(GS_regulation, by = c("tID", "reaction", "modtype")) %>%
    mutate(is_GS = ifelse(is.na(is_GS), F, is_GS))
  
  #
  GS_contingency <- table(SIG = tested_regulation$is_sig, GS = tested_regulation$is_GS)
  
  # GS annotated regulation is more likely to be consistent
  #(GS_contingency[2,2]/sum(GS_contingency[,2])) / (GS_contingency[2,1]/sum(GS_contingency[,1]))
  #set.seed(1234)
  #chisq.test(GS_contingency, simulate.p.value = T, B = 1e7)
  
  ### Collapse highly correlated entries ###
  corr_penalty <- 0.6 # subtract from correlation s.t. grouping is only favored when cor > corr_penalty
  
  regulator_trends <- lapply(tested_regulation$rMech, function(x){
    
    regs <- rxnList_form[[x]]$rxnFormData %>% filter(Type %in% c("inh", "act")) %>% dplyr::select(SubstrateID, Type)
    regs <- regs %>% left_join(data.frame(SubstrateID = names(rxnList_form[[x]]$metNames), commonName = unname(rxnList_form[[x]]$metNames)), by = "SubstrateID")
    
    rxnList_form[[x]]$rxnMet[,colnames(rxnList_form[[x]]$rxnMet) %in% regs$SubstrateID, drop = F] %>%
      mutate(condition = rownames(.)) %>% gather("SubstrateID", "RA", -condition, convert = T) %>%
      left_join(regs, by = "SubstrateID") %>% mutate(rMech = x)
    
  })
  regulator_trends <- do.call("rbind", regulator_trends)
  regulator_trends <- regulator_trends %>% left_join(tested_regulation %>% dplyr::select(rMech, reaction), by = "rMech") %>%
    dplyr::select(condition, tID = SubstrateID, RA) %>% unique()
  
  met_cluster <- NULL
  
  for(a_reaction in unique(tested_regulation$reaction)){
    for(a_reg_type in c("inh", "act")){
      for(signif in c(T, F)){
        
        reg_subset <- tested_regulation %>% filter(reaction == a_reaction, modtype == a_reg_type, is_sig == signif)
        
        if(nrow(reg_subset) == 0){
          # no regulators in class
          next 
        }else if(nrow(reg_subset) == 1){
          # only one metabolite
          met_cluster <- rbind(met_cluster, data.frame(reaction = a_reaction, modtype = a_reg_type, tID = reg_subset$tID, group = 1))
          next
          
        }
        
        regulator_RA <- regulator_trends %>% filter(tID %in% reg_subset$tID) %>%
          spread(tID, value = RA)
        
        regulator_corr <- regulator_RA %>% dplyr::select(-condition) %>% as.matrix() %>% cor()
        regulator_corr <- regulator_corr - corr_penalty
        
        all_species <- rownames(regulator_corr)
        nclus <- length(all_species)
        clus_assignments <- data.frame(ID = all_species, cluster = 1:nclus)
        
        has_converged = F
        iter <- 0
        old_score <- -Inf
        
        while(!has_converged){
          
          iter <- iter + 1
          
          for(i in 1:nclus){
            current_assignments <- clus_assignments %>% dplyr::slice(-i)
            
            cluster_score <- rbind(current_assignments, data.frame(ID = all_species[i], cluster = 1:nclus)) %>%
              group_by(cluster) %>% summarize(score = 
                                                mean(regulator_corr[rownames(regulator_corr) %in% ID, colnames(regulator_corr) %in% ID][
                                                  lower.tri(regulator_corr[rownames(regulator_corr) %in% ID, colnames(regulator_corr) %in% ID], diag = F)]))
            cluster_score$score[is.na(cluster_score$score)] <- 0
            
            clus_assignments$cluster[i] <- cluster_score$cluster[which.max(cluster_score$score)]
          }
          
          current_score <- clus_assignments %>% group_by(cluster) %>%
            summarize(score = mean(regulator_corr[rownames(regulator_corr) %in% ID, colnames(regulator_corr) %in% ID][
              lower.tri(regulator_corr[rownames(regulator_corr) %in% ID, colnames(regulator_corr) %in% ID], diag = F)]))
          current_score$score[is.na(current_score$score)] <- 0
          
          current_score <- current_score %>% dplyr::select(score) %>% unlist() %>% unname() %>% sum()
          
          if(old_score >= current_score){
            has_converged <- T
          }else{
            old_score <- current_score
          }
        }
        
        met_cluster <- rbind(met_cluster, data.frame(reaction = a_reaction, modtype = a_reg_type, tID = clus_assignments$ID, group = clus_assignments$cluster))
      }
    }
  }
  
  tested_regulation <- tested_regulation %>% left_join(met_cluster, by = c("reaction", "modtype", "tID"))
  
  tested_regulation_aggregate <- tested_regulation %>% group_by(reaction, modtype, is_sig, group) %>%
    dplyr::summarize(N = n(), is_GS = sum(is_GS), sce = sum(sce), other = sum(other))
  
  tested_regulation_aggregate <- tested_regulation_aggregate %>% dplyr::select(-group) %>%
    rowwise() %>% mutate(sce_mean = sce/N, other_mean = other/N) %>%
    group_by(reaction) %>% mutate(totalCitations = sum(sce_mean) + sum(other_mean), 
                                  totalYeast = sum(sce_mean),
                                  regTests = sum(N)) %>%
    rowwise() %>% mutate(sce_frac = ifelse(totalYeast != 0, sce_mean/totalYeast, 0),
                         other_frac = other_mean/totalCitations)
  
  tested_regulation_aggregate <- tested_regulation_aggregate %>%
    mutate(sce_factor = ifelse(sce == 0, "None", "NA"),
           sce_factor = ifelse(sce > 3, "High", sce_factor),
           sce_factor = ifelse(sce_factor == "NA", "Noted", sce_factor),
           sce_factor = factor(sce_factor, levels = c("None", "Noted", "High"), ordered = T)) %>%
    mutate(other_factor = ifelse(other > 20 | other_frac > 0.3, "High", "NA"),
           other_factor = ifelse(other <= 2 | other_frac < 0.04, "Low", other_factor),
           other_factor = ifelse(other_factor == "NA", "Mid", other_factor),
           other_factor = factor(other_factor, levels = c("Low", "Mid", "High"), ordered = T))
  
  ### Hypothesis testing
  
  weightedReg_withGS = glm(is_sig ~ modtype + sce_frac + sce_mean + other_frac + other_mean + is_GS,
                           weight = N, family=binomial(logit), data = tested_regulation_aggregate)
  
  weightedReg_withGS = glm(is_sig ~ modtype + regTests + sce_mean + other_mean,
                           weight = N, family=binomial(logit), data = tested_regulation_aggregate)
  
  
  
  
  
  
  
  
  sig_assoc_norm <- sig_association %>% group_by(reaction) %>% mutate(totalObs = sum(sce + other), totalYeast = sum(sce)) %>%
    rowwise() %>% mutate(sce_frac = ifelse(totalYeast != 0, sce/totalYeast, 0), other_frac = other/totalObs)
  
  sig_assoc_norm <- sig_assoc_norm %>%
    mutate(sce_factor = ifelse(sce == 0, "None", "NA"),
           sce_factor = ifelse(sce > 3 & sce_frac > 0.3, "High", sce_factor),
           sce_factor = ifelse(sce_factor == "NA", "Noted", sce_factor),
           sce_factor = factor(sce_factor, levels = c("None", "Noted", "High"))) %>%
    mutate(other_factor = ifelse(other > 20 | other_frac > 0.3, "High", "NA"),
           other_factor = ifelse(other <= 2 | other_frac < 0.07, "Low", other_factor),
           other_factor = ifelse(other_factor == "NA", "Mid", other_factor),
           other_factor = factor(other_factor, levels = c("Low", "Mid", "High")))
  
  # table(sig_assoc_norm$is_sig, sig_assoc_norm$sce != 0)
  
  
  # standard models
  #glm(is_sig ~ modtype + sce_frac + other_frac, family=binomial(logit), data=sig_assoc_norm)
  #glm(is_sig ~ modtype + sce_frac + other_frac + totalObs, family=binomial(logit), data=sig_assoc_norm)
  #glm(is_sig ~ modtype + sce + other + totalObs, family=binomial(logit), data=sig_assoc_norm)
  #glm(is_sig ~ modtype*sce_frac + modtype*other_frac + modtype*totalObs, family=binomial(logit), data=sig_assoc_norm)
  
  # weighted models for significance
  
  #weightedReg_noGS = glm(is_sig ~ modtype + sce_frac + other_frac + totalObs,
  #    weights = (other_frac + sce_frac)*totalObs, family=binomial(logit), data=sig_assoc_norm)
  
  weightedReg_withGS = glm(is_sig ~ modtype + sce_frac + other_frac,
                           family=binomial(logit), data=sig_assoc_norm)
  
  weightedReg_withGS = glm(is_sig ~ modtype + sce_factor + other_factor,
                           family=binomial(logit), data=sig_assoc_norm)
  summary(weightedReg_withGS)
  
  weightedReg_withGS = glm(is_sig ~ modtype + sce_factor + other_factor + is_GS,
                           family=binomial(logit), data=sig_assoc_norm)
  summary(weightedReg_withGS)
  
  weightedReg_withGS = glm(is_sig ~ modtype + sce_factor + other_factor + is_GS,
                           family=binomial(logit), data=sig_assoc_norm)
  
  
  weightedReg_withGS = glm(is_sig ~ modtype + sce_frac + other_frac + is_GS + totalObs,
                           family=binomial(logit), data=sig_assoc_norm)
  
  sig_assoc_norm <- sig_assoc_norm %>% ungroup() %>% mutate(prob_sig = unname(weightedReg_withGS$fitted))
  
  #ggplot(sig_assoc_norm, aes(x = prob_sig)) + facet_grid(is_sig ~., scale = "free_y") + geom_bar(binwidth = 0.02)
  #ggplot(sig_assoc_norm, aes(x = is_sig, y = log(prob_sig/(1 - prob_sig)))) + geom_boxplot(notch = T)
  
  return_list[['GS_contingency']] = GS_contingency
  return_list[['fitted_model']] = weightedReg_withGS
  return_list[['prob_reg']] = sig_assoc_norm
  
  return(return_list)
  
}

filter_rMech_by_prior <- function(relevant_rxns, reactionInfo, literature_support, rxnList_form){
  
  # Compare the relative support of models for each reaction
  
  require(tidyr)
  require(dplyr)
  
  # For each reaction, compare rMM with models of regulation penalizing regulatory models
  # based on the BRENDA-informed plausibility
  
  rMechs_considered <- reactionInfo %>% filter(reaction %in% relevant_rxns, modelType %in% c("rMM", "regulator", "cooperativity", "2+ regulators")) %>% tbl_df()
  
  regulators <- lapply(1:nrow(rMechs_considered), function(i){
    
    regs <- rxnList_form[[rMechs_considered$rMech[i]]]$rxnFormData %>% filter(Type %in% c("inh", "act")) %>% dplyr::select(SubstrateID, Type)
    regs <- regs %>% left_join(data.frame(SubstrateID = names(rxnList_form[[rMechs_considered$rMech[i]]]$metNames), commonName = unname(rxnList_form[[rMechs_considered$rMech[i]]]$metNames)), by = "SubstrateID")
    
    if(nrow(regs) != 0){
      return(data.frame(rMech = rMechs_considered$rMech[i], reaction = rMechs_considered$reaction[i], regs))
    }
    
  })
  regulators <- do.call("rbind", regulators) %>% tbl_df()
  
  regulators <- regulators %>% mutate(tID_type = paste(SubstrateID, Type, sep = "-"),
           common_type = paste(commonName, ifelse(Type == "act", "+", "-"))) %>%
    dplyr::select(rMech, reaction, tID_type, common_type)
  
  regulator_prior <- literature_support %>% mutate(tID_type = paste(tID, modtype, sep = "-")) %>%
    dplyr::select(reaction, tID_type, prob_sig)
  
  # join regulators with literature annotations and model fit (and model d.o.f.)
  regulators <- regulators %>% left_join(regulator_prior, by = c("reaction", "tID_type")) %>%
    left_join(reactionInfo %>% dplyr::select(rMech, modelType, ncond, npar, ML), by = "rMech")
  
  # If multiple types of inhibitor or activation are included, collapse them to the best form
  # for combinatorial regulation this involves taking a single rMech combination
  
  similar_rMech_matrix <- regulators %>% dplyr::select(rMech, reaction, modelType, ncond, tID_type, ML) %>% mutate(present = 1) %>%
    spread(key = tID_type, value = present, fill = 0) 
  
  similar_rMech_matrix <- similar_rMech_matrix %>% dplyr::select(rMech, reaction, modelType, ncond, ML) %>%
    mutate(regulators = apply(as.data.frame(similar_rMech_matrix[,!(colnames(similar_rMech_matrix) %in% c("rMech", "reaction", "ncond", "ML"))]), 1, function(x){
      paste(x, collapse = "") 
    }))
  
  # the most likely rMech of each type of similar regulation
  similar_rMech_matrix <- similar_rMech_matrix %>% group_by(reaction, ncond, regulators) %>% filter(ML == max(ML))
  
  # for each regulatory mechanism create a prior relative to the null (i.e. p/(1-p))
  # if multiple regulators exist take the product over all species
  regulatory_model_relative_prior <- regulators %>% filter(rMech %in% similar_rMech_matrix$rMech) %>% group_by(rMech) %>%
    dplyr::summarize(relative_prior = prod(prob_sig/(1-prob_sig))) %>%
    mutate(relative_prior = ifelse(relative_prior > 1, 1, relative_prior))
  
  # Summarize p(model)*p(data|model) for each rMech
  # look at AICc
  
  support_summary <- rMechs_considered %>% filter(modelType == "rMM" | rMech %in% regulatory_model_relative_prior$rMech) %>%
    left_join(regulatory_model_relative_prior, by = "rMech") %>%
    mutate(relative_prior = ifelse(is.na(relative_prior), 1, relative_prior)) %>%
    dplyr::select(reaction, rMech, ncond, npar, ML, relative_prior)
  
  support_summary <- support_summary %>% mutate(ML_prior = ML + log(relative_prior),
                                                AICc = 2*npar - 2*ML_prior + 2*npar*(npar + 1)/(ncond - npar - 1))
  
  support_summary <- support_summary %>% group_by(reaction, ncond) %>%
    dplyr::mutate(minAIC = min(AICc),
                  rel_prob = exp((minAIC - AICc)/2),
                  AIC_prob = rel_prob/sum(rel_prob))
  
  return(support_summary)
  
  # Find the most likely model
  #support_summary %>% group_by(reaction) %>% filter(AIC_prob == max(AIC_prob)) %>% View()
  
  }

mode_of_regulation <- function(rMech_support, rxnList_form){
  
  ### For each reaction with 1+ significant regulators ###
  # use reaction stoichiometry and flux carried to assign regulation at feed-back, feed-forward or cross-pathway
  # return classify, number of steps, 
  
  require(dplyr)
  require(tidyr)
  require(igraph)
  
  # for each regualted rMech, point to its 
  
  regulators <- lapply(1:nrow(rMech_support), function(i){
    
    regs <- rxnList_form[[rMech_support$rMech[i]]]$rxnFormData %>% filter(Type %in% c("inh", "act")) %>% dplyr::select(SubstrateID, Type)
    regs <- regs %>% left_join(data.frame(SubstrateID = names(rxnList_form[[rMech_support$rMech[i]]]$metNames), commonName = unname(rxnList_form[[rMech_support$rMech[i]]]$metNames)), by = "SubstrateID")
    
    if(nrow(regs) != 0){
      return(data.frame(rMech = rMech_support$rMech[i], reaction = rMech_support$reaction[i], regs))
    }
    
  })
  regulators <- do.call("rbind", regulators) %>% tbl_df()
  
  unique_met_rxn_pairs <- regulators %>% dplyr::select(reaction, SubstrateID, commonName) %>% unique()
  
  ### Form a bipartite graph using reaction stoichiometry and measured flux ###
  
  load("flux_cache/yeast_stoi_directed.Rdata") # S: metabolism stoichiometry
  load("flux_cache/reconstructionWithCustom.Rdata") # metabolic reconstruction files
  carried_flux <- read.table("flux_cache/fluxCarriedSimple.tsv", header = T, sep = "\t")
  
  # reactions which carry flux and direction
  Used_reactions <- data.frame(reaction = rownames(carried_flux), F = rowSums(carried_flux > 0), R = rowSums(carried_flux < 0)) %>% gather("Dir", "N", -reaction) %>% tbl_df() %>% 
    filter(N != 0) %>% mutate(name = paste(reaction, Dir, sep = "_")) %>%
    filter(grepl('r_[0-9]{4}', reaction))
  
  S_carried <- S[,colnames(S) %in% Used_reactions$name]
  colnames(S_carried) <- sub('_[FR]$', '', colnames(S_carried))
  S_carried <- S_carried[rowSums(S_carried != 0) != 0,] # metabolites which are utilized
  
  
  # common species and cofactors should only be included when they have a biosynthetic (not-energetic) role 
  # remove protons, water, NAD(P)(H) and ammonium
  S_carried <- S_carried[!(rownames(S_carried) %in% corrFile$SpeciesID[grep('^H\\+|^H2O|^NAD|^ammonium|carbon dioxide|bicarbonate', corrFile$SpeciesName)]),] 
  
  rxnFile_fluxcarried <- rxnFile %>% filter(ReactionID %in% colnames(S_carried)) %>% tbl_df()
  
  s_ATP <- corrFile$SpeciesID[grep('^ATP ', corrFile$SpeciesName)]
  s_ADP <- corrFile$SpeciesID[grep('^ADP ', corrFile$SpeciesName)]
  s_AMP <- corrFile$SpeciesID[grep('^AMP ', corrFile$SpeciesName)]
  s_Pi <- corrFile$SpeciesID[grep('^phosphate ', corrFile$SpeciesName)]  
  s_PPi <- corrFile$SpeciesID[grep('^diphosphate ', corrFile$SpeciesName)]  
  
  # Cofactor use is limited to reactions where either the adenylate group is donated or the primary role is interconverting nucleotides
  
  # Reactions including ATP + ADP|AMP
  ATP_cofactors <- rxnFile_fluxcarried %>% group_by(ReactionID) %>% filter(any(s_ATP %in% Metabolite) &
                                                            (any(s_ADP %in% Metabolite)|any(s_AMP %in% Metabolite))) %>%
    dplyr::select(ReactionID) %>% unlist() %>% unname() %>% unique()
  
  # Rescue reaction with more 3+ N[MDT]P
  ATP_convert <- rxnFile_fluxcarried %>% group_by(ReactionID) %>% dplyr::summarize(N_Nuc = length(grep('[ACGTU][MDT]P ', MetName))) %>% filter(N_Nuc >= 3) %>%
    dplyr::select(ReactionID) %>% unlist() %>% unname()
  
  ATP_cofactors <- setdiff(ATP_cofactors, ATP_convert)
  
  # For reactions in ATP_cofactors, zero out ATP, ADP, AMP, Pi, PPi
  
  S_carried[rownames(S_carried) %in% c(s_ATP, s_ADP, s_AMP, s_Pi, s_PPi), colnames(S_carried) %in% ATP_cofactors] <- 0
  
  # remove phosphate as a substrate from reactions, will allow it to serve as a feedback inhibitor w/o
  # forming a complete path
  # pyrophosphate should behave as it is directed into phosphate
  
  S_carried[rownames(S_carried) %in% s_Pi,][S_carried[rownames(S_carried) %in% s_Pi,] < 0] <- 0
  
  # Also deal with transamination reactions - look for reactions where Gln -> Glu or Glu -> aKG
  # rescue nitrogen fixing reactions GS and GD based on inclusion of ammonium
  
  s_Glu <- corrFile$SpeciesID[grep('^L-glutamate ', corrFile$SpeciesName)]
  s_Gln <- corrFile$SpeciesID[grep('^L-glutamine ', corrFile$SpeciesName)]  
  s_aKG <- corrFile$SpeciesID[grep('^2-oxoglutarate', corrFile$SpeciesName)]  
  s_amm <- corrFile$SpeciesID[grep('^ammonium', corrFile$SpeciesName)]  
  
  transamination_reactions <- rxnFile_fluxcarried %>% group_by(ReactionID) %>%
    filter(any(s_Glu %in% Metabolite) &
             (any(s_Gln %in% Metabolite)|any(s_aKG %in% Metabolite)) &
             !any(s_amm %in% Metabolite)) %>% 
    dplyr::select(ReactionID) %>% unlist() %>% unname() %>% unique()
  
  S_carried[rownames(S_carried) %in% c(s_Glu, s_Gln, s_aKG), transamination_reactions] <- 0
  
  # Look at annotation of pathways, to help detect some distant interactions
  # Some reactions are not included or have a blank annotation
  
  rxn_pathways <- read.delim("./flux_cache/reactionPathways.tsv") %>%
    filter(reactionID %in% colnames(S_carried))
  
  rxn_pathways <- lapply(1:nrow(rxn_pathways), function(i){
    pathway <- strsplit(rxn_pathways$pathway[i], split = '__')[[1]]
    if(length(pathway) != 0){
      data.frame(reactionID = rxn_pathways$reactionID[i], pathway)
    }
  })
  rxn_pathways <- do.call("rbind", rxn_pathways)
  
  # Counts by pathway
  pathway_summary <- rxn_pathways %>% group_by(pathway) %>% dplyr::summarize(N = n()) %>% arrange(desc(N)) %>%
    filter(N > 5 & N < 30)
  
  rxn_pathways <- rxn_pathways %>% filter(pathway %in% pathway_summary$pathway)
  
  
  ### convert stoichiometric matrix to a directed bipartite graph ###
  
  S_graph <- melt(t(S_carried), varnames = c("source", "sink"))
  S_graph <- S_graph[S_graph$value != 0,] # reactions - metabolite links
  S_graph$source <- as.character(S_graph$source); S_graph$sink <- as.character(S_graph$sink) # class coercion
  
  S_graph[S_graph$value < 0,] <- S_graph[S_graph$value < 0,][,c(2,1,3)] # for consumed metabolites, invert direction
  S_graph <- S_graph[,-3]
  
  S_igraph <- graph.data.frame(S_graph)
  V(S_igraph)$type <- ifelse(substr(V(S_igraph)$name, 1, 1) == "r", "reaction", "metabolite")
  
  #sort(betweenness(S_igraph), decreasing = T)[1:4] # check which species have the highest betweeness to see if they should be dealt with
  #rxnFile_fluxcarried %>% group_by(ReactionID) %>% filter("s_0066" %in% Metabolite) %>% View()

  # Look at length of shortest path from all regulated reactions to all measured metabolites
  
  # Calculate paths from reaction to regulator and from regulator to reaction to define feedbacks and feedforward
  
  all_reg_paths <- lapply(1:nrow(unique_met_rxn_pairs), function(i){
    
    # because tIDs map to metabolites in multiple compartments, while reactions are unique
    # each compartment-specific sID is seperately linked to a reaction
    
    met_sID <- corrFile$SpeciesID[corrFile$SpeciesType == unique_met_rxn_pairs$SubstrateID[i]]
    met_sID <- met_sID[met_sID %in% V(S_igraph)$name]
    
    if(length(met_sID) == 0){
      print("external met")
      return(list())
      
    }
    if(unique_met_rxn_pairs$reaction[i] %in% V(S_igraph)$name == F){
      print("excluded reaction - shouldnt happen")
      return(list())
    }
    
    regulation_path <- list()
    
    # distance rID to sID
    feedbacks <- get.all.shortest.paths(S_igraph,
                                        which(V(S_igraph)$name == unique_met_rxn_pairs$reaction[i]),
                                        which(V(S_igraph)$name %in% met_sID))$res
    if(length(feedbacks) != 0){
      feedbacks <- melt(feedbacks)
      colnames(feedbacks) <- c("node", "path")
      feedbacks$pairnum <- i
      feedbacks$type <- "feedback"
      regulation_path <- rbind(regulation_path, feedbacks)
    }
    
    # distance sID to rID
    feedforward <- get.all.shortest.paths(S_igraph,
                                          which(V(S_igraph)$name %in% met_sID),
                                          which(V(S_igraph)$name == unique_met_rxn_pairs$reaction[i]))$res
    if(length(feedforward) != 0){
      feedforward <- melt(feedforward)
      colnames(feedforward) <- c("node", "path")
      feedforward$pairnum <- i
      feedforward$type <- "feedforward"
      regulation_path <- rbind(regulation_path, feedforward)
    }
   
    return(regulation_path)
    
  })
  all_reg_paths <- do.call("rbind", all_reg_paths) %>% tbl_df()
  
  vertex_info <- data.frame(node = 1:length(V(S_igraph)$name), name = V(S_igraph)$name, spec_type = V(S_igraph)$type, in_degree = unname(igraph::degree(S_igraph, mode = "in")),
                            out_degree = unname(igraph::degree(S_igraph, mode = "out")),
                            betweenness = unname(betweenness(S_igraph)))
  
  # Do not classify distant regulation as feedback/feedfoward, if there are many bifurcations
  vertex_info[,c('in_degree', 'out_degree')][vertex_info[,c('in_degree', 'out_degree')] == 0] <- 1
  
  pathway_class_call <- all_reg_paths %>% left_join(vertex_info, by = "node") %>% filter(spec_type == "reaction") %>%
    dplyr::select(path, pairnum, type, reactionID = name) %>% left_join(rxn_pathways, by = "reactionID")
  
  pathway_class_call <- pathway_class_call %>% group_by(path, pairnum, type, pathway) %>% dplyr::summarize(Npathway = n()) %>%
    left_join(pathway_class_call %>% group_by(path, pairnum, type) %>% dplyr::summarize(N = length(unique(reactionID))), by = c("path", "pairnum", "type"))
  
  # For each pathway, look for overrepresentation of a single pathway term
  pathway_class_call <- pathway_class_call %>% filter(!is.na(pathway)) 
  
  pathway_class_call <- pathway_class_call %>% mutate(pathway_match = ifelse((Npathway >= 0.5*N & N >= 5) | (Npathway == N & N >= 3), T, F))
  pathway_class_call <- pathway_class_call %>% group_by(path, pairnum, type) %>% dplyr::summarize(pathway_match = ifelse(any(pathway_match), T, F))
  
  # add back reactions with no pathway annotation for any reactions
  
  pathway_class_call <- rbind(pathway_class_call,  
                              all_reg_paths %>% dplyr::select(path, pairnum, type) %>% unique() %>% anti_join(pathway_class_call, by = c("path", "pairnum", "type")) %>% mutate(pathway_match = F)
  )
  
  # Also call paths where the reaction is very close to the regulator in the pruned network #
  
  regulator_betweenness <- all_reg_paths %>% left_join(vertex_info, by = "node") %>% filter(spec_type == "metabolite") %>%
    group_by(path, pairnum, type) %>% dplyr::summarize(met_betweenness = ifelse(type[1] == "feedback", last(betweenness), first(betweenness)))
  
  reg_class_call <- all_reg_paths %>% left_join(vertex_info, by = "node") %>% filter(spec_type == "reaction") %>%
    group_by(path, pairnum, type) %>% dplyr::summarize(N_steps = n(),
                                                       split_prod = ifelse(type[1] == "feedback", prod(1/in_degree),  prod(1/out_degree))) %>%
    left_join(regulator_betweenness, by = c("path", "pairnum", "type")) %>%
    left_join(pathway_class_call, by = c("path", "pairnum", "type")) %>%
    group_by(pairnum) %>% filter(N_steps == min(N_steps), split_prod == max(split_prod)) %>% dplyr::slice(1) %>%
    mutate(type = ifelse(split_prod < 0.2 & !pathway_match, "cross-pathway", type))
    
  # By regulation mechanism, regulation type call:
  all_regulation_type <- unique_met_rxn_pairs %>% mutate(pairnum = 1:n()) %>% left_join(reg_class_call %>% dplyr::select(pairnum, type, N_steps), by = "pairnum") %>%
    mutate(type = ifelse(is.na(type), 'cross-pathway', type))
  
  # For combinatorial regulation, collapse multiple entries into one-per-rMech
  rMech_class <- regulators %>% left_join(all_regulation_type %>% dplyr::select(reaction, SubstrateID, type, N_steps), by = c("reaction", "SubstrateID")) %>%
    group_by(rMech) %>% dplyr::summarize(type = paste(sort(unique(type)), collapse = "/"), N_steps = paste(sort(N_steps), collapse = ","))
  
  # add back unregulated reactions
  
  rMech_class <- rMech_class %>% mutate(N_steps = ifelse(N_steps == "", NA, N_steps)) %>% rbind(
    rMech_support %>% ungroup() %>% dplyr::select(rMech) %>% anti_join(rMech_class, by = "rMech") %>% mutate(type = "unregulated", N_steps = NA)
  )
  
  return(rMech_class)
  
  #reaction_regulation_type %>% spread("type", "relLik", fill = 0) %>% View()
  
  # look at distance between regulator and metabolites
  #reg_dist <- reg_class_call %>% group_by(type, N_steps) %>% dplyr::summarize(N = n())
  #max_reg_dist_density <- reg_dist %>% group_by(N_steps) %>% summarize(Nmax = sum(N))
  
  #S_dist <- data.frame(type = "overall", N_steps = 1:length(path.length.hist(S_igraph)$res), N = path.length.hist(S_igraph)$res) %>%
  #  dplyr::mutate(N_norm = N * max(max_reg_dist_density$Nmax)/max(N))
  
  #ggplot(reg_dist, aes(x = N_steps, y = N, fill = type)) + geom_bar(stat = "identity") +
  #  geom_line(data = S_dist, aes(x = N_steps, y = N_norm))
  
  # look at centrality of regulators
  
  #ggplot(reg_class_call, aes(x = met_betweenness, fill = type)) + geom_bar()
  #qplot(betweenness(S_igraph))
  
}







enzyme_control_source <- function(){
  
  ##### Associating enzyme metabolic leverage with transcription factors, thermodynamics ... #####
  # Look at the ML of enzymes to identify cases where:
  # A) Potential inducibility is high verus low
  # B) Potential inducibility varies based upon nutrient condition
  
  ML_inducibility <- MLdata[reaction %in% adequate_fit_optimal_rxn_form,,]
  ML_inducibility <- ML_inducibility[Type == "enzyme",,]
  ML_inducibility <- ML_inducibility[,list(ML = sum(get("0.5")), nenzyme = length(get("0.5"))), by = c("reaction", "condition")]
  
  # collapse across conditions to mean(ML) and SD(ML)
  ML_inducibility_summary <- ML_inducibility[,list(ML_mean = mean(ML), ML_sd = sd(ML), ML_min = min(ML), ML_max = max(ML), ML_range = max(ML)-min(ML), nenzyme = nenzyme[1]), by = "reaction"]
  ML_inducibility_summary[,CV := ML_sd/ML_mean] 
  ML_inducibility_summary$rxn <- substr(ML_inducibility_summary$reaction, 1, 6)
  ML_inducibility_summary[,ML_logit := log(ML_mean/(1-ML_mean))]
  
  rxn_meta_info <- read.delim("flux_cache/rxnParYeast.tsv")
  aligned_ML_meta <- rxn_meta_info[chmatch(ML_inducibility_summary$rxn, rxn_meta_info$ReactionID),]
  
  # Determine pathway E.C. numbers
  # If multiple E.C. identifiers exist, take the first one (this will be the identifier associated with the kegg R ID if there is one)
  aligned_ML_meta$EC <- sapply(aligned_ML_meta$EC, function(x){
    strsplit(x, split = ",")[[1]][1]
  })
  
  aligned_ML_meta$EC1 <- sapply(aligned_ML_meta$EC, function(x){strsplit(x, split = "\\.")[[1]][1]})
  aligned_ML_meta$EC2 <- sapply(aligned_ML_meta$EC, function(x){paste(strsplit(x, split = "\\.")[[1]][1:2], collapse = ".")})
  
  # Determine pathway-by-reaction
  
  aligned_ML_meta$pathname <- sapply(aligned_ML_meta$pathname, function(x){strsplit(x, split = "__")[[1]][1]})
  
  rxn_meta_path_prune <- melt(lapply(aligned_ML_meta$pathname, function(x){
    strsplit(x, split = "__")
  }))
  
  rxn_meta_path_pruned <- rxn_meta_path_prune[rxn_meta_path_prune$value %in% names(table(rxn_meta_path_prune$value))[table(rxn_meta_path_prune$value) >= 4 & table(rxn_meta_path_prune$value) <= 20],]
  rxn_meta_path_pruned$rxn <- aligned_ML_meta$ReactionID[rxn_meta_path_pruned$L1]
  rxn_meta_path_pruned <- rbind(rxn_meta_path_pruned, data.frame(value = "Misc", L2 = 1, L1 = NA, rxn = aligned_ML_meta$ReactionID[!(aligned_ML_meta$ReactionID %in% rxn_meta_path_pruned$rxn)]))
  rxn_pathways_cast <- acast(rxn_meta_path_pruned, formula = rxn ~ value, value.var = "L2", fill = 0)
  
  # Determine reaction reversibility
  aligned_ML_meta$reversibility <- ifelse(reversibleRx$modelBound[chmatch(aligned_ML_meta$ReactionID, reversibleRx$rx)] == "greaterEqual", "F", "T")
  
  # Regression of EC and pathway on metabolic leverage
  cat("\nEnzyme control vs. EC\n")
  print(anova(lm(ML_inducibility_summary$ML_logit ~ aligned_ML_meta$EC1))) # associate with E.C. number
  cat("\nEnzyme control vs. pathway\n")
  print(anova(lm(ML_inducibility_summary$ML_logit ~ rxn_pathways_cast))) # assocaite with pathway
  cat("\nEnzyme control vs. rxn reversibility\n")
  print(anova(lm(ML_inducibility_summary$ML_logit ~ aligned_ML_meta$reversibility))) # associate with reversibility
  
  regulatory_contingency <- table(regulated = c(1:nrow(ML_inducibility_summary)) %in% grep('act|inh', ML_inducibility_summary$reaction), reversible = aligned_ML_meta$reversibility)
  cat("\nRegulatory reaction vs. rxn reversibility\n")
  print(chisq.test(regulatory_contingency, simulate.p.value = T))
  #regulatory_contingency[2,2]/sum(regulatory_contingency[,2])
  #regulatory_contingency[2,1]/sum(regulatory_contingency[,1])
  
  ML_inducibility_summary$reversibility <- aligned_ML_meta$reversibility
  setkey(ML_inducibility_summary, 'ML_mean')
  ML_inducibility_summary$rxn <- factor(ML_inducibility_summary$rxn, levels = ML_inducibility_summary$rxn)
  
  # Association with trancription factor targets
  
  # convert reaction enzymes to common names
  rxn_enzyme_groups <- read.delim("./flux_cache/rxn_enzyme_groups.tsv")
  library("org.Sc.sgd.db")
  c2o <- toTable(org.Sc.sgdCOMMON2ORF)
  
  # convert from systematic names to common names and then generate a compact summary of genes with consecutive numbers
  ML_reaction_enzymes <- sapply(ML_inducibility_summary$rxn, function(x){
    commonSubset <- sort(c2o$gene_name[chmatch(unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction == x]), c2o$systematic_name)])
    commonSubsetDF <- data.frame(a = regmatches(commonSubset, regexpr('^[A-Z]{3}', commonSubset)), n = regmatches(commonSubset, regexpr('[0-9]+', commonSubset)))
    
    gene_name_compact = NULL
    for(an_a in unique(commonSubsetDF$a)){
      if(length(commonSubsetDF$n[commonSubsetDF$a == an_a]) == 1){
        gene_name_compact <- c(gene_name_compact, paste(an_a, commonSubsetDF$n[commonSubsetDF$a == an_a], sep = ""))
      }else{
        q_seq <- as.numeric(commonSubsetDF$n[commonSubsetDF$a == an_a])
        q_group <- rep(1:length(q_seq))
        for(q_el in 1:(length(q_seq)-1)){
          if(q_seq[q_el] + 1 == q_seq[q_el + 1]){
            q_group[q_el + 1] <- q_group[q_el]
          }
        }
        group_track <- NULL
        for(a_group in unique(q_group)){
          if(length(q_seq[q_group == a_group]) == 1){
            group_track <- c(group_track, q_seq[q_group == a_group])
          }else{
            group_track <- c(group_track, paste(q_seq[q_group == a_group][1], q_seq[q_group == a_group][length(q_seq[q_group == a_group])], sep = "-"))
          }
        }
        gene_name_compact <- c(gene_name_compact, paste(an_a, paste(group_track, collapse = ", ") , sep = ""))
      }
    }
    return(c(expanded = paste(commonSubset, collapse = ", "), collapsed = paste(gene_name_compact, collapse = ", ")))
  })
  ML_reaction_enzymes <- as.data.frame(t(ML_reaction_enzymes))
  
  ML_inducibility_summary$genes <- ML_reaction_enzymes$collapsed
  
  # Determine transcription factors regulating reaction subset
  
  ML_gene_summaries <- data.frame(systematic = unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% ML_inducibility_summary$rxn]),
                                  common = c2o$gene_name[chmatch(unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% ML_inducibility_summary$rxn]), c2o$systematic_name)])
  
  # Import matrix relating transcription factors to their targets (http://www.yeastract.com/generateregulationmatrix.php)
  # Interaction based on "DNA binding and expression evidence", all genes considered #
  # expression or affinity : 10% non-zero
  # affinity : 2% non-zero
  
  # connect every gene involved in a reaction with whether it is a TF target
  TF_indirect <- as.matrix(read.delim("./companionFiles/yeast_TF_regulation.csv", sep = ";", row.names = 1)) # direct and indirect
  TF_indirect <- TF_indirect[,colnames(TF_indirect) %in% ML_gene_summaries$common]
  
  TF_direct <- as.matrix(read.delim("./companionFiles/Yeast_TF_affinity.csv", sep = ";", row.names = 1)) # direct targets
  TF_direct <- TF_direct[,colnames(TF_direct) %in% ML_gene_summaries$common]
  
  if(all(colnames(TF_indirect) == colnames(TF_direct))){
    
    rxn2gene <- matrix(0, nrow = nrow(ML_inducibility_summary), ncol = ncol(TF_direct))
    rownames(rxn2gene) <- ML_inducibility_summary$rxn; colnames(rxn2gene) <- colnames(TF_direct)
    for(i in 1:nrow(ML_reaction_enzymes)){
      rxn2gene[i,colnames(rxn2gene) %in% strsplit(ML_reaction_enzymes$expanded[i], split = ", ")[[1]]] <- 1
    }  
    
  }else{
    stop("TF_direct and TF_indirect gene complements differ -> rxn2gene transformation needs to be modified")
  }
  
  # convert matrix from TF ~ Gene to TF ~ Rxn
  
  TF_indirect_byrxn <- TF_indirect %*% t(rxn2gene); TF_indirect_byrxn[TF_indirect_byrxn != 0] <- 1
  TF_direct_byrxn <- TF_direct %*% t(rxn2gene); TF_direct_byrxn[TF_direct_byrxn != 0] <- 1
  
  # Reduce the number of transcription factors to those with a role in steady-state metabolism #
  # from FIRE, generate a subset of TFs shaping transcription across these conditions (based upon Brauer data)
  FIRE_TFs <- c("Msn2p", "Msn4p", "Gcn4p", "Bas1p", "Cbf1p", "Mbp1p", "Swi4p")
  
  TF_indirect_byrxn <- TF_indirect_byrxn[rownames(TF_indirect_byrxn) %in% FIRE_TFs,]
  TF_direct_byrxn <- TF_direct_byrxn[rownames(TF_direct_byrxn) %in% FIRE_TFs,]
  
  TF_ML_assoc <- data.frame(TF = c(rownames(TF_indirect_byrxn), rownames(TF_direct_byrxn)), Effect = c(rep("indirect", times = nrow(TF_indirect_byrxn)),rep("direct", times = nrow(TF_direct_byrxn))) , p = NA)
  
  for(i in 1:nrow(TF_ML_assoc)){
    if(TF_ML_assoc$Effect[i] == "indirect"){
      refVec <- TF_indirect_byrxn[rownames(TF_indirect_byrxn) == TF_ML_assoc$TF[i],]
    }else{
      refVec <- TF_direct_byrxn[rownames(TF_direct_byrxn) == TF_ML_assoc$TF[i],]
    }
    if(length(unique(refVec)) == 1){next}
    
    TF_ML_assoc$p[i] <- wilcox.test(ML_inducibility_summary$ML_mean[refVec == 1], ML_inducibility_summary$ML_mean[refVec == 0], alternative = "two.sided")$p.value
  }
  
  library(qvalue)
  
  TF_ML_assoc$q <- qvalue(TF_ML_assoc$p)$q
  TF_ML_assoc <- TF_ML_assoc[TF_ML_assoc$q < 0.1,] # look at TFs with an FDR of less than 0.1
  TF_ML_assoc$label <- paste(sub('p$', '', TF_ML_assoc$TF), TF_ML_assoc$Effect, sep = "-")
  
  TFsigSubset <- rbind(TF_indirect_byrxn[chmatch(TF_ML_assoc$TF[TF_ML_assoc$Effect == "indirect"], rownames(TF_indirect_byrxn)),],
                       TF_direct_byrxn[chmatch(TF_ML_assoc$TF[TF_ML_assoc$Effect == "direct"], rownames(TF_direct_byrxn)),])
  rownames(TFsigSubset) <- TF_ML_assoc$label
  
  TF_effect_melt <- data.table(melt(TFsigSubset))
  setnames(TF_effect_melt, colnames(TF_effect_melt), c("TF", "rxn", "target"))
  TF_effect_melt$rxn <- factor(TF_effect_melt$rxn, levels = levels(TF_effect_melt$rxn))
  TF_effect_melt$TF <- factor(TF_effect_melt$TF, levels = rev(sort(unique(as.character(TF_effect_melt$TF)))))
  TF_effect_melt$ypos <- max(ML_inducibility_summary$ML_max) + as.numeric(TF_effect_melt$TF)/25
  TF_effect_melt$fillCol <- ifelse(TF_effect_melt$target == 1, "chocolate1", "aliceblue")  
  
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                         panel.background = element_rect(fill = "gray90"), legend.position = "top", 
                         axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                         axis.text = element_text(color = "black"), axis.text.x = element_text(size = 18, angle = 75, hjust = 1, vjust = 1),
                         panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                         axis.line = element_line(color = "black", size = 1), legend.title=element_blank()
  )
  
  write.table(ML_inducibility_summary, file = "flux_cache/ML_inducibility_summary.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
  
  ggplot() + geom_pointrange(data = ML_inducibility_summary, aes(x = rxn, y = ML_mean, ymin = ML_min, ymax = ML_max), size = 2, color = "blue3") +
    scale_x_discrete("Reactions", breaks = ML_inducibility_summary$rxn, labels = ML_inducibility_summary$genes) + scale_y_continuous("Enyzme Metabolic Leverage", breaks = seq(0,0.8, by = 0.2), expand = c(0,0)) +
    geom_raster(data = TF_effect_melt, aes(x = rxn, y = ypos, fill = fillCol)) + 
    geom_text(data = TF_effect_melt[TF_effect_melt$rxn == ML_inducibility_summary$rxn[1],], aes(x = rxn, y = ypos, label = TF), hjust = 0, size = 7) +
    barplot_theme + scale_fill_identity() + scale_color_identity() + expand_limits(y = 0)
  
  ggsave("Figures/MLstrength.pdf", height = 10, width = 14)
  
  return(ML_inducibility_summary)
  
}

group_regulation <- function(rMech_support, rxnList_form, GS_regulation_support, print_plots = F){
  
  # core functions
  require(dplyr)
  require(igraph)
  # color scheme
  require(colorRamps)
  require(wesanderson)
  
  # two tuning parameters
  corr_penalty <- 0.5 # subtract from correlation s.t. grouping is only favored when cor > 0.5
  interact_penalty <- -2 # penalty for having an interaction within a group of consistent metabolites
  
  # for each reaction with significant regulation (or combinatorial regulation) attempt to group
  # regulators into sets with a similar role
  
  all_reg <- rMech_support %>% filter(type != "unregulated") %>% dplyr::select(reaction, rMech, ncond, AIC_prob, spearman)
  
  # for each reaction consider a metabolite as an activator and inhibitor seperately
  
  regulator_trends <- lapply(all_reg$rMech, function(x){
    
    regs <- rxnList_form[[x]]$rxnFormData %>% filter(Type %in% c("inh", "act")) %>% dplyr::select(SubstrateID, Type)
    regs <- regs %>% left_join(data.frame(SubstrateID = names(rxnList_form[[x]]$metNames), commonName = unname(rxnList_form[[x]]$metNames)), by = "SubstrateID")
    
    rxnList_form[[x]]$rxnMet[,colnames(rxnList_form[[x]]$rxnMet) %in% regs$SubstrateID, drop = F] %>%
      mutate(condition = rownames(.)) %>% gather("SubstrateID", "RA", -condition, convert = T) %>%
      left_join(regs, by = "SubstrateID") %>% mutate(rMech = x)
    
  })
  regulator_trends <- do.call("rbind", regulator_trends)
  
  regulator_trends <- regulator_trends %>% left_join(all_reg, by = "rMech")
  
  reg_ID <- regulator_trends %>% dplyr::select(SubstrateID, Type, commonName) %>% unique() %>%
    mutate(tID_type = paste(SubstrateID, Type, sep = "-"),
           common_type = paste(commonName, ifelse(Type == "act", "+", "-")))
  
  # prepare outputs
  # plot showing the seperation of regulators into sets with a similar role
  reaction_graphs <- list()
  
  for(a_rxn in unique(all_reg$reaction)){
    
    reaction_regulators <- regulator_trends %>% filter(reaction %in% a_rxn) %>% tbl_df()
    
    for(a_cond_subset in unique(sort(reaction_regulators$ncond))){ # allow full and reduced sets of conditions to be evalutated
      
      #warning(paste0(a_rxn, "-", a_cond_subset))
      
      # All metabolite trends for the relevent reaction
      relevant_regulators <- reaction_regulators %>% filter(ncond == a_cond_subset)
      
      # pull down name of reaction for plotting
      
      reaction_name <- ifelse(length(unique(reaction_regulators$ncond)) == 1,
                              paste0(a_rxn, ": ", rxnList_form[[relevant_regulators$rMech[1]]]$reaction),
                              paste0(a_rxn, ": ", rxnList_form[[relevant_regulators$rMech[1]]]$reaction, " (", a_cond_subset, " conditions)")
      )
      
      # look at correlation of regulators
      
      regulator_info <- relevant_regulators %>% dplyr::select(SubstrateID, Type) %>% unique() %>% mutate(tID_type = paste(SubstrateID, Type, sep = "-"))
      
      if(nrow(regulator_info) == 1){
        reaction_graphs$regulation[[paste(a_rxn, a_cond_subset, sep = "-")]] <- paste(relevant_regulators$commonName[1], ifelse(relevant_regulators$Type[1] == "inh", "-", "+"))
        next
      }
      
      # flip inhibitors, so that anticorrelated activators and inhibitors are grouped
      reg_RA <- relevant_regulators %>% left_join(regulator_info, by = c("SubstrateID", "Type")) %>%
        dplyr::select(tID_type, Type, condition, RA) %>%
        unique() %>% mutate(RA = ifelse(Type == "inh", -1*RA, RA)) %>%
        spread("condition", "RA") %>% as.data.frame(.)
      rownames(reg_RA) <- reg_RA$tID_type
      reg_RA <- reg_RA[,!(colnames(reg_RA) %in% c('tID_type', 'Type'))]
      
      reg_corr <- as.matrix(reg_RA) %>% t() %>% cor(.)
      
      # look at interaction of pairwise regulators
      
      regulator_overlap <- relevant_regulators %>% dplyr::select(rMech, SubstrateID, Type) %>% unique() %>%
        group_by(rMech) %>% filter(n() != 1) %>% mutate(overlap = T) %>% spread("rMech", "overlap", fill = F) %>%
        as.data.frame()
      
      rownames(regulator_overlap) <- paste(regulator_overlap$SubstrateID, regulator_overlap$Type, sep = "-")
      regulator_overlap <- regulator_overlap[,-c(1:2), drop = F]
      
      regulator_overlap <- as.matrix(regulator_overlap*1) %*% as.matrix(t(regulator_overlap*1))
      regulator_overlap <- regulator_overlap > 0
      
      # two types of edges between metabolites
      # - joint regulation
      # - equivalent regulation
      
      if(length(regulator_overlap) == 0){
        # no pairwise interaction
        interaction_edges <- NULL
      }else{
        interaction_edges <- regulator_overlap
        interaction_edges[upper.tri(interaction_edges, diag = T)] <- NA
        interaction_edges <- interaction_edges %>% as.data.frame() %>% mutate(tID_type_1 = rownames(.)) %>%
          gather("tID_type_2", "interact", -tID_type_1, convert = T) %>% filter(interact & !is.na(interact)) %>%
          mutate(interact = interact_penalty * interact, type = "pairwise")
      }
      
      corr_edges <- reg_corr
      corr_edges[upper.tri(reg_corr, diag = T)] <- NA
      corr_edges <- corr_edges %>% as.data.frame() %>% mutate(tID_type_1 = rownames(.)) %>%
        gather("tID_type_2", "interact", -tID_type_1, convert = T) %>% filter(!is.na(interact)) %>%
        mutate(interact = interact - corr_penalty, type = "corr")
      
      edge_weights <- rbind(corr_edges, interaction_edges) %>%
        group_by(tID_type_1, tID_type_2) %>% dplyr::summarize(weight = sum(interact))
      
      edge_color_key <- data.frame(weight = seq(-1 - corr_penalty, 1 - corr_penalty, length.out = 1000),
                                   color = green2red(1000))
      
      edge_weights <- edge_weights %>% rowwise() %>% mutate(color = edge_color_key$color[which.min(abs(edge_color_key$weight - weight))])
      
      ### We want to seperate the fully-connected graph into subgraphs with high weight
      # Each metabolites is initially grouped into its own cluster
      # One at a time, each metabolite is tested with every cluster and then
      # assigned to the cluster with the greatest score (+ correlation, - interaction)
      
      all_species <- rownames(reg_corr)
      nclus <- length(all_species)
      clus_assignments <- data.frame(ID = all_species, cluster = 1:nclus)
      
      has_converged = F
      iter <- 0
      old_score <- -Inf
      
      while(!has_converged){
        
        iter <- iter + 1
        
        for(i in 1:nclus){
          current_assignments <- clus_assignments %>% dplyr::slice(-i)
          cluster_score <- rbind(current_assignments, data.frame(ID = all_species[i], cluster = 1:nclus)) %>%
            group_by(cluster) %>% summarize(score = sum(edge_weights$weight[edge_weights$tID_type_1 %in% ID & edge_weights$tID_type_2 %in% ID]) + 0)
          clus_assignments$cluster[i] <- cluster_score$cluster[which.max(cluster_score$score)]
        }
        
        current_score = clus_assignments %>% group_by(cluster) %>% summarize(score = sum(edge_weights$weight[edge_weights$tID_type_1 %in% ID & edge_weights$tID_type_2 %in% ID]) + 0) %>%
          dplyr::select(score) %>% unlist() %>% unname() %>% sum()
        
        if(old_score >= current_score){
          has_converged <- T
        }else{
          old_score <- current_score
        }
      }
      
      if(length(unique(clus_assignments$cluster)) > 1){
        
        if(is.null(interaction_edges)){
          
          interaction_freq <- t(combn(unique(clus_assignments$cluster), 2)) %>% as.data.frame()
          colnames(interaction_freq) <- c("C1", "C2")
          interaction_freq <- interaction_freq %>% mutate(freq = 0, logic_connection = "OR")
          
        }else{
          
          # If multiple clusters exist, how are they related
          
          inter_cluster_interactions <- interaction_edges %>% rowwise() %>% mutate(C1 = clus_assignments$cluster[clus_assignments$ID == tID_type_1],
                                                                                   C2 = clus_assignments$cluster[clus_assignments$ID == tID_type_2]) %>%
            dplyr::select(C1, C2) %>% ungroup() %>% group_by(C1, C2) %>%
            dplyr::summarize(count = n())
          
          inter_cluster_interactions <- rbind(inter_cluster_interactions,
                                              expand.grid(C1 = unique(clus_assignments$cluster), C2 = unique(clus_assignments$cluster), count = 0) %>% apply(c(1,2), as.numeric) %>%
                                                as.data.frame() %>% anti_join(inter_cluster_interactions, by = c("C1", "C2")))
          
          inter_cluster_interactions <- inter_cluster_interactions %>% spread(C1, count) %>% as.data.frame()
          rownames(inter_cluster_interactions) <- inter_cluster_interactions[,1]
          inter_cluster_interactions <- inter_cluster_interactions[,-1]
          
          inter_cluster_interactions[lower.tri(inter_cluster_interactions)] <- inter_cluster_interactions[lower.tri(inter_cluster_interactions)] + t(inter_cluster_interactions)[lower.tri(inter_cluster_interactions)]
          inter_cluster_interactions[upper.tri(inter_cluster_interactions)] <- NA
          
          possible_interactions <- clus_assignments %>% group_by(cluster) %>% summarize(N = n()) %>%
            as.data.frame()
          rownames(possible_interactions) <- possible_interactions$cluster
          possible_interactions <- possible_interactions[,-1, drop = F]
          possible_interactions <- as.matrix(possible_interactions)
          
          possible_pairs <- possible_interactions %*% t(possible_interactions)
          diag(possible_pairs) <- possible_interactions * (possible_interactions-1)
          possible_pairs[upper.tri(possible_pairs)] <- NA
          
          interaction_freq <- inter_cluster_interactions/possible_pairs %>% as.data.frame()
          interaction_freq$C1 <- rownames(interaction_freq)
          interaction_freq <- interaction_freq %>% gather(C2, freq, -C1) %>% filter(!is.na(freq), C1 != C2) %>%
            mutate(logic_connection = ifelse(freq > 0.1, "AND", "OR")) %>%
            mutate(C1 = as.numeric(as.character(C1)), C2 = as.numeric(as.character(C2)))
        }
      }
      ### Optional printing of graph layout ###
      
      pruned_subnetworks <- edge_weights %>% ungroup() %>% group_by(tID_type_1) %>%
        filter(tID_type_2 %in% clus_assignments$ID[clus_assignments$cluster == clus_assignments$cluster[clus_assignments$ID == tID_type_1[1]]])
      
      if(is.null(interaction_edges)){
        
        pruned_subnetworks <- pruned_subnetworks %>% mutate(width = 1.2)
        
      }else{
        
        pruned_subnetworks <- pruned_subnetworks %>% full_join(
          interaction_edges %>% mutate(width = 1) %>% dplyr::select(-interact),
          by = c("tID_type_1", "tID_type_2")) %>%
          mutate(lty = ifelse(is.na(color), 3, 1)) %>%
          mutate(width = ifelse(is.na(width), 1.2, width),
                 weight = ifelse(is.na(weight), 0, weight),
                 color = ifelse(is.na(color), "grey50", color))
        
      }
      
      # Initialize undirected graph from edge list (adding edge attributes too)
      
      reg_connection_graph <- graph.data.frame(pruned_subnetworks, directed = F)
      
      # Setup vertex attributes
      
      reg_significance <- relevant_regulators %>% dplyr::select(rMech, SubstrateID, Type, AIC_prob) %>% unique() %>%
        group_by(rMech) %>% filter(n() == 1) %>% group_by(SubstrateID, Type) %>%
        filter(AIC_prob == max(AIC_prob)) %>%
        mutate(size = log10(AIC_prob) + 5) %>% mutate(size = ifelse(size < 0.5, 0.5, size)) %>%
        ungroup() %>% mutate(tID_type = paste(SubstrateID, Type, sep = "-"))
      
      # add nodes with no edges
      reg_connection_graph <- reg_connection_graph + vertices(reg_significance$tID_type[!(reg_significance$tID_type %in% V(reg_connection_graph)$name)])
      
      # color nodes based on cluster
      color_cluster <- data.frame(cluster = unique(clus_assignments$cluster))
      color_cluster$color <- c(wes_palette("Moonrise3", 5), wes_palette("Zissou", 5))[1:nrow(color_cluster)]
      
      vertex_match <- data.frame(tID_type = V(reg_connection_graph)$name) %>% left_join(reg_ID, by = "tID_type") %>%
        left_join(reg_significance, by = "tID_type") %>%
        left_join(clus_assignments, by = c("tID_type" = "ID")) %>%
        left_join(color_cluster, by = "cluster") %>%
        mutate(size = ifelse(is.na(size), 0.5, size))
      
      V(reg_connection_graph)$label <- vertex_match$common_type
      V(reg_connection_graph)$size <- vertex_match$size*10
      V(reg_connection_graph)$color <- vertex_match$color
      
      # All regulation that improves simpler models
      
      l <- layout.kamada.kawai(reg_connection_graph, weights = E(reg_connection_graph)$weight, initemp = 100, coolexp = 0.995, niter = 5000)
      reg_connection_graph$layout <- l
      reg_connection_graph$main = reaction_name
      
      reaction_graphs$complete[[paste(a_rxn, a_cond_subset, sep = "-")]] <- reg_connection_graph
      if(print_plots){
        plot(reg_connection_graph)
      }
      
      # Filter regulators if they have low support (high AIC) and a substantially lower spearman correlation than the best form (~0.05)
      supported_rMech <- rMech_support %>% filter(reaction == a_rxn, ncond == a_cond_subset) %>% ungroup() %>%
        filter(spearman > max(spearman)-0.05 | AIC_prob > 0.001) %>% arrange(desc(AIC_prob)) %>% dplyr::slice(1:20)
      
      sig_GS_regulation <- GS_regulation_support %>% filter(reaction == a_rxn, is_sig) %>% dplyr::select(SubstrateID = tID, Type = modtype) %>%
        mutate(tID_type = paste(SubstrateID, Type, sep = "-"), rMech = NA)
      
      reduced_regulators <- relevant_regulators %>% dplyr::select(rMech, SubstrateID, Type) %>% unique() %>% filter(rMech %in% supported_rMech$rMech) %>%
        mutate(tID_type = paste(SubstrateID, Type, sep = "-"))
      
      if(nrow(sig_GS_regulation) != 0){
        shedded_reg <- sig_GS_regulation %>% anti_join(reduced_regulators, by = "tID_type")
        shedded_reg <- shedded_reg %>% filter(tID_type %in% regulator_info$tID_type)
        
        if(nrow(shedded_reg) != 0){
          reduced_regulators <- rbind(reduced_regulators, 
                                      sig_GS_regulation %>% anti_join(reduced_regulators, by = "tID_type")
          )
          #print(paste(a_rxn))
        }
      }
      
      # delete irrelevant regulators
      reg_connection_graph <- delete.vertices(reg_connection_graph, V(reg_connection_graph)$name[!(V(reg_connection_graph)$name %in% unique(reduced_regulators$tID_type))])
      
      reduced_edges <- reduced_regulators %>% group_by(rMech) %>% filter(n() != 1)
      
      if(nrow(reduced_edges) != 0){
        reduced_edge_matches <- which(apply(get.edgelist(reg_connection_graph), 1, function(x){
          if(clus_assignments$cluster[clus_assignments$ID == x[1]] == clus_assignments$cluster[clus_assignments$ID == x[2]]){
            # same cluster
            T
          }else{
            # different cluster 
            node_matches <- reduced_edges %>% filter(tID_type %in% x) %>% dplyr::summarize(edge_ind = n())
            if(nrow(node_matches) == 0){
              F 
            }else if(any(node_matches$edge_ind == 2)){
              T}else{F}
          }
        }))
        
        # delete irrelevant pairwise interaction
        reg_connection_graph <- delete.edges(reg_connection_graph, which(!(c(1:length(E(reg_connection_graph))) %in% reduced_edge_matches)))
      }
      
      # Best supported regulation
      
      l <- layout.kamada.kawai(reg_connection_graph, weights = E(reg_connection_graph)$weight, initemp = 100, coolexp = 0.995, niter = 5000)
      reg_connection_graph$layout <- l
     
      reaction_graphs$reduced[[paste(a_rxn, a_cond_subset, sep = "-")]] <- reg_connection_graph
      
      if(print_plots){
        plot(reg_connection_graph)
      }
      
      # From parsimonious regulators, specify combinatorial regulation
      
      reduced_clus_assignments <- clus_assignments %>% filter(ID %in% reduced_regulators$tID_type)
      reduced_clusters <- interaction_freq %>% filter(C1 %in% unique(reduced_clus_assignments$cluster) & C2 %in% unique(reduced_clus_assignments$cluster))
      
      cluster_AIC <- reduced_regulators %>% left_join(supported_rMech %>% dplyr::select(rMech, AIC_prob), by = "rMech") %>%
        left_join(reduced_clus_assignments, by = c("tID_type" = "ID"))
      
      interaction_AIC <- matrix(0, nrow = length(unique(cluster_AIC$cluster)), ncol = length(unique(cluster_AIC$cluster)))
      rownames(interaction_AIC) <- colnames(interaction_AIC) <- sort(unique(cluster_AIC$cluster))
      for(a_rMech in unique(cluster_AIC$rMech)){
        rMech_data <- cluster_AIC %>% filter(rMech == a_rMech) 
        if(nrow(rMech_data) == 1){
          interaction_AIC[rownames(interaction_AIC) == rMech_data$cluster, colnames(interaction_AIC) == rMech_data$cluster] <- 
            interaction_AIC[rownames(interaction_AIC) == rMech_data$cluster, colnames(interaction_AIC) == rMech_data$cluster] + rMech_data$AIC_prob
        }else{
          chunk_store <- interaction_AIC[rownames(interaction_AIC) %in% rMech_data$cluster, colnames(interaction_AIC) %in% rMech_data$cluster]
          chunk_store[lower.tri(chunk_store)] <- chunk_store[lower.tri(chunk_store)] + rMech_data$AIC_prob[1]
          interaction_AIC[rownames(interaction_AIC) %in% rMech_data$cluster, colnames(interaction_AIC) %in% rMech_data$cluster] <- chunk_store
        }
      }
      
      # Summarize the clusters or combinations that need be represented
      all_summary_combinations <- interaction_AIC %>% as.data.frame() %>% mutate(C1 = rownames(.)) %>% gather(C2, "AIC", -C1) %>%
        filter(AIC != 0) %>% mutate(C1 = as.numeric(C1), C2 = as.numeric(as.character(C2))) %>%
        left_join(interaction_freq, by = c("C1", "C2")) %>%
        filter(is.na(logic_connection) | logic_connection == "AND") %>%
        mutate(combinatorial = ifelse(C1 == C2, F, T)) %>%
        dplyr::select(C1, C2, AIC, combinatorial) %>% 
        mutate(summary_order = NA, primary_cluster = NA)
      
      order_counter <- 1
      while(any(is.na(all_summary_combinations$summary_order))){
        
        single_reg <- all_summary_combinations %>% filter(!combinatorial & is.na(summary_order))
        combo_reg <- all_summary_combinations %>% filter(combinatorial & is.na(summary_order))
        
        if(nrow(combo_reg) != 0){
          clusters_remaining <- data.frame(query_cluster = unique(c(combo_reg$C1, combo_reg$C2)), AIC = NA)
          for(a_cluster_n in 1:nrow(clusters_remaining)){
            clus_subset <- combo_reg %>% filter(C1 == clusters_remaining$query_cluster[a_cluster_n] | C2 == clusters_remaining$query_cluster[a_cluster_n])
            clusters_remaining$AIC[a_cluster_n] <- sum(clus_subset$AIC)
          }
          
          if(nrow(single_reg) == 0){
            combo <- T
          }else{
            if(max(clusters_remaining$AIC) > max(single_reg$AIC)){
              combo <- T 
            }else{
              combo <- F
            }}
        }else{
          combo <- F
        }
        
        if(combo){
          # add combinatorial with a primary regulator
          all_summary_combinations$summary_order[all_summary_combinations$combinatorial & (all_summary_combinations$C1 == clusters_remaining$query_cluster[which.max(clusters_remaining$AIC)] |
                                                                                             all_summary_combinations$C2 == clusters_remaining$query_cluster[which.max(clusters_remaining$AIC)])] <- order_counter
          all_summary_combinations$primary_cluster[all_summary_combinations$combinatorial & (all_summary_combinations$C1 == clusters_remaining$query_cluster[which.max(clusters_remaining$AIC)] |
                                                                                               all_summary_combinations$C2 == clusters_remaining$query_cluster[which.max(clusters_remaining$AIC)])] <- clusters_remaining$query_cluster[which.max(clusters_remaining$AIC)]
        }else{
          # add single regulation
          all_summary_combinations$summary_order[is.na(all_summary_combinations$summary_order) & !all_summary_combinations$combinatorial & all_summary_combinations$C1 == single_reg$C1[which.max(single_reg$AIC)]] <- order_counter
        }
        order_counter <- order_counter + 1
      }
      
      # Layout the cluster combinations using a cluster identifier as a standin for metabolites
      
      cluster_printout <- c()
      for(i in 1:max(all_summary_combinations$summary_order)){
        
        if(all_summary_combinations$combinatorial[all_summary_combinations$summary_order == i][1]){
          primary_clus <- all_summary_combinations$primary_cluster[all_summary_combinations$summary_order == i][1]
          secondary_clus <- all_summary_combinations %>% filter(primary_cluster == primary_clus) %>% arrange(desc(AIC)) %>% mutate(secondary_cluster = ifelse(C1 == primary_cluster, C2, C1))
          
          if(nrow(secondary_clus) == 1){
            cluster_printout[i] <- paste0(paste0("C", primary_clus), " AND ",
                                          paste0("C", secondary_clus$secondary_cluster))
          }else{
            cluster_printout[i] <- paste0(paste0("C", primary_clus), " AND (",
                                          paste(paste0("C", secondary_clus$secondary_cluster), collapse = " OR "),
                                          ")")
          }
          
        }else{
          cluster_printout[i] <- paste0("C", all_summary_combinations$C1[all_summary_combinations$summary_order == i])
        }
      }
      
      cluster_members <- cluster_AIC %>% group_by(tID_type, cluster) %>% dplyr::summarize(AIC = sum(AIC_prob)) %>%
        group_by(cluster) %>% arrange(desc(AIC)) %>% mutate(clusterID = paste0("C", cluster)) %>%
        left_join(vertex_match %>% dplyr::select(tID_type, common_type), by = "tID_type")
      
      cluster_members <- cluster_members %>% group_by(clusterID) %>% dplyr::summarize(print = paste(common_type, collapse = " / "), N = n()) %>%
        mutate(print = ifelse(N != 1, sub('$', ')', sub('^', '(', print)), print))
      
      # sub cluster members into cluster placeholders
      for(i in 1:nrow(cluster_members)){
        cluster_printout <- sub(paste0(cluster_members$clusterID[i], '([ )]|$)'), paste0(cluster_members$print[i],'\\1'), cluster_printout)
      }
      cluster_printout <- ifelse(grepl('AND', cluster_printout), sub('$', ']', sub('^', '[', cluster_printout)), cluster_printout)
      
      reaction_graphs$regulation[[paste(a_rxn, a_cond_subset, sep = "-")]] <- paste(cluster_printout, collapse = " OR ")
      
    }
  }
  return(reaction_graphs)
}

allostery_affinity <- function(){
  
  # plot physiological concentrations of alanine
  
  parSubset <- param_set_list[parSetInfo$rx == "r_0816-rm-t_0461-inh-comp_rmCond"]
  load(paste(c("FBGA_files/paramSets/", param_run_info$file[parSubset[[1]]$name$index]), collapse = ""))
  run_rxn <- run_summary[["r_0816-rm-t_0461-inh-comp_rmCond"]]
  
  alanine_conc <- 2^run_rxn$metabolites[,names(run_rxn$rxnSummary$metNames)[run_rxn$rxnSummary$metNames == "L-alanine"], drop = F]
  alanine_conc <- alanine_conc %>% as.data.frame() %>% mutate(condition = rownames(.)) %>% dplyr::select(alanine = t_0461, condition) %>%
    mutate(limitation = substr(condition, 1, 1), GR = substr(condition, 2, length(condition))) %>%
    mutate(limitation = factor(limitation, levels = c("P", "C", "N", "U")), condition = factor(condition, levels = condition))
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 20, face = "bold"), 
                         panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                         axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5),
                         panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                         axis.line = element_line(color = "black", size = 1),
                         axis.ticks.y = element_line(size = 3, color = "black"), axis.ticks.x = element_blank(), axis.ticks.length = unit(0.6, "lines")
  )
  
  ggplot(alanine_conc, aes(x = condition, y = alanine * 1000, group = limitation)) + geom_path(size = 4, alpha = 1, col = "RED",) + geom_point(size = 10, col = "RED") +
    barplot_theme + scale_color_identity(guide = "none") + scale_size_identity() + expand_limits(y = 0) + geom_blank(aes(y = alanine*1000*1.05)) +
    scale_y_continuous("Alanine Concentration (mM)", breaks = c(0, 25, 50, 75, 100, 125), expand = c(0,0)) + scale_x_discrete("Conditions") +
    geom_vline(x = 0, size = 3) + geom_hline(y = 0, size = 3)
  
  ggsave("Figures/alanineConcentration.pdf", height = 10, width = 10)
  
  
  # plot marginal distribution of alanine ki
  
  par_likelihood <- NULL
  par_markov_chain <- NULL
  
  parSubset <- param_set_list[parSetInfo$rx == "r_0816-rm-t_0461-inh-uncomp_rmCond"]
  
  for(i in 1:length(parSubset)){
    
    par_likelihood <- rbind(par_likelihood, data.frame(sample = 1:param_run_info$n_samples[parSubset[[i]]$name$index], likelihood = parSubset[[i]]$lik, index = parSubset[[i]]$name$index))
    par_markov_chain <- rbind(par_markov_chain, parSubset[[i]]$MC)
  }
  
  
  #load(paste(c("FBGA_files/paramSets/", param_run_info$file[param_run_info$index == par_likelihood$index[1]]), collapse = ""))
  #run_rxn <- run_summary[["r_0816-rm-t_0461-inh-uncomp_rmCond"]]
  
  otcase_alanine_ki <- data.frame(ki = 2^par_markov_chain[,'t_0461'])
  otcase_alanine_ki_mle <- otcase_alanine_ki[which.max(par_likelihood$likelihood),]
  otcase_alanine_ki_exp <- 14.8e-3
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 20, face = "bold"), 
                         panel.background = element_rect(fill = "gray80"), legend.position = "right", 
                         axis.ticks.x = element_line(color = "black", size = 1), axis.ticks.y = element_line(color = "black", size = 1),
                         axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5),
                         panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                         axis.line = element_line(color = "black", size = 1)
  )
  
  ggplot(otcase_alanine_ki, aes(x = ki)) + geom_bar(binwidth = 0.05) + 
    geom_vline(xintercept = otcase_alanine_ki_mle, color = "RED", size = 2) + geom_vline(xintercept = otcase_alanine_ki_exp, color = "BLUE", size = 2) +
    scale_x_log10("Affinity", breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2), labels = c("10uM", "100uM", "1mM", "10mM", "100mM", "1M", "10M", "100M"), expand = c(0,0)) + barplot_theme + 
    scale_y_continuous("Counts", expand = c(0,0)) + geom_blank(aes(y=1.1*..count..), binwidth = 0.05, stat="bin") +
    geom_text(data = data.frame(label = c("MLE", "Measured"), x = 8, y = c(87.5, 80), color = c("RED", "BLUE")), aes(x = x, y = y, label = label, color = color), size = 10) +
    scale_color_identity()
  ggsave("Figures/OTCaseAla_dist.pdf", width = 7, height = 7)
  
}


flux_mass_action_regression <- function(){
  
  ##### predict flux using linear regression using either metabolites and enzymes or their log - partition variance explained into that explained by metabolite and by enzymes ####
  
  reaction_pred_log <- data.frame(nmetab = rep(NA, length(grep('rm$', names(rxnList)))), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, varExplainedJointly = NA, TSS = NA)
  reaction_pred_linear <- data.frame(nmetab = rep(NA, length(grep('rm$', names(rxnList)))), nenz = NA, nCond = NA, Fmetab = NA, Fenz = NA, varExplainedTotal = NA, varExplainedMetab = NA, varExplainedEnzy = NA, varExplainedEither = NA, varExplainedJointly = NA, TSS = NA)
  
  for(rxN in grep('rm$', names(rxnList))){
    reaction_pred_linear$nenz[rxN] <- reaction_pred_log$nenz[rxN] <- length(rxnList[[rxN]]$enzymeComplexes[,1])
    reaction_pred_linear$nmetab[rxN] <- reaction_pred_log$nmetab[rxN] <- sum(!is.na(rxnList[[rxN]]$rxnMet))/25
    
    if(all(rxnList[[rxN]]$flux$standardQP >= 0)){
      rxFlux <- log2(rxnList[[rxN]]$flux$standardQP)
    } else if(all(rxnList[[rxN]]$flux$standardQP <= 0)){
      rxFlux <- log2(-1*rxnList[[rxN]]$flux$standardQP)
    } else{
      next
    } #if only forward flux consider its log, if only backwards flux consider the log of -1*flux, if the directionality changes then skip this rxn.
    rxFlux[!is.finite(rxFlux)] <- NA
    
    reaction_pred_log$nCond[rxN] <- reaction_pred_linear$nCond[rxN] <- sum(!is.na(rxFlux))
    
    if(reaction_pred_linear$nenz[rxN] != 0){
      enzyme_df = rxnList[[rxN]]$enzymeComplexes
      
      lenzymes <- as.matrix(data.frame(t(enzyme_df)))
      enzymes <- 2^lenzymes
    }else{lenzymes <- enzymes <- NULL}
    if(reaction_pred_linear$nmetab[rxN] != 0){
      metab_df = rxnList[[rxN]]$rxnMet
      
      lmetabs <- as.matrix(data.frame(metab_df[,colSums(!is.na(metab_df) != 0) != 0])); colnames(lmetabs) <- colnames(rxnList[[rxN]]$rxnMet)[colSums(!is.na(rxnList[[rxN]]$rxnMet) != 0) != 0]
      metabs <- 2^lmetabs
    }else{lmetabs <- metabs <- NULL}
    
    if(reaction_pred_linear$nenz[rxN] == 0 & reaction_pred_linear$nmetab[rxN] == 0){
      next
    }
    
    if(reaction_pred_linear$nenz[rxN] != 0){
      ### only enzymes ###
      ### prediction using log measures ###
      reaction_pred_log$Fenz[rxN] <- anova(lm(rxFlux ~ lenzymes))$F[1]
      reaction_pred_log$varExplainedEnzy[rxN] <- anova(lm(rxFlux ~ lenzymes))$Sum[1]
      reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lenzymes))$Sum)
      
      ### prediction using linear measures ###
      reaction_pred_linear$Fenz[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes))$F[1]
      reaction_pred_linear$varExplainedEnzy[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes))$Sum[1]
      reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes))$Sum)  
    }  
    
    if(reaction_pred_linear$nmetab[rxN] != 0){
      ### only metabolites ###
      ### prediction using log measures ###
      reaction_pred_log$Fmetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$F[1]
      reaction_pred_log$varExplainedMetab[rxN] <- anova(lm(rxFlux ~ lmetabs))$Sum[1]
      reaction_pred_log$TSS[rxN] <- sum(anova(lm(rxFlux ~ lmetabs))$Sum)
      
      ### prediction using linear measures ###
      reaction_pred_linear$Fmetab[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ metabs))$F[1]
      reaction_pred_linear$varExplainedMetab[rxN] <- anova(lm(rxnList[[rxN]]$flux$standardQP ~ metabs))$Sum[1]
      reaction_pred_linear$TSS[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux$standardQP ~ metabs))$Sum)
    }  
    
    if(reaction_pred_linear$nmetab[rxN] != 0 & reaction_pred_linear$nenz[rxN] != 0){
      ### both metabolites and enzymes ###
      reaction_pred_linear$varExplainedTotal[rxN] <- sum(anova(lm(rxnList[[rxN]]$flux$standardQP ~ enzymes + metabs))$Sum[1:2])
      if(reaction_pred_linear$varExplainedTotal[rxN] < max(reaction_pred_linear$varExplainedMetab[rxN], reaction_pred_linear$varExplainedEnzy[rxN])){
        reaction_pred_linear$varExplainedTotal[rxN] <- max(reaction_pred_linear$varExplainedMetab[rxN], reaction_pred_linear$varExplainedEnzy[rxN])  
      }
      
      ### Variance explained in a full model versus reduced ones
      
      if(reaction_pred_linear$varExplainedTotal[rxN] > sum(reaction_pred_linear$varExplainedMetab[rxN] + reaction_pred_linear$varExplainedEnzy[rxN])){
        ### add variance explained in complete model to new class - jointly described
        reaction_pred_linear$varExplainedJointly[rxN] <- reaction_pred_linear$varExplainedTotal[rxN] - sum(reaction_pred_linear$varExplainedMetab[rxN] + reaction_pred_linear$varExplainedEnzy[rxN])
      }else{
        ### some variance is equally accounted for by metabolites or enzymes and should be pulled out and removed from both enzymes and mets
        reaction_pred_linear$varExplainedEither[rxN] <- sum(reaction_pred_linear$varExplainedMetab[rxN] + reaction_pred_linear$varExplainedEnzy[rxN]) - reaction_pred_linear$varExplainedTotal[rxN]
        reaction_pred_linear$varExplainedMetab[rxN] <- reaction_pred_linear$varExplainedMetab[rxN] - reaction_pred_linear$varExplainedEither[rxN]
        reaction_pred_linear$varExplainedEnzy[rxN] <- reaction_pred_linear$varExplainedEnzy[rxN] - reaction_pred_linear$varExplainedEither[rxN]
      }
      
      reaction_pred_log$varExplainedTotal[rxN] <- sum(anova(lm(rxFlux ~ lenzymes + lmetabs))$Sum[1:2])
      if(is.na(reaction_pred_log$varExplainedTotal[rxN])){next}
      
      if(reaction_pred_log$varExplainedTotal[rxN] < max(reaction_pred_log$varExplainedMetab[rxN], reaction_pred_log$varExplainedEnzy[rxN])){
        reaction_pred_log$varExplainedTotal[rxN] <- max(reaction_pred_log$varExplainedMetab[rxN], reaction_pred_log$varExplainedEnzy[rxN])  
      }
      
      if(reaction_pred_log$varExplainedTotal[rxN] > sum(reaction_pred_log$varExplainedMetab[rxN] + reaction_pred_log$varExplainedEnzy[rxN])){
        ### add variance explained in complete model to new class - jointly described
        reaction_pred_log$varExplainedJointly[rxN] <- reaction_pred_log$varExplainedTotal[rxN] - sum(reaction_pred_log$varExplainedMetab[rxN] + reaction_pred_log$varExplainedEnzy[rxN])
      }else{
        ### some variance is equally accounted for by metabolites or enzymes and should be pulled out and removed from both enzymes and mets
        reaction_pred_log$varExplainedEither[rxN] <- sum(reaction_pred_log$varExplainedMetab[rxN] + reaction_pred_log$varExplainedEnzy[rxN]) - reaction_pred_log$varExplainedTotal[rxN]
        reaction_pred_log$varExplainedMetab[rxN] <- reaction_pred_log$varExplainedMetab[rxN] - reaction_pred_log$varExplainedEither[rxN]
        reaction_pred_log$varExplainedEnzy[rxN] <- reaction_pred_log$varExplainedEnzy[rxN] - reaction_pred_log$varExplainedEither[rxN]
      }
      
    }
  }
  
  
  #hist(reaction_pred_linear$Fmetab)
  #hist(reaction_pred_linear$Fenz)
  
  reaction_pred_summary_log <- data.frame(N = reaction_pred_log$nCond, metaboliteVarianceExplained = reaction_pred_log$varExplainedMetab/reaction_pred_log$TSS, enzymeVarianceExplained = reaction_pred_log$varExplainedEnzy/reaction_pred_log$TSS, 
                                          varianceAmbiguouslyExplained = reaction_pred_log$varExplainedEither/reaction_pred_log$TSS, varianceJointlyExplained = reaction_pred_log$varExplainedJointly/reaction_pred_log$TSS)
  
  reaction_pred_summary_linear <- data.frame(N = reaction_pred_linear$nCond, metaboliteVarianceExplained = reaction_pred_linear$varExplainedMetab/reaction_pred_linear$TSS, enzymeVarianceExplained = reaction_pred_linear$varExplainedEnzy/reaction_pred_linear$TSS,
                                             varianceAmbiguouslyExplained = reaction_pred_linear$varExplainedEither/reaction_pred_linear$TSS, varianceJointlyExplained = reaction_pred_linear$varExplainedJointly/reaction_pred_linear$TSS)
  
  rownames(reaction_pred_summary_log) <- rownames(reaction_pred_summary_linear) <- names(rxnList)[grep('rm$', names(rxnList))]
  
  
  residualDF <- reaction_pred_log$nCond - reaction_pred_log$nmetab - reaction_pred_log$nenz
  residualDF <- residualDF[reaction_pred_summary_log$N >= 15]
  
  reaction_pred_summary_log <- reaction_pred_summary_log[reaction_pred_summary_log$N >= 15,-1]
  reaction_pred_summary_linear <- reaction_pred_summary_linear[reaction_pred_summary_linear$N >= 15,-1]
  
  reaction_pred_summary_log <- reaction_pred_summary_log[order(apply(reaction_pred_summary_log, 1, sum, na.rm = TRUE)),]
  reaction_pred_summary_linear <- reaction_pred_summary_linear[order(apply(reaction_pred_summary_linear, 1, sum, na.rm = TRUE)),]
  
  
  
  
  reaction_pred_summary_plotter <- rbind(data.frame(modelType = "logFlux ~ logMetab + logEnzyme", melt(data.frame(index = c(1:length(reaction_pred_summary_log[,1])), reaction_pred_summary_log), id.vars = "index")),
                                         data.frame(modelType = "Flux ~ Metab + Enzyme", melt(data.frame(index = c(1:length(reaction_pred_summary_linear[,1])), reaction_pred_summary_linear), id.vars = "index")))
  
  reaction_pred_summary_plotter <- reaction_pred_summary_plotter[!is.na(reaction_pred_summary_plotter$value),]
  
  barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "top", 
                         panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line = element_blank(), panel.margin = unit(1.5, "lines"), axis.text = element_text(color = "BLACK"),
                         strip.background = element_rect(fill = "chocolate1"), strip.text = element_text(vjust = 0.6)) 
  
  # For 75 reactions, linear regression was used to determine the fraction of variation in flux that could be explained by linear combinations of metabolite and enzyme concentrations.
  # This comparison was also made to relate the log of flux carried to the log abundances of enzymes and metabolites
  
  rxnPredictionPlot <- ggplot(reaction_pred_summary_plotter, aes(x = factor(index), y = value, fill = as.factor(variable), color = "black")) + facet_grid(modelType ~ .)
  rxnPredictionPlot + geom_bar(stat ="identity", width=0.75) + barplot_theme + geom_vline(aes(xintercept = 0), size = 0.5) + geom_hline(aes(yintercept = 0), size = 0.5) + 
    scale_x_discrete(name = "Reactions", expand = c(0,0)) + scale_y_continuous(name = "% Variance Explained", expand = c(0,0), limits = c(0,1), label = percent_format()) +
    scale_fill_brewer("Prediction Method", palette = "Set2") + scale_color_identity()
  
  #scale_fill_manual(values = c("enzymeVarianceExplained" = "sienna1", "metaboliteVarianceExplained" = "steelblue1", varianceAmbiguouslyExplained = "olivedrab3", varianceJointlyExplained = "red")) 
  
  ggsave("varianceExplained.pdf", width = 20, height = 12)
  
}

#res = 200
#control_co = 0
#discretize = T
#dbin = 6

color_simplex <- function(res = 100, control_co = 0.1, discretize = FALSE, dbin = 8){
  
  # generate a color profile that can be used to express a consensus color based upon the relative
  # value of 3 measures summing to one.
  # res - how many bins divide the space into
  # contorl_co - values less than contorl co will display as white
  
  require(reshape2)
  require(data.table)
  require(scales)
  
  color_simplex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), 
                               legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
                               panel.grid.major = element_blank(), axis.line = element_blank(), axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black", size = 2)) 
  
  
  color_summary <- list()
  
  control_lattice <- matrix(NA, ncol = res+1, nrow = res+1)
  rownames(control_lattice) <- seq(0,1,by = 1/res) # parts enzyme - red
  colnames(control_lattice) <- seq(0,1,by = 1/res) # parts allostery - blue
  
  control_colors <- melt(control_lattice)
  colnames(control_colors) <- c("enzyme", "allostery", "color")
  control_colors <- data.table(control_colors)
  
  control_colors <- control_colors[enzyme + allostery <= 1,]
  control_colors[,residual := round(1-allostery-enzyme, 4)]
  
  control_colors[,angle := ifelse(enzyme + allostery != 0, 120*allostery/(allostery+enzyme), 360)]
  #control_colors[,angle := ifelse(enzyme + allostery != 0, 0 + 240*allostery/(allostery+enzyme), 360)]
  #control_colors[,angle := ifelse(enzyme + allostery != 0, 230 + 140*allostery/(allostery+enzyme), 360)]
  control_colors[,angle := ifelse(enzyme + allostery != 0, 0 + 240*allostery/(allostery+enzyme), 360)]
  
  
  if(discretize == F){
    control_colors$residual[control_colors$residual < 0.2] <- 0.2
    control_colors[,color := hcl(h = angle, c = 100, l = 100*residual),]
    
    control_colors$color[control_colors$residual > 1-control_co] <- "#FFFFFF"
    
    color_summary$Figure <- ggplot(control_colors, aes(x = enzyme, y = allostery, fill = color)) + geom_tile() + scale_fill_identity() +
      color_simplex_theme + scale_x_continuous("Enzyme Control", label = percent, expand = c(0,0)) +
      scale_y_continuous("Allosteric Control", label = percent, expand = c(0,0)) +
      geom_hline(yintercept = 0, size = 3) + geom_vline(xintercept = 0, size = 3) +
      geom_abline(intercept = 1, slope = -1, size = 3) + geom_abline(intercept = control_co, slope = -1, size = 3)
    
  }else{
    
    # first determine control magnitude [0,1] -> mag_bin
    # this is the distance projected from the origin, split up into dbin intervals
    # divide angle within a magnitude bin into mag_bin pieces e.g. lowest level is a single bin, outer-most is eight
    
    control_colors[,mag_bin := ceiling((1-residual) * dbin),]
    control_colors <- control_colors[mag_bin != 0,,]
    control_colors$angle_bin <- NA
    
    angle_range <- range(control_colors$angle)
    
    for(a_mag_bin in unique(control_colors$mag_bin)){
      angle_bin <- ceiling((control_colors[mag_bin == a_mag_bin,angle] - min(angle_range)) / (angle_range[2]-angle_range[1]) * a_mag_bin)
      angle_bin[angle_bin == 0] <- 1
      control_colors$angle_bin[control_colors$mag_bin == a_mag_bin] <- angle_bin
    }
    
    control_colors[,consensus_angle := mean(angle), by = c("mag_bin", "angle_bin")]
    control_colors[,consensus_residual := mean(residual), by = c("mag_bin", "angle_bin")]
    
    #control_colors[,color := hcl(h = consensus_angle, c = 100, l = 100*sqrt(consensus_residual)),]
    control_colors$consensus_residual[control_colors$consensus_residual < 0.2] <- 0.2
    control_colors[,color := hcl(h = consensus_angle, c = 100, l = 100*consensus_residual),]
    
    control_colors$color[control_colors$mag_bin == 1] <- "#FFFFFF"
    
    color_summary$Figure <- ggplot(control_colors, aes(x = enzyme, y = allostery, fill = color)) + geom_tile() + scale_fill_identity() +
      color_simplex_theme + scale_x_continuous("Enzyme Control", label = percent, expand = c(0,0)) +
      scale_y_continuous("Allosteric Control", label = percent, expand = c(0,0)) +
      geom_hline(yintercept = 0, size = 3) + geom_vline(xintercept = 0, size = 3) +
      geom_abline(data = data.frame(int = seq(1:dbin)/dbin, slope = rep(-1, dbin)), aes(intercept = seq(1:dbin)/dbin, slope = rep(-1, dbin)), size = 2) 
    
    color_summary$Figure
    
   }
  
  
  color_summary$Table <- control_colors[,list(enzyme, allostery, residual, color)]
  
  return(color_summary)
}





color_ternary <- function(res = 100){
  
  # generate a color profile that can be used to express a consensus color based upon the relative
  # value of 3 measures summing to one.
  # res - how many bins divide the space into
  # contorl_co - values less than contorl co will display as white
  
  
  require(reshape2)
  require(dplyr)
  require(scales)
  require(ggtern)
  
  angle2rad <- function(angle){angle / 180 * pi}
  rad2angle <- function(rad){rad/pi * 180}
  
  
  color_summary <- list()
  
  control_lattice <- matrix(NA, ncol = res+1, nrow = res+1)
  rownames(control_lattice) <- seq(0,1,by = 1/res) # parts enzyme - red
  colnames(control_lattice) <- seq(0,1,by = 1/res) # parts allostery - blue
  
  control_colors <- melt(control_lattice)
  colnames(control_colors) <- c("enzyme", "allostery", "substrates")
  
  control_colors <- control_colors %>% tbl_df() %>% filter(enzyme + allostery <= 1) %>% dplyr::mutate(substrates = round(1 - (allostery + enzyme), 4))
  
  # pair every point with the one above and to the right of it - these will translate into triangles in the ternary plot
  
  x_levels <- sort(unique(control_colors$allostery))
  y_levels <- sort(unique(control_colors$enzyme))
  
  find_trios <- lapply(1:nrow(control_colors), function(i){
    
    # to generate a filled triangle - tile the big triangle with little triangle to upper R /\ and a little triangle to lower L \/
    
    query <- control_colors %>% dplyr::slice(i) %>% as.data.frame()
    query[,'position'] = "L"
    query[,'group'] = i
    
    # define the point to the right and left of the queried point - one or both may not exist
    
    if(query[,'allostery'] == last(x_levels)){right = NULL}else{
      
      right = data.frame(enzyme = query$enzyme, allostery = x_levels[which(query[,'allostery'] == x_levels) + 1])
      right[,'substrates'] = 1 - sum(right)
      right[,'position'] = "R"
      right[,'group'] = i
      if(right[,'substrates'] < 0){right = NULL}
      
    }
    
    if(query[,'allostery'] == first(x_levels)){left = NULL}else{
      
      left = data.frame(enzyme = query$enzyme, allostery = x_levels[which(query[,'allostery'] == x_levels) - 1])
      left[,'substrates'] = 1 - sum(left)
      left[,'position'] = "L"
      left[,'group'] = i
      if(left[,'substrates'] < 0){left = NULL}
      
    }
    
    # define the point above and below the queried - one or both may not exist
    
    if(query[,'enzyme'] == last(y_levels)){top = NULL}else{
      
      top <- data.frame(enzyme = y_levels[which(query[,'enzyme'] == y_levels) + 1], allostery = query$allostery)
      top[,'substrates'] = 1 - sum(top)
      top[,'position'] = "T"
      top[,'group'] = i
      if(top[,'substrates'] < 0){top = NULL}
      
    }
    
    if(query[,'enzyme'] == first(y_levels)){bottom = NULL}else{
      
      bottom <- data.frame(enzyme = y_levels[which(query[,'enzyme'] == y_levels) - 1], allostery = query$allostery)
      bottom[,'substrates'] = 1 - sum(bottom)
      bottom[,'position'] = "B"
      bottom[,'group'] = i
      if(bottom[,'substrates'] < 0){bottom = NULL}
      
    }
    
    # define top triangle
    
    if(!is.null(top) & !is.null(right)){
      
      center = data.frame(enzyme = mean(c(top[,'enzyme'], query[,'enzyme'])), 
                          allostery = mean(c(right[,'allostery'], query[,'allostery'])))
      center[,'substrates'] = 1 - sum(center)
      center[,'position'] = "CT"
      center[,'group'] = i
      
      top_triangle <- rbind(query, top, right, center)
      top_triangle[,'group'] = paste0(top_triangle[,'group'], 'T')
      
    }else{
      top_triangle <- NULL
    }
      
    
    # define bottom triangle
    
    if(!is.null(bottom) & !is.null(left)){
      
      center = data.frame(enzyme = mean(c(bottom[,'enzyme'], query[,'enzyme'])), 
                          allostery = mean(c(left[,'allostery'], query[,'allostery'])))
      center[,'substrates'] = 1 - sum(center)
      center[,'position'] = "CB"
      center[,'group'] = i
      
      bottom_triangle <- rbind(query, bottom, left, center)
      bottom_triangle[,'group'] = paste0(bottom_triangle[,'group'], 'B')
      
    }else{
      bottom_triangle <- NULL
    }
    
    return(rbind(top_triangle, bottom_triangle))
    
  })

  control_colors <- do.call("rbind", find_trios) %>% tbl_df()
  
  if(!all(round(control_colors$enzyme + control_colors$allostery + control_colors$substrates, 3) == 1)){stop("all entries don't sum to unity")}
  
  
  # translate to ternary coordinate system - equilateral triangle
  control_colors <- control_colors %>% dplyr::mutate(x = (1/2)*(2*allostery + enzyme) / (allostery + enzyme + substrates),
                                              y = sqrt(3)/2 * enzyme*(allostery + enzyme + substrates))
  
  # Use centers to define color-scheme
  center_colors <- control_colors %>% group_by(group) %>% filter(position %in% c("CT", "CB"))
  #center_colors <- center_colors %>% rowwise() %>% dplyr::mutate(color_purity = max(enzyme, allostery, substrates) - min(enzyme, allostery, substrates)) %>%
  #  dplyr::select(group, x, y, color_purity)
  center_colors <- center_colors %>% rowwise() %>% dplyr::mutate(color_purity = max(enzyme, allostery, substrates) - min(enzyme, allostery, substrates))
 
  
  # determine angle relative to center for hue
  ternary_center <- c(x = 0.5, y = sin(angle2rad(30)) * 0.5)
  
  center_colors <- center_colors %>% dplyr::mutate(x_centered = x - ternary_center['x'], y_centered = y - ternary_center['y'])
  center_colors <- center_colors %>% dplyr::mutate(angle = rad2angle(atan2(x = x_centered, y = y_centered )))
  center_colors <- center_colors %>% dplyr::mutate(angle = ifelse(angle < 0, angle + 360, angle))
  
  # establish color scheme
  #center_colors <- center_colors %>% dplyr::mutate(color = hcl(h = 240 - (angle - 90), c = 0 + 140*color_purity, l = 70 - 45 * color_purity, fixup = T))
  center_colors <- center_colors %>% dplyr::mutate(color = hcl(h = 240 - (angle - 90), c = 0 + 180*color_purity, l = 90 - 55 * color_purity, fixup = T))
  
  
  # add colors back to defining polygons
  control_colors <- control_colors %>% filter(!(position %in% c("CT", "CB"))) %>% left_join((center_colors %>% dplyr::select(group, color)))
  
  
  color_ternary_theme <- theme(text = element_blank(), panel.background = element_rect(fill = "white"), 
                               legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
                               panel.grid.major = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
  
  vertex_label = data.frame(x = c(-0.15, 0.5, 1.15), y = c(-0.125, sqrt(3)/2 + 0.1, -0.08), label = c('Substrates &\nProducts', 'Enzymes', 'Allostery'), 
                            color = c(control_colors$color[control_colors$substrates == 1], control_colors$color[control_colors$enzyme == 1], control_colors$color[control_colors$allostery == 1]))
  
  ternary_boundary <- data.frame(x = c(0, 0.5, 1), y = c(0, sqrt(3)/2, 0))
  
  text_dist <- 0.08
  tick_dist <- 0.03
  
  ternary_text <- data.frame(x = c(0 - text_dist*cos(angle2rad(30)), 0.5, 1 + text_dist*cos(angle2rad(30)), 
                             cos(angle2rad(60))*0.5 - text_dist*cos(angle2rad(30)), 0.5, 1 - ( cos(angle2rad(60))*0.5 - text_dist*cos(angle2rad(30)))),
                             y = c(0 - text_dist*sin(angle2rad(30)), -text_dist, 0 - text_dist*sin(angle2rad(30)),
                             sin(angle2rad(60))*0.5 + text_dist*sin(angle2rad(30)), sin(angle2rad(60)) + text_dist, sin(angle2rad(60))*0.5 + text_dist*sin(angle2rad(30))),
                             angle = c(0, 0, 0, 60, 0, -60),
                             text = c("100%", "50:50", "100%", "50:50", "100%", "50:50"))
                             
  
  ternary_ticks <- data.frame(x1 = c(0, 0.5, 1, cos(angle2rad(60))/2, 0.5, 1 - cos(angle2rad(60))/2),
                              y1 = c(0, 0, 0, sin(angle2rad(60))/2, sin(angle2rad(60)), sin(angle2rad(60))/2),
                              x2 = c(0 - tick_dist*cos(angle2rad(30)), 0.5, 1 + tick_dist*cos(angle2rad(30)), cos(angle2rad(60))/2 - tick_dist*cos(angle2rad(30)),
                                     0.5, 1 - (cos(angle2rad(60))/2 - tick_dist*sin(angle2rad(60)))),
                              y2 = c(-tick_dist*sin(angle2rad(30)), -tick_dist, -tick_dist*sin(angle2rad(30)), sin(angle2rad(60))/2 + tick_dist*sin(angle2rad(30)),
                                     sin(angle2rad(60)) + tick_dist, sin(angle2rad(60))/2 + tick_dist*sin(angle2rad(30))))
  
  color_summary$Figure_BW <- ggplot() + geom_polygon(data = control_colors, aes(x = x, y = y, group = group), color = "gray80", fill = "gray80") +  color_ternary_theme +
    geom_text(data = vertex_label, aes(x = x, y = y, label = label, color = color), size = 9) +
    geom_polygon(data = ternary_boundary, aes(x = x, y = y), color = "BLACK", fill = NA, size = 1) +
    geom_text(data = ternary_text, aes(x = x, y = y, label = text, angle = angle), color = "BLACK", size = 7) + 
    geom_segment(data = ternary_ticks, aes(x = x1, y = y1, xend = x2, yend = y2), size = 1) +
    scale_fill_identity() + scale_color_identity() + scale_size_identity() 
  
  color_summary$Figure_Color <- ggplot() + geom_polygon(data = control_colors, aes(x = x, y = y, group = group, fill = color, color = color)) +  color_ternary_theme +
    geom_text(data = vertex_label, aes(x = x, y = y, label = label, color = color), size = 9) +
    geom_polygon(data = ternary_boundary, aes(x = x, y = y), color = "BLACK", fill = NA, size = 1) +
    geom_text(data = ternary_text, aes(x = x, y = y, label = text, angle = angle), color = "BLACK", size = 7) + 
    geom_segment(data = ternary_ticks, aes(x = x1, y = y1, xend = x2, yend = y2), size = 1) +
    scale_fill_identity() + scale_color_identity() + scale_size_identity()
  
  color_summary$Table <- control_colors %>% dplyr::select(enzyme, allostery, substrates, color)
  
  return(color_summary)
}


  