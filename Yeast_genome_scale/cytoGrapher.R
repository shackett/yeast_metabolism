#Visualizing flux results using Rcytoscape bridge to cytoscape
#creating a graphical model 

library("RCytoscape")
library(gplots)
library(combinat)

### Functions ###
  				
met_assigner <- function(head_node, rxn_angle, sub_posn, newlen){
						
	#assign metabolites to angles radiating from a node such that minimize the sum of squared angle adjustment
	#returns a filled in version of sub_posn
					
	nmets <- length(sub_posn[,1])
	if(odd(nmets)){
		met_angles <- angle_set_odd[1:nmets]
		}else{
			met_angles <- angle_set_even[1:nmets]
			}
					
	ideal_met_angles <- rxn_angle + met_angles
	actual_met_angles <- apply(sub_posn, 1, function(x){
	if(is.na(x[1])){
		NA
		}else{
	ifelse((x - head_node)[1] >= 0, atan((head_node - x)[2]/(head_node - x)[1]), atan((head_node - x)[2]/(head_node - x)[1]) + pi)}})
	angle_permutations <- permn(ideal_met_angles)
	new_angles <- angle_permutations[which.min(lapply(angle_permutations, function(x){sum((x - actual_met_angles)^2, na.rm = TRUE)}))][[1]]
					
	for(metab in c(1:length(sub_posn[,1]))){
		if(is.na(sub_posn[metab,][1])){
			sub_posn[metab,] <- head_node + newlen*c(cos(new_angles[metab]), sin(new_angles[metab]))
			}
		}
	sub_posn
	}
	
color_by_flux <- function(flux_mat, col_num, edge_sf){
  	
	nonzeroRx <- rownames(flux_mat)[flux_mat[,col_num] != 0]
	rxnFlux <- flux_mat[flux_mat[,col_num] != 0,col_num]
	
	
	medFlux <- summary(abs(rxnFlux))[names(summary(abs(rxnFlux))) == "Median"]
	edgeWeight <- abs(rxnFlux)
	edgeWidth <- sapply(edgeWeight/(10*medFlux), function(x){max(x, 1)})*edge_sf
	edgeStyle <- rep("SOLID", times = length(edgeWidth))
	
	number.col = 1001
	colorz <- blue2red(number.col)
	col_index <- round(edgeWeight/(max(edgeWeight))*(number.col-1))+1
	
	dyn_pop_frame_total <- NULL	
		
	for(rx in 1:length(nonzeroRx)){
		rxname <- nonzeroRx[rx]
		rxn_stoi <- stoiActive[stoiActive[,colnames(stoiActive) == rxname] != 0,colnames(stoiActive) == rxname]
		cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
		if(length(cofactor_change) != 0){
			cofactor_change <- cofactor_change[names(cofactor_change) %in% cofactor.rxns$cofactor[cofactor.rxns$cofactor %in% names(cofactor_change)][sapply(cofactor.list[cofactor.rxns$cofactor %in% names(cofactor_change)], function(x){!(rxname %in% x)})]]
			}
		
		principal_change <-  rxn_stoi[!(names(rxn_stoi) %in% names(cofactor_change))]
			 
		if(sum(names(principal_change) %in% split.metab[,1]) != 0){
			
			meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){rxname %in% x}),]
		
			for(i in 1:length(meta_switch[,1])){
				names(principal_change)[names(principal_change) == meta_switch[i,1]] <- meta_switch$new_name[i]
				}
				}
		
		dyn_pop_frame <- NULL
		if(length(principal_change[principal_change < 0]) > 0){
			dyn_pop_frame <- rbind(dyn_pop_frame, cbind(nonzeroRx[rx], names(principal_change[principal_change < 0]), paste(rxname, "sub", sep = "_"), edgeWidth[rx]*abs(principal_change[principal_change < 0]), colorz[col_index][rx], edgeStyle[rx], ifelse(rxnFlux[rx] > 0, 1, 0)))
			}
		if(length(principal_change[principal_change > 0]) > 0){
			dyn_pop_frame <- rbind(dyn_pop_frame, cbind(nonzeroRx[rx], paste(rxname, "prod", sep = "_"), names(principal_change[principal_change > 0]), edgeWidth[rx]*abs(principal_change[principal_change > 0]), colorz[col_index][rx], edgeStyle[rx], ifelse(rxnFlux[rx] > 0, 1, 0)))
			}
		dyn_pop_frame <- rbind(dyn_pop_frame, cbind(nonzeroRx[rx], paste(rxname, "sub", sep = "_"), paste(rxname, "prod", sep = "_"), edgeWidth[rx]*1, colorz[col_index][rx], edgeStyle[rx], ifelse(rxnFlux[rx] > 0, 1, 0)))
		
	
		dyn_pop_frame_total <- rbind(dyn_pop_frame_total, dyn_pop_frame)
	}
	rownames(dyn_pop_frame_total) <- NULL
	colnames(dyn_pop_frame_total) <- c("rxn", "source", "dest", "width", "color", "line_style", "flip")
	dyn_pop_frame_total
}			

####





max.add.new <- 2000

# read in reaction which are going to be laid out - those which carried flux under some condition generated in FBA_run_full_reco.R.  Also load fluxes across each condition.

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

load("totalStoiAux.Rdata")
load("../ChemicalSpeciesQuant/boundaryFluxes.Rdata")

### replace composition reactions (where for instance all amino acids are bundled together into AA flux composition because this is the appropriate measured pool)
### using species-specific boundary fluxes - i.e. valine incorporated into protein with the concaminant hydrolysis of ATP/GTP

stoisub <- Stotal
stoisub <- stoisub[,grep('composition', colnames(stoisub), invert = TRUE)]


# find the sIDs of each metabolite in the composition function
biomassSpecies <- unique(c(comp_by_cond$compositionFile$AltName, colnames(comp_by_cond$biomassExtensionE)[-1]))
resource_matches <- lapply(biomassSpecies, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)
biomassMet <- NULL
for(x in 1:length(biomassSpecies)){
biomassMet <- rbind(biomassMet, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "cytoplasm"],])
}
biomassMet$stoisubIndex <- chmatch(biomassMet$SpeciesID, rownames(stoisub))


### matrix with composition stoichiometry which will be added to reactions
stoiAppended <- NULL

# maintenance ATP hydrolysis
rxStoi <- data.frame(specie = comp_by_cond$compositionFile$AltName[comp_by_cond$compositionFile$varCategory == "Maintenance ATP hydrolysis"], sID = NA, rowIndex = NA, coef = NA)
rxStoi$sID <- biomassMet$SpeciesID[chmatch(rxStoi$specie, biomassMet$SpeciesName)]
rxStoi$rowIndex <- biomassMet$stoisubIndex[chmatch(rxStoi$specie, biomassMet$SpeciesName)]
rxStoi$coef <- ifelse(rxStoi$specie == "ATP", -1, 1)

atpOut <- matrix(0, ncol = 1, nrow = nrow(stoisub))
atpOut[rxStoi$rowIndex,] <- rxStoi$coef
colnames(atpOut) <- "Maintenance ATP hydrolysis"
stoiAppended <- cbind(stoiAppended, atpOut)

# all other species

for(specie in comp_by_cond$biomassExtensionE$name){
  rxStoi <- data.table(specie = c(specie, colnames(comp_by_cond$biomassExtensionE)[-1]), sID = NA, rowIndex = NA, coef = unlist(c(-1, comp_by_cond$biomassExtensionE[comp_by_cond$biomassExtensionE$name == specie, -1])))
  rxStoi$sID <- biomassMet$SpeciesID[chmatch(rxStoi$specie, biomassMet$SpeciesName)]
  rxStoi$rowIndex <- biomassMet$stoisubIndex[chmatch(rxStoi$specie, biomassMet$SpeciesName)]

  
  specOut <- matrix(0, ncol = 1, nrow = nrow(stoisub))
  specOut[rxStoi[,sum(coef), by = rowIndex]$rowIndex,] <- rxStoi[,sum(coef), by = rowIndex]$V1
  colnames(specOut) <- paste(specie, "to biomass") 
  stoiAppended <- cbind(stoiAppended, specOut)
  }

stoisub <- cbind(stoisub, stoiAppended)




### reorder metSty and rxnSty to reflect the row and column names of S
metSty = read.delim("metSty.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metSty$y <- metSty$y*-1

rxnSty = read.delim("rxnSty.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rxnSty <- rxnSty[sapply(c(1:length(rxnSty[,1])), function(x){c(1:length(rxnSty[,1]))[rxnSty$ReactionID == colnames(stoisub)[x]]}),]
rownames(rxnSty) <- NULL

metab.coord <- metSty[!is.na(metSty$x),]

nodeOver <- rxnSty[!is.na(rxnSty$xsub),]


cofactor.rxns <- read.delim("Layout/cofactor_exceptions.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cofactor.list <- list()
for(i in 1:length(cofactor.rxns[,1])){
  cofactor.list[[i]] <- strsplit(cofactor.rxns$reaction[i], split = ", ")[[1]]
  temp_inv <- strsplit(cofactor.rxns$inverse_rxn[i], split = ", ")[[1]]
  if(length(temp_inv) != 0){
    cofactor.list[[i]] <- union(cofactor.list[[i]], names(stoisub[rownames(stoisub) == cofactor.rxns$cofactor[i],][stoisub[rownames(stoisub) == cofactor.rxns$cofactor[i],] != 0])[!(names(stoisub[rownames(stoisub) == cofactor.rxns$cofactor[i],][stoisub[rownames(stoisub) == cofactor.rxns$cofactor[i],] != 0]) %in% temp_inv)])
  }
}


#fix split.metab
split.metab <- read.delim("Layout/met_split.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rxn.list <- list()
for(i in 1:length(split.metab[,1])){
  rxn.list[[i]] <- strsplit(split.metab$reaction[i], split = ", ")[[1]]
}



#plot(metab.coord$y ~ metab.coord$x, pch = 16)

arm_lengths <- 2
arm_ratio <- 1
spread_angle <- 60/360*2*pi
angle_set_odd <- spread_angle * c(0, rep(1:10, each = 2)*rep(c(1,-1), times = 10))
angle_set_even <- angle_set_odd - spread_angle/2

cofactors <- cofactor.rxns$cofactor

metab_names <- c(rownames(stoisub)[!(rownames(stoisub) %in% split.metab$metabolite)], split.metab$new_name)
#rbind(stoisub[,!(apply(!is.na(rxn_nodes), 1, sum) != 0)], c(1:length(!(apply(!is.na(rxn_nodes), 1, sum) != 0)))[!(apply(!is.na(rxn_nodes), 1, sum) != 0)])

graph_center <- c(0,0)
met_pos <- data.frame(x = rep(NA, times = length(metab_names)), y = rep(NA, times = length(metab_names))); rownames(met_pos) <- metab_names
rxn_added <- rep(FALSE, times = length(stoisub[1,]))
rxn_nodes <- matrix(NA, ncol = 4, nrow = length(stoisub[1,])); colnames(rxn_nodes) <- c("rn_x", "rn_y", "pn_x", "pn_y"); rownames(rxn_nodes) <- colnames(stoisub)
cof_nodes <- NULL

#cof_nodes <- data.frame(cofactor = NA, xpos = X, ypos = X, stoi = X, stringAsFactors = FALSE)

for(met in 1:length(metab.coord[,1])){
  met_pos[rownames(met_pos) == metab.coord[met,1],] <- c(metab.coord$x[met], metab.coord$y[met])
}



#reactions to be evaluated must have at least one already positioned substrate or product

def_splits <- rownames(met_pos)[!is.na(met_pos$x)][sapply(c(1:length(rownames(met_pos)[!is.na(met_pos[,1])])), function(x){rownames(met_pos)[!is.na(met_pos[,1])][x] %in% split.metab$new_name})]
all_semi_defined <- union(rownames(met_pos)[apply(is.na(met_pos), 1, sum) == 0], split.metab$metabolite[split.metab$new_name %in% def_splits])

tmp_mat <- matrix(stoisub[rownames(stoisub) %in% all_semi_defined,], ncol = length(stoisub[1,])); rownames(tmp_mat) <- rownames(stoisub[rownames(stoisub) %in% all_semi_defined,])
tmp_mat[rownames(tmp_mat) %in% split.metab$metabolite[split.metab$new_name %in% def_splits],] <- 0
colnames(tmp_mat) <- colnames(stoisub)

for(k in 1:length(def_splits)){
  tmp_mat[rownames(tmp_mat) == split.metab$metabolite[split.metab$new_name %in% def_splits][k], colnames(tmp_mat) %in% rxn.list[[c(1:length(split.metab[,1]))[split.metab$new_name == def_splits[k]]]]] <- 1
}

def_cof <- rownames(tmp_mat) %in% cofactor.rxns$cofactor
def_cof_index <- c(1:length(def_cof))[def_cof]
for(met in c(1:sum(def_cof))){
  tmp <- 	tmp_mat[def_cof_index[met],] 
  tmp[!(colnames(stoisub) %in% (cofactor.list[[c(1:length(cofactor.rxns[,1]))[cofactor.rxns$cofactor %in% rownames(tmp_mat)[def_cof_index[met]]]]]))] <- 0
  tmp_mat[def_cof_index[met],] <- tmp
}


rxn_to_do <- c(1:ncol(stoisub))[rxn_added == FALSE & apply(tmp_mat!= 0, 2, sum) != 0]

rxn_added[rxn_to_do] <- TRUE


#within each iteration determine the reactions that haven't already been layed out and are connected to at least one already defined specie
#tmp <- c(1:length(stoisub[1,]))[colnames(stoisub) == "r_0510"]

while(length(rxn_to_do) != 0){
  
  #loop through reactions that are going to be defined and 
  rxns.added <- 1
  for(rx in rxn_to_do){
    rxns.added <- rxns.added + 1
    if(rxns.added > max.add.new){break}
    
    #determine which species are products, substrates and which should be treated as cofactors
    
    #specify the reaction products and reactants
    rxn_stoi <- stoisub[,rx][stoisub[,rx] != 0]
    
    #determine whether a metabolite should be a fully connected node or one that is only associated with this single reaction (a cofactor)
    
    cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
    
    if(length(cofactor_change) != 0){
      cofactor_change <- cofactor_change[names(cofactor_change) %in% cofactor.rxns$cofactor[cofactor.rxns$cofactor %in% names(cofactor_change)][sapply(cofactor.list[cofactor.rxns$cofactor %in% names(cofactor_change)], function(x){!(colnames(stoisub)[rx] %in% x)})]]
    }
    principal_change <-  rxn_stoi[!(names(rxn_stoi) %in% names(cofactor_change))]
    
    #determine whether a metabolite should be renamed (to allow for better layout...)
    
    if(sum(names(principal_change) %in% split.metab[,1]) != 0){
      
      meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){colnames(stoisub)[rx] %in% x}),]
      
      for(i in 1:length(meta_switch[,1])){
        names(principal_change)[names(principal_change) == meta_switch[i,1]] <- meta_switch$new_name[i]
      }
    }
    
    odd_react <- odd(length(principal_change[principal_change < 0]))
    odd_prod <- odd(length(principal_change[principal_change > 0]))
    n_react <- length(principal_change[principal_change < 0])
    n_prod <- length(principal_change[principal_change > 0])
    ndefined_react <- c(1:length(principal_change))[names(principal_change) %in% names((apply(!is.na(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),]), 1, sum) != 0)[(apply(!is.na(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),]), 1, sum) != 0) == TRUE])]
    ndefined_prod <- c(1:length(principal_change))[names(principal_change) %in% names((apply(!is.na(met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),]), 1, sum) != 0)[(apply(!is.na(met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),]), 1, sum) != 0) == TRUE])]
    
    if(sum(!(names(principal_change) %in% rownames(met_pos))) != 0){
      print(paste("some species are undefined - probably because they were split:", names(principal_change)[!(names(principal_change) %in% rownames(met_pos))]))
      print(paste("check reaction", colnames(stoisub)[rx]))
      next
    }
    
    if(n_react == 0 | n_prod == 0){
      if(n_react == 0 & n_prod == 0){
        print(paste(colnames(stoisub)[rx], "has zero non-cofactor species"))
        rxn_nodes[rx,1:2] <- c(1,1)
        rxn_nodes[rx,3:4] <- c(-1,-1)
        next
      }
      
      # if a reaction does not have principal reactants and products then lay out reaction nodes
      defined_mets <- met_pos[rownames(met_pos) %in% names(principal_change),] 
      
      if(length(ndefined_react) == 0 & length(ndefined_prod) == 0){
        
        met_pos[rownames(met_pos) %in% rownames(defined_mets)[is.na(defined_mets$x)],] <- 1
        rxn_nodes[rx,1:2] <- c(2,0)
        rxn_nodes[rx,3:4] <- c(2,2)
        
      }else{
        defined_mets$x[is.na(defined_mets$x)] <- mean(defined_mets$x[!is.na(defined_mets$x)])
        defined_mets$y[is.na(defined_mets$y)] <- mean(defined_mets$y[!is.na(defined_mets$y)])
        met_pos[chmatch(rownames(defined_mets), rownames(met_pos)),] <- defined_mets
        rxn_nodes[rx,1:2] <- c(mean(defined_mets$x) + 1, mean(defined_mets$y) - 1)
        rxn_nodes[rx,3:4] <- c(mean(defined_mets$x) + 1, mean(defined_mets$y) + 1)
      }
    next
    }
    
    ### position metabolites and reactions where either only substrates or only reactants have already been laid out (or in the case of neither lay arbitrarely) ###
    
    #if only either a subset of products or reactants is defined, but not both, the principal direction vector (going from reactants to products) is determined by the center of the graph.  Otherwise this vector is determined by the position of the defined metabolites, making adjustments to account for whether the number of principal products and reactants is odd or even
    #if there are already reactions attached to a metabolite polarize the new reaction in the opposite direction from the mean angle
    
    if(length(ndefined_react) == 0 | length(ndefined_prod) == 0){
      
      if(length(ndefined_react) == 0 & length(ndefined_prod) == 0){
        met_pos[rownames(met_pos) %in% rownames(defined_mets)[is.na(defined_mets$x)],] <- 1
        rxn_nodes[rx,1:2] <- c(2,0)
        rxn_nodes[rx,3:4] <- c(2,2)
        next  
      }
      
      #use ifelse to indicate whether reactants or products were defined
      
      reactDef = ifelse(length(ndefined_react) != 0, TRUE, FALSE)
      
      
      if(reactDef){changing <- principal_change < 0}else{changing <- principal_change > 0}
      #center_pos <- apply(met_pos[rownames(met_pos) %in% names(principal_change[changing]),][ifelse(reactDef, ndefined_react, ndefined_prod),], 2, mean)
      
      center_pos <- apply(met_pos[rownames(met_pos) %in% names(principal_change[ifelse(reactDef, ndefined_react, ndefined_prod)]),], 2, mean)
      
      #test_exist <- matrix(rxn_nodes[stoisub[rownames(stoisub) %in% rownames(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),][ndefined_react,]),] != 0,], ncol = 4, byrow = FALSE)
      
      test_exist <- matrix(rxn_nodes[stoisub[rownames(stoisub) %in% rownames(met_pos[rownames(met_pos) %in% names(principal_change[ifelse(reactDef, ndefined_react, ndefined_prod)]),]),] != 0,], ncol = 4, byrow = FALSE)
      
      if(sum(!is.na(test_exist[,1])) != 0){
        graph_center <- apply(matrix(test_exist[c(1:length(test_exist[,1]))[!is.na(test_exist[,1])],], ncol = 4, byrow = FALSE), 2, mean)[1:2]
      }else{
        graph_center <- c(0,0)
      }
      
      angle_det <- center_pos - graph_center; angle_det <- ifelse((center_pos - graph_center)[1] >= 0, atan(angle_det[2]/angle_det[1]), pi + atan(angle_det[2]/angle_det[1]))
      
      midpoint <- c(center_pos[1] + cos(angle_det)*(arm_lengths*3 + 2*arm_lengths*arm_ratio)/2, center_pos[2] + sin(angle_det)*(arm_lengths*3 + 2*arm_lengths*arm_ratio)/2)
      
      #define the center of mass for the side of the reaction with defined species
      if((odd_react == TRUE & reactDef) | (odd_prod == TRUE & !reactDef)){
        anglez <- ifelse((midpoint - center_pos)[1] >= 0, atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1]), pi + atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1])); anglez <- c(cos(anglez), sin(anglez))
        
        if(reactDef){cell_choice <- c(1:2)}else{cell_choice <- c(3:4)}
        rxn_nodes[rx, cell_choice] <- center_pos + arm_lengths*arm_ratio*anglez
        
      }else{
        anglez <- ifelse((midpoint - center_pos)[1] >= 0, atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1]), pi + atan((midpoint - center_pos)[2]/(midpoint - center_pos)[1])) + spread_angle*arm_lengths*arm_ratio/(arm_lengths*1.5 + arm_lengths*arm_ratio); anglez <- c(cos(anglez), sin(anglez))
        
        if(reactDef){cell_choice <- c(1:2)}else{cell_choice <- c(3:4)}
        rxn_nodes[rx, cell_choice] <- center_pos + arm_lengths*arm_ratio*anglez					
      }
      
      if(reactDef){changing <- principal_change < 0}else{changing <- principal_change > 0}
      
      sub_posn <- met_pos[rownames(met_pos) %in% names(changing[changing == TRUE]),]
      met_pos[rownames(met_pos) %in% names(changing[changing == TRUE]),] <- met_assigner(rxn_nodes[rx, cell_choice], anglez, sub_posn, arm_lengths*arm_ratio)
      
      #define the center of mass for the previously undefined side of the reaction		
      anglez <- ifelse(c(rxn_nodes[rx,cell_choice] - graph_center)[1] >= 0, atan((rxn_nodes[rx,cell_choice] - graph_center)[2]/(rxn_nodes[rx,cell_choice] - graph_center)[1]), pi + atan((rxn_nodes[rx,cell_choice] - graph_center)[2]/(rxn_nodes[rx,cell_choice] - graph_center)[1])); anglez <- c(cos(anglez), sin(anglez))
      
      cell_choice <- list()
      if(reactDef){cell_choice[[1]] <- c(3:4); cell_choice[[2]] <- c(1:2)}else{cell_choice[[1]] <- c(1:2); cell_choice[[2]] <- c(3:4)}
      rxn_nodes[rx, cell_choice[[1]]] <- rxn_nodes[rx, cell_choice[[2]]] + (arm_lengths*3 + arm_lengths*arm_ratio)*anglez
      
      anglez <- ifelse((rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[1] >= 0, atan((rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[2]/(rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[1]), pi + atan((rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[2]/(rxn_nodes[rx, cell_choice[[1]]] - rxn_nodes[rx, cell_choice[[2]]])[1]))
      
      if(odd_react == TRUE){splayed_angle = angle_set_odd[c(1:n_react)]}else{splayed_angle = angle_set_even[c(1:n_react)]}
      
      met_posn <- lapply(anglez + splayed_angle, function(x){rxn_nodes[rx, cell_choice[[1]]] + c(cos(x), sin(x))*arm_lengths*arm_ratio})
      if(reactDef){changing <- principal_change > 0}else{changing <- principal_change < 0}
      
      for(i in 1:length(met_posn)){
        met_pos[rownames(met_pos) %in% names(principal_change[changing]),] <- met_posn[[i]]
      }
    }
    
    ### position metabolites and reactions where some substrates and reactants have already been laid out ###
    
    if(length(ndefined_react) != 0 & length(ndefined_prod) != 0){
      
      #if there are defined reactants and products then the direction of the reaction edge linking them will be determined by their position with some offset to account for whether there is an odd or even number of reactants/products
      #get the position of defined reactants and products
      
      sub_pos <- met_pos[rownames(met_pos) %in% names(principal_change)[ndefined_react],]
      prod_pos <- met_pos[rownames(met_pos) %in% names(principal_change)[ndefined_prod],]
      sub_pivot <- apply(sub_pos, 2, mean)
      prod_pivot <- apply(prod_pos, 2, mean)
      
      #scale the pivot length to equal the minimum distance between a 'reactant or product pivot node', initially defined as the average position of defined products or reactants, and each defined product/reactant.
      len_scale <- min(sqrt(apply((matrix(sub_pivot, ncol = 2, nrow = length(ndefined_prod), byrow = TRUE) - prod_pos)^2, 1, sum)))/sqrt(sum((sub_pivot - prod_pivot)^2))
      sub_prod_diff <- prod_pivot - sub_pivot
      sub_prod_angle <- ifelse(sub_prod_diff[1] >= 0, atan(sub_prod_diff[2]/sub_prod_diff[1]), pi + atan(sub_prod_diff[2]/sub_prod_diff[1]))
      
      prod_pivot <- sub_pivot + len_scale*sqrt(sum((sub_pivot - prod_pivot)^2))*c(cos(sub_prod_angle), sin(sub_prod_angle))
      
      len_scale <- min(sqrt(apply((matrix(prod_pivot, ncol = 2, nrow = length(ndefined_react), byrow = TRUE) - sub_pos)^2, 1, sum)))/sqrt(sum((sub_pivot - prod_pivot)^2))
      sub_prod_diff <- prod_pivot - sub_pivot
      sub_prod_angle <- ifelse(sub_prod_diff[1] >= 0, atan(sub_prod_diff[2]/sub_prod_diff[1]), pi + atan(sub_prod_diff[2]/sub_prod_diff[1])) + pi
      
      sub_pivot <- prod_pivot + len_scale*sqrt(sum((sub_pivot - prod_pivot)^2))*c(cos(sub_prod_angle), sin(sub_prod_angle))
      
      #rotate when the number of defined metabolites on one end is even and the total is odd or the number of defined species is odd and the total is even.
      
      sub_prod_diff_l <- sqrt(sum((prod_pivot - sub_pivot)^2))
      
      if(!((odd(n_prod) & odd(length(ndefined_prod))) | (even(n_prod) & even(length(ndefined_prod))))){
        sub_prod_angle <- ifelse(sub_prod_diff[1] >= 0, atan(sub_prod_diff[2]/sub_prod_diff[1]), pi + atan(sub_prod_diff[2]/sub_prod_diff[1]))
        
        angle_adj <- (sub_prod_angle + atan((sin(spread_angle/2)*(sub_prod_diff_l*arm_lengths*arm_ratio/(2*arm_lengths*arm_ratio + 3*arm_lengths)))/(2*arm_lengths*arm_ratio + 3*arm_lengths)))
        prod_pivot <- sub_pivot + sub_prod_diff_l*c(cos(angle_adj), sin(angle_adj))
      }
      
      if(!((odd(n_react) & odd(length(ndefined_react))) | (even(n_react) & even(length(ndefined_react))))){
        sub_prod_angle <- ifelse(sub_prod_diff[1] >= 0, atan(sub_prod_diff[2]/sub_prod_diff[1]), pi + atan(sub_prod_diff[2]/sub_prod_diff[1])) + pi
        
        angle_adj <- (sub_prod_angle + atan((sin(spread_angle/2)*(sub_prod_diff_l*arm_lengths*arm_ratio/(2*arm_lengths*arm_ratio + 3*arm_lengths)))/(2*arm_lengths*arm_ratio + 3*arm_lengths)))
        sub_pivot <- prod_pivot + sub_prod_diff_l*c(cos(angle_adj), sin(angle_adj))
      }	
      
      
      sub_prod_angle <- ifelse((prod_pivot - sub_pivot)[1] >= 0, atan((prod_pivot - sub_pivot)[2]/(prod_pivot - sub_pivot)[1]), pi + atan((prod_pivot - sub_pivot)[2]/(prod_pivot - sub_pivot)[1]))
      
      
      rxn_nodes[rx,1:2] <- sub_pivot + sub_prod_diff_l*arm_lengths*arm_ratio/(2*arm_lengths*arm_ratio + 3*arm_lengths)*c(cos(sub_prod_angle), sin(sub_prod_angle))
      rxn_nodes[rx,3:4] <- prod_pivot + sub_prod_diff_l*arm_lengths*arm_ratio/(2*arm_lengths*arm_ratio + 3*arm_lengths)*c(cos(sub_prod_angle + pi), sin(sub_prod_angle + pi))
      
      #assign metabolites
      
      newlen <- sub_prod_diff_l*arm_lengths*arm_ratio/(2*arm_lengths*arm_ratio + 3*arm_lengths)
      sub_posn <- met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),]
      
      met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),] <- met_assigner(rxn_nodes[rx,1:2], sub_prod_angle, sub_posn, newlen)
      
      sub_posn <- met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),]
      
      met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),] <- met_assigner(rxn_nodes[rx,3:4], sub_prod_angle + pi, sub_posn, newlen)
      
    }
    
  }
  
  #if(!all(!is.nan(rxn_nodes[rx,])) | !all(!is.nan(unlist(met_pos[rownames(met_pos) %in% names(principal_change),])))){print("die"); die}
     
     
  def_splits <- rownames(met_pos)[!is.na(met_pos[,1])][sapply(c(1:length(rownames(met_pos)[!is.na(met_pos[,1])])), function(x){rownames(met_pos)[!is.na(met_pos[,1])][x] %in% split.metab$new_name})]
  all_semi_defined <- union(rownames(met_pos)[apply(is.na(met_pos), 1, sum) == 0], split.metab$metabolite[split.metab$new_name %in% def_splits])
  
  tmp_mat <- matrix(stoisub[rownames(stoisub) %in% all_semi_defined,], ncol = length(stoisub[1,])); rownames(tmp_mat) <- rownames(stoisub[rownames(stoisub) %in% all_semi_defined,])
  tmp_mat[rownames(tmp_mat) %in% split.metab$metabolite[split.metab$new_name %in% def_splits],] <- 0
  colnames(tmp_mat) <- colnames(stoisub)
  
  for(k in 1:length(def_splits)){
    tmp_mat[rownames(tmp_mat) == split.metab$metabolite[split.metab$new_name %in% def_splits][k], colnames(tmp_mat) %in% rxn.list[[c(1:length(split.metab[,1]))[split.metab$new_name == def_splits[k]]]]] <- 1
  }
  
  def_cof <- rownames(tmp_mat) %in% cofactor.rxns$cofactor
  def_cof_index <- c(1:length(def_cof))[def_cof]
  for(met in c(1:sum(def_cof))){
    tmp <- 	tmp_mat[def_cof_index[met],] 
    tmp[!(colnames(stoisub) %in% (cofactor.list[[c(1:length(cofactor.rxns[,1]))[cofactor.rxns$cofactor %in% rownames(tmp_mat)[def_cof_index[met]]]]]))] <- 0
    tmp_mat[def_cof_index[met],] <- tmp
  }
  
  
  rxn_to_do <- c(1:length(stoisub[1,]))[rxn_added == FALSE & apply(tmp_mat!= 0, 2, sum) != 0]
  if(length(rxn_to_do) ==0){
    rxn_to_do <- c(1:length(stoisub[1,]))[rxn_added == FALSE]
  }
  
  rxn_added[rxn_to_do] <- TRUE
  
}

#x <- "r_1436"
#reaction_info(x)
#rxnFile$Compartment[rxnFile$ReactionID == x]
#stoisub[stoisub[,colnames(stoisub) == x] != 0,colnames(stoisub) == x]



#overwrite some nodes with pre-specified coordinates
for(i in 1:length(nodeOver[,1])){
	rxn_nodes[nodeOver$reaction[i],] <- unlist(nodeOver[i, colnames(nodeOver) %in% c("xsub", "ysub", "xprod", "yprod")])
	}
	




#### reduce reactions (rxn_nodes) and metabolites (met_pos) to the set of reactions which carry flux under some condition ####

fluxMat <- read.delim("Flux_analysis/fluxCarriedSimple.tsv")


allLaidRx <- data.frame(index = grep('^r_', rownames(rxn_nodes)), rxName = grep('^r_', rownames(rxn_nodes), value = T))
validRx <- allLaidRx[allLaidRx$rxName %in% rownames(fluxMat),]
validRx <- rbind(validRx, data.frame(index = grep('^r_', rownames(rxn_nodes), invert = T), rxName = grep('^r_', rownames(rxn_nodes), invert = T, value = T)))

stoiActive <- stoisub[,validRx$index]
stoiActive <- stoiActive[rowSums(stoiActive != 0) != 0,]

rxn_nodes <- rxn_nodes[rownames(rxn_nodes) %in% validRx$rxName,] #remove reactions which don't carry flux in any condition


valid_split <- split.metab$new_name[sapply(rxn.list, function(x){length(x[x %in% validRx$rxName])}) != 0] # split metabolites which still exist after reactions are pruned
met_pos <- met_pos[rownames(met_pos) %in% c(valid_split, rownames(stoiActive)[!(rownames(stoiActive) %in% split.metab$metabolite)]),] # all metabolites with some valid associated rxns






### optionally add edges and nodes for cofactors ###


########## Add the created reaction diagram to cytoscape and specify attributes allowing for flexible visualization ######

library(RCytoscape)
# install cytoscapeRPC plugin > can be done through "manage plugins" in Cytoscape

### initialize some graph properties ###

metabGraph <- new("graphNEL", edgemode = "directed")
metabGraph <- initNodeAttribute(graph = metabGraph, attribute.name = "moleculeType", attribute.type = "char", default.value = "undefined")
metabGraph <- initNodeAttribute(graph = metabGraph, attribute.name = "name", attribute.type = "char", default.value = "undefined")
metabGraph <- initNodeAttribute(graph = metabGraph, attribute.name = "compartment", attribute.type = "char", default.value = "undefined")

metabGraph <- initEdgeAttribute(graph = metabGraph, attribute.name = "edgeType", attribute.type = "char", default.value = "produces")
metabGraph <- initEdgeAttribute(graph = metabGraph, attribute.name = "weights", attribute.type = "numeric", default.value = 1)
metabGraph <- initEdgeAttribute(graph = metabGraph, attribute.name = "reaction", attribute.type = "char", default.value = "undefined")
metabGraph <- initEdgeAttribute(graph = metabGraph, attribute.name = "reactionName", attribute.type = "char", default.value = "undefined")

### display names and compartments of reactionse and metabolites ###

met_pos <- data.frame(met_pos, display_name = rownames(met_pos), stringsAsFactors = FALSE)
rxn_nodes <- data.frame(rxn_nodes, display_name = rownames(rxn_nodes), stringsAsFactors = FALSE)

for(i in 1:nrow(met_pos)){
	if(rownames(met_pos)[i] %in% metSty$SpeciesID){
		met_pos$display_name[i] <- paste(metSty[metSty$SpeciesID == rownames(met_pos)[i],c(2:3)], collapse = "_")
		met_pos$comp[i] <- metSty$Compartment[metSty$SpeciesID == rownames(met_pos)[i]]
		}}
		
for(i in 1:nrow(rxn_nodes)){
	if(rownames(rxn_nodes)[i] %in% rxnSty$ReactionID[!is.na(rxnSty$Reaction)]){
		rxn_nodes$display_name[i] <- paste(rxnSty[rxnSty$ReactionID == rownames(rxn_nodes)[i],c(2:3)], collapse = "_")
		}}
	
## Specify Nodes ###

# Metabolite nodes
for(mets in 1:nrow(met_pos)){
	if(is.nan(met_pos[mets,1])){print(paste(mets, " NaN", collapse = " "))}
	if(!is.na(met_pos[mets,1])){
		metabGraph <- graph::addNode(rownames(met_pos)[mets], metabGraph)
		nodeData(metabGraph, rownames(met_pos)[mets], "moleculeType") <- "primaryMet"
		nodeData(metabGraph, rownames(met_pos)[mets], "name") <- met_pos$display_name[mets]
		nodeData(metabGraph, rownames(met_pos)[mets], "compartment") <- met_pos$comp[mets]
		}
	}

# Substrate/product reaction linker nodes 
for(rxns in 1:nrow(rxn_nodes)){
	if(is.nan(rxn_nodes[rxns,1])){print(paste(rxns, " NaN", collapse = " "))}
	if(!is.na(rxn_nodes[rxns,1])){
		for(x in c("sub", "prod")){
			metabGraph <- graph::addNode(paste(rownames(rxn_nodes)[rxns], x, sep = "_"), metabGraph)
			nodeData(metabGraph, paste(rownames(rxn_nodes)[rxns], x, sep = "_"), "moleculeType") <- "spNode"
			nodeData(metabGraph, paste(rownames(rxn_nodes)[rxns], x, sep = "_"), "name") <- rxn_nodes$display_name[rxns]
			}
		}
	}

### Specify edges ###
cofactor_layout <- NULL
cof_node_dist <- 3
preferred_angles <- c(pi/4, -pi/4, pi/3, -pi/3)

for(rx in 1:ncol(stoiActive)){
#for(rx in 622:length(stoisub[1,])){
		
	rxn_stoi <- stoiActive[,rx][stoiActive[,rx] != 0]
	cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
	
  if(length(cofactor_change) != 0){
    cofactor_change <- cofactor_change[names(cofactor_change) %in% cofactor.rxns$cofactor[cofactor.rxns$cofactor %in% names(cofactor_change)][sapply(cofactor.list[cofactor.rxns$cofactor %in% names(cofactor_change)], function(x){!(colnames(stoiActive)[rx] %in% x)})]]
    }
  	
  if(length(cofactor_change) != 0){
	  
	  ### Add cofactor nodes ###
	  parent_pos <- rxn_nodes[rownames(rxn_nodes) == colnames(stoiActive)[rx],]
	  rxAngle <- atan((parent_pos$pn_y - parent_pos$rn_y)/(parent_pos$pn_x - parent_pos$rn_x))
	  
	  rxCoef <- data.frame(reaction = colnames(stoiActive)[rx], ID = names(cofactor_change), name = paste(colnames(stoiActive)[rx], names(cofactor_change), sep = "/"), LongcommonName = metIDtoSpec(names(cofactor_change), T), Stoi = unname(cofactor_change), x = NA, y = NA)
	  rxCoef$commonName <- sapply(rxCoef$LongcommonName, function(x){strsplit(x, split = '-c_[0-9]{1,2}$')[[1]]})
    rxCoef$compartment <- regmatches(rxCoef$LongcommonName, regexpr('c_[0-9]{1,2}$', rxCoef$LongcommonName))
    for(cof in 1:length(cofactor_change)){
	    
	    if(cofactor_change[cof] < 0){
	      coefAngle <- pi + rxAngle - preferred_angles[c(1:length(cofactor_change[cofactor_change < 0]))[names(cofactor_change[cofactor_change < 0]) == rxCoef$ID[cof]]]
	      rxCoef$y[cof] <- parent_pos$rn_y + sin(coefAngle)*cof_node_dist
	      rxCoef$x[cof] <- parent_pos$rn_x + cos(coefAngle)*cof_node_dist
	    }else{
	      coefAngle <- rxAngle + preferred_angles[c(1:length(cofactor_change[cofactor_change > 0]))[names(cofactor_change[cofactor_change > 0]) == rxCoef$ID[cof]]]
	      rxCoef$y[cof] <- parent_pos$pn_y + sin(coefAngle)*cof_node_dist
	      rxCoef$x[cof] <- parent_pos$pn_x + cos(coefAngle)*cof_node_dist
	    }
	  }
    
    #plot(parent_pos$rn_x, parent_pos$rn_y, xlim = c(-30, 15), ylim = c(-15, 30))
    #points(parent_pos$pn_x, parent_pos$pn_y)
    #points(rxCoef$y ~ rxCoef$x, col = "RED")
    
	  metabGraph <- graph::addNode(rxCoef$name, metabGraph)
	  nodeData(metabGraph, rxCoef$name, "moleculeType") <- "Cofactor"
	  nodeData(metabGraph, rxCoef$name, "name") <- rxCoef$commonName
	  nodeData(metabGraph, rxCoef$name, "compartment") <- rxCoef$compartment
	  
	  if(sum(rxCoef$Stoi < 0) != 0){
	    
	    metabGraph <- graph::addEdge(rxCoef$name[rxCoef$Stoi < 0], paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), metabGraph)
	    edgeData(metabGraph, rxCoef$name[rxCoef$Stoi < 0], paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "weights") <- abs(rxCoef$Stoi[rxCoef$Stoi < 0])
	    edgeData(metabGraph, rxCoef$name[rxCoef$Stoi < 0], paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "edgeType") <- "Consumed"
	    edgeData(metabGraph, rxCoef$name[rxCoef$Stoi < 0], paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "reaction") <- rownames(rxn_nodes)[rx]
	    
	  }
	  if(sum(rxCoef$Stoi > 0) != 0){
	    
	    metabGraph <- graph::addEdge(paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), rxCoef$name[rxCoef$Stoi > 0], metabGraph)
	    edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), rxCoef$name[rxCoef$Stoi > 0], "weights") <- abs(rxCoef$Stoi[rxCoef$Stoi > 0])
	    edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), rxCoef$name[rxCoef$Stoi > 0], "edgeType") <- "Produced"
	    edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), rxCoef$name[rxCoef$Stoi > 0], "reaction") <- rownames(rxn_nodes)[rx]
	    
	  }
    
    cofactor_layout <- rbind(cofactor_layout, rxCoef)
    
	  }
	
	principal_change <-  rxn_stoi[!(names(rxn_stoi) %in% names(cofactor_change))]
		 
	if(sum(names(principal_change) %in% split.metab[,1]) != 0){
		
		meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){colnames(stoiActive)[rx] %in% x}),]
		
		for(i in 1:length(meta_switch[,1])){
			names(principal_change)[names(principal_change) == meta_switch[i,1]] <- meta_switch$new_name[i]
			}
			}
	
  if(sum(is.na(met_pos[rownames(met_pos) %in% names(principal_change),1])) != 0){
    print(paste("in reaction", colnames(stoiActive)[rx], ":", paste(rownames(met_pos)[rownames(met_pos) %in% names(principal_change)][is.na(met_pos[rownames(met_pos) %in% names(principal_change),1])], collapse = " & "), "not positioned"))
    next
    }
  
  
	if(!is.na(rxn_nodes[rx,1])){
	
	if(length(principal_change[principal_change < 0]) > 0){	
		metabGraph <- graph::addEdge(names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), metabGraph)
		edgeData(metabGraph, names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "weights") <- unname(abs(principal_change[principal_change < 0]))
		edgeData(metabGraph, names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "edgeType") <- "Consumed"
		edgeData(metabGraph, names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "reaction") <- rownames(rxn_nodes)[rx]
		}
		
	if(length(principal_change[principal_change > 0]) > 0){
	 	metabGraph <- graph::addEdge(paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), metabGraph)
	 	edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), "weights") <- unname(abs(principal_change[principal_change > 0]))
	 	edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), "edgeType") <- "Produced"
	 	edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), "reaction") <- rownames(rxn_nodes)[rx]
	 	}
	 	
	 	metabGraph <- addEdge(paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), metabGraph)
		edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "weights") <- 1
		edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "edgeType") <- "reacts"
		edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "reaction") <- rownames(rxn_nodes)[rx]
		edgeData(metabGraph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "reactionName") <- rxn_nodes$display_name[rx]
	}}

#eda.names(metabGraph)
#eda(metabGraph, "weights")

plotter = new.CytoscapeWindow("yeastie4", graph = metabGraph)
#specify node positioning
#options(error = recover)
#options(help.ports=2120)
displayGraph(plotter)

setNodePosition(plotter, rownames(met_pos)[!is.na(met_pos[,1])], met_pos$x[!is.na(met_pos[,1])], -1*met_pos$y[!is.na(met_pos[,1])])
setNodePosition(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "sub", sep = "_"),rxn_nodes[,1][!is.na(rxn_nodes[,1])], -1*rxn_nodes[,2][!is.na(rxn_nodes[,1])])
setNodePosition(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,3])], "prod", sep = "_"),rxn_nodes[,3][!is.na(rxn_nodes[,3])], -1*rxn_nodes[,4][!is.na(rxn_nodes[,3])])

setNodePosition(plotter, cofactor_layout$name, cofactor_layout$x, -1*cofactor_layout$y)
setNodeSizeDirect(plotter, cofactor_layout$name, 2)
setNodeShapeDirect(plotter, cofactor_layout$name, "diamond")


#hide all of the reaction nodes
setNodeFillOpacityDirect(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "sub", sep = "_"), 0)
setNodeFillOpacityDirect(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "prod", sep = "_"), 0)
setNodeBorderOpacityDirect(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "sub", sep = "_"), 0)
setNodeBorderOpacityDirect(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "prod", sep = "_"), 0)
setNodeSizeDirect(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "sub", sep = "_"), 0.01)
setNodeSizeDirect(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "prod", sep = "_"), 0.01)

setDefaultNodeSize(plotter, 2)
setDefaultNodeFontSize(plotter, 0.5)
setNodeLabelRule(plotter, "name")

setEdgeLabelRule(plotter, "reaction")
setDefaultEdgeFontSize(plotter, 0.1)
setDefaultBackgroundColor(plotter, '#000000') # black background
setDefaultEdgeColor(plotter, '#FF0033') # red edges
setDefaultNodeLabelColor(plotter, '#FFFFFF') # white node labels
showGraphicsDetails(plotter, TRUE) #make it so that cytoscape doesn't suppress labels when zoomed out

redraw(plotter)

#setEdgeLabelRule(obj, edge.attribute.name)
#setEdgeLabelWidthDirect(obj, edge.names, new.value)

### color nodes by compartment ###

col_mat <- matrix(col2rgb(rich.colors(length(unique(met_pos$comp)))), ncol = 3, byrow = TRUE)
color_index <- data.frame(compartment = sort(unique(met_pos$comp)), color = mapply(rgb, red = col_mat[,1]/255, green = col_mat[,2]/255, blue = col_mat[,3]/255), stringsAsFactors = FALSE)

for(mets in 1:length(met_pos[,1])){
  if(!is.na(met_pos[mets,1])){
    setNodeColorDirect(plotter, rownames(met_pos)[mets], color_index$color[color_index$compartment == met_pos$comp[mets]])
  }
}

flux_mat <- fluxMat
col_num <- 6
edge_sf <- 0.1
library(colorRamps)

flux_att <- color_by_flux(flux_mat, col_num, edge_sf)
edge_names <- sapply(c(1:length(flux_att[,1])), function(x){paste(flux_att[,colnames(flux_att) == "source"][x], flux_att[,colnames(flux_att) == "dest"][x], sep = "~")})
all_edges <- cy2.edge.names(metabGraph)
edge_names2 <- unname(unlist(sapply(edge_names, function(x){all_edges[names(all_edges) == x]})))
edge_names2 <- unname(unlist(sapply(edge_names, function(x){
  if(x %in% names(all_edges)){
    all_edges[names(all_edges) == x]
  }else{
    NA
  }
})))

valid_edge <- edge_names2 %in% all_edges

noflux_edge <- !(names(all_edges) %in% edge_names)
noflux_edgename <- unname(all_edges[noflux_edge])



frac_match <- rep(NA, times = length(unique(flux_att[,colnames(flux_att) == "rxn"])))
for(i in 1:length(unique(flux_att[,colnames(flux_att) == "rxn"]))){
  submatch <- flux_att[flux_att[,colnames(flux_att) == "rxn"] %in% unique(flux_att[,colnames(flux_att) == "rxn"])[i],]
  if(is.vector(submatch) == TRUE){
    sub_names  <- paste(submatch[names(submatch) == "source"], submatch[names(submatch) == "dest"], sep = "~")
  }else{
    sub_names <- sapply(c(1:length(submatch[,1])), function(x){paste(submatch[,colnames(submatch) == "source"][x], submatch[,colnames(submatch) == "dest"][x], sep = "~")})
  }
  frac_match[i] <- sum(sub_names %in% names(all_edges))/length(sub_names %in% names(all_edges))
}


flux_att_red <- flux_att[valid_edge,]
edge_names2_red <- edge_names2[valid_edge]

for(edge in c(1:length(edge_names2_red))){	
  setEdgeLineWidthDirect(plotter, edge_names2_red[edge], as.numeric(flux_att_red[,colnames(flux_att_red) == "width"][edge]))
  setEdgeColorDirect(plotter, edge_names2_red[edge], flux_att_red[,colnames(flux_att_red) == "color"][edge])
}



for(edge in c(1:length(noflux_edgename))){	
  setEdgeLineWidthDirect(plotter, noflux_edgename[edge], edge_sf)
  setEdgeColorDirect(plotter, noflux_edgename[edge], rgb(0,0,0))
}

redraw(plotter)

#setEdge
#grep(matchRE, cy2.edge.names(metabGraph))
#grep(flux_att[1,1], cy2.edge.names(metabGraph))
#grep(flux_att[2,1], cy2.edge.names(metabGraph))


#setEdgeLineStyleDirect
#setEdgeLineWidthDirect
#setEdgeColorDirect



rx <- 73



			
		









nodePos <- matrix(unlist(getNodePosition(plotter, rownames(met_pos)[!is.na(met_pos$x)])), ncol = 2, byrow = TRUE)
nodeName <- rownames(met_pos)[!is.na(met_pos$x)]
allNodes <- sort(union(nodeName,  metSty[,1]))

new_mat <- as.data.frame(matrix(NA, nrow = length(allNodes), ncol = length(metSty[1,])+2), stringsAsFactors = FALSE)
 rownames(new_mat) <- sort(allNodes); colnames(new_mat) <- c(colnames(metSty), (paste("x", (length(metSty[1,]) - 3)/2, sep = "")), (paste("y", (length(metSty[1,]) - 3)/2, sep = "")))

for(n in 1:length(allNodes)){
	if(allNodes[n] %in% nodeName){
		new_mat[n,(length(metSty[1,])+1):(length(metSty[1,])+2)] <- nodePos[nodeName == allNodes[n],]
		if(allNodes[n] %in% metSty$SpeciesID){
			new_mat[n,-((length(metSty[1,])+1):(length(metSty[1,])+2))] <- metSty[metSty$SpeciesID == allNodes[n],]
			}else{
				new_mat[n,1] <- allNodes[n]
				new_mat[n,c(2:3)] <- metSty[metSty$SpeciesID == split.metab[split.metab$new_name == allNodes[n],]$metabolite,][c(2,3)]
				}
			}else{
				new_mat[n,] <- c(metSty[metSty$SpeciesID == allNodes[n],], metSty[metSty$SpeciesID == allNodes[n],][(length(metSty[metSty$SpeciesID == allNodes[n],]) - 1):length(metSty[metSty$SpeciesID == allNodes[n],])])
				}}
    
joint.table <- new_mat	


#all_mets <- metSty[,1]
#bind_frame <- data.frame(temp1 = rep(NA, times = length(all_mets)), temp2 = rep(NA, times = length(all_mets)))
#colnames(bind_frame) <- c((paste("x", (length(metSty[1,]) - 3)/2, sep = "")), (paste("y", (length(metSty[1,]) - 3)/2, sep = "")))
#for(i in 1:length(nodeName)){
#	bind_frame[all_mets == nodeName[i],] <- nodePos[i,]
#		}
#joint.table <- data.frame(metSty, bind_frame)
write.table(joint.table, file = "joint.table.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


#get a node position
getNodePosition(obj, node.names)


getNodePosition(plotter, rownames(met_pos))

 



	
	
plot(met_pos[,2] ~ met_pos[,1], col = "RED", pch = 16)
segments(rxn_nodes[,1], rxn_nodes[,2], rxn_nodes[,3], rxn_nodes[,4])	
	
for(rx in 1:length(stoiActive[1,])){
	
	rxn_stoi <- stoiActive[,rx][stoiActive[,rx] != 0]
	cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
	if(length(cofactor_change) != 0){
	cofactor_change <- cofactor_change[names(cofactor_change) %in% cofactor.rxns$cofactor[cofactor.rxns$cofactor %in% names(cofactor_change)][sapply(cofactor.list[cofactor.rxns$cofactor %in% names(cofactor_change)], function(x){!(colnames(stoiActive)[rx] %in% x)})]]
		}
	
	principal_change <-  rxn_stoi[!(names(rxn_stoi) %in% names(cofactor_change))]
		 
	if(sum(names(principal_change) %in% split.metab[,1]) != 0){
		
		meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){colnames(stoiActive)[rx] %in% x}),]
		#meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){rx %in% x}),]
		for(i in 1:length(meta_switch[,1])){
			names(principal_change)[names(principal_change) == meta_switch[i,1]] <- meta_switch$new_name[i]
			}
			#if(length(meta_switch) != 0){print(rx)}
		}
	
	if(length(principal_change[principal_change < 0]) > 0){	
	segments(met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),][,1], met_pos[rownames(met_pos) %in% names(principal_change[principal_change < 0]),][,2], rxn_nodes[rx,1], rxn_nodes[rx,2], col = "GREEN")
		}
	if(length(principal_change[principal_change > 0]) > 0){
	 segments(met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),][,1], met_pos[rownames(met_pos) %in% names(principal_change[principal_change > 0]),][,2], rxn_nodes[rx,3], rxn_nodes[rx,4], col = "GREEN")
		}
	}
	
	
					
					
					
					
					
	
	
	
	
	
	
	
	