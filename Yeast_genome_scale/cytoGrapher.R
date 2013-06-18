#Visualizing flux results using Rcytoscape bridge to cytoscape
#creating a graphical model 

library("RCytoscape")
library(gplots)
library(combinat)

max.add.new <- 2000

#Read in reaction which are going to be laid out - those which carried flux under some condition generated in FBA_run_full_reco.R.  Also load fluxes across each condition.

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")

load("totalStoiAux.Rdata")
stoisub <- Stotal

#reorder metSty and rxnSty to reflect the row and column names of S
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

#read in flux values
flux_vals <- read.delim("carriedFlux.tsv", sep = "\t")



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
        
        met_pos[rownames(met_pos) %in% rownames(defined_mets)[is.na(defined_mets$x)]] <- 1
        rxn_nodes[rx,1:2] <- c(2,0)
        rxn_nodes[rx,3:4] <- c(2,2)
        
      }else{
        defined_mets$x[is.na(defined_mets$x)] <- mean(defined_mets$x[!is.na(defined_mets$x)])
        defined_mets$y[is.na(defined_mets$y)] <- mean(defined_mets$y[!is.na(defined_mets$y)])
        met_pos[chmatch(rownames(defined_mets), rownames(met_pos)),] <- defined_mets
        rxn_nodes[rx,1:2] <- c(mean(defined_mets$x) + 1, mean(defined_mets$y) - 1)
        rxn_nodes[rx,3:4] <- c(mean(defined_mets$x) + 1, mean(defined_mets$y) + 1)
      }
    }
    
    #if only either a subset of products or reactants is defined, but not both, the principal direction vector (going from reactants to products) is determined by the center of the graph.  Otherwise this vector is determined by the position of the defined metabolites, making adjustments to account for whether the number of principal products and reactants is odd or even
    #if there are already reactions attached to a metabolite polarize the new reaction in the opposite direction from the mean angle
    
    if(length(ndefined_react) == 0 | length(ndefined_prod) == 0){
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



#overwrite some nodes with pre-specified coordinates
for(i in 1:length(nodeOver[,1])){
	rxn_nodes[nodeOver$reaction[i],] <- unlist(nodeOver[i, colnames(nodeOver) %in% c("xsub", "ysub", "xprod", "yprod")])
	}
	
	







########## Add the created reaction diagram to cytoscape and specify attributes allowing for flexible visualization ######

library(RCytoscape)

mel_graph <- new("graphNEL", edgemode = "directed")
mel_graph <- initNodeAttribute(graph = mel_graph, attribute.name = "moleculeType", attribute.type = "char", default.value = "undefined")
mel_graph <- initNodeAttribute(graph = mel_graph, attribute.name = "name", attribute.type = "char", default.value = "undefined")
mel_graph <- initNodeAttribute(graph = mel_graph, attribute.name = "compartment", attribute.type = "char", default.value = "undefined")

mel_graph <- initEdgeAttribute(graph = mel_graph, attribute.name = "edgeType", attribute.type = "char", default.value = "produces")
mel_graph <- initEdgeAttribute(graph = mel_graph, attribute.name = "weights", attribute.type = "numeric", default.value = 1)
mel_graph <- initEdgeAttribute(graph = mel_graph, attribute.name = "reaction", attribute.type = "char", default.value = "undefined")
mel_graph <- initEdgeAttribute(graph = mel_graph, attribute.name = "reactionName", attribute.type = "char", default.value = "undefined")


met_pos <- data.frame(met_pos, display_name = rownames(met_pos), stringsAsFactors = FALSE)
rxn_nodes <- data.frame(rxn_nodes, display_name = rownames(rxn_nodes), stringsAsFactors = FALSE)

if(organism == "yeast"){
for(i in 1:length(met_pos[,1])){
	if(rownames(met_pos)[i] %in% metSty$SpeciesID){
		met_pos$display_name[i] <- paste(metSty[metSty$SpeciesID == rownames(met_pos)[i],c(2:3)], collapse = "_")
		met_pos$comp[i] <- metSty$Compartment[metSty$SpeciesID == rownames(met_pos)[i]]
		}}
		
for(i in 1:length(rxn_nodes[,1])){
	if(rownames(rxn_nodes)[i] %in% rxnSty$ReactionID[!is.na(rxnSty$Reaction)]){
		rxn_nodes$display_name[i] <- paste(rxnSty[rxnSty$ReactionID == rownames(rxn_nodes)[i],c(2:3)], collapse = "_")
		}}
	}

#specify metabolite nodes
for(mets in 1:length(met_pos[,1])){
	if(is.nan(met_pos[mets,1])){print(paste(mets, " NaN", collapse = " "))}
	if(!is.na(met_pos[mets,1])){
		mel_graph <- graph::addNode(rownames(met_pos)[mets], mel_graph)
		nodeData(mel_graph, rownames(met_pos)[mets], "moleculeType") <- "primaryMet"
		nodeData(mel_graph, rownames(met_pos)[mets], "name") <- met_pos$display_name[mets]
		nodeData(mel_graph, rownames(met_pos)[mets], "compartment") <- met_pos$comp[mets]
		}
	}
#setNodePosition(mel_graph, rownames(met_pos)[!is.na(met_pos[,1])], met_pos$x[!is.na(met_pos[,1])], met_pos$y[!is.na(met_pos[,1])])
#specify substrate/product linker nodes
for(rxns in 1:length(rxn_nodes[,1])){
	if(is.nan(rxn_nodes[rxns,1])){print(paste(rxns, " NaN", collapse = " "))}
	if(!is.na(rxn_nodes[rxns,1])){
		for(x in c("sub", "prod")){
			mel_graph <- graph::addNode(paste(rownames(rxn_nodes)[rxns], x, sep = "_"), mel_graph)
			nodeData(mel_graph, paste(rownames(rxn_nodes)[rxns], x, sep = "_"), "moleculeType") <- "spNode"
			nodeData(mel_graph, paste(rownames(rxn_nodes)[rxns], x, sep = "_"), "name") <- rxn_nodes$display_name[rxns]
			}
		}
	}
#specify edges
for(rx in 1:length(stoisub[1,])){
#for(rx in 622:length(stoisub[1,])){
		
	rxn_stoi <- stoisub[,rx][stoisub[,rx] != 0]
	cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
	if(length(cofactor_change) != 0){
	cofactor_change <- cofactor_change[names(cofactor_change) %in% cofactor.rxns$cofactor[cofactor.rxns$cofactor %in% names(cofactor_change)][sapply(cofactor.list[cofactor.rxns$cofactor %in% names(cofactor_change)], function(x){!(colnames(stoisub)[rx] %in% x)})]]
		}
	
	principal_change <-  rxn_stoi[!(names(rxn_stoi) %in% names(cofactor_change))]
		 
	if(sum(names(principal_change) %in% split.metab[,1]) != 0){
		
		meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){colnames(stoisub)[rx] %in% x}),]
		
		for(i in 1:length(meta_switch[,1])){
			names(principal_change)[names(principal_change) == meta_switch[i,1]] <- meta_switch$new_name[i]
			}
			}
	
	if(!is.na(rxn_nodes[rx,1])){
	
	if(length(principal_change[principal_change < 0]) > 0){	
		mel_graph <- graph::addEdge(names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), mel_graph)
		edgeData(mel_graph, names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "weights") <- unname(abs(principal_change[principal_change < 0]))
		edgeData(mel_graph, names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "edgeType") <- "produced"
		edgeData(mel_graph, names(principal_change[principal_change < 0]), paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), "reaction") <- rownames(rxn_nodes)[rx]
		}
		
	if(length(principal_change[principal_change > 0]) > 0){
	 	mel_graph <- graph::addEdge(paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), mel_graph)
	 	edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), "weights") <- unname(abs(principal_change[principal_change > 0]))
	 	edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), "edgeType") <- "consumed"
	 	edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), names(principal_change[principal_change > 0]), "reaction") <- rownames(rxn_nodes)[rx]
	 	}
	 	
	 	mel_graph <- addEdge(paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), mel_graph)
		edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "weights") <- 1
		edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "edgeType") <- "reacts"
		edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "reaction") <- rownames(rxn_nodes)[rx]
		edgeData(mel_graph, paste(rownames(rxn_nodes)[rx], "sub", sep = "_"), paste(rownames(rxn_nodes)[rx], "prod", sep = "_"), "reactionName") <- rxn_nodes$display_name[rx]
	}}

#eda.names(mel_graph)
#eda(mel_graph, "weights")

plotter = new.CytoscapeWindow("yeastie4", graph = mel_graph)
#specify node positioning
#options(error = recover)
#options(help.ports=2120)
displayGraph(plotter)

setNodePosition(plotter, rownames(met_pos)[!is.na(met_pos[,1])], met_pos$x[!is.na(met_pos[,1])], -1*met_pos$y[!is.na(met_pos[,1])])
setNodePosition(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,1])], "sub", sep = "_"),rxn_nodes[,1][!is.na(rxn_nodes[,1])], -1*rxn_nodes[,2][!is.na(rxn_nodes[,1])])
setNodePosition(plotter, paste(rownames(rxn_nodes)[!is.na(rxn_nodes[,3])], "prod", sep = "_"),rxn_nodes[,3][!is.na(rxn_nodes[,3])], -1*rxn_nodes[,4][!is.na(rxn_nodes[,3])])

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
#setNodeLabelRule()


#setEdgeLabelRule(obj, edge.attribute.name)
#setEdgeLabelWidthDirect(obj, edge.names, new.value)
#add in 2 more nodes between the sub/prod node and cofactors emerging from them




if(organism == "yeast"){
	col_mat <- matrix(col2rgb(rich.colors(length(unique(met_pos$comp)))), ncol = 3, byrow = TRUE)
	#col_mat <- matrix(col2rgb(rainbow(7)), ncol = 3, byrow = TRUE)
	color_index <- data.frame(compartment = sort(unique(met_pos$comp)), color = mapply(rgb, red = col_mat[,1]/255, green = col_mat[,2]/255, blue = col_mat[,3]/255), stringsAsFactors = FALSE)
	
	for(mets in 1:length(met_pos[,1])){
	if(!is.na(met_pos[mets,1])){
		setNodeColorDirect(plotter, rownames(met_pos)[mets], color_index$color[color_index$compartment == met_pos$comp[mets]])
		}
	}

	flux_mat <- flux_vals
	col_num <- 2
	edge_sf <- 0.1
	library(colorRamps)
	
	flux_att <- color_by_flux(flux_vals, col_num, edge_sf)
	edge_names <- sapply(c(1:length(flux_att[,1])), function(x){paste(flux_att[,colnames(flux_att) == "source"][x], flux_att[,colnames(flux_att) == "dest"][x], sep = "~")})
	all_edges <- cy2.edge.names(mel_graph)
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


	
	}

redraw(plotter)

#setEdge
#grep(matchRE, cy2.edge.names(mel_graph))
#grep(flux_att[1,1], cy2.edge.names(mel_graph))
#grep(flux_att[2,1], cy2.edge.names(mel_graph))


#setEdgeLineStyleDirect
#setEdgeLineWidthDirect
#setEdgeColorDirect



rx <- 73

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
		rxn_stoi <- stoisub[stoisub[,colnames(stoisub) == rxname] != 0,colnames(stoisub) == rxname]
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

			
		



#edgeVals <- sort(unique(eda(mel_graph, "weights")))
#names(eda(mel_graph, "weights")), unname(eda(mel_graph, "weights"))


#setEdgeLineWidthRule(obj, edge.attribute.name, attribute.values, line.widths, default.width)
#getArrowShapes(plotter)
#setEdgeTargetArrowRule
#add flux value - either effects edge width or color

if(organism == "mel"){

load("LHPC_resp.Rdata")
library(colorRamps)
lh_pc[,1] <- lh_pc[,1]*-1


pc_num <- 1
number.col = 1001
colorz <- blue2red(number.col)
col_index <- round((lh_pc + max(abs(range(lh_pc[,pc_num]))))/(2*max(abs(range(lh_pc))))*(number.col-1))


edgeSF <- 0.3
for(x in c(1:length(eda(mel_graph, "weights")))){
	setEdgeLineWidthDirect(plotter, cy2.edge.names(mel_graph)[x], (unname(eda(mel_graph, "weights"))*edgeSF)[x])
	rxn = unname(eda(mel_graph, "reaction"))[x]
	if(!(rxn %in% rownames(col_index))){
		setEdgeColorDirect(plotter, cy2.edge.names(mel_graph)[x], rgb(0,0,0))
		}else{
			setEdgeColorDirect(plotter, cy2.edge.names(mel_graph)[x], colorz[col_index[rownames(col_index) %in% rxn,pc_num]+1])
			}}

	}

#source the options in this script and save the model positioning part
redraw(plotter)








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
				new_mat[n,-((length(metSty[1,])+1):(length(metSty[1,])+2))] <- metSty[metSty$SpeciesID == allNodes[n],]
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
	
for(rx in 1:length(stoisub[1,])){
	
	rxn_stoi <- stoisub[,rx][stoisub[,rx] != 0]
	cofactor_change <- rxn_stoi[names(rxn_stoi) %in% cofactors]
	if(length(cofactor_change) != 0){
	cofactor_change <- cofactor_change[names(cofactor_change) %in% cofactor.rxns$cofactor[cofactor.rxns$cofactor %in% names(cofactor_change)][sapply(cofactor.list[cofactor.rxns$cofactor %in% names(cofactor_change)], function(x){!(colnames(stoisub)[rx] %in% x)})]]
		}
	
	principal_change <-  rxn_stoi[!(names(rxn_stoi) %in% names(cofactor_change))]
		 
	if(sum(names(principal_change) %in% split.metab[,1]) != 0){
		
		meta_switch <- split.metab[split.metab$metabolite %in% names(principal_change),][sapply(rxn.list[split.metab$metabolite %in% names(principal_change)], function(x){colnames(stoisub)[rx] %in% x}),]
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
	
	
	
	
	
	
	
	
	
	