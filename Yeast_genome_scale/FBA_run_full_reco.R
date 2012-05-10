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
		
resource_matches <- lapply(sources, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

boundary_met <- NULL
for(x in 1:length(sources)){
boundary_met <- rbind(boundary_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}

#extrate the IDs of excreted metabolites

excreted <- c("acetate", "ethanol", "succinate(2-)", "(R)-lactate", "L-alanine", "L-glutamate")
#excreted <- c("acetate", "ethanol")



resource_matches <- lapply(excreted, perfect.match, query = corrFile$SpeciesName, corrFile = corrFile)

excreted_met <- NULL
for(x in 1:length(excreted)){
excreted_met <- rbind(excreted_met, resource_matches[[x]][resource_matches[[x]]$Compartment %in% compFile$compID[compFile$compName == "extracellular"],])
}



#tmp <- stoiMat[metIDtoSpec(rownames(stoiMat)) == "(S)-lactate",]
#tmp2 <- tmp[,apply(tmp != 0, 2, sum) != 0]
#rownames(tmp2) <- metIDtoSpec(rownames(tmp2)); colnames(tmp2) <- rxnIDtoEnz(colnames(tmp2))
#grep("lactate", rxnIDtoEnz(colnames(stoiMat)), value = TRUE)





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

#stoiMat[colnames(stoiMat) %in% "r_1812"][stoiMat[colnames(stoiMat) %in% "r_1812"] != 0]
#metIDtoSpec(rownames(stoiMat)[stoiMat[colnames(stoiMat) %in% "r_1812"] != 0])












growth_rate <- data.frame(limit = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[1]}), dr = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[2]}), growth = NA)

flux_vectors <- list()
######################## Set up the linear equations for FBA #######################

for(treatment in 1:length(names(treatment_par))){

#remove reactions which are defective (such 
S_rxns = stoiMat[,!(colnames(stoiMat) %in% c(treatment_par[[treatment]]$auxotrophies, rem.unbalanced))]

#added reactions for boundary fluxes

##### Nutrient Influx #####
##### Unconstrained Chemical Influx #####
##### Excreted Metabolite Efflux #####

influx_rxns <- c(1:length(metabolites))[metabolites %in% c(boundary_met$SpeciesID, freeExchange_met$SpeciesID)]

influxS <- matrix(0, ncol = length(influx_rxns), nrow = length(metabolites))
for(i in 1:length(influx_rxns)){
	influxS[influx_rxns[i],i] <- -1
	}

efflux_rxns <- c(1:length(metabolites))[metabolites %in% excreted_met$SpeciesID]

effluxS <- matrix(0, ncol = length(efflux_rxns), nrow = length(metabolites))
for(i in 1:length(efflux_rxns)){
	effluxS[efflux_rxns[i],i] <- 1
	}

	
	
##### Composition fxn #####

compVec <- rep(0, times = length(metabolites))
for(i in 1:length(comp_met$SpeciesID)){
	compVec[rownames(stoiMat) == comp_met$SpeciesID[i]] <- as.numeric(compositionFile$StoiCoef)[i]
	}

S <- cbind(S_rxns, influxS, effluxS, compVec)		
#colnames(S) <- c(colnames(S_rxns), sapply(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName), function(x){paste(x, "boundary")}), "composition")
colnames(S) <- c(colnames(S_rxns), sapply(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName, excreted_met$SpeciesName), function(x){paste(x, "boundary")}), "composition")


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


effluxG <- matrix(0, ncol = length(S[1,]), nrow = length(excreted_met$SpeciesID))
effluxh <- NULL

for(i in 1:length(excreted_met$SpeciesID)){
effluxG[i, c(1:length(S[1,]))[colnames(S) %in% paste(as.character(excreted_met$SpeciesName)[i], "boundary")]] <- 1
effluxh	<- c(effluxh, 0)
	
	}	
	
#Gtot <- rbind(influxG)#, effluxG)	
#htot <- c(influxh)#), effluxh) 	

Gtot <- rbind(influxG, effluxG)	
htot <- c(influxh, effluxh) 	


############### costFxn - indicates the final rxn in S ######

costFxn = c(rep(0, times = length(S[1,]) -1), 1)

######## use linear programming to maximize biomass #######

linp_solution <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = costFxn, ispos = FALSE)
invert_fluxes <- ifelse(colnames(S) %in% c(sapply(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName), function(x){paste(x, "boundary")}), "composition"), -1, 1)


flux_vectors[[names(treatment_par)[treatment]]] <- linp_solution$X*invert_fluxes

growth_rate$growth[treatment] <- linp_solution$solutionNorm*-1

}


#v <- rep(1, times = length(S[1,]))
#S %*%v
Gtot %*% v >= htot



rxnIDtoEnz("r_1384")

#leucine 0.3
2.515818e-05/-0.2964
#uracil 0.3
6.959310e-06/-0.0599


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

std.reduced_flux_mat <- (reduced_flux_mat - apply(reduced_flux_mat, 1, mean, na.rm = TRUE))/apply(reduced_flux_mat, 1, sd, na.rm = TRUE)

heatmap.2(std.reduced_flux_mat, Colv = FALSE, trace = "n", dendrogram = "row")



######### compare boundary fluxes with their limit #########

limiting_fluxes <- matrix(nrow = length(names(treatment_par)), ncol = length(treatment_par[[1]]$nutrients[,1])); rownames(limiting_fluxes) <- names(treatment_par); colnames(limiting_fluxes) <- treatment_par[[1]]$nutrients[,1]

for(treatment in 1:length(names(treatment_par))){

limiting_fluxes[treatment,]	<-(reduced_flux_mat[sapply(sapply(treatment_par[[treatment]]$nutrients$nutrient, function(x){paste(x, "boundary")}), function(x){c(1:length(reduced_flux_mat[,1]))[rownames(reduced_flux_mat) %in% x]}), treatment]*-1)/treatment_par[[treatment]]$nutrients$conc_per_t
	
	}




########## evaluating rxns ########
"isa acyl-CoA", "isa fatty acid", 

query_rxns <- c("protein production")

eval_mat <- stoiMat[apply(stoiMat[,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns] != 0, 1, sum) != 0,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))


named_stoi <- stoiMat; rownames(named_stoi) <- metIDtoSpec(rownames(named_stoi)); colnames(named_stoi) <- rxnIDtoEnz(colnames(named_stoi))


query_mets <- c("Ser-tRNA(Ser)")

x <- rxn_search(named_stoi, query_mets, is_rxn = FALSE)



matches <- metIDtoSpec(rownames(stoiMat)) %in% query_mets
if(sum(matches) >


stoiMat[metIDtoSpec(rownames(stoiMat)) %in% query_mets,] != 0










#fxns

perfect.match <- function(source, query, corrFile){
all_char <- "[[:graph:][:space:]]"
	
tmp <- corrFile[grep(source, query)[!(grep(source, query) %in% union(grep(paste(all_char, source, sep = ""), query), grep(paste(source, all_char, sep = ""), query)))],]
if(length(tmp[,1]) == 0){tmp <- corrFile[grep(source, query, fixed = TRUE),]}
tmp
	}


rxn_search = function(stoiMat, search_string, is_rxn = TRUE){
	#search by metabolite or reactant and return all reactions and nonzero metabolites.
	if (is_rxn == TRUE){
		colz = grep(search_string, colnames(stoiMat), ignore.case = TRUE)
		} else {
		met = grep(search_string, rownames(stoiMat), ignore.case = TRUE)
		if (length(met) == 1){
			colz = c(1:length(stoiMat[1,]))[stoiMat[met,] != 0]
			} else {
			colz = c(1:length(stoiMat[1,]))[apply(stoiMat[met,], 2, is.not.zero)]
		}}
	
	if(length(colz) == 0){
		print("no hits")
		} else {
			rxns = stoiMat[,colz]
			if(is.vector(rxns)){
				c(colz, rxns[rxns != 0])
				} else {
					output <- rbind(colz, rxns[apply(rxns, 1, is.not.zero),])
					colnames(output) = colnames(stoiMat)[colz]
					output
					}
		}
	}

