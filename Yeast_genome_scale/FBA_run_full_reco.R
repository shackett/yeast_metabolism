library(lpSolve)
library(limSolve)
library(gplots)

setwd("/Users/seanhackett/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale")
source("FBA_lib.R")

options(stringsAsFactors = FALSE)

#inputFilebase = "yeast_GS"
inputFilebase = "yeast"


rxnFile = read.delim(paste("rxn_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
rxnparFile = read.delim(paste("species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
corrFile = read.delim(paste("spec_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
compFile <- read.delim(paste("comp_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)

metComp <- read.delim("METeleComp.tsv", stringsAsFactors = FALSE)
compositionFile <- read.csv2("../Yeast_comp.csv", sep = ",", stringsAsFactors = FALSE)
nutrientFile <- read.delim("Boer_nutrients.txt")[1:6,1:6]
rownames(nutrientFile) <- nutrientFile[,1]; nutrientFile <- nutrientFile[,-1]
reversibleRx <- read.delim("../EcoliYeastMatch/revRxns.tsv", sep = "\t", header = TRUE)

load("checkRev.R")
reversibleRx[,2][reversibleRx[,1] %in% misclass] <- 0
reversibleRx[,2][reversibleRx[,1] %in% c("r_0241", "r_0862", "r_0938", "r_0246", "r_0247", "r_0248", "r_0478", "r_0277")] <- 1

#ATPsynthase: c("r_0246", "r_0247", "r_0248")
#r_0478
#r_0277 - may be reversible

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]
metabolites <- unique(rxnStoi$Metabolite)


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

#rxnIDtoEnz(unique(rxnFile[grep("orotidine", rxnFile$Reaction),]$ReactionID))
rxchoose <- "r_0653"
x <- stoiMat[stoiMat[,colnames(stoiMat) == rxchoose] != 0,colnames(stoiMat) == rxchoose]
names(x) <- metIDtoSpec(names(x))
x

metchoose <- "s_1160"
x <- rxn_search(Stotal, metchoose, is_rxn = FALSE)
rownames(x) <- metIDtoSpec(rownames(x)); #colnames(x) <- rxnIDtoEnz(colnames(x))

#metIDtoSpec(c("s_1063", "s_1059"))
#rxnIDtoEnz(c("r_0885", "r_0707"))

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


skip_me = TRUE

if(skip_me == FALSE){

######### Confirm mass balance of reactions ############

met_chebi <- unlist(sapply(metabolites, metToCHEBI))
met_chebi_comp <- metComp[metComp$ID %in% met_chebi,]
met_chebi_comp <- met_chebi_comp[,c(TRUE, TRUE, apply(met_chebi_comp[,-c(1,2)], 2, sum) != 0)]

ele_comp <- lapply(met_chebi, function(x){
	if(length(met_chebi_comp[met_chebi_comp$ID %in% x,]) != 0){
	met_chebi_comp[met_chebi_comp$ID %in% x,]}
	})

#matrix of elemental abundance of species corresponding to rows of the stoichiometric matrix
ele_comp_mat <- matrix(NA, nrow = length(stoiMat[,1]), ncol = length(ele_comp[[1]])-2); rownames(ele_comp_mat) <- rownames(stoiMat); colnames(ele_comp_mat) <- names(ele_comp[[1]])[-c(1:2)]
for(el in 1:length(stoiMat[1,])){
	if(length(unlist(ele_comp[[el]][-c(1:2)])) != 0){
		ele_comp_mat[el,] <- unlist(ele_comp[[el]][-c(1:2)])
		}
	}

missed_compounds <- met_chebi[!is.na(met_chebi)][names(met_chebi[!is.na(met_chebi)]) %in% rownames(ele_comp_mat)[apply(is.na(ele_comp_mat), 1, sum) != 0]]

#create a data.frame of species with a chebi ID, but without a chemical formula match in my composition file
#missed_df <- data.frame(metIDtoSpec(names(missed_compounds)), names(missed_compounds), missed_compounds)
#unique_missed_df <- as.data.frame(matrix(NA, nrow = length(unique(missed_compounds)), ncol = 3))
#for(i in 1:length(unique(missed_compounds))){
#	unique_missed_df[i,] <- missed_df[missed_df[,3] == sort(unique(missed_compounds))[i],][1,]
#	}

#unique_missed_df[,3] %in% chem_form$COMPOUND_ID


#out <- as.data.frame(matrix(0, ncol = length(colnames(ele_comp_mat)), nrow = length(unique_missed_df[,1])))
#colnames(out) <- c(colnames(ele_comp_mat))

#write.table(cbind(unique_missed_df, out), row.names = FALSE, col.names = TRUE, file = "woot", sep = "\t")
add.chebi.comp <- read.delim("chebi_srh_curated.tsv", stringsAsFactors = FALSE)
add.chebi.comp <- add.chebi.comp[is.na(add.chebi.comp$generic),]
add.chebi.comp$internal_ID <- sapply(add.chebi.comp$internal_ID, function(x){corrFile$SpeciesType[corrFile$SpeciesID == x]})

for(i in 1:length(ele_comp_mat[,1])){
	
	if(corrFile$SpeciesType[corrFile$SpeciesID == rownames(ele_comp_mat)[i]] %in% add.chebi.comp$internal_ID){
		ele_comp_mat[i,] <- unlist(add.chebi.comp[add.chebi.comp$internal_ID %in% corrFile$SpeciesType[corrFile$SpeciesID == rownames(ele_comp_mat)[i]], -c(1:3, length(add.chebi.comp[1,]))])
		}}
		
		
		




mass_balanced <- matrix(NA, nrow = length(stoiMat[1,]), ncol = length(ele_comp_mat[1,]) + 1); rownames(mass_balanced) <- colnames(stoiMat); colnames(mass_balanced) <- c("missingIDs", colnames(ele_comp_mat))
mass_balanced <- as.data.frame(mass_balanced, stringsAsFactors = FALSE)

for(colnum in c(1:length(stoiMat[1,]))){
	
	if(length(stoiMat[,colnum][stoiMat[,colnum] != 0]) > 1){
	if(sum(is.na(ele_comp_mat[stoiMat[,colnum] != 0,][,1])) != 0){
		#some species are missing IDs
		#ignore those species and determine for which species the reaction is balanced
		mass_balanced$missingIDs[colnum] <- TRUE
	
		defined_spec <- !is.na(ele_comp_mat[stoiMat[,colnum] != 0,][,1])
		if(sum(defined_spec) > 1){
		mass_balanced[colnum, -1] <- t(ele_comp_mat[stoiMat[,colnum] != 0,][defined_spec,])%*% (stoiMat[,colnum][stoiMat[,colnum] != 0][defined_spec]) == 0
		}
		}else{
			
			mass_balanced$missingIDs[colnum] <- FALSE
			
			if(is.vector(ele_comp_mat[stoiMat[,colnum] != 0,]) == FALSE){
			mass_balanced[colnum, -1] <- t(ele_comp_mat[stoiMat[,colnum] != 0,]) %*% (stoiMat[,colnum][stoiMat[,colnum] != 0]) == 0
			}
			}}}
	
#reactions to use: those that have transformed metabolites
good_rxns <- colnames(stoiMat)[!is.na(mass_balanced[,2])]	
	
add_rxns <- mass_balanced[colnames(stoiMat) %in% good_rxns,][mass_balanced[colnames(stoiMat) %in% good_rxns,]$missingIDs == TRUE,]	
	
	
	
	
#reactions carrying flux with elemental composition information	
non_na_MB <- mass_balanced[rownames(mass_balanced) %in% rownames(reduced_flux_mat),][!is.na(mass_balanced[rownames(mass_balanced) %in% rownames(reduced_flux_mat),]$P),]

#reactions carrying flux that are not mass balanced for an element

non_na_MB_stoi <- stoiMat[apply(stoiMat[,colnames(stoiMat) %in% rownames(non_na_MB[non_na_MB$P == FALSE,])] != 0, 1, sum) != 0,colnames(stoiMat) %in% rownames(non_na_MB[non_na_MB$P == FALSE,])]

#for met
rxnparFile[rxnparFile[,1] == corrFile$SpeciesType[corrFile$SpeciesID == "s_0504"],]


mb.ids <- sapply(rownames(non_na_MB_stoi)[is.na(ele_comp_mat[rownames(ele_comp_mat) %in% rownames(non_na_MB_stoi),][,1])], function(x){corrFile$SpeciesType[corrFile$SpeciesID == x]})
#for rxn 
rxnFile[rxnFile$ReactionID %in% "r_1279",]
rxnFile[rxnFile$ReactionID %in% "r_1281",]
rxnFile[rxnFile$ReactionID %in% "r_1448",]
#destroy glycine-cleavage complex (lipoylprotein)

stoiMat <- stoiMat[,!is.na(mass_balanced[,2])]

rxnIDtoEnz(colnames(stoiMat)[is.na(mass_balanced[,2])])




rxnparFile[rxnparFile[,1] %in% unique(mb.ids),]

metIDtoSpec(rownames(non_na_MB_stoi)[is.na(ele_comp_mat[rownames(ele_comp_mat) %in% rownames(non_na_MB_stoi),][,1])])
metIDtoSpec("s_0504")

	
met_dict <- metIDtoSpec(rownames(non_na_MB_stoi))
rxn_dict <- rxnIDtoEnz(colnames(non_na_MB_stoi))

rownames(non_na_MB_stoi) <- sapply(c(1:length(non_na_MB_stoi[,1])), function(x){met_dict[x][[1]]})
colnames(non_na_MB_stoi) <- sapply(c(1:length(non_na_MB_stoi[1,])), function(x){rxn_dict[x][[1]]})

rxnparFile[rxnparFile[,1] %in% unique(mb.ids),]


}

#metIDtoSpec("s_0334")

#rxnparFile[,3][rxnparFile[,1] %in% corrFile$SpeciesType[corrFile$SpeciesID %in% "r_0267"]]

#taken from chebi
#chem_form <- read.delim("../Yeast_reconstruction/Sequences/chemical_data.tsv", stringsAsFactors = FALSE)
#chem_form <- chem_form[chem_form$SOURCE == "ChEBI",]


#chem_form[chem_form$ID %in% met_chebi,]
#rxnparFile[grep("7814", rxnparFile[,3]),]




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


#boundary_put <- stoiMat[,is.unbalanced][apply(abs(stoiMat[,is.unbalanced]), 1, sum) != 0,]

#rownames(boundary_put) <- metIDtoSpec(rownames(boundary_put))
#colnames(boundary_put) <- rxnIDtoEnz(colnames(boundary_put))

#all_rxns <- rxnIDtoEnz(colnames(stoiMat))



##### generic reactions #####

named_stoi <- stoiMat
met_dict <- metIDtoSpec(rownames(named_stoi))
met_dict <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
rxn_dict <- rxnIDtoEnz(colnames(named_stoi))
rxn_dict <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})

rownames(named_stoi) <- met_dict
colnames(named_stoi) <- rxn_dict

labelz <- c("isa", "protein production", "biomass production", "growth", "lipid production", "IPC synthase")
aggregate_rxns <- NULL

for(l in 1:length(labelz)){
	aggregate_rxns <- union(aggregate_rxns, rxn_search(named_stoi, labelz[l], is_rxn = TRUE, index = TRUE))
	}

#grep(labelz[l], colnames(named_stoi), fixed = TRUE)
#aggregate_rxns <- c(colnames(rxn_search(named_stoi, "isa", is_rxn = TRUE)), colnames(rxn_search(named_stoi, "protein production", is_rxn = TRUE)))

rem.aggregate <- colnames(stoiMat)[aggregate_rxns]



#look for rxns that produce CO2 and make them irreversible

carb_match <- rxn_search(named_stoi, "carbon dioxide", is_rxn = FALSE, index = TRUE)

co_two_producing_rx <- apply(stoiMat[met_dict == "carbon dioxide",carb_match] < 0, 2, sum) == 0
co_two_producing_rx <- names(co_two_producing_rx)[co_two_producing_rx]

reversibleRx[,2][reversibleRx[,1] %in% co_two_producing_rx] <- 1






growth_rate <- data.frame(limit = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[1]}), dr = sapply(names(treatment_par), function(x){unlist(strsplit(x, " "))[2]}), growth = NA)

flux_vectors <- list()
######################## Set up the linear equations for FBA #######################

for(treatment in 1:length(names(treatment_par))){

#reversibleRx[reversibleRx[,1] %in% c("r_0246", "r_0328", "r_0329", "r_0629", "r_0630", "r_0631", "r_0632", "r_0691", "r_0692", "r_0693", "r_0694", "r_0938", "r_0939", "r_0940"),2] <- 1

#remove reactions which are defective 
S_rxns = stoiMat[,!(colnames(stoiMat) %in% c(treatment_par[[treatment]]$auxotrophies, rem.unbalanced, rem.aggregate))]


#added reactions for boundary fluxes

##### Nutrient Influx #####
##### Unconstrained Chemical Influx #####
##### Excreted Metabolite Efflux #####

influx_rxns <- sapply(c(boundary_met$SpeciesID, freeExchange_met$SpeciesID), function(id){c(1:length(metabolites))[metabolites == id]})

influxS <- matrix(0, ncol = length(influx_rxns), nrow = length(metabolites))
for(i in 1:length(influx_rxns)){
	influxS[influx_rxns[i],i] <- 1
	}

efflux_rxns <- sapply(excreted_met$SpeciesID, function(id){c(1:length(metabolites))[metabolites == id]})

effluxS <- matrix(0, ncol = length(efflux_rxns), nrow = length(metabolites))
for(i in 1:length(efflux_rxns)){
	effluxS[efflux_rxns[i],i] <- -1
	}

	
	
##### Composition fxn #####

compVec <- rep(0, times = length(metabolites))
for(i in 1:length(comp_met$SpeciesID)){
	compVec[rownames(stoiMat) == comp_met$SpeciesID[i]] <- as.numeric(compositionFile$StoiCoef)[i]
	}

S <- cbind(S_rxns, influxS, effluxS, compVec)		
colnames(S) <- c(colnames(S_rxns), sapply(c(boundary_met$SpeciesName, freeExchange_met$SpeciesName, excreted_met$SpeciesName), function(x){paste(x, "boundary")}), "composition")
#if(is.na(treatment_par[[treatment]]$auxotrophies)){save(S, file = "totalStoi.Rdata")}


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
	
#thermodynamic flux-constraints

thermoG <- NULL
thermoh <- NULL
rxstack <- NULL
cnstr_rxns <- c(1:length(S[1,]))[colnames(S) %in% reversibleRx[,1]]
for(i in cnstr_rxns){
	revC <- reversibleRx[,2][reversibleRx[,1] == colnames(S)[i]]
	if(revC != 0){
		pointVec <- rep(0, times = length(S[1,]))
		pointVec[i] <- 1*revC
		thermoG <- rbind(thermoG, pointVec)
		thermoh <- c(thermoh, 0)
		rxstack <- c(rxstack, colnames(S)[colnames(S) %in% reversibleRx[,1]][i])
		}
	}

	

Gtot <- rbind(influxG, effluxG, thermoG)
htot <- c(influxh, effluxh, thermoh) 	

#misclass <- "r_0885"
#save(misclass, file = "checkRev.R")

#Gtot <- rbind(influxG, effluxG, thermoG[167,])	
#htot <- c(influxh, effluxh, thermoh[167]) 

#Gtot <- rbind(influxG, effluxG, thermoG[c(1:166,168:204),])	
#htot <- c(influxh, effluxh, thermoh[c(1:166,168:204)]) 


############### costFxn - indicates the final rxn in S ######
#lp loss fxn
costFxn = c(rep(0, times = length(S[1,]) -1), -1)

#qp loss fxn
#costFxn = c(rep(0, times = length(S[1,]) -1), 1)





######## use linear programming to maximize biomass #######

linp_solution <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = costFxn, ispos = FALSE)

#linp_solution$X[names(linp_solution$X) %in% c("r_0486", "r_0488", "r_0689", "r_1003", "r_0246")]


flux_vectors[[names(treatment_par)[treatment]]] <- linp_solution$X

growth_rate$growth[treatment] <- linp_solution$solutionNorm*-1

}



#determine whether reactions can carry flux under the current boundary conditons

rxnGR <- data.frame(rxn = colnames(S), fluxF = NA, boundedF = NA, fluxR = NA, boundedR = NA)
x <- t(sapply(rxnGR[1:5,1], optimize_rx_flux))
rxnGR[,2:5] <- t(sapply(rxnGR[,1], optimize_rx_flux))

optimize_rx_flux <- function(rx){
	
	cost1 <- ifelse(colnames(S) %in% rx, -1, 0)
	linp_solution1 <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = cost1, ispos = FALSE)
	
	cost2 <- ifelse(colnames(S) %in% rx, 1, 0)
	linp_solution2 <- linp(E = S, F = Fzero, G = Gtot, H = htot, Cost = cost2, ispos = FALSE)
	
	c(linp_solution1$solutionNorm*-1, linp_solution1$IsError, linp_solution2$solutionNorm*-1, linp_solution2$IsError)
	
	}


no_flux <- rxnGR$rxn[!((rxnGR$fluxF != 0 | rxnGR$boundedF == 1) | (rxnGR$fluxR != 0 | rxnGR$boundedR == 1))]
#rxns that can't carry flux 
rxnIDtoEnz(no_flux)


rxnGR$rxn[rxnGR$growth == 0 & rxnGR$bounded == 0]




#colorz <- rep(c(1:5), each = 6)
#plot(growth_rate$growth, col = colorz, , xlab = "condition", ylab = "growth rate", pch = 16)
#legend("topleft", unique(growth_rate$limit), text.col = c(1:5))


#uncooperative rxns

#rxnSet <- reversibleRx[,1][reversibleRx[,1] %in% colnames(S)[apply((thermoG[c(1:length(thermoG[,1]))[!(c(1:length(thermoG[,1])) %in% c(1:28,30:32,34:35,37:44,47:58,60:64,66:71,73:74,76:77,79:82))],]) != 0, 2, sum) != 0]]

#rxnum <- 7
#rxnstoi <- S[S[,colnames(S) == rxnSet[rxnum]] != 0 ,colnames(S) == rxnSet[rxnum]]
#names(rxnstoi) <- metIDtoSpec(names(rxnstoi))
#rxnstoi

save(flux_vectors, growth_rate, treatment_par, file = "Flux_analysis/SweaveFluxFilez.Rdata")
save(rxnFile, rxnparFile, corrFile, compFile, metComp, compositionFile, nutrientFile, reversibleRx, file = "Flux_analysis/SweaveNetFilez.Rdata")





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

std.reduced_flux_mat <- (renamed_reduced_flux - apply(renamed_reduced_flux, 1, mean, na.rm = TRUE))/apply(renamed_reduced_flux, 1, sd, na.rm = TRUE)

#write.table(reduced_flux_mat, file = "carriedFlux.tsv", sep = "\t", row.names = TRUE, col.names = TRUE)

heatmap.2(std.reduced_flux_mat, Colv = FALSE, trace = "n", dendrogram = "row", cexRow = 0.05)

#look at the subset of rxns carrying flux in ura and leu

renamed_reduced_flux[apply(renamed_reduced_flux[,!(growth_rate$limit %in% c("Leucine", "Uracil"))] != 0, 1, sum) == 0,]


#heatmap of fluxes-per-unit growth

flux_per_gr <- reduced_flux_mat/matrix(growth_rate$growth, ncol = length(reduced_flux_mat[1,]), nrow = length(reduced_flux_mat[,1]), byrow = TRUE)
flux_per_gr <- flux_per_gr[apply(flux_per_gr, 1, sd) != 0,]

heatmap.2(t(scale(t(flux_per_gr), TRUE, TRUE)), trace = "n")

####### BRIDGE TO NETWORK LAYOUT ########



# generate the stoichiometry matrix from rxns carrying flux
load("totalStoi.Rdata")
first_write = FALSE

if(first_write == TRUE){
#writing a position file from scratch
Stotal <- S[,colnames(S) %in% rownames(reduced_flux_mat)]
Stotal <- Stotal[apply(Stotal != 0, 1, sum) != 0,]

metSty <- data.frame(SpeciesID = rep(NA, times = length(Stotal[,1])), SpeciesName = rep(NA, times = length(Stotal[,1])), Compartment = rep(NA, times = length(Stotal[,1])), x = rep(NA, times = length(Stotal[,1])), y = rep(NA, times = length(Stotal[,1]))) 
rxnSty <- data.frame(ReactionID = rep(NA, times = length(Stotal[1,])), Reaction = rep(NA, times = length(Stotal[1,])), Compartment = rep(NA, times = length(Stotal[1,])), x1 = rep(NA, times = length(Stotal[1,])), y1 = rep(NA, times = length(Stotal[1,])), x2 = rep(NA, times = length(Stotal[1,])), y2 = rep(NA, times = length(Stotal[1,]))) 

for(i in 1:length(Stotal[,1])){
	metSty[i,c(1:3)] <- corrFile[corrFile$SpeciesID %in% rownames(Stotal)[i],][,c(1,2,4)]
	}
for(i in 1:length(Stotal[1,])){
	if(colnames(Stotal)[i] %in% rxnFile$ReactionID){
		rxnSty[i,c(1:3)] <- rxnFile[rxnFile$ReactionID %in% colnames(Stotal)[i],][1,][c(2,1,3)]
		}else{
			rxnSty$ReactionID[i] <- colnames(Stotal)[i]
			}}
}else{
	#updating an existing position file
	metSty.old <- read.delim("metSty.tsv", sep = "\t", header = TRUE)
	rxnSty.old <- read.delim("rxnSty.tsv", sep = "\t", header = TRUE)
	
	Stotal <- S[,colnames(S) %in% union(rownames(reduced_flux_mat), rxnSty.old$ReactionID)]
	Stotal <- Stotal[apply(Stotal != 0, 1, sum) != 0,]
	
	#new rxns
	new_rxns <- rownames(reduced_flux_mat)[!(rownames(reduced_flux_mat) %in% rxnSty.old$ReactionID)]
	new_mets <- rownames(Stotal)[!(rownames(Stotal) %in% metSty.old$SpeciesID)]
	
	metSty_bind <- as.data.frame(matrix(NA, ncol = length(metSty.old[1,]), nrow = length(new_mets)))
	colnames(metSty_bind) <- colnames(metSty.old)
	for(i in 1:length(new_mets)){
		metSty_bind[i,c(1:3)] <- corrFile[corrFile$SpeciesID %in% new_mets[i],][,c(1,2,4)]
		}
	
	rxnSty_bind <- as.data.frame(matrix(NA, ncol = length(rxnSty.old[1,]), nrow = length(new_rxns)))
	colnames(rxnSty_bind) <- colnames(rxnSty.old)
	for(i in 1:length(new_rxns)){
	if(new_rxns[i] %in% rxnFile$ReactionID){
		rxnSty_bind[i,c(1:3)] <- rxnFile[rxnFile$ReactionID %in% new_rxns[i],][1,][c(2,1,3)]
		}else{
			rxnSty_bind$ReactionID[i] <- new_rxns[i]
			}}
	
	metSty = rbind(metSty.old, metSty_bind)
	rxnSty = rbind(rxnSty.old, rxnSty_bind)
	
	}

save(Stotal, metSty, rxnSty, reversibleRx, file = "totalStoiAux.Rdata")						
write.table(metSty, file = "metSty.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(rxnSty, file = "rxnSty.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)






######### compare boundary fluxes with their limit #########

limiting_fluxes <- matrix(nrow = length(names(treatment_par)), ncol = length(treatment_par[[1]]$nutrients[,1])); rownames(limiting_fluxes) <- names(treatment_par); colnames(limiting_fluxes) <- treatment_par[[1]]$nutrients[,1]

for(treatment in 1:length(names(treatment_par))){

limiting_fluxes[treatment,]	<-(reduced_flux_mat[sapply(sapply(treatment_par[[treatment]]$nutrients$nutrient, function(x){paste(x, "boundary")}), function(x){c(1:length(reduced_flux_mat[,1]))[rownames(reduced_flux_mat) %in% x]}), treatment]*-1)/treatment_par[[treatment]]$nutrients$conc_per_t
	
	}
#library(xtable)
#xtable(limiting_fluxes)


########## evaluating rxns ########

#query_rxns <- c("protein production")

#eval_mat <- stoiMat[apply(stoiMat[,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns] != 0, 1, sum) != 0,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))

######### evaluate metabolite #########

#eval_mat <- stoiMat[,rxnIDtoEnz(colnames(stoiMat)) %in% query_rxns]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))

reduced_flux_mat[rownames(reduced_flux_mat) %in% c("r_0393", "r_0394"),]


query_met <- c("stearoyl")

eval_mets(query_met, TRUE)

metToCHEBI("s_1002")
#
########






named_stoi <- stoiMat
met_dict <- metIDtoSpec(rownames(named_stoi))
rxn_dict <- rxnIDtoEnz(colnames(named_stoi))

rownames(named_stoi) <- sapply(c(1:length(named_stoi[,1])), function(x){met_dict[x][[1]]})
colnames(named_stoi) <- sapply(c(1:length(named_stoi[1,])), function(x){rxn_dict[x][[1]]})

aggregate_rxns <- c(colnames(rxn_search(named_stoi, "isa", is_rxn = TRUE)), colnames(rxn_search(named_stoi, "protein production", is_rxn = TRUE)))

##########


#flux through regexed rxn
renamed_reduced_flux[rownames(renamed_reduced_flux) %in% names(rxn_search(named_stoi, "g", is_rxn = FALSE)[1,]),]



growth_rate$growth

"uracil boundary"

growth_rate$growth[30]/(treatment_par$'Uracil 0.3'$nutrients$conc_per_t[treatment_par$'Uracil 0.3'$nutrients$nutrient == "uracil"]/as.numeric(compositionFile[compositionFile$MetName %in% "UMP",]$StoiCoef))



colnames(stoiMat)



#### looking at extracellular rxns

extra <- stoiMat[,compartment == "c_05"]
extra_stoi <- extra[apply(extra != 0, 1, sum) != 0,]

met_dict <- metIDtoSpec(rownames(extra_stoi))
rxn_dict <- rxnIDtoEnz(colnames(extra_stoi))

rownames(extra_stoi) <- sapply(c(1:length(extra_stoi[,1])), function(x){met_dict[x][[1]]})
colnames(extra_stoi) <- sapply(c(1:length(extra_stoi[1,])), function(x){rxn_dict[x][[1]]})

#look at all reactions involving a metabolite from a particular compartment, that actually carry flux

extra <- stoiMat[rownames(stoiMat) %in% corrFile[corrFile$Compartment == "c_02",]$SpeciesID,apply(stoiMat[rownames(stoiMat) %in% corrFile[corrFile$Compartment == "c_02",]$SpeciesID,] != 0, 2, sum) != 0]

extra_red <- stoiMat[,colnames(stoiMat) %in% colnames(extra)[colnames(extra) %in% rownames(reduced_flux_mat)]]
extra_red <- extra_red[apply(extra_red != 0, 1, sum) != 0,]

met_dict <- metIDtoSpec(rownames(extra_red))
rxn_dict <- rxnIDtoEnz(colnames(extra_red))

rownames(extra_red) <- sapply(c(1:length(extra_red[,1])), function(x){met_dict[x][[1]]})
colnames(extra_red) <- sapply(c(1:length(extra_red[1,])), function(x){rxn_dict[x][[1]]})

 
 

#fxns


	
	