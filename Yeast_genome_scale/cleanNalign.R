options(stringsAsFactors = FALSE)

setwd("~/Desktop/Rabinowitz/FBA_SRH/")
source("Yeast_genome_scale/FBA_lib.R")

inputFilebase = "yeast"
run_full = FALSE

rxnFile = read.delim(paste("Yeast_genome_scale/companionFiles/rxn_", inputFilebase, ".tsv", sep = ""))
rxnparFile = read.delim(paste("Yeast_genome_scale/companionFiles/species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE)
corrFile = read.delim(paste("Yeast_genome_scale/companionFiles/spec_", inputFilebase, ".tsv", sep = ""))
compFile <- read.delim(paste("Yeast_genome_scale/companionFiles/comp_", inputFilebase, ".tsv", sep = ""))

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]

if(file.exists("Yeast_genome_scale/yeast_stoi.R")){
	load("Yeast_genome_scale/yeast_stoi.R")
	} else {write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}

metabolites <- rownames(stoiMat)

#Determine the ChEBI IDs of all metabolites in the stoichiometric matrix
	
add.chebi.comp <- read.delim("Yeast_genome_scale/chebi_srh_curated.tsv")
add.chebi.comp <- add.chebi.comp[is.na(add.chebi.comp$generic),]
add.chebi.comp$internal_ID <- sapply(add.chebi.comp$internal_ID, function(x){corrFile$SpeciesType[corrFile$SpeciesID == x]})

yeastChEBI <- cbind(metabolites, as.numeric(unlist(sapply(metabolites, metToCHEBI))))

for(el in 1:length(add.chebi.comp[,1])){
	new_ids <- corrFile$SpeciesID[corrFile$SpeciesType == add.chebi.comp$internal_ID[el]]
	if(length(new_ids) == 1){
		yeastChEBI[c(1:length(yeastChEBI[,1]))[yeastChEBI[,1] == new_ids],2] <- add.chebi.comp$CHEBI_ID[el]
		}else{
			yeastChEBI[sapply(new_ids, function(x){c(1:length(yeastChEBI[,1]))[yeastChEBI[,1] == x]}),][,2] <- add.chebi.comp$CHEBI_ID[el]
			}
		}

#read in Yifan's absolute metabolite abundance data
yeastAbsConc <- read.delim("ChemicalSpeciesQuant/Metabolites/yeast_absolute_concentration_yifan.txt", sep = "\t")

#read in the metabolite abundance dataframe from Boer 2010
yeastChemoConc <- read.delim("ChemicalSpeciesQuant/brauer-microarray/BoerMetLevels.txt", header = FALSE, sep = "\t")
chemoMets <- unlist(yeastChemoConc[1,-c(1:7)])

#list constructed from chemoMets + manual ChEBI ID specification
yeastChemoIDs = read.delim("ChemicalSpeciesQuant/brauer-microarray/BoerMetIDs.csv", sep = ",", fill = TRUE)
chemoID.list <- list()
chemoMatch.list <- list()
	for(i in 1:length(yeastChemoIDs[,1])){
		chemoID.list[[i]] <- strsplit(yeastChemoIDs$ChEBI[i], split = ", ")[[1]]
		chemoMatch.list[[i]] <- yeastChEBI[,1][yeastChEBI[,2] %in% chemoID.list[[i]]]
	}
n.chemoMatch <- sapply(chemoMatch.list, function(x){length(x)})
mathedMets <- sort(unique(unlist(chemoMatch.list)))


#match protein ascertained by proteomics to the reactions they catalyze

yeastChemoProt <- read.delim("ChemicalSpeciesQuant/Proteomics/20120715_Scaffold_Files/ProteinReport20120714ChemostatYeastP.txt", skip = 42)
yeastChemoProt <- yeastChemoProt[!(yeastChemoProt[,1] == "END OF FILE"),]
samplez <- unique(yeastChemoProt$Biological.sample.name[yeastChemoProt$Biological.sample.name != ""])
orfz <- unique(yeastChemoProt$Protein.accession.numbers)
coverageMat <- matrix(NA, ncol = length(samplez), nrow = length(orfz))
colnames(coverageMat) <- samplez; rownames(coverageMat) <- orfz

for(i in 1:length(yeastChemoProt[,1])){
	coverageMat[rownames(coverageMat) == yeastChemoProt$Protein.accession.numbers[i], colnames(coverageMat) == yeastChemoProt$Biological.sample.name[i]] <- as.numeric(strsplit(yeastChemoProt$Percentage.sequence.coverage[i], "%")[[1]])/100
	}	
coverageMat[is.na(coverageMat)] <- 0


yeastRxnCoverage <- data.frame(rxn = colnames(stoiMat), coverageDepth = 0, coverageBreadth = 0)
for(ID in colnames(stoiMat)){
	geneNamez <- rxnFile[rxnFile$ReactionID == ID,]$MetName[is.na(rxnFile[rxnFile$ReactionID == ID,]$StoiCoef)]
	geneNamez <- unique(unlist(strsplit(geneNamez, split = ":")))
	if(length(geneNamez) == 0){
		next
		}
	rxnMatches <- matrix(sapply(geneNamez, function(ID){
		if(ID %in% rownames(coverageMat)){
		coverageMat[grep(ID, rownames(coverageMat)),]
			}else{
				rep(0, times = length(coverageMat[1,]))
				}}), ncol = length(geneNamez))
	
	bestMatch <- rxnMatches[,which.max(apply(rxnMatches, 2, mean))]
	yeastRxnCoverage$coverageBreadth[yeastRxnCoverage$rxn == ID] <- sum(bestMatch != 0)
	yeastRxnCoverage$coverageDepth[yeastRxnCoverage$rxn == ID] <- mean(bestMatch)
	}

save(mathedMets, yeastRxnCoverage, file = "ChemicalSpeciesQuant/speciesMatch.Rdata")




