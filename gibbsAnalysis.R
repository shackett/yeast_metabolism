options(stringsAsFactors = FALSE)

setwd("~/Desktop/Rabinowitz/FBA_SRH/")
source("Yeast_genome_scale/FBA_lib.R")


inputFilebase = "yeast"
run_full = FALSE

rxnFile = read.delim(paste("Yeast_genome_scale/rxn_", inputFilebase, ".tsv", sep = ""))
rxnparFile = read.delim(paste("Yeast_genome_scale/species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE)
corrFile = read.delim(paste("Yeast_genome_scale/spec_", inputFilebase, ".tsv", sep = ""))
compFile <- read.delim(paste("Yeast_genome_scale/comp_", inputFilebase, ".tsv", sep = ""))

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]

if(file.exists("Yeast_genome_scale/yeast_stoi.R")){
	load("Yeast_genome_scale/yeast_stoi.R")
	} else {write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}




metabolites <- rownames(stoiMat)

e_coli_mets <- read.delim("EcoliYeastMatch/gibbs_mets.txt", header = TRUE, sep = "\t")
e_coli_mets[,1:3] <- apply(e_coli_mets[,1:3], c(1,2), function(x){
	if(x == "''"){NA}else{
	strsplit(x, split = "'")[[1]][2]
	}})
e_coli_rxns <- read.delim("EcoliYeastMatch/gibbs_rxns.txt", header = TRUE, sep = "\t")
e_coli_rxns[,c(1,3)] <- apply(e_coli_rxns[,c(1,3)], c(1,2), function(x){
	if(x == "''"){NA}else{
	strsplit(x, split = "'")[[1]][2]
	}})


e_coli_S <- as.matrix(read.delim("EcoliYeastMatch/e_coliStoi", header = FALSE, sep = ",", stringsAsFactors = FALSE))
rownames(e_coli_S) <- e_coli_mets$ID
colnames(e_coli_S) <- e_coli_rxns$ID



if(run_full == TRUE){
	#write out the ecoli KEGG IDs 
	
	write.table(e_coli_mets[,c(1,3)], file = "ecoli_kegg_id.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

	#write out the yeast ChEBI IDs
	
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
	
	write.table(yeastChEBI, file = "yeast_chebi_id.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)
	
	#use the script ecoliYeastMatch.py to match the yeast species ids to ecoli species ids
	}	
	






#corrFile[,2][corrFile[,1] %in% c("s_1098", "s_1094")]
#YE_match[YE_match[,1] %in% c(corrFile[grep("water", corrFile[,2]),]$SpeciesID),]
#corrFile[grep("water", corrFile[,2]),]$SpeciesID


#match based on shared stoichiometry

YE_match <- read.delim("EcoliYeastMatch/correspondence.dict.csv", sep = "\t", header = TRUE)
YE_match <- YE_match[sapply(c(1:length(YE_match[,1])), function(x){c(1:length(YE_match[,1]))[YE_match[,1] == rownames(stoiMat)[x]]}),]

#correct some obvious mismatches

YE_match <- match_correction(YE_match, "water", "H2O")
YE_match <- match_correction(YE_match, "^ATP", "^ATP")
YE_match <- match_correction(YE_match, "^ADP$", "^ADP$")
YE_match <- match_correction(YE_match, "FAD$", "Flavin-adenine-dinucleotide-oxidized")
YE_match <- match_correction(YE_match, "FADH2", "Flavin-adenine-dinucleotide-reduced")
YE_match <- match_correction(YE_match, "^NAD\\(", "Nicotinamide-adenine-dinucleotide$")
YE_match <- match_correction(YE_match, "NADH", "Nicotinamide-adenine-dinucleotide--reduced")
YE_match <- match_correction(YE_match, "^NADP\\(", "Nicotinamide-adenine-dinucleotide-phosphate$")
YE_match <- match_correction(YE_match, "NADPH", "Nicotinamide-adenine-dinucleotide-phosphate--reduced")
#glucose
#f16bp
#D-ribose
#UTP

#metIDtoSpec(YE_match$yeastID)[grep("^NADP\\(", metIDtoSpec(YE_match$yeastID))]
#YE_match[grep("^NADP\\(", metIDtoSpec(YE_match$yeastID)),]
#e_coli_mets[grep("NADH", e_coli_mets$Name),]


#yeastReg <- "^ATP"
#coliReg <- "^ATP"

YE_match_list <- list()
for(spec in c(1:length(YE_match[,1]))){
	YE_match_list[[YE_match[spec,1]]] <- unlist(strsplit(YE_match[spec,2], split = ","))
	}



all_yeast_spec <- unique(corrFile$SpeciesType)
all_yeast_names <- sapply(all_yeast_spec, function(ID){
	corrFile$SpeciesName[corrFile$SpeciesType == ID][1]
	})


yeastMet_possible_match <- data.frame(specID = all_yeast_spec, specName = unname(all_yeast_names), ecoliMatchSugg = NA, ecoliMatchSuggManual = NA, carryFlux = NA)

for(i in 1:length(all_yeast_spec)){
	matched_spec <- corrFile$SpeciesID[corrFile$SpeciesType == all_yeast_spec[i]]
	matched_spec <- matched_spec[matched_spec %in% names(YE_match_list)]
	ecolimatches <- unique(unlist(YE_match_list[matched_spec]))
	if(length(ecolimatches) != 0){
		yeastMet_possible_match$ecoliMatchSugg[i] <- paste(ecolimatches, collapse = ", ")
		}
	}

#determine which mets are involved in rxns that carry flux
carriedFlux <- read.delim("Yeast_genome_scale/carriedFlux.tsv", sep = "\t")

flux_involved_spec <- unique(corrFile$SpeciesType[corrFile$SpeciesID %in% rownames(stoiMat)[apply(stoiMat[,colnames(stoiMat) %in% rownames(carriedFlux)] != 0, 1, sum) != 0]])

yeastMet_possible_match$carryFlux <- ifelse(yeastMet_possible_match$specID %in% flux_involved_spec, 1, 0)
yeastMet_possible_match <- yeastMet_possible_match[order(yeastMet_possible_match$carryFlux, decreasing = TRUE),]

#write.table(yeastMet_possible_match, file = "yeastMet_possible_match.tsv", row.names = FALSE, col.names = TRUE, sep = "\t")

yeastMet_possible_match <- read.delim("EcoliYeastMatch/yeastMet_possible_match.tsv", header = TRUE, sep = "\t")

yeastMet_possible_match$ecoliMatchSuggManual <- gsub("'", "", yeastMet_possible_match$ecoliMatchSuggManual)
yeastMet_possible_match$ecoliMatchSuggManual[is.na(yeastMet_possible_match$ecoliMatchSuggManual)] <- ""


for(overwrite in yeastMet_possible_match$specID[yeastMet_possible_match$carryFlux == 1 | yeastMet_possible_match$ecoliMatchSuggManual != ""]){
	
	YE_match$ecoliID[YE_match$yeastID %in% corrFile$SpeciesID[corrFile$SpeciesType == overwrite]] <- 
	yeastMet_possible_match$ecoliMatchSuggManual[yeastMet_possible_match$specID == overwrite]
	}

YE_match_list <- list()
for(spec in c(1:length(YE_match[,1]))){
	YE_match_list[[YE_match[spec,1]]] <- unlist(strsplit(YE_match[spec,2], split = ","))
	}









#for each yeast rxn: read the rows of the ecoli stoichiometry matrix corresponding to the species matched with yeast.  Look for matching stoichiometry

perfect_matches <- list()
near_matches <- list()
is_perfect <- rep(NA, times = length(stoiMat[1,]))

for(yRx in 1:length(stoiMat[1,])){
	
	rxn_stoi <- stoiMat[stoiMat[,yRx] != 0, yRx]
	if(length(rxn_stoi) == 1){
		names(rxn_stoi) <- rownames(stoiMat)[stoiMat[,yRx] != 0]
		}
			
	YE_match_list[names(rxn_stoi)]
	
	spec_col_id <- NULL
	matching_mat <- NULL
	involved_spec <- NULL
	for(spec in 1:length(rxn_stoi)){
		if(length(YE_match_list[names(rxn_stoi)][[spec]]) == 0){
			print(paste("no match for", names(rxn_stoi)[spec], ":", metIDtoSpec(names(rxn_stoi)[spec])))
			}else{
				matching_mat <- rbind(matching_mat, matrix(e_coli_S[rownames(e_coli_S) %in% YE_match_list[names(rxn_stoi)][[spec]],], nrow = length(YE_match_list[names(rxn_stoi)][[spec]])))
				spec_col_id <- c(spec_col_id, rep(spec, times = length(YE_match_list[names(rxn_stoi)][[spec]])))
				involved_spec <- c(involved_spec, rownames(e_coli_S)[rownames(e_coli_S) %in% YE_match_list[names(rxn_stoi)][[spec]]])
			}
		}
	
	
	
	
	if(is.null(spec_col_id)){
		next
		}
	
	unity_mat <- matching_mat * 1/matrix(sapply(spec_col_id, function(x){rxn_stoi[x]}), nrow = length(spec_col_id), ncol = length(e_coli_S[1,]))
	
	pos_match <- NULL
	neg_match <- NULL
	
	for(spec in 1:length(rxn_stoi)){
		if(length(spec_col_id[spec_col_id == spec]) == 0){
			pos_match <- rbind(pos_match, rep(FALSE, times = length(e_coli_S[1,])))
			neg_match <- rbind(neg_match, rep(FALSE, times = length(e_coli_S[1,])))
			}else{
				if(length(spec_col_id[spec_col_id == spec]) == 1){
					pos_match <- rbind(pos_match, unity_mat[c(1:length(spec_col_id))[spec_col_id == spec],] == 1)
					neg_match <- rbind(neg_match, unity_mat[c(1:length(spec_col_id))[spec_col_id == spec],] == -1)
					}else{
						pos_match <- rbind(pos_match, apply(unity_mat[c(1:length(spec_col_id))[spec_col_id == spec],] == 1, 2, sum) != 0)
						neg_match <- rbind(neg_match, apply(unity_mat[c(1:length(spec_col_id))[spec_col_id == spec],] == -1, 2, sum) != 0)
						
						}
					}
			}
	
	#match score = #of yeast species matched to an ecoli rxn - #of ecoli species not included
	missing_elements <- apply(e_coli_S[!(rownames(e_coli_S) %in% involved_spec),] != 0, 2, sum)
	
	near_matches[[colnames(stoiMat)[yRx]]] <- list()
	perfect_matches[[colnames(stoiMat)[yRx]]] <- list()
							
	
	if(length(rxn_stoi) > 1){
		pos_match <- apply(pos_match, 2, sum) - missing_elements
		neg_match <- apply(neg_match, 2, sum) - missing_elements
		
		most_similar <- max(c(pos_match, neg_match))
		
		if(most_similar == length(rxn_stoi)){
			is_perfect[yRx] <- TRUE
			if(length(c(1:length(pos_match))[pos_match == most_similar]) != 0){
				perfect_matches[[c(1:length(names(perfect_matches)))[names(perfect_matches) == colnames(stoiMat)[yRx]]]]$"F" <- c(1:length(pos_match))[pos_match == most_similar]
				}else{
					perfect_matches[[c(1:length(names(perfect_matches)))[names(perfect_matches) == colnames(stoiMat)[yRx]]]]$"F" <- NA
					}
			if(length(c(1:length(neg_match))[neg_match == most_similar]) != 0){
				perfect_matches[[c(1:length(names(perfect_matches)))[names(perfect_matches) == colnames(stoiMat)[yRx]]]]$"R" <- c(1:length(neg_match))[neg_match == most_similar]
				}else{
					perfect_matches[[c(1:length(names(perfect_matches)))[names(perfect_matches) == colnames(stoiMat)[yRx]]]]$"R" <- NA
					}
			
			}else{
				is_perfect[yRx] <- FALSE
				}
		
		if(most_similar >= length(rxn_stoi)/2){
			
			n_matches <- table(c(pos_match, neg_match))
			most_similar <- max(names(n_matches)[names(n_matches) != length(rxn_stoi)])
			
			if(length(c(1:length(pos_match))[pos_match == most_similar]) != 0){
				near_matches[[c(1:length(names(near_matches)))[names(near_matches) == colnames(stoiMat)[yRx]]]]$"F" <- c(1:length(pos_match))[pos_match == most_similar]
				}else{
					near_matches[[c(1:length(names(near_matches)))[names(near_matches) == colnames(stoiMat)[yRx]]]]$"F" <- NA
					}
			if(length(c(1:length(neg_match))[neg_match == most_similar]) != 0){
				near_matches[[c(1:length(names(near_matches)))[names(near_matches) == colnames(stoiMat)[yRx]]]]$"R" <- c(1:length(neg_match))[neg_match == most_similar]
				}else{
					near_matches[[c(1:length(names(near_matches)))[names(near_matches) == colnames(stoiMat)[yRx]]]]$"R" <- NA
					}
			
			}else{
				near_matches[[c(1:length(names(near_matches)))[names(near_matches) == colnames(stoiMat)[yRx]]]]$"F" <- NA
				near_matches[[c(1:length(names(near_matches)))[names(near_matches) == colnames(stoiMat)[yRx]]]]$"R" <- NA
				}
		}
	}
		


#breakdown of assignment
c(table(is_perfect), nNA = table(is.na(is_perfect))[names(table(is.na(is_perfect))) == "TRUE"])

#determine whether perfect matches are unique		


SM_gibbs <- data.frame(rxn = colnames(stoiMat), rxName = rxnIDtoEnz(colnames(stoiMat)), ecoliRxName = NA, delta_G_Kj_mol = NA, G_form_sd = NA)

perfect <- c(1:length(is_perfect))[!is.na(is_perfect)][is_perfect[!is.na(is_perfect)]]
n_matches <- rep(NA, times = length(perfect))
rxn_free_e <- data.frame(delta_G_Kj_mol = rep(NA, times = length(perfect)), G_form_sd = rep(NA, times = length(perfect)))
for(match in perfect){
	n_matches[c(1:length(perfect))[perfect == match]] <- sum(!is.na(unlist(perfect_matches[[colnames(stoiMat)[match]]])))
	
	
	if(!is.na(perfect_matches[[colnames(stoiMat)[match]]]$F)[1] & !is.na(perfect_matches[[colnames(stoiMat)[match]]]$R)[1]){
		#some reactions match with forward & reverse polarity
		F_gibbs <- e_coli_rxns[perfect_matches[[colnames(stoiMat)[match]]]$F,][,c(2,4:5)]
		R_gibbs <- e_coli_rxns[perfect_matches[[colnames(stoiMat)[match]]]$R,][,c(2,4:5)]
		R_gibbs[,2] <- R_gibbs[,2]*-1
		freeE <- rbind(F_gibbs, R_gibbs)
		}else{
			if(!is.na(perfect_matches[[colnames(stoiMat)[match]]]$F)[1]){
				#only forward matches
				freeE <- e_coli_rxns[perfect_matches[[colnames(stoiMat)[match]]]$F,][,c(2,4:5)]
				
				}else{
					#only reverse matches
					R_gibbs <- e_coli_rxns[perfect_matches[[colnames(stoiMat)[match]]]$R,][,c(2,4:5)]
					R_gibbs[,2] <- R_gibbs[,2]*-1
					freeE <- R_gibbs
					}
			}
	
	if(length(freeE[,1]) == 1){
		SM_gibbs[match,3:5] <- c(e_coli_rxns$"Name"[as.numeric(rownames(freeE))], freeE$delta_G_Kj_mol, freeE$G_form_sd)
		next
		}
	
	if(sd(freeE$delta_G_Kj_mol) == 0){
		SM_gibbs[match,3:5] <- c(e_coli_rxns$"Name"[as.numeric(rownames(freeE)[1])], mean(freeE$delta_G_Kj_mol), max(freeE$G_form_sd))
		next
		}
		
	if(abs(sd(freeE$delta_G_Kj_mol)/mean(freeE$delta_G_Kj_mol)) < 0.2){
		SM_gibbs[match,3:5] <- c(e_coli_rxns$"Name"[as.numeric(rownames(freeE)[1])], mean(freeE$delta_G_Kj_mol), max(freeE$G_form_sd))
		}
	}
	
	
	table(!is.na(SM_gibbs$delta_G_Kj_mol))


	
		
#match based on shared E.C. numbers

e_coli_genes <- read.table("EcoliYeastMatch/gibbs_genes.txt", sep = "\t")
e_coli_rxn_corr <- read.table("EcoliYeastMatch/rxnGeneMat", sep = ",")
#ecoli_dict <- read.table("EcoliYeastMatch/ecoliNameDict.txt", sep = "\t", header = TRUE)
ecoli_dict <- read.table("EcoliYeastMatch/gene-links.dat", skip = 16, sep = "\t")

e_coli_gene_names <- unlist(sapply(e_coli_genes[,1], function(ID){ifelse(ID %in% ecoli_dict[,3], ecoli_dict[,7][ecoli_dict[,3] == ID], NA)}))

#source("http://bioconductor.org/biocLite.R")
#biocLite("org.EcK12.eg.db")

library("org.EcK12.eg.db")
# Convert the object to a list
name_to_entrez_list <- as.list(org.EcK12.egALIAS2EG)

name_to_entrez <- sapply(e_coli_gene_names, function(ID){
	if(is.na(ID)){
		NA
		}else{
		if(ID %in% names(name_to_entrez_list)){
			name_to_entrez_list[ID]
			}else{
				NA
				}
		}
	})

entrez_to_EC_list <- as.list(org.EcK12.egENZYME)

entrez_to_EC <- sapply(name_to_entrez, function(ID){
	if(is.na(ID[1])){
		NA
		}else{
			ECnums <- NULL
			for(i in 1:length(ID)){
				if(ID[i] %in% names(entrez_to_EC_list)){
					ECnums <- union(ECnums, entrez_to_EC_list[names(entrez_to_EC_list) == ID[i]])
					}
				}
			if(length(ECnums) == 0){
				NA}else{
					ECnums
					}
			}
	})

### Pass free energy to E.C. numbers

unique_ECs <- unique(unlist(entrez_to_EC)[!is.na(unlist(entrez_to_EC))])
e_coli_EC_gibbs <- list()
for(ID in unique_ECs){
	e_coli_EC_gibbs[ID] <- NA
	}


for(gene in 1:length(e_coli_rxn_corr[1,])){
	
	ECvals <- unlist(entrez_to_EC[[gene]])
	if(!is.na(ECvals[1])){
		for(ECval in ECvals){
			if(is.vector(e_coli_EC_gibbs[[ECval]])){
				e_coli_EC_gibbs[[ECval]] <- e_coli_rxns[e_coli_rxn_corr[,gene] == 1,]
				}else{
					e_coli_EC_gibbs[[ECval]] <- rbind(e_coli_EC_gibbs[[ECval]], e_coli_rxns[e_coli_rxn_corr[,gene] == 1,])
					}
			}
		}
	}

#determine the EC numbers of yeast genes

library("org.Sc.sgd.db")
yeast_gene_to_EC_list <- as.list(org.Sc.sgdENZYME)

yeastEC <- list()
for(ID in colnames(stoiMat)){
	geneNamez <- rxnFile[rxnFile$ReactionID == ID,]$MetName[is.na(rxnFile[rxnFile$ReactionID == ID,]$StoiCoef)]
	ECnums <- unlist(yeast_gene_to_EC_list[unique(unlist(strsplit(geneNamez, split = ":")))])
	if(is.null(ECnums)){yeastEC[[ID]] <- NA; next}
	if(length(ECnums[!is.na(ECnums)]) != 0){
		yeastEC[[ID]] <- ECnums[!is.na(ECnums)]
		}else{yeastEC[[ID]] <- NA}
	}

#take the reaction free energies attributed to EC numbers and pass them to the yeast rxns with the same EC number

yeast_EC_pass <- list()
for(ID in colnames(stoiMat)){
	if(is.na(yeastEC[[ID]][1])){
		yeast_EC_pass[[ID]] <- NA; next
		}
	if(sum(yeastEC[[ID]] %in% names(e_coli_EC_gibbs)) == 0){
		yeast_EC_pass[[ID]] <- NA; next
		}
		
	matchez <- yeastEC[[ID]][yeastEC[[ID]] %in% names(e_coli_EC_gibbs)]
	out <- NULL
	
	for(match in matchez){out <- rbind(out, e_coli_EC_gibbs[[match]])}
	
	yeast_EC_pass[[ID]] <- out
		
	}

#make a table of the high confidence gibbs free energies drawn by E.C. comparison

EC_gibbs <- data.frame(rxn = names(yeast_EC_pass), rxName = rxnIDtoEnz(names(yeast_EC_pass)), ecoliRxName = NA, delta_G_Kj_mol = NA, G_form_sd = NA)

for(rx in 1:length(EC_gibbs[,1])){
	if(is.vector(yeast_EC_pass[[EC_gibbs$rxn[rx]]])){
		next
		}
	matchez <- yeast_EC_pass[[EC_gibbs$rxn[rx]]]
	if(length(matchez$delta_G_Kj_mol) == 1){
		EC_gibbs[rx,3:5] <- matchez[,3:5]
		next
		}
	if(sd(matchez$delta_G_Kj_mol) == 0){
		EC_gibbs[rx,3:5] <- c(matchez$Name[1], mean(matchez$delta_G_Kj_mol), max(matchez$G_form_sd))
		next
		}
	if(abs(sd(matchez$delta_G_Kj_mol)/mean(matchez$delta_G_Kj_mol)) < 0.2){
		EC_gibbs[rx,3:5] <- c(matchez$Name[1], mean(matchez$delta_G_Kj_mol), max(matchez$G_form_sd))
		}
	}


dim(EC_gibbs[!is.na(EC_gibbs$ecoliRxName),])
#288 1-1 matches by E.C. alignment
#cbind(EC_gibbs$rxName, EC_gibbs$ecoliRxName)[!is.na(EC_gibbs$ecoliRxName),][1:20,]

table(apply(!is.na(cbind(SM_gibbs[,4], EC_gibbs[,4])), 1, sum))
#about 700 IDs total

EC_gibbs[,4:5] <- apply(EC_gibbs[,4:5], c(1,2), as.numeric)
SM_gibbs[,4:5] <- apply(SM_gibbs[,4:5], c(1,2), as.numeric)

gibbsFused <- data.frame(rxn = names(yeast_EC_pass), rxName = rxnIDtoEnz(names(yeast_EC_pass)), delta_G_Kj_mol = NA, G_form_sd = NA)
for(i in 1:length(gibbsFused[,1])){
	if(!is.na(SM_gibbs[i,4]) & !is.na(EC_gibbs[i,4])){
		gibbsFused[i,3:4] <- c(mean(SM_gibbs[i,4], EC_gibbs[i,4]), max(SM_gibbs[i,5], EC_gibbs[i,5]))
		next
		}
	if(!is.na(SM_gibbs[i,4])){
		gibbsFused[i,3:4] <- SM_gibbs[i,4:5]
		next
		}	
	if(!is.na(EC_gibbs[i,4])){
		gibbsFused[i,3:4] <- EC_gibbs[i,4:5]
		next
		}
	}

	
table(abs(gibbsFused[,3]) > 30)

gibbsCO = 30
reversible <- unlist(sapply(c(1:length(gibbsFused[,1])), function(x){
	if(!is.na(gibbsFused[x,3])){
	if(gibbsFused[x,3] - 2*gibbsFused[x,4] > gibbsCO){
		-1
		}else{
	if(gibbsFused[x,3] + 2*gibbsFused[x,4] < -1*gibbsCO){
		1
		}else{
			0
			}	
		}}else{
			0
			}}))

write.table(data.frame(rx = gibbsFused[,1], reversible = reversible), file = "revRxns.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)




match_correction <- function(YE_match, yeastReg, coliReg, report = FALSE){
	
	#remove matches in other rxns
	for(id_match in e_coli_mets[grep(coliReg, e_coli_mets$Name),]$ID){
		id_match <- sub("\\[", "\\\\[", id_match)
		id_match <- sub("\\]", "\\\\]", id_match)
		YE_match$ecoliID <- sub(paste("^", id_match, ",", sep = ""), "", YE_match$ecoliID)
		YE_match$ecoliID <- sub(paste(",", id_match, sep = ""), "", YE_match$ecoliID)
		YE_match$ecoliID <- sub(paste("^", id_match, "$", sep = ""), "", YE_match$ecoliID)
		}
		
		
	#add appropriate mathches
	for(match in grep(yeastReg, metIDtoSpec(YE_match$yeastID))){
		YE_match$ecoliID[match] <- paste(e_coli_mets[grep(coliReg, e_coli_mets$Name),]$ID, collapse = ",")
		}
	
	print(paste("Replaced ", yeastReg, " ", length(grep(yeastReg, metIDtoSpec(YE_match$yeastID))), " times", sep = ""))
	
	if(report == TRUE){
		YE_match[c(grep(paste("^", id_match, ",", sep = ""), YE_match$ecoliID), grep(paste(",", id_match, sep = ""), YE_match$ecoliID), grep(paste("^", id_match, "$", sep = ""), YE_match$ecoliID)),]
		}
	YE_match
	}








