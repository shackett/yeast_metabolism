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
	












#match based on shared stoichiometry

YE_match <- read.delim("EcoliYeastMatch/correspondence.dict.csv", sep = "\t", header = TRUE)
YE_match <- YE_match[sapply(c(1:length(YE_match[,1])), function(x){c(1:length(YE_match[,1]))[YE_match[,1] == rownames(stoiMat)[x]]}),]

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
	for(spec in 1:length(rxn_stoi)){
		if(length(YE_match_list[names(rxn_stoi)][[spec]]) == 0){
			print(paste("no match for", names(rxn_stoi)[spec], ":", metIDtoSpec(names(rxn_stoi)[spec])))
			}else{
				matching_mat <- rbind(matching_mat, matrix(e_coli_S[rownames(e_coli_S) %in% YE_match_list[names(rxn_stoi)][[spec]],], nrow = length(YE_match_list[names(rxn_stoi)][[spec]])))
				spec_col_id <- c(spec_col_id, rep(spec, times = length(YE_match_list[names(rxn_stoi)][[spec]])))
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
	
	near_matches[[colnames(stoiMat)[yRx]]] <- list()
	perfect_matches[[colnames(stoiMat)[yRx]]] <- list()
							
	
	if(length(rxn_stoi) > 1){
		pos_match <- apply(pos_match, 2, sum)
		neg_match <- apply(neg_match, 2, sum)
		
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
perfect <- c(1:length(is_perfect))[!is.na(is_perfect)][is_perfect[!is.na(is_perfect)]]
n_matches <- rep(NA, times = length(perfect))
rxn_free_e <- data.frame(delta_G_Kj_mol = rep(NA, times = length(perfect)), G_form_sd = rep(NA, times = length(perfect)))
for(match in perfect){
	n_matches[c(1:length(perfect))[perfect == match]] <- length(unlist(perfect_matches[[colnames(stoiMat)[match]]]))
	
	rxn_free_e <- 
	
	
	perfect_matches[[colnames(stoiMat)[match]]]
	
	e_coli_rxns[2047,]
	
	}

imperfect <- perfect[n_matches != 1]





	
		
		
	
	
		
	
	
	
	
	

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

e_coli_rxns

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




rxnFile












#all of the mapped genes in yeast
mapped_genes

library("org.Sc.sgd.db")
x <- org.Sc.sgdENZYME
    # Get the ORF identifiers that are mapped to an EC number 
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
  
    }    

table(unique(unlist(xx)) %in% unique(unlist(yy)))




org.Sc.sgdENZYME()







