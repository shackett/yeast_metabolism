#source("http://bioconductor.org/biocLite.R")
#biocLite("org.Sc.sgd.db")
source("http://bioconductor.org/biocLite.R")
biocLite("org.EcK12.eg.db")

#all of the mapped genes in yeast
mapped_genes

library("org.Sc.sgd.db")
x <- org.Sc.sgdENZYME
    # Get the ORF identifiers that are mapped to an EC number 
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
  
    }    

library("org.EcK12.eg.db")

library("org.EcK12.eg.db")
y <- org.EcK12.egENZYME
	# Get the ORF identifiers that are mapped to an EC number 
    mapped_genes <- mappedkeys(y)
    # Convert to a list
    yy <- as.list(y[mapped_genes])
    if(length(yy) > 0) {
      # Get the ENZYME id for the first five genes
      yy[1:5]
      # Get the first one
      yy[[1]]
    }    

table(unique(unlist(yy)) %in% unique(unlist(xx)))
table(unique(unlist(xx)) %in% unique(unlist(yy)))




org.Sc.sgdENZYME()

setwd("~/Desktop/Rabinowitz/FBA_SRH/")
source("Yeast_genome_scale/FBA_lib.R")


options(stringsAsFactors = FALSE)

inputFilebase = "yeast"

rxnFile = read.delim(paste("Yeast_genome_scale/rxn_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
rxnparFile = read.delim(paste("Yeast_genome_scale/species_par_", inputFilebase, ".tsv", sep = ""), header = FALSE, stringsAsFactors = FALSE)
corrFile = read.delim(paste("Yeast_genome_scale/spec_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)
compFile <- read.delim(paste("Yeast_genome_scale/comp_", inputFilebase, ".tsv", sep = ""), stringsAsFactors = FALSE)

reactions = unique(rxnFile$ReactionID)
rxnStoi <- rxnFile[is.na(rxnFile$StoiCoef) == FALSE,]

if(file.exists("Yeast_genome_scale/yeast_stoi.R")){
	load("Yeast_genome_scale/yeast_stoi.R")
	} else {write_stoiMat(metabolites, reactions, corrFile, rxnFile, internal_names = TRUE)}

metabolites <- rownames(stoiMat)
#write out the yeast ChEBI IDs
yeastChEBI <- cbind(metabolites, as.numeric(unlist(sapply(metabolites, metToCHEBI))))
write.table(yeastChEBI, file = "yeast_chebi_id.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)

#write out the ecoli KEGG IDs 
e_coli_mets <- read.delim("gibbs_mets.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
e_coli_rxns <- read.delim("gibbs_rxns.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
e_coli_S <- as.matrix(read.delim("e_coliStoi", header = FALSE, sep = ",", stringsAsFactors = FALSE))

e_coli_mets[,1:3] <- apply(e_coli_mets[,1:3], c(1,2), function(x){
	if(x == "''"){NA}else{
	strsplit(x, split = "'")[[1]][2]
	}})

write.table(e_coli_mets[,c(1,3)], file = "ecoli_kegg_id.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)













#read in to get a list of yeast chebi:KEGG & chemical species synonyms
yeastCTS_1 <- read.table("yeast_CTS_output.csv", sep = ";", header = TRUE)

ChEBIsyn <- NULL
for(i in 1:length(yeastCTS_1[,1])){
	if(!(yeastCTS_1$Synonyms[i] %in% c("None", ""))){ 
		i_syn <- unlist(strsplit(yeastCTS_1$Synonyms[i], ", "))
		i_syn <- sub('[^:space:+]', "", i_syn)
		i_syn <- unname(sapply(i_syn, function(x){unlist(strsplit(x, " - "))[1]}))
		ChEBIsyn <- rbind(ChEBIsyn, data.frame(yChEBI = yeastCTS_1$From[i], ySyn = i_syn))
		}
	}
	
ChEBIsyn <- ChEBIsyn[ChEBIsyn[,2] != "",]
	
	
joint.ids <- sapply(c(1:length(ChEBIsyn[,1])), function(x){
	paste(ChEBIsyn[x,1], ChEBIsyn[x,2], sep = "fubar")
	})

joint.mat <- matrix(unlist(sapply(unique(joint.ids), function(x){strsplit(x, split = "fubar")})), ncol = 2, byrow = TRUE)
	
		
write.table(joint.mat[1:6000,2], file = "yeast_syn1.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(joint.mat[6001:length(joint.mat[,1]), 2], file = "yeast_syn2.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE)



#get a kegg id for each metabolite



keggChebiCorr <- read.table("Yeast_reconstruction/Sequences/database_accession.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
keggChebiCorr <- keggChebiCorr[grep("KEGG", keggChebiCorr$TYPE),]

met_dict <- data.frame("yeastID" = metabolites, "ChEBI" = as.numeric(unlist(sapply(metabolites, metToCHEBI))), "KEGG" = NA)

matched_mets <- met_dict$ChEBI[c(1:length(met_dict[,1]))[(met_dict$ChEBI %in% keggChebiCorr$COMPOUND_ID)]]
#keggChebiCorr[keggChebiCorr$COMPOUND_ID %in% matched_mets,]

met_dict$KEGG <- sapply(met_dict$ChEBI, function(x){
	if(x %in% keggChebiCorr$COMPOUND_ID){
		keggChebiCorr$ACCESSION_NUMBER[keggChebiCorr$COMPOUND_ID %in% x]
		}else{
			NA
			}})
		
	
	 

#cross-reference unique chemical names with ChEBI and KEGG using CTS to create a degenerate correspondence between these species identifiers
CTS_degen_files <- list.files(pattern = "yCTS")
unlist(strsplit(CTS_degen_files, split = "yCTS_"))
files_numz <- as.numeric(sapply(CTS_degen_files, function(x){unlist(strsplit(unlist(strsplit(x, split = "yCTS_"))[2], split = ".csv"))[1]}))

import_order <- CTS_degen_files[order(files_numz)]
degen_library <- NULL

for(a_file in import_order){
	extending_lib <- read.delim(a_file, sep = ";", header = TRUE, stringsAsFactors = FALSE)
	extending_lib <- extending_lib[,colnames(extending_lib) %in% c("From", "Name", "CAS", "CHEBI", "KEGG")]
	
	
	}






#crossreference the Ecoli IDs with ChEBI and CAS using CTS






e_coli_mets <- 

unique_rxns <- 

#reducing the size of the stoichiometric matrix by consolidating reactions with identical stoichiometry

n.rx = length(e_coli_rxns[,1]); n.met = length(e_coli_mets[,1])




#u.mat <- matrix(NA, ncol = n.rx, nrow = n.rx)
#for(i in 1:n.rx){
#	tmp <- matrix(e_coli_S[,i], ncol = n.rx, n.met) - e_coli_S
#	u.mat[i,] <- apply(tmp != 0, 2, sum) == 0
#	}
