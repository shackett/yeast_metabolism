setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
load('run_Deg/protLMfile.Rdata')
library(qvalue)
source("pep_library.R")

n_bs_t <- 5000
n_bs_f <- 1000
FDR_desired <- 0.05

cond_info <- data.frame(t(sapply(rownames(prot_abund_final), function(cond){
	if(strsplit(cond, "")[[1]][1] == "L"){
		tmp <- cond
		}else{
			tmp <- paste(c(tolower(strsplit(cond, "")[[1]][1]), strsplit(cond, "")[[1]][-1]), collapse = "")
			}
		t(conditions[conditions$condition == tmp,])
	})), stringsAsFactors = FALSE)
colnames(cond_info) <- colnames(conditions)	
	
cond_info$limitation <- factor(cond_info$limitation)
cond_info$actualDR <- as.numeric(cond_info$actualDR)
cond_info$logSF <- as.numeric(cond_info$logSF)

####### determine estimates of coefficients for all proteins ########
specAbund <- prot_abund[,1:n_prot]
specPrec <- prot_prec[,1:n_prot]
specAbund[as.matrix(specPrec == 0)] <- NA

prot_tstat <- bs_tstat(specAbund, specPrec, n_bs_t, FDR_desired)
prot_Fstat <- bs_Fstat(specAbund, specPrec, n_bs_f, FDR_desired)

######### Use gene-set enrichment analysis to determine if genes with a significant coefficient are enriched for ontology categories #######
######## Consider genes with a positive or negative test statistic seperately #######
 

library(GSEABase)
library(org.Sc.sgd.db)
library(methods)
library(GOstats)

frame = toTable(org.Sc.sgdGO)

goframeData = data.frame(frame$go_id, frame$Evidence, frame$systematic_name)
goFrame=GOFrame(goframeData,organism="Saccharomyces cerevisiae")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())




#MF - molecular function, BP - biological process, CC cellular component
GSEA.test = function(genes, category) {
	universe = sapply(colnames(specAbund), function(id){strsplit(id, split = '/')[[1]][1]})
	
	params <- GSEAGOHyperGParams(name="GSEA for proteomics data",
		geneSetCollection=gsc,
		geneIds = genes,
		universeGeneIds = universe,
		ontology = category,
		pvalueCutoff = 1,
		conditional = FALSE,
		testDirection = "over")
	hyperGTest(params)
	}

 


for(test in 1:length(prot_tstat$FDR_stats[,1])){
	if(prot_tstat$FDR_stats$discoveries[test] < 20){
		next
		}
	#split by +/- test statistic
	high_genes <- sapply(prot_tstat$Discoveries[[test]]$protein[prot_tstat$Discoveries[[test]]$change > 0], function(id){strsplit(id, split = '/')[[1]][1]})
	low_genes <- sapply(prot_tstat$Discoveries[[test]]$protein[prot_tstat$Discoveries[[test]]$change < 0], function(id){strsplit(id, split = '/')[[1]][1]})
	
	sum_table <- NULL
	
	for(ont in c("MF", "BP", "CC")){
		go_sum <- summary(GSEA.test(high_genes, ont))
		bonfp <- go_sum$Pvalue*length(go_sum$Pvalue)
		lowp <- go_sum[bonfp < 0.05,]
		if(length(lowp[,1]) != 0){
			colnames(lowp)[1] <- "GOID"
			sum_table <- rbind(sum_table, cbind(category = ont, stat_p = lowp, bonferroni_p = bonfp[bonfp < 0.05], dir = "UP"))
			}
		go_sum <- summary(GSEA.test(low_genes, ont))
		bonfp <- go_sum$Pvalue*length(go_sum$Pvalue)
		lowp <- go_sum[bonfp < 0.05,]
		if(length(lowp[,1]) != 0){
			colnames(lowp)[1] <- "GOID"
			sum_table <- rbind(sum_table, cbind(category = ont, stat_p = lowp, bonferroni_p = bonfp[bonfp < 0.05], dir = "DOWN"))
			}
		
	}
	prot_tstat$GO[[test]] <- sum_table
}


for(test in 1:length(prot_Fstat$FDR_stats[,1])){
	if(prot_Fstat$FDR_stats$discoveries[test] < 20){
		next
		}
	#split by +/- test statistic
	high_genes <- sapply(prot_Fstat$Discoveries[[test]]$protein[prot_Fstat$Discoveries[[test]]$change > 0], function(id){strsplit(id, split = '/')[[1]][1]})
	low_genes <- sapply(prot_Fstat$Discoveries[[test]]$protein[prot_Fstat$Discoveries[[test]]$change < 0], function(id){strsplit(id, split = '/')[[1]][1]})
	
	sum_table <- NULL
	
	for(ont in c("MF", "BP", "CC")){
		go_sum <- summary(GSEA.test(high_genes, ont))
		bonfp <- go_sum$Pvalue*length(go_sum$Pvalue)
		lowp <- go_sum[bonfp < 0.05,]
		if(length(lowp[,1]) != 0){
			colnames(lowp)[1] <- "GOID"
			sum_table <- rbind(sum_table, cbind(category = ont, stat_p = lowp, bonferroni_p = bonfp[bonfp < 0.05], dir = "UP"))
			}
		go_sum <- summary(GSEA.test(low_genes, ont))
		bonfp <- go_sum$Pvalue*length(go_sum$Pvalue)
		lowp <- go_sum[bonfp < 0.05,]
		if(length(lowp[,1]) != 0){
			colnames(lowp)[1] <- "GOID"
			sum_table <- rbind(sum_table, cbind(category = ont, stat_p = lowp, bonferroni_p = bonfp[bonfp < 0.05], dir = "DOWN"))
			}
		
	}
	prot_Fstat$GO[[test]] <- sum_table
}






####### determine estimates of coefficients for all divergent peptides ######
specAbund <- prot_abund[,(n_prot+1):length(prot_abund[1,])][,apply(prot_abund[,(n_prot+1):length(prot_abund[1,])] != 0, 2, sum) != 0]
specPrec <- prot_prec[,(n_prot+1):length(prot_abund[1,])][,apply(prot_abund[,(n_prot+1):length(prot_abund[1,])] != 0, 2, sum) != 0]
specAbund[as.matrix(specPrec == 0)] <- NA

divpep_tstat <- bs_tstat(specAbund, specPrec, n_bs_t, FDR_desired)
divpep_Fstat <- bs_Fstat(specAbund, specPrec, n_bs_f, FDR_desired)

protRegressions <- list(prot_tstat = prot_tstat, prot_Fstat = prot_Fstat, divpep_tstat = divpep_tstat, divpep_Fstat = divpep_Fstat)

save(protRegressions, n_bs_t, n_bs_f, FDR_desired, file = "protRegressions.Rdata")





