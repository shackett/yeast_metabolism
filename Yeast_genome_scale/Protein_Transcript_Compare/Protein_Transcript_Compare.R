#### This script is a pared down version of Pep_Prot.Rnw focusing on direct comparisons of proteomic and transcriptional data ####

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/Protein_Transcript_Compare")

options(stringsAsFactors = F)

### Packages ####
library(impute)
library(missMDA)
library(grid)
library(data.table)
library(ggplot2)
library(reshape2)
library(scales)
library(gridExtra) # for combining ggplot objects

load("../../ChemicalSpeciesQuant/Proteomics/EMoutputDeg.Rdata")
load("../../ChemicalSpeciesQuant/Proteomics/EMimport.Rdata")
source("../../ChemicalSpeciesQuant/Proteomics/pep_library.R")


prot_abund_final <- prot_abund
mixing_fract[likdiff_df$matched == "Divergent-trend",] <- 0

#remove proteins that were subsumable or poorly fit
prot_abund_final <- prot_abund_final[,!(colSums(mixing_fract) == 0)]
prot_prec <- prot_prec[,!(colSums(mixing_fract) == 0)]
possibleMap <- Matrix(unique_mappingMat[,alpha_pres == 1] == 1)
mixing_fract <- mixing_fract[,!(colSums(mixing_fract) == 0)]

prot_abund_final[prot_abund_final == 0] <- NA

n_prot <- ncol(prot_abund_final)

impute_abund  <- impute.knn(as.matrix(t(prot_abund_final)), k = 10, rowmax = 1)$data

impute_svd <- svd(impute_abund, nu = 10, nv = 10)
resid_wrt_svd <- impute_svd$u[,1:10] %*% diag(impute_svd$d[1:10]) %*% t(impute_svd$v[,1:10])
squared_resid <- (impute_abund - resid_wrt_svd)^2
squared_resid[as.matrix(t(prot_prec) == 0)] <- NA
realized_prec <- (1 / apply(squared_resid, 1, max, na.rm = T))

# for imputed abundances, set the precision equal to the maximal reconstruction residual (actual - UDt(V)) for that protein #
impute_prec <- t(prot_prec)
for(i in 1:nrow(impute_prec)){
  impute_prec[i,][impute_prec[i,] == 0] <- min(min(impute_prec[i,][impute_prec[i,] != 0]), realized_prec[i])
  }

write.output(impute_abund, "Output/proteinAbundance.tsv")
write.table(as.data.frame(as.matrix(impute_prec)), "Output/proteinPrecision.tsv", col.names = T, row.names = T, quote = F)

#### row-Center relative abundances ####

impute_abund <- impute_abund - apply(impute_abund, 1, mean)

#### Load Transcriptional data ####

transcript_brauer <- read.delim("dilution_rate_00_raw.cdt") # supplemental data in Brauer et al. 2008
rownames(transcript_brauer) <- transcript_brauer$YORF
transcript_brauer <- transcript_brauer[-1,-c(1:5)]

#remove genes with many missing values
#impute missing values with KNN imputation (as per original Brauer paper)

transcript_brauer <- transcript_brauer[rowSums(!is.na(transcript_brauer)) > 0.5*length(transcript_brauer[1,]),]
transcript_brauer <- impute.knn(as.matrix(transcript_brauer), k = 10)$data

transcript_brauer <- transcript_brauer[,grep('^S', colnames(transcript_brauer), invert = T)] # disregard sulfur limitation

transcript.condition <- read.delim("../../ChemicalSpeciesQuant/Proteomics/TranscriptomicsChemostatsActualDR.tsv") # import transcriptomics chemostat dilution rates
load('~/Desktop/Composition/RNA_abundance/RNAabundance.R') # import proteomics chemostat dilution rates


#missing mappings
missing_on_array <- rownames(impute_abund)[!(rownames(impute_abund) %in% rownames(transcript_brauer))]
length(missing_on_array)
genes_to_compare <- rownames(impute_abund)[rownames(impute_abund) %in% rownames(transcript_brauer)]

#determine the relative abundance of array data by looking at the corresponding condition and imputing unobserved growth rates by drawing a line between the abundance of flanking conditions

proteomics_reduced <- impute_abund[chmatch(genes_to_compare, rownames(impute_abund)),]
transcript_brauer_reduced <- transcript_brauer[chmatch(genes_to_compare, rownames(transcript_brauer)),]
all(rownames(proteomics_reduced) == rownames(transcript_brauer_reduced))

#load true dilution rates for each of the proteomics samples

prot_cond <- as.data.frame(matrix(NA, ncol = n_c, nrow = 2))
prot_cond <- sapply(colnames(impute_abund), function(cond){
	c(unlist(strsplit(cond, '[0-9]+'))[1], unlist(strsplit(cond, '[a-zA-Z]'))[2])
	})
colnames(prot_cond) <- colnames(impute_abund)
rownames(prot_cond) <- c("limitation", "GR")
prot_cond <- as.data.frame(cbind(t(prot_cond), realDR = NA), stringsAsFactors = FALSE)
prot_cond$realDR <- sapply(1:length(prot_cond[,1]), function(cond){
	RNAabund$actual.dr[grep(paste(prot_cond$limitation[cond], prot_cond$GR[cond], sep = ""), RNAabund$condition, ignore.case = TRUE)[1]]
	})

DR_change_mat <- matrix(0, nrow = nrow(transcript.condition), ncol = nrow(prot_cond))
colnames(DR_change_mat) <- rownames(prot_cond); rownames(DR_change_mat) <- colnames(transcript_brauer)
for(cond in 1:length(prot_cond[,1])){
	#find the 2 closest DR within the same limitation
	c_match <- c(1:nrow(transcript.condition))[transcript.condition$Nutrient %in% toupper(prot_cond$limitation[cond])]
	flanking_match <- c_match[order(abs(as.numeric(transcript.condition$GR[c_match]) - prot_cond[cond,]$realDR))[1:2]]
	lb_diff <- (prot_cond$realDR[cond] - as.numeric(transcript.condition$GR[flanking_match])[1])/diff(as.numeric(transcript.condition$GR[flanking_match]))
	DR_change_mat[flanking_match,cond] <- c((1-lb_diff), lb_diff)
	}
	
remapped_transc <- transcript_brauer_reduced %*% DR_change_mat

transcript_brauer_reduced <- remapped_transc - apply(remapped_transc, 1, mean) # mean center each genes expresion



#### Compare protein RA to transcript RA across all genes ####

pool_thresh <- 4
shared_prot_pool <- proteomics_reduced; shared_prot_pool[shared_prot_pool < -1*pool_thresh] <- -1*pool_thresh; shared_prot_pool[shared_prot_pool > pool_thresh] <- pool_thresh
shared_trans_pool <- transcript_brauer_reduced; shared_trans_pool[shared_trans_pool < -1*pool_thresh] <- -1*pool_thresh; shared_trans_pool[shared_trans_pool > pool_thresh] <- pool_thresh

PTcor <- signif(c(cor(c(proteomics_reduced), c(transcript_brauer_reduced), method = "pearson"), cor(c(proteomics_reduced), c(transcript_brauer_reduced), method = "spearman")), 3)

PTscatterDF <- data.table(rbind(data.frame(melt(shared_trans_pool), type = "Transcript"), data.frame(melt(shared_prot_pool), type = "Protein")))
colnames(PTscatterDF) <- c("Gene", "Condition", "RA", "Type")
PTscatterDF <- PTscatterDF[,list(Protein = RA[Type == "Protein"], Transcript = RA[Type == "Transcript"]), by = c("Gene", "Condition")]


hex_bin_max <- 500
n_hex_breaks <- 9
hex_breaks <- c(1,2,4,8,16,32,64,125,250,500)
#hex_breaks <- round(exp(seq(0, log(hex_bin_max), by = log(hex_bin_max)/n_hex_breaks)))

hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "gray90"), legend.position = "top", 
                   panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line"), axis.text = element_text(color = "black")) 

PTscatter_plotter <- ggplot(cbind(PTscatterDF, z = 1), aes(x = Protein, y = Transcript, z = z)) + scale_x_continuous("Relative protein abundance", expand = c(0.02,0.02)) + scale_y_continuous("Relative transcript abundance", expand = c(0.02,0.02)) + 
  scale_fill_gradient(name = "Counts", low = "black", high = "firebrick1", trans = "log", breaks = hex_breaks, labels = hex_breaks) + hex_theme
PTscatter_plotter <- PTscatter_plotter + geom_hex(bins = 60)
ggsave(plot = PTscatter_plotter, "Output/globalPTcomp.pdf", width = 10, height = 10)

Ptaggreement_all <- table(PTscatterDF$Protein > 0, PTscatterDF$Transcript > 0)
sum(diag(Ptaggreement_all)) / sum(Ptaggreement_all)
PTagreement <- table(PTscatterDF$Protein[abs(PTscatterDF$Transcript) > 1] > 0,
      PTscatterDF$Transcript[abs(PTscatterDF$Transcript) > 1] > 0)
sum(diag(PTagreement)) / sum(PTagreement)

#### Compare protein RA to transcript RA on a gene-by-gene basis ####

pt_corrs <- sapply(c(1:nrow(proteomics_reduced)), function(row){
  cor(proteomics_reduced[row,], transcript_brauer_reduced[row,], method = "spearman")
  }); pt_corrs <- data.frame(correlation = pt_corrs)


hist_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), legend.position = "none", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), axis.title.x = element_text(vjust = -0.1), axis.title.y = element_text(vjust = 0.3),
  axis.text = element_text(color = "black")) 

plotted_count_max <- max(table(cut(pt_corrs$correlation,seq(min(pt_corrs$correlation),max(pt_corrs$correlation),0.04)))) 
                           
cor_plot <- ggplot(pt_corrs, aes(x = correlation)) + hist_theme
cor_plot <- cor_plot + scale_x_continuous(name = "Spearman Correlation", expand = c(0,0), limits = c(-1,1), breaks = seq(-1, 1, 0.2)) + scale_y_discrete(name = "Counts", expand = c(0,0), breaks = (0:ceiling(plotted_count_max/5))*5) + geom_histogram(colour = "aliceblue", fill = "limegreen", binwidth = 0.04)
ggsave(plot = cor_plot, "Output/PTspearmans.pdf", width = 12, height = 12)


#generate a null distribution of spearman correlations via permutation testing
perms <- 1000
#generate a matrix of permutation indices
spearCorrNULL <- matrix(NA, ncol = perms, nrow = nrow(proteomics_reduced))
for(gene in 1:nrow(proteomics_reduced)){
  #permute protein pattern w.r.t. transcript pattern
  sampleIndices <- sapply(1:perms, function(perm){sample(1:ncol(proteomics_reduced), ncol(proteomics_reduced))})
  permutedData <- apply(sampleIndices, 2, function(permDat){proteomics_reduced[gene,permDat]})
  spearCorrNULL[gene,] <- apply(permutedData, 2, function(nullCorr){cor(nullCorr, transcript_brauer_reduced[gene,], method = "spearman")})
}

nullSpears <- melt(spearCorrNULL); colnames(nullSpears) <- c("gene", "permutation", "correlation")
cor_plotNULL <- ggplot(nullSpears, aes(x = correlation)) + hist_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
cor_plotNULL <- cor_plotNULL + scale_x_continuous(name = "Spearman Correlation", expand = c(0,0), limits = c(-1,1), breaks = seq(-1, 1, 0.2)) + scale_y_discrete(name = "Counts", expand = c(0,0)) + geom_histogram(colour = "aliceblue", fill = "chocolate1", binwidth = 0.04)
ggsave(plot = cor_plotNULL, "Output/PTspearmansNULL.pdf", width = 12, height = 12)


#two tailed test using permutations to establish which correlations are discoveries
library(qvalue)

spearPvalue <- sapply(1:nrow(pt_corrs), function(spearp){
  sum(pt_corrs$correlation[spearp] <= spearCorrNULL[spearp,])/perms
  })
spearPvalue_twotailed <- 1 - abs((0.5 - spearPvalue)*2)
corrQval <- qvalue(spearPvalue_twotailed) 
corrQval$pi0  

pt_corrs$significant <- corrQval$qvalues <= 0.05

cor_plot <- ggplot(pt_corrs, aes(x = correlation, fill = significant)) + hist_theme
cor_plot <- cor_plot + scale_x_continuous(name = "Spearman Correlation", expand = c(0,0), limits = c(-1,1), breaks = seq(-1, 1, 0.2)) + scale_y_discrete(name = "Counts", expand = c(0,0), breaks = (0:ceiling(plotted_count_max/5))*5) + geom_histogram(colour = "aliceblue", binwidth = 0.04) + scale_fill_manual(name = "Correlation Significance", values = c("TRUE" = "limegreen", "FALSE" = "orange"))
ggsave(plot = cor_plot, "Output/PTspearmansSig.pdf", width = 12, height = 12)


#### Investigating systematic departures between proteins and transcripts


cond_reorder <- c(16:20, 1:5, 11:15, 6:10, 21:25)
ptDiff <- proteomics_reduced[,cond_reorder] - transcript_brauer_reduced[,cond_reorder]
ptDiff <- ptDiff - apply(ptDiff, 1, mean)

### Save so that dendrogram of PTdiff can be generated and each matrix can be kmeans clustered so that shared motifs can be found

write.output(ptDiff, "Output/coclust_ptdiff.tsv")
write.output(proteomics_reduced[,cond_reorder], "Output/proteinRA.tsv")
write.output(transcript_brauer_reduced[,cond_reorder], "Output/transcriptRA.tsv")

#### Generate a list of by-gene comparisons of transcripts and proteins ####

PTcomparison <- list()

for(i in 1:nrow(proteomics_reduced)){
  Gene = rownames(proteomics_reduced)[i]
  
  PTdt <- data.table(Condition = toupper(colnames(proteomics_reduced)[cond_reorder]), Protein = proteomics_reduced[i,cond_reorder], Transcript = transcript_brauer_reduced[i,cond_reorder])
  PTdt[,Proteins_per_Transcript := Protein - Transcript,]
  
  corrected_corr <- cor(PTdt$Protein, PTdt$Transcript, method = "pearson")/sqrt(var(PTdt$Protein)/(var(PTdt$Protein) + median(squared_resid[rownames(squared_resid) == Gene,], na.rm = T)))
  corrected_corr <- ifelse(corrected_corr > 0, min(corrected_corr, 1), max(corrected_corr, -1))
  
  PTcomparison[[Gene]]$PTabund <- PTdt
  PTcomparison[[Gene]]$PTcorr <- corrected_corr
  
  }

save(PTcomparison, file = "../companionFiles/PTcomparison_list.Rdata")



####

spearCounts <- 3

gene_index <- c(order(pt_corrs$correlation)[1:spearCounts],
order(abs(pt_corrs$correlation))[1:spearCounts],
order(pt_corrs$correlation)[(length(pt_corrs[,1]) - (spearCounts - 1)):length(pt_corrs[,1])])
                             
dispartecors <- data.frame(Class = rep(c("anticorrelated", "uncorrelated", "correlated"), each = spearCounts), Corr = pt_corrs$correlation[gene_index], Gene = rownames(proteomics_reduced)[gene_index])

reducedPTscatter <- PTscatterDF[Gene %in% dispartecors$Gene,,]
reducedPTscatter[,Limitation := toupper(substr(Condition, 1, 1))]
reducedPTscatter$Limitation <- factor(reducedPTscatter$Limitation, levels = c("P", "C", "N", "L", "U"))
reducedPTscatter$Gene <- factor(reducedPTscatter$Gene, levels = dispartecors$Gene)
reducedPTscatter$GeneCommon <- orf2common(as.character(reducedPTscatter$Gene))
setkeyv(reducedPTscatter, c("Gene", "Limitation"))
reducedPTscatter$GeneCommon <- factor(reducedPTscatter$GeneCommon, levels = unique(reducedPTscatter$GeneCommon))


scatter_facet_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), 
                             legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major = element_line(size = 0.25, color = "black", linetype = "longdash"), 
                             axis.text.x = element_text(angle = 90),strip.background = element_rect(fill = "cadetblue2"), axis.text = element_text(color = "black")) 
                                                      
divScatterPlot <- ggplot(reducedPTscatter, aes(x = Transcript, y = Protein, col = Limitation, fill = Limitation)) + scatter_facet_theme  + facet_wrap(~ GeneCommon, ncol = spearCounts, scales = "free") +
  geom_point(aes(size = 3)) + scale_x_continuous(name = expression(log[2] ~ "Transcript relative abundance")) + scale_y_continuous(name = expression(log[2] ~ "Protein relative abundance")) +
  scale_colour_brewer(guide = "none", palette = "Set1") + scale_fill_brewer(name = "Limitation", palette = "Set1") + scale_size_continuous(guide = "none") + scale_alpha_continuous(guide = "none", range = c(0,0)) +
  geom_smooth(method = lm, aes(fill = NULL, alpha = as.numeric(0)), size = 1.5)
divScatterPlot
ggsave("Output/PTgoodBadCompare.pdf", height = 10, width = 10)





#### Determine how many principal components are prominent in the proteomics data based upon cross-validation and reconstruction error ####

pcrange <- c(2,18)
npc.compare <- estim_ncpPCA(prot_abund_final, ncp.min = pcrange[1], ncp.max = pcrange[2], method.cv = 'Kfold', pNA = 0.10, nbsim = 10)
npc <- (pcrange[1]:pcrange[2])[npc.compare$criterion < (max(npc.compare$criterion) - min(npc.compare$criterion))*0.01 + min(npc.compare$criterion)][1]


ScreePlots <- data.frame(PC = 1:25, T = svd(transcript_brauer_reduced)$d^2 / sum(svd(transcript_brauer_reduced)$d^2), P = svd(proteomics_reduced)$d^2 / sum(svd(proteomics_reduced)$d^2), PTratio = svd(ptDiff)$d^2 / sum(svd(ptDiff)$d^2))
ScreePlots <- melt(ScreePlots, id.vars = "PC", variable.name = "DataType", value.name = "VarianceFraction")

scatter_theme <- theme(text = element_text(size = 23), title = element_text(size = 25), panel.background = element_rect(fill = "white"), 
                       legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "navy"), axis.ticks = element_blank(),
                       axis.text = element_text(color = "black")) 

ggplot(ScreePlots, aes(x = PC, y = VarianceFraction*100, color = DataType)) + geom_point(size = 5) + scatter_theme + scale_color_brewer(palette = "Set2") + expand_limits(y = 0)

### PC-based summary of proteomics data

# Scree
# PCA scores plot: PC1 ~ PC2
# Major PCs

proteomicsPC_plots <- list()

proteomics_svd <- svd(proteomics_reduced[,cond_reorder])
proteomics_PCs <-proteomics_svd$v[,1:6]
colnames(proteomics_PCs) <- paste0("PC", 1:6)
rownames(proteomics_PCs) <- colnames(proteomics_reduced)[cond_reorder]

scatterPlotTheme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(color = "black", fill = "white"), legend.position = "none",
                          axis.title.x = element_text(vjust = -0.3), axis.title.y = element_text(vjust = 0.25),
                          panel.grid.minor = element_blank(), legend.key.width = unit(6, "line"), panel.grid.major = element_line(colour = "black"), axis.ticks = element_line(colour = "black"), strip.background = element_rect(fill = "cyan"),
                          axis.text = element_text(color = "blacK")) 


ScreePlot <- data.frame(PC = 1:25, P = proteomics_svd$d^2 / sum(proteomics_svd$d^2))
ScreePlot$color <- ifelse(ScreePlot$PC <= npc, "RED", "BLACK")

proteomicsPC_plots$Scree <- ggplot(ScreePlot, aes(x = PC, y = P, color = color)) + geom_point(size = 5) + scatterPlotTheme + scale_color_identity() + expand_limits(y = c(0, max(ScreePlot$P)*1.04)) +
  scale_y_continuous("Fraction of variance explained", labels = percent, breaks = seq(0,0.25, by = 0.05), expand = c(0,0)) +
  scale_x_discrete("Principal component", breaks = c(1,5,10,15,20,25))

# PCA scores
proteomics_design <- data.frame(condition = rownames(proteomics_PCs), limitation = toupper(substr(rownames(proteomics_PCs), 1, 1)), DR = substr(rownames(proteomics_PCs), 2, 5))
proteomics_design$size <- sqrt(as.numeric(proteomics_design$DR)*1000)
proteomics_design$limitation <- factor(proteomics_design$limitation, levels = unique(proteomics_design$limitation))

proteomics_PC_aug <- data.frame(proteomics_PCs, proteomics_design)

proteomicsPC_plots$PCA_scores <- ggplot(proteomics_PC_aug, aes(x = PC2, y = PC1, color = limitation, size = size)) + geom_point(shape = 20, alpha = 0.7) + scale_size_identity() + 
  scale_color_brewer(palette = "Set2") + scatterPlotTheme +
  scale_x_continuous("Principal component 2") + scale_y_continuous("Principal component 1")

#library(gplots)
#heatmap.2(cor(proteomics_reduced[,cond_reorder]), trace = "none", Rowv = F, Colv = F)

# Proteomics PCs

scatterPlotTheme_facet <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(color = "black", fill = "cornsilk1"), legend.position = "none",
                          axis.title.x = element_text(vjust = -0.3), axis.title.y = element_text(vjust = 0.25), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), axis.ticks = element_blank(), strip.background = element_rect(fill = "burlywood1"),
                          axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 90)) 


proteomics_PC_melt <- melt(proteomics_PC_aug, id.vars = colnames(proteomics_design))
proteomics_PC_melt$condition <- factor(proteomics_PC_melt$condition, levels = unique(proteomics_PC_melt$condition))
proteomics_PC_melt$limitation <- factor(proteomics_PC_melt$limitation, levels = unique(proteomics_PC_melt$limitation))

proteomicsPC_plots$protPCs <- ggplot(proteomics_PC_melt, aes(x = condition, y = value, group = limitation, col = limitation)) + facet_grid(variable ~ ., scales = "free_y") + 
  geom_hline(yintercept = 0, size = 1) +
  geom_point(size = 4) + geom_line(size = 2) + scatterPlotTheme_facet + scale_color_brewer(palette = "Set2") +
  scale_size_identity()
  
proteomicsPC_plots$ncol <- 3

pdf(file = "Output/ProteomicsPCsummary.pdf", height = 7, width = 21)
do.call(grid.arrange,  proteomicsPC_plots)
dev.off()





#### Plot proteomics principal components ####

write.output.preclustered <- function(tab, output, gap = NA){
  #write an output table which can be read directly into Java treeview (without a dendrogram)
  
  if(all(!is.na(gap))){
    
  	for(i in 1:length(gap)){
  		if(i == 1){
  			tab_save <- cbind(tab[,1:gap[1]], gap = NA)
  			}else if(i == length(gap)){
  				tab_save <- cbind(tab_save, tab[,(gap[i-1]+1):gap[i]], gap = NA,tab[,(gap[length(gap)]+1):ncol(tab)])
  				}else{
  					tab_save <- cbind(tab_save, tab[,(gap[i-1]+1):gap[i]], gap = NA)
  					}
  	}
  	tab <- tab_save
  	}
  
  tab <- data.frame(GID = rownames(tab), Gene = rownames(tab), NAME = rownames(tab), GWEIGHT = 1, tab, stringsAsFactors = FALSE)
  Eweight <- c("Eweight", rep("", 3), rep(1, 25))
  tab <- rbind(Eweight, tab)
  
  write.table(tab, file = output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


npc <- 12
impSVD <- svd(impute_abund, nu = npc, nv = npc)
impSVD_pcs <- impSVD$v
colnames(impSVD_pcs) <- paste("PC", c(1:npc))
rownames(impSVD_pcs) <- toupper(colnames(impute_abund))
impSVD_pcs <- impSVD_pcs[c(16:20,1:5,11:15,6:10,21:25),]

write.output.preclustered(t(impSVD_pcs), "Output/proteomicPCs.cdt", gap = c(5,10,15,20))

#### plot proteomics / transcripts PCs ####
npc <- 12
impSVD <- svd(ptDiff, nu = npc, nv = npc)
impSVD_pcs <- impSVD$v
colnames(impSVD_pcs) <- paste("PC", c(1:npc))
rownames(impSVD_pcs) <- toupper(colnames(impute_abund))
impSVD_pcs <- impSVD_pcs[c(16:20,1:5,11:15,6:10,21:25),]

write.output.preclustered(t(impSVD_pcs), "Output/ptDiffPCs.cdt", gap = c(5,10,15,20))


######### Compare omics datasets in terms of how much variation is explained by experimental design ########

library(reshape2)
library(data.table)
library(ggplot2)
#library(eigenR2)
library(impute)
library(scales)

options(stringsAsFactor = F)

boer_metab <- read.delim('~/Desktop/Rabinowitz/Presentations/CommonImages/CombinedHeatmap/BoerMetabolites.txt')[-1,-2]
boer_metab <- data.table(melt(boer_metab, id.vars = "Metabolite", value.name = "RA"))
boer_metab$variable <- as.character(boer_metab$variable)
boer_metab[,Lim := substr(variable, 1, 3), by = "variable"]
boer_metab[,DR := substr(variable, 5, 1000000L), by = "variable"]
boer_metab <- boer_metab[,colnames(boer_metab) != "variable", with = F]
setnames(boer_metab, 'Metabolite', 'Feature')

chemo_DR <- read.delim("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/chemostatInformation/ListOfChemostats_augmented.txt")

protein_RA <- read.delim('~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/Protein_Transcript_Compare/Output/proteinRA.tsv')
rownames(protein_RA) <- protein_RA$Gene
protein_RA <- protein_RA[,-1]
protein_RA <- data.table(melt(data.frame(chemo_DR[chemo_DR$include_this, c('Limitation', 'actualDR')], t(protein_RA)), id.vars = c("Limitation", "actualDR")))
setnames(protein_RA, c('Limitation', 'actualDR', 'variable', 'value'), c('Lim', 'DR', 'Feature', 'RA'))

flux_rel <- read.delim('~/Desktop/Rabinowitz/Presentations/CommonImages/CombinedHeatmap/fluxCarried.tsv')
rownames(flux_rel) <- flux_rel$Gene
flux_rel <- flux_rel[,-1]
flux_rel <- data.table(melt(data.frame(chemo_DR[chemo_DR$include_this, c('Limitation', 'actualDR')], t(flux_rel)), id.vars = c("Limitation", "actualDR")))
setnames(flux_rel, c('Limitation', 'actualDR', 'variable', 'value'), c('Lim', 'DR', 'Feature', 'RA'))

pt_diff_rel <- read.delim('~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/Protein_Transcript_Compare/Output/coclust_ptdiff.tsv')
rownames(pt_diff_rel) <- pt_diff_rel$Gene
pt_diff_rel <- pt_diff_rel[,-1]
pt_diff_rel <- data.table(melt(data.frame(chemo_DR[chemo_DR$include_this, c('Limitation', 'actualDR')], t(pt_diff_rel)), id.vars = c("Limitation", "actualDR")))
setnames(pt_diff_rel, c('Limitation', 'actualDR', 'variable', 'value'), c('Lim', 'DR', 'Feature', 'RA'))

brauer_chemo_DR <- read.delim("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics/TranscriptomicsChemostatsActualDR.tsv") # import transcriptomics chemostat dilution rates

transcript_RA <- read.delim("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/Protein_Transcript_Compare/dilution_rate_00_raw.cdt")[-1,-c(1,3:5)]
rownames(transcript_RA) <- transcript_RA$YORF
transcript_RA <- transcript_RA[,-1]
colnames(transcript_RA) <- sub('([0-9.]{3})[0-9.]{2}$', '\\1', sub('limD', '0', colnames(transcript_RA)))
colnames(transcript_RA) <- sub('C', 'G', colnames(transcript_RA))
transcript_RA <- transcript_RA[,colnames(transcript_RA) %in% brauer_chemo_DR$SampleID]

if(!all(as.character(brauer_chemo_DR$SampleID) == colnames(transcript_RA))){
  stop("out of order")
  }

transcript_RA <- data.table(melt(data.frame(brauer_chemo_DR[, c('Nutrient', 'GR')], t(transcript_RA)), id.vars = c('Nutrient', 'GR')))
setnames(transcript_RA, c('Nutrient', 'GR', 'variable', 'value'), c('Lim', 'DR', 'Feature', 'RA'))

transcript_subset <- transcript_RA
transcript_subset <- transcript_subset[Feature %in% protein_RA$Feature,]

transcript_RA[,Experiment := 'Transcripts']
transcript_subset[,Experiment := 'Transcripts\n[proteomics subset]']
protein_RA[,Experiment := 'Proteins']
boer_metab[,Experiment := 'Metabolites'] 
flux_rel[,Experiment := 'Fluxes']
pt_diff_rel[,Experiment := 'Proteins\nper\ntranscript']

aggregated_omics <- rbind(transcript_RA, protein_RA, transcript_subset, pt_diff_rel, boer_metab, flux_rel, use.names = T)


eigenR2summary <- data.frame(Experiment = unique(aggregated_omics$Experiment), GR = NA, Lim = NA, GR_Lim = NA)
eigenR2_per_df <- data.frame(Experiment = unique(aggregated_omics$Experiment), GR = NA, Lim = NA, GR_Lim = NA)

for(exp in unique(aggregated_omics$Experiment)){
  
  omics_subset <- acast(aggregated_omics[Experiment == exp,], Feature ~ Lim + DR, value.var = "RA")
  omics_subset <- omics_subset[rowSums(!is.na(omics_subset)) >= ncol(omics_subset)/2,]
  omics_subset <- impute.knn(omics_subset)$data # this will only effect the brauer transcriptomics dataset
  
  # center each row
  omics_subset <- omics_subset - rowMeans(omics_subset)
  
  omics_covariates <- as.data.frame(t(sapply(colnames(omics_subset), function(x){strsplit(x, '_')[[1]]})))
  colnames(omics_covariates) <- c("Lim", "DR")
  omics_covariates$DR <- as.numeric(as.character(omics_covariates$DR))
  
  # calculate explained sum of squares explained by sets of covariates & explained mean square
  
  TSS <- sum(rowSums(omics_subset^2))
  
  # Test effect of GR as DR + 1 vs. 1
  
  DReffect <- model.matrix(~ omics_covariates$DR + 1) # 2 df
  
  DR_RSS <- sum(apply(omics_subset, 1, function(x){sum(lm(x ~ DReffect + 0)$resid^2)}))
  DR_Rsquared <- 1 - DR_RSS/TSS #(DR_RSS/(ncol(omics_subset)-ncol(DReffect)) / (TSS/(ncol(omics_subset)-1)))
  
  eigenR2summary$GR[eigenR2summary$Experiment == exp] <- DR_Rsquared
  eigenR2_per_df$GR[eigenR2summary$Experiment == exp] <- DR_Rsquared/(ncol(DReffect) - 1)
  
  
  # Test effect of Lim as two nested comparisons:
  # Lim + GR vs. GR + 1
  # Lim * GR vs. Lim + GR
  
  # add in limitation-specific intercept
  
  LimDReffect <- model.matrix(~ omics_covariates$Lim + omics_covariates$DR + 1) # 6 df
  
  LimDR_RSS <- sum(apply(omics_subset, 1, function(x){sum(lm(x ~ LimDReffect + 0)$resid^2)}))
  LimDR_Rsquared <- 1 - LimDR_RSS/TSS #(LimDR_RSS/(25-ncol(LimDReffect)) / (TSS/24))
  
  eigenR2summary$Lim[eigenR2summary$Experiment == exp] <- LimDR_Rsquared - DR_Rsquared
  eigenR2_per_df$Lim[eigenR2summary$Experiment == exp] <- (LimDR_Rsquared - DR_Rsquared)/(ncol(LimDReffect) - ncol(DReffect))
  
  # add in limitation-specific slope
  
  LimDRInteffect <- model.matrix(~ omics_covariates$Lim * omics_covariates$DR + 1) # 10 df
  
  LimDRInt_RSS <- sum(apply(omics_subset, 1, function(x){sum(lm(x ~ LimDRInteffect + 0)$resid^2)}))
  LimDRInt_Rsquared <- 1 - LimDRInt_RSS/TSS#(LimDRInt_RSS/(25-ncol(LimDRInteffect)) / (TSS/24))
  
  eigenR2summary$GR_Lim[eigenR2summary$Experiment == exp] <- LimDRInt_Rsquared - LimDR_Rsquared
  eigenR2_per_df$GR_Lim[eigenR2summary$Experiment == exp] <- (LimDRInt_Rsquared - LimDR_Rsquared)/(ncol(LimDRInteffect) - ncol(LimDReffect))
  
  }


eigen_melt <- melt(eigenR2summary)
eigen_melt <- eigen_melt[eigen_melt$Experiment != "Fluxes",]

eigen_melt$Experiment <- factor(eigen_melt$Experiment, levels = c("Transcripts", "Transcripts\n[proteomics subset]",
                                                                  "Proteins", "Proteins\nper\ntranscript", "Metabolites"))

eigen_melt <- data.table(eigen_melt)
eigen_melt[,scaled := value/sum(value), by = "Experiment"]

eigen_bar_label <- eigen_melt[eigen_melt$variable != "GR_Lim",]
eigen_bar_label$value[eigen_bar_label$variable == "Lim"] <- eigen_bar_label$value[eigen_bar_label$variable == "Lim"] + eigen_bar_label$value[eigen_bar_label$variable == "GR"]
eigen_bar_label$scaled[eigen_bar_label$variable == "Lim"] <- eigen_bar_label$scaled[eigen_bar_label$variable == "Lim"] + eigen_bar_label$scaled[eigen_bar_label$variable == "GR"]
eigen_bar_label$label = paste0(round(eigen_bar_label$scaled*100, 1), "%")

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), 
                       panel.background = element_rect(fill = "gray80"), legend.position = "top", 
                       axis.ticks.x = element_blank(), axis.ticks.y = element_line(color = "black"),
                       axis.text = element_text(color = "black"), axis.text.x = element_text(size = 20, angle = 60, hjust = 0.5, vjust = 0.5),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(size = 1),
                       axis.line = element_line(color = "black", size = 1), legend.title=element_blank(), legend.key = element_rect(color = "white")
                       )

ggplot(eigen_melt, aes(x = Experiment, y = value, fill = variable)) + 
  geom_bar(stat = "identity") +
  geom_text(data = eigen_bar_label, aes(label = label, y = value), vjust = -0.2) + 
  scale_fill_brewer(palette = "Set2", breaks = c("GR", "Lim", "GR_Lim"), labels = c("Growth-rate dependence", "Limitation dependence", "Limitation specific growth-rate dependence")) + expand_limits(y = c(0,1)) +
  scale_y_continuous("Variance Explained by\nExperimental Design", labels = percent_format(), expand = c(0,0)) +
  scale_x_discrete("") + barplot_theme + guides(fill = guide_legend(nrow = 2, byrow = T))

ggsave("Output/omicsVarianceExplained.pdf", width = 9, height = 7)



omics_comparison <- list()

# Summarize amount of variability in each data set
# bootstrapped quantiles of CV #

omics_CV <- aggregated_omics[, list(SD = sd(RA, na.rm = T)), by = c("Feature", "Experiment")]

omics_CV <- omics_CV[Experiment != "Fluxes",]
omics_CV$Experiment <- factor(omics_CV$Experiment, levels = c("Transcripts", "Transcripts\n[proteomics subset]",
                                                              "Proteins", "Proteins\nper\ntranscript", "Metabolites"))

omics_quantile_bs <- NULL

nbs <- 1000
for(exp in levels(omics_CV$Experiment)){
  
  exp_data <- omics_CV$SD[omics_CV$Experiment == exp]
  SD_bs_quantiles <- sapply(1:nbs, function(x){quantile(sample(exp_data, size = length(exp_data), replace = T), probs = c(0.5, 0.95))})
  
  quantile_CI <- apply(SD_bs_quantiles, 1, function(x){quantile(x, probs = c(0.025, 0.975))})
  rownames(quantile_CI) <- c("min", "max")
  m_quantile_CI <- melt(quantile_CI); colnames(m_quantile_CI) <- c("bound", "quantile", "value")
  
  omics_quantile_bs <- rbind(omics_quantile_bs, data.frame(Experiment = exp, m_quantile_CI))
  
  }

omics_quantile_bs <- dcast(omics_quantile_bs, "Experiment + quantile ~ bound", value.var = "value")

omics_quantile_bs$Experiment <- factor(omics_quantile_bs$Experiment, levels = c("Transcripts", "Transcripts\n[proteomics subset]",
                                                                                "Proteins", "Proteins\nper\ntranscript", "Metabolites"))
omics_comparison$omics_variability <- 
  ggplot() + geom_violin(data = omics_CV, aes(x = Experiment, y = SD), fill = "turquoise3") +
  #geom_errorbar(data = omics_quantile_bs, aes(x = Experiment, ymin = min, ymax = max, group = quantile), size = 1, color = "purple4") +
  scale_y_continuous(expression('Standard Deviation of' ~ log[2] ~ 'abundances') , expand = c(0,0)) +
  scale_x_discrete("") + barplot_theme + expand_limits(y = 0)


omics_varianceExplained <-  eigen_melt[, list(value = sum(value)), by = "Experiment"]

omics_comparison$variance_explained <- ggplot(omics_varianceExplained, aes(x = Experiment, y = value)) + 
  geom_bar(stat = "identity", fill = "chocolate1") +
  expand_limits(y = c(0,1)) +
  scale_y_continuous("Fraction of total variance explained\nby experimental design", labels = percent_format(), expand = c(0,0)) +
  scale_x_discrete("") + barplot_theme

omics_comparison$variance_partition <- ggplot(eigen_melt[variable == "GR",,], aes(x = Experiment, y = scaled, fill = variable)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2", breaks = c("GR", "Lim", "GR_Lim"), labels = c("Growth-rate dependence", "Limitation dependence", "Limitation specific growth-rate dependence"), guide = F) +
  scale_y_continuous("Fraction of explained variance\ndue to growth-rate", labels = percent_format(), expand = c(0,0)) +
  scale_x_discrete("") + barplot_theme


library(gridExtra)

pdf(file = "Output/omics_layout.pdf", height = 25, width = 9)
do.call(grid.arrange,  omics_comparison)
dev.off()

for(a_plot in names(omics_comparison)){
  
  ggsave(omics_comparison[[a_plot]], file = paste0("Output/", a_plot, ".pdf"), width = 9, height = 7)

}


#### PCs of multi-omics ####

omic_PCcomp <- NULL

for(exp in unique(aggregated_omics$Experiment)){
  
  omics_subset <- acast(aggregated_omics[Experiment == exp,], Feature ~ Lim + DR, value.var = "RA")
  omics_subset <- omics_subset[rowSums(!is.na(omics_subset)) >= ncol(omics_subset)/2,]
  omics_subset <- impute.knn(omics_subset)$data # this will only effect the brauer transcriptomics dataset
  
  # center each row
  omics_subset <- omics_subset - rowMeans(omics_subset)
  
  omics_covariates <- as.data.frame(t(sapply(colnames(omics_subset), function(x){strsplit(x, '_')[[1]]})))
  colnames(omics_covariates) <- c("Lim", "DR")
  omics_covariates$DR <- as.numeric(as.character(omics_covariates$DR))
  
  omic_svd <- as.data.frame(svd(omics_subset)$v[,1:3])
  colnames(omic_svd) <- c("PC1", "PC2", "PC3")
  
  omics_summary <- data.frame(omic_svd, omics_covariates)
  
  omics_summary$Lim <- factor(omics_summary$Lim, levels = c("P","C","N","L","U"))

  omic_PCcomp <- rbind(omic_PCcomp, data.frame(omic = exp, omics_summary))
  
  }

ggplot(data = omic_PCcomp, aes(x = PC1, y = PC2, color = Lim, size = DR)) + geom_point() + facet_wrap(~ omic)
  



