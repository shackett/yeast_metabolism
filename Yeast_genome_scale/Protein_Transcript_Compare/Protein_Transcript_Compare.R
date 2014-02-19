#### This script is a pared down version of Pep_Prot.Rnw focusing on direct comparisons of proteomic and transcriptional data ####

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/Protein_Transcript_Compare")
load("../../ChemicalSpeciesQuant/Proteomics/EMoutputDeg.Rdata")
load("../../ChemicalSpeciesQuant/Proteomics/EMimport.Rdata")
source("../../ChemicalSpeciesQuant/Proteomics/pep_library.R")

### Packages ####
library(impute)
library(missMDA)

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

transcript_brauer <- read.delim("../../ChemicalSpeciesQuant/brauer-microarray/Brauer_2008.pcl") # supplemental data in Brauer et al. 2008
rownames(transcript_brauer) <- transcript_brauer$SYSTEMATIC_NAME
transcript_brauer <- transcript_brauer[-1,-c(1:3)]

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
PTscatter_plotter
ggsave(plot = PTscatter_plotter, "Output/globalPTcomp.pdf", width = 10, height = 10)

table(PTscatterDF$Protein > 0, PTscatterDF$Transcript > 0)


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
ptDiff <- proteomics_reduced[,cond_reorder] - shared_trans_pool[,cond_reorder]
ptDiff <- ptDiff - apply(ptDiff, 1, mean)

### Save so that dendrogram of PTdiff can be generated and each matrix can be kmeans clustered so that shared motifs can be found

write.output(ptDiff, "Output/coclust_ptdiff.tsv")
write.output(proteomics_reduced, "Output/proteinRA.tsv")
write.output(transcript_brauer_reduced, "Output/transcriptRA.tsv")


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
reducedPTscatter[,GeneCommon := orf2common(as.character(Gene)), by = "Gene"]
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

ggplot(ScreePlots, aes(x = PC, y = VarianceFraction*100, color = DataType)) + geom_point(size = 5) + scatter_theme + scale_color_brewer(palette = "Set2")


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
