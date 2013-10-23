#source("http://bioconductor.org/biocLite.R")
#biocLite("sva")

library(gplots)
library(ggplot2)
library(reshape2)
library(data.table)
library(colorRamps)

source("../../Yeast_genome_scale/FBA_lib.R")

boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), legend.position = "top", 
  panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 12, angle = 90, hjust = 1), axis.line = element_blank(),
  axis.text = element_text(color = "black"))

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Lipids")
options(stringsAsFactors = F)

### Fatty acid quanitification experiment 1 ###
FAquant <- read.delim(file = "yeastFA_10_14_2013/FAgroups_10_14_13.csv", sep = ",")
FAheader <- read.delim(file = "yeastFA_10_14_2013/FAgroups_10_14_13.csv", sep = ",", header = F, nrows = 1)

FAmeta <- FAquant[,grep('[A-Z]{1,2}[0-9]{1,2}', FAheader, invert = T)]; colnames(FAmeta) <- FAheader[grep('[A-Z]{1,2}[0-9]{1,2}', FAheader, invert = T)]
FAmeta <- FAmeta[,c('medMz', 'medRt', 'compound')]
FAmeta$NC[grep('C[0-9]{1,2}', FAmeta$compound)] <- regmatches(FAmeta$compound, regexpr('C[0-9]{1,2}', FAmeta$compound))
FAmeta$NC <- sapply(FAmeta$NC, function(x){strsplit(x, split = 'C')[[1]][2]})
FAmeta$label <- ifelse(!is.na(FAmeta$NC), F, T)

# determine which compound matches unannotated labelled compounds

IDfa <- FAmeta[FAmeta$compound != "",]

for(i in c(1:nrow(FAmeta))[FAmeta$compound == ""]){
  
  mzDiff <- abs(IDfa$medMz + as.numeric(IDfa$NC) - FAmeta$medMz[i])
  rtDiff <- abs(FAmeta$medRt[i] - IDfa$medRt)
  
  FAmeta$compound[i] <- IDfa$compound[which.min(100 * mzDiff + rtDiff)]
  
  }
FAmeta$name <- sapply(1:nrow(FAmeta), function(x){
  paste(FAmeta$compound[x], ifelse(FAmeta$label[x], 'lab', 'unlab'), sep = "_")
  })
FAmeta$mgLipid <- NA
# standard is 1.6mg/mL - 20uL added
FAmeta$mgLipid[FAmeta$label] <- 1.6 * 0.02

### Sample specific information ###

wellinfo <- read.delim(file = "../BulkComposition/DRYWELL.txt")
wellinfo <- wellinfo[wellinfo$Type == "sample",]

wellinfo$Column <- as.character(wellinfo$Column)
wellinfo$Column <- sapply(wellinfo$Column, function(x){
  vsplit <- strsplit(x, split = "")[[1]]
  if(length(vsplit) == 1){paste("0", vsplit, sep = "")}else{paste(vsplit, collapse = "")}
  })

wellinfo$Well <- mapply(function(x,y){paste(x, y, sep = "")}, x = wellinfo$Row, y = wellinfo$Column)

FAmatrix <- FAquant[,grep('[A-Z]{1,2}[0-9]{1,2}', FAheader)]; colnames(FAmatrix) <- FAheader[grep('[A-Z]{1,2}[0-9]{1,2}', FAheader)]

FAheader <- data.frame(sample_name_raw = colnames(FAmatrix), sample_name = regmatches(colnames(FAmatrix), regexpr('[A-Z]{1,2}[0-9]{1,2}', colnames(FAmatrix))))

FAheader$sampleClass <- "sample"
FAheader$sampleClass[grep('M', FAheader$sample_name)] <- "mock"
FAheader$sampleClass[grep('BL', FAheader$sample_name)] <- "blank"

FAheader$DR = FAheader$Limitation = FAheader$mgDry= NA
FAheader[FAheader$sampleClass == "sample",c('Limitation', 'DR', 'mgDry')] <- wellinfo[chmatch(FAheader$sample_name[FAheader$sampleClass == "sample"], wellinfo$Well),c('Limitation', 'DR', 'mg_per_well')]

FAheader$condition <- FAheader$sampleClass
FAheader$condition[FAheader$sampleClass == "sample"] <- mapply(function(x, y){paste(x, y, sep = "")}, x = FAheader$Limitation[FAheader$sampleClass == "sample"], y = sprintf("%.2f", FAheader$DR[FAheader$sampleClass == "sample"]))
FAheader$condition <- factor(FAheader$condition, levels = c("blank", "mock", sort(unique(FAheader$condition[FAheader$sampleClass == "sample"]))))


### FA data ###

rownames(FAmatrix) <- FAmeta$name
FAmatrix[FAmatrix == 0] <- min(FAmatrix[FAmatrix != 0])
lFA <- log2(FAmatrix)
lFAcent <- as.matrix(lFA - apply(lFA, 1, median) %*% t(rep(1, ncol(lFA))))
heatmap.2(lFAcent[,FAheader$sampleClass == "sample"], trace = "none")  

augDF <- cbind(FAheader, t(lFA))

augDF <- augDF[FAheader$sampleClass != "blank",]

meltDF <- melt(augDF, id.vars = colnames(FAheader))

ggplot(meltDF, aes(x = condition, y = value, color = factor(Limitation))) + geom_boxplot() + facet_wrap(~ variable, ncol = 2, scales = "free_y") + boxplot_theme
ggsave("FAquant.pdf", height = 18, width = 18)



### Subtract off mock samples in absolute space ###

mockAbund <- apply(FAmatrix[,FAheader$sampleClass == "mock"], 1, mean) %*% t(rep(1, ncol(FAmatrix)))
mockAbund[FAmeta$label,] <- 0

FA_mock_subtracted <- FAmatrix - mockAbund
sample_header <- FAheader[FAheader$sampleClass == "sample",]
sample_abund <- FA_mock_subtracted[,FAheader$sampleClass == "sample"]
sample_abund[sample_abund < min(FAmatrix[FAmatrix != 0])] <- min(FAmatrix[FAmatrix != 0]) #floor small values to the minimum ascertained IC

### align unlabelled samples and their labeled standards to convert them into fraction DW ###

rownames(sample_abund) <- FAmeta$name

labeled_compounds <- FAmeta$compound[FAmeta$label]

labeled_unlabeled_index = chmatch(FAmeta$compound[FAmeta$compound %in% labeled_compounds & FAmeta$label == T], FAmeta$compound[FAmeta$compound %in% labeled_compounds & FAmeta$label == F])

sample_abund[FAmeta$compound %in% labeled_compounds & FAmeta$label == F,] <- sample_abund[FAmeta$compound %in% labeled_compounds & FAmeta$label == F,]/sample_abund[FAmeta$compound %in% labeled_compounds & FAmeta$label == T,][labeled_unlabeled_index,] * FAmeta$mgLipid[FAmeta$compound %in% labeled_compounds & FAmeta$label == T][labeled_unlabeled_index] %*% t(rep(1, ncol(sample_abund)))

sample_abund <- sample_abund / rep(1, nrow(sample_abund)) %*% t(sample_header$mgDry)

### remove labelled peaks and junk peaks ###

sample_abund <- sample_abund[FAmeta$label == F & !(FAmeta$compound %in% c("C22:1", "C26:1")),]
good_FAmeta <- FAmeta[FAmeta$label == F & !(FAmeta$compound %in% c("C22:1", "C26:1")),]
rownames(sample_abund) <- good_FAmeta$compound

### Absolute sample abundances corrected for 60% recovery during lipid extraction ###

sample_abund[good_FAmeta$compound %in% labeled_compounds,] <- sample_abund[good_FAmeta$compound %in% labeled_compounds,] * 10/6


### since FA with absolute quantification are a ratio, they should be analyzed in linear space ###

FAquant_append <- cbind(sample_header, t(sample_abund[good_FAmeta$compound %in% labeled_compounds,]))
FAquant_append$"All FA" <- apply(FAquant_append[,!colnames(FAquant_append) %in% colnames(sample_header)], 1, sum)

FAquant_melt <- data.table(melt(FAquant_append, id.vars = colnames(sample_header)))

FAquant_summary <- FAquant_melt[,list(mean = mean(value), se = sd(value)/sqrt(length(value))), by = c("condition", "variable", "Limitation")]
FAquant_summary[,CV := se/mean]

ggplot(FAquant_summary, aes(x = condition, y = mean, ymin = mean - 2*se, ymax = mean + 2*se, color = factor(Limitation))) + geom_pointrange(size = 1) + facet_wrap(~ variable, ncol = 1, scales = "free_y") + boxplot_theme +
  scale_color_brewer(palette = "Set2") + scale_y_continuous("% DW") + expand_limits(y = 0)
ggsave("FAdw.pdf", height = 10, width = 10)

write.table(FAquant_summary, quote = F, row.names = F, col.names = T, file = "FAabsoluteQuant.tsv", sep = "\t")

### combine all absolute peaks into a single pooled weight ###

sample_abund <- rbind(sample_abund, apply(sample_abund[good_FAmeta$compound %in% labeled_compounds,], 2, sum))
rownames(sample_abund)[nrow(sample_abund)] <- "All FA"


### look at the log2 data ###

l_sample_abund <- log2(sample_abund)
  
sample_info_append <- cbind(sample_header, t(l_sample_abund))
samplemeltDF <- melt(sample_info_append, id.vars = colnames(FAheader))


ggplot(samplemeltDF, aes(x = condition, y = value, fill = factor(Limitation))) + geom_boxplot() + facet_wrap(~ variable, ncol = 2, scales = "free_y") + boxplot_theme +
  scale_fill_brewer(palette = "Set1") + scale_y_continuous("log2 IC/mgDW or log2 mgFA/mgDW")
ggsave("FAquant2.pdf", height = 18, width = 18)




condSummary <- data.table(samplemeltDF)
condSummary <- condSummary[, list(mean = mean(value)), by = c("Limitation", "DR", "variable")]

condMatrix <- acast(condSummary, variable ~ Limitation + DR)

l_sample_cent <- condMatrix - apply(condMatrix, 1, mean) %*% t(rep(1, ncol(condMatrix)))
l_sample_cent <- as.matrix(l_sample_cent)

write.output(l_sample_cent, "FAabund.tsv")
  
heatmap.2(l_sample_cent, trace = "none", col = blue2yellow(100), Colv = F)






