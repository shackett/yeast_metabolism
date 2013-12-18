##### Combine experimental data on yeast composition, nutrient intake and excretion and expected composition into a set of condition-specific boundary fluxes #####

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(gplots)
library(reshape2)
library(data.table)

scatter_facet_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), 
      legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major = element_line(size = 0.5), axis.text.x = element_text(angle = 90),
      strip.background = element_rect(fill = "cadetblue2")) 

### Default composition function ####

compositionFile <- read.csv2("../Yeast_comp.csv", sep  = ",", stringsAsFactors = FALSE)

modelMetComp <- read.delim("../Yeast_genome_scale/flux_cache/stoiMetsComp.tsv", header = TRUE) # This file needs to be generated in FBA_run_full_reco as it includes additional metabolites which supplement the SBML model

compComp <- data.frame(t(sapply(compositionFile$AltName, function(x){
  incorporated <- sub('TP', 'MP', x)
  output <- unlist(modelMetComp[modelMetComp[,2] == incorporated,][1,])
  output[2] <- x
  output
  }))) #elemental composition each macromolecule in biomass function
  # for NTPs and dNTPs the relevent stoichiometry for determining the fractional contribution to dry-weight is (d)NMPs

compComp[compComp$name == "(1->3)-beta-D-glucan [cytoplasm]", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5) #one oxygen shared bc of condensation
compComp[compComp$name == "glycogen [cytoplasm]", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5)
compComp[compComp$name == "mannan [cytoplasm]", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5)
compComp[compComp$name == "polyphosphate [cytoplasm]", colnames(compComp) %in% c("O", "P")] <- c(4, 1)
compComp[compComp$name == "complex sphingolipid [cytoplasm]", colnames(compComp) %in% c("C", "H", "N", "O")] <- c(16, 35, 1, 2)  # this is primarily a begnign standin so the MW is not zero resulting in conversion problems

compComp[is.na(compComp)] <- 0
compComp[,-c(1:2)] <- apply(compComp[,-c(1:2)], c(1,2), as.numeric)
compComp <- compComp[,c(TRUE, TRUE, colSums(compComp[,-c(1:2)]) != 0)]

atomicMasses <- data.frame(element = c("C", "H", "N", "O", "P", "R", "S"), mass = c(12.0107, 1.00794, 14.00674, 15.9994, 30.973761, 0, 32.066))
atomicMasses <- atomicMasses[atomicMasses$element %in% colnames(compComp[,-c(1,2)]),]

compositionFile$MW <- t(t(compComp[,-c(1,2)])) %*% t(t(c(atomicMasses[,2])))
compositionFile$weightPerUn <- as.numeric(compositionFile$StoiCoef) * compositionFile$MW
#compositionFile$weight_per_t[compositionFile$Class == "Energy Balance"] <- NA

if(any(compositionFile$MW == 0)){
  print("incompatible MW")  
  }


### What is the dry-weight fraction of macromolecules in a conventional objective function ###

class_composition <- sapply(unique(compositionFile$Class), function(x){sum(compositionFile$weightPerUn[compositionFile$Class == x])}) * -1
class_composition <- data.frame(Category = names(class_composition[names(class_composition) != "Energy Balance"]), Abundance = unname(class_composition[names(class_composition) != "Energy Balance"]))
class_composition$CellComposition <- "Default"
class_composition <- class_composition[!is.na(class_composition$Abundance),]

pie_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), legend.position = "right", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(3, "line"), axis.text = element_blank(), axis.title = element_blank()) 


class_comp_plot <- ggplot(class_composition, aes(y = Abundance, factor(CellComposition), fill = factor(Category))) + pie_theme
class_comp_plot + geom_bar(stat = "identity") + coord_polar(theta = "y") + scale_fill_discrete(name = "Class")

class_comp_sep <- compositionFile[compositionFile$Class != "Energy Balance" & !is.na(compositionFile$weightPerUn),]
class_comp_sep$Abundance <- class_comp_sep$weightPerUn/sum(class_comp_sep$weightPerUn)
class_comp_sep$CellComposition <- "Default"
class_comp_sep <- class_comp_sep[order(class_comp_sep$Class),]
for(i in unique(class_comp_sep$Class)){
  class_comp_sep[class_comp_sep$Class == i,] <- class_comp_sep[c(1:length(class_comp_sep[,1]))[class_comp_sep$Class == i][order(class_comp_sep[class_comp_sep$Class == i, colnames(class_comp_sep) == "Abundance"], decreasing = TRUE)],]
  }

pie_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), legend.position = "right", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(3, "line"), axis.text = element_blank(), axis.title = element_blank(), line = element_line(size = 1)) 

class_comp_plot <- ggplot(class_comp_sep, aes(y = Abundance, factor(CellComposition), fill = factor(Class))) + pie_theme + scale_fill_brewer(name = "Class", palette = "Set1")
class_comp_plot + geom_bar(colour = "black", stat = "identity") + coord_polar(theta = "y")
ggsave(file = "default_composition.pdf", height = 15, width = 15)




########### Chemostat specific information ###########


chemostatInfo_bu <- read.table("BulkComposition/chemostatDRmisc.tsv", sep = "\t", header = TRUE) #VolFrac_mean - uL cellular vol per mL media

# Load all chemostat conditions
chemostatInfo <- read.table("chemostatInformation/ListOfChemostats_augmented.txt", sep = "\t", header = T)

sampleInfo <- read.table("chemostatInformation/Sample_information/sampleInfo.txt", sep = "\t", header = TRUE) ## dry weight from culture
sampleInfo <- sampleInfo[sampleInfo$Sample %in% chemostatInfo$composition_data,]
sampleInfo$DWperML <- sampleInfo$DryWeight/sampleInfo$cultureV #mg dry weight per mL of culture
sampleInfo$ChemostatID <- sampleInfo$Sample

sampleInfoMerge <- merge(chemostatInfo, sampleInfo, by = "ChemostatID")
# mg dry weight per uL of cells - this direct experimental ratio can be scaled the observed chemostat density to infer the expected dry weight density
sampleInfoMerge$mg_per_uL_cells <- sampleInfoMerge$DWperML/sampleInfoMerge$VolFrac_mean
   
# Switch to the subset which will be analyzed
chemostatInfo <- chemostatInfo[chemostatInfo$include_this,]
chemostatInfo$DWperML <- chemostatInfo$VolFrac_mean * sampleInfoMerge$mg_per_uL_cells[chmatch(chemostatInfo$composition_data, sampleInfoMerge$ChemostatID)]


# Initialize
n_c <- nrow(chemostatInfo)
comp_by_cond <- list()
tmp <- matrix(NA, ncol = n_c, nrow = length(compositionFile[,1]))
colnames(tmp) <- chemostatInfo$ChemostatCond; rownames(tmp) <- compositionFile$MetName
comp_by_cond$moles_per_cell <- comp_by_cond$grams_per_cell <- tmp

## coefficient of variation of experimental composition measurements ##
CV_table <- data.table(condition = chemostatInfo$ChemostatCond)
comp_by_cond$CV_table <- CV_table

chemostatInfo$DWperCell = chemostatInfo$DWperML * (1/1000) / chemostatInfo$VolFrac_mean * (chemostatInfo$medcellVol * 10^-9) #g dry weight per cell


##### For protein and RNA fraction, combine high-quality external experimental data with condition-specific #####
##### measurements (with some under-estimation) #######

# Reference
external_data <- read.delim("chemostatInformation/Sample_information/LiteratureComp.csv", sep = ",")

# Condition-specific protein fraction
proteinFrac <- read.table('BulkComposition/Protein_data.txt',sep='\t', header = T)
proteinFrac <- proteinFrac[chmatch(chemostatInfo$composition_data, proteinFrac$condition),]

protein_summary <- cbind(chemostatInfo[,c("ChemostatCond", "Limitation", "actualDR")], Specie = "Protein", Fraction = proteinFrac$fraction, Article = "Hackett")
colnames(protein_summary) <- colnames(external_data)
protein_summary$Fraction <- protein_summary$Fraction * 100

# Condition-specific RNA fraction
RNA_file <- data.table(read.delim("BulkComposition/RNA_abundance/RNAabund.csv", sep = ",", header = TRUE))
RNAFrac = RNA_file[,list(mean = mean(RNAconc/1000), SD = sd(RNAconc/1000, na.rm = T)/sqrt(length(RNAconc))), by = "condition"] #ug of RNA per mL of culture
RNAFrac$mg_per_mL_DRY <- chemostatInfo$DWperML[chmatch(RNAFrac$condition, chemostatInfo$composition_data)] #mg dry per mL culture
RNAFrac[,fraction := mean/mg_per_mL_DRY,]

RNA_summary <- cbind(chemostatInfo[,c("ChemostatCond", "Limitation", "actualDR")], Specie = "RNA", Fraction = RNAFrac$fraction[chmatch(chemostatInfo$composition_data, RNAFrac$condition)], Article = "Hackett")
colnames(RNA_summary) <- colnames(external_data)
RNA_summary$Fraction <- RNA_summary$Fraction * 100

all_prot_RNA_data <- rbind(external_data, protein_summary, RNA_summary)

## Use regression to find a consistent estimate of % composition treating Lange et al. data as best reference ##

fitted_species <- NULL

for(spec in unique(all_prot_RNA_data$Specie)){
  
  # Contrasts used are:
  # Limitation (5 levels) + DR (continuous) + Reference (2 - offsets relative to Lange)
  
  spec_subset <- all_prot_RNA_data[all_prot_RNA_data$Specie == spec,]
  
  limDesign <- cbind(ifelse(spec_subset$Limitation == "C", 1, 0), ifelse(spec_subset$Limitation == "N", 1, 0),
  ifelse(spec_subset$Limitation == "P", 1, 0), ifelse(spec_subset$Limitation == "L", 1, 0), ifelse(spec_subset$Limitation == "U", 1, 0))
  colnames(limDesign) <- c("C", "N", "P", "L", "U")
  
  studyDesign <- cbind(ifelse(spec_subset$Article == "Schulze 1995", 1, 0), ifelse(spec_subset$Article == "Hackett", 1, 0))
  colnames(studyDesign) <- c("Schulze offset", "Hackett offset")
  
  speciesDesignMat <- as.matrix(cbind(limDesign, DR = spec_subset$DR, studyDesign))
  
  speciesRegression <- lm(spec_subset$Fraction ~ speciesDesignMat + 0)
  
  # prediction of chemostat abundance, se
  
  speciesPredict <- predict.lm(speciesRegression, newdata = as.data.frame(speciesDesignMat), se.fit = T)
  speciesPredict$fit[spec_subset$Article == "Schulze 1995"] <- speciesPredict$fit[spec_subset$Article == "Schulze 1995"] - speciesRegression$coef[names(speciesRegression$coef) == "speciesDesignMatSchulze offset"]
  speciesPredict$fit[spec_subset$Article == "Hackett"] <- speciesPredict$fit[spec_subset$Article == "Hackett"] - speciesRegression$coef[names(speciesRegression$coef) == "speciesDesignMatHackett offset"]
  
  spec_subset$Fitted_Fraction <- speciesPredict$fit
  spec_subset$Fitted_SE <- speciesPredict$se
  
  fitted_species <- rbind(fitted_species, spec_subset)
  
  
}

ggplot(fitted_species, aes(x = DR)) + facet_grid(Specie ~ Limitation, scale = "free_y") + 
  geom_errorbar(aes(ymax = Fitted_Fraction + 2*Fitted_SE, ymin = Fitted_Fraction - 2*Fitted_SE, color = factor(Article)), size = 1) +
  geom_point(aes(y = Fraction), size = 4, color = "black", shape = 16 ) +
  geom_point(aes(y = Fraction, color = factor(Article)), size = 3, shape = 16 ) + expand_limits(y = 0) + scatter_facet_theme +
  ggtitle('RNA and Protein Dry-weight data \n from 3 sources and concensus trends')

ggsave("LiteratureConsensusComp.pdf", height = 12, width = 12)


## Consensus composition measurements added to condition-specific composition ##
## Protein ##

protein_subset <- fitted_species[fitted_species$Specie == "Protein" & fitted_species$Article == "Hackett",]
protein_subset <- protein_subset[chmatch(chemostatInfo$ChemostatCond, protein_subset$Condition),]
protein_subset$DWperCell <- chemostatInfo$DWperCell

if(!(all(protein_subset$Condition == chemostatInfo$condition))){
 print("Species mismatch") 
}

CV_table[,"AA flux" := protein_subset$Fitted_SE/protein_subset$Fitted_Fraction,]

aaRelAbunds <- compositionFile[compositionFile$Class == "Amino Acid",]
aaWeightFrac <- (aaRelAbunds$weightPerUn*-1)/sum(aaRelAbunds$weightPerUn*-1)

comp_by_cond$moles_per_cell[compositionFile$Class == "Amino Acid",] <- t(t(t((protein_subset$Fitted_Fraction/100 * chemostatInfo$DWperCell))) %*% t(aaWeightFrac/aaRelAbunds$MW))

## RNA ##

RNA_subset <- fitted_species[fitted_species$Specie == "RNA" & fitted_species$Article == "Hackett",]
chemostatInfo_subset <- chemostatInfo[chemostatInfo$ChemostatCond %in% RNA_subset$Condition,]
RNA_subset <- RNA_subset[chmatch(chemostatInfo_subset$ChemostatCond, RNA_subset$Condition),]
RNA_subset$DWperCell <- chemostatInfo_subset$DWperCell

if(!(all(RNA_subset$Condition == chemostatInfo_subset$condition))){
 print("Species mismatch") 
}

CV_table[chmatch(RNA_subset$Condition, CV_table$condition), "NTP flux" := RNA_subset$Fitted_SE/RNA_subset$Fitted_Fraction]

rnaRelAbunds <- compositionFile[compositionFile$MetName %in% c("ATP", "UTP", "CTP", "GTP"),]
rnaWeightFrac <- (rnaRelAbunds$weightPerUn*-1)/sum(rnaRelAbunds$weightPerUn*-1)

comp_by_cond$moles_per_cell[compositionFile$MetName %in% c("ATP", "UTP", "CTP", "GTP"), chmatch(RNA_subset$Condition, CV_table$condition)] <-
  t(t(t((RNA_subset$Fitted_Fraction/100 * chemostatInfo_subset$DWperCell))) %*% t(rnaWeightFrac/rnaRelAbunds$MW))

##### Stand-alone abundance of fatty-acids, carbohydrate and glycerol ######
## expressed as fraction of dry weight


## carbohydrate

carbFrac <- read.table('BulkComposition/TCfrac_data.txt',sep='\t', header = T)
carbFrac <- carbFrac[chmatch(chemostatInfo$composition_data, carbFrac$conditions),]
CV_table[,"sugar polymer flux" := sqrt(carbFrac$assayVariance)/carbFrac$fraction,]

carbRelAbunds <- compositionFile[compositionFile$Class %in% c("Cell Wall Carbohydrates", "Storage Carbohydrates"),]
carbWeightFrac <- (carbRelAbunds$weightPerUn*-1)/sum(carbRelAbunds$weightPerUn*-1)

comp_by_cond$moles_per_cell[compositionFile$Class %in% c("Cell Wall Carbohydrates", "Storage Carbohydrates"),] <- t(t(t((carbFrac$fraction * chemostatInfo$DWperCell))) %*% t(carbWeightFrac/carbRelAbunds$MW))

## glycerol

glycerolFrac <- read.table('BulkComposition/glycfrac_data.txt',sep='\t', header = T)
glycerolFrac <- glycerolFrac[chmatch(chemostatInfo$composition_data, glycerolFrac$condition),]
CV_table[,"glycerol washout" := sqrt(glycerolFrac$assayVariance)/glycerolFrac$fraction,]

comp_by_cond$moles_per_cell[compositionFile$MetName == "glycerol",] <- t(t(t((glycerolFrac$fraction * chemostatInfo$DWperCell))) %*% t(1/compositionFile$MW[compositionFile$MetName == "glycerol"]))

## polyphosphate production ##

polyphosphateFrac <- read.table('BulkComposition/PPfrac_data.txt',sep='\t', header = T)
polyphosphateFrac <- polyphosphateFrac[chmatch(chemostatInfo$composition_data, polyphosphateFrac$condition),]
CV_table[,"polyphosphate flux" := sqrt(polyphosphateFrac$assayVariance)/polyphosphateFrac$fraction,]

comp_by_cond$moles_per_cell[compositionFile$MetName == "polyphosphate",] <- t(t(t((polyphosphateFrac$fraction * chemostatInfo$DWperCell))) %*% t(1/compositionFile$MW[compositionFile$MetName == "polyphosphate"]))


## Total DNA ##
## genome size + genome size * fraction of cells not in G1 # from Brauer

buddingFrac = 0.936 - 1.971*chemostatInfo$actualDR #brauer 2008 relationship between unbudded fraction and growth rate
genomeLength = 12157105 # yeast genome length from SGD http://www.yeastgenome.org/cache/genomeSnapshot.html
GContent = 0.383 #http://bionumbers.hms.harvard.edu//bionumber.aspx?id=102126&ver=0
avogadros = 6.02214e23
 
comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) %in% c("dATP", "dTTP"),] <- rbind((genomeLength * (1-GContent) / avogadros)*(1+(1-buddingFrac)), (genomeLength * (1-GContent) / avogadros)*(1+(1-buddingFrac)))
comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) %in% c("dGTP", "dCTP"),] <- rbind((genomeLength * GContent / avogadros)*(1+(1-buddingFrac)), (genomeLength * GContent / avogadros)*(1+(1-buddingFrac)))


## Fatty acid content ##
## Absolute quantification for each of 4 fatty acids C16:0,1 & C18:0,1 ##

FAFrac <- read.table('Lipids/FAabsoluteQuant.tsv', sep = "\t", header = T)

for(FA in compositionFile[compositionFile$Class == "Fatty Acids",]$MetName){
  
  a_FA <- FAFrac[FAFrac$variable == FA,]
  a_FA <- a_FA[chmatch(chemostatInfo$composition_data, a_FA$condition),]
  
  comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) == FA,] <- t(t(t((a_FA$mean * chemostatInfo$DWperCell))) %*% t(1/compositionFile$MW[compositionFile$MetName == FA]))

  CV_table[,eval(FA) := a_FA$CV]
  
  }



##### Combining observed abundances with total dry weight per cell and inferring contributions of non-measured elements based upon the assumption that they remain in constant proportions ####
# scaling energy usage by cell weight relative to default weight.

comp_by_cond$grams_per_cell <- comp_by_cond$moles_per_cell * (t(t(compositionFile$MW)) %*% t(rep(1,n_c)))
cond_dryweight <- chemostatInfo$DWperCell * 10^12#pg per cell

weightSummary <- comp_by_cond$grams_per_cell[rowSums(!is.na(comp_by_cond$grams_per_cell)) != 0,]; weightSummary[is.nan(weightSummary)] <- NA
weightSummary <- weightSummary * 10^12 #convert from moles to pmoles
summaryInfo <- compositionFile[rowSums(!is.na(comp_by_cond$grams_per_cell)) != 0,]

weightSummary <- rbind(weightSummary, sapply(cond_dryweight - colSums(weightSummary, na.rm = TRUE), function(x){max(x, 0)})) # determine residual dry weight by subtracting measured components
rownames(weightSummary)[length(weightSummary[,1])] <- "Residual"
summaryInfo <- summaryInfo[,colnames(summaryInfo) %in% c("MetName", "Class", "Abbreviated")]
summaryInfo <- rbind(summaryInfo, c("Residual", "Residual dry weight", "Residual"))

summaryInfo$Class[grep('Carbohydrates', summaryInfo$Class)] <- "Carbohydrates"

weightAndSum <- cbind(summaryInfo, weightSummary)
weightAndSum <- weightAndSum[order(weightAndSum$Class),] # arrange by metabolite class

# arrange within metabolite class by abundance
for(i in unique(weightAndSum$Class)){
  weightAndSum[weightAndSum$Class == i,] <- weightAndSum[c(1:length(weightAndSum[,1]))[weightAndSum$Class == i][order(weightAndSum$P0.05[weightAndSum$Class == i], decreasing = TRUE)],]
  }


weightStackDF <- melt(weightAndSum)
colnames(weightStackDF)[4:5] <- c("Condition", "Abundance")


weightStackDF$fraction <- sapply(1:length(weightStackDF[,1]), function(x){
  weightStackDF$Abundance[x]/sum(weightStackDF$Abundance[weightStackDF$Condition == weightStackDF$Condition[x]], na.rm = TRUE)
  })

weightStackDF$Class <- factor(weightStackDF$Class)

#### Confidence intervals for pg per cell ###

weightStackDF <- weightStackDF[grep('p0.05H', weightStackDF$Condition, invert = T),]
weightStackDF$Limitation <- factor(toupper(sapply(as.character(weightStackDF$Condition), function(x){strsplit(x, split = '')[[1]][1]})), levels = c("P", "C", "N", "L", "U"))
weightStackDF$DR <- sapply(as.character(weightStackDF$Condition), function(x){paste(strsplit(x, split = '')[[1]][-1], collapse = "")})

class_sum_wpc <- acast(weightStackDF, formula = Class ~ Condition, sum, value.var = "Abundance")
cumsum_height <- apply(class_sum_wpc, 2, cumsum)
class_sum_melt <- data.frame(melt(class_sum_wpc), melt(cumsum_height)[,3])
colnames(class_sum_melt) <- c("Class", "Condition", "Abundance", "cumsum")
class_sum_melt$SE <- NA

class_sum_melt$SE[class_sum_melt$Class == "Amino Acid"] <- CV_table$"AA flux"[chmatch(as.character(class_sum_melt$Condition[class_sum_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_melt$SE[class_sum_melt$Class == "Carbohydrates"] <- CV_table$"sugar polymer flux"[chmatch(as.character(class_sum_melt$Condition[class_sum_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_melt$SE[class_sum_melt$Class == "Glycerol"] <- CV_table$"glycerol washout"[chmatch(as.character(class_sum_melt$Condition[class_sum_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_melt$SE[class_sum_melt$Class == "Nucleic Acid"] <- CV_table$"NTP flux"[chmatch(as.character(class_sum_melt$Condition[class_sum_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_melt$SE[class_sum_melt$Class == "Polyphosphates"] <- CV_table$"polyphosphate flux"[chmatch(as.character(class_sum_melt$Condition[class_sum_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_melt$SE[class_sum_melt$Class == "Fatty Acids"] <- apply(CV_table[,compositionFile[compositionFile$Class == "Fatty Acids",]$MetName,with = F], 1, mean)[chmatch(as.character(class_sum_melt$Condition[class_sum_melt$Class == "Fatty Acids"]), CV_table$condition)]

class_sum_melt <- class_sum_melt[!is.na(class_sum_melt$SE),]
class_sum_melt$lb <- class_sum_melt$cumsum - class_sum_melt$SE*class_sum_melt$Abundance
class_sum_melt$ub <- class_sum_melt$cumsum + class_sum_melt$SE*class_sum_melt$Abundance
class_sum_melt$Limitation <- factor(toupper(sapply(as.character(class_sum_melt$Condition), function(x){strsplit(x, split = '')[[1]][1]})), levels = c("P", "C", "N", "L", "U"))
class_sum_melt$DR <- sapply(as.character(class_sum_melt$Condition), function(x){paste(strsplit(x, split = '')[[1]][-1], collapse = "")})


#### Confidence intervals for % comp ####

class_sum_fraction <- acast(weightStackDF, formula = Class ~ Condition, sum, value.var = "fraction")
cumsum_height <- apply(class_sum_fraction, 2, cumsum)
class_sum_fraction_melt <- data.frame(melt(class_sum_fraction), melt(cumsum_height)[,3])
colnames(class_sum_fraction_melt) <- c("Class", "Condition", "Abundance", "cumsum")
class_sum_fraction_melt$SE <- NA

class_sum_fraction_melt$SE[class_sum_fraction_melt$Class == "Amino Acid"] <- CV_table$"AA flux"[chmatch(as.character(class_sum_fraction_melt$Condition[class_sum_fraction_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_fraction_melt$SE[class_sum_fraction_melt$Class == "Carbohydrates"] <- CV_table$"sugar polymer flux"[chmatch(as.character(class_sum_fraction_melt$Condition[class_sum_fraction_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_fraction_melt$SE[class_sum_fraction_melt$Class == "Glycerol"] <- CV_table$"glycerol washout"[chmatch(as.character(class_sum_fraction_melt$Condition[class_sum_fraction_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_fraction_melt$SE[class_sum_fraction_melt$Class == "Nucleic Acid"] <- CV_table$"NTP flux"[chmatch(as.character(class_sum_fraction_melt$Condition[class_sum_fraction_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_fraction_melt$SE[class_sum_fraction_melt$Class == "Polyphosphates"] <- CV_table$"polyphosphate flux"[chmatch(as.character(class_sum_fraction_melt$Condition[class_sum_fraction_melt$Class == "Amino Acid"]), CV_table$condition)]
class_sum_fraction_melt$SE[class_sum_fraction_melt$Class == "Fatty Acids"] <- apply(CV_table[,compositionFile[compositionFile$Class == "Fatty Acids",]$MetName,with = F], 1, mean)[chmatch(as.character(class_sum_fraction_melt$Condition[class_sum_fraction_melt$Class == "Fatty Acids"]), CV_table$condition)]

class_sum_fraction_melt <- class_sum_fraction_melt[!is.na(class_sum_fraction_melt$SE),]
class_sum_fraction_melt$lb <- class_sum_fraction_melt$cumsum - class_sum_fraction_melt$SE*class_sum_fraction_melt$Abundance
class_sum_fraction_melt$ub <- class_sum_fraction_melt$cumsum + class_sum_fraction_melt$SE*class_sum_fraction_melt$Abundance
class_sum_fraction_melt$lbnonstack <- class_sum_fraction_melt$Abundance - class_sum_fraction_melt$SE*class_sum_fraction_melt$Abundance
class_sum_fraction_melt$ubnonstack <- class_sum_fraction_melt$Abundance + class_sum_fraction_melt$SE*class_sum_fraction_melt$Abundance

class_sum_fraction_melt$Limitation <- factor(toupper(sapply(as.character(class_sum_fraction_melt$Condition), function(x){strsplit(x, split = '')[[1]][1]})), levels = c("P", "C", "N", "L", "U"))
class_sum_fraction_melt$DR <- sapply(as.character(class_sum_fraction_melt$Condition), function(x){paste(strsplit(x, split = '')[[1]][-1], collapse = "")})

#plot condition specific composition from experimental measurements

barplot_theme <- theme(text = element_text(size = 40, face = "bold"), title = element_text(size = 40, face = "bold"), panel.background = element_blank(), legend.position = "top", 
                       panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text = element_text(size = 30, color = "black"), 
                       axis.text.x = element_text(size = 20, angle = 70, vjust = 0.5, colour = "BLACK"), legend.key.width = unit(3, "line"),
                       strip.background = element_rect(fill = "darkgoldenrod1"), panel.margin = unit(2, "lines")) 

cbPalette <- c("#56B4E9", "#D55E00", "#E69F00", "#009E73", "#F0E442", "#999999", "#CC79A7")

ggplot() + geom_bar(data = weightStackDF, aes(y = Abundance, x = DR, fill = factor(Class)), colour = "black", stat = "identity") + scale_y_continuous("pg per cell +/- se", expand = c(0,0)) + scale_x_discrete("Dilution Rate") + scale_fill_manual(name = "Class", values = cbPalette, guide = guide_legend(nrow = 2)) +
  geom_errorbar(data = class_sum_melt, aes(x = DR, ymin = lb, ymax = ub), color = "WHITE") + barplot_theme + facet_grid(~ Limitation) + ggtitle("Total macromolecules and composition greatly varies across growth conditions")
ggsave(file = "condition_composition.pdf", height = 12, width = 20)


ggplot() + geom_bar(data = weightStackDF, aes(y = fraction, x = DR, fill = factor(Class)), stat = "identity") + scale_y_continuous("Fraction of dry weight +/- se", expand = c(0,0)) + scale_x_discrete("Dilution Rate") + scale_fill_manual(name = "Class", values = cbPalette, guide = guide_legend(nrow = 2)) +
  geom_errorbar(data = class_sum_fraction_melt, aes(x = DR, ymin = lbnonstack, ymax = ubnonstack), color = "black") + barplot_theme + facet_grid(Class ~ Limitation, scale = "free") + ggtitle("Relative abundance of macromolecules varies across growth conditions")
ggsave(file = "facet_composition.pdf", height = 32, width = 20)


ggplot() + geom_bar(data = weightStackDF[!is.na(weightStackDF$Abundance),], aes(y = fraction, x = DR, fill = factor(Class), position = "fill"), stat = "identity") + scale_y_continuous("Fraction of dry weight +/- se", expand = c(0,0)) + scale_x_discrete("Dilution Rate") + scale_fill_manual(name = "Class", values = cbPalette, guide = guide_legend(nrow = 2)) +
  geom_errorbar(data = class_sum_fraction_melt, aes(x = DR, ymin = lb, ymax = ub), color = "white") + barplot_theme + facet_grid( ~ Limitation, scale = "free") + ggtitle("Fractional composition of macromolecules across growth conditions")
ggsave(file = "fractional_composition.pdf", height = 14, width = 25)




#fill in missing components relative to measured species

comp_by_cond$grams_per_cell[rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c,] <- t(t(compositionFile$weightPerUn[rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c]/sum(compositionFile$weightPerUn, na.rm = TRUE))) %*% colSums(comp_by_cond$grams_per_cell, na.rm = TRUE)

comp_by_cond$moles_per_cell[compositionFile$Class != "Energy Balance",] <- comp_by_cond$grams_per_cell[compositionFile$Class != "Energy Balance",]/t(t(compositionFile$MW[compositionFile$Class != "Energy Balance"])) %*% rep(1, n_c)

# calculate the maintenance flux from dry weight and polymerization costs from molarity

# maintenance flux - relative to DW
maintFlux <- 0.001 * colSums(comp_by_cond$grams_per_cell) #1 mmole per gDW/hr
maintVec <- c(1, -1, -1); names(maintVec) <- c("ATP [cytoplasm]", "ADP [cytoplasm]", "phosphate [cytoplasm]")

maintMatrix <- maintVec %*% t(maintFlux/chemostatInfo$actualDR) #divide out the dilution rate so that when culture concentration is multiplied by DR in FBA_run_full_reco.R, this adjustment is canceled leaving flux per hr
rownames(maintMatrix) <- names(maintVec)

# energetic flux in order to polymerization monomers
composition_part <- NULL
for(i in 1:length(compositionFile[,1])){
  
  poly_stoi <- compositionFile$Polymerization_Cost[i]
  poly_stoi <- strsplit(poly_stoi, '->')[[1]]
  reactants <- strsplit(poly_stoi[1], '\\+')[[1]]
  products <- strsplit(poly_stoi[2], '\\+')[[1]]
  for(reactant in reactants){
    reactant <- sub('^ ', '', reactant)
    reactant <- sub(' $', '', reactant)
    
    reactantSpec <- ifelse(length(grep('[0-9]\\.*[0-9]*', reactant)) == 0, reactant, regmatches(reactant, regexpr('[0-9]\\.*[0-9]*', reactant), invert = TRUE)[[1]][2])
    reactantSpec <- sub('^ ', '', reactantSpec)
    reactantStoi <- ifelse(length(grep('[0-9]\\.*[0-9]*', reactant)) == 0, 1, regmatches(reactant, regexpr('[0-9]\\.*[0-9]*', reactant)))
    
    composition_part <- rbind(composition_part, data.frame(name = compositionFile$AltName[i], compound = reactantSpec, stoi = as.numeric(reactantStoi)))
    }
  for(product in products){
    product <- sub('^ ', '', product)
    product <- sub(' $', '', product)
    
    productSpec <- ifelse(length(grep('[0-9]\\.*[0-9]*', product)) == 0, product, regmatches(product, regexpr('[0-9]\\.*[0-9]*', product), invert = TRUE)[[1]][2])
    productSpec <- sub('^ ', '', productSpec)
    productStoi <- ifelse(length(grep('[0-9]\\.*[0-9]*', product)) == 0, 1, regmatches(product, regexpr('[0-9]\\.*[0-9]*', product)))
    
    composition_part <- rbind(composition_part, data.frame(name = compositionFile$AltName[i], compound = productSpec, stoi = -1*as.numeric(productStoi)))
    }
  }

composition_part <- dcast(composition_part, formula = name ~ compound, sum, value.var = "stoi")
composition_part <- composition_part[,colnames(composition_part) != "NA"]

composition_part$name <- factor(composition_part$name, levels = compositionFile$AltName)
composition_part <- composition_part[order(composition_part$name),]
composition_part$name <- as.character(composition_part$name)

### integrate these costs of polymerization directly into FBA rather than smooshing them together
### save composition_part to comp_by_cond before output


comp_by_cond$moles_per_cell <- rbind(comp_by_cond$moles_per_cell, maintMatrix)

### rewrite compositionFile so that added species in polymerization cost have appropriate names in FBA model
added_cmpds <- data.frame(MetName = rownames(maintMatrix), StoiCoef = NA, AltName = rownames(maintMatrix), Class = "Maintenance ATP hydrolysis", Abbreviated = rownames(maintMatrix), Polymerization_Cost = NA, MW = NA, weightPerUn = NA)
compositionFile_Extended = data.frame(t(cbind(t(compositionFile), t(added_cmpds))))
#write.table(compositionFile_Extended, file = "../Yeast_comp_energy.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





# generate weight per cellular volume
colSums(comp_by_cond$grams_per_cell, na.rm = TRUE) / (chemostatInfo$medcellVol * 10^-15) / 10  #grams per 100mL of macromolecules - chunky but reasonable

# determine the fluxes into macromolecule biosynthesis and energy per volume
comp_by_cond$intacellularMolarity = (comp_by_cond$moles_per_cell)/(t(t(rep(1, length(compositionFile_Extended[,1])))) %*% chemostatInfo$medcellVol * 10^-15)



comp_by_cond$anabolicFlux = comp_by_cond$intacellularMolarity * (t(t(rep(1, length(compositionFile_Extended[,1])))) %*% chemostatInfo$actualDR) #moles/L-hr 

comp_by_cond$cultureMolarity <- comp_by_cond$intacellularMolarity * t(t(rep(1, length(compositionFile_Extended[,1])))) %*% chemostatInfo$VolFrac_mean/1000 #moles/L culture








########### Boundary fluxes from nutrient uptake or metabolite excretion ############

SPM <- read.delim("Media/measuredNutrients.tsv")
NMR <- read.delim("Media/mediaComposition_NMR.tsv")
nutrientFile <- read.delim("../Yeast_genome_scale/Boer_nutrients.txt")[1:6,1:6] #molar concentration of nutrients
nutrientFileConv <- data.frame(mediaName = colnames(nutrientFile)[-1], shortName = c("N", "P", "C", "L", "U"))

mediaSummary <- NULL

### Nutrient not present here will be solely bounded by mediaConcentration*DR

### phosphate uptake
SPM <- SPM[chmatch(chemostatInfo$composition_data, SPM$condition),]
SPM$SDconcensus <- median(SPM$phosphateSD/(SPM$ceiling - SPM$phosphateUptake))*(SPM$ceiling - SPM$phosphateUptake) #use a fixed coefficient of variation 13%

SPM$maxSD <- mapply(function(x,y){max(x, y)}, x = SPM$phosphateSD, y = SPM$SDconcensus)


# take a consensus SD for p-lim and !p-lim

SPM_out <- data.frame(condition = chemostatInfo$ChemostatCond, specie = "phosphate [extracellular]", change = SPM$phosphateUptake, sd = SPM$maxSD, lb = 0, ub = SPM$ceiling, type = "uptake")
mediaSummary <- rbind(mediaSummary, SPM_out)

### NMR data - glucose -> EtOH, Ac, glycerol 
speciesOfInterest <- data.frame(NMRname = c("Glucose", "Ethanol", "Acetate", "Lactate", "Glycerol"), mediaName = c("D-glucose", NA, NA, NA, NA), 
                                modelName = c("D-glucose [extracellular]", "ethanol [extracellular]", "acetate [extracellular]", "(R)-lactate [extracellular]", "glycerol [extracellular]"), 
                                type = c("uptake", rep("excretion", 4)))

for(a_specie_n in 1:nrow(speciesOfInterest)){
  
  # match NMR data by associated NMR_name and NMR_date
  data_index <- mapply(function(x,y){c(1:nrow(NMR))[NMR$peak == speciesOfInterest$NMRname[a_specie_n] & NMR$condition == x & NMR$date == y]
  }, x = chemostatInfo$NMR_name, y = chemostatInfo$NMR_date)
  
  NMR_data_subset <- NMR[data_index,]
  NMR_data_subset$condition <- factor(NMR_data_subset$condition, levels = chemostatInfo$NMR_name)
  NMR_data_subset <- NMR_data_subset[order(NMR_data_subset$condition),]
  NMR_data_subset$consensusSD <- NMR_data_subset$estimate * median(NMR_data_subset$se/NMR_data_subset$estimate, na.rm = T) #take a consensus SD using a fixed CV
  NMR_data_subset$maxSD <-  mapply(function(x,y){max(x, y, na.rm = T)}, x = NMR_data_subset$se, y = NMR_data_subset$consensusSD)
  
  if(speciesOfInterest$type[a_specie_n] == "uptake"){
    
    mediaComp <- nutrientFile[nutrientFile$X == speciesOfInterest$mediaName[a_specie_n],-1]
    upper_bound <- unname(unlist(sapply(as.character(NMR_data_subset$condition), function(x){mediaComp[names(mediaComp) == nutrientFileConv$mediaName[nutrientFileConv$shortName == strsplit(x, '')[[1]][1]]]})))
    NMR_out <- data.frame(condition = NMR_data_subset$condition, specie = speciesOfInterest$modelName[a_specie_n], change = (upper_bound - NMR_data_subset$estimate/1000), sd = NMR_data_subset$maxSD/1000, lb = 0, ub = upper_bound, type = "uptake")
    
  }else{
    NMR_out <- data.frame(condition = NMR_data_subset$condition, specie = speciesOfInterest$modelName[a_specie_n], change = NMR_data_subset$estimate/1000, sd = NMR_data_subset$maxSD/1000, lb = 0, ub = Inf, type = "excretion")
    }
  
  mediaSummary <- rbind(mediaSummary, NMR_out)
  }
mediaSummary <- data.table(mediaSummary)

boundary_ele_comp <- mediaSummary[,modelMetComp[modelMetComp$name == specie,][1,],by = specie]
boundary_ele_comp[is.na(boundary_ele_comp)] <- 0
boundary_ele_comp <- boundary_ele_comp[,-grep('Fe|K|Na', colnames(boundary_ele_comp)),with = F]

boundary_ele_comp$MW <- rowSums(boundary_ele_comp[,atomicMasses$element,with = F] * rep(1, nrow(boundary_ele_comp)) %*% t(atomicMasses$mass))

mediaSummary$density <- sapply(1:nrow(mediaSummary), function(x){
  mediaSummary$change[x] * boundary_ele_comp$MW[boundary_ele_comp$specie == mediaSummary$specie[x]]
  })




#grams per L of anabolic components
cellularComponents <- comp_by_cond$cultureMolarity[compositionFile_Extended$Class != "Maintenance ATP hydrolysis",] * as.numeric(compositionFile_Extended$MW[compositionFile_Extended$Class != "Maintenance ATP hydrolysis"]) %*% t(rep(1, ncol(comp_by_cond$cultureMolarity[compositionFile_Extended$Class != "Maintenance ATP hydrolysis",])))
cellularDF <- melt(t(cellularComponents))
colnames(cellularDF) <- c("condition", "metabolite", "density")
cellularDF$specie <- sapply(cellularDF$metabolite, function(x){
  compositionFile$Class[compositionFile$MetName == x]
  })
cellularDF <- cellularDF[!is.nan(cellularDF$density),]

cellularDFmelt <- melt(acast(cellularDF, condition ~ specie, value.var  = "density", fun.aggregate = sum))
colnames(cellularDFmelt) <- c("condition", "specie", "density")
cellularDFmelt$type <- "biomass"
cellularDFmelt <- data.table(cellularDFmelt)
cellularDFmelt$condition <- as.character(cellularDFmelt$condition); cellularDFmelt$specie <- as.character(cellularDFmelt$specie)

jointCompSummary <- rbind(mediaSummary[,list(condition, specie, density, type),], cellularDFmelt)
jointCompSummary[,Limitation := chemostatInfo$Limitation[chemostatInfo$NMR_name == condition], by = condition]
jointCompSummary$Limitation <- factor(jointCompSummary$Limitation, level = unique(chemostatInfo$Limitation))
jointCompSummary[,DR := chemostatInfo$DRgoal[chemostatInfo$NMR_name == condition], by = condition]

jointCompSummary$type <- factor(jointCompSummary$type, levels = c("uptake", "excretion", "biomass"))

barplot_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), legend.position = "bottom", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.key.width = unit(3, "line"), axis.text.x = element_text(angle = 90), strip.background = element_rect(fill = "darkgoldenrod1"))


class_comp_plot <- ggplot(jointCompSummary, aes(y = density, x = type, fill = factor(specie))) + barplot_theme + facet_grid(Limitation ~ DR, scales = "free_y") + scale_fill_discrete(name = "Class", guide = guide_legend(nrow = 4))
class_comp_plot + geom_bar(colour = "black", stat = "identity") + scale_y_continuous("Density (g/L)", expand = c(0.1, 0.1)) + scale_x_discrete("Type")
ggsave("speciesUtilization.pdf", height = 20, width = 14)


#class_comp_plot <- ggplot(mediaSummary, aes(y = density, factor(type), fill = factor(specie))) + pie_theme + scale_fill_brewer(name = "Class", palette = "Set1") + facet_wrap(~ condition)
#class_comp_plot + geom_bar(colour = "black", stat = "identity")

### for matching biomass fluxes which species should be combined into a single variance category ####
compositionFile_Extended$varCategory <- NA
compositionFile_Extended$varCategory[compositionFile_Extended$Class == "Amino Acid"] <- "AA flux"
compositionFile_Extended$varCategory[grep('Carbohydrates', compositionFile_Extended$Class)] <- "sugar polymer flux"
compositionFile_Extended$varCategory[grep('^d[A-Z]TP', compositionFile_Extended$MetName)] <- "dNTP flux"
compositionFile_Extended$varCategory[grep('^[A-Z]TP', compositionFile_Extended$MetName)] <- "NTP flux"
compositionFile_Extended$varCategory[grep('Sulfate', compositionFile_Extended$MetName)] <- "sulfate flux"
compositionFile_Extended$varCategory[grep('Maintenance ATP hydrolysis', compositionFile_Extended$Class)] <- "Maintenance ATP hydrolysis"
compositionFile_Extended$varCategory[compositionFile_Extended$MetName == "glycerol"] <- "glycerol washout"
compositionFile_Extended$varCategory[compositionFile_Extended$MetName == "polyphosphate"] <- "polyphosphate flux"
compositionFile_Extended$varCategory[compositionFile_Extended$Class == "Fatty Acids"] <- compositionFile_Extended$MetName[compositionFile_Extended$Class == "Fatty Acids"]

# Specify whether fluxes are drained off of the periphery of metabolism or  tracking an internal reaction or combination of reactions
# fatty acid incorporation was previously specified the latter way but 
compositionFile_Extended$FluxType <- "Boundary"


### Final output -> FBA_run ####

comp_by_cond$compositionFile <- compositionFile_Extended
comp_by_cond$biomassExtensionE <- composition_part
save(comp_by_cond, chemostatInfo, mediaSummary, file = "boundaryFluxes.Rdata")

# Table of measured % DW

write.table(class_sum_fraction_melt, file = "percentComposition.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


