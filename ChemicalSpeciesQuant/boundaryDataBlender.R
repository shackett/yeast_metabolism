##### Combine experimental data on yeast composition, nutrient intake and excretion and expected composition into a set of condition-specific boundary fluxes #####

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(reshape)

### Default composition function ####

compositionFile <- read.csv2("../Yeast_comp.csv", sep = ",", stringsAsFactors = FALSE)

modelMetComp <- read.table("stoiMetsComp.tsv", header = TRUE)

compComp <- data.frame(t(sapply(compositionFile$AltName, function(x){
  unlist(modelMetComp[modelMetComp[,2] == x,][1,])
  }))) #elemental composition each macromolecule in biomass function

compComp[compComp$name == "(1->3)-beta-D-glucan", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5) #one oxygen shared bc of condensation
compComp[compComp$name == "glycogen", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5)
compComp[compComp$name == "mannan", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5)
compComp[is.na(compComp)] <- 0
compComp[,-c(1:2)] <- apply(compComp[,-c(1:2)], c(1,2), as.numeric)
compComp <- compComp[,c(TRUE, TRUE, colSums(compComp[,-c(1:2)]) != 0)]

atomicMasses <- data.frame(element = c("C", "H", "N", "O", "P", "S"), mass = c(12.0107, 1.00794, 14.00674, 15.9994, 30.973761, 32.066))

compositionFile$MW <- t(t(compComp[,-c(1,2)])) %*% t(t(c(atomicMasses[,2])))
compositionFile$weightPerUn <- as.numeric(compositionFile$StoiCoef) * compositionFile$MW
#compositionFile$weight_per_t[compositionFile$Class == "Energy Balance"] <- NA

class_composition <- sapply(unique(compositionFile$Class), function(x){sum(compositionFile$weightPerUn[compositionFile$Class == x])}) * -1
class_composition <- data.frame(Category = names(class_composition[names(class_composition) != "Energy Balance"]), Abundance = unname(class_composition[names(class_composition) != "Energy Balance"]))
class_composition$CellComposition <- "Default"

pie_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), legend.position = "right", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(3, "line"), axis.text = element_blank(), axis.title = element_blank()) 


class_comp_plot <- ggplot(class_composition, aes(y = Abundance, factor(CellComposition), fill = factor(Category))) + pie_theme
class_comp_plot + geom_bar(stat = "identity") + coord_polar(theta = "y") + scale_fill_discrete(name = "Class")

class_comp_sep <- compositionFile[compositionFile$Class != "Energy Balance",]
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

chemostatInfo <- read.table("BulkComposition/chemostatDRmisc.tsv", sep = "\t", header = TRUE) #VolFrac_mean - uL cellular vol per mL media

comp_by_cond <- list()
tmp <- matrix(NA, ncol = length(chemostatInfo[,1]), nrow = length(compositionFile[,1]))
colnames(tmp) <- chemostatInfo$condition; rownames(tmp) <- compositionFile$MetName
comp_by_cond$moles_per_cell <- comp_by_cond$grams_per_cell <- tmp

n_c <- length(chemostatInfo[,1])

##### Total Protein ######

load('BulkComposition/protSpecQuantOut.Rdata')

proteinWeight <- dry_weight_perCell_DFmelt[dry_weight_perCell_DFmelt$variable == "Protein",] #pg protein per cell

aaRelAbunds <- compositionFile[compositionFile$Class == "Amino Acid",]
aaWeightFrac <- (aaRelAbunds$weightPerUn*-1)/sum(aaRelAbunds$weightPerUn*-1)

comp_by_cond$moles_per_cell[compositionFile$Class == "Amino Acid",] <- t(t(t(proteinWeight$value)) %*% t(aaWeightFrac/aaRelAbunds$MW))

dry_weight_perCell_plot <- ggplot(dry_weight_perCell_DFmelt, aes(x = factor(Condition), y = value, fill = variable)) + theme(axis.text.x = element_text(size = 4, face = "bold"), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + scale_y_continuous("pg protein/dry-material per cell", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(palette = "Set1")
dry_weight_perCell_plot + geom_bar(stat = "identity", position = "stack")

##### Total RNA ######

RNA_file <- read.delim("BulkComposition/RNA_abundance/RNAabund.csv", sep = ",", header = TRUE)
RNAconc <- sapply(chemostatInfo$condition, function(x){mean(RNA_file$RNAconc[RNA_file$condition == x])})
RNA_peruLcellVol <- RNAconc/chemostatInfo$VolFrac_mean #numerator is ug of RNA per mL of cells, denominator is uL of cells per mL of cells
totalRNApercell <- RNA_peruLcellVol/((10^-6)/(chemostatInfo$medcellVol * 10^-15)) * 10^6 #pg per cell
RNAdefaultcomp <- compositionFile[compositionFile$MetName %in% c("AMP", "UMP", "CMP", "GMP"),]
comp_by_cond$moles_per_cell[compositionFile$MetName %in% c("AMP", "UMP", "CMP", "GMP"),] <- t(t(RNAdefaultcomp$weightPerUn/sum(RNAdefaultcomp$weightPerUn)/RNAdefaultcomp$MW)) %*% t(totalRNApercell)

##### Total DNA #####

## genome size + genome size * fraction of cells not in G1 # from Brauer

buddingFrac = 0.936 - 1.971*chemostatInfo$actualDR #brauer 2008 relationship between unbudded fraction and growth rate
genomeLength = 12157105 # yeast genome length from SGD http://www.yeastgenome.org/cache/genomeSnapshot.html
GContent = 0.383 #http://bionumbers.hms.harvard.edu//bionumber.aspx?id=102126&ver=0
avogadros = 6.02214e23
 
comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) %in% c("dAMP", "dTMP"),] <- rbind((genomeLength * (1-GContent) / avogadros)*(1+(1-buddingFrac)), (genomeLength * (1-GContent) / avogadros)*(1+(1-buddingFrac))) * 10^12
comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) %in% c("dGMP", "dCMP"),] <- rbind((genomeLength * GContent / avogadros)*(1+(1-buddingFrac)), (genomeLength * GContent / avogadros)*(1+(1-buddingFrac))) * 10^12

##### Combining observed abundances with total dry weight per cell and inferring contributions of non-measured elements based upon the assumption that they remain in constant proportions, filling the rest of the dry weight ####

#scaling energy usage by cell weight relative to default weight.

comp_by_cond$grams_per_cell <- comp_by_cond$moles_per_cell * (t(t(compositionFile$MW)) %*% t(rep(1,n_c)))

cond_dryweight <- sapply(unique(dry_weight_perCell_DFmelt$Condition), function(x){
  sum(dry_weight_perCell_DFmelt$value[dry_weight_perCell_DFmelt$Condition == x])
  })#pg per cell

weightSummary <- comp_by_cond$grams_per_cell[rowSums(!is.na(comp_by_cond$grams_per_cell)) != 0,]; weightSummary[is.nan(weightSummary)] <- NA
summaryInfo <- compositionFile[rowSums(!is.na(comp_by_cond$grams_per_cell)) != 0,]
weightSummary <- rbind(weightSummary, cond_dryweight - colSums(weightSummary, na.rm = TRUE))
rownames(weightSummary)[length(weightSummary[,1])] <- "Residual"
summaryInfo <- summaryInfo[,colnames(summaryInfo) %in% c("MetName", "Class", "Abbreviated")]
summaryInfo <- rbind(summaryInfo, c("Residual", "Residual dry weight", "Residual"))

weightAndSum <- cbind(summaryInfo, weightSummary)
weightAndSum <- weightAndSum[order(weightAndSum$Class),]

for(i in unique(weightAndSum$Class)){
  weightAndSum[weightAndSum$Class == i,] <- weightAndSum[c(1:length(weightAndSum[,1]))[weightAndSum$Class == i][order(weightAndSum$p0.05[weightAndSum$Class == i], decreasing = TRUE)],]
  }


weightStackDF <- melt(weightAndSum)
colnames(weightStackDF)[4:5] <- c("Condition", "Abundance")

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_blank(), legend.position = "top", 
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(size = 20, angle = 70, vjust = 0.5, colour = "BLACK"), legend.key.width = unit(3, "line")) 

#plot condition specific composition from experimental measurements

class_comp_plot <- ggplot(weightStackDF, aes(y = Abundance, factor(Condition), fill = factor(Class))) + barplot_theme
class_comp_plot + geom_bar(colour = "black", stat = "identity") + scale_y_continuous("pg per cell", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(name = "Class", palette = "Set1")
ggsave(file = "condition_composition.pdf", height = 12, width = 20)

class_comp_plot <- ggplot(weightStackDF[!is.na(weightStackDF$Abundance),], aes(y = Abundance, factor(Condition), fill = factor(Class))) + barplot_theme
class_comp_plot + geom_bar(stat = "identity", position = "fill") + scale_y_continuous("Fraction of cellular material", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(name = "Class", palette = "Set1")
ggsave(file = "fractional_composition.pdf", height = 12, width = 20)


#fill in structural composition components based on residual dry weight
comp_by_cond$grams_per_cell[compositionFile$Class != "Energy Balance" & rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c,] <- t(t(compositionFile$weightPerUn[compositionFile$Class != "Energy Balance" & rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c]/sum(compositionFile$weightPerUn[compositionFile$Class != "Energy Balance" & rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c])
)) %*% t(cond_dryweight - colSums(comp_by_cond$grams_per_cell[compositionFile$Class == "Amino Acid",]))

#initially assume that ATP -> ADP + Pi flux is proportional to cellular dry weight

comp_by_cond$moles_per_cell[compositionFile$Class != "Energy Balance",] <- comp_by_cond$grams_per_cell[compositionFile$Class != "Energy Balance",]/t(t(compositionFile$MW[compositionFile$Class != "Energy Balance"])) %*% rep(1, n_c)
comp_by_cond$moles_per_cell[compositionFile$Class == "Energy Balance",] <- -1*t(t(as.numeric(compositionFile$StoiCoef[compositionFile$Class == "Energy Balance"]))) %*% t(cond_dryweight/sum(compositionFile$weightPerUn[compositionFile$Class != "Energy Balance"] * -1))



# generate weight per volume
apply((comp_by_cond$grams_per_cell * 10^-12)/(t(t(rep(1, length(compositionFile[,1])))) %*% chemostatInfo$medcellVol * 10^-15), 2, sum, na.rm = TRUE)/10 #grams per 100mL of macromolecules - chunky but reasonable

# determine the fluxes into macromolecule biosynthesis and energy per volume
comp_by_cond$intacellularMolarity = (comp_by_cond$moles_per_cell * 10^-12)/(t(t(rep(1, length(compositionFile[,1])))) %*% chemostatInfo$medcellVol * 10^-15)
comp_by_cond$anabolicFlux = comp_by_cond$intacellularMolarity * (t(t(rep(1, length(compositionFile[,1])))) %*% chemostatInfo$actualDR) #moles/L-hr 

comp_by_cond$cultureMolarity <- comp_by_cond$intacellularMolarity * t(t(rep(1, length(compositionFile[,1])))) %*% chemostatInfo$VolFrac_mean/1000 #moles/L culture

save(comp_by_cond, chemostatInfo, file = "boundaryFluxes.Rdata")

