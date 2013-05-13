##### Combine experimental data on yeast composition, nutrient intake and excretion and expected composition into a set of condition-specific boundary fluxes #####

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(gplots)
library(reshape2)
library(data.table)


### Default composition function ####

compositionFile <- read.csv2("../Yeast_comp.csv", sep = ",", stringsAsFactors = FALSE)

modelMetComp <- read.table("stoiMetsComp.tsv", header = TRUE)

compComp <- data.frame(t(sapply(compositionFile$AltName, function(x){
  incorporated <- sub('TP', 'MP', x)
  output <- unlist(modelMetComp[modelMetComp[,2] == incorporated,][1,])
  output[2] <- x
  output
  }))) #elemental composition each macromolecule in biomass function
  # for NTPs and dNTPs the relevent stoichiometry for determining the fractional contribution to dry-weight is (d)NMPs

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
RNAdefaultcomp <- compositionFile[compositionFile$MetName %in% c("ATP", "UTP", "CTP", "GTP"),]
comp_by_cond$moles_per_cell[compositionFile$MetName %in% c("ATP", "UTP", "CTP", "GTP"),] <- t(t(RNAdefaultcomp$weightPerUn/sum(RNAdefaultcomp$weightPerUn)/RNAdefaultcomp$MW)) %*% t(totalRNApercell)

##### Total DNA #####

## genome size + genome size * fraction of cells not in G1 # from Brauer

buddingFrac = 0.936 - 1.971*chemostatInfo$actualDR #brauer 2008 relationship between unbudded fraction and growth rate
genomeLength = 12157105 # yeast genome length from SGD http://www.yeastgenome.org/cache/genomeSnapshot.html
GContent = 0.383 #http://bionumbers.hms.harvard.edu//bionumber.aspx?id=102126&ver=0
avogadros = 6.02214e23
 
comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) %in% c("dATP", "dTTP"),] <- rbind((genomeLength * (1-GContent) / avogadros)*(1+(1-buddingFrac)), (genomeLength * (1-GContent) / avogadros)*(1+(1-buddingFrac))) * 10^12
comp_by_cond$moles_per_cell[rownames(comp_by_cond$grams_per_cell) %in% c("dGTP", "dCTP"),] <- rbind((genomeLength * GContent / avogadros)*(1+(1-buddingFrac)), (genomeLength * GContent / avogadros)*(1+(1-buddingFrac))) * 10^12

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
comp_by_cond$grams_per_cell[rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c,] <- t(t(compositionFile$weightPerUn[rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c]/sum(compositionFile$weightPerUn[compositionFile$Class != "Energy Balance" & rowSums(is.na(comp_by_cond$grams_per_cell)) == n_c])
)) %*% t(cond_dryweight - colSums(comp_by_cond$grams_per_cell[compositionFile$Class == "Amino Acid",]))

comp_by_cond$moles_per_cell[compositionFile$Class != "Energy Balance",] <- comp_by_cond$grams_per_cell[compositionFile$Class != "Energy Balance",]/t(t(compositionFile$MW[compositionFile$Class != "Energy Balance"])) %*% rep(1, n_c)

# calculate the maintenance flux from dry weight and polymerization costs from molarity

# maintenance flux - relative to DW
maintFlux <- 0.0001 * colSums(comp_by_cond$grams_per_cell)
maintVec <- c(1, -1, -1); names(maintVec) <- c("ATP", "ADP", "phosphate")

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

composition_part <- cast(composition_part, formula = name ~ compound, sum, value = "stoi")
composition_part <- composition_part[,colnames(composition_part) != "NA"]

composition_part$name <- factor(composition_part$name, levels = compositionFile$AltName)
composition_part <- composition_part[order(composition_part$name),]

energy_balance <- t(comp_by_cond$moles_per_cell) %*% as.matrix(composition_part[,-1])
colnames(energy_balance) <- colnames(composition_part)[-1]

for(maint_cmpd in rownames(maintMatrix)){
  if(maint_cmpd %in% colnames(energy_balance)){
    energy_balance[,colnames(energy_balance) == maint_cmpd] <- maintMatrix[rownames(maintMatrix) == maint_cmpd,] + energy_balance[,colnames(energy_balance) == maint_cmpd]
    }
}

added_cmpds <- data.frame(maintMatrix[!(rownames(maintMatrix) %in% colnames(energy_balance)),])
colnames(added_cmpds) <- rownames(maintMatrix)[!(rownames(maintMatrix) %in% colnames(energy_balance))]

energy_balance <- cbind(energy_balance, added_cmpds)

### add energetic balances back into composition fxn
comp_by_cond$moles_per_cell[rownames(comp_by_cond$moles_per_cell) %in% colnames(energy_balance),] <- comp_by_cond$moles_per_cell[rownames(comp_by_cond$moles_per_cell) %in% colnames(energy_balance),] + t(energy_balance[,colnames(energy_balance) %in% compositionFile$AltName])

comp_by_cond$moles_per_cell <- rbind(comp_by_cond$moles_per_cell, t(energy_balance[,!(colnames(energy_balance) %in% compositionFile$AltName)]))



### rewrite compositionFile so that added species in polymerization cost have appropriate names in FBA model
added_cmpds <- data.frame(MetName = colnames(energy_balance[,!(colnames(energy_balance) %in% compositionFile$AltName)]), StoiCoef = NA, AltName = colnames(energy_balance[,!(colnames(energy_balance) %in% compositionFile$AltName)]), Class = "Energy Balance", Abbreviated = colnames(energy_balance[,!(colnames(energy_balance) %in% compositionFile$AltName)]), Polymerization_Cost = NA, MW = NA, weightPerUn = NA)
compositionFile_Extended = data.frame(t(cbind(t(compositionFile), t(added_cmpds))))
write.table(compositionFile_Extended, file = "../Yeast_comp_energy.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)





# generate weight per cellular volume
apply((comp_by_cond$grams_per_cell * 10^-12)/(t(t(rep(1, length(compositionFile[,1])))) %*% chemostatInfo$medcellVol * 10^-15), 2, sum, na.rm = TRUE)/10 #grams per 100mL of macromolecules - chunky but reasonable

# determine the fluxes into macromolecule biosynthesis and energy per volume
comp_by_cond$intacellularMolarity = (comp_by_cond$moles_per_cell * 10^-12)/(t(t(rep(1, length(compositionFile_Extended[,1])))) %*% chemostatInfo$medcellVol * 10^-15)



comp_by_cond$anabolicFlux = comp_by_cond$intacellularMolarity * (t(t(rep(1, length(compositionFile_Extended[,1])))) %*% chemostatInfo$actualDR) #moles/L-hr 

comp_by_cond$cultureMolarity <- comp_by_cond$intacellularMolarity * t(t(rep(1, length(compositionFile_Extended[,1])))) %*% chemostatInfo$VolFrac_mean/1000 #moles/L culture








########### Boundary fluxes from nutrient uptake or metabolite excretion ############

SPM <- read.delim("Media/measuredNutrients.tsv")
NMR <- read.delim("Media/mediaComposition_NMR.tsv")
nutrientFile <- read.delim("../Yeast_genome_scale/Boer_nutrients.txt")[1:6,1:6]; nutrientCode <- data.frame(nutrient = colnames(nutrientFile)[-1], shorthand = c("n", "p", "c", "L", "u")) #molar concentration of nutrients
nutrientFileConv <- data.frame(mediaName = colnames(nutrientFile)[-1], shortName = c("n", "p", "c", "L", "u"))

mediaSummary <- NULL

### Nutrient not present here will be solely bounded by mediaConcentration*DR

### phosphate uptake
SPM$condition <- factor(SPMtmp$condition, levels = chemostatInfo$condition)
SPM <- SPM[order(SPM$condition),]
SPM_out <- data.frame(condition = SPM$condition, specie = "phosphate", change = SPM$phosphateUptake, sd = SPM$phosphateSD, lb = 0, ub = SPM$ceiling, type = "uptake")
mediaSummary <- rbind(mediaSummary, SPM_out)

### NMR data - glucose -> EtOH, Ac, glycerol 
speciesOfInterest <- data.frame(NMRname = c("Glucose", "Ethanol", "Acetate", "Glycerol"), modelName = c("D-glucose", "ethanol", "acetate", "glycerol"), type = c("uptake", rep("excretion", 3)))

for(a_specie_n in 1:length(speciesOfInterest[,1])){
  NMR_data_subset <- NMR[(NMR$peak == speciesOfInterest$NMRname[a_specie_n]) & NMR$condition %in% chemostatInfo$condition,]
  NMR_data_subset$condition <- factor(NMR_data_subset$condition, levels = chemostatInfo$condition)
  NMR_data_subset <- NMR_data_subset[order(NMR_data_subset$condition),]
  
  if(speciesOfInterest$type[a_specie_n] == "uptake"){
    mediaComp <- nutrientFile[nutrientFile$X == speciesOfInterest$modelName[a_specie_n],-1]
    upper_bound <- unname(unlist(sapply(as.character(NMR_data_subset$condition), function(x){mediaComp[names(mediaComp) == nutrientFileConv$mediaName[nutrientFileConv$shortName == strsplit(x, '')[[1]][1]]]})))
    NMR_out <- data.frame(condition = NMR_data_subset$condition, specie = speciesOfInterest$modelName[a_specie_n], change = (upper_bound - NMR_data_subset$estimate/1000), sd = NMR_data_subset$se/1000, lb = 0, ub = upper_bound, type = "uptake")
    
  }else{
    NMR_out <- data.frame(condition = NMR_data_subset$condition, specie = speciesOfInterest$modelName[a_specie_n], change = NMR_data_subset$estimate/1000, sd = NMR_data_subset$se/1000, lb = 0, ub = Inf, type = "excretion")
    }
  
  mediaSummary <- rbind(mediaSummary, NMR_out)
  }
mediaSummary <- data.table(mediaSummary)

boundary_ele_comp <- mediaSummary[,modelMetComp[modelMetComp$name == specie,][1,],by = specie]
boundary_ele_comp[boundary_ele_comp$specie == "D-glucose",]$C <- 6; boundary_ele_comp[boundary_ele_comp$specie == "D-glucose",]$H <- 12; boundary_ele_comp[boundary_ele_comp$specie == "D-glucose",]$O <- 6
boundary_ele_comp[is.na(boundary_ele_comp)] <- 0
boundary_ele_comp <- boundary_ele_comp[,-grep('Fe|K|Na|R', colnames(boundary_ele_comp)),with = F]

boundary_ele_comp$MW <- rowSums(boundary_ele_comp[,atomicMasses$element,with = F] * rep(1, nrow(boundary_ele_comp)) %*% t(atomicMasses$mass))

mediaSummary$density <- sapply(1:nrow(mediaSummary), function(x){
  mediaSummary$change[x] * boundary_ele_comp$MW[boundary_ele_comp$specie == mediaSummary$specie[x]]
  })




#grams per L of anabolic components
cellularComponents <- comp_by_cond$cultureMolarity[compositionFile_Extended$Class != "Energy Balance",] * as.numeric(compositionFile_Extended$MW[compositionFile_Extended$Class != "Energy Balance"]) %*% t(rep(1, ncol(comp_by_cond$cultureMolarity[compositionFile_Extended$Class != "Energy Balance",])))
cellularComponents <- cellularComponents[!(rownames(cellularComponents) %in% c("ATP", "GTP")),]
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
jointCompSummary[,Limitation := chemostatInfo$limitation[chemostatInfo$condition == condition], by = condition]
jointCompSummary[,DR := chemostatInfo$DRgoal[chemostatInfo$condition == condition], by = condition]
jointCompSummary <- jointCompSummary[!(jointCompSummary$condition %in% c("p0.05H1", "p0.05H2")),]

jointCompSummary$type <- factor(jointCompSummary$type, levels = c("uptake", "excretion", "biomass"))

barplot_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), legend.position = "bottom", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.key.width = unit(3, "line"), axis.text.x = element_text(angle = 90), strip.background = element_rect(fill = "darkgoldenrod1"))


class_comp_plot <- ggplot(jointCompSummary, aes(y = density, x = type, fill = factor(specie))) + barplot_theme + facet_grid(Limitation ~ DR, scales = "free_y") + scale_fill_discrete(name = "Class", guide = guide_legend(nrow = 3))
class_comp_plot + geom_bar(colour = "black", stat = "identity") + scale_y_continuous("Density (g/L)", expand = c(0.1, 0.1)) + scale_x_discrete("Type")
ggsave("speciesUtilization.pdf", height = 20, width = 14)


#class_comp_plot <- ggplot(mediaSummary, aes(y = density, factor(type), fill = factor(specie))) + pie_theme + scale_fill_brewer(name = "Class", palette = "Set1") + facet_wrap(~ condition)
#class_comp_plot + geom_bar(colour = "black", stat = "identity")

save(comp_by_cond, chemostatInfo, mediaSummary, file = "boundaryFluxes.Rdata")


