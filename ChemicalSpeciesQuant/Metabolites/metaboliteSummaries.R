# Combine metabolomics datasets to return a consensus metabolite relative abundance
# and a conversion to absolute abundances where applicable

setwd("/Users/Sean/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Metabolites")
options(stringsAsFactors = F)

library(reshape2)
library(nlme)
library(data.table)
library(corpcor)
library(colorRamps)
library(dplyr)

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute) #impute missing values using knn imputation

### 4 datasets ###
#1) Reanalysis of original Boer data for metabolite relative abundances (and precision)
#2) Conversion to align chemostat and batch samples to use Yifan concentrations
#3) Yifan absolute quantification of metabolite abundance in batch culture
#4) Measurement of amino acid concentrations in N-labelled chemostats

# Import relationship between metabolite names, model name, and KEGG/ChEBI IDs
source("../../Yeast_genome_scale/cleanNalign.R", chdir = T)
# each metabolomics dataset is seperately matched to model IDs and then different metabolite measurements
# such as relative and absolute can be combined by virtue of their model IDs
# Dataset 2 & 3 are directly matched to each other before they are matched to the genome scale model

##### 1) Reanalysis of original Boer data for metabolite relative abundances (and precision) #####

tab_allBoer <- read.delim('./boerquant/boer_data_2.txt',header=T,sep='\t')
tab_allBoerNames <- read.delim('./boerquant/boer_data_2.txt',sep='\t', header = F, nrows = 1)
boerMeta <- read.delim("BoerMetabolites.txt")[-1,1:2]
boerHeatmap <- read.delim("BoerMetabolites.txt")[-1,-c(2:3)]

## Approach 1: get relative abundances by dividing through the reference condition
refTab <- tab_allBoer[ tab_allBoer$Exp.Ref == 'ref',-c(1,3,4,5,6)]
for (month in unique(refTab$Month)){
  monfil <- refTab$Month == month
  for (met in unique(refTab$Method)){
    metfil <- refTab$Method == met
    refTab[monfil & metfil,3:ncol(refTab)][1,] <- apply(refTab[monfil & metfil,3:ncol(refTab)],2,mean)
  }
}
refTab <- refTab[!duplicated(refTab[,1:2]),]
row.names(refTab) <- apply(refTab[,1:2],1,function(x){paste(x,collapse='')})

rIDBoer <- apply(tab_allBoer[,c('Month','Method')],1,function(x){paste(x,collapse='')})

### normalized tab_boer
ntab_allBoer <- tab_allBoer
ntab_allBoer[,8:113] <- ntab_allBoer[,8:113]/as.matrix(refTab[rIDBoer,3:108])

#heatmap.2(as.matrix(log2(ntab_allBoer[,8:113])), trace = "none")

boerRAmat <- as.matrix(log2(ntab_allBoer[ntab_allBoer$Exp.Ref == "exp",8:113]))
colnames(boerRAmat) <- tab_allBoerNames[8:113]
boerRAmat[boerRAmat == 0] <- NA

### determine whether the experimental sample was measured (> 300 IC)
boerExpMissing = tab_allBoer[,8:113][ntab_allBoer$Exp.Ref == "exp",] == 300

boerSampleInfo <- ntab_allBoer[ntab_allBoer$Exp.Ref == "exp",1:7]
boerSampleInfo$Condition <- paste(boerSampleInfo$Nutr, boerSampleInfo$Gr, sep = ".")

#plot(apply(boerRAmat, 1, median, na.rm = T) ~ factor(boerSampleInfo$Condition))

### relative abundances ###

### robust median polish - scaling factor of met w.r.t. median abundance

boerRAmatExpNA <- boerRAmat
boerRAmatExpNA[boerExpMissing] <- NA

medianDF <- data.table(cond = boerSampleInfo$Condition, method = boerSampleInfo$Method, medRA = apply(boerRAmatExpNA, 1, median, na.rm = T),
      medScaler = apply(boerRAmatExpNA - t(t(rep(1, nrow(boerRAmatExpNA)))) %*% apply(boerRAmatExpNA, 2, median, na.rm = T), 1, median, na.rm = T))
medianDF[,medianConsensus := mean(medRA), by = cond]
medianDF[,scalerConsensus := mean(medScaler), by = cond]
medianDF[,medianChange := medRA - medianConsensus]
medianDF[,scalerChange := medScaler - scalerConsensus]

plot(medianDF$medRA ~ medianDF$medScaler)

### Save the standard deviation by metabolite, the SD by metabolite*condition will be stored seperately
SDbyMet <- as.data.frame(matrix(ncol = 3, nrow = ncol(boerRAmat)))
colnames(SDbyMet) <- c("standard", "medianPolish", "medianScale")

### Save point estimates of metabolite relative abundance and SD using alternative normalizations
MetRAestimate <- NULL

### Save residuals of median scaled point estimates to look at covariance
residualStacks <- NULL

method_measurement <- data.frame(species = colnames(boerRAmat), chisquare_pvalue = NA) # pvalue for whether prep method effect p(missing)
for(i in 1:ncol(boerRAmat)){
  
  # determine whether a metabolite is not measured using one of the prep methods
  method_measurement$chisquare_pvalue[i] <- chisq.test(table(is.na(boerRAmat[,i]), boerSampleInfo$Method))$p.value
  
  reducedTable = data.table(cond = boerSampleInfo$Condition, method = boerSampleInfo$Method, RA = boerRAmat[,i], expMissing = boerExpMissing[,i])
  
  if(method_measurement$chisquare_pvalue[i] > 0.01){
    reducedTable$RA[is.na(reducedTable$RA)] <- 0
    }
  
  fitted_params = length(unique(reducedTable$cond[!is.na(reducedTable$RA)])) + nlevels(as.factor(boerSampleInfo$Method[!is.na(reducedTable$RA)])) - 1
  sd_infl <- sum(!is.na(reducedTable$RA))/(sum(!is.na(reducedTable$RA)) - fitted_params) # inflate the residual sd by this fraction to account for degrees of freedom consumed by fitting
  
  ####
  
  reducedTable$std_fitted = NA
  if(nlevels(as.factor(boerSampleInfo$Method[!is.na(reducedTable$RA)])) == 1){
    reducedTable$std_fitted[!is.na(reducedTable$RA)] = lm(data = reducedTable, RA ~ factor(cond))$fitted
    }else{
      reducedTable$std_fitted[!is.na(reducedTable$RA)] = lm(data = reducedTable, RA ~ factor(cond) + factor(method))$fitted
      }
  
  std_summary = reducedTable[,list(RA = mean(std_fitted, na.rm = T), SD = sd(RA - std_fitted, na.rm = T)*sd_infl, n = length(RA[!is.na(RA)])), by = cond]
  std_summary$metabolite = colnames(boerRAmat)[i]; std_summary$normalization = "standard"
  SDbyMet$standard[i] <- reducedTable[,sd(RA - std_fitted, na.rm = T)*sd_infl]
  std_summary$SD[std_summary$SD < SDbyMet$standard[i]/3] <- SDbyMet$standard[i]# overwrite within condition SD using metabolite SD if within met SD is less than average / 3
  
  ### normlization of RA
  
  ### median polish ###
  
  reducedTable$RAmedPol = boerRAmat[,i] - medianDF$medianChange
  reducedTable$RAmedPol[reducedTable$expMissing] <- reducedTable$RA[reducedTable$expMissing] # median polish is not used if experimental sample was below minimal IC
  
  reducedTable$std_fitted = NA
  if(nlevels(as.factor(boerSampleInfo$Method[!is.na(reducedTable$RAmedPol)])) == 1){
    reducedTable$std_fitted[!is.na(reducedTable$RAmedPol)] = lm(data = reducedTable, RAmedPol ~ factor(cond))$fitted
    }else{
      reducedTable$std_fitted[!is.na(reducedTable$RAmedPol)] = lm(data = reducedTable, RAmedPol ~ factor(cond) + factor(method))$fitted
      }
  
  median_summary = reducedTable[,list(RA = mean(std_fitted, na.rm = T), SD = sd(RAmedPol - std_fitted, na.rm = T)*sd_infl, n = length(RA[!is.na(RA)])), by = cond]
  median_summary$metabolite = colnames(boerRAmat)[i]; median_summary$normalization = "medianPolish"
  SDbyMet$medianPolish[i] <- reducedTable[,sd(RAmedPol - std_fitted, na.rm = T)*sd_infl]
  median_summary$SD[median_summary$SD < SDbyMet$medianPolish[i]/3] <- SDbyMet$medianPolish[i]# overwrite within condition SD using metabolite SD if within met SD is less than average / 3
  
  ### median scaling factors - deals better with missing values
  
  reducedTable$RAscaled = boerRAmat[,i] - medianDF$scalerChange
  reducedTable$RAscaled[reducedTable$expMissing] <- reducedTable$RA[reducedTable$expMissing] # median polish is not used if experimental sample was below minimal IC
  
  reducedTable$std_fitted = NA
  if(nlevels(as.factor(boerSampleInfo$Method[!is.na(reducedTable$RAscaled)])) == 1){
    reducedTable$std_fitted[!is.na(reducedTable$RAscaled)] = lm(data = reducedTable, RAscaled ~ factor(cond))$fitted
    }else{
      reducedTable$std_fitted[!is.na(reducedTable$RAscaled)] = lm(data = reducedTable, RAscaled ~ factor(cond) + factor(method))$fitted
      }
  
  
  scale_summary = reducedTable[,list(RA = mean(std_fitted, na.rm = T), SD = sd(RAscaled - std_fitted, na.rm = T)*sd_infl, n = length(RA[!is.na(RA)])), by = cond]
  scale_summary$metabolite = colnames(boerRAmat)[i]; scale_summary$normalization = "medianScale"
  SDbyMet$medianScale[i] <- reducedTable[,sd(RAscaled - std_fitted, na.rm = T)*sd_infl]
  scale_summary$SD[scale_summary$SD < SDbyMet$medianScale[i]/3] <- SDbyMet$medianScale[i]# overwrite within condition SD using metabolite SD if within met SD is less than average / 3
  
  MetRAestimate <- rbind(MetRAestimate, std_summary, median_summary, scale_summary)
  
  residualStacks <- rbind(residualStacks, data.frame(metabolite = colnames(boerRAmat)[i], sampleNum = 1:nrow(boerRAmat), residual = reducedTable$RA - reducedTable$std_fitted))
  }

### robust median scaling outperforms skipping normalization or just using medians
table(apply(SDbyMet, 1, which.min))
apply(SDbyMet^2, 2, sum)

### Calculate the relative abundance and standard error (labelled as SD) for each met * condition ###
metRA = acast(MetRAestimate[MetRAestimate[,normalization == "medianScale",],], formula = metabolite ~ cond, value.var = 'RA')
metSD = acast(MetRAestimate[MetRAestimate[,normalization == "medianScale",],], formula = metabolite ~ cond, value.var = 'SD')/sqrt(acast(MetRAestimate[MetRAestimate[,normalization == "medianScale",],], formula = metabolite ~ cond, value.var = 'n'))

metMatch <- chmatch(rownames(metRA), boerHeatmap$Metabolite)
sampleMatch <- chmatch(colnames(metRA), colnames(boerHeatmap))

reanalyzedBoer <- metRA[!is.na(metMatch), !is.na(sampleMatch)]
originalBoer <- boerHeatmap[metMatch[!is.na(metMatch)],sampleMatch[!is.na(sampleMatch)]] # this heatmap is mean-centered

reanalyzedBoer_meanCent <- reanalyzedBoer - t(t(apply(reanalyzedBoer, 1, mean, na.rm = T))) %*% rep(1, ncol(reanalyzedBoer))
#reanalyzedBoer[1:10,1:10]
#originalBoer[1:10,1:10]

combinedData <- data.frame(melt(reanalyzedBoer_meanCent), melt(originalBoer)[,2])
colnames(combinedData) <- c("Metabolite", "Condition", "Reanalyzed", "Boer")
ggplot(combinedData, aes(x = Reanalyzed, y = Boer)) + geom_hex(bins = 100) # the offset in the boer heatmap is slow DR c-lim


## determine the correlation of all residuals across metabolites ##
residualStackDF <- acast(residualStacks, formula = metabolite ~ sampleNum, value.var = "residual")
residualCorr <- cor(t(residualStackDF), use = "pairwise.complete.obs")

heatmap.2(residualCorr, trace = "none") # correlation of residuals

residualStackDF <- impute.knn(residualStackDF, rowmax = 0.7)$data

shrunkCorr <- cor.shrink(t(residualStackDF)) # shrink residual correlation matrix towards identity
pdf("metaboliteSummaries/boerCorr.pdf", height = 10, width = 10)
heatmap.2(shrunkCorr, trace = "none", symbreaks = T, col = blue2yellow(100)) 
dev.off()

partialCorrs <- cor2pcor(shrunkCorr) # partial correlation matrix
diag(partialCorrs) <- NA
rownames(partialCorrs) <- colnames(partialCorrs) <- rownames(shrunkCorr)
pdf("metaboliteSummaries/boerPartialcorr.pdf", height = 10, width = 10)
heatmap.2(partialCorrs, trace = "none", symbreaks = T, col = blue2yellow(100))
dev.off()
####

metRA = impute.knn(metRA)$data

for(a_row in 1:nrow(metSD)){ # use metabolite specific standard deviation as the sd for non-determined conditions
  metSD[a_row, is.na(metSD[a_row,])] <- SDbyMet$medianScale[a_row]
  }

write.table(metRA, "metaboliteSummaries/boerMean.tsv", quote = F, col.names = T, row.names = T, sep = "\t")
write.table(metSD, "metaboliteSummaries/boerSD.tsv", quote = F, col.names = T, row.names = T, sep = "\t")
write.table(shrunkCorr, "metaboliteSummaries/boerCorr.tsv", quote = F, col.names = T, row.names = T, sep = "\t")

boerSummary <- melt(metRA, value.name = "log2_RA") %>% tbl_df() %>% inner_join(
  melt(metSD, value.name = "log2_CV") %>% tbl_df()) %>% select(compound = Var1, condition = Var2, log2_RA, log2_CV) %>%
  mutate(compound = as.character(compound), condition = as.character(condition))

### Now relate these species to the model based on their name to attribute a model name and other systematic IDs ###

boerMeta_expanded <- lapply(1:nrow(boerMeta), function(i){
  data.frame(boerRow = i, expand.grid(Metabolite = strsplit(boerMeta$Metabolite[i], split = '/')[[1]], KEGG = strsplit(boerMeta$KEGG[i], split = '/')[[1]]))
})
boerMeta_expanded <- do.call("rbind", boerMeta_expanded)
boerMeta_expanded[,2:3] <- apply(boerMeta_expanded[,2:3], c(1,2), as.character)

boerMeta_expanded$CHEBI <- sapply(boerMeta_expanded$KEGG,Kegg2chebi)

# match the few yet unmatched by name
boerMeta_expanded$CHEBI[is.na(boerMeta_expanded$CHEBI)] <- sapply(boerMeta_expanded$Metabolite[is.na(boerMeta_expanded$CHEBI)],MatchName2Chebi)
boerMeta_expanded$fuzCHEBI <- sapply(boerMeta_expanded$CHEBI,Chebi2fuzchebi) # the primary chebi synonym corresponding to each CHEBI
boerMeta_expanded$fuzmCHEBI <- sapply(boerMeta_expanded$Metabolite,MatchName2Chebi) # best string match of metabolite names to a chebi metabolite name

# match IDs by KEGG-
boerMeta_expanded_KEGGmatch <- boerMeta_expanded %>% tbl_df() %>% left_join(listTID, by = "KEGG") %>%
  filter(!is.na(SpeciesName))
boerMeta_expanded_KEGGmatch %>% View()

# match IDs by CHEBI
boerMeta_resid <- boerMeta_expanded %>% filter(!(KEGG %in% boerMeta_expanded_KEGGmatch$KEGG))
boerMeta_expanded_CHEBImatch <- left_join(boerMeta_resid, listTID, by = "fuzCHEBI") %>% filter(!is.na(SpeciesName))
boerMeta_expanded_CHEBImatch %>% View()

# manually annotate remainder or ignore them if they are not in the model
# D-gluconate, glycerate, N-acetyl-glutamine, & dimethylglycine are not currently used

boerMeta_resid <- boerMeta_expanded %>% filter(!((KEGG %in% boerMeta_expanded_KEGGmatch$KEGG) | (fuzCHEBI %in% boerMeta_expanded_CHEBImatch$fuzCHEBI)))

# Combine metabolites back together
boerMeta_annotated <- rbind(
  boerMeta_expanded_KEGGmatch %>% select(boerRow, Metabolite, KEGG, CHEBI = CHEBI.y, SpeciesName, SpeciesType),
  boerMeta_expanded_CHEBImatch %>% select(boerRow, Metabolite, KEGG=KEGG.y, CHEBI = CHEBI.y, SpeciesName, SpeciesType),
  boerMeta_resid %>% select(boerRow:CHEBI) %>% mutate(SpeciesName = NA, SpeciesType = NA) 
) %>% arrange(boerRow)




##### 2) Conversion to align chemostat and batch samples to use Yifan concentrations #####

### setup exactive data and formatting ###

Exactive <- read.table("absolute_quant_compare/exactiveSummary.csv", sep = ",", header = F)
ExactiveDesign <- data.table(sample = as.character(Exactive[1,-c(1:13)]))
ExactiveDesign <- ExactiveDesign[ExactiveDesign[,grep('chemo|batchL|batchH', sample),],]
ExactiveDesign[,sampleType := regmatches(sample, regexpr('chemo|batchL|batchH', sample)),]

ExactiveDesign$dilution <- sapply(ExactiveDesign$sample, function(x){
  tail_vec <- strsplit(x, split = 'chemo|batchL|batchH')[[1]][2]
  regmatches(tail_vec, regexpr('[0-9.]+', tail_vec))
  })

colnames(Exactive) <- Exactive[1,]; Exactive <- Exactive[-1,]
Exactive <- Exactive[Exactive$note != "C12 PARENT",]
ExactiveMatrix <- Exactive[,-c(1:13)]
ExactiveMatrix <- ExactiveMatrix[,order(ExactiveDesign$sample)]
ExactiveDesign <- ExactiveDesign[order(ExactiveDesign$sample),]

ExactiveMatrix <- as.matrix(apply(ExactiveMatrix, c(1,2), as.numeric))
rownames(ExactiveMatrix) <- Exactive$compound
ExactiveMatrix[ExactiveMatrix < 32] <- NA
ExactiveMatrix <- ExactiveMatrix[,colnames(ExactiveMatrix) %in% ExactiveDesign$sample]

### setup max data and formatting ###

Max <- read.table("absolute_quant_compare/Max_2013_05_03-highestPeaks.csv", sep = ",", header = F)
MaxDesign <- data.table(sample = as.character(Max[1,-(1:5)]))
MaxDesign <- MaxDesign[MaxDesign[,grep('chemo|batchL|batchH', sample),],]
MaxDesign[,sampleType := regmatches(sample, regexpr('chemo|batchL|batchH', sample)),]

MaxDesign$dilution <- sapply(MaxDesign$sample, function(x){
  tail_vec <- strsplit(x, split = 'chemo|batchL|batchH')[[1]][2]
  regmatches(tail_vec, regexpr('[0-9.]+', tail_vec))
  })

colnames(Max) <- Max[1,]; Max <- Max[-1,]
Max <- Max[Max$chosen == TRUE,]

MaxMatrix <- Max[,-c(1:5)]
MaxMatrix <- MaxMatrix[,colnames(MaxMatrix) %in% MaxDesign$sample][,order(MaxDesign$sample)]
MaxDesign <- MaxDesign[order(MaxDesign$sample),]
MaxMatrix <- as.matrix(apply(MaxMatrix, c(1,2), as.numeric))
rownames(MaxMatrix) <- Max$Xcompound
MaxMatrix[MaxMatrix < 32] <- NA

### Combine max and exactive data into a single matrix - point estimation of
# chemostat logIC and regression of dilution on batch culture to determine
# ratio of chemostat to full concentration batch ###

totalMatrix <- log2(rbind(ExactiveMatrix, MaxMatrix))

metScaling <- data.table(compound = rownames(totalMatrix), chemostatLogAbund = apply(totalMatrix[,ExactiveDesign$sampleType == "chemo"], 1, mean, na.rm = TRUE),
  chemostatLogSD = apply(totalMatrix[,ExactiveDesign$sampleType == "chemo"], 1, sd, na.rm = TRUE))

heatmap.2(totalMatrix - metScaling[,chemostatLogAbund,], trace = "none")

met_scale_lm <- data.frame(intercept = rep(NA, length(totalMatrix[,1])), slope = NA)
for(metN in 1:length(totalMatrix[,1])){
  #met_scale_lm[metN,] <- lm(totalMatrix[,MaxDesign$sampleType == "batchL"][metN,] ~ log2(as.numeric(unlist(MaxDesign[MaxDesign$sampleType == "batchL",]$dilution))))$coef
  met_scale_lm[metN,] <- lm(totalMatrix[,MaxDesign$sampleType %in% c("batchL", "batchH")][metN,] ~ log2(as.numeric(unlist(MaxDesign[MaxDesign$sampleType %in% c("batchL", "batchH"),]$dilution))) + ifelse(MaxDesign[MaxDesign$sampleType %in% c("batchL", "batchH"),]$sampleType == "batchH", 1, 0))$coef[1:2]
  }
met_scale_lm$slope[met_scale_lm$slope < 0.6] <- 0.6
met_scale_lm$slope[met_scale_lm$slope > 1] <- 1

metScaling$batchInt <- met_scale_lm$int
metScaling$batchSlope <- met_scale_lm$slope

metScaling[,logScaling := (chemostatLogAbund - met_scale_lm$intercept)/met_scale_lm$slope - log2(0.73),]

##### 3) Yifan absolute quantification of metabolite abundance in batch culture #####

yifanConc <- read.delim("yeast absolute concentration yifan 1.1.txt")

equivalent_compounds <- c("2,3-diphosphoglycerate" = "2_3-Diphosphoglyceric acid", "Acetyl-CoA" = "acetyl-CoA", "Alanine" = "alanine", "D-Gluconic acid" = "D-gluconate", 
  "dihydroxyacetone-phosphate" = "dihydroxy-acetone-phosphate", "fructose-1-6-bisphosphate" = "fructose bisphosphate", "fumarate" = "fumarate 3-3", "leucine" = "leucine/isoleucine",
  "isoleucine" = "leucine/isoleucine", "N-acetyl-glucosamine-1-phosphate" = "N-acetyl-glucosamine-1/6-phosphate", "N-acetyl-glucosamine-6-phosphate" = "N-acetyl-glucosamine-1/6-phosphate",
  "NADP" = "NADP+", "Proline" = "proline", "PRPP" = "5-phosphoribosyl-1-pyrophosphate", "ribose 1-phosphate" = "ribose-phosphate", "ribose 5-phosphate" = "ribose-phosphate",
  "xylulose 5-phosphate" = "ribose-phosphate", "xylose 5-phosphate" = "ribose-phosphate", "sedoheptulose-7-phosphate" = "D-sedoheptulose-1/7-phosphate", "Succinate" = "succinate",
  "trehalose" = "trehalose/sucrose", "Valine" = "valine")

yifanConc$c_lim_scaling <- NA

for(i in 1:length(yifanConc[,1])){
  
  if(yifanConc$Compound[i] %in% metScaling$compound){
    yifanConc$c_lim_scaling[i] <- metScaling$logScaling[metScaling$compound == yifanConc$Compound[i]]
    } else if(yifanConc$Compound[i] %in% names(equivalent_compounds)){
      yifanConc$c_lim_scaling[i] <- metScaling$logScaling[metScaling$compound == equivalent_compounds[names(equivalent_compounds) == yifanConc$Compound[i]]]
      }
  }

yifanConc$c_lim_scaling[yifanConc$Compound == "Alanine"] <- 0 # low-ball estimate because measured difference is not physiological - makes rate of metabolic washout feasible

yifanConc$c_lim_conc <- yifanConc$Glucose * 2^yifanConc$c_lim_scaling / 1000 #concentation in M

yifanConc <- yifanConc[!is.na(yifanConc$c_lim_conc),]

yifanConc <- yifanConc %>% tbl_df() %>% mutate(condition = "C0.30") %>% select(compound = Compound, Instrument, Mixed, ChEBI, KEGG, condition, concentration = c_lim_conc)

#write.table(yifanConc, file = "metaboliteSummaries/yeast_absolute_concentration_chemo.txt", sep = "\t", quote = F, col.names = TRUE, row.names = F)

# Pass model IDs to Yifan compounds

yifan_meta <- yifanConc %>% select(compound, Instrument, Mixed, CHEBI = ChEBI, KEGG)

# Match each metabolites KEGG compounds to the model and find those IDs that match
yifan_meta$KEGGmatch <- sapply(yifan_meta$KEGG, function(x){
  KEGGs <- strsplit(x, split = ', ')[[1]]
  KEGGreturn <- KEGGs[KEGGs %in% listTID$KEGG]
  if(length(KEGGreturn) > 1){
   paste(KEGGreturn, collapse = '_') 
  }else if(length(KEGGreturn) == 0){
    NA
  }else{
    KEGGreturn
  }
})

# For metabolites without a KEGG match, try the same with CHEBI
yifan_meta$CHEBImatch <- sapply(yifan_meta$CHEBI, function(x){
  CHEBIs <- strsplit(x, split = ', ')[[1]]
  CHEBIreturn <- CHEBIs[CHEBIs %in% c(listTID$CHEBI, listTID$fuzCHEBI)]
  if(length(CHEBIreturn) > 1){
   paste(CHEBIreturn, collapse = '_') 
  }else if(length(CHEBIreturn) == 0){
    NA
  }else{
    CHEBIreturn
  }
})

# orphaned metabolites - Gluconate (gluconic acid), Hydroxyisocaproic acid, N-acetyl-glutamine, N-Acetyl-L-alanine, xylose 5P

yifan_meta %>% filter(is.na(KEGGmatch) & is.na(CHEBImatch))
table(!is.na(yifan_meta$KEGGmatch), !is.na(yifan_meta$CHEBImatch))

# add model IDs based on matched KEGG and CHEBI IDs

# KEGG
yifan_meta_KEGG <- yifan_meta %>% filter(!is.na(KEGGmatch)) %>% select(compound, Instrument, Mixed, KEGG = KEGGmatch)
yifan_meta_KEGG_unfold <- lapply(1:nrow(yifan_meta_KEGG), function(i){
  data.frame(yifan_meta_KEGG[i,colnames(yifan_meta_KEGG) != 'KEGG'], KEGG = strsplit(yifan_meta_KEGG[i,'KEGG'] %>% unlist(), split = '_')[[1]])
})
yifan_meta_KEGG_unfold <- do.call("rbind", yifan_meta_KEGG_unfold)
yifan_meta_KEGG_unfold <- yifan_meta_KEGG_unfold %>% left_join(listTID)

yifan_meta_KEGG_unfold <- yifan_meta_KEGG_unfold %>% select(compound, Instrument, SpeciesType)

# CHEBI where KEGG is not matched
yifan_meta_CHEBI <- yifan_meta %>% filter(is.na(KEGGmatch), !is.na(CHEBImatch)) %>% select(compound, Instrument, Mixed, CHEBI = CHEBImatch)

# pull down the CHEBI or fuzCHEBI matches
model_yifanCHEBImatches <- listTID %>% filter(CHEBI %in% yifan_meta_CHEBI$CHEBI | fuzCHEBI %in% yifan_meta_CHEBI$CHEBI)

yifan_meta_CHEBI_unfold <- lapply(1:nrow(yifan_meta_CHEBI), function(i){
  CHEBIs <- strsplit(yifan_meta_CHEBI$CHEBI[i], split = "_")[[1]]
  cbind(model_yifanCHEBImatches %>% filter(CHEBI %in% CHEBIs | fuzCHEBI %in% CHEBIs), 
        yifan_meta_CHEBI[i,] %>% select(-CHEBI))
})
yifan_meta_CHEBI_unfold <-  do.call("rbind", yifan_meta_CHEBI_unfold)

yifan_meta_CHEBI_unfold <- yifan_meta_CHEBI_unfold %>% select(compound, Instrument, SpeciesType)

rbind(yifan_meta_KEGG_unfold, yifan_meta_CHEBI_unfold)


##### 4) Measurement of amino acid concentrations in N-labelled chemostats #####

AA_absolute_quant <- read.delim('./absolute_quant_compare/N15chemoAA/N15absoluteQuant.txt') %>% tbl_df()

# separate leucine and isoleucine into discrete metabolites
AA_absolute_quant <- AA_absolute_quant %>% filter(compound != "leucine-isoleucine") %>% rbind(
  AA_absolute_quant[AA_absolute_quant$compound == "leucine-isoleucine",] %>% mutate(compound = "leucine"),
  AA_absolute_quant[AA_absolute_quant$compound == "leucine-isoleucine",] %>% mutate(compound = "isoleucine"))


# Pass model IDs to N15-AA quant experiment (based on pre-specified KEGG IDs)
N15AA_mets_meta <- AA_absolute_quant %>% select(compound, KEGG) %>% unique()

N15AA_mets_meta_expanded <- lapply(1:nrow(N15AA_mets_meta), function(i){
  data.frame(row = i, expand.grid(Metabolite = strsplit(N15AA_mets_meta$compound[i], split = '/')[[1]], KEGG = strsplit(N15AA_mets_meta$KEGG[i], split = '/')[[1]]))
})
N15AA_mets_meta_expanded <- do.call("rbind", N15AA_mets_meta_expanded)

N15AA_mets_meta_expanded <- N15AA_mets_meta_expanded %>% left_join(listTID) 
N15AA_mets_meta_expanded %>% View()

N15AA_mets_meta_expanded %>% select(row, KEGG, CHEBI, SpeciesName, SpeciesType)


##### Source cleanNalign so that tIDs corresponding to each metabolite can be found #### 



boerSummary %>% rowwise() %>% mutate(chebi = MatchName2Chebi(unique(compound),1))

boerSummary %>% group_by(compound) %>% mutate(chebi = MatchName2Chebi(unique(compound),1))
name <- "1,3-diphopshateglycerate"

sapply(AA_absolute_quant$compound, function(x){MatchName2Chebi(x,1)})

listTID

yifanConc %>% tbl_df() %>% rowwise() %>% mutate(chebi = MatchName2Chebi(Compound,1))

