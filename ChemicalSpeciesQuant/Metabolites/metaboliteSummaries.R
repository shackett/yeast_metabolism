# Combine metabolomics datasets to return a consensus metabolite relative abundance
# and a conversion to absolute abundances where applicable

setwd("/Users/Sean/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Metabolites")
options(stringsAsFactors = F)

# data manipulation
library(dplyr)
library(reshape2)
library(data.table)

library(nlme)
library(corpcor) # estimation of residual covariance matrix (with shrinkage)
library(gplots) # for generating heatmaps
library(colorRamps) # alternative heatmap color gradients
library(missMDA) # for doing cross-validation based estimation of number of significant PCs

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute) # impute missing values using knn imputation

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

# Raw data from boer et al. 2010
tab_allBoer <- read.delim('./boerquant/boer_data_2.txt',header=T,sep='\t')
tab_allBoerNames <- unname(unlist(read.delim('./boerquant/boer_data_2.txt',sep='\t', header = F, nrows = 1)))
# Published boer et al. 2010 heatmap and KEGG IDs for metabolites
# annoyingly these names don't perfectly match the raw data
boerMeta <- read.delim("BoerMetabolites.txt")[-1,1:2]
boerHeatmap <- read.delim("BoerMetabolites.txt")[-1,-c(2:3)]

# take care of that name mismatching 
raw_boer_names <- tab_allBoerNames[8:113]

unmatched_raw_names <- raw_boer_names[!(raw_boer_names %in% boerMeta$Metabolite)]
unmatched_hm_names <- boerMeta$Metabolite[!(boerMeta$Metabolite %in% raw_boer_names)]
# mismatched names when ordered, are correctly aligned
# double check that the mismatched compounds are the same following ordering
raw_hm_matches <- data.frame(Raw = sort(unmatched_raw_names), HM = sort(unmatched_hm_names)) 

raw_boer_names[!(raw_boer_names %in% boerMeta$Metabolite)] <- raw_hm_matches$HM[chmatch(unmatched_raw_names, raw_hm_matches$Raw)]

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
colnames(boerRAmat) <- raw_boer_names
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

### Reserve space for metabolities which are inferred in other ways ###

boerMeta <- boerMeta %>% rbind(data.frame(Metabolite = c("phosphate", "phosphoenolpyruvate"), KEGG = c("C00009", "C00074")))

# filter one slow-growth p-limited chemostat
metRA <- metRA[,!(colnames(metRA) %in% "PO4.0.061")]
metSD <- metSD[,!(colnames(metSD) %in% "PO4.0.061")]
boerSummary <- boerSummary %>% filter(condition != "PO4.0.061")

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
) %>% arrange(boerRow) %>% filter(!is.na(SpeciesName))

if(!all(table(boerMeta_annotated$SpeciesType) == 1)){
  degenerate_compounds <- names(which(table(boerMeta_annotated$SpeciesType) > 1))
  best_matched_compounds <- NULL
  for(a_cmpd in degenerate_compounds){
    # find the best match based upon string similarity
    sim_subset <- boerMeta_annotated %>% filter(SpeciesType == a_cmpd) %>% mutate(stringSim = levenshteinDist(Metabolite, SpeciesName)) %>%
      filter(stringSim == min(stringSim)) %>% select(-stringSim)
    if(nrow(sim_subset) != 1){
      # take the peak that is uniquely identified in case one compound is a mixture on a subset of instruments
      sim_subset <- sim_subset %>% group_by(boerRow) %>% mutate(nShared = length(boerMeta_annotated$boerRow[boerMeta_annotated$boerRow == boerRow])) %>%
        group_by() %>% filter(nShared == min(nShared)) %>% select(-nShared)
      if(nrow(sim_subset) != 1){
        warning("Multiple metabolite measurements must be combined before this point to generate a dataset level consensus")
        }
      }
    best_matched_compounds <- rbind(best_matched_compounds, sim_subset)
  }
  
  boerMeta_annotated <- boerMeta_annotated %>% filter(!(SpeciesType %in% degenerate_compounds)) %>%
    rbind(best_matched_compounds)
}


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

yifan_meta <- yifanConc %>% mutate(row = 1:nrow(yifanConc)) %>% select(row, compound, Instrument, Mixed, CHEBI = ChEBI, KEGG)

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
yifan_meta_KEGG <- yifan_meta %>% filter(!is.na(KEGGmatch)) %>% select(row, compound, Instrument, Mixed, KEGG = KEGGmatch)
yifan_meta_KEGG_unfold <- lapply(1:nrow(yifan_meta_KEGG), function(i){
  data.frame(yifan_meta_KEGG[i,colnames(yifan_meta_KEGG) != 'KEGG'], KEGG = strsplit(yifan_meta_KEGG[i,'KEGG'] %>% unlist(), split = '_')[[1]])
})
yifan_meta_KEGG_unfold <- do.call("rbind", yifan_meta_KEGG_unfold)
yifan_meta_KEGG_unfold <- yifan_meta_KEGG_unfold %>% left_join(listTID)

yifan_meta_KEGG_unfold <- yifan_meta_KEGG_unfold %>% select(row, Metabolite = compound, Instrument, SpeciesName, SpeciesType)

# CHEBI where KEGG is not matched
yifan_meta_CHEBI <- yifan_meta %>% filter(is.na(KEGGmatch), !is.na(CHEBImatch)) %>% select(row, compound, Instrument, Mixed, CHEBI = CHEBImatch)

# pull down the CHEBI or fuzCHEBI matches
model_yifanCHEBImatches <- listTID %>% filter(CHEBI %in% yifan_meta_CHEBI$CHEBI | fuzCHEBI %in% yifan_meta_CHEBI$CHEBI)

yifan_meta_CHEBI_unfold <- lapply(1:nrow(yifan_meta_CHEBI), function(i){
  CHEBIs <- strsplit(yifan_meta_CHEBI$CHEBI[i], split = "_")[[1]]
  cbind(model_yifanCHEBImatches %>% filter(CHEBI %in% CHEBIs | fuzCHEBI %in% CHEBIs), 
        yifan_meta_CHEBI[i,] %>% select(-CHEBI))
})
yifan_meta_CHEBI_unfold <-  do.call("rbind", yifan_meta_CHEBI_unfold)

yifan_meta_CHEBI_unfold <- yifan_meta_CHEBI_unfold %>% select(row, Metabolite = compound, Instrument, SpeciesName, SpeciesType)

yifan_meta_correspondence <- rbind(yifan_meta_KEGG_unfold, yifan_meta_CHEBI_unfold)

# Ensure that reverse mapping of t_IDs to metabolites is unique
if(!all(table(yifan_meta_correspondence$SpeciesType) == 1)){
  stop("Each tID must be matched to a single compound in this dataset")
  }


##### 4) Measurement of amino acid concentrations in N-labelled chemostats #####

AA_absolute_quant <- read.delim('./absolute_quant_compare/N15chemoAA/N15absoluteQuant.txt') %>% tbl_df()

# separate leucine and isoleucine into discrete metabolites
AA_absolute_quant <- AA_absolute_quant %>% filter(compound != "leucine-isoleucine") %>% rbind(
  AA_absolute_quant[AA_absolute_quant$compound == "leucine-isoleucine",] %>% mutate(compound = "leucine"),
  AA_absolute_quant[AA_absolute_quant$compound == "leucine-isoleucine",] %>% mutate(compound = "isoleucine"))

AA_absolute_quant <- AA_absolute_quant %>% mutate(concentration = concentration_mM/1000)

# Take the mean of replicates
AA_absolute_quant <- AA_absolute_quant %>% group_by(condition, compound, KEGG, DR) %>%
  summarize(concentration = mean(concentration)) %>% group_by()

# Pass model IDs to N15-AA quant experiment (based on pre-specified KEGG IDs)
N15AA_mets_meta <- AA_absolute_quant %>% select(compound, KEGG) %>% unique()

N15AA_mets_meta_expanded <- lapply(1:nrow(N15AA_mets_meta), function(i){
  data.frame(row = i, expand.grid(Metabolite = strsplit(N15AA_mets_meta$compound[i], split = '/')[[1]], KEGG = strsplit(N15AA_mets_meta$KEGG[i], split = '/')[[1]]))
})
N15AA_mets_meta_expanded <- do.call("rbind", N15AA_mets_meta_expanded)

N15AA_mets_meta_expanded <- N15AA_mets_meta_expanded %>% left_join(listTID) 
N15AA_mets_meta_expanded %>% View()

N15AA_mets_meta_expanded <- N15AA_mets_meta_expanded %>% select(row, Metabolite, KEGG, CHEBI, SpeciesName, SpeciesType)


if(!all(table(N15AA_mets_meta_expanded$SpeciesType) == 1)){
  degenerate_compounds <- names(which(table(N15AA_mets_meta_expanded$SpeciesType) > 1))
  best_matched_compounds <- NULL
  for(a_cmpd in degenerate_compounds){
    # find the best match based upon string similarity
    sim_subset <- N15AA_mets_meta_expanded %>% filter(SpeciesType == a_cmpd) %>% mutate(stringSim = levenshteinDist(Metabolite, SpeciesName)) %>%
      filter(stringSim == min(stringSim)) %>% select(-stringSim)
    if(nrow(sim_subset) != 1){
      # take the peak that is uniquely identified in case one compound is a mixture on a subset of instruments
      sim_subset <- sim_subset %>% group_by(boerRow) %>% mutate(nShared = length(boerMeta_annotated$boerRow[boerMeta_annotated$boerRow == boerRow])) %>%
        group_by() %>% filter(nShared == min(nShared)) %>% select(-nShared)
      if(nrow(sim_subset) != 1){
        warning("Multiple metabolite measurements must be combined before this point to generate a dataset level consensus")
        }
      }
    best_matched_compounds <- rbind(best_matched_compounds, sim_subset)
  }
  
  N15AA_mets_meta_expanded <- N15AA_mets_meta_expanded %>% filter(!(SpeciesType %in% degenerate_compounds)) %>%
    rbind(best_matched_compounds)
  
}

##### Compare model matched metabolites for each dataset #####
# Validate the model matches the correct metabolite to each dataset

# data.frames to be compared
metabolite_datasets_meta <- c("boerMeta_annotated", "yifan_meta_correspondence", "N15AA_mets_meta_expanded")

metabolite_dataset_summary <- lapply(metabolite_datasets_meta, function(x){
  get(x) %>% select(Metabolite, SpeciesName, SpeciesType) %>% mutate(dataset = x)
})
metabolite_dataset_summary <- do.call("rbind", metabolite_dataset_summary)

metabolite_dataset_summary_table <- dcast(metabolite_dataset_summary, SpeciesName ~ dataset, value.var = "Metabolite")
View(metabolite_dataset_summary_table)

tIDtoMet <- metabolite_dataset_summary %>% filter(dataset == "boerMeta_annotated") 
if(any(table(tIDtoMet$SpeciesType) != 1)){
 stop("Find a unique mapping between tIDs and metabolites")
}

##### Match experimental designs of each dataset #####
# so that a consensus relative abundance and conversion to absolute abundance can be made

met_cond_match <- data.frame(standard = c("P", "C", "N", "L", "U"), boer = c("PO4", "GLU", "NH4", "LEU", "URA"))

chemostatInfo <- read.table("../chemostatInformation/ListOfChemostats_augmented.txt", sep = "\t", header = T)
chemostatInfo <- chemostatInfo[chemostatInfo$include_this,]  

# pulldown all dataset conditions

metabolite_datasets <- c("boerSummary", "AA_absolute_quant", "yifanConc")

metabolite_conditions <- lapply(metabolite_datasets, function(x){
  get(x) %>% select(condition) %>% unique() %>% mutate(dataset = x)
})
metabolite_conditions <- do.call("rbind", metabolite_conditions)

# Populate dataset specific limitation and DR
metabolite_conditions <- rbind(
  # Boer
  metabolite_conditions %>% filter(dataset == "boerSummary") %>% mutate(boer_Limitation = substr(condition, 1, 3)) %>%
  group_by(boer_Limitation) %>% mutate(Limitation = met_cond_match$standard[met_cond_match$boer == unique(boer_Limitation)]) %>%
  rowwise() %>% mutate(actualDR = as.numeric(substr(condition, 5, 20))) %>% select(-boer_Limitation),
  # N15
  metabolite_conditions %>% filter(dataset == "AA_absolute_quant") %>% rowwise() %>% mutate(Limitation = "P", actualDR = AA_absolute_quant$DR[AA_absolute_quant$condition == condition][1]),
  # YifanAlign
  metabolite_conditions %>% filter(dataset == "yifanConc") %>% mutate(Limitation = "C", actualDR = 0.30)
)

##### Use metabolite abundances to inform a common set of conditions #####

metabolomicsMatrix <- boerSummary %>% left_join(metabolite_conditions %>% filter(dataset == "boerSummary")) %>%
  mutate(condition = paste0(Limitation, actualDR)) %>% acast(formula = "compound ~ condition", value.var = "log2_RA")
metabolomicsSD <- boerSummary %>% left_join(metabolite_conditions %>% filter(dataset == "boerSummary")) %>%
  mutate(condition = paste0(Limitation, actualDR)) %>% acast(formula = "compound ~ condition", value.var = "log2_CV")

if(!all(colnames(metabolomicsMatrix) == colnames(metabolomicsSD))){stop("Metabolite relative abundances and variances don't match")}

n_c <- ncol(metabolomicsMatrix)
n_m <- nrow(metabolomicsMatrix)
  
### Determine how many significant principal components exist in the metabolomics matrix ###
  
matrix_svd <- svd(metabolomicsMatrix)
plot(matrix_svd$d^2/sum(matrix_svd$d^2)) #scree-plot - fraction of variance explained by each PC
  
### determine how many significant principal components should be included based on repeated random sub-sampling validation ###
pcrange <- c(2,18)
npc.compare <- estim_ncpPCA(metabolomicsMatrix, ncp.min = pcrange[1], ncp.max = pcrange[2], method.cv = 'Kfold', pNA = 0.10, nbsim = 100)

# take the value of # PCs which minimizes the 
npc.cons <- (pcrange[1]:pcrange[2])[npc.compare$criterion < (max(npc.compare$criterion) - min(npc.compare$criterion))*0.3 + min(npc.compare$criterion)][1]
npc.min <- as.numeric(names(which.min(npc.compare$criterion)))

### metabolomic PC summary ###
  
pdf(file = "metaboliteSummaries/metPCnum.pdf", height = 6, width = 6)
plot(npc.compare$criterion ~ c(pcrange[1]:pcrange[2]), pch = 16, ylab = "MS error of prediction", xlab = "number of PCs", main = "Optimal number of metabolomic principal components")
abline(v = npc.cons, col = "RED", lwd = 2)
abline(v = npc.min, col = "BLUE", lwd = 2)
plot((matrix_svd$d)^2 / sum((matrix_svd$d)^2) ~ c(1:length(matrix_svd$d)), pch = 16, cex = 2, col = "RED", xlab = "PC Number", ylab = "Fraction of variance explained")
dev.off()
  
metPCs <- matrix_svd$v[,1:npc.min]; rownames(metPCs) <- colnames(metabolomicsMatrix); colnames(metPCs) <- paste("PC", c(1:npc.min))
pc_plot_df <- melt(metPCs)
colnames(pc_plot_df) <- c("condition", "PC", "value")
pc_plot_df$cond <- factor(sapply(as.character(pc_plot_df$condition), function(x){unlist(strsplit(x, ""))[1]}))
pc_plot_df$PC <- factor(pc_plot_df$PC, levels = paste("PC", c(1:npc.min)))

factor_plot <- ggplot(pc_plot_df, aes(x = condition, y = value, group = cond, col = PC)) + facet_wrap(~ PC, ncol = 2, scales = "free_y") + scale_x_discrete("Experimental condition") + scale_y_continuous("Principal Component Value") + theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = "aliceblue"), strip.background = element_rect(fill = "cadetblue1"), text = element_text(size = 15), axis.text.x = element_text(angle = 90), title = element_text(size = 25, face = "bold"))
factor_plot + geom_line() + ggtitle("Metabolomic principal components")
ggsave(file = "metaboliteSummaries/metPCs.pdf", height = 10, width = 12)

### summary of SVD of metabolomics matrix ###
reorgMets <- metabolomicsMatrix[hclust(d = dist(metabolomicsMatrix))$order,]
reorgMetSVD <- svd(reorgMets, nu = npc.min, nv = npc.min)
  
pdf("metaboliteSummaries/PCsummary.pdf", height = 10, width = 10)
  
heatmap.2(reorgMets, Rowv = F, Colv = F, trace = "none", symkey = T, col = greenred(50), main = "Raw")
heatmap.2(reorgMetSVD$u, Rowv = F, Colv = F, trace = "none", symkey = T, col = greenred(50), main = "U: Principal Component Loadings")
heatmap.2(diag(reorgMetSVD$d[1:npc]), Rowv = F, Colv = F, trace = "none", symkey = T, col = greenred(50), main = "D: Principal Component Eigenvalues")
heatmap.2(t(reorgMetSVD$v), Rowv = F, Colv = F, trace = "none", symkey = T, col = greenred(50), main = "t(V): Principal Components")
  
svd_projection <- reorgMetSVD$u %*% diag(reorgMetSVD$d[1:npc.min]) %*% t(reorgMetSVD$v)
  
heatmap.2(svd_projection, Rowv = F, Colv = F, trace = "none", symkey = T, col = greenred(50), main = paste("UDt(V) - ", npc, " dimensional summary", sep =""))
heatmap.2(reorgMets - svd_projection, trace = "none", symkey = T, col = greenred(50), main = "Residual Variation - Reclustered")
  
dev.off()
  
##### Project boer data onto signficant principal components and then interpolate values at Lim/DR ####

metSVD <- svd(metabolomicsMatrix)
metMatrixProj <- metSVD$u[,1:npc.min] %*% diag(metSVD$d[1:npc.min]) %*% t(metSVD$v[,1:npc.min])
rownames(metMatrixProj) <- rownames(metabolomicsMatrix); colnames(metMatrixProj) <- colnames(metabolomicsMatrix)

goal_conditions <- chemostatInfo %>% select(ChemostatCond, Limitation, actualDR)
boer_conditions <- metabolite_conditions %>% filter(dataset == "boerSummary") %>% mutate(condition = paste0(Limitation, actualDR))
if(nrow(goal_conditions) != nrow(boer_conditions)){
  stop("Number of design conditions does not match the number of metabolomics conditions")
  }

# Rows are boer conditions, columns are goal conditions
DR_change_mat <- matrix(0, nrow = n_c, ncol = n_c)
for(cond in 1:n_c){
  #find the 2 closest DR within the same limitation
  c_match <- c(1:n_c)[boer_conditions$Limitation == goal_conditions$Limitation[cond]]
  flanking_match <- c_match[order(abs(boer_conditions[c_match,]$actualDR - goal_conditions$actualDR[cond]))[1:2]]
  
  lb_diff <- (goal_conditions$actualDR[cond] - boer_conditions$actualDR[flanking_match][1])/diff(boer_conditions$actualDR[flanking_match])
  DR_change_mat[flanking_match,cond] <- c((1-lb_diff), lb_diff)
}

if(!(all((boer_conditions$actualDR %*% DR_change_mat) - goal_conditions$actualDR < 10^-15))){
 stop("Matching between conditions seems to have failed")
}

remapped_metabolites <- metMatrixProj %*% DR_change_mat
remapped_SD <- as.matrix(metabolomicsSD) %*% DR_change_mat
colnames(remapped_metabolites) <- colnames(remapped_SD) <- goal_conditions$ChemostatCond

# SVD of remapped metabolite -> save so that the principal components of met relative abundance can be used
# using a slightly conservative number of principal components
metSVD <- svd(remapped_metabolites, nu = npc.cons, nv = npc.cons)
rownames(metSVD$v) <- colnames(remapped_metabolites)
  
# Convert boer dataset names to model names, rows will now correspond to boerMeta_annotated

# add slots of inferred metabolites to quantitative data so that dimensionality matches boerMeta

inferredMets <- matrix(NA, ncol = ncol(remapped_metabolites), nrow = 2)
rownames(inferredMets) <- c("phosphate", "phosphoenolpyruvate")

if(!all(rownames(remapped_metabolites) %in% boerMeta[!(boerMeta$Metabolite %in% rownames(inferredMets)),'Metabolite'])){
  stop("Boer metabolites cannot be matched to their meta data, check when these names were aligned at the beginning of this script")
  }

remapped_metabolites <- remapped_metabolites[chmatch(boerMeta[!(boerMeta$Metabolite %in% rownames(inferredMets)),'Metabolite'], rownames(remapped_metabolites)),]
remapped_SD <- remapped_SD[chmatch(boerMeta[!(boerMeta$Metabolite %in% rownames(inferredMets)),'Metabolite'], rownames(remapped_SD)),]

remapped_metabolites <- rbind(remapped_metabolites, inferredMets)
remapped_SD <- rbind(remapped_SD, inferredMets)

remapped_corr <- shrunkCorr[chmatch(boerMeta[!(boerMeta$Metabolite %in% rownames(inferredMets)),'Metabolite'], rownames(shrunkCorr)),
                            chmatch(boerMeta[!(boerMeta$Metabolite %in% rownames(inferredMets)),'Metabolite'], colnames(shrunkCorr))]
expanded_met_correlations <- matrix(NA, ncol = nrow(boerMeta), nrow = nrow(boerMeta))
rownames(expanded_met_correlations) <- colnames(expanded_met_correlations) <- boerMeta$Metabolite
  
expanded_met_correlations[chmatch(rownames(shrunkCorr), boerMeta$Metabolite), chmatch(colnames(shrunkCorr), boerMeta$Metabolite)] <- shrunkCorr
  

if(!(all(rownames(remapped_metabolites) == boerMeta[,'Metabolite']))){
  stop("Names were misaligned")
  }

# Now that names are matched, the model conversions between the boer meta data and the model can be used
remapped_metabolites <- remapped_metabolites[boerMeta_annotated$boerRow,]
remapped_SD <- remapped_SD[boerMeta_annotated$boerRow,]
expanded_met_correlations <- expanded_met_correlations[boerMeta_annotated$boerRow,boerMeta_annotated$boerRow]

rownames(remapped_metabolites) <- rownames(remapped_SD) <- boerMeta_annotated$SpeciesName
rownames(expanded_met_correlations) <- colnames(expanded_met_correlations) <- boerMeta_annotated$SpeciesName

#### Convert the relative to absolute abundances (where available) ####

absolute_quant_datasets <- data.frame(quant_data = c("AA_absolute_quant", "yifanConc"), meta_data = c("N15AA_mets_meta_expanded", "yifan_meta_correspondence"))
absolute_met_dataset_subset <- metabolite_dataset_summary_table[!is.na(metabolite_dataset_summary_table %>% select(boerMeta_annotated)) & rowSums(!is.na(metabolite_dataset_summary_table[,colnames(metabolite_dataset_summary_table) %in% absolute_quant_datasets$meta_data])) != 0,]

# find C = [absolute met]/[2^relative met]
# so 2^relative met * C = [absolute met]

absolute_conditions <- metabolite_conditions %>% filter(dataset %in% absolute_quant_datasets$quant_data)

DR_change_mat <- matrix(0, nrow = n_c, ncol = nrow(absolute_conditions))
for(cond in 1:nrow(absolute_conditions)){
  #find the 2 closest DR within the same limitation
  c_match <- c(1:n_c)[goal_conditions$Limitation == absolute_conditions$Limitation[cond]]
  flanking_match <- c_match[order(abs(goal_conditions[c_match,]$actualDR - absolute_conditions$actualDR[cond]))[1:2]]
  
  lb_diff <- (absolute_conditions$actualDR[cond] - goal_conditions$actualDR[flanking_match][1])/diff(goal_conditions$actualDR[flanking_match])
  DR_change_mat[flanking_match,cond] <- c((1-lb_diff), lb_diff)
}

# determine metabolite relative abundance at conditions where absolute measurements were made
absolute_boer_converts <- remapped_metabolites %*% DR_change_mat
absolute_boer_converts <- absolute_boer_converts[rownames(absolute_boer_converts) %in% absolute_met_dataset_subset$SpeciesName,]
colnames(absolute_boer_converts) <- absolute_conditions$condition

absolute_abund_combined <- rbind(
  AA_absolute_quant %>% mutate(Metabolite = compound) %>% left_join(N15AA_mets_meta_expanded, by = "Metabolite") %>% 
    select(Metabolite, SpeciesName, condition, concentration) %>% mutate(dataset = "N15AA_mets_meta_expanded"),
  yifanConc[yifan_meta_correspondence$row,] %>% mutate(Metabolite = compound) %>% left_join(yifan_meta_correspondence) %>% 
    select(Metabolite, SpeciesName, condition, concentration) %>% mutate(dataset = "yifan_meta_correspondence")
)

absolute_rel_comp <- absolute_abund_combined %>% left_join(melt(absolute_boer_converts) %>% select(SpeciesName = Var1, condition = Var2, relConc = value))

absolute_rel_comp <- absolute_rel_comp %>% group_by(SpeciesName) %>% filter(!(dataset == "yifan_meta_correspondence" & SpeciesName %in% unique(absolute_abund_combined$SpeciesName[absolute_abund_combined$dataset == "N15AA_mets_meta_expanded"]))) %>%
  summarize(concConv = mean(concentration / 2^relConc)) %>% filter(!is.na(concConv))


##### Infer concentrations for a subset of species #####

# Recalculate the phosphate concentration with the absolute ATP/ADP values

atp <- 2^remapped_metabolites[rownames(remapped_metabolites) == "ATP",]*absolute_rel_comp$concConv[absolute_rel_comp$SpeciesName == "ATP"]
adp <- 2^remapped_metabolites[rownames(remapped_metabolites) == "ADP",]*absolute_rel_comp$concConv[absolute_rel_comp$SpeciesName == "ADP"]
h2o <- 55 # (assumption pure water)

dG0 = 37.9 # kJ / mol,eQuilibrator
dG = 57 # kJ/mol Bionumbers http://bionumbers.hms.harvard.edu//bionumber.aspx?id=100775&ver=0

R = 8.3144621*10^(-3) # kJ/(mol*K)
Tm= 310.15 # K

p = (atp*h2o)/adp * exp((dG0-dG)/(R*Tm))

# calculate the standard deviation of the phosphate concentration

peqtn <- "log(((2^ATP)*h2o)/(2^ADP) * 2.71828^((dG0-dG)/(R*Tm)))/log(2)"
eq <- eval(parse(text = paste('expression(',peqtn,')')))

derivList <- list()
for(spec in c('ATP', 'ADP')){
  derivList[[spec]] <- D(eq, spec)
}

for(cond in 1:n_c){
  standard_devs <- remapped_SD[chmatch(c('ATP', 'ADP'), rownames(remapped_SD)), cond]
  species_corr <- shrunkCorr[chmatch(c('ATP', 'ADP'), rownames(shrunkCorr)), chmatch(c('ATP', 'ADP'), rownames(shrunkCorr))]
  species_cov <- species_corr*(standard_devs %*% t(rep(1, 2)))*(rep(1, 2) %*% t(standard_devs))
  cond_values <- remapped_metabolites[chmatch(c('ATP', 'ADP'), rownames(remapped_metabolites)),cond]; names(cond_values) <- c('ATP', 'ADP')
  cond_values <- data.frame(ATP = cond_values[1], ADP = cond_values[2], h2o, dG0, dG, R, Tm)
  
  partial_deriv <- rep(NA, 2)
  for(n_der in 1:length(derivList)){
    partial_deriv[n_der] <- with(cond_values , eval(derivList[[n_der]]))
  }
  
  remapped_SD[rownames(remapped_SD) == "phosphate",cond] <- sqrt(t(partial_deriv) %*% species_cov %*% partial_deriv) # save the compute log-space standard deviation
}

remapped_metabolites[rownames(remapped_metabolites) == "phosphate",] <- log2(p) - log2(p)[names(p) == "P0.05"]
absolute_rel_comp <- rbind(absolute_rel_comp, data.frame(SpeciesName = "phosphate", concConv = p['P0.05']/2^remapped_metabolites[rownames(remapped_metabolites) == "phosphate",'P0.05']))

# add 3pg as phosphoenolpyruvate

remapped_metabolites[rownames(remapped_metabolites) == "phosphoenolpyruvate",] <- remapped_metabolites[rownames(remapped_metabolites) == "3-phosphoglycerate",]
remapped_SD[rownames(remapped_SD) == "phosphoenolpyruvate",] <- remapped_SD[rownames(remapped_SD) == "3-phosphoglycerate",]

  
# determine the correlation of added components
expanded_met_correlations[,rownames(expanded_met_correlations) == "phosphoenolpyruvate"] <- expanded_met_correlations[rownames(expanded_met_correlations) == "phosphoenolpyruvate",] <- expanded_met_correlations[rownames(expanded_met_correlations) == "3-phosphoglycerate",]
expanded_met_correlations[,rownames(expanded_met_correlations) == "phosphate"] <- expanded_met_correlations[rownames(expanded_met_correlations) == "phosphate",] <- 0 # Im not sure how to propagate the correlations of ATP and ADP through here
expanded_met_correlations[rownames(expanded_met_correlations) %in% c("phosphoenolpyruvate", "phosphate"),colnames(expanded_met_correlations) %in% c("phosphoenolpyruvate", "phosphate")] <- diag(2)
  
#### Save outputs ####

# Metabolite relative abundances at chemostat conditions with model name and ID

tab_boer <- boerMeta_annotated %>% select(SpeciesName, SpeciesType) %>% left_join(absolute_rel_comp) %>% cbind(remapped_metabolites)
write.table(tab_boer, '../../Yeast_genome_scale/flux_cache/tab_boer.txt', sep='\t', row.names = F, col.names = T, quote = F)

# Additional metabolite information - Coefficient of variation, residual correlations, and SVD of original metabolite matrix

save(tab_boer,metSVD,remapped_SD,expanded_met_correlations,file='../../Yeast_genome_scale/flux_cache/metaboliteTables.RData')

