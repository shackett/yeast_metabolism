library(reshape2)
library(nlme)
library(data.table)

setwd("/Users/Sean/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Metabolites/boerquant")
options(stringsAsFactors = FALSE)

tab_allBoer <- read.delim('./boer_data_2.txt',header=T,sep='\t')
tab_allBoerNames <- read.delim('./boer_data_2.txt',sep='\t', header = F, nrows = 1)
boerHeatmap <- read.delim("../BoerMetabolites.txt")[-1,-c(2:3)]

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

heatmap.2(as.matrix(log2(ntab_allBoer[,8:113])), trace = "none")

boerRAmat <- as.matrix(log2(ntab_allBoer[ntab_allBoer$Exp.Ref == "exp",8:113]))
colnames(boerRAmat) <- tab_allBoerNames[8:113]
boerRAmat[boerRAmat == 0] <- NA

### determine whether the experimental sample was measured (> 300 IC)
boerExpMissing = tab_allBoer[,8:113][ntab_allBoer$Exp.Ref == "exp",] == 300


boerSampleInfo <- ntab_allBoer[ntab_allBoer$Exp.Ref == "exp",1:7]
boerSampleInfo$Condition <- paste(boerSampleInfo$Nutr, boerSampleInfo$Gr, sep = ".")


plot(apply(boerRAmat, 1, median, na.rm = T) ~ factor(boerSampleInfo$Condition))
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


###### determine the correlation of all residuals across metabolites #####
residualStackDF <- acast(residualStacks, formula = metabolite ~ sampleNum, value.var = "residual")
residualCorr <- cor(t(residualStackDF), use = "pairwise.complete.obs")

heatmap.2(residualCorr, trace = "none") # correlation of residuals

library(corpcor)
library(impute) #impute missing values using knn imputation

residualStackDF <- impute.knn(residualStackDF, rowmax = 0.7)$data

library(colorRamps)
shrunkCorr <- cor.shrink(t(residualStackDF)) # shrink residual correlation matrix towards identity
pdf("boerCorr.pdf", height = 10, width = 10)
heatmap.2(shrunkCorr, trace = "none", symbreaks = T, col = blue2yellow(100)) 
dev.off()

partialCorrs <- cor2pcor(shrunkCorr) # partial correlation matrix
diag(partialCorrs) <- NA
rownames(partialCorrs) <- colnames(partialCorrs) <- rownames(shrunkCorr)
pdf("boerPartialcorr.pdf", height = 10, width = 10)
heatmap.2(partialCorrs, trace = "none", symbreaks = T, col = blue2yellow(100))
dev.off()
####

library(impute) #impute missing values using knn imputation
metRA = impute.knn(metRA)$data

for(a_row in 1:nrow(metSD)){ # use metabolite specific standard deviation as the sd for non-determined conditions
  metSD[a_row, is.na(metSD[a_row,])] <- SDbyMet$medianScale[a_row]
  }

write.table(metRA, "boerMean.tsv", quote = F, col.names = T, row.names = T, sep = "\t")
write.table(metSD, "boerSD.tsv", quote = F, col.names = T, row.names = T, sep = "\t")
write.table(shrunkCorr, "boerCorr.tsv", quote = F, col.names = T, row.names = T, sep = "\t")



