#library("rNMR")

### Convert bruker files to UCSF format - select the highest level folder 
#cf()
### open UCSF files
#fs()
### open region of interest (roi) tab
#roi()
#write.table(roiTable, "RoiTable.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = F)

setwd("/Users/Sean/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Media")

options(stringsAsFactors = FALSE)
library(reshape2)
library(gplots)
library(ggplot2)
library(data.table)
library(colorRamps)

##### Analysis #######

#### Describe the NMR samples ####

NMRsummary <- read.delim("effluentQuant_roiHeight.txt")

NMRsample_table <- data.frame(file = NMRsummary$File, date = NA, condition = NA, SampleType = NA, Rep = NA, is_sample = NA, limitation = NA, DR = NA)
for(i in 1:nrow(NMRsample_table)){
  dateSplit <- '[0-9]{4}_[0-9]{2}'
  NMRsample_table$date[i] <- regmatches(NMRsample_table$file[i], regexpr(dateSplit, NMRsample_table$file[i]))

  dateSplit <- '[0-9]{4}_[0-9]{2}_[0-9]{2,}'
  longMatch <- regmatches(NMRsample_table$file[i], regexpr(dateSplit, NMRsample_table$file[i]))
  if(length(longMatch) != 0){
    NMRsample_table$date[i] <- longMatch
    }
  
  condName <- strsplit(NMRsample_table$file[i], split = NMRsample_table$date[i])[[1]][2]
  NMRsample_table$condition[i] <- strsplit(condName, split = '^_|.ucsf')[[1]][2]

  }

NMRsample_table$condition <- unname(sapply(NMRsample_table$condition, function(x){sub('ref_', 'p0.05H', x)}))
NMRsample_table$Rep <- sapply(NMRsample_table$condition, function(x){regmatches(x, regexpr('[0-9A-Z]$', x))})
NMRsample_table$SampleType <- sapply(NMRsample_table$condition, function(x){strsplit(x, '_*[0-9A-Z]$', x)[[1]][1]})
NMRsample_table$SampleType <- toupper(NMRsample_table$SampleType)

NMRsample_table$is_sample <- TRUE; NMRsample_table$is_sample[grep('STD', NMRsample_table$SampleType)] <- FALSE

NMRsample_table$limitation[NMRsample_table$is_sample] <- sapply(NMRsample_table$SampleType[NMRsample_table$is_sample], function(x){strsplit(x, '[0-9.]+')[[1]][1]})
NMRsample_table$DR[NMRsample_table$is_sample] <- sapply(NMRsample_table$SampleType[NMRsample_table$is_sample], function(x){regmatches(x, regexpr('[0-9.]{4}', x))})


#### Create a matrix of NMR data from peak height ####

NMRmatrix <- as.matrix(NMRsummary[,-1])
rownames(NMRmatrix) <- NMRsample_table$condition

### use the NMR peak with the lowest coefficient of variation ###

for(a_compound in c('Glucose', 'Ethanol', 'Glycerol', 'Acetate', 'Lactate')){
  
  if(a_compound %in% 'Ethanol' & length(grep(a_compound, colnames(NMRmatrix))) != 1){ # if specified take the biggest peak
    colnames(NMRmatrix)[grep(a_compound, colnames(NMRmatrix))][which.max(apply(NMRmatrix[,grep(a_compound, colnames(NMRmatrix))], 2, mean))] <- a_compound
  } else if(length(grep(a_compound, colnames(NMRmatrix))) == 1){
    
    colnames(NMRmatrix)[grep(a_compound, colnames(NMRmatrix))] <- a_compound
    
  } else if(length(grep(a_compound, colnames(NMRmatrix))) > 1){
    
    subMat <- data.table(NMRsample_table[,c('date', 'SampleType')], NMRmatrix[,grep(a_compound, colnames(NMRmatrix))])
    subMatmelt <- data.table(melt(subMat, id.vars = c('date', 'SampleType')))
    
    peakCV <- subMatmelt[,list(CV = sd(value)/mean(value)), by = c('variable', 'date', 'SampleType')]
    peakAvgCV <- peakCV[, list(medCV = median(CV)), by = variable]
    
    colnames(NMRmatrix)[colnames(NMRmatrix) == peakAvgCV$variable[which.min(peakAvgCV$medCV)]] <- a_compound
  } else {
    print(paste(a_compound, "not found"))
    
  }
}

heatmap.2(cor(NMRmatrix[,colnames(NMRmatrix) != "DSS_0"]), symkey = T, col = blue2yellow(100))
heatmap.2(scale(log2(NMRmatrix[grep('std', rownames(NMRmatrix), invert = T), grep('DSS', colnames(NMRmatrix), invert = T)]), center = T, scale = F), trace = "none", symkey = T, col = blue2yellow(100))

#### Define NMR standards ####

NMRstandards <- NMRsample_table[!NMRsample_table$is_sample,]
NMRstandards$relative_conc = as.numeric(sub('([0-9.]{1,5})([xX])', '\\1', regmatches(NMRstandards$condition, regexpr('^([0-9.]{1,5})[xX]', NMRstandards$condition))))

# standards specific to each date
standard_conc <- rbind(data.frame(date = "2013_04_27", compound = c("Glucose", "Ethanol", "Acetate", "Glycerol"), concentration = c(122, 244, 10, 10)), 
                       data.frame(date = "2013_12", compound = c("Glucose", "Ethanol", "Acetate", "Glycerol", "Lactate"), concentration = c(122, 244, 10, 10, 10))) #concentrations of 1x standard (mM)



### Regression of standard peak areas onto known concentration ####
# mapping standard abundance to concentration
# Confirm that concentrations are linear on a log-log plot

NMRstandardConc <- NMRmatrix[chmatch(NMRstandards$condition, rownames(NMRmatrix)), colnames(NMRmatrix) %in% standard_conc$compound]

LMlist <- list()

for(i in 1:ncol(NMRstandardConc)){
  
  specieConc <- data.table(standard = rownames(NMRstandardConc), stdConc = NMRstandards$relative_conc, peakHeight = NMRstandardConc[,i])
  dayConc <- standard_conc[standard_conc$compound == colnames(NMRstandardConc)[i],]
  specieConc$MasterConc <- dayConc$concentration[chmatch(NMRstandards$date[chmatch(specieConc$standard, NMRstandards$condition)], dayConc$date)]
  specieConc <- specieConc[!is.na(specieConc$MasterConc),] # remove standards which do not contain the relevent metabolite
  specieConc[,standardConc := stdConc * MasterConc] 
  
  log_peakSize = log2(specieConc$peakHeight)
  log_Conc = log2(specieConc$standardConc)
  
  concLM <- lm(log_Conc ~ log_peakSize)
  LMlist[[colnames(NMRstandardConc)[i]]] <- concLM
  
  print(plot(log_Conc ~ log_peakSize, xlab = "log_peak", ylab = "log_conc", main = colnames(NMRstandardConc)[i], pch = 16))
  abline(a = concLM$coef[1], b = concLM$coef[2], col = "red")
  print(plot(concLM, which = 1))
  
}

### Quantifying identified metabolite concentrations in each experimental sample

NMRsample_data <- NMRmatrix[NMRsample_table$is_sample, grep('DSS', colnames(NMRmatrix), invert = T)]
NMRsamples <- NMRsample_table[NMRsample_table$is_sample,]

NMRsample_conc <- NMRsample_data; NMRsample_conc[!is.na(NMRsample_conc)] <- NA

for(i in 1:ncol(NMRsample_data)){
  if(colnames(NMRsample_data)[i] %in% standard_conc$compound){
    ### scale peak height to standards ###
    NMRsample_conc[,i] <- 2^predict(LMlist[[colnames(NMRsample_data)[i]]], data.frame(log_peakSize = log2(NMRsample_data[,i])))
    
    }else{
      ### scale peak height to DSS molarity (5mM*(1/9)) (although proton number of a peak affects this)
      NMRsample_conc[,i] <- sapply(NMRsample_data[,i], function(x){max(x, 0)}) * 5/9
      }
  }


#### Deconvolve overlapping peaks of Ethanol-Lactate and Lactate-Glycerol using standards w & w/o lactate ####

# A # from standards peakArea(concentration)
# B # from pure peak, quantify concentration
# C # for mixed peak subtract E[conc] from peak area

NMRarea <- read.delim("effluentQuant_roiArea.txt")

if(!all(NMRarea$File == NMRsummary$File)){
  print("roiSummary (height) and roiArea do not match")
  }

NMRarea <- as.matrix(NMRarea[,-1])
rownames(NMRarea) <- NMRsample_table$condition

NMRstandardArea <- NMRarea[chmatch(NMRstandards$condition, rownames(NMRmatrix)),]

### Ethanol-Lactate ###

specieConc <- data.table(standard = rownames(NMRstandardConc), stdConc = NMRstandards$relative_conc, peakHeight = NMRstandardArea[,'Lactate_13'])
specieConc$Lactate <- sapply(NMRstandards$date, function(x){conc = standard_conc$conc[standard_conc$compound == "Lactate" & standard_conc$date == x]; ifelse(length(conc) == 1, conc, 0)})*specieConc$stdConc
specieConc$Ethanol <- sapply(NMRstandards$date, function(x){conc = standard_conc$conc[standard_conc$compound == "Ethanol" & standard_conc$date == x]; ifelse(length(conc) == 1, conc, 0)})*specieConc$stdConc

plot(specieConc$peakHeight ~ specieConc$stdConc, col = factor(specieConc$Lactate))

LEfit <- lm(data = specieConc, formula = peakHeight ~ Lactate + Ethanol + 0)
lactate_conc <- (NMRarea[NMRsample_table$is_sample,'Lactate_13'] - (NMRsample_conc[,c('Ethanol')] * LEfit$coef[names(LEfit$coef) == "Ethanol"]))/LEfit$coef[names(LEfit$coef) == "Lactate"]

NMRsample_conc[,c('Lactate')] <- lactate_conc

### Lactate - Unknown ###

#specieConc <- data.table(standard = rownames(NMRstandardConc), stdConc = NMRstandards$relative_conc, peakHeight = NMRstandardArea[,'Lactate_41'])
#specieConc$Lactate <- sapply(NMRstandards$date, function(x){conc = standard_conc$conc[standard_conc$compound == "Lactate" & standard_conc$date == x]; ifelse(length(conc) == 1, conc, 0)})*specieConc$stdConc
#specieConc$Glycerol <- sapply(NMRstandards$date, function(x){conc = standard_conc$conc[standard_conc$compound == "Glycerol" & standard_conc$date == x]; ifelse(length(conc) == 1, conc, 0)})*specieConc$stdConc

#plot(specieConc$peakHeight ~ specieConc$stdConc, col = factor(specieConc$Lactate))

#LGfit <- lm(data = specieConc, formula = peakHeight ~ Lactate + Glycerol + 0)
#glycerol_conc <- (NMRarea[NMRsample_table$is_sample,'Lactate_41'] - (NMRsample_conc[,c('lactate_inferred')] * LGfit$coef[names(LGfit$coef) == "Lactate"]))/LGfit$coef[names(LGfit$coef) == "Glycerol"]

#NMRsample_conc <- cbind(NMRsample_conc, data.frame(glycerol_inferred = glycerol_conc))

#plot(NMRsample_conc$Glucose ~ NMRsample_conc$glycerol_inferred)


### Remove some outlier measurements (ethanol contamination) and samples ###
# determine the studentized-residuals for the linear model: concentration ~ condition + slope*condition + residual

NMR_resids <- NMRsample_conc; NMR_resids[!is.na(NMR_resids)] <- NA

condInfo <- NMRsamples[,c('limitation', 'DR')]
condInfo[,'limitation'] <- factor(condInfo[,'limitation'])
condInfo[,'DR'] <- as.numeric(condInfo[,'DR'])

regModel <- model.matrix(~ limitation + limitation*DR, condInfo)
condInfo$leverage <- diag(regModel %*% solve(t(regModel) %*% regModel) %*% t(regModel))

for(a_specie_n in 1:ncol(NMRsample_conc)){
  
  specie_lm <- lm(NMRsample_conc[,a_specie_n] ~ regModel + 0)
  NMR_resids[,a_specie_n] <- specie_lm$resid / (summary(specie_lm)$sigma * sqrt(1 - condInfo$leverage))
  
}

sampleDeviation <- apply(NMR_resids, 1, function(x){sum(abs(x))}) 

NMRsample_conc[NMRsamples$condition == 'L0.16_2',] <- NA # broadly inconsistent
NMRsample_conc[NMRsamples$condition == 'u0.30_1' & NMRsamples$date == "2013_12",] <- NA # noted turbidity
NMRsample_conc[NMRsamples$condition == 'c0.30_2','Ethanol'] <- NA # more C in ethanol measured than in provided glucose
NMRsample_conc[NMRsamples$condition == 'u0.05_1','Ethanol'] <- NA
NMRsample_conc[NMRsamples$condition == 'u0.30_1' & NMRsamples$date == "2013_04_27",'Ethanol'] <- NA
NMRsample_conc[NMRsamples$condition == 'u0.30_2' & NMRsamples$date == "2013_12",'Glucose'] <- NA


#### Summarize concentration by chemostat and generate summary figures of knowns + unknowns and just knowns ####


sampleAug <- data.frame(NMRsamples[,c('SampleType', 'date', 'limitation', 'DR')], NMRsample_conc)
sampleMelt <- data.table(melt(sampleAug, id.vars = c('SampleType', 'date', 'limitation', 'DR'), variable.name = "peak", value.name = 'estimate'))

NMR_point_estimate <- sampleMelt[,list(estimate = mean(estimate, na.rm = T), se = sd(estimate, na.rm = T)/sqrt(length(estimate[!is.na(estimate)]))), by = c("SampleType", "peak", "date", "limitation", "DR")]
NMR_point_estimate[,lb := estimate - 1.96*se]
NMR_point_estimate[,ub := estimate + 1.96*se]

setnames(NMR_point_estimate, c('SampleType'), c('condition'))




barplot_theme <- theme(text = element_text(size = 60, face = "bold"), title = element_text(size = 50, face = "bold"), panel.background = element_blank(), legend.position = "none", 
                       panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(), axis.text = element_text(size = 30, color = "black"), axis.text.x = element_text(angle = 90)) 

#NMR_point_estimate  <- NMR_point_estimate[NMR_point_estimate$peak %in% c("Ethanol", "Acetate", "Glucose", "Glycerol"),]



NMR_point_estimate_plot <- NMR_point_estimate
NMR_point_estimate_plot$peak <- factor(NMR_point_estimate_plot$peak, levels = c("Glucose", "Ethanol", "Acetate", "Lactate", "Glycerol", as.character(sort(unique(NMR_point_estimate_plot$peak)[!(unique(NMR_point_estimate_plot$peak) %in%  c("Glucose", "Ethanol", "Acetate", "Lactate", "Glycerol"))]))))
NMR_point_estimate_plot$limitation <- factor(NMR_point_estimate_plot$limitation, levels = c("P", "C", "N", "L", "U"))

setkeyv(NMR_point_estimate_plot, c('limitation', 'condition', 'date', 'peak'))

NMR_point_estimate_plot[order(limitation, condition, date, peak)]
NMR_point_estimate_plot[,cond_date := paste(condition, date, sep = "_")]
NMR_point_estimate_plot$cond_date <- factor(NMR_point_estimate_plot$cond_date, levels = unique(NMR_point_estimate_plot$cond_date))

unique_conditions <- NMR_point_estimate_plot[,list(condition = unique(condition)), by = cond_date]


NMRbarplot <- ggplot(NMR_point_estimate_plot, aes(x = cond_date, y = estimate, fill = factor(limitation))) + facet_wrap( ~ peak, scale = "free_y", ncol = 2) + barplot_theme
NMRbarplot + geom_bar(stat = "identity") + ggtitle('Concentration of metabolites (and unknowns) \n in chemostat effluent') + scale_fill_brewer(palette = "Set2") + 
  geom_errorbar(aes(ymin = lb, ymax = ub)) + scale_x_discrete("Chemostat condition", labels = NMR_point_estimate_plot$condition[chmatch(levels(NMR_point_estimate_plot$cond_date), as.character(NMR_point_estimate_plot$cond_date))]) + scale_y_continuous("Concentration (mM)")

ggsave("mediaComposition_all.pdf", width = 24, height = 22)

NMRbarplot <- ggplot(NMR_point_estimate_plot[NMR_point_estimate_plot$peak %in% standard_conc$compound,], aes(x = cond_date, y = estimate, fill = factor(limitation))) + facet_grid(peak ~ ., scale = "free_y") + barplot_theme
NMRbarplot + geom_bar(stat = "identity") + ggtitle("Concentration of metabolites in chemostat effluent") + scale_fill_brewer(palette = "Set2") + 
  geom_errorbar(aes(ymin = lb, ymax = ub)) + scale_x_discrete("Chemostat condition", labels = NMR_point_estimate_plot$condition[chmatch(levels(NMR_point_estimate_plot$cond_date), as.character(NMR_point_estimate_plot$cond_date))]) + scale_y_continuous("Concentration (mM)")

ggsave("mediaComposition.pdf", width = 14, height = 22)


write.table(NMR_point_estimate, "mediaComposition_NMR.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)



