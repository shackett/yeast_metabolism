setwd("/Users/Sean/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Metabolites/absolute_quant_compare")
options(stringsAsFactors = F)

### setup exactive data and formatting ####

Exactive <- read.table("exactiveSummary.csv", sep = ",", header = F)
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

### setup max data and formatting ####

Max <- read.table("Max_2013_05_03-highestPeaks.csv", sep = ",", header = F)
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

### Combine max and exactive data into a single matrix - point estimation of chemostat logIC and regression of dilution on batch culture to determine ratio of chemostat to full concentration batch ####

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
# chemostat difference relative to batch samples, log2



#####
yifanConc <- read.delim("../yeast_absolute_concentration_yifan.txt")

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

write.table(yifanConc[!is.na(yifanConc$c_lim_conc),], file = "../yeast_absolute_concentration_chemo.txt", sep = "\t", quote = F, col.names = TRUE, row.names = F)

#metScaling$compound[!(metScaling$compound %in% yifanConc$Compound)]
#yifanConc$Compound[!(yifanConc$Compound %in% metScaling$compound)]
