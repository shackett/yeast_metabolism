setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Media")

nutrientUptakeRates <- list()

#### import DR and 

chemostatInfo <- read.table("../BulkComposition/chemostatDRmisc.tsv", sep = "\t", header = TRUE)

nutrientFile <- read.delim("../../Yeast_genome_scale/Boer_nutrients.txt")[1:6,1:6]

highp <- max(nutrientFile[nutrientFile[,1] == "phosphate",][-1]) * 10^-6 * 10^9 #20ul * 1/20 dilution = 1ul of media per assay #nmoles max per assay
lowp <- min(nutrientFile[nutrientFile[,1] == "phosphate",][-1]) * 0.0002 * 10^9 #200ul per assay

################################################################

#### Import absorbance data and fit phosphate concentration for non-phosphate limited conditions ###

kineticDat <- read.delim("phosphateSpec/phosphateLow.txt", sep = "\t", header = FALSE)
plateIDs <- unique(kineticDat[,14])[unique(kineticDat[,14]) != ""]
kineticDat <- kineticDat[,-c(1,14)]

#the plate was measured 6 times, so split these into sepearte data.frames
plate_split <- list()
for(i in 1:length(plateIDs)){
  kinData <- kineticDat[c(2:9) + 9*(i-1),]
  kinData[kinData < 0 | kinData == "?????"] <- NA
  plate_split[[i]] <- kinData
}
#collapse measurement replicates into a single plate
PLcorrect <- FALSE#use corrected or not
if(PLcorrect){plateIndices <- grep('Corrected', plateIDs)}else{plateIndices <- grep('Corrected', plateIDs, invert = TRUE)}

absorb_plate <- matrix(NA, ncol = 12, nrow = 8)
absorb_plate_sd <- matrix(NA, ncol = 12, nrow = 8)
for(i in 1:8){
  for(j in 1:12){
    absorbMeasure <- as.numeric(sapply(plateIndices, function(x){plate_split[[x]][i,j]}))
    absorb_plate[i,j] <- mean(absorbMeasure, na.rm = TRUE)
    absorb_plate_sd[i,j] <- sd(absorbMeasure, na.rm = TRUE)
  }}
#absorb_plate_sd/absorb_plate: indicates that the means of all samples should be fine

#check the order
experimental_samples <- data.frame(condition = c(rep(rep(c("p", "n", "c", "u", "L"), each = 5), times = 2), rep("p", 4)), DR = c(rep(c("0.05", "0.11", "0.16", "0.22", "0.30"), times = 10), rep(c("0.05H1", "0.05H2"), times =2)))

absorbDF <- melt(absorb_plate)
colnames(absorbDF) <- c("row", "column", "absorbance")
absorbDF$conc <- NA
absorbDF$conc[absorbDF$column == 8] <- absorbDF$conc[absorbDF$column == 9] <- c(5:0, 0, NA)

absorbDF$limitation <- c(experimental_samples$condition, rep(NA, length(absorbDF[,1]) - length(experimental_samples[,1])))
absorbDF$DR <- c(experimental_samples$DR, rep(NA, length(absorbDF[,1]) - length(experimental_samples[,1])))
absorbDF <- absorbDF[!is.na(absorbDF$limitation) | !is.na(absorbDF$conc),]

conc_fit <- lm(data = absorbDF[!is.na(absorbDF$conc),], conc ~ absorbance)$coef

absorbDF$conc[is.na(absorbDF$conc)] <- conc_fit[1] + conc_fit[2]*absorbDF$absorbance[is.na(absorbDF$conc)]
plot(absorbDF[!is.na(absorbDF$conc),]$conc ~ absorbDF[!is.na(absorbDF$conc),]$absorbance)

absorbDF$collapsedName <- apply(absorbDF, 1, function(namepaste){paste(namepaste[colnames(absorbDF) %in% c("limitation", "DR")], collapse = "")})

plot(absorbDF$conc[!is.na(absorbDF$limitation) & !(absorbDF$limitation %in% c("p", "Pctrl"))] ~ factor(absorbDF$collapsedName[!is.na(absorbDF$limitation) & !(absorbDF$limitation %in% c("p", "Pctrl"))])) 
abline(h = highp, col = "RED") #this is the amount of phosphate in the original media

relChemos <- chemostatInfo[chemostatInfo$condition %in% absorbDF$collapsedName,]
relChemos <- relChemos[unlist(sapply(absorbDF$collapsedName, function(x){c(1:length(relChemos[,1]))[relChemos$condition == x]})),]
relChemos <- relChemos[relChemos$limitation != "p",]

conditionDF <- relChemos[!is.na(absorbDF[absorbDF$collapsedName %in% relChemos$condition,]$absorbance),]
measuredVal <- absorbDF[absorbDF$collapsedName %in% relChemos$condition,][!is.na(absorbDF[absorbDF$collapsedName %in% relChemos$condition,]$absorbance),]

shrunkConcEst <- shrinkageRegression(conditionDF, measuredVal)
plot(shrunkConcEst[1,] ~ factor(relChemos$condition)) 
plot(shrunkConcEst[length(shrunkConcEst[,1]),] ~ factor(relChemos$condition)) 
abline(h = highp, col = "RED") #this is the amount of phosphate in the original media
measuredVal$shrunkConc = shrunkConcEst[length(shrunkConcEst[,1]),]

repeatedCondDF <- measuredVal[measuredVal$collapsedName %in% names(table(measuredVal$collapsedName))[table(measuredVal$collapsedName) >= 2],]
#commonSD <- sqrt(sum((repeatedCondDF$conc - sapply(repeatedCondDF$collapsedName, function(x){mean(repeatedCondDF$conc[repeatedCondDF$collapsedName == x])}))^2)/(length(repeatedCondDF[,1]) - 1))

phosphateUptakeNL = data.frame(condition = unique(measuredVal$collapsedName), phosphateUptakeShrunk = highp - sapply(unique(measuredVal$collapsedName), function(x){mean(measuredVal$shrunkConc[measuredVal$collapsedName == x])}), 
                               phosphateUptakeRaw = highp - sapply(unique(measuredVal$collapsedName), function(x){mean(measuredVal$conc[measuredVal$collapsedName == x])}),
                               phosphateSD = sapply(unique(measuredVal$collapsedName), function(x){sd(measuredVal$conc[measuredVal$collapsedName == x])}))
#phosphateUptakeNL$phosphateSD[is.na(phosphateUptakeNL$phosphateSD) | phosphateUptakeNL$phosphateSD < commonSD] <- commonSD
phosphateUptakeNL$ceiling <- highp
phosphateUptakeNL[,-1] <- phosphateUptakeNL[,-1]/(10^-6 * 10^9)



################################################################

#### Import absorbance data and fit phosphate concentration for non-phosphate limited conditions ##

kineticDat <- read.delim("phosphateSpec/phosphateHighestRedo.txt", sep = "\t", header = FALSE)
plateIDs <- unique(kineticDat[,14])[unique(kineticDat[,14]) != ""]
kineticDat <- kineticDat[,-c(1,14)]

#the plate was measured 6 times, so split these into sepearte data.frames
plate_split <- list()
for(i in 1:length(plateIDs)){
  kinData <- kineticDat[c(2:9) + 9*(i-1),]
  kinData[kinData < 0 | kinData == "?????"] <- NA
  plate_split[[i]] <- kinData
}
#collapse measurement replicates into a single plate
PLcorrect <- FALSE#use corrected or not
if(PLcorrect){plateIndices <- grep('Corrected', plateIDs)}else{plateIndices <- grep('Corrected', plateIDs, invert = TRUE)}

absorb_plate <- matrix(NA, ncol = 12, nrow = 8)
absorb_plate_sd <- matrix(NA, ncol = 12, nrow = 8)
for(i in 1:8){
  for(j in 1:12){
    absorbMeasure <- as.numeric(sapply(plateIndices, function(x){plate_split[[x]][i,j]}))
    absorb_plate[i,j] <- mean(absorbMeasure, na.rm = TRUE)
    absorb_plate_sd[i,j] <- sd(absorbMeasure, na.rm = TRUE)
  }}

absorbDFlow <- melt(absorb_plate)
colnames(absorbDFlow) <- c("row", "column", "absorbance")
absorbDFlow$conc <- NA
absorbDFlow$conc[absorbDFlow$column == 8] <- absorbDFlow$conc[absorbDFlow$column == 9] <- c(rep(NA, 5), 5:3)
absorbDFlow$conc[absorbDFlow$column == 10] <- absorbDFlow$conc[absorbDFlow$column == 11] <- c(2:0, 0, rep(NA, 4))

absorbDFlow$limitation <- c(experimental_samples$condition, rep(NA, length(absorbDFlow[,1]) - length(experimental_samples[,1])))
absorbDFlow$DR <- c(experimental_samples$DR, rep(NA, length(absorbDFlow[,1]) - length(experimental_samples[,1])))
absorbDFlow <- absorbDFlow[!is.na(absorbDFlow$limitation) | !is.na(absorbDFlow$conc),]
absorbDFlow <- absorbDFlow[is.na(absorbDFlow$limitation) | (absorbDFlow$limitation %in% c("p", "Pctrl")),]
  
absorbDFlow$collapsedName <- apply(absorbDFlow, 1, function(namepaste){paste(namepaste[colnames(absorbDFlow) %in% c("limitation", "DR")], collapse = '')})

plot(absorbDFlow$absorbance[!is.na(absorbDFlow$limitation)] ~ factor(absorbDFlow$collapsedName[!is.na(absorbDFlow$limitation)]))
#there are two obvious outliers (one of which was misfiltered (the elevated p0.05 sample)
absorbDFlow$absorbance[!is.na(absorbDFlow$limitation)][absorbDFlow$absorbance[!is.na(absorbDFlow$limitation)] > 0.4] <- NA
plot(absorbDFlow$absorbance[!is.na(absorbDFlow$limitation)] ~ factor(absorbDFlow$collapsedName[!is.na(absorbDFlow$limitation)]))
#strangely high density low DR chemostats have a moderate level of residual phosphate

conc_fit <- lm(data = absorbDFlow[!is.na(absorbDFlow$conc),], conc ~ absorbance)$coef

absorbDFlow$conc[is.na(absorbDFlow$conc)] <- conc_fit[1] + conc_fit[2]*absorbDFlow$absorbance[is.na(absorbDFlow$conc)]
plot(absorbDF[!is.na(absorbDF$conc),]$conc ~ absorbDF[!is.na(absorbDF$conc),]$absorbance)

relChemos <- chemostatInfo[chemostatInfo$condition %in% absorbDFlow$collapsedName,]
relChemos <- relChemos[unlist(sapply(absorbDFlow$collapsedName, function(x){c(1:length(relChemos[,1]))[relChemos$condition == x]})),]

plot(absorbDFlow$conc[!is.na(absorbDFlow$limitation)] ~ factor(absorbDFlow$collapsedName[!is.na(absorbDFlow$limitation)]))


conditionDF <- relChemos[!is.na(absorbDFlow[absorbDFlow$collapsedName %in% relChemos$condition,]$absorbance),]
measuredVal <- absorbDFlow[absorbDFlow$collapsedName %in% relChemos$condition,][!is.na(absorbDFlow[absorbDFlow$collapsedName %in% relChemos$condition,]$absorbance),]

#use a the higher of a condition-specific variance and a common variance to err on conservative estimation
repeatedCondDF <- measuredVal[measuredVal$collapsedName %in% names(table(measuredVal$collapsedName))[table(measuredVal$collapsedName) >= 2],]
commonSD <- sqrt(sum((repeatedCondDF$conc - sapply(repeatedCondDF$collapsedName, function(x){mean(repeatedCondDF$conc[repeatedCondDF$collapsedName == x])}))^2)/(length(repeatedCondDF[,1]) - 1))

phosphateUptake = data.frame(condition = unique(measuredVal$collapsedName), phosphateUptakeShrunk = NA, phosphateUptakeRaw = lowp - sapply(unique(measuredVal$collapsedName), function(x){mean(measuredVal$conc[measuredVal$collapsedName == x])}),
                             phosphateSD = sapply(unique(measuredVal$collapsedName), function(x){sd(measuredVal$conc[measuredVal$collapsedName == x])}))
phosphateUptake$phosphateSD[is.na(phosphateUptake$phosphateSD) | phosphateUptake$phosphateSD < commonSD] <- commonSD
phosphateUptake$ceiling <- lowp

phosphateUptake[,-1] <- phosphateUptake[,-1]/(0.0002 * 10^9)

nutrientUptakeRates$phosphate <- rbind(phosphateUptakeNL, phosphateUptake)


################################################################

#### Import absorbance data and fit ammonium concentration for all cultures (which have different dilution factors) ###
## We want to measure the change in absorbance associated with the enzymatic depletion of NADPH (in an ammount equal to NH4) ##
## Difference between pre-enzyme absorbance (preread) and final incubated absorbance (after 70 minutes)

## Preread ##

kineticDat <- read.delim("ammoniumSpec/ammonium_preread.txt", sep = "\t", header = FALSE)
plateIDs <- unique(kineticDat[,14])[unique(kineticDat[,14]) != ""]
kineticDat <- kineticDat[,-c(1,14)]

wellpl <- read.delim("ammoniumSpec/ammonium_pl.txt", sep = "\t", header = FALSE)
wellpl <- melt(wellpl[-1,-c(1,14)])[,2] #in cm

#the plate was measured 6 times, so split these into sepearte data.frames
plate_split <- list()
for(i in 1:length(plateIDs)){
  kinData <- kineticDat[c(2:9) + 9*(i-1),]
  kinData[kinData < 0 | kinData == "?????"] <- NA
  plate_split[[i]] <- kinData
}
#collapse measurement replicates into a single plate
PLcorrect <- TRUE#use corrected or not
if(PLcorrect){plateIndices <- grep('Corrected', plateIDs)}else{plateIndices <- grep('Corrected', plateIDs, invert = TRUE)}

absorb_plate_preread <- matrix(NA, ncol = 12, nrow = 8)
for(i in 1:8){
  for(j in 1:12){
    absorbMeasure <- as.numeric(sapply(plateIndices, function(x){plate_split[[x]][i,j]}))
    absorb_plate_preread[i,j] <- mean(absorbMeasure, na.rm = TRUE)
  }}

## post-incubation ###

kineticDat <- read.delim("ammoniumSpec/ammoniumConc30min.txt", sep = "\t", header = FALSE)
plateIDs <- unique(kineticDat[,14])[unique(kineticDat[,14]) != ""]
kineticDat <- kineticDat[,-c(1,14)]

#the plate was measured 6 times, so split these into sepearte data.frames
plate_split <- list()
for(i in 1:length(plateIDs)){
  kinData <- kineticDat[c(2:9) + 9*(i-1),]
  kinData[kinData < 0 | kinData == "?????"] <- NA
  plate_split[[i]] <- kinData
}
#collapse measurement replicates into a single plate
if(PLcorrect){plateIndices <- grep('Corrected', plateIDs)}else{plateIndices <- grep('Corrected', plateIDs, invert = TRUE)}

absorb_plate_inc <- matrix(NA, ncol = 12, nrow = 8)
for(i in 1:8){
  for(j in 1:12){
    absorbMeasure <- as.numeric(sapply(plateIndices, function(x){plate_split[[x]][i,j]}))
    absorb_plate_inc[i,j] <- mean(absorbMeasure, na.rm = TRUE)
  }}

### aggregate multiple reads of the same plate to try to measure the change in absorbance due to NADPH competitive reaction, not related to ammonium

kineticDat <- read.delim("ammoniumSpec/ammoniumConcTimecourse.txt", sep = "\t", header = FALSE)
plateIDs <- data.frame(identifier = 0:7, time = c(ISOdatetime(2013, 2, 25, hour = 15, min = 49, sec = 33), ISOdatetime(2013, 2, 25, 16, 16, 12), ISOdatetime(2013, 2, 25, 16, 20, 48), 
                                                  ISOdatetime(2013, 2, 25, 16, 42, 56), ISOdatetime(2013, 2, 25, 16, 50, 54), ISOdatetime(2013, 2, 25, 17, 12, 56), ISOdatetime(2013, 2, 25, 23, 56, 28), ISOdatetime(2013, 2, 26, 11, 04, 14)))
plateIDs$timeDiff <- plateIDs$time - plateIDs$time[1]
kineticDat <- kineticDat[,-c(1,14)]

#the plate was measured 6 times, so split these into sepearte data.frames
plate_split <- list()
for(i in 1:length(plateIDs[,1])){
  kinData <- kineticDat[c(2:9) + 9*(i-1),]
  kinData[kinData < 0 | kinData == "?????"] <- NA
  plate_split[[i]] <- kinData
}

#collapse measurement replicates into a single plate
if(PLcorrect){plateIndices <- grep('Corrected', plateIDs)}else{plateIndices <- grep('Corrected', plateIDs, invert = TRUE)}

absorb_timecourse <- matrix(NA, nrow = 96, ncol = length(plateIDs[,1]))
for(c in 1:length(plateIDs[,1])){
  absorb_timecourse[,c] <- as.numeric(melt(as.matrix(plate_split[[c]]))[,3])
}

### determine how much of the change in absorbance is from the competing reaction from change between 10 minutes and 50 minutes

concAdj <- NA;
concAdj[rowSums(is.na(absorb_timecourse)) == 0] <- apply(absorb_timecourse[rowSums(is.na(absorb_timecourse)) == 0,], 1, function(x){lm(x[2:5] ~ as.numeric(plateIDs$timeDiff)[2:5])$coef[2]}) * as.numeric(plateIDs$timeDiff[3])
      

image(absorb_plate_preread - absorb_plate_inc)

absorb_change <- melt(absorb_plate_preread - absorb_plate_inc)
colnames(absorb_change) <- c("row", "column", "absorbance")
absorb_change$absorbance <- absorb_change$absorbance + concAdj


#reagent difference starting at H6 results in large differences in absorbance (taking the difference between pre-read and post-incubation plates takes care of most of this effect
#but this will also be used as an adjustment variable in the difference
reagCO <- ifelse(PLcorrect, 1.2, 0.8)
absorb_change$reagent <- ifelse(melt(absorb_plate_preread)[,3] > reagCO, "B", "A")
absorb_change$absorbance[is.nan(absorb_change$absorbance)] <- NA

absorb_change$conc <- NA
absorb_change$conc[absorb_change$column == 8] <- absorb_change$conc[absorb_change$column == 9] <- c(40, 20, 10, 5, 0, 0, NA, NA)


absorb_change$limitation <- c(experimental_samples$condition, rep(NA, length(absorb_change[,1]) - length(experimental_samples[,1])))
absorb_change$DR <- c(experimental_samples$DR, rep(NA, length(absorb_change[,1]) - length(experimental_samples[,1])))
absorb_change$collapsedName <- apply(absorb_change, 1, function(namepaste){paste(namepaste[colnames(absorb_change) %in% c("limitation", "DR")], collapse = '')})
absorb_change$pl <- wellpl

absorb_timecourse <- absorb_timecourse[!is.na(absorb_change$limitation) | !is.na(absorb_change$conc),]
absorb_change <- absorb_change[!is.na(absorb_change$limitation) | !is.na(absorb_change$conc),]

plot(absorb_timecourse[1,]  ~ plateIDs$timeDiff)
       
absorb_change$microliterMedia <- NA;
absorb_change$microliterMedia[absorb_change$limitation %in% c("p", "c", "u", "L")] <- 10*(1/20)*(2/19)
absorb_change$microliterMedia[absorb_change$limitation %in% c("n")] <- 40

plot(absorb_change$absorbance[!is.na(absorb_change$limitation)] ~ factor(absorb_change$collapsedName[!is.na(absorb_change$limitation)]))


absorb_change$absorbanceReagentAdj <- absorb_change$absorbance - ifelse(absorb_change$reagent == "B", lm(absorb_change, formula = absorbance ~ factor(reagent) + factor(limitation)*factor(DR))$coef[2], 0) #regress out the effect of reagent B

conc_fit <- lm(data = absorb_change[!is.na(absorb_change$conc),], conc ~ absorbanceReagentAdj)$coef #predicting equivalent uL of standard from OD change
plot(absorb_change[!is.na(absorb_change$conc),]$conc ~ absorb_change[!is.na(absorb_change$conc),]$absorbanceReagentAdj)
abline(a = conc_fit[1], b = conc_fit[2], lwd = 2, col = "RED")

#confirm the concentration of ammonium sulphate standard
invreg <- lm(data = absorb_change[!is.na(absorb_change$conc),], absorbanceReagentAdj ~ conc)$coef
#invreg[2]/(6.22 * mean(absorb_change$pl[is.na(absorb_change$limitation)])) #gives mM NADPH per uL of standard
#invreg[2]/(6.22 * mean(absorb_change$pl[is.na(absorb_change$limitation)])) * 200 * 10^-6 #mmoles per uL of media
#stdconc <- invreg[2]/(6.22 * mean(absorb_change$pl[is.na(absorb_change$limitation)])) * 200 * 10^-6 * 10^6 * 1000 #umoles per L of media #stated concentration is 588 uM (not bad)
stdconc <- invreg[2]/(6.22) * 200 * 10^-6 * 10^6 * 1000 #umoles per L of media #stated concentration is 588 uM (not bad) #use if pl is accounted for


#concentration from molar absorptivity
(absorb_change$absorbance[is.na(absorb_change$conc)] - invreg[1]) * 0.2 * (40/absorb_change$microliterMedia[is.na(absorb_change$conc)]) / 0.04 * 0.00273 / 17
invreg[2] * 0.2 / 0.001 * 0.00273 * 1000 #ug/ml


absorb_change$absorbanceReagentAdj[is.na(absorb_change$conc)]/(6.22 * absorb_change$pl[is.na(absorb_change$conc)]) * 200 * 10^-6 # mmoles in well
absorb_change$absorbanceReagentAdj[is.na(absorb_change$conc)]/(6.22 * absorb_change$pl[is.na(absorb_change$conc)]) * 200 * 10^-6 / absorb_change$microliterMedia[is.na(absorb_change$conc)] #mmoles per uL media
absorb_change$absorbanceReagentAdj[is.na(absorb_change$conc)]/(6.22 * absorb_change$pl[is.na(absorb_change$conc)]) * 200 * 10^-6 / absorb_change$microliterMedia[is.na(absorb_change$conc)] * 10^3 #moles per L media


absorb_change$conc[is.na(absorb_change$conc)] <- conc_fit[1] + conc_fit[2]*absorb_change$absorbanceReagentAdj[is.na(absorb_change$conc)]
plot(absorb_change$conc[!is.na(absorb_change$limitation)] ~ factor(absorb_change$collapsedName[!is.na(absorb_change$limitation)]))

relChemos <- chemostatInfo[chemostatInfo$condition %in% absorb_change$collapsedName,]
relChemos <- relChemos[unlist(sapply(absorb_change$collapsedName, function(x){c(1:length(relChemos[,1]))[relChemos$condition == x]})),]

conditionDF <- relChemos[!is.na(absorb_change[absorb_change$collapsedName %in% relChemos$condition,]$absorbance),]
measuredVal <- absorb_change[absorb_change$collapsedName %in% relChemos$condition,][!is.na(absorb_change[absorbDF$collapsedName %in% relChemos$condition,]$absorbance),]

shrunkConcEst <- shrinkageRegression(conditionDF, measuredVal)
plot(shrunkConcEst[1,] ~ factor(relChemos$condition)) 
plot(shrunkConcEst[length(shrunkConcEst[,1]),] ~ factor(relChemos$condition)) 
measuredVal$shrunkConc <- shrunkConcEst[length(shrunkConcEst[,1]),]

#concentration is in uL

# 588 * 10^-12 #concentration in uL of 588 uM ammonium sulfate

measuredVal$ceiling <- ifelse(measuredVal$limitation == "n", nutrientFile[nutrientFile$X == "ammonium",-1]$Ammonium, nutrientFile[nutrientFile$X == "ammonium",-1]$Glucose)


measuredVal$conc * stdconc * 10^-12 * 1/measuredVal$microliterMedia * 10^6


measuredVal$shrunkConc * stdconc * 10^-12 * 1/measuredVal$microliterMedia * 10^6






nutrientFile[nutrientFile$X == "ammonium",-1]$Ammonium #N-limited

nutrientFile[nutrientFile$X == "ammonium",-1]$Glucose #not N-limited


ammoniumUptake = data.frame(condition = unique(measuredVal$collapsedName), ammoniumUptakeShrunk = highp - sapply(unique(measuredVal$collapsedName), function(x){mean(measuredVal$shrunkConc[measuredVal$collapsedName == x])}), 
                               ammoniumUptakeRaw = highp - sapply(unique(measuredVal$collapsedName), function(x){mean(measuredVal$conc[measuredVal$collapsedName == x])}),
                               phosphateSD = sapply(unique(measuredVal$collapsedName), function(x){sd(measuredVal$conc[measuredVal$collapsedName == x])}))



write.table(rbind(phosphateUptakeNL, phosphateUptake), file = "measuredNutrients.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)#phosphate concentration



######### Functions ##########

shrinkageRegression <- function(conditionDF, measuredVal){
  
  refittedconc <- matrix(measuredVal$conc, ncol = length(conditionDF[,1]))
  SS_track <- Inf
  SS_change <- Inf
  it_continue <- TRUE
  while(it_continue){
    
    if(length(unique(conditionDF$limitation)) == 1){
      regModel <- model.matrix(data = conditionDF, ~ actualDR)  
    }else{
      regModel <- model.matrix(data = conditionDF, ~ factor(limitation) + actualDR + factor(limitation)*actualDR)
    }
    
    fit_linMod <- lm(refittedconc[length(refittedconc[,1]),] ~ regModel + 0)
    
    MSE <- sum((fit_linMod$fitted - measuredVal$conc)^2)/(length(conditionDF[,1]) - length(regModel[1,]))
    
    withinCondVar <- sapply(conditionDF$condition, function(x){
      var(measuredVal$conc[measuredVal$collapsedName == x])
    })
    
    #shrink by within-treatment MSE relative to average RMSE error
    shrinkFrac <- sapply(measuredVal$collapsedName, function(x){withinCondVar[conditionDF$condition == x][1]})/(sapply(measuredVal$collapsedName, function(x){withinCondVar[conditionDF$condition == x][1]}) + MSE)
    
    updated_conc <- fit_linMod$fitted*shrinkFrac + (1-shrinkFrac)*measuredVal$conc
    updated_conc[is.na(updated_conc)] <- refittedconc[length(refittedconc[,1]),is.na(updated_conc)]
    refittedconc <- rbind(refittedconc, updated_conc)
    
    SS_change <- SS_track[length(SS_track)] - sum(fit_linMod$resid^2)
    SS_track <- c(SS_track, sum(fit_linMod$resid^2))
    if(SS_change < 10^-6){
      it_continue <- FALSE
    }else{
      it_continue <- TRUE
    }
  }
  refittedconc
}




