setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/BulkComposition")
options(stringsAsFactors = FALSE)
library(reshape)
library(ggplot2)

chemoStat <- read.table("releventDR.csv", sep = ",", header = TRUE)
chemoStat <- chemoStat[apply(is.na(chemoStat), 1, sum) < 2,]


#plot(chemoStat$colterVol.uL.mL ~ chemoStat$PCMVol.mL)
#the coulter counter under-estimates volume by ~20%
coulter_underest <- median(chemoStat$colterVol.uL.mL/chemoStat$PCMVol.mL, na.rm = TRUE)

#determine whether the coulter counter, packed cell volume (PCV) or klett have lower relative error

klettTab <- chemoStat[!is.na(chemoStat$Klett),]
coultTab <- chemoStat[!is.na(chemoStat$colterVol.uL.mL),]
pcvTab <- chemoStat[!is.na(chemoStat$PCMVol.mL),]

anova(lm(klettTab$Klett ~ factor(klettTab$Chemostat)))
anova(lm(coultTab$colterVol.uL.mL ~ factor(coultTab$Chemostat)))
anova(lm(pcvTab$PCMVol.mL ~ factor(pcvTab$Chemostat)))

#klett has the greatest signal to noise and should be used as the normalization factor, but it needs to be scaled using the coulter-counter to give cellular volumes

lm(data = chemoStat, formula = colterVol.uL.mL ~ Klett)$coef

coultReg <- lm(data = chemoStat, formula = colterVol.uL.mL ~ Klett)
pcvReg <- lm(data = chemoStat, formula = PCMVol.mL ~ Klett)

coult_pred <- predict(coultReg, newdata = chemoStat)/coulter_underest
pcv_pred <- predict(pcvReg, newdata = chemoStat)

plot(chemoStat$colterVol.uL.mL ~ chemoStat$Klett); lines(coultReg$coef[1] + coultReg$coef[2]*seq(0, 200, by = 1) ~ seq(0, 200, by = 1))
plot(chemoStat$PCMVol.mL ~ chemoStat$Klett); lines(pcvReg$coef[1] + pcvReg$coef[2]*seq(0, 200, by = 1) ~ seq(0, 200, by = 1))

#### Klett has the lowest relative error, coulter counter has the best linear relationship with klett.  Klett will be scaled to coulter abundance via a lm and then this intracellular volume will be adjusted to account for the underestimation of volume by coulter relative to the packed cell volume

chemo.cols <- c("Chemostat", "DR", "coulterMed")
chemoStat <- chemoStat[,colnames(chemoStat) %in% chemo.cols]
chemoStat$cellularVolFrac = coult_pred


conditions <- data.frame(condition = rep(NA, times = 27), limitation = c(rep(c("p", "c", "n", "L", "u"), each = 5), "p", "p"), DRgoal = c(rep(c("0.05", "0.11", "0.16", "0.22", "0.30"), times = 5), c(0.05, 0.05)), actualDR = rep(NA, times = 27), VolFrac_mean = rep(NA, times = 27), VolFrac_SD = rep(NA, times = 27), medcellVol = rep(NA, times = 27))
conditions[,2] <- as.character(conditions[,2])
conditions[,3] <- as.character(conditions[,3])
conditions$condition <- apply(conditions[,c(2,3)], 1, paste, collapse = "")
conditions$condition[c(26,27)] <- c("ctrl1_0.05", "ctrl2_0.05")

for(i in 1:length(conditions[,1])){
	
	condSubset <- chemoStat[chemoStat$Chemostat == conditions$condition[i],]
	
	conditions$VolFrac_mean[i] <- mean(condSubset$cellularVolFrac[!is.na(condSubset$cellularVolFrac)])
	conditions$VolFrac_SD[i] <- sd(condSubset$cellularVolFrac[!is.na(condSubset$cellularVolFrac)])
	conditions$actualDR[i] <- mean(condSubset$DR[!is.na(condSubset$DR)])
	conditions$medcellVol[i] <- mean(condSubset$coulterMed)
	
	}
conditions$condition[c(26,27)] <- c("p0.05H1", "p0.05H2")

#### Import absorbance data and fit protein concentration ###

kineticDat <- read.delim("chemo_absoluteProtquant.txt", sep = "\t", header = FALSE)
kineticDat <- kineticDat[,-c(1,14)]

#the plate was measured 6 times, so split these into sepearte data.frames
plate_split <- list()
for(i in 1:6){
	plate_split[[i]] <- kineticDat[c(2:9) + 9*(i-1),]
	}
#collapse measurement replicates into a single plate
absorb_plate <- matrix(NA, ncol = 12, nrow = 8)
absorb_plate_sd <- matrix(NA, ncol = 12, nrow = 8)
for(i in 1:8){
for(j in 1:12){
	absorb_plate[i,j] <- mean(sapply(1:6, function(x){plate_split[[x]][i,j]}))
	absorb_plate_sd[i,j] <- sd(sapply(1:6, function(x){plate_split[[x]][i,j]}))
	}}
#absorb_plate_sd/absorb_plate: indicates that the means of all samples should be fine

standards <- data.frame(std.row = rep(c(4:8), each = 3), std.col = rep(c(10:12), times = 5), conc = rep(c(2000, 1500, 1000, 500, 0), each = 3), abs = NA)

for(i in 1:length(standards[,1])){
	standards$abs[i] <- absorb_plate[standards$std.row[i], standards$std.col[i]]
	}

#fit of absorption to protein concentration in ug/mL
fit_all <- lm(standards$conc ~ standards$abs)
fit_lt2k <- lm(standards$conc[standards$conc != 2000] ~ standards$abs[standards$conc != 2000])

plot(standards$conc ~ standards$abs, pch = 16, col = "RED", xlab = "absorption", ylab = "protein conc (ug/mL)")
lines(fit_all$coef[1] + fit_all$coef[2]*seq(from = 0, to = 1.5, by = 0.1) ~ seq(from = 0, to = 1.5, by = 0.1), col = "ORANGE", lwd = 2)
lines(fit_lt2k$coef[1] + fit_lt2k$coef[2]*seq(from = 0, to = 1.5, by = 0.1) ~ seq(from = 0, to = 1.5, by = 0.1), col = "PINK", lwd = 2)

standards$fitted_conc <- fit_lt2k$coef[1] + fit_lt2k$coef[2]*standards$abs
standards <- standards[!(standards$conc %in% 2000),]


#####

wellInfo <- read.table("wellInfo.txt", sep = "\t", header = TRUE)
wellInfo$abs <- NA; wellInfo$conc <- NA
strsplit(wellInfo$Well, "")

wellPos <- t(sapply(wellInfo$Well, function(x){c(unlist(strsplit(x, ""))[1], unlist(strsplit(x, "[A-Z]"))[2])}))
wellPos <- cbind(cell.row = sapply(wellPos[,1], function(x){c(1:8)[LETTERS[1:8] == x]}), cell.col = as.numeric(wellPos[,2]))

for(i in 1:length(wellInfo[,1])){
	wellInfo$abs[i] <- absorb_plate[wellPos[i,1], wellPos[i,2]]
	wellInfo$conc[i] <- fit_lt2k$coef[1] + fit_lt2k$coef[2]*wellInfo$abs[i]
	}

###### Remove bad wells ####
#negative concentration
wellInfo <- wellInfo[wellInfo$conc > 0,]
#pre-recorecorded underloaded
wellInfo <- wellInfo[!(wellInfo$Well %in% c("B04", "E09", "B01")),]

plot(wellInfo$conc ~ factor(wellInfo$Sample), col = "RED", xlab = "absorption", ylab = "protein conc (ug/mL)")

conditions$assayConc_mean <- NA; conditions$assayConc_SD <- NA

for(i in 1:length(conditions[,1])){
	conditions$assayConc_mean[i] <- mean(wellInfo$conc[wellInfo$Sample %in% conditions$condition[i]])
	conditions$assayConc_SD[i] <- sd(wellInfo$conc[wellInfo$Sample %in% conditions$condition[i]])
	}



sampleInfo <- read.table("sampleInfo.txt", sep = "\t", header = TRUE)
sampleInfo <- sampleInfo[sapply(sampleInfo$Sample, function(x){c(1:length(sampleInfo[,1]))[conditions$condition %in% x]}),]

#plot(conditions$assayConc_mean/sampleInfo$homog_weight)

#ug in assay
#conditions$assayConc_mean * 9/1000
#mg in original plate
#conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000
#% Dry-weight
#conditions$protein_DW_frac <- (conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)/sampleInfo$homog_weight
#mg protein per ml culture
#(conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)*(sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV
#concentration mg protein per uL cellular volume - check
#(conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)*(sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV * conditions$VolFrac_mean
#grams protein per 100mL cellular volume
#(conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)*(sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV * conditions$VolFrac_mean/1000 * 10^5


#use shrinkage towards fitted values to improve point-estimate

wellInfo$limitation <- NA; wellInfo$actualDR <- NA; wellInfo$cellVol <- NA; wellInfo$DWcorrPermL <- NA; wellInfo$DWperCultmL <- NA
for(well in 1:length(wellInfo[,1])){
	
	wellInfo$limitation[well] <- conditions$limitation[conditions$condition == wellInfo$Sample[well]]
	wellInfo$actualDR[well] <- conditions$actualDR[conditions$condition == wellInfo$Sample[well]]
	#correction for fraction of dried material homogenized and culture volume dried
	wellInfo$DWcorrPermL[well] <- ((sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV)[sampleInfo$Sample == wellInfo$Sample[well]]
	#coulter counter-measured uL of cellular volume per mL culture
	wellInfo$cellVol[well] <- conditions$VolFrac_mean[conditions$condition == wellInfo$Sample[well]]
	#mg dry-weight per mL culture
	wellInfo$DWperCultmL[well] <- (sampleInfo$DryWeight/sampleInfo$cultureV)[sampleInfo$Sample == wellInfo$Sample[well]]
	}
# 9 uL loaded * 200/9 fraction of dilution plate loaded * 900/10 fraction of homogenate in dilution plate - gives mg protein dry weight
wellInfo$protConc <- wellInfo$conc * 9/1000 * 200/9 * 90 * 1/1000 * wellInfo$DWcorrPermL * 1/wellInfo$cellVol

plot(wellInfo$protConc ~ factor(wellInfo$Sample))
plot(log(wellInfo$protConc) ~ factor(wellInfo$Sample))

#fitting a linear model predicting protein concentrations using limitation and limitation*DR
conclm1 <- lm(data = wellInfo, formula = protConc ~ factor(limitation) + actualDR + factor(limitation)*actualDR)
plot(conclm1, which = 2)
plot(conclm1, which = 3)
#violates the assumption of homoscedasticity

#fitting a linear model predicting log(protein) concentrations using limitation and limitation*DR
conclm2 <- lm(data = wellInfo, formula = log(protConc) ~ factor(limitation) + actualDR + factor(limitation)*actualDR)
plot(conclm2, which = 2)
plot(conclm2, which = 3)
#conforms to the assumption of equivariance gaussian residuals fairly well

#fitting a glm treating protConc as gamma and using the "inverse" link
conclm3 <- glm(data = wellInfo, formula = protConc ~ factor(limitation) + actualDR + factor(limitation)*actualDR, family = "Gamma")
plot(conclm3, which = 2)
plot(conclm3, which = 3)



refittedconc <- matrix(wellInfo$protConc, ncol = length(wellInfo[,1]))
SS_track <- Inf
SS_change <- Inf
it_continue <- TRUE
while(it_continue){
	
	fit_linMod <- lm(data = wellInfo, formula = refittedconc[length(refittedconc[,1]),] ~ factor(limitation) + actualDR + factor(limitation)*actualDR)
	
	MSE <- sum((fit_linMod$fitted - wellInfo$protConc)^2)/(length(wellInfo[,1]) - 10)
	
	withinCondVar <- sapply(conditions$condition, function(x){
		var(wellInfo$protConc[wellInfo$Sample == x])
		})
	
	#shrink by within-treatment MSE relative to average RMSE error
	shrinkFrac <- sapply(wellInfo$Sample, function(x){withinCondVar[conditions$condition == x]})/(sapply(wellInfo$Sample, function(x){withinCondVar[conditions$condition == x]}) + MSE)
	
	updated_conc <- fit_linMod$fitted*shrinkFrac + (1-shrinkFrac)*wellInfo$protConc
	refittedconc <- rbind(refittedconc, updated_conc)
	
	SS_change <- SS_track[length(SS_track)] - sum(fit_linMod$resid^2)
	SS_track <- c(SS_track, sum(fit_linMod$resid^2))
	if(SS_change < 10^-6){
		it_continue <- FALSE
		}else{
			it_continue <- TRUE
			}
	
	
	}

plot(refittedconc[length(refittedconc[,1]),] ~ factor(wellInfo$Sample)) 
plot(refittedconc[1,] ~ factor(wellInfo$Sample))


rownames(refittedconc) <- paste("It", 0:(length(refittedconc[,1]) - 1)); colnames(refittedconc) <- wellInfo$Sample
#choosing which iterations to display: intial, 1/3, 2/3 and final
if(length(refittedconc[,1]) > 4){
	col_det <- sapply(1 + c((length(refittedconc[,1])-2)/3, 2*(length(refittedconc[,1])-2)/3), function(x){
		c(2:(length(refittedconc[,1])-1))[which.min(abs(c(2:(length(refittedconc[,1])-1)) - x))]
		})
	refittedconc <- refittedconc[c(1, col_det, length(refittedconc[,1])),]
	}
refittedconc <- refittedconc[c(1, length(refittedconc[,1])),]

plotting_DF <- melt(refittedconc)
levels(plotting_DF[,1]) <- levels(plotting_DF[,1])[order(as.numeric(sapply(levels(plotting_DF[,1]), function(x){unlist(strsplit(x, split = " "))[2]})))]
colnames(plotting_DF) <- c("Iteration", "Condition", "Concentration")
plotting_DF$limitation <- sapply(as.character(plotting_DF$Condition), function(cond){
	wellInfo$limitation[wellInfo$Sample == cond][1]
	})

pconc_plot <- ggplot(plotting_DF, aes(x = factor(Condition), y = Concentration, col = limitation)) + facet_wrap(~ Iteration, ncol = 1) + theme(axis.text.x = element_text(size = 4, face = "bold")) + scale_y_continuous("mg protein/uL cellular volume") + scale_x_discrete("Experimental condition")
pconc_plot + geom_boxplot() 

conditions$prot_conc <- sapply(conditions$condition, function(cond){
	mean(refittedconc[length(refittedconc[,1]),][names(refittedconc[length(refittedconc[,1]),]) == cond])
	})

conditions$logSF <- log2(conditions$prot_conc) - mean(log2(conditions$prot_conc)[conditions$condition %in% c("p0.05H1", "p0.05H2")])
conditions$logSF[conditions$condition %in% c("p0.05H1", "p0.05H2")] <- NA


#backcalculating fraction of dry-weight from shrunk protein concentration values
dry_weight_frac <- refittedconc[length(refittedconc[,1]),] * wellInfo$cellVol/wellInfo$DWperCultmL
dry_weight_DF <- melt(dry_weight_frac); colnames(dry_weight_DF) <- "ProteinDW_fraction"; dry_weight_DF$condition <- names(dry_weight_frac); dry_weight_DF$limitation <-  sapply(as.character(dry_weight_DF$condition), function(cond){wellInfo$limitation[wellInfo$Sample == cond][1]})

dw_plot <- ggplot(dry_weight_DF, aes(x = factor(condition), y = ProteinDW_fraction, col = limitation)) + theme(axis.text.x = element_text(size = 4, face = "bold")) + scale_y_continuous("Protein fraction of dry-material") + scale_x_discrete("Experimental condition")
dw_plot + geom_boxplot()

#look at the yield in terms of dry-weight per culture volume

#mass per mL
dry_weight_turnover_DF <- data.frame(Condition = conditions$condition, DryWeight_conc = (sampleInfo$DryWeight/sampleInfo$cultureV), Protein = sapply(conditions$condition, function(cond){mean(((refittedconc[length(refittedconc[,1]),] * wellInfo$cellVol))[wellInfo$Sample == cond])}), DR = conditions$actualDR)

dry_weight_turnover_DF$non_Protein <- dry_weight_turnover_DF$DryWeight_conc - dry_weight_turnover_DF$Protein
dry_weight_turnover_DFmelt <- melt(dry_weight_turnover_DF, id.vars = c("DR", "Condition"))

#make a stacked barplot # this is the amount of material in a chemostat, ignoring dif
dry_weight_turnover_DFmelt <- dry_weight_turnover_DFmelt[!(dry_weight_turnover_DFmelt$variable == "DryWeight_conc"),]
dry_weight_turnover_plot1 <- ggplot(dry_weight_turnover_DFmelt, aes(x = factor(Condition), y = value, fill = variable)) + theme(axis.text.x = element_text(size = 4, face = "bold")) + scale_y_continuous("Dry weight(mg) per mL culture", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(palette = "Set1")
dry_weight_turnover_plot1 + geom_bar(stat = "identity", position = "stack")

#dry weight conc * DR - optimal growth, in terms of biomass production at intermediate GRs
dry_weight_turnover_DFmelt2 <- dry_weight_turnover_DFmelt; dry_weight_turnover_DFmelt2$value <- dry_weight_turnover_DFmelt$value * dry_weight_turnover_DFmelt$DR
dry_weight_turnover_plot2 <- ggplot(dry_weight_turnover_DFmelt2, aes(x = factor(Condition), y = value, fill = variable)) + theme(axis.text.x = element_text(size = 4, face = "bold")) + scale_y_continuous("Dry weight production (mg/hr) per mL culture", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(palette = "Set1")
dry_weight_turnover_plot2 + geom_bar(stat = "identity", position = "stack")


####
#divide by cellVolume (uL) per mL media to get mg protein/dry-matter per uL of cellular volume
dry_weight_turnover_DF_cellular <- dry_weight_turnover_DF
dry_weight_turnover_DF_cellular$DryWeight_conc <- dry_weight_turnover_DF_cellular$DryWeight_conc/conditions$VolFrac_mean
dry_weight_turnover_DF_cellular$Protein <- dry_weight_turnover_DF_cellular$Protein/conditions$VolFrac_mean
dry_weight_turnover_DF_cellular$non_Protein <- dry_weight_turnover_DF_cellular$non_Protein/conditions$VolFrac_mean

dry_weight_turnover_D_cellular_DFmelt <- melt(dry_weight_turnover_DF_cellular, id.vars = c("DR", "Condition"))
dry_weight_turnover_D_cellular_DFmelt <- dry_weight_turnover_D_cellular_DFmelt[!(dry_weight_turnover_D_cellular_DFmelt$variable == "DryWeight_conc"),]
dry_weight_turnover_plot3 <- ggplot(dry_weight_turnover_D_cellular_DFmelt, aes(x = factor(Condition), y = value, fill = variable)) + theme(axis.text.x = element_text(size = 4, face = "bold")) + scale_y_continuous("Dry weight (mg) per uL of cellular volume", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(palette = "Set1")
dry_weight_turnover_plot3 + geom_bar(stat = "identity", position = "stack")

dry_weight_perCell_DF <- dry_weight_turnover_DF_cellular
dry_weight_perCell_DF$cellVolume <- conditions$medcellVol
dry_weight_perCell_DFmelt <- melt(dry_weight_perCell_DF, id.vars = c("DR", "Condition", "cellVolume"))
dry_weight_perCell_DFmelt$value <- dry_weight_perCell_DFmelt$value * dry_weight_perCell_DFmelt$cellVolume
dry_weight_perCell_DFmelt <- dry_weight_perCell_DFmelt[!(dry_weight_perCell_DFmelt$variable == "DryWeight_conc"),]
dry_weight_perCell_plot <- ggplot(dry_weight_perCell_DFmelt, aes(x = factor(Condition), y = value, fill = variable)) + theme(axis.text.x = element_text(size = 4, face = "bold"), panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + scale_y_continuous("pg protein/dry-material per cell", expand = c(0,0)) + scale_x_discrete("Experimental condition") + scale_fill_brewer(palette = "Set1")
dry_weight_perCell_plot + geom_bar(stat = "identity", position = "stack")





save(conditions, plotting_DF, dry_weight_DF, dry_weight_turnover_DFmelt, dry_weight_turnover_DFmelt2, dry_weight_turnover_D_cellular_DFmelt, dry_weight_perCell_DFmelt, file = "protSpecQuantOut.Rdata")