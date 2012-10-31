setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/BulkComposition")
options(stringsAsFactors = FALSE)
library(reshape)
library(ggplot2)

chemoStat <- read.table("releventDR.csv", sep = ",", header = TRUE)
chemo.cols <- c("Chemostat", "DR", "colterVol.uL.mL")
chemoStat <- chemoStat[,colnames(chemoStat) %in% chemo.cols]
chemoStat <- chemoStat[apply(is.na(chemoStat), 1, sum) < 2,]

conditions <- data.frame(condition = rep(NA, times = 27), limitation = c(rep(c("p", "c", "n", "L", "u"), each = 5), "p", "p"), DRgoal = c(rep(c("0.05", "0.11", "0.16", "0.22", "0.30"), times = 5), c(0.05, 0.05)), actualDR = rep(NA, times = 27), cellVol_mean = rep(NA, times = 27), cellVol_SD = rep(NA, times = 27))
conditions[,2] <- as.character(conditions[,2])
conditions[,3] <- as.character(conditions[,3])
conditions$condition <- apply(conditions[,c(2,3)], 1, paste, collapse = "")
conditions$condition[c(26,27)] <- c("ctrl1_0.05", "ctrl2_0.05")

#chemoStat$Chemostat[!(chemoStat$Chemostat %in% conditions$condition)]

for(i in 1:27){
	
	condSubset <- chemoStat[chemoStat$Chemostat == conditions$condition[i],]
	
	conditions$cellVol_mean[i] <- mean(condSubset$colterVol.uL.mL[!is.na(condSubset$colterVol.uL.mL)])
	conditions$cellVol_SD[i] <- sd(condSubset$colterVol.uL.mL[!is.na(condSubset$colterVol.uL.mL)])
	conditions$actualDR[i] <- mean(condSubset$DR[!is.na(condSubset$DR)])
	
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

plot(conditions$assayConc_mean/sampleInfo$homog_weight)



#ug in assay
conditions$assayConc_mean * 9/1000
#mg in original plate
conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000
#% Dry-weight
conditions$protein_DW_frac <- (conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)/sampleInfo$homog_weight
#mg protein per ml culture
(conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)*(sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV
#concentration mg protein per uL cellular volume - check
(conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)*(sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV * conditions$cellVol_mean
#grams protein per 100mL cellular volume
(conditions$assayConc_mean * 9/1000 * 200/9 * 90 * 1/1000)*(sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV * conditions$cellVol_mean/1000 * 10^5


#use shrinkage towards fitted values to improve point-estimate

conditions

wellInfo$limitation <- NA; wellInfo$actualDR <- NA; wellInfo$cellVol <- NA; wellInfo$DWperCultmL <- NA
for(well in 1:length(wellInfo[,1])){
	
	wellInfo$limitation[well] <- conditions$limitation[conditions$condition == wellInfo$Sample[well]]
	wellInfo$actualDR[well] <- conditions$actualDR[conditions$condition == wellInfo$Sample[well]]
	wellInfo$DWperCultmL[well] <- ((sampleInfo$DryWeight / sampleInfo$homog_weight)/sampleInfo$cultureV)[sampleInfo$Sample == wellInfo$Sample[well]]
	wellInfo$cellVol[well] <- conditions$cellVol_mean[conditions$condition == wellInfo$Sample[well]]
	
	}

wellInfo$protConc <- wellInfo$conc * 9/1000 * 200/9 * 90 * 1/1000 * wellInfo$DWperCultmL * wellInfo$cellVol

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

#fitting a glm treating protConc as gamma and using the "inverse" link
conclm4 <- glm(data = wellInfo, formula = protConc ~ factor(limitation) + actualDR + factor(limitation)*actualDR, family = gaussian(link = "log"))
plot(conclm3, which = 2)
plot(conclm3, which = 3)

rmse <- sqrt(sum(conclm2$resid^2)/(length(wellInfo[,1]) - 10))

shrinkFrac <- 1 - 1/(abs(conclm2$resid)^2*(length(wellInfo[,1])/(length(wellInfo[,1]) - 10))/(rmse^2))
shrinkFrac <- sapply(shrinkFrac, function(x){max(0, x)})




refittedconc <- matrix(wellInfo$protConc, ncol = length(wellInfo[,1]))
SS_track <- Inf#sum((log(wellInfo$protConc) - mean(log(wellInfo$protConc)))^2)
SS_change <- Inf
it_continue <- TRUE
while(it_continue){
	
	fit_linMod <- lm(data = wellInfo, formula = log(refittedconc[length(refittedconc[,1]),]) ~ factor(limitation) + actualDR + factor(limitation)*actualDR)
	#fit_linMod <- glm(data = wellInfo, formula = refittedconc[length(refittedconc[,1]),] ~ factor(limitation) + actualDR + factor(limitation)*actualDR, family = gaussian(link = "log"))

	#rmse <- sqrt(sum(fit_linMod$resid^2)/(length(wellInfo[,1]) - 10))
	rmse <- sqrt(sum((fit_linMod$fitted - log(wellInfo$protConc))^2)/(length(wellInfo[,1]) - 10))
	
	
	#MS logCond - mean(logCond) #should this be about the mean or the fitted value
	MSlog <- sapply(conditions$condition, function(x){
		sum((log(wellInfo$protConc[wellInfo$Sample == x]) - mean(log(wellInfo$protConc[wellInfo$Sample == x])))^2)/sum(wellInfo$Sample == x)
		})
	#MSlog <- sapply(conditions$condition, function(x){
	#	sum((log(wellInfo$protConc[wellInfo$Sample == x]) - fit_linMod$fitted[wellInfo$Sample == x])^2)/sum(wellInfo$Sample == x)
	#	})
	
	
	shrinkFrac <- sapply(wellInfo$Sample, function(x){MSlog[conditions$condition == x]})/(sapply(wellInfo$Sample, function(x){MSlog[conditions$condition == x]}) + rmse^2)
	
	
	#shrinkFrac <- 1 - 1/(abs(conclm2$resid)^2*(length(wellInfo[,1])/(length(wellInfo[,1]) - 10))/(rmse^2))
	#shrinkFrac <- sapply(shrinkFrac, function(x){max(0, x)})
	
	#shrinkFrac <- 1 - 1/(abs(conclm2$resid)^2*(length(wellInfo[,1])/(length(wellInfo[,1]) - 10))/(rmse^2))
	#shrinkFrac <- sapply(shrinkFrac, function(x){max(0, x)})

	updated_conc <- exp(fit_linMod$fitted)*shrinkFrac + (1-shrinkFrac)*wellInfo$protConc
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
#the log
plot(log(refittedconc[length(refittedconc[,1]),]) ~ factor(wellInfo$Sample)) 
plot(log(refittedconc[1,]) ~ factor(wellInfo$Sample))


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

plotting_plot <- ggplot(plotting_DF, aes(x = factor(Condition), y = Concentration, col = limitation)) + facet_wrap(~ Iteration, ncol = 1)
plotting_plot + geom_boxplot()

