library(data.table)

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/chemostatInformation")

#### Load chemostat information ####
# original block of chemostats (relevantDR.txt) and follow-up experiments (followup_chemostats.txt)

chemostat_list <- read.delim('ListOfChemostats.txt')
primary_chemo <- data.table(read.delim('releventDR.csv', sep = ','))
followup_chemo <- data.table(read.delim('followup_chemostats.txt'))

# For each chemostat we want a functional summary of density and treatment:
# ChemostatCond : the chemostat treatment/nutrient condition
# ChemostatID : a unique identifier for each chemostat
# limitation : what type of chemostat media used
# DRgoal : the dilution rate that was sought
# actualDR : the dilution rate achieved
# VolFrac_mean : uL cell volume per mL culture
# VolFrac_SE : standard error of VolFrac_mean
# medcellVol : median cell volume (fL)
# NMRdate : point to which date is associated with the chemostat's NMR data

lm_coult_to_density <- lm(primary_chemo$PCMVol.mL ~ primary_chemo$colterVol.uL.mL + 0)$coef
plot(primary_chemo$PCMVol.mL ~ primary_chemo$colterVol.uL.mL, xlab = "coulter volume", ylab = "PCV") # coulter countVolume under-estimates culture volume relative to packed cell volume
abline(a = 0, b = lm_coult_to_density)

lm_klett_to_density <- lm(primary_chemo$PCMVol.mL ~ primary_chemo$Klett + 0)$coef
plot(primary_chemo$PCMVol.mL ~ primary_chemo$Klett, xlab = "Klett", ylab = "PCV") # coulter countVolume under-estimates culture volume relative to packed cell volume
abline(a = 0, b = lm_klett_to_density)

# use lm_coult_to_density to relate directly measured coulter volumes to VolFrac
primary_chemo$VolFrac = primary_chemo$colterVol.uL.mL * lm_coult_to_density
followup_chemo$VolFrac = followup_chemo$colterVol.uL.mL * lm_coult_to_density

### Clean-up the primary chemostat information ###

primary_summary <- primary_chemo[,list(actualDR = mean(DR, na.rm = T), VolFrac_mean = mean(VolFrac, na.rm = T),
                                       VolFrac_SD = sd(VolFrac, na.rm = T), medcellVol = median(coulterMed, na.rm = T),
                                       NMRdate = "2013_04_27"), by = "Chemostat"]
primary_summary[,ChemostatCond := toupper(Chemostat),]
primary_summary$ChemostatCond <- sub('CTRL[0-9]_', 'P', primary_summary$ChemostatCond)
primary_summary[,ChemostatID := paste(ChemostatCond, "2011", sep = "_"),]
primary_summary$ChemostatID[primary_summary$Chemostat == "ctrl1_0.05"] <- "P0.05H1_2011"
primary_summary$ChemostatID[primary_summary$Chemostat == "ctrl2_0.05"] <- "P0.05H2_2011"


### Clean-up followup chemostat information ###

followup_summary <- followup_chemo[,list(actualDR = mean(DR, na.rm = T), VolFrac_mean = mean(VolFrac, na.rm = T),
                                         VolFrac_SD = sd(VolFrac, na.rm = T), medcellVol = median(coulterMed, na.rm = T),
                                         NMRdate = NMR_date[NMR_date != ""]), by = "ChemostatCond"]
followup_chemo$

followup_chemo[,list(actualDR = mean(DR, na.rm = T), 



