library(data.table)
library(ggplot2)

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/chemostatInformation")

#### Load chemostat information ####
# original block of chemostats (relevantDR.txt) and follow-up experiments (followup_chemostats.txt)

chemostat_list <- read.delim('ListOfChemostats.txt')
primary_chemo <- data.table(read.delim('releventDR.csv', sep = ','))
followup_chemo <- data.table(read.delim('followup_chemostats.txt'))

if(!all(chemostat_list$Sheet %in% c("releventDR.csv", "followup_chemostats.txt"))){
  stop("add chemostat information to an existing sheet or update code")
  }

print(paste(length(chemostat_list$ChemostatCond[chemostat_list$include_this]), "chemostats chosen for analysis"))

redundant_conditions <- names(table(chemostat_list$ChemostatCond[chemostat_list$include_this]))[table(chemostat_list$ChemostatCond[chemostat_list$include_this]) != 1]
if(length(redundant_conditions) != 0){
  warning(paste(length(redundant_conditions), "conditions appear to be redundant, ensure that they are informed by seperate experimental data"))
  }


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

### Compare measures of culture density ###

# abundMetrics <- primary_chemo[,c('ChemostatID', 'Klett', 'colterVol.uL.mL', 'PCMVol.mL'),with = F]
# abundMetrics <- data.table(melt(abundMetrics, id.vars = "ChemostatID"))
# abundCV <- abundMetrics[,list(CV = sd(value, na.rm = T)/mean(value, na.rm = T)), by = c("ChemostatID", "variable")]
# ggplot(abundCV, aes(x = factor(variable), y = CV)) + geom_boxplot()

# Klett is the least noisy abundance measurement, while PCV probably gives the best measure of biomass

lm_klett_to_density <- lm(primary_chemo$PCMVol.mL ~ primary_chemo$Klett + 0)$coef
plot(primary_chemo$PCMVol.mL ~ primary_chemo$Klett, xlab = "Klett", ylab = "PCV") # coulter countVolume under-estimates culture volume relative to packed cell volume
abline(a = 0, b = lm_klett_to_density)

# use lm_klett_to_density to relate directly measured Klett to VolFrac

primary_chemo$VolFrac = primary_chemo$Klett * lm_klett_to_density
followup_chemo$VolFrac = followup_chemo$Klett * lm_klett_to_density


#### Clean-up chemostat information ####

primary_summary <- primary_chemo[,list(actualDR = mean(DR, na.rm = T), VolFrac_mean = mean(VolFrac, na.rm = T),
                                       VolFrac_SE = sd(VolFrac, na.rm = T)/sqrt(length(VolFrac[!is.na(VolFrac)])), medcellVol = median(coulterMed, na.rm = T)), by = "ChemostatID"]
followup_summary <- followup_chemo[,list(actualDR = mean(DR, na.rm = T), VolFrac_mean = mean(VolFrac, na.rm = T),
                                         VolFrac_SE = sd(VolFrac, na.rm = T)/sqrt(length(VolFrac[!is.na(VolFrac)])), medcellVol = median(coulterMed, na.rm = T)), by = "ChemostatID"]
run_info <- rbind(primary_summary, followup_summary)

if(!all(chemostat_list$ChemostatID[chemostat_list$include_this] %in% run_info$ChemostatID)){
  stop("You are missing chemostat running information (DR...)")
  }

#### Write chemostat information file containing cell density ####

ListOfChemostats_augmented <- merge(chemostat_list, run_info)
ListOfChemostats_augmented$Limitation <- factor(ListOfChemostats_augmented$Limitation, levels = c("P", "C", "N", "L", "U"))
ListOfChemostats_augmented <- ListOfChemostats_augmented[order(ListOfChemostats_augmented$Limitation, ListOfChemostats_augmented$DRgoal),]

write.table(ListOfChemostats_augmented, file = "ListOfChemostats_augmented.txt", sep = "\t", col.names = T,  row.names = F, quote = F)