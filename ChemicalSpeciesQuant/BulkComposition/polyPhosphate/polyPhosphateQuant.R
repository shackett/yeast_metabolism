setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/BulkComposition/polyPhosphate")

library(reshape2)
library(ggplot2)
library(nlme)

options(stringsAsFactors = FALSE)

PP_dataframe <- read.delim("phosphate_chemoQuant.txt", header = FALSE, sep = "\t")
PP_wellinfo <- read.delim("PP_wellinfo.txt", header = TRUE, sep = "\t")

for(i in 1:2){
  
    one_tech_rep <- PP_dataframe[c(1:9)+9*(i-1), -14]
    rownames(one_tech_rep) <- one_tech_rep[,1]; colnames(one_tech_rep) <- one_tech_rep[1,]
    one_tech_rep <- as.matrix(one_tech_rep[-1,-1])
    
    if(i == 1){
      PP_plate_set <- melt(one_tech_rep)
      }else{
        PP_plate_set <- cbind(PP_plate_set, melt(one_tech_rep)[,3])
        }
  }

colnames(PP_plate_set) <- c("Row", "Column", paste("R", 1:2, sep = ""))
PP_plate_matrix <- as.matrix(PP_plate_set[,-c(1:2)]) 


PP_standards <- NULL
for(i in 1:(ncol(PP_plate_set)-2)){
  PP_standards <- rbind(PP_standards, data.frame(plateNum = i,
    absorbance = PP_plate_set[PP_wellinfo$Type == "standard",][,i+2],
    abundance = PP_wellinfo[PP_wellinfo$Type == "standard",]$mg_PP,
    condition = PP_wellinfo[PP_wellinfo$Type == "standard",]$Limitation))
  }
std_plot <- ggplot(PP_standards, aes(y = absorbance, x = abundance, col = factor(condition))) + facet_wrap(~ plateNum, ncol = 2)
std_plot + geom_point()


for(i in 1:2){
  subset_plot <- PP_standards[PP_standards$plateNum == i & PP_standards$condition %in% c("Pi", "blank"),]
  print(plot(subset_plot$absorbance ~ subset_plot$abundance, pch = 16, main = i))
  abline(lm(subset_plot$absorbance ~ subset_plot$abundance)[1], lm(subset_plot$absorbance ~ subset_plot$abundance)[2], col = "RED")
  
  subset_plot <- subset_plot[subset_plot$abundance <= 0.0008 ,]
  abline(lm(subset_plot$absorbance ~ subset_plot$abundance)[1], lm(subset_plot$absorbance ~ subset_plot$abundance)[2], col = "GREEN")
  }

### remove high standard concentration, and lowest 3 concentrations cover experimental range and are highly linear
PP_plate_set <- PP_plate_set[is.na(PP_wellinfo$mg_PP) | PP_wellinfo$mg_PP <= 0.0008,]
PP_wellinfo <- PP_wellinfo[is.na(PP_wellinfo$mg_PP) | PP_wellinfo$mg_PP <= 0.0008,]


#### look at the subset of plates where the standards were highly linear

samples <- PP_wellinfo[PP_wellinfo$Type == "sample",]
samples$cond <- mapply(function(x,y){
  paste(c(x, y), collapse = "")
  }, x = samples$Limitation, y = samples$DR)
samples$cond <- sub('0.3', '0.30', samples$cond)

weight_frac <- NULL
for(i in 1:2){
  
  #data.frame(PP_wellinfo$mg_PP[PP_wellinfo$Type == "standard"], PP_plate_set[PP_wellinfo$Type == "standard",i+2])
  
  std_lm <- lm(PP_wellinfo$mg_PP[PP_wellinfo$Type == "standard"] ~ PP_plate_set[PP_wellinfo$Type == "standard",i+2])
  
  mgP <- std_lm$coef[1] + std_lm$coef[2]*PP_plate_set[PP_wellinfo$Type == "sample",i+2]
  weight_frac <- cbind(weight_frac, mgP/PP_wellinfo$mg_dry_well[PP_wellinfo$Type == "sample"])
  #plot(mgC/TC_wellinfo$mg_dry_well[TC_wellinfo$Type == "sample"] ~ factor(samples$cond))
  }

weight_frac <- data.frame(samples$cond, weight_frac)
colnames(weight_frac) <- c("condition", paste("R", 1:2, sep = ""))

all_plates <- melt(data.frame(weight_frac), id.vars = "condition")
all_plates$cond_plate <- factor(mapply(function(x,y){paste(c(x,y), collapse = "_")}, x = all_plates$condition, y = all_plates$variable))
all_plates$condition <- factor(all_plates$condition)
all_plates$variable <- factor(all_plates$variable)
all_plates$value <- as.numeric(all_plates$value)
all_plates$limitation <- toupper(substr(all_plates$condition, 1, 1))

ggplot(all_plates, aes(x = condition, y = value, col = limitation)) + geom_boxplot() + facet_grid (variable ~ .) + ylim(0,0.15)
ggsave("polyphosphateAbund.pdf", height = 8, width = 10)

all_plates <- data.table(all_plates)
all_plates[, pointEstimate := mean(value), by = condition]
all_plates[, residual := value,]


all_plates[condition == "p0.05", mean(value), by = condition]

PPsummary <- all_plates[, unique(pointEstimate), by = condition]
colnames(PPsummary) <- c("condition", "fraction")

PPsummary$indVariance <- all_plates[, var(residual - pointEstimate), by = condition]$V1
PPsummary$assayVariance <- all_plates[, var(residual - pointEstimate),]

write.table(PPsummary, file = "../PPfrac_data.txt", sep = "\t", col.names = TRUE, row.names = F, quote = F)

#### mL equivalents of media corresponding to phosphate content of cells ######
mgPPdry <- PPsummary$fraction * sapply(PPsummary$condition, function(x){PP_wellinfo$mg_dry[paste(PP_wellinfo$Limitation, PP_wellinfo$DR, sep = "") == x][1]})
### ~ 0.6 g/L non-P-lim
### ~ 6 mg/L P-lim
#mgPPdry/600 - 2 - 5mL of replete media worth of Pi in cell pellet.  This is much greater than what would be left after centrifugation and washing.
