##### Combine experimental data on yeast composition, nutrient intake and excretion and expected composition into a set of condition-specific boundary fluxes #####

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant")
options(stringsAsFactors = FALSE)

library(ggplot2)

### Default composition function ####

compositionFile <- read.csv2("../Yeast_comp.csv", sep = ",", stringsAsFactors = FALSE)

modelMetComp <- read.table("stoiMetsComp.tsv", header = TRUE)

compComp <- data.frame(t(sapply(compositionFile$AltName, function(x){
  unlist(modelMetComp[modelMetComp[,2] == x,][1,])
  }))) #elemental composition each macromolecule in biomass function

compComp[compComp$name == "(1->3)-beta-D-glucan", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5) #one oxygen shared bc of condensation
compComp[compComp$name == "glycogen", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5)
compComp[compComp$name == "mannan", colnames(compComp) %in% c("C", "H", "O")] <- c(6, 10, 5)
compComp[is.na(compComp)] <- 0
compComp[,-c(1:2)] <- apply(compComp[,-c(1:2)], c(1,2), as.numeric)
compComp <- compComp[,c(TRUE, TRUE, colSums(compComp[,-c(1:2)]) != 0)]

atomicMasses <- data.frame(element = c("C", "H", "N", "O", "P", "S"), mass = c(12.0107, 1.00794, 14.00674, 15.9994, 30.973761, 32.066))

compositionFile$MW <- t(t(compComp[,-c(1,2)])) %*% t(t(c(atomicMasses[,2])))
compositionFile$weight_per_t <- as.numeric(compositionFile$StoiCoef) * compositionFile$MW
#compositionFile$weight_per_t[compositionFile$Class == "Energy Balance"] <- NA

class_composition <- sapply(unique(compositionFile$Class), function(x){sum(compositionFile$weight_per_t[compositionFile$Class == x])}) * -1
class_composition <- data.frame(Category = names(class_composition[names(class_composition) != "Energy Balance"]), Abundance = unname(class_composition[names(class_composition) != "Energy Balance"]))
class_composition$CellComposition <- "Default"


pie_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "white"), legend.position = "right", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(3, "line"), axis.text = element_blank(), axis.title = element_blank()) 


class_comp_plot <- ggplot(class_composition, aes(y = Abundance, factor(CellComposition), fill = factor(Category))) + pie_theme
class_comp_plot + geom_bar(stat = "identity") + coord_polar(theta = "y") + scale_fill_discrete(name = "Class")


########### Chemostat specific information ###########

chemostatInfo <- read.table("BulkComposition/chemostatDRmisc.tsv", sep = "\t", header = TRUE)







##### Total Protein ######

load('BulkComposition/protSpecQuantOut.Rdata')

added_conds <- conditions[conditions$condition == "p0.30",]; added_conds$condition <- "p0.30_ORBI"
conditions <- rbind(conditions, added_conds)

cond_info <- sapply(colnames(good_light), function(cond){
  if(strsplit(cond, "")[[1]][1] == "L"){
		tmp <- cond
		}else{
			tmp <- paste(c(tolower(strsplit(cond, "")[[1]][1]), strsplit(cond, "")[[1]][-1]), collapse = "")
			}
		c(1:length(conditions[,1]))[conditions$condition == tmp]
		})


##### Total RNA ######





##### Total DNA #####

## genome size + genome size * fraction of cells not in G1 # from Brauer
0.936–1.971D #brauer 2008 relationship between unbudded fraction and growth rate



