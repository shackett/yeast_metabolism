setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/FTIR")

library(gplots)
library(ggplot2)
library(reshape)

FTIR_abund_list <- list()
FTIR_abund_tallDF <- NULL

directories <- NULL
dates <- list.files(paste(getwd(),"/FTIR_samples/",sep = ""))

FTIR_file_info <- read.delim("FTIR_samples.txt")
FTIR_file_info$Valid[is.na(FTIR_file_info$Valid) == TRUE] <- TRUE

compound_type <- data.frame(class = c(rep("protein", 3), rep("carbohydrate", 5), rep("D/RNA", 2), rep("lipid-related", 3), rep("misc", 3)) , compound = c("BSA", "chick_albumin", "chymotrypsin", "glucan", "chitin", "glycogen", "mannan", "trehalose", "DNA", "RNA", "glycerol", "triton-X", "palmitate", "KBr", "KPO4-EDTA", "poly-phosphate"))

unique_compounds <- unique(FTIR_file_info$Compound)[unique(FTIR_file_info$Compound) != "sample"]	

wavelengths <- NULL
for(cmpd in unique(FTIR_file_info$Compound)){
	
	cmpd_info <- FTIR_file_info[FTIR_file_info$Compound == cmpd,]
	cmpd_info <- cmpd_info[cmpd_info$Valid,]
	
	absorbance <- NULL
	
	if(length(cmpd_info[,1]) == 0){
		print(paste(cmpd, "has no valid samples"))
		next
		}
		
	for(runN in 1:length(cmpd_info[,1])){
		run_data_file <- read.table(file.path("FTIR_samples", cmpd_info$Date_Folder[runN], paste(cmpd_info$Sample_Name[runN], ".csv", sep = "")), sep = ",", stringsAsFactors = FALSE)
		
		if(length(wavelengths) == 0){
			wavelengths <- run_data_file[,1]
			}
		if(sum(wavelengths != run_data_file[,1]) != 0){
			print(paste("wavelengths don't match for", runN))
			}
		
		absorbance <- cbind(absorbance, run_data_file[,2])
		
		}
	
  
  
	rownames(absorbance) <- wavelengths[wavelengths %in% run_data_file[,1]]; colnames(absorbance) <- paste("R", c(1:length(cmpd_info[,1])), sep = "")
	
	#remove effect of negative baseline
	absorbance <- absorbance - apply(absorbance, 2, min)
	#orthonormalize
	absorbance <- absorbance/apply(absorbance, 2, sum)
	
	FTIR_abund_list[[cmpd]] <- absorbance
	
	FTIR_abund_tallDF <- rbind(FTIR_abund_tallDF, cbind(cmpd, melt(absorbance)))
		
	}	


colnames(FTIR_abund_tallDF) <- c("Compound", "Wavenumber", "Replicate", "Absorbance")
FTIR_abund_tallDF$Compound <- factor(FTIR_abund_tallDF$Compound)

#plot by compound
plot_abs <- ggplot(data = FTIR_abund_tallDF, aes(x = Wavenumber, y = Absorbance))
plot_abs <- plot_abs + facet_wrap( ~ Compound, nrow = 4)
plot_abs + geom_area()

unique(FTIR_abund_tallDF$Compound)[!unique(FTIR_abund_tallDF$Compound) %in% compound_type$compound]

#plot by type and color by compound
FTIR_abund_tallDF$Class <- sapply(FTIR_abund_tallDF$Compound, function(x){compound_type$class[compound_type$compound %in% x]})
FTIR_abund_tallDF$Samples <- as.factor(sapply(c(1:length(FTIR_abund_tallDF[,1])), function(a_row){paste(c(as.character(FTIR_abund_tallDF$Compound[a_row]), as.character(FTIR_abund_tallDF$Replicate[a_row])), collapse = "-")}))

plot_abs_class <- ggplot(data = FTIR_abund_tallDF, aes(x = Wavenumber, y = Absorbance, group = Samples, col = Compound))
plot_abs_class <- plot_abs_class + facet_wrap( ~ Class, ncol = 1)
plot_abs_class + geom_line(size = 1)

	
	
