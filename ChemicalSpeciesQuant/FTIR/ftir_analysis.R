setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/FTIR")

library(gplots)
library(ggplot2)
library(reshape)

FTIR_abund_list <- list()

directories <- NULL
dates <- list.files(paste(getwd(),"/FTIR_samples/",sep = ""))

for(i in 1:length(dates)){
	all_files <- list.files(paste(getwd(),"/FTIR_samples/",dates[i],sep = ""))
	csv_files <- all_files[grep(".CSV", all_files)]
	
	
	FTIR_abund_list[[dates[i]]]$compounds <- t(sapply(csv_files, function(ID){
		prefix <- unlist(strsplit(ID, "_"))[1]
		if(length(grep('[0-9]', prefix)) >= 1){
			c(unlist(strsplit(prefix, '[0-9]'))[1], unlist(strsplit(prefix, '[-A-Za-z]+'))[2])
			}else{
				c(prefix, 1)
				}}))
		
	colnames(FTIR_abund_list[[dates[i]]]$compounds) <- c("Compound", "Measurement")
	
	wavelengths <- NULL
	absorbance <- NULL
	
	for(run in rownames(FTIR_abund_list[[dates[i]]]$compounds)){
		run_data_file <- read.table(file.path("FTIR_samples", dates[i], run), sep = ",", stringsAsFactors = FALSE)
		if(length(wavelengths) == 0){
			wavelengths <- run_data_file[,1]
			}
		
		if(sum(wavelengths != run_data_file[,1]) != 0){
			print(paste("wavelengths don't match for", run))
			}
		absorbance <- cbind(absorbance, run_data_file[,2])
		}
	rownames(absorbance) <- wavelengths; colnames(absorbance) <- rownames(FTIR_abund_list[[dates[i]]]$compounds)
	
	valid_samples <- sapply(unique(FTIR_abund_list[[dates[i]]]$compounds[,colnames(FTIR_abund_list[[dates[i]]]$compounds) == "Compound"]), function(unique_cmpd){
		c(1:length(FTIR_abund_list[[dates[i]]]$compounds[,1]))[FTIR_abund_list[[dates[i]]]$compounds[,colnames(FTIR_abund_list[[dates[i]]]$compounds) == "Compound"] == unique_cmpd][FTIR_abund_list[[dates[i]]]$compounds[,colnames(FTIR_abund_list[[dates[i]]]$compounds) == "Measurement"][FTIR_abund_list[[dates[i]]]$compounds[,colnames(FTIR_abund_list[[dates[i]]]$compounds) == "Compound"] == unique_cmpd] == max(FTIR_abund_list[[dates[i]]]$compounds[,colnames(FTIR_abund_list[[dates[i]]]$compounds) == "Measurement"][FTIR_abund_list[[dates[i]]]$compounds[,colnames(FTIR_abund_list[[dates[i]]]$compounds) == "Compound"] == unique_cmpd])]
		})
	
	compound_type <- data.frame(class = c(rep("protein", 3), rep("carbohydrate", 5), rep("D/RNA", 2), rep("lipid-related", 3), rep("misc", 3)) , compound = c("BSA", "chick-albumin", "chymotrypsin", "Glucan", "Chitin", "Glycogen", "Mannan", "trehalose", "DNA", "RNA", "glycerol(KBr)", "Triton-X(KBr)", "Palmitate", "KBr", "pottasium-phosphate-EDTA", "poly-phosphate"))
	
	
	
	valid_samples <- valid_samples[!(names(valid_samples) %in% c("KBr"))]
	
	valid_sample_Abs <- absorbance[,valid_samples]
	colnames(valid_sample_Abs) <- names(valid_samples)
	
	valid_sampleNorm <- valid_sample_Abs*1/(t(t(rep(1, times = length(wavelengths)))) %*% t(apply(valid_sample_Abs, 2, sum)))
	
	plotting_frame <- melt(valid_sampleNorm)
	colnames(plotting_frame) <- c("Wavelength", "Compound", "Absorbance")
	
	plot_abs <- ggplot(data = plotting_frame, aes(x = Wavelength, y = Absorbance))
	plot_abs <- plot_abs + facet_wrap( ~ Compound, nrow = 4)
	plot_abs + geom_area()
	
	}
	
	
