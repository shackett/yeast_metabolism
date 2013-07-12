#library("rNMR")
setwd("/Users/Sean/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Media")

options(stringsAsFactors = FALSE)
library(reshape2)
library(gplots)
library(ggplot2)

### Convert bruker files to UCSF format - select the highest level folder 
cf()

### open UCSF files
fs()

### open region of interest (roi) tab
roi()

write.table(roiTable, "RoiTable.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = F)

##### Analysis #######

### Describe the NMR samples ####

NMRsummary <- read.delim("effluentQuant_vz_roiSummary.txt")
NMRsamples <- NMRsummary$File
NMRsamples <- unname(sapply(sapply(NMRsamples, function(x){strsplit(x, '2013_04_27_')[[1]][2]}), function(x){strsplit(x, '.ucsf')[[1]][1]}))
NMRsamples <- unname(sapply(NMRsamples, function(x){sub('ref_', 'p0.05H', x)}))
NMRsample_table <- data.frame(t(sapply(NMRsamples, function(x){strsplit(x, '_')[[1]]})))
rownames(NMRsample_table) <- unname(apply(NMRsample_table, 1, function(x){paste(x, collapse = "_")}))
colnames(NMRsample_table) <- c("SampleType", "Rep")

NMRsample_table$is_sample <- TRUE; NMRsample_table$is_sample[grep('std', NMRsample_table$SampleType)] <- FALSE
NMRsample_table$DR <- NMRsample_table$limitation <- NA

NMRsample_table$limitation[NMRsample_table$is_sample] <- sapply(NMRsample_table$SampleType[NMRsample_table$is_sample], function(x){strsplit(x, '[0-9.]+')[[1]][1]})
NMRsample_table$DR[NMRsample_table$is_sample] <- sapply(NMRsample_table$SampleType[NMRsample_table$is_sample], function(x){regmatches(x, regexpr('[0-9.]{4}', x))})



### Create a matrix of NMR data ####

NMRmatrix <- as.matrix(NMRsummary[,-1])
rownames(NMRmatrix) <- rownames(NMRsample_table)
colnames(NMRmatrix)[grep('Glucose',colnames(NMRmatrix))] = 'Glucose'
colnames(NMRmatrix)[grep('Ethanol',colnames(NMRmatrix))] = 'Ethanol'
colnames(NMRmatrix)[grep('Glycerol',colnames(NMRmatrix))] = 'Glycerol'


heatmap.2(cor(NMRmatrix[,-1]))
heatmap.2(log2(NMRmatrix[,-1]), Rowv = FALSE)

### Regression of standard peak areas onto known concentration ####

NMRstandards <- data.frame(name = NMRsample_table[grep('std', NMRsample_table[,1]),], grep('std', NMRsample_table[,1]))[,c(1,2,6)]
colnames(NMRstandards) <- c("SampleType", "Rep", "Index")
NMRstandards$relative_conc = as.numeric(sapply(NMRstandards$SampleType, function(x){strsplit(x, 'xstd')[[1]]}))

standard_conc <- c(Glucose = 122, Ethanol = 244, Acetate = 10, Glycerol = 10) #concentrations of 1x standard (mM)

NMRstandardConc <- NMRmatrix[NMRstandards$Index,colnames(NMRmatrix) %in% names(standard_conc)]
plot(NMRmatrix[c(1:length(NMRmatrix[,1]))[-NMRstandards$Index],colnames(NMRmatrix) %in% names(standard_conc)][,3])

### mapping standard abundance to concentration ####

for(i in 1:length(NMRstandardConc[1,])){
  print(plot(log(NMRstandards$relative_conc) ~ log(NMRstandardConc[,i]), xlab = "log_peak", ylab = "log_conc", main = colnames(NMRstandardConc)[i], pch = 16))
  print(plot(lm(log(NMRstandards$relative_conc) ~ log(NMRstandardConc[,i])), which = 1))
}

### for ethanol a linear model is sufficient, while for acetate and glucose a quadratic fit is needed

standard_fit_degree <- c(Glucose = 3, Ethanol = 2, Acetate = 3, Glycerol = 3)
standard_fit_lm <- list()

for(i in colnames(NMRstandardConc)){
  
  log_pred_peak <- log2(NMRstandardConc[,colnames(NMRstandardConc) == i])
  log_ref_conc <- log2(standard_conc[names(standard_conc) == i] * NMRstandards$relative_conc)
  
  if(standard_fit_degree[names(standard_fit_degree) == i] == 2){
    standard_fit_lm[[i]] <- lm(log_ref_conc ~ log_pred_peak)
  }else if(standard_fit_degree[names(standard_fit_degree) == i] == 3){
    standard_fit_lm[[i]] <- lm(log_ref_conc ~ log_pred_peak + I(log_pred_peak^2))
    standard_fit_lm[[i]] <- lm(log_ref_conc ~ poly(log_pred_peak, 2, raw = TRUE))
  }else{print("check yo self")}
  
  #predict(standard_fit_lm[[i]], data.frame(log_pred_peak = log_pred_peak), interval = "confidence")
}

### Determining sample concentrations from either comparison to standard regression or DSS ####

NMRmatrix_conc <- NMRmatrix

for(i in 1:length(NMRmatrix_conc[1,])){
  if(colnames(NMRmatrix_conc)[i] %in% names(standard_conc)){
    ### scale peak height to standards ###
    NMRmatrix_conc[,i] <- 2^predict(standard_fit_lm[[colnames(NMRmatrix_conc)[i]]], data.frame(log_pred_peak = log2(NMRmatrix_conc[,i])))
    
    }else{
      ### scale peak height to DSS molarity (5mM*(1/9)) (although proton number of a peak affects this)
      NMRmatrix_conc[,i] <- sapply(NMRmatrix_conc[,i], function(x){max(x, 0)}) * 5/9
      }
  }

NMRmatrix_HM <- NMRmatrix_conc[c(1:length(NMRmatrix[,1]))[-NMRstandards$Index],grep('DSS', colnames(NMRmatrix_conc), invert = TRUE) ]
NMRdesign <- NMRsample_table[NMRsample_table$is_sample,]

heatmap.2(NMRmatrix_HM/t(t(rep(1, length(NMRmatrix_HM[,1])))) %*% apply(NMRmatrix_HM, 2, sd), Rowv = FALSE)

### Point estimation (and variance) of concentration in each condition ####

NMR_point_estimate <- NULL

### contrasts are relative to zero rather than intercept

treatment_mat <- matrix(0, ncol = length(unique(NMRdesign$SampleType)), nrow = length(NMRdesign[,1]))
colnames(treatment_mat) <- unique(NMRdesign$SampleType)
for(t in 1:length(treatment_mat[1,])){
  treatment_mat[,t][NMRdesign$SampleType == colnames(treatment_mat)[t]] <- 1
  }

for(j in 1:length(treatment_mat[1,])){
  
  outputDF <- data.frame(condition = colnames(treatment_mat)[j], peak = colnames(NMRmatrix_HM), estimate = apply(NMRmatrix_HM[treatment_mat[,j] == 1,], 2, mean), se = apply(NMRmatrix_HM[treatment_mat[,j] == 1,], 2, sd)/sqrt(colSums(!is.na(NMRmatrix_HM[treatment_mat[,j] == 1,]))))
  outputDF$lb <- outputDF$estimate - 1.96*outputDF$se
  outputDF$ub <- outputDF$estimate + 1.96*outputDF$se
  rownames(outputDF) <- NULL
  
  NMR_point_estimate <- rbind(NMR_point_estimate, outputDF)
  
  }

cond_finder <- t(sapply(NMR_point_estimate$condition, function(x){NMRdesign[NMRdesign$SampleType == x,][1,]}))[,4:5]
NMR_point_estimate <- cbind(NMR_point_estimate, data.frame(limitation = unlist(cond_finder[,1]), DR = unlist(cond_finder[,2])))

barplot_theme <- theme(text = element_text(size = 20, face = "bold"), title = element_text(size = 20, face = "bold"), panel.background = element_blank(), legend.position = "none", 
            panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.ticks.x = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle = 90)) 

#NMR_point_estimate  <- NMR_point_estimate[NMR_point_estimate$peak %in% c("Ethanol", "Acetate", "Glucose", "Glycerol"),]
  
NMRbarplot <- ggplot(NMR_point_estimate, aes(x = factor(condition), y = estimate, fill = factor(limitation))) + facet_grid(peak ~ ., scale = "free_y") + barplot_theme
NMRbarplot + geom_bar(stat = "identity") + ggtitle("Concentration of metabolites (and unknowns) in chemostat effluent") + scale_fill_brewer(palette = "Set2") + 
  geom_errorbar(aes(ymin = lb, ymax = ub)) + scale_x_discrete("Chemostat condition") + scale_y_continuous("Concentration (mM)")
ggsave("mediaComposition.pdf", width = 14, height = 22)
write.table(NMR_point_estimate, "mediaComposition_NMR.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = F)
