#library("rNMR")
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

##### analysis ####

NMRsummary <- read.delim("effluentQuant_roiSummary.txt")
NMRsamples <- NMRsummary$File
NMRsamples <- unname(sapply(sapply(NMRsamples, function(x){strsplit(x, '2013_04_27_')[[1]][2]}), function(x){strsplit(x, '.ucsf')[[1]][1]}))
NMRsamples <- unname(sapply(NMRsamples, function(x){sub('ref_', 'p0.05H', x)}))
NMRsample_table <- data.frame(t(sapply(NMRsamples, function(x){strsplit(x, '_')[[1]]})))
rownames(NMRsample_table) <- unname(apply(NMRsample_table, 1, function(x){paste(x, collapse = "_")}))
colnames(NMRsample_table) <- c("SampleType", "Rep")
NMRsample_table$is_sample <- grep(NMRsample_table

NMRmatrix <- as.matrix(NMRsummary[,-1])
rownames(NMRmatrix) <- rownames(NMRsample_table)



heatmap.2(cor(NMRmatrix[,-2]))
heatmap.2(log2(NMRmatrix), Rowv = FALSE)

NMRstandards <- data.frame(name = NMRsample_table[grep('std', NMRsample_table[,1]),], grep('std', NMRsample_table[,1]))
colnames(NMRstandards) <- c("SampleType", "Rep", "Index")
NMRstandards$relative_conc = as.numeric(sapply(NMRstandards$SampleType, function(x){strsplit(x, 'xstd')[[1]]}))

standard_conc <- c(Glucose = 122, Ethanol = 244, Acetate = 10) #concentrations of 1x standard (mM)

NMRstandardConc <- NMRmatrix[NMRstandards$Index,colnames(NMRmatrix) %in% names(standard_conc)]
plot(NMRmatrix[c(1:length(NMRmatrix[,1]))[-NMRstandards$Index],colnames(NMRmatrix) %in% names(standard_conc)][,3])

### mapping standard abundance to concentration ####

for(i in 1:length(NMRstandardConc[1,])){
  print(plot(log(NMRstandards$relative_conc) ~ log(NMRstandardConc[,i]), xlab = "log_peak", ylab = "log_conc", main = colnames(NMRstandardConc)[i], pch = 16))
  print(plot(lm(log(NMRstandards$relative_conc) ~ log(NMRstandardConc[,i])), which = 1))
}

### for glucose a linear model is sufficient, while for acetate and ethanol a quadratic fit is needed

standard_fit_degree <- c(Glucose = 3, Ethanol = 2, Acetate = 3)
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

heatmap.2(NMRmatrix_HM/t(t(rep(1, length(NMRmatrix_HM[,1])))) %*% apply(NMRmatrix_HM, 2, sd), Rowv = FALSE)

melt(NMRmatrix_HM)
  


