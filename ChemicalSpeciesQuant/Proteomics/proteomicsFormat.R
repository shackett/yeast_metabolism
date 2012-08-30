####### An EM to go from peptides to proteins #######

setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')

load('20120807ProtPepMatrices.Rdata')
source("pep_library.R")
plotting_fxn()
matrix_fxn()

# Filter to only look at peptides with less than quality_frac fraction of missing values

quality_frac <- 0.8
ICthreshold <- 2^15

lightIC[lightIC < ICthreshold] <- NA; heavyIC[heavyIC < ICthreshold] <- NA

good_samples <- rowSums(is.finite(PepMatrix) & !is.na(lightIC) & !is.na(heavyIC)) >= (length(PepMatrix[1,])*quality_frac)

#the relative abundance of a peptide across conditions w.r.t a common reference
abundMat <- PepMatrix[good_samples,]

#possible mappings between a protein and all matching peptides
mappingMat <- ProtPepMatrix[good_samples,]
mappingMat <- mappingMat[,colSums(mappingMat) != 0]
good_light <- lightIC[good_samples,]
good_heavy <- heavyIC[good_samples,]

#map measured peaks to unique peptide sequences

pepNames <- rownames(mappingMat)
unique_NameCorr <- sapply(pepNames, function(name){
	unlist(strsplit(name, "\\."))[1]
	})
unique_pepNames <- unique(unique_NameCorr)
unique_NameCorrCol <- sapply(unique_NameCorr, function(name){
	c(1:length(unique_pepNames))[name == unique_pepNames]
	})

pepToUniq <- matrix(0, nrow = length(pepNames), ncol = length(unique_pepNames))
rownames(pepToUniq) <- pepNames; colnames(pepToUniq) <- unique_pepNames
for(i in 1:length(pepNames)){
	pepToUniq[i,unique_NameCorrCol[i]] <- 1
	}
pepToUniq <- Matrix(pepToUniq)

#determine the expectation of the standard deviation as a heteroschedastic fxn of IC using p0.05 light v. p0.05 heavy

avgSignalSTD <- apply(cbind(heavyIC[,colnames(heavyIC) == "P0.05"], lightIC[,colnames(lightIC) == "P0.05"]), 1, mean)
logLight <- log2(lightIC[is.finite(avgSignalSTD),colnames(lightIC) == "P0.05"])
logHeavy <- log2(heavyIC[is.finite(avgSignalSTD),colnames(heavyIC) == "P0.05"])
avgSignalSTD <- avgSignalSTD[is.finite(avgSignalSTD)]
logLight <- logLight + optimize(normFactor, c(-1, 1), logLight = logLight, logHeavy = logHeavy)$minimum

#variance of the replicate differences accounting for small sample (scaling factor of 2)
STDvar <- (logLight - logHeavy)^2*2
STDvar_fit <- lm(log2(STDvar) ~ log2(avgSignalSTD))$coef
STDvar_fit_initial <- STDvar_fit

gplot.hexbin(hexbin(logLight, logHeavy, xbins = 200), colramp = rainbow)

sd_reg_plots <- function(STDvar, avgSignalSTD, STDvar_fit){

	print(plot(log2(STDvar) ~ log2(avgSignalSTD), pch = 16, cex = 0.3))
	abline(STDvar_fit, col = "RED")
	abline(STDvar_fit_initial, col = "GREEN")
	print(gplot.hexbin(hexbin(log2(avgSignalSTD), log2(STDvar), xbins = 80), colramp = rainbow))
	
	}

####

function(generate_plots = FALSE){

n_tp <- length(good_light[,1])
n_req_members <- 3

allAttribPep <- pepToUniq %*% (mixing_fract == 1)[,1:n_prot]
valid_prot <- c(1:n_prot)[colSums(allAttribPep) >= n_req_members]
valid_pep <- c(1:n_tp)[rowSums(allAttribPep[,valid_prot]) == 1]

nsharedpeps <- colSums(allAttribPep[,valid_prot]) %*% t(mixing_fract[,valid_prot]) %*% t(pepToUniq[valid_pep,])
residuals <- as.matrix(prot_abund[,valid_prot]  %*% t(mixing_fract[,valid_prot]) %*% t(pepToUniq[valid_pep,]) - t(abundMat[valid_pep,]))
sqresid_corrected <- as.vector(t(residuals^2 * (t(t(rep(1, times = n_c))) %*% (nsharedpeps/(nsharedpeps-1)))))

sampleIC <- mapply(x = c(good_light[valid_pep,]), y = c(good_heavy[valid_pep,]), function(x,y){
	if(!is.na(x) & !is.na(y)){mean(x,y)}else{NA}
	})

zero_thresh <- 0.00001
STDvar_fit <- lm(log2(sqresid_corrected[!is.na(sampleIC) & sqresid_corrected > zero_thresh]) ~ log2(sampleIC[!is.na(sampleIC) & sqresid_corrected > zero_thresh]))

if(generate_plots == TRUE){
	sd_reg_plots(sqresid_corrected[!is.na(sampleIC) & sqresid_corrected > zero_thresh], sampleIC[!is.na(sampleIC) & sqresid_corrected > zero_thresh], STDvar_fit)
	}
	
z <- cbind(log2(sqresid_corrected[!is.na(sampleIC) & sqresid_corrected > zero_thresh]), log2(sampleIC[!is.na(sampleIC) & sqresid_corrected > zero_thresh]))
plot(z[,2] ~ z[,1])
abline(a = mean(z[,1], b	



#calculate the expected sampling variance of the heavy-low diff for experimental measurement

lightHeavyCellmean <- sapply(c(1:length(heavyIC[,1]))[good_samples], function(row){
	mapply(FUN = function(a,b){
		if(!is.na(a) & !is.na(b)){mean(a, b)}else{NA}
		}, lightIC[row,], heavyIC[row,])
	})

fittedVar <- var_calc(lightHeavyCellmean, STDvar_fit)
fittedPrec <- fittedVar^-1

#for each unique peptide, combine the multiple ionization states to produce a single point estimate, using integrated likelihood

#set the precision of missing values to 0; equivalent to no impact, infinite variance

fittedPrec[t(is.na(abundMat))] <- 0
fittedPrec[is.na(fittedPrec)] <- 0
abundMat[t(fittedPrec) == 0] <- 0

n_p = length(unique_pepNames) #4042
n_pp = length(unique_mappingMat[,1]) + length(unique_mappingMat[1,]) #4964
n_prot <- n_pp - n_p
n_c = length(uniquePepMean[,1]) #15


uniquePepMean <- matrix((((t(abundMat) * fittedPrec) %*% pepToUniq)/(fittedPrec %*% pepToUniq)), ncol = n_p, nrow = n_c)
uniquePepPrecision <- Matrix(fittedPrec %*% pepToUniq)
uniquePepMean[is.nan(uniquePepMean)] <- NA

#number of non-missing values for peptides
Nmissing_val <- table(rowSums(!is.na(t(uniquePepMean))))
barplot(Nmissing_val, lwd = 5)
uniquePepMean[is.na(uniquePepMean)] <- 0


#change mapping from peptides to unique peptides (averaging over ionization states)
unique_mappingMat <- t(pepToUniq) %*% mappingMat
unique_mappingMat[!(unique_mappingMat %in% c(0,1))] <- 1


######## EM ######

save(unique_pepNames, uniquePepMean, unique_mappingMat, uniquePepPrecision, file = "peptide.files.Rdata")
rm(list = ls()); gc()
load("peptide.files.Rdata")
source("pep_library.R")
matrix_fxn()
prerun_fixed_mat <- TRUE

### Initalization ###

n_p = length(unique_pepNames) #4042
n_pp = length(unique_mappingMat[,1]) + length(unique_mappingMat[1,]) #4964
n_prot <- n_pp - n_p
n_c = length(uniquePepMean[,1]) #15

prior_p_div <- exp(-1*qchisq(0.999, n_c))


tmp <- diag(rep(prior_p_div, times = n_p)); colnames(tmp) <- paste(unique_pepNames, "divergent", sep = "_")
mixing_fract <- unique_mappingMat/((1-prior_p_div)^-1*rowSums(unique_mappingMat))
mixing_fract <- cbind(mixing_fract, tmp)
mixing_fract <- Matrix(mixing_fract)
prior_mat_likadj <- unique_mappingMat/((1-prior_p_div)^-1)
prior_mat_likadj <- cbind(prior_mat_likadj, tmp)

tmp[!(tmp %in% c(0,1))] <- 1 
prior_mat <- cbind(unique_mappingMat, tmp) 
prior_mat_logical <- prior_mat; prior_mat_logical <- prior_mat_logical == 1
prior_mat_sparse <- Matrix(prior_mat_logical)
prior_mat_likadj <- prior_mat_likadj*prior_mat_sparse
prior_mat_likadj[prior_mat_sparse] <- log(prior_mat_likadj[prior_mat_sparse])
rm(tmp, prior_mat, prior_mat_logical); gc(); gc()

#preprocess sparse data and precision matrices

if(prerun_fixed_mat == TRUE){
#lower processor usage, higher memory
sampleEst_list <- list()
samplePrec_list <- list()
for(c in 1:n_c){
	sampleEst_list[[c]] <- uniquePepMean[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	samplePrec_list[[c]] <- 0.5*(uniquePepPrecision[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse)
	print(c)
	}
	};gc()
	
#fractSVD <- svd(t(mixing_fract))
#pseudoinv <- fractSVD$v%*%solve(diag(fractSVD$d))%*%t(fractSVD$u)
#prot_abund <- uniquePepMean%*%pseudoinv
#rownames(prot_abund) <- rownames(uniquePepMean); colnames(prot_abund) <- colnames(mixing_fract)
#prot_prec <- uniquePepPrecision %*% mixing_fract
	

### Iteration ###

whole_data_logL <- NULL
previous_it <- -Inf
continue_it <- TRUE

t <- 1
while(continue_it){
	
	#update protein abundances: using integrated likelihood
	#update protein means
	prot_abund = ((uniquePepMean*uniquePepPrecision) %*% mixing_fract)/(uniquePepPrecision %*% mixing_fract)
	prot_abund <- as.matrix(Matrix(prot_abund, sparse = FALSE))
	prot_abund[is.nan(prot_abund)] <- 0
	
	#update protein precision
	prot_prec <- uniquePepPrecision %*% mixing_fract
	
	#update mixing fract
	
	#evaluate the log-likelihood of sample abundances about the inferred protein abundances with a precision-specific to the signal strength of the samples
	if(prerun_fixed_mat == TRUE){
	sampleLik <- -1*((sampleEst_list[[1]] - (t(prot_abund[1,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec_list[[1]])
	for(c in 2:n_c){
		sampleLik <- sampleLik - ((sampleEst_list[[c]] - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec_list[[c]])
			}
		}else{
			for(c in 1:n_c){
				sampleEst <- uniquePepMean[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse
				samplePrec <- 0.5*(uniquePepPrecision[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse)
				if(c == 1){
					sampleLik <- -1*((sampleEst - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec)
					}else{
						sampleLik <- sampleLik - ((sampleEst - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec)
						}
				}
			}
		
	#adjust for prior	
	sampleLik <- sampleLik + prior_mat_likadj
	
	relLik <- sampleLik - apply(sampleLik, 1, max_non_zero) %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	relLik[prior_mat_sparse] <- exp(relLik[prior_mat_sparse])
	
	liksums <- apply(relLik, 1, sum)
	
	mixing_fract <- 1/liksums %*% t(rep(1, times = n_pp)) * prior_mat_sparse * relLik
	
	#update complete data log-likelihood
	
	new_log_lik <- logL(prot_abund, uniquePepMean, mixing_fract, uniquePepPrecision)
	
	#check for convergence
	
	if(abs(new_log_lik - previous_it) < 0.1){
		continue_it <- FALSE
		print("done")
		}else{
			whole_data_logL <- c(whole_data_logL, new_log_lik)
			print(paste("iteration", t, "log likelihood", round(new_log_lik, 3)))
			t <- t + 1
			previous_it <- new_log_lik
			}
	}

initial_convergence <- t
save(whole_data_logL, initial_convergence,  mixing_fract, prot_abund, prot_prec, file = "EM1_results.Rdata")
load("EM1_results.Rdata")
#rewrite prior sparse matrix of allowable states to either fully attribute a peptide to canonical protein trends or to the divergent trends

max_state <- apply(mixing_fract, 1, which.max)
div_max <- n_prot < max_state

#number of peptides not-conforming to some general protein trend
table(div_max)

#how many peptides are informing a protein trend, with only the largest effect considered
table(table(max_state))

for(p in 1:n_p){
	if(div_max[p] == TRUE){
		prior_mat_sparse[p, -(p + n_prot)] <- FALSE
		prior_mat_likadj[p, -(p + n_prot)] <- 0
		}else{
			prior_mat_sparse[p, p + n_prot] <- FALSE
			prior_mat_likadj[p, p + n_prot] <- 0
			}
	}



##### Run the EM again with the new sparse prior
print("EM 2 running")

if(prerun_fixed_mat == TRUE){
#lower processor usage, higher memory
sampleEst_list <- list()
samplePrec_list <- list()
for(c in 1:n_c){
	sampleEst_list[[c]] <- uniquePepMean[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	samplePrec_list[[c]] <- 0.5*(uniquePepPrecision[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse)
	print(c)
	}
	}

previous_it <- -Inf
continue_it <- TRUE

while(continue_it){
	
	#update protein abundances: using integrated likelihood
	#update protein means
	prot_abund = ((uniquePepMean*uniquePepPrecision) %*% mixing_fract)/(uniquePepPrecision %*% mixing_fract)
	prot_abund <- as.matrix(Matrix(prot_abund, sparse = FALSE))
	prot_abund[is.nan(prot_abund)] <- 0
	
	#update protein precision
	prot_prec <- uniquePepPrecision %*% mixing_fract
	
	#update mixing fract
	
	#evaluate the log-likelihood of sample abundances about the inferred protein abundances with a precision-specific to the signal strength of the samples
	if(prerun_fixed_mat == TRUE){
	sampleLik <- -1*((sampleEst_list[[1]] - (t(prot_abund[1,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec_list[[1]])
	for(c in 2:n_c){
		sampleLik <- sampleLik - ((sampleEst_list[[c]] - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec_list[[c]])
			}
		}else{
			for(c in 1:n_c){
				sampleEst <- uniquePepMean[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse
				samplePrec <- 0.5*(uniquePepPrecision[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse)
				if(c == 1){
					sampleLik <- -1*((sampleEst - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec)
					}else{
						sampleLik <- sampleLik - ((sampleEst - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec)
						}
				}
			}
	#adjust likelihood 0 values so that they aren't consumed by the sparse matrix
	sampleLik[sampleLik == 0 & prior_mat_sparse] <- -0.001
	
	#adjust for prior	
	sampleLik <- sampleLik + prior_mat_likadj
	
	relLik <- sampleLik - apply(sampleLik, 1, max_non_zero) %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	relLik[prior_mat_sparse] <- exp(relLik[prior_mat_sparse])
	
	liksums <- apply(relLik, 1, sum)
	
	mixing_fract <- 1/liksums %*% t(rep(1, times = n_pp)) * prior_mat_sparse * relLik
	
	#update complete data log-likelihood
	
	new_log_lik <- logL(prot_abund, uniquePepMean, mixing_fract, uniquePepPrecision)
	
	#check for convergence
	
	if(abs(new_log_lik - previous_it) < 0.01){
		continue_it <- FALSE
		print("done")
		}else{
			whole_data_logL <- c(whole_data_logL, new_log_lik)
			print(paste("iteration", t, "log likelihood", round(new_log_lik, 3)))
			t <- t + 1
			previous_it <- new_log_lik
			}
	}


plot(whole_data_logL, col = ifelse(c(1:length(whole_data_logL)) <= initial_convergence, "ORANGE", "BLACK"), pch = 16)
plot(whole_data_logL[-1], col = ifelse(c(2:length(whole_data_logL)) <= initial_convergence, "ORANGE", "BLACK"), pch = 16)

prot_abund_final <- prot_abund[,1:n_prot]
prot_abund_final[prot_abund_final == 0] <- NA

save(prot_abund_final, prot_prec, mixing_fract, whole_data_logL, initial_convergence, file = "EMoutput.Rdata")

pdf(file = "proteinHeat.pdf")
heatmap.2(t(prot_abund_final[,1:n_prot]), trace = "none", Colv = NULL, dendrogram = "row", na.color = "white", col = blue2red(500), labRow = FALSE)
dev.off()


library(missMDA)
#determine how many significant principal components should be included based on repeated random sub-sampling validation
pcrange <- c(2,12)
npc.compare <- estim_ncpPCA(prot_abund_final, ncp.min = pcrange[1], ncp.max = pcrange[2], method.cv = "Kfold", pNA = 0.10, nbsim = 10)
npc <- npc.compare$ncp
plot(npc.compare$criterion ~ c(pcrange[1]:pcrange[2]), pch = 16, ylab = "MS error of prediction", xlab = "number of PCs")
abline(v = npc, col = "RED", lwd = 2)

#determine the most likely values of the missing data
impute_abund <- imputePCA(prot_abund_final, npc, scale = FALSE)$completeObs
impute_abund_thresh <- impute_abund
impute_abund_thresh[impute_abund_thresh > 5] <- 5; impute_abund_thresh[impute_abund_thresh < -5] <- -5
heatmap.2(t(impute_abund_thresh), trace = "none", Colv = NULL, dendrogram = "row", na.color = "white", col = blue2red(500), labRow = FALSE, symkey = TRUE, scale = "none", denscol = "black", breaks = seq(-1*max(range(impute_abund)), max(range(impute_abund)), by = max(range(impute_abund))/250))
pdf(file = "proteinHM.pdf")
heatmap.2(t(impute_abund_thresh), trace = "none", Colv = NULL, dendrogram = "row", na.color = "white", col = blue2red(500), labRow = FALSE, symkey = TRUE, scale = "none", denscol = "black", breaks = seq(-1*max(range(impute_abund)), max(range(impute_abund)), by = max(range(impute_abund))/250))
dev.off()



#determne distributions of missing values and factor parameters
multi_impute_abund <- MIPCA(prot_abund_final, npc, nboot = 100)
plot(multi_impute_abund)





library(ggplot2)

plot_protein_add <- function(prot, possibleMap, prot_abund, uniquePepMean, uniquePepPrecision){
row_match <- c(1:n_p)[possibleMap[,prot] != 0]
div_match <- row_match[diag(mixing_fract[row_match, n_prot + row_match]) == 1]
mixed_match <- row_match[!(row_match %in% div_match)]

prot_abundM <- cbind(prot_abund[,prot], uniquePepMean[,c(mixed_match,div_match)])
prot_precM <- cbind(prot_prec[,prot], uniquePepPrecision[,c(mixed_match, div_match)])
prot_abundM[prot_abundM == 0] <- NA; prot_precM[is.na(prot_abundM)] <- NA
prot_max <- prot_abundM + 2*sqrt(1/prot_precM)
prot_min <- prot_abundM - 2*sqrt(1/prot_precM)

n_pp_plot <- length(prot_abundM[1,])
n_c <- length(prot_abundM[,1])

ncolors <- 100
color_fac <- round(mixing_fract[mixed_match, prot], length(ncolors)+1)
match_cols <- c(rainbow(ncolors + 1, start = 0, end = 0.8)[unique(color_fac)*ncolors + 1], "black", "darkgray")
match_fact <- c(paste(unique(color_fac), "x", " match", sep = ""), "protein", "non-match")
all_cols <- c("protein", paste(color_fac, "x", " match", sep = ""), rep("non-match", times = length(div_match)))
prot_cols <- apply((rep(1, times = n_c) %*% t(c(1:n_pp_plot))), c(1,2), function(col){all_cols[col]})

xpos <- (1:n_c) %*% t(rep(1, n_pp_plot)) + rep(1, n_c) %*% t(((1:n_pp_plot)-1)*n_c + ((1:n_pp_plot)-1)*2)
all_data <- data.frame(abund = c(t(prot_abundM)), abundMax = c(t(prot_max)), abundMin = c(t(prot_min)), xpos = c(xpos), color_fact = c(t(prot_cols)))

plotter <- ggplot(all_data, aes(xpos, abund, ymin = abundMin, ymax = abundMax, colour = color_fact))
plotter <- plotter + opts(panel_background = theme_rect(colour = "pink")) + xlab("sample * peptide") + ylab("relative abundance") 
plotter <- plotter + labs(colour = paste(colnames(possibleMap)[prot], ": protein ", prot, sep = ""))
print(plotter + geom_linerange()
+ scale_colour_manual(values = match_cols, limits = match_fact)
+ opts(axis.text.y = theme_text(colour = "red")))
}


pdf("prot_output_plots.pdf")
for(p in 1:100){
	plot_protein_add(p, possibleMap, prot_abund, uniquePepMean, uniquePepPrecision)
	}
dev.off()



#determine the abundance of transcripts corresponding to ascertained proteins


colnames(impute_abund_thresh)


transcript_brauer <- read.delim("../brauer-microarray/Brauer_2008.pcl")
rownames(transcript_brauer) <- transcript_brauer$SYSTEMATIC_NAME
transcript_brauer <- transcript_brauer[-1,-c(1:3)]

transcript.condition <- as.data.frame(matrix(NA, ncol = 36, nrow = 2))
colnames(transcript.condition) <- colnames(transcript_brauer)
rownames(transcript.condition) <- c("limitation", "GR")
limitations <- c("C", "N", "P", "S", "L", "U")
transcript.condition[1,] <- rep(limitations, each = 6)
transcript.condition[2,] <- rep(c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30), times = 6)

barplot(table(rowSums(is.na(transcript_brauer))))

#filter genes with many missing values
transcript_brauer <- transcript_brauer[rowSums(is.na(transcript_brauer)) < length(transcript.condition[1,])*0.2,]

pcrange <- c(16, 22)
npc.compare.trans <- estim_ncpPCA(transcript_brauer, ncp.min = pcrange[1], ncp.max = pcrange[2], method.cv = "Kfold", pNA = 0.10, nbsim = 5)
npc.trans <- npc.compare$ncp
plot(npc.compare.trans$criterion ~ c(pcrange[1]:pcrange[2]), pch = 16, ylab = "MS error of prediction", xlab = "number of PCs")
abline(v = npc.trans, col = "RED", lwd = 2)


#impute missing values (very few)


imputePCA(prot_abund_final, npc, scale = FALSE)


#missing mappings
missing_on_array <- colnames(impute_abund_thresh)[!(colnames(impute_abund_thresh) %in% rownames(transcript_brauer))]
length(missing_on_array)
genes_to_compare <- colnames(impute_abund)[colnames(impute_abund) %in% rownames(transcript_brauer)]

#determine the relative abundance of array data by looking at the corresponding condition and imputing unobserved growth rates by drawing a line between the abundance of flanking conditions

tmp <- sapply(genes_to_compare, function(match){
	transcript_brauer[rownames(transcript_brauer) == match,]
	})
transcript_brauer_reduced <- matrix(unlist(tmp), ncol = length(tmp[1,]), nrow = length(tmp[,1]))
colnames(transcript_brauer_reduced) <- colnames(tmp); rownames(transcript_brauer_reduced) <- rownames(tmp)


#load true dilution rates for each of the proteomics samples
load('~/Desktop/Composition/RNA_abundance/RNAabundance.R')

prot_cond <- as.data.frame(matrix(NA, ncol = n_c, nrow = 2))
prot_cond <- sapply(rownames(impute_abund), function(cond){
	c(unlist(strsplit(cond, '[0-9]+'))[1], unlist(strsplit(cond, '[A-Z]'))[2])
	})
colnames(prot_cond) <- rownames(impute_abund)
rownames(prot_cond) <- c("limitation", "GR")
prot_cond <- as.data.frame(cbind(t(prot_cond), realDR = NA), stringsAsFactors = FALSE)
prot_cond$realDR <- sapply(1:length(prot_cond[,1]), function(cond){
	RNAabund$actual.dr[grep(paste(prot_cond$limitation[cond], prot_cond$GR[cond], sep = ""), RNAabund$condition, ignore.case = TRUE)[1]]
	})

DR_change_mat <- matrix(0, nrow = length(transcript.condition[1,]), ncol = length(prot_cond[,1]))
colnames(DR_change_mat) <- rownames(prot_cond); rownames(DR_change_mat) <- colnames(transcript.condition)
for(cond in 1:length(prot_cond[,1])){
	#find the 2 closest DR within the same limitation
	c_match <- c(1:length(transcript.condition[1,]))[transcript.condition[1,] %in% prot_cond[cond,]$limitation]
	flanking_match <- c_match[order(abs(as.numeric(transcript.condition[2,c_match]) - prot_cond[cond,]$realDR))[1:2]]
	lb_diff <- (prot_cond$realDR[cond] - as.numeric(transcript.condition[2,flanking_match])[1])/diff(as.numeric(transcript.condition[2,flanking_match]))
	DR_change_mat[flanking_match,cond] <- c((1-lb_diff), lb_diff)
	}
	
	
remapped_transc <- t(transcript_brauer_reduced) %*% DR_change_mat

heatmap.2(remapped_transc, trace = "none", Colv = NULL, dendrogram = "row", na.rm = TRUE, na.color = "white", col = blue2red(500), labRow = FALSE, symkey = TRUE, scale = "none", denscol = "black", breaks = seq(-1*max(range(impute_abund)), max(range(impute_abund)), by = max(range(impute_abund))/250))



impute_abund[,(colnames(impute_abund) %in% rownames(transcript_brauer))] 







