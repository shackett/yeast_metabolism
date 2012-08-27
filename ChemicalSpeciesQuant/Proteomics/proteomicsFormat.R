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

#determine the expectation of the standard deviation as a heteroschedastic fxn of IC using p0.05 light v. p0.05 heavy

avgSignalSTD <- apply(cbind(heavyIC[,colnames(heavyIC) == "P0.05"], lightIC[,colnames(lightIC) == "P0.05"]), 1, mean)
logLight <- log2(lightIC[is.finite(avgSignalSTD),colnames(lightIC) == "P0.05"])
logHeavy <- log2(heavyIC[is.finite(avgSignalSTD),colnames(heavyIC) == "P0.05"])
avgSignalSTD <- avgSignalSTD[is.finite(avgSignalSTD)]


logLight <- logLight + optimize(normFactor, c(-1, 1), logLight = logLight, logHeavy = logHeavy)$minimum

#variance of the replicate differences accounting for small sample (scaling factor of 2)
STDvar <- (logLight - logHeavy)^2*2
STDvar_fit <- lm(log2(STDvar) ~ log2(avgSignalSTD))$coef


plot(log2(STDvar) ~ log2(avgSignalSTD), pch = 16, cex = 0.3)
abline(STDvar_fit, col = "RED")

gplot.hexbin(hexbin(log2(avgSignalSTD), log2(STDvar), xbins = 80), colramp = rainbow)

gplot.hexbin(hexbin(logLight, logHeavy, xbins = 200), colramp = rainbow)


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

uniquePepMean <- ((t(abundMat) * fittedPrec) %*% pepToUniq)/(fittedPrec %*% pepToUniq)
uniquePepPrecision <- fittedPrec %*% pepToUniq

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


prot_abund_final <- prot_abund[,1:n_prot]
prot_abund_final[prot_abund_final == 0] <- NA

save(prot_abund_final, prot_prec, mixing_fract, initial_convergence, final_convergence, file = "EMoutput.Rdata")

#impute missing values based on the inverted-Wishart model of the norm package - use EM, also could do multiple imputations via MCMC

#tmp <- prelim.norm(t(prot_abund_final))
#theta_hat <- em.norm(tmp, criterion = 10^-8)
#rngseed(1234567)	#set random number generator seed 
#par_est2 <- getparam.norm(tmp, par_est)
#imp_data <- imp.norm(tmp, theta_hat, t(prot_abund_final))

pdf(file = "proteinHeat.pdf")
heatmap.2(t(prot_abund_final[,1:n_prot]), trace = "none", Colv = NULL, dendrogram = "row", na.color = "white", col = blue2red(500))
dev.off()


#data(mdata)
#s <- prelim.norm(mdata)	#do preliminary manipulations
#thetahat <- em.norm(s)	#compute mle 
#getparam.norm(s,thetahat,corr=TRUE)$r #look at estimated correlations

#data(mdata)
#s <- prelim.norm(mdata)
#thetahat <- em.norm(s, criterion = 10^-8)
#rngseed(1234567)	#set random number generator seed 
#ximp <- imp.norm(s,thetahat,mdata) #impute missing data under the MLE














