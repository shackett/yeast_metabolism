
setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
load('20120807ProtPepMatrices.Rdata')

####### An EM to go from peptides to proteins #######


# Filter to only look at peptides with no missing values

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


library(hexbin)
plot(log2(STDvar) ~ log2(avgSignalSTD), pch = 16, cex = 0.3)
abline(STDvar_fit, col = "RED")

gplot.hexbin(hexbin(log2(avgSignalSTD), log2(STDvar), xbins = 80), colramp = rainbow)

plot(logLight ~ logHeavy, pch = 16, cex = 0.3)
abline(0, 1, col = "RED")

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

#set the precision of missing values to 0

fittedPrec[t(is.na(abundMat))] <- 0
fittedPrec[is.na(fittedPrec)] <- 0
abundMat[t(fittedPrec) == 0] <- 0

uniquePepMean <- ((t(abundMat) * fittedPrec) %*% pepToUniq)/(fittedPrec %*% pepToUniq)
uniquePepPrecision <- fittedPrec %*% pepToUniq

library(gplots)
uniquePepMean[is.nan(uniquePepMean)] <- NA
table(rowSums(!is.na(t(uniquePepMean))))

#change mapping from peptides to unique peptides (averaging over ionization states)
unique_mappingMat <- t(pepToUniq) %*% mappingMat
unique_mappingMat[!(unique_mappingMat %in% c(0,1))] <- 1


######## EM ######

### Initalization ###

n_p = length(unique_pepNames) #4042
n_pp = length(unique_mappingMat[,1]) + length(unique_mappingMat[1,]) #4964
n_c = length(uniquePepMean[,1]) #15

prior_p_div <- 0.01

tmp <- diag(rep(prior_p_div, times = n_p)); colnames(tmp) <- paste(unique_pepNames, "divergent", sep = "_")

mixing_fract <- unique_mappingMat/((1-prior_p_div)^-1*rowSums(unique_mappingMat))
mixing_fract <- cbind(mixing_fract, tmp) 
prior_mat <- cbind(unique_mappingMat, tmp) 


#prot_abund <- matrix(0, ncol = length(mixing_fract[1,]), nrow = length(uniquePepMean[,1]))
#rownames(prot_abund) <- rownames(uniquePepMean); colnames(prot_abund) <- colnames(mixing_fract)

fractSVD <- svd(t(mixing_fract))
pseudoinv <- fractSVD$v%*%solve(diag(fractSVD$d))%*%t(fractSVD$u)
uniquePepMean[is.na(uniquePepMean)] <- 0
prot_abund <- uniquePepMean%*%pseudoinv
rownames(prot_abund) <- rownames(uniquePepMean); colnames(prot_abund) <- colnames(mixing_fract)
prot_prec <- uniquePepPrecision %*% mixing_fract
	

### Iteration ###




for(t in 1:20){
	
	#update protein means
	prot_abund = ((uniquePepMean*uniquePepPrecision) %*% mixing_fract)/(uniquePepPrecision %*% mixing_fract)
	prot_abund[is.nan(prot_abund)] <- 0
	
	#update protein precision
	prot_prec <- uniquePepPrecision %*% mixing_fract
	
	#update mixing fract
	
	
	library(Matrix)
	?sparseMatrix
	sparsePrior <- Matrix(prior_mat)
	
	
	
	totalLik <- matrix(0, ncol = n_p, nrow = n_pp)
	for(c in 1:length(prot_abund[,1])){
		tmp1 = ((prot_abund[c,] %*% t(rep(1, times = n_p)))*t(prior_mat) - (t(uniquePepMean[c,] %*% t(rep(1, times = n_pp))))*t(prior_mat))^2
		tmp2 = 0.5*(t(uniquePepPrecision[c,] %*% t(rep(1, times = n_pp)))*t(prior_mat))
		sampleLik = tmp1*tmp2
		totalLik <- sampleLik + totalLik
		}
	
	
	
	#update complete data log-likelihood
	
	
	
	#check for convergence
	
	
	}

logL <- function(prot_abund, uniquePepMean, mixing_fract, uniquePepPrecision){
	
	#calculate the complete data log-likelihood
	p_est <- prot_abund %*% t(mixing_fract)
	
	}




#heatmap.2(t(uniquePepMean), trace = "none")







var_calc <- function(sampleIC, STDvar_fit){
	2^(STDvar_fit[1] + STDvar_fit[2]*log2(sampleIC))
	}



normFactor <- function(alpha, logLight, logHeavy){
	#determine a factor that will be added to each logLight abundance such that the global sum of squares is minimized between logLight and logHeavy
	sum(((logLight + alpha) - logHeavy)^2)
	}


