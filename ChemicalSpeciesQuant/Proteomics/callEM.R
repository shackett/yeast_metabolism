setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
source("pep_library.R")
matrix_fxn()
load("EMimport.Rdata")

#print(sessionInfo())

tmp <- Matrix(diag(rep(prior_p_div, times = n_p))); colnames(tmp) <- paste(unique_pepNames, "divergent", sep = "_")
mixing_fract <- unique_mappingMat * (1/(t(t((1-prior_p_div)^-1*rowSums(unique_mappingMat)))) %*% rep(1, n_prot))
mixing_fract <- cBind(mixing_fract, tmp)
prior_mat_likadj <- unique_mappingMat/((1-prior_p_div)^-1)
prior_mat_likadj <- cBind(prior_mat_likadj, tmp)

tmp <- as.matrix(tmp)
tmp[!(tmp %in% c(0,1))] <- 1 
prior_mat_logical <- cbind(as.matrix(unique_mappingMat), tmp) 
prior_mat_logical <- prior_mat_logical == 1
prior_mat_sparse <- Matrix(prior_mat_logical)
#prior_mat_likadj2 <- prior_mat_likadj*prior_mat_sparse
prior_mat_likadj[prior_mat_sparse] <- log(prior_mat_likadj[prior_mat_sparse])
colnames(prior_mat_sparse) <- colnames(mixing_fract) <- colnames(prot_abund)
rm(tmp, prior_mat_logical); gc(); gc()


#preprocess sparse data and precision matrices

if(prerun_fixed_mat == TRUE){
#lower processor usage, higher memory
sampleEst_list <- list()
samplePrec_list <- list()
for(c in 1:n_c){
	sampleEst_list[[c]] <- uniquePepMean[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	samplePrec_list[[c]] <- 0.5*(uniquePepPrecision[c,] %*% t(rep(1, times = n_pp)) * prior_mat_sparse)
	colnames(sampleEst_list[[c]]) <- colnames(prot_abund)
	colnames(samplePrec_list[[c]]) <- colnames(prot_abund)
	print(c)
	}
	};gc()
		

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
	
	#correct for proteins that only have peptides that are also found in other proteins: penalzie the log-likelihood of the peptides associated with these proteins evenly
	weight_store <- list()	
	for(ap in ambigprots){
		tmp <- sampleLik[,ap][prior_mat_sparse[,ap]]
		weights <- 1/length(tmp)
		weight_store[[ap]] <- weights
		sampleLik[,ap][prior_mat_sparse[,ap]] <- tmp + log(prior_p_div)*weights
		}
	
	#adjust for prior
	colnames(sampleLik) <- colnames(prior_mat_likadj)
	sampleLik <- sampleLik + prior_mat_likadj
	
	#if an ambiguous protein has the highest likelihood (besides divergent peptides), than set the divergent peptide to log(prior_p_div) + penalty applied to sampleLIk
	
	if(!unq_matches_only){
	# record the likelihood of all peptide-protein matches, set all non-matches to NA
	tmp <- sampleLik[,1:n_prot]; tmp[!prior_mat_sparse[,1:n_prot]] <- NA
	
	# does the ambiguous protein have the highest likelihood match to some peptides
	ambig_good_fit <- is.finite(apply(tmp[,ambigprots], 1, function(x){if(length(x[!is.na(x)]) == 0){-Inf}else{max(x, na.rm = TRUE)}})) & (apply(tmp[,ambigprots], 1, function(x){if(length(x[!is.na(x)]) == 0){-Inf}else{max(x, na.rm = TRUE)}}) == apply(tmp, 1, max, na.rm = TRUE))
	good_fit_peps <- c(1:n_p)[ambig_good_fit]
	good_fit_match <- ambigprots[apply(tmp[ambig_good_fit,ambigprots], 1, which.max)]
	
	#penalize the divergent pep priors
	diag(sampleLik[good_fit_peps, n_prot + good_fit_peps]) <- diag(sampleLik[good_fit_peps, n_prot + good_fit_peps]) + log(prior_p_div)*sapply(good_fit_match, function(x){weight_store[[x]]})
	
		}
	
	
	# a peptide is attributed to a protein proportionally to its likelihood of that match - ideally it would fit the residuals in an iterative fashion (gibbs sampling), but it shouldn't be crucial.
	
	sampleLik_NA <- sampleLik
	sampleLik_NA[!prior_mat_sparse] <- NA
	relLik <- sampleLik - apply(sampleLik_NA, 1, max_non_NA) %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	relLik[prior_mat_sparse] <- exp(relLik[prior_mat_sparse])
	
	liksums <- apply(relLik, 1, sum)
	
	mixing_fract <- 1/liksums %*% t(rep(1, times = n_pp)) * prior_mat_sparse * relLik
	
	#update complete data log-likelihood
	
	new_log_lik <- logL(prot_abund, uniquePepMean, mixing_fract, uniquePepPrecision)
	
	#check for convergence
	
	if(abs(new_log_lik - previous_it) < 0.01 | (t > 30)){
		continue_it <- FALSE
		whole_data_logL <- c(whole_data_logL, new_log_lik)
		t <- t + 1
		print("done")
		}else{
			whole_data_logL <- c(whole_data_logL, new_log_lik)
			print(paste("iteration", t, "log likelihood", round(new_log_lik, 3)))
			t <- t + 1
			previous_it <- new_log_lik
			}
	}

#for each peptide compare a model where peptide patterns are equivalent to protein patterns with one where a peptide doesn't match - the divergent model will fit perfectly vs. the alternative where the protein fits

likPmatch <- -1*uniquePepPrecision*(uniquePepMean - prot_abund %*% t(mixing_fract))^2
likPmatch <- apply(likPmatch, 2, sum)

divergentPep_summary <- data.frame(peptide = names(likPmatch)[likPmatch < log(prior_p_div)], prot_matches = NA, likelihood = unname(likPmatch[likPmatch < log(prior_p_div)]), nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)

poor_match_redMat <- prior_mat_sparse[likPmatch < log(prior_p_div), 1:n_prot]

for(pep in 1:length(poor_match_redMat[,1])){
	divergentPep_summary$prot_matches[pep] <- paste(colnames(poor_match_redMat)[poor_match_redMat[pep,]], collapse = "/")
	divergentPep_summary$nS[pep] <- length(grep("S", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	divergentPep_summary$nR[pep] <- length(grep("R", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	divergentPep_summary$nY[pep] <- length(grep("Y", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	}
divergentPep_summary$phosphoSite <- divergentPep_summary$nS != 0 | divergentPep_summary$nR != 0 | divergentPep_summary$nY != 0


background_SRYfreq <- data.frame(peptide = names(likPmatch)[likPmatch >= log(prior_p_div)], prot_matches = NA, likelihood = unname(likPmatch[likPmatch >= log(prior_p_div)]), nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)

good_match_redMat <- prior_mat_sparse[likPmatch >= log(prior_p_div), 1:n_prot]

for(pep in 1:length(good_match_redMat[,1])){
	background_SRYfreq$prot_matches[pep] <- paste(colnames(good_match_redMat)[good_match_redMat[pep,]], collapse = "/")
	background_SRYfreq$nS[pep] <- length(grep("S", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	background_SRYfreq$nR[pep] <- length(grep("R", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	background_SRYfreq$nY[pep] <- length(grep("Y", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	}
background_SRYfreq$phosphoSite <- background_SRYfreq$nS != 0 | background_SRYfreq$nR != 0 | background_SRYfreq$nY != 0


#background_SRYfreq <- data.frame(nS = rep(0, 15), nR = rep(0, 15), nY = rep(0, 15), all = rep(0, 15))
#rownames(background_SRYfreq) <- 0:14

#for(pep in names(likPmatch)){
#	background_SRYfreq$nS[length(grep("S", unlist(strsplit(pep, ""))))] <- background_SRYfreq$nS[length(grep("S", unlist(strsplit(pep, ""))))] + 1
#	background_SRYfreq$nR[length(grep("R", unlist(strsplit(pep, ""))))] <- background_SRYfreq$nR[length(grep("R", unlist(strsplit(pep, ""))))] + 1
#	background_SRYfreq$nY[length(grep("Y", unlist(strsplit(pep, ""))))] <- background_SRYfreq$nY[length(grep("Y", unlist(strsplit(pep, ""))))] + 1
#	background_SRYfreq$all[length(grep("S", unlist(strsplit(pep, "")))) + length(grep("R", unlist(strsplit(pep, "")))) + length(grep("Y", unlist(strsplit(pep, ""))))] <- background_SRYfreq$all[length(grep("S", unlist(strsplit(pep, "")))) + length(grep("R", unlist(strsplit(pep, "")))) + length(grep("Y", unlist(strsplit(pep, ""))))] + 1
#	}


hist(likPmatch, xlab = "Peptide log-likelihood", col = "orangered")
abline(v = log(prior_p_div), lwd = 3)


initial_convergence <- t - 1
#rewrite prior sparse matrix of allowable states to either fully attribute a peptide to canonical protein trends or to the divergent trends

max_state <- apply(mixing_fract, 1, which.max)
div_max <- apply(mixing_fract[,1:n_prot], 1, sum) <= 0.5
#div_max <- n_prot < max_state

#number of peptides not-conforming to some general protein trend
table(div_max)

#how many peptides are informing a protein trend, with only the largest effect considered
barplot(table(table(max_state[max_state <= n_prot])))

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
	prot_prec[prot_prec < 1] <- 0
	
	#update mixing fract
	
	#evaluate the log-likelihood of sample abundances about the inferred protein abundances with a precision-specific to the signal strength of the samples
	if(prerun_fixed_mat == TRUE){
	sampleLik <- - ((sampleEst_list[[1]] - (t(prot_abund[1,] %*% t(rep(1, times = n_p)))*prior_mat_sparse))^2*samplePrec_list[[1]])
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
	
	#correct for proteins that only have peptides that are also found in other proteins: penalzie the log-likelihood of the peptides associated with these proteins evenly
	weight_store <- list()	
	for(ap in ambigprots){
		tmp <- sampleLik[,ap][prior_mat_sparse[,ap]]
		weights <- 1/length(tmp)
		weight_store[[ap]] <- weights
		sampleLik[,ap][prior_mat_sparse[,ap]] <- tmp + log(prior_p_div)*weights
		}
	
	#adjust for prior	
	sampleLik <- prior_mat_likadj + sampleLik
	
	sampleLik_NA <- sampleLik
	sampleLik_NA[!prior_mat_sparse] <- NA
	relLik <- sampleLik - apply(sampleLik_NA, 1, max_non_NA) %*% t(rep(1, times = n_pp)) * prior_mat_sparse
	relLik[prior_mat_sparse] <- exp(relLik[prior_mat_sparse])
	
	liksums <- apply(relLik, 1, sum)
	
	mixing_fract <- 1/liksums %*% t(rep(1, times = n_pp)) * prior_mat_sparse * relLik
	
	#update complete data log-likelihood
	
	new_log_lik <- logL(prot_abund, uniquePepMean, mixing_fract, uniquePepPrecision)
	
	#check for convergence
	
	if(abs(new_log_lik - previous_it) < 0.01){
		continue_it <- FALSE
		whole_data_logL <- c(whole_data_logL, new_log_lik)
		print("done")
		}else{
			whole_data_logL <- c(whole_data_logL, new_log_lik)
			print(paste("iteration", t, "log likelihood", round(new_log_lik, 3)))
			t <- t + 1
			previous_it <- new_log_lik
			}
	}


#prot_abund_final <- prot_abund[,1:n_prot]
prot_abund_final[prot_abund_final == 0] <- NA

#Remove proteins with less than 0.5 peptides attributed to them
n_prot_nonmatched <- sum(apply(mixing_fract[,1:n_prot], 2, sum) <= 0.5)

prot_abund_final <- prot_abund[,c(1:n_prot)[apply(mixing_fract[,1:n_prot], 2, sum) > 0.5]]
prot_abund_final[prot_abund_final == 0] <- NA


if(unq_matches_only == TRUE){
	save(prot_abund_final, prot_prec, mixing_fract, whole_data_logL, initial_convergence, max_state, div_max, n_prot_nonmatched, div_max, divergentPep_summary, background_SRYfreq, file = "EMoutputUnq.Rdata")
	}else{
		save(prot_abund_final, prot_prec, mixing_fract, whole_data_logL, initial_convergence, max_state, div_max, n_prot_nonmatched, div_max, divergentPep_summary, background_SRYfreq, file = "EMoutputDeg.Rdata")
	}


