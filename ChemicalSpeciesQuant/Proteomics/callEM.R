setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
source("pep_library.R")
matrix_fxn()
load("EMimport.Rdata")

#print(sessionInfo())

#initialize theta (mixing_fraction) and pi (prob non-divergent peptide)

mixing_fract <- unique_mappingMat/rowSums(unique_mappingMat)
sparse_mapping <- Matrix(unique_mappingMat)

pi_fit <- rep(1, times = n_p)
#record the log-likelihood of a match versus non-match
pi_lik_comp <- data.frame(match = rep(NA, times = n_p), divergent = rep(NA, times = n_p))

#preprocess sparse data and precision matrices

sparsePrecision <- Matrix(uniquePepPrecision != 0)

if(prerun_fixed_mat == TRUE){
#lower processor usage, higher memory
sampleEst_list <- list()
samplePrec_list <- list()
for(c in 1:n_c){
	sampleEst_list[[c]] <- uniquePepMean[c,] %*% t(rep(1, times = n_prot)) * sparse_mapping
	samplePrec_list[[c]] <- 0.5*(uniquePepPrecision[c,] %*% t(rep(1, times = n_prot)) * sparse_mapping)
	colnames(sampleEst_list[[c]]) <- colnames(mixing_fract)
	colnames(samplePrec_list[[c]]) <- colnames(mixing_fract)
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
	#update protein means (c x j)
	prot_abund = ((uniquePepMean*uniquePepPrecision) %*% mixing_fract)/(uniquePepPrecision %*% mixing_fract)
	prot_abund <- as.matrix(Matrix(prot_abund, sparse = FALSE))
	prot_abund[is.nan(prot_abund)] <- 0
	
	#update protein precision
	prot_prec <- uniquePepPrecision %*% mixing_fract
	
	#update pi - prob peptide match
  
  prot_matchLik <- -(1/2)*uniquePepPrecision*(((prot_abund %*% t(mixing_fract)) - uniquePepMean)^2)
	pi_lik_comp$match <- colSums(prot_matchLik); pi_lik_comp$divergent <- log(prior_p_div)
  #table(apply(pi_lik_comp, 1, diff) < 0)
	pi_fit <- ifelse(apply(pi_lik_comp, 1, diff) < 0, 1, 0)
  
  #update mixing fract
	
	#evaluate the log-likelihood of sample abundances about the inferred protein abundances with a precision-specific to the signal strength of the samples
	if(prerun_fixed_mat == TRUE){
	sampleLik <- -1*((sampleEst_list[[1]] - (t(prot_abund[1,] %*% t(rep(1, times = n_p)))*sparse_mapping))^2*samplePrec_list[[1]])
	for(c in 2:n_c){
		sampleLik <- sampleLik - ((sampleEst_list[[c]] - (t(prot_abund[c,] %*% t(rep(1, times = n_p)))*sparse_mapping))^2*samplePrec_list[[c]])
      }
		}
	
  
  
  
	#correct for subsumable proteins (only have peptides that are also found in other proteins): penalize the log-likelihood of the peptides associated with these proteins evenly
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
	
	new_log_lik <- sum(apply(log(mixing_fract) + sampleLik_NA, 1, max, na.rm = TRUE))
	
	
	#check for convergence
	
	if(new_log_lik - previous_it < 0.01 | (t > 30)){
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

initial_convergence <- t - 1
#rewrite prior sparse matrix of allowable states to either fully attribute a peptide to canonical protein trends or to the divergent trends

max_state <- apply(mixing_fract, 1, which.max)
div_max <- apply(mixing_fract[,1:n_prot], 1, sum) <= 0.5

#number of peptides not-conforming to some general protein trend
table(div_max)


# for each model compare the best peptide-protein match with the non-match likelihood*prior

likdiff <- apply(sampleLik_NA[,1:n_prot], 1, max, na.rm = TRUE) - diag(sampleLik_NA[1:n_p,(n_prot+1):n_pp])

likdiff_df <- data.frame(likelihood = likdiff, matched = ifelse(likdiff >= 0, "Protein-match", "Divergent-trend"))
likdiff_df$likelihood[likdiff_df$likelihood <= -400] <- -400

likdiff_plot <- ggplot(likdiff_df, aes(x = likelihood, fill = matched))
likdiff_plot + geom_histogram()


#for each peptide compare a model where peptide patterns are equivalent to protein patterns with one where a peptide doesn't match - the divergent model will fit perfectly vs. the alternative where the protein fits


divergentPep_summary <- data.frame(peptide = names(likdiff)[div_max], prot_matches = NA, likelihood = unname(likdiff[div_max]), nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)

poor_match_redMat <- prior_mat_sparse[div_max, 1:n_prot]

for(pep in 1:length(poor_match_redMat[,1])){
	divergentPep_summary$prot_matches[pep] <- paste(colnames(poor_match_redMat)[poor_match_redMat[pep,]], collapse = "/")
	divergentPep_summary$nS[pep] <- length(grep("S", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	divergentPep_summary$nR[pep] <- length(grep("R", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	divergentPep_summary$nY[pep] <- length(grep("Y", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	}
divergentPep_summary$possible_phosphoSite <- divergentPep_summary$nS != 0 | divergentPep_summary$nR != 0 | divergentPep_summary$nY != 0


background_SRYfreq <- data.frame(peptide = names(likdiff)[!div_max], prot_matches = NA, likelihood = unname(likdiff[!div_max]), nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)

good_match_redMat <- prior_mat_sparse[!div_max, 1:n_prot]

for(pep in 1:length(good_match_redMat[,1])){
	background_SRYfreq$prot_matches[pep] <- paste(colnames(good_match_redMat)[good_match_redMat[pep,]], collapse = "/")
	background_SRYfreq$nS[pep] <- length(grep("S", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	background_SRYfreq$nR[pep] <- length(grep("R", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	background_SRYfreq$nY[pep] <- length(grep("Y", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	}
background_SRYfreq$possible_phosphoSite <- background_SRYfreq$nS != 0 | background_SRYfreq$nR != 0 | background_SRYfreq$nY != 0





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
	
	new_log_lik <- sum(apply(log(mixing_fract) + sampleLik_NA, 1, max, na.rm = TRUE))
	
	#check for convergence
	
	if(new_log_lik - previous_it < 0.01){
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


if(unq_matches_only == TRUE){
	save(prot_abund, prot_prec, mixing_fract, whole_data_logL, initial_convergence, max_state, div_max, div_max, divergentPep_summary, background_SRYfreq, likdiff_df, file = "EMoutputUnq.Rdata")
	}else{
		save(prot_abund, prot_prec, mixing_fract, whole_data_logL, initial_convergence, max_state, div_max, div_max, divergentPep_summary, background_SRYfreq, likdiff_df, file = "EMoutputDeg.Rdata")
	}


