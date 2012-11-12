setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
load('../BulkComposition/protSpecQuantOut.Rdata')
load('tmp.Rdata')

#prot_abund_final
#prot_prec
n_prot <- 859

cond_info <- data.frame(t(sapply(rownames(prot_abund_final), function(cond){
	if(strsplit(cond, "")[[1]][1] == "L"){
		tmp <- cond
		}else{
			tmp <- paste(c(tolower(strsplit(cond, "")[[1]][1]), strsplit(cond, "")[[1]][-1]), collapse = "")
			}
		t(conditions[conditions$condition == tmp,])
	})), stringsAsFactors = FALSE)
colnames(cond_info) <- colnames(conditions)	
	
cond_info$limitation <- factor(cond_info$limitation)
cond_info$actualDR <- as.numeric(cond_info$actualDR)
cond_info$logSF <- as.numeric(cond_info$logSF)

#set 
specAbund <- prot_abund_final[,1:n_prot]
specPrec <- prot_prec[,1:n_prot]
specAbund[as.matrix(specPrec == 0)] <- NA

n_bs <- 5000

p_valMat <- matrix(NA, nrow = n_prot, ncol = length(levels(cond_info$limitation))*2 + 1)
rownames(p_valMat) <- colnames(specAbund)[1:n_prot]
colnames(p_valMat) <- c("prot_abund", unlist(sapply(c("int", "slope"), function(variable){
	paste(levels(cond_info$limitation), variable, sep = "_")
	})))
t_valMat <- p_valMat


for(pr in 1:n_prot){
	
	valid.conds <- !is.na(specAbund[,pr])
	prot_val <- specAbund[,pr][valid.conds]
	prot_W <- specPrec[,pr][valid.conds]
	n_samples <- sum(valid.conds)
	
	design_mat <- cbind(cond_info$logSF, sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}), sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)})*cond_info$actualDR)[valid.conds,]
	
	prot_val_weighted <- prot_val * sqrt(prot_W)
	design_weighted <- design_mat * t(t(sqrt(prot_W))) %*% rep(1, times = length(p_valMat[1,]))
	
	beta_val <- solve(t(design_weighted) %*% design_weighted) %*% t(design_weighted) %*%  prot_val_weighted
	residuals <- (prot_val_weighted - design_weighted %*% beta_val)
	model_degree <- length(beta_val)
	RSE <- sqrt(sum(residuals^2)/(n_samples - length(beta_val)))
	beta_se <- RSE^2 * solve(t(design_weighted) %*% design_weighted) 
	
	#remove weights on residuals for now
	infl_residuals <- residuals/sqrt(prot_W) * sqrt(n_samples/(n_samples - model_degree))
	
	bs_resids <- matrix(sample(infl_residuals, n_samples*n_bs, replace = TRUE), ncol = n_samples)
	bs_resids <- bs_resids * (t(t(rep(1, times = n_bs))) %*% sqrt(prot_W))
	
	bs_response <- t(t(rep(1, times = n_bs))) %*% t(design_weighted %*% beta_val) + bs_resids
	
	bs_coefs <- sapply(c(1:n_bs), function(bsN){
	#calculate the bootstrap-beta and se(bootstrap-beta) for the bs_response and offset by the fitted beta.
	beta_bs <- solve(t(design_weighted)%*% design_weighted) %*% t(design_weighted) %*% bs_response[bsN,]
	(beta_val - beta_bs)/sqrt(diag(sum((bs_response[bsN,] - (design_weighted %*% beta_bs))^2)/(n_samples - length(beta_val)) * solve(t(design_weighted) %*% design_weighted)))
	#lm(bs_response[bsN,] ~ design_weighted + 0)
	})
	
	#beta_val %*% rep(1, times = n_bs) + bs_coefs*sqrt(diag(beta_se))
	
	### check ###
	p_valMat[pr,] <- 1 - (abs(apply((beta_val %*% rep(1, times = n_bs) + bs_coefs*sqrt(diag(beta_se))) >= 0, 1, sum) - n_bs/2)/(n_bs/2)) 
	t_valMat[pr,] <- beta_val/sqrt(diag(beta_se))
	}

FDR_desired <- 0.05
source('~/Desktop/Rabinowitz/Cancer_metabolomics/cancer_lib.R')


condDiff_FDR_stats <- data.frame(fractionNull = rep(NA, times = length(colnames(p_valMat))), testStatCO = rep(NA, times = length(colnames(p_valMat))), FDRattained = rep(NA, times = length(colnames(p_valMat))), discoveries = rep(NA, times = length(colnames(p_valMat))))
rownames(condDiff_FDR_stats) <- colnames(p_valMat)
Discoveries <- list()


for(ncomp in 1:length(condDiff_FDR_stats[,1])){
	
	pvals <- p_valMat[,ncomp]
	lambda.vals <- seq(0, 1, by = 0.001)
	num.p <- length(pvals)
	
	pi <- rep(NA, times = length(lambda.vals))
	for (l in 1:length(lambda.vals)){
		pi[l] <- sum(pvals > lambda.vals[l])/(num.p * (1 - lambda.vals[l]))
		}
	#print(plot(pi ~ lambda.vals, main = paste("estimate of fraction of TPs for", rownames(condDiff_FDR_stats)[ncomp]), ylim = c(0,1), pch = 16))}
	
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "prot_abund"] <- 0.2
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "c_int"] <- 0.5
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "n_int"] <- 0.4
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "p_int"] <- 0.9
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "c_slope"] <- 0.5
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "n_slope"] <- 0.4
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "p_slope"] <- 0.6
	
	condDiff_FDR_stats$testStatCO[ncomp] <- seq(0.0001, 1, by = 0.0001)[which.min(unlist(lapply(seq(0.0001, 1, by = 0.0001), FDR.eval, pi = condDiff_FDR_stats$fractionNull[ncomp], num.p = num.p, test.stats = pvals, FDRlev = FDR_desired, use.unif = TRUE)))]

	condDiff_FDR_stats$FDRattained[ncomp] <- FDR.calc(pi = condDiff_FDR_stats$fractionNull[ncomp], lambda = condDiff_FDR_stats$testStatCO[ncomp], num.p = num.p, test.stats = pvals, use.unif = TRUE)
	
	
	discoveries <- data.frame(protein = rownames(p_valMat)[p_valMat[,ncomp] < condDiff_FDR_stats$testStatCO[ncomp]], change = t_valMat[,ncomp][p_valMat[,ncomp] < condDiff_FDR_stats$testStatCO[ncomp]], p_value = p_valMat[,ncomp][p_valMat[,ncomp] < condDiff_FDR_stats$testStatCO[ncomp]], q_value = NA, stringsAsFactors = FALSE)
	
	discoveries$q_value <- sapply(c(1:length(discoveries[,1])), function(i){
		(condDiff_FDR_stats$fractionNull[ncomp]*discoveries$p_value[i]*num.p) / (sum(pvals <= discoveries$p_value[i]))
		})
	
	condDiff_FDR_stats$discoveries[ncomp] <- length(discoveries[,1])
	Discoveries[[ncomp]] <- discoveries
				
	}


###############

#recompute with removing the protein concentration offset to give changes in proteome composition

p_compMat <- matrix(NA, nrow = n_prot, ncol = length(levels(cond_info$limitation))*2)
rownames(p_compMat) <- colnames(specAbund)[1:n_prot]
colnames(p_compMat) <- c(unlist(sapply(c("int", "slope"), function(variable){
	paste(levels(cond_info$limitation), variable, sep = "_")
	})))
t_compMat <- p_compMat


for(pr in 1:n_prot){
	
	valid.conds <- !is.na(specAbund[,pr])
	prot_val <- (specAbund[,pr] - cond_info$logSF)[valid.conds]
	prot_W <- specPrec[,pr][valid.conds]
	n_samples <- sum(valid.conds)
	
	design_mat <- cbind(sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}), sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)})*cond_info$actualDR)[valid.conds,]
	
	prot_val_weighted <- prot_val * sqrt(prot_W)
	design_weighted <- design_mat * t(t(sqrt(prot_W))) %*% rep(1, times = length(p_compMat[1,]))
	
	beta_val <- solve(t(design_weighted) %*% design_weighted) %*% t(design_weighted) %*%  prot_val_weighted
	residuals <- (prot_val_weighted - design_weighted %*% beta_val)
	model_degree <- length(beta_val)
	RSE <- sqrt(sum(residuals^2)/(n_samples - length(beta_val)))
	beta_se <- RSE^2 * solve(t(design_weighted) %*% design_weighted) 
	
	#remove weights on residuals for now
	infl_residuals <- residuals/sqrt(prot_W) * sqrt(n_samples/(n_samples - model_degree))
	
	bs_resids <- matrix(sample(infl_residuals, n_samples*n_bs, replace = TRUE), ncol = n_samples)
	bs_resids <- bs_resids * (t(t(rep(1, times = n_bs))) %*% sqrt(prot_W))
	
	bs_response <- t(t(rep(1, times = n_bs))) %*% t(design_weighted %*% beta_val) + bs_resids
	
	bs_coefs <- sapply(c(1:n_bs), function(bsN){
	#calculate the bootstrap-beta and se(bootstrap-beta) for the bs_response and offset by the fitted beta.
	beta_bs <- solve(t(design_weighted)%*% design_weighted) %*% t(design_weighted) %*% bs_response[bsN,]
	(beta_val - beta_bs)/sqrt(diag(sum((bs_response[bsN,] - (design_weighted %*% beta_bs))^2)/(n_samples - length(beta_val)) * solve(t(design_weighted) %*% design_weighted)))
	#lm(bs_response[bsN,] ~ design_weighted + 0)
	})
	
	### check ###
	p_compMat[pr,] <- 1 - (abs(apply((beta_val %*% rep(1, times = n_bs) + bs_coefs*sqrt(diag(beta_se))) >= 0, 1, sum) - n_bs/2)/(n_bs/2)) 
	t_compMat[pr,] <- beta_val/sqrt(diag(beta_se))
	}

FDR_desired <- 0.05
source('~/Desktop/Rabinowitz/Cancer_metabolomics/cancer_lib.R')


condDiffcomp_FDR_stats <- data.frame(fractionNull = rep(NA, times = length(colnames(p_compMat))), testStatCO = rep(NA, times = length(colnames(p_compMat))), FDRattained = rep(NA, times = length(colnames(p_compMat))), discoveries = rep(NA, times = length(colnames(p_compMat))))
rownames(condDiffcomp_FDR_stats) <- colnames(p_compMat)
Discoveries <- list()


for(ncomp in 1:length(condDiffcomp_FDR_stats[,1])){
	
	pvals <- p_compMat[,ncomp]
	lambda.vals <- seq(0, 1, by = 0.001)
	num.p <- length(pvals)
	
	pi <- rep(NA, times = length(lambda.vals))
	for (l in 1:length(lambda.vals)){
		pi[l] <- sum(pvals > lambda.vals[l])/(num.p * (1 - lambda.vals[l]))
		}
	#print(plot(pi ~ lambda.vals, main = paste("estimate of fraction of TPs for", rownames(condDiffcomp_FDR_stats)[ncomp]), ylim = c(0,1), pch = 16))}
	
	condDiffcomp_FDR_stats$fractionNull[rownames(condDiffcomp_FDR_stats) == "c_int"] <- 0.3
	condDiffcomp_FDR_stats$fractionNull[rownames(condDiffcomp_FDR_stats) == "n_int"] <- 0.2
	condDiffcomp_FDR_stats$fractionNull[rownames(condDiffcomp_FDR_stats) == "p_int"] <- 0.7
	condDiffcomp_FDR_stats$fractionNull[rownames(condDiffcomp_FDR_stats) == "c_slope"] <- 0.35
	condDiffcomp_FDR_stats$fractionNull[rownames(condDiffcomp_FDR_stats) == "n_slope"] <- 0.35
	condDiffcomp_FDR_stats$fractionNull[rownames(condDiffcomp_FDR_stats) == "p_slope"] <- 0.35
	
	condDiffcomp_FDR_stats$testStatCO[ncomp] <- seq(0.0001, 1, by = 0.0001)[which.min(unlist(lapply(seq(0.0001, 1, by = 0.0001), FDR.eval, pi = condDiffcomp_FDR_stats$fractionNull[ncomp], num.p = num.p, test.stats = pvals, FDRlev = FDR_desired, use.unif = TRUE)))]

	condDiffcomp_FDR_stats$FDRattained[ncomp] <- FDR.calc(pi = condDiffcomp_FDR_stats$fractionNull[ncomp], lambda = condDiffcomp_FDR_stats$testStatCO[ncomp], num.p = num.p, test.stats = pvals, use.unif = TRUE)
	
	
	discoveries <- data.frame(protein = rownames(p_compMat)[p_compMat[,ncomp] < condDiffcomp_FDR_stats$testStatCO[ncomp]], change = t_compMat[,ncomp][p_compMat[,ncomp] < condDiffcomp_FDR_stats$testStatCO[ncomp]], p_value = p_compMat[,ncomp][p_compMat[,ncomp] < condDiffcomp_FDR_stats$testStatCO[ncomp]], q_value = NA, stringsAsFactors = FALSE)
	
	discoveries$q_value <- sapply(c(1:length(discoveries[,1])), function(i){
		(condDiffcomp_FDR_stats$fractionNull[ncomp]*discoveries$p_value[i]*num.p) / (sum(pvals <= discoveries$p_value[i]))
		})
	
	condDiffcomp_FDR_stats$discoveries[ncomp] <- length(discoveries[,1])
	Discoveries[[ncomp]] <- discoveries
				
	}









##############

beta_val <- solve(t(design_mat) %*% diag(prot_W) %*% design_mat) %*% t(design_mat) %*% diag(prot_W) %*%  prot_val
	
	lm(prot_val ~ design_mat + 0, weights = prot_W)
	lm(prot_val_weighted ~ design_weighted + 0)
	
	
	residuals <- (design_mat %*% beta_val - prot_val)
	
	
	sum(residuals^2)/(n_samples - model_degree) * solve(t(design_mat) %*% diag(prot_W) %*% design_mat)
	
	infl_residuals <- residuals * sqrt(n_samples/(n_samples - model_degree))
	
	bs_resids <- matrix(sample(infl_residuals, n_samples*n_bs, replace = TRUE), ncol = n_samples)
	bs_response <- t(t(rep(1, times = n_bs))) %*% t(design_mat %*% beta_val) + bs_resids
		
	bs_coefs <- sapply(c(1:n_bs), function(bsN){
	solve(t(design_mat) %*% diag(prot_W) %*% design_mat) %*% t(design_mat) %*% diag(prot_W) %*% bs_response[bsN,]		
	})	
			
	p_valMat[pr,] <- pnorm(abs(beta_val), 0, apply(bs_coefs, 1, sd), lower.tail = FALSE)*2
	t_valMat[pr,] <- beta_val/diag(beta_se)
	



	
	












p_valMat <- matrix(NA, nrow = n_prot, ncol = length(levels(cond_info$limitation))*2)
rownames(p_valMat) <- colnames(specAbund)[1:n_prot]; 
colnames(p_valMat) <- unlist(sapply(c("int", "slope"), function(variable){
	paste(levels(cond_info$limitation), variable, sep = "_")
	}))
dir_valMat <- p_valMat

for(pr in 1:n_prot){
	
	valid.conds <- !is.na(specAbund[,pr])
	prot_val <- specAbund[,pr][valid.conds]
	prot_W <- specPrec[,pr][valid.conds]
	n_samples <- sum(valid.conds)
	
	
	for(cond in levels(cond_info$limitation)){
		######### test for a condition-specific abundance ##########
		
		#alternative model
		#alt_design <- cbind(ifelse(cond_info$limitation == cond, 1, 0), ifelse(cond_info$limitation != cond, 1, 0), cond_info$actualDR)[valid.conds,]
		alt_design <- cbind(cond_info$logSF, ifelse(cond_info$limitation == cond, 1, 0), ifelse(cond_info$limitation != cond, 1, 0), cond_info$actualDR)[valid.conds,]
		
		
		beta_val_alt <- solve(t(alt_design) %*% diag(prot_W) %*% alt_design) %*% t(alt_design) %*% diag(prot_W) %*%  prot_val
		#equivalent to, but faster than: lm(prot_val ~ alt_design + 0, weights = prot_W)
		
		alt_resids <- alt_design %*% beta_val_alt - prot_val
		altRSS <- sum(alt_resids^2)
		alt_degree <- length(beta_val_alt)
		
		#null model
		#null_design <- cbind(1, cond_info$actualDR)[valid.conds,]
		null_design <- cbind(cond_info$logSF, 1, cond_info$actualDR)[valid.conds,]
		
		beta_val_null <- solve(t(null_design) %*% diag(prot_W) %*% null_design) %*% t(null_design) %*% diag(prot_W) %*% prot_val
		null_fit <- null_design %*% beta_val_null
		nullRSS <- sum((null_fit - prot_val)^2)
		null_degree <- length(beta_val_null)
		
		test_f <- ((nullRSS - altRSS)/(alt_degree - null_degree))/(altRSS/(n_samples - alt_degree))
		
		alt_bs_resids <- alt_resids * sqrt((n_samples - null_degree)/(n_samples - alt_degree))		
		bs_resids <- matrix(sample(alt_bs_resids, n_samples*n_bs, replace = TRUE), ncol = n_samples)
		bs_response <- t(t(rep(1, times = n_bs))) %*% t(null_fit) + bs_resids
		
		#altRSS
		bs_altRSS <- sapply(c(1:n_bs), function(bsN){
			sum((alt_design %*% (solve(t(alt_design) %*% diag(prot_W) %*% alt_design) %*% t(alt_design) %*% diag(prot_W) %*% bs_response[bsN,]) - bs_response[bsN,])^2)
			})
		
		bs_nullRSS <- sapply(c(1:n_bs), function(bsN){
			sum((null_design %*% (solve(t(null_design) %*% diag(prot_W) %*% null_design) %*% t(null_design) %*% diag(prot_W) %*% bs_response[bsN,]) - bs_response[bsN,])^2)
			})
		
		F_bs <- ((bs_nullRSS - bs_altRSS)/(alt_degree - null_degree))/(altRSS/(n_samples - alt_degree))
		
		p_valMat[pr, paste(cond, "int", sep = "_")] <- sum(F_bs >= test_f)/n_bs
		#difference between the condition offset and all other offsets
		dir_valMat[pr, paste(cond, "int", sep = "_")] <- offset_diff <- beta_val_alt[1] - beta_val_alt[2]
		
		
		######### test for a condition-specific slope ##########
		
		#alternative model
		#alt_design <- cbind(sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}),  cond_info$actualDR*ifelse(cond_info$limitation == cond, 1, 0), cond_info$actualDR*ifelse(cond_info$limitation != cond, 1, 0))[valid.conds,]
		alt_design <- cbind(cond_info$logSF, sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}),  cond_info$actualDR*ifelse(cond_info$limitation == cond, 1, 0), cond_info$actualDR*ifelse(cond_info$limitation != cond, 1, 0))[valid.conds,]
		
		beta_val_alt <- solve(t(alt_design) %*% diag(prot_W) %*% alt_design) %*% t(alt_design) %*% diag(prot_W) %*% prot_val
		
		alt_resids <- alt_design %*% beta_val_alt - prot_val
		altRSS <- sum(alt_resids^2)
		alt_degree <- length(beta_val_alt)
		
		#null model
		#null_design <- cbind(sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}), cond_info$actualDR)[valid.conds,]
		null_design <- cbind(cond_info$logSF, sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}), cond_info$actualDR)[valid.conds,]
		
		beta_val_null <- solve(t(null_design) %*% diag(prot_W) %*% null_design) %*% t(null_design) %*% diag(prot_W) %*% prot_val
		null_fit <- null_design %*% beta_val_null
		nullRSS <- sum((null_fit - prot_val)^2)
		null_degree <- length(beta_val_null)
		
		
		test_f <- ((nullRSS - altRSS)/(alt_degree - null_degree))/(altRSS/(n_samples - alt_degree))
		
		alt_bs_resids <- alt_resids * sqrt((n_samples - null_degree)/(n_samples - alt_degree))		
		bs_resids <- matrix(sample(alt_bs_resids, n_samples*n_bs, replace = TRUE), ncol = n_samples)
		bs_response <- t(t(rep(1, times = n_bs))) %*% t(null_fit) + bs_resids
		
		#altRSS
		bs_altRSS <- sapply(c(1:n_bs), function(bsN){
			sum((alt_design %*% (solve(t(alt_design) %*% diag(prot_W) %*% alt_design) %*% t(alt_design) %*% diag(prot_W) %*% bs_response[bsN,]) - bs_response[bsN,])^2)
			})
		
		bs_nullRSS <- sapply(c(1:n_bs), function(bsN){
			sum((null_design %*% (solve(t(null_design) %*% diag(prot_W) %*% null_design) %*% t(null_design) %*% diag(prot_W) %*% bs_response[bsN,]) - bs_response[bsN,])^2)
			})
		
		F_bs <- ((bs_nullRSS - bs_altRSS)/(alt_degree - null_degree))/(altRSS/(n_samples - alt_degree))
		
		p_valMat[pr, paste(cond, "slope", sep = "_")] <- sum(F_bs >= test_f)/n_bs
		dir_valMat[pr, paste(cond, "slope", sep = "_")] <- offset_diff <- beta_val_alt[alt_degree - 1] - beta_val_alt[alt_degree]
		
		
		
		}
	}



FDR_desired <- 0.01
source('~/Desktop/Rabinowitz/Cancer_metabolomics/cancer_lib.R')


condDiff_FDR_stats <- data.frame(fractionNull = rep(NA, times = length(colnames(p_valMat))), testStatCO = rep(NA, times = length(colnames(p_valMat))), FDRattained = rep(NA, times = length(colnames(p_valMat))), discoveries = rep(NA, times = length(colnames(p_valMat))))
rownames(condDiff_FDR_stats) <- colnames(p_valMat)
Discoveries <- list()


for(ncomp in 1:length(condDiff_FDR_stats[,1])){
	
	pvals <- p_valMat[,ncomp]
	lambda.vals <- seq(0, 1, by = 0.001)
	num.p <- length(pvals)
	
	pi <- rep(NA, times = length(lambda.vals))
	for (l in 1:length(lambda.vals)){
		pi[l] <- sum(pvals > lambda.vals[l])/(num.p * (1 - lambda.vals[l]))
		}
	#plot(pi ~ lambda.vals, main = paste("estimate of fraction of TPs for", rownames(condDiff_FDR_stats)[ncomp]), ylim = c(0,1), pch = 16)
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "c_int"] <- 0.15
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "n_int"] <- 0.6
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "p_int"] <- 0.3
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "c_slope"] <- 0.85
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "n_slope"] <- 0.2
	condDiff_FDR_stats$fractionNull[rownames(condDiff_FDR_stats) == "p_slope"] <- 0.35
	
	condDiff_FDR_stats$testStatCO[ncomp] <- seq(0.0001, 1, by = 0.0001)[which.min(unlist(lapply(seq(0.0001, 1, by = 0.0001), FDR.eval, pi = condDiff_FDR_stats$fractionNull[ncomp], num.p = num.p, test.stats = pvals, FDRlev = FDR_desired, use.unif = TRUE)))]

	condDiff_FDR_stats$FDRattained[ncomp] <- FDR.calc(pi = condDiff_FDR_stats$fractionNull[ncomp], lambda = condDiff_FDR_stats$testStatCO[ncomp], num.p = num.p, test.stats = pvals, use.unif = TRUE)
	
	
	discoveries <- data.frame(protein = rownames(p_valMat)[p_valMat[,ncomp] < condDiff_FDR_stats$testStatCO[ncomp]], change = dir_valMat[,ncomp][p_valMat[,ncomp] < condDiff_FDR_stats$testStatCO[ncomp]], p_value = p_valMat[,ncomp][p_valMat[,ncomp] < condDiff_FDR_stats$testStatCO[ncomp]], q_value = NA, stringsAsFactors = FALSE)
	
	discoveries$q_value <- sapply(c(1:length(discoveries[,1])), function(i){
		(condDiff_FDR_stats$fractionNull[ncomp]*discoveries$p_value[i]*num.p) / (sum(pvals <= discoveries$p_value[i]))
		})
	
	condDiff_FDR_stats$discoveries[ncomp] <- length(discoveries[,1])
	Discoveries[[ncomp]] <- discoveries
				
	}
	
	


