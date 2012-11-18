
#library(norm)

# Functions

plotting_fxn <- function(){
	library(hexbin)
	library(gplots)
	library(ggplot2)
	library(colorRamps)	
	}

matrix_fxn <- function(){
	library(Matrix)
	library(flexclust)
	}

write.output <- function(tab, output){
  #write an output table in a way where all columns will have a column name
  #rownames go to 1st column, 1st column name is ``Gene''
  
  tab <- data.frame(Gene = rownames(tab), tab, stringsAsFactors = FALSE)
  write.table(tab, file = output, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


var_calc <- function(sampleIC, STDvar_fit){
	2^(STDvar_fit[1] + STDvar_fit[2]*log2(sampleIC))
	}

normFactor <- function(alpha, logLight, logHeavy){
	#determine a factor that will be added to each logLight abundance such that the global sum of squares is minimized between logLight and logHeavy
	sum(((logLight + alpha) - logHeavy)^2)
	}

max_non_zero <- function(vec){
	max(vec[vec != 0])
	}


max_non_NA <- function(vec){
	max(vec[!is.na(vec)])
	}

orf2desc <- function(x) {
	
	require(org.Sc.sgd.db)
	
	 d = org.Sc.sgdDESCRIPTION
	 mapped_probes <- mappedkeys(d)
	 dd = as.list(d[mapped_probes])
	 out = sapply(x, function(x) {dd[as.character(x)][[1]]})
	 name_frame <- data.frame(geneName = names(out), geneDescrip = out, stringsAsFactors = FALSE)
	 name_out <- apply(name_frame, 1, function(gene){
	 	paste(gene, collapse = ": ")
	 	})
	 name_out <- unname(name_out)
	 name_out
	}


residPrec <- function(generate_plots = FALSE){

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
	
STDvar_fit	
}

plot_protein_lattice <- function(prots, possibleMap, prot_abund, prot_prec, uniquePepMean, uniquePepPrecision, mixing_fract, conditions, num.cols = 5){

#Plot a lattice of protein relative abundances as well as the peptides that match it

mix_frac_set <- c()
all_data_bind <- NULL
abline_df <- NULL
for(prot in prots){
	row_match <- c(1:n_p)[possibleMap[,prot]]
	div_match <- row_match[diag(mixing_fract[row_match, n_prot + row_match]) == 1]
	mixed_match <- row_match[!(row_match %in% div_match)]
	
	prot_abundM <- cbind(prot_abund[,prot], uniquePepMean[,c(mixed_match,div_match)])
	prot_precM <- cbind(prot_prec[,prot], as.matrix(uniquePepPrecision[,c(mixed_match, div_match)]))
	prot_abundM[prot_abundM == 0] <- NA; prot_precM[is.na(prot_abundM)] <- NA
	prot_max <- prot_abundM + 2*sqrt(1/prot_precM)
	prot_min <- prot_abundM - 2*sqrt(1/prot_precM)
	
	n_pp_plot <- length(prot_abundM[1,])
	n_c <- length(prot_abundM[,1])
	
	ncolors <- 100
	color_fac <- round(mixing_fract[mixed_match, prot], length(ncolors)+1)
	mix_frac_set <- union(mix_frac_set, color_fac)
	match_fact <- c(paste(unique(color_fac), "x", " match", sep = ""), "protein", "non-match")
	all_cols <- c("protein", paste(color_fac, "x", " match", sep = ""), rep("non-match", times = length(div_match)))
	prot_cols <- apply((rep(1, times = n_c) %*% t(c(1:n_pp_plot))), c(1,2), function(col){all_cols[col]})
	
  cond_offset <- data.frame(limitation = unique(conditions$limitation), offset = (0:(length(unique(conditions$limitation))-1))*2)
	
  xpos <- (1:n_c) %*% t(rep(1, n_pp_plot)) + rep(1, n_c) %*% t(((1:n_pp_plot)-1)*n_c) + sapply(rep(conditions$limitation, each = n_pp_plot), function(offset){cond_offset$offset[cond_offset$limitation == offset]})
  all_data <- data.frame(abund = c(t(prot_abundM)), abundMax = c(t(prot_max)), abundMin = c(t(prot_min)), xpos = c(xpos), color_fact = c(t(prot_cols)), stringsAsFactors = FALSE)
	all_data_bind <- rbind(all_data_bind, cbind(proteinNumber = colnames(prot_abund)[prot], all_data))
  abline_df <- rbind(abline_df, data.frame(xpos = xpos[cumsum(rle(rep(conditions$limitation, each = n_pp_plot))$lengths[-length(cond_offset[,1])])]+1, proteinNumber = colnames(prot_abund)[prot], stringsAsFactors = FALSE))
	}
all_data_bind$proteinNumber <- as.factor(all_data_bind$proteinNumber)
abline_df$proteinNumber <- as.factor(abline_df$proteinNumber)

match_cols <- c(rainbow(ncolors + 1, start = 0.7, end = 0)[unique(sort(mix_frac_set))*ncolors + 1], "black", "darkgray")
match_fact <- c(paste(sort(mix_frac_set), "x", " match", sep = ""), "protein", "non-match")

all_data_bind$color_fact <- as.factor(all_data_bind$color_fact)


plotter <- ggplot(all_data_bind, aes(xpos, abund, ymin = abundMin, ymax = abundMax, colour = color_fact))
plotter <- plotter + xlab("Sample * Peptide") + ylab("Relative Abundance") + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())
plotter <- plotter + labs(colour = "Color label")
plotter <- plotter + facet_wrap( ~ proteinNumber, ncol = num.cols, scales = "free")
print(plotter + geom_linerange() +  geom_vline(data = abline_df, aes(xintercept = xpos, proteinNumber = proteinNumber))
+ scale_colour_manual(values = match_cols, limits = match_fact)
+ opts(axis.text.y = theme_text(colour = "red")))

}







plot_protein_add <- function(prot, possibleMap, prot_abund, uniquePepMean, uniquePepPrecision){
row_match <- c(1:n_p)[possibleMap[,prot]]
#row_match <- c(1:n_p)[possibleMap[,prot] != 0]
div_match <- row_match[diag(mixing_fract[row_match, n_prot + row_match]) == 1]
mixed_match <- row_match[!(row_match %in% div_match)]

prot_abundM <- cbind(prot_abund[,prot], uniquePepMean[,c(mixed_match,div_match)])
prot_precM <- cbind(prot_prec[,prot], as.matrix(uniquePepPrecision[,c(mixed_match, div_match)]))
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
plotter <- plotter + xlab("sample * peptide") + ylab("relative abundance") 
plotter <- plotter + labs(colour = paste(colnames(possibleMap)[prot], ": protein ", prot, sep = ""))
print(plotter + geom_linerange()
+ scale_colour_manual(values = match_cols, limits = match_fact)
+ opts(axis.text.y = theme_text(colour = "red")))
}



psite_table <- function(pepSummary, phospho_sites_red){

	phospho_match <- sapply(c(1:length(pepSummary[,1])), function(pep){
		#for(pep in c(1:length(pepSummary[,1]))){
		
		known_prot_mods <- phospho_sites_red[phospho_sites_red$A %in% unlist(strsplit(pepSummary$prot_matches[pep], "/")),]
		if(length(known_prot_mods[,1]) != 0){
			match_ev <- sapply(1:length(known_prot_mods[,1]), function(p_site){
				#find the position of the relevent peptide in the matched protein: compare this interval ot the phosphosites
				pep_pos <- str_locate(known_prot_mods[p_site,3], tolower(pepSummary$peptide[pep]))
				if(is.na(pep_pos[1,1])){
					NA
					}else{
						pep_interval <- str_locate(known_prot_mods[p_site,3], tolower(pepSummary$peptide[pep]))
						if(sub('[A-Za-z]+', '', known_prot_mods[p_site,2]) >= pep_interval[1] & sub('[A-Za-z]+', '', known_prot_mods[p_site,2]) <= pep_interval[2]){
							TRUE
							}else{
								FALSE
								}
						}
			})
			data.frame(nTRUE = sum(match_ev[!is.na(match_ev)]), nFALSE = sum(match_ev[!is.na(match_ev)] == FALSE), nNA = sum(is.na(match_ev)))
			}else{
				data.frame(nTRUE = 0, nFALSE = 0, nNA = 0)
				}
		})
	
	}


knownPsites <- function(pepSummary, phospho_sites_red){

	#For a data.frame with the protein matches to a given peptide, go through peptides 1 by 1 and determine all of the phospho-sites for proteins matching that peptide.  If the peptide is known to be phosphorylated return TRUE and then count up the number of matches
	
	phospho_match <- psite_table(pepSummary, phospho_sites_red)
	
	phospho_match <- matrix(unlist(t(phospho_match)), ncol = 3)	
	table(phospho_match[,1] != 0)
	}


bs_tstat <- function(specAbund, specPrec, n_bs = 100, FDR_desired = 0.05){

output_list <- list()

p_valMat <- matrix(NA, nrow = length(specAbund[1,]), ncol = length(levels(cond_info$limitation))*2 + 1)
rownames(p_valMat) <- colnames(specAbund)
colnames(p_valMat) <- c("prot_abund", unlist(sapply(c("int", "slope"), function(variable){
	paste(levels(cond_info$limitation), variable, sep = "_")
	})))
t_valMat <- p_valMat


for(pr in 1:length(specAbund[1,])){
	
	valid.conds <- !is.na(specAbund[,pr])
	prot_val <- specAbund[,pr][valid.conds]
	prot_W <- specPrec[,pr][valid.conds]
	n_samples <- sum(valid.conds)
	
	design_mat <- cbind(cond_info$logSF, sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}), sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)})*cond_info$actualDR)[valid.conds,]
	
	prot_val_weighted <- prot_val * sqrt(prot_W)
	design_weighted <- design_mat * t(t(sqrt(prot_W))) %*% rep(1, times = length(p_valMat[1,]))
	
	#see if any coefficients can't be estimated and change design matrix accordingly
	evaluatable_coefs <- !is.na(lm(prot_val_weighted ~ design_weighted + 0)$coef)
	design_weighted <- design_weighted[,evaluatable_coefs]
	
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
	p_valMat[pr,evaluatable_coefs] <- 1 - (abs(apply((beta_val %*% rep(1, times = n_bs) + bs_coefs*sqrt(diag(beta_se))) >= 0, 1, sum) - n_bs/2)/(n_bs/2)) 
	t_valMat[pr,evaluatable_coefs] <- beta_val/sqrt(diag(beta_se))
	
	if(floor(pr/100) == pr/100){
		print(paste(pr, "proteins down"))
		}
	
	
	}

condDiff_FDR_stats <- data.frame(fractionNull = rep(NA, times = length(colnames(p_valMat))), testStatCO = rep(NA, times = length(colnames(p_valMat))), discoveries = rep(NA, times = length(colnames(p_valMat))))
rownames(condDiff_FDR_stats) <- colnames(p_valMat)
Discoveries <- list()

for(ncomp in 1:length(condDiff_FDR_stats[,1])){
	
	pvals <- p_valMat[,ncomp][!is.na(p_valMat[,ncomp])]
	qval_obj <- qvalue(pvals, lambda=seq(0,0.90,0.01), pi0.method="bootstrap", fdr.level=FDR_desired)
	condDiff_FDR_stats$fractionNull[ncomp] <- qval_obj$pi0
	condDiff_FDR_stats$testStatCO[ncomp] <- max(qval_obj$pvalues[qval_obj$qvalues <= FDR_desired])
	condDiff_FDR_stats$testStatCO[ncomp] <- max(condDiff_FDR_stats$testStatCO[ncomp], 0)
	condDiff_FDR_stats$discoveries[ncomp] <- sum(qval_obj$qvalues <= FDR_desired)
	
	discoveries <- data.frame(protein = rownames(p_valMat)[!is.na(p_valMat[,ncomp])][qval_obj$qvalues <= FDR_desired], change = sapply(rownames(p_valMat)[!is.na(p_valMat[,ncomp])][qval_obj$qvalues <= FDR_desired], function(name_match){t_valMat[,ncomp][names(t_valMat[,ncomp]) == name_match]}), p_value = qval_obj$pvalues[qval_obj$qvalues <= FDR_desired], q_value = qval_obj$qvalues[qval_obj$qvalues <= FDR_desired], stringsAsFactors = FALSE)
	rownames(discoveries) <- NULL
	
	Discoveries[[ncomp]] <- discoveries
	}

output_list$FDR_stats <- condDiff_FDR_stats
output_list$Discoveries <- Discoveries
output_list
}

	




bs_Fstat <- function(specAbund, specPrec, n_bs = 10, FDR_desired = 0.05){

output_list <- list()

p_valMat <- matrix(NA, nrow = length(specAbund[1,]), ncol = length(levels(cond_info$limitation))*2)
rownames(p_valMat) <- colnames(specAbund)
colnames(p_valMat) <- unlist(sapply(c("int", "slope"), function(variable){
	paste(levels(cond_info$limitation), variable, sep = "_")
	}))
dir_valMat <- p_valMat

for(pr in 1:length(specAbund[1,])){
	
	valid.conds <- !is.na(specAbund[,pr])
	prot_val <- specAbund[,pr][valid.conds]
	prot_W <- specPrec[,pr][valid.conds]
	n_samples <- sum(valid.conds)
	
	
	for(cond in levels(cond_info$limitation)){
		######### test for a condition-specific abundance ##########
		
		#alternative model
		
		alt_design <- cbind(cond_info$logSF, ifelse(cond_info$limitation == cond, 1, 0), ifelse(cond_info$limitation != cond, 1, 0), cond_info$actualDR)[valid.conds,]
		evaluatable_coef_alt <- !is.na(lm(prot_val ~ alt_design + 0, weights = prot_W)$coef)
		alt_design <- alt_design[,evaluatable_coef_alt]
	
		beta_val_alt <- solve(t(alt_design) %*% diag(prot_W) %*% alt_design) %*% t(alt_design) %*% diag(prot_W) %*%  prot_val
		#equivalent to, but faster than: lm(prot_val ~ alt_design + 0, weights = prot_W)
		
		alt_resids <- alt_design %*% beta_val_alt - prot_val
		altRSS <- sum(alt_resids^2)
		alt_degree <- length(beta_val_alt)
		
		#null model
		
		null_design <- cbind(cond_info$logSF, 1, cond_info$actualDR)[valid.conds,]
		evaluatable_coef_null <- !is.na(lm(prot_val ~ null_design + 0, weights = prot_W)$coef)
		null_design <- null_design[,evaluatable_coef_null]
		
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
		
		alt_design <- cbind(cond_info$logSF, sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}),  cond_info$actualDR*ifelse(cond_info$limitation == cond, 1, 0), cond_info$actualDR*ifelse(cond_info$limitation != cond, 1, 0))[valid.conds,]
		evaluatable_coef_alt <- !is.na(lm(prot_val ~ alt_design + 0, weights = prot_W)$coef)
		alt_design <- alt_design[,evaluatable_coef_alt]
		
		beta_val_alt <- solve(t(alt_design) %*% diag(prot_W) %*% alt_design) %*% t(alt_design) %*% diag(prot_W) %*% prot_val
		
		alt_resids <- alt_design %*% beta_val_alt - prot_val
		altRSS <- sum(alt_resids^2)
		alt_degree <- length(beta_val_alt)
		
		#null model
		null_design <- cbind(cond_info$logSF, sapply(levels(cond_info$limitation), function(cond){ifelse(cond_info$limitation == cond, 1, 0)}), cond_info$actualDR)[valid.conds,]
		evaluatable_coef_null <- !is.na(lm(prot_val ~ null_design + 0, weights = prot_W)$coef)
		null_design <- null_design[,evaluatable_coef_null]
		
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
	
	if(floor(pr/10) == pr/10){
		print(paste(pr, "proteins down"))
		}
	
	}

condDiff_FDR_stats <- data.frame(fractionNull = rep(NA, times = length(colnames(p_valMat))), testStatCO = rep(NA, times = length(colnames(p_valMat))), discoveries = rep(NA, times = length(colnames(p_valMat))))
rownames(condDiff_FDR_stats) <- colnames(p_valMat)
Discoveries <- list()

for(ncomp in 1:length(condDiff_FDR_stats[,1])){
	
	pvals <- p_valMat[,ncomp][!is.na(p_valMat[,ncomp])]
	qval_obj <- qvalue(pvals, lambda=seq(0,0.90,0.01), pi0.method="bootstrap", fdr.level=FDR_desired)
	condDiff_FDR_stats$fractionNull[ncomp] <- qval_obj$pi0
	condDiff_FDR_stats$testStatCO[ncomp] <- max(qval_obj$pvalues[qval_obj$qvalues <= FDR_desired])
	condDiff_FDR_stats$testStatCO[ncomp] <- max(condDiff_FDR_stats$testStatCO[ncomp], 0)
	condDiff_FDR_stats$discoveries[ncomp] <- sum(qval_obj$qvalues <= FDR_desired)
	
	discoveries <- data.frame(protein = rownames(p_valMat)[!is.na(p_valMat[,ncomp])][qval_obj$qvalues <= FDR_desired], change = sapply(rownames(p_valMat)[!is.na(p_valMat[,ncomp])][qval_obj$qvalues <= FDR_desired], function(name_match){dir_valMat[,ncomp][names(dir_valMat[,ncomp]) == name_match]}), p_value = qval_obj$pvalues[qval_obj$qvalues <= FDR_desired], q_value = qval_obj$qvalues[qval_obj$qvalues <= FDR_desired], stringsAsFactors = FALSE)
	rownames(discoveries) <- NULL
	
	Discoveries[[ncomp]] <- discoveries
	}

output_list$FDR_stats <- condDiff_FDR_stats
output_list$Discoveries <- Discoveries
output_list
}





