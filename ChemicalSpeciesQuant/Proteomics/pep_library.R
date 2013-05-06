

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


orf2common <- function(x){
  
  require(org.Sc.sgd.db)
  
  d = org.Sc.sgdGENENAME
  mapped_genes = mappedkeys(d)
  dd = as.list(d[mapped_genes])
  
  sapply(x, function(orfs){
  orfIDs <- strsplit(orfs, split = '/')[[1]]
  matchedIDs <- orfIDs %in% names(dd)
  commonNames <- sapply(c(1:length(orfIDs)), function(NID){
      if(matchedIDs[NID]){
        unname(dd[orfIDs[NID]])[[1]]
        }else{
          orfIDs[NID]
          }
      })
  paste(commonNames, collapse = '/')    
  })
}
  
orf2entrez <- function(x){
  d <- org.Sc.sgdENTREZID
  mapped_genes <- mappedkeys(d)
  dd <- as.list(d[mapped_genes])
  
  matched_genes <- x %in% names(dd)
  
  entrez_names <- sapply(c(1:length(x)), function(el){
  if(matched_genes[el]){
    dd[[x[el]]]
    }else{
    NA 
      }
  })
  entrez_names
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

#prot_abund <- prot_abund_final; prots <- 1:6; num.cols <- 2; label_phospho = FALSE
#phospho_lib = phospho_sites_red; label_phospho = TRUE
#gene_samples[1:5], possibleMap, prot_abund, prot_prec, uniquePepMean, uniquePepPrecision, mixing_fract, conditions, label_phospho = TRUE, phospho_lib = phospho_sites_red, plotting_names = plotting_names, 
#num.cols = 1


plot_protein_lattice <- function(prots, possibleMap, prot_abund, prot_prec, uniquePepMean, uniquePepPrecision, mixing_fract, conditions, label_phospho = FALSE, phospho_lib = NULL, plotting_names = NULL, num.cols = 5){

#Plot a lattice of protein relative abundances as well as the peptides that match it
  
require(stringr)
require(gplots)
  
mix_frac_set <- c()
all_data_bind <- NULL
abline_df <- NULL
nutlab_df <- NULL
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
	
  if(label_phospho){
    #search for phospho-sites and label
  	psite_test_mat <- data.frame(peptide = colnames(uniquePepPrecision)[c(mixed_match,div_match)], prot_matches = colnames(prot_abund)[prot], nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)
    	for(pep in 1:length(psite_test_mat[,1])){
    	  psite_test_mat$nS[pep] <- length(grep("S", unlist(strsplit(psite_test_mat$peptide[pep], ""))))
    	  psite_test_mat$nR[pep] <- length(grep("R", unlist(strsplit(psite_test_mat$peptide[pep], ""))))
    	  psite_test_mat$nY[pep] <- length(grep("Y", unlist(strsplit(psite_test_mat$peptide[pep], ""))))
      }
  	
  	Pmods <- c("No P-site known", ifelse(knownPtest(psite_test_mat, phospho_lib), "Known phospho-site", "No P-site known"))  
  	PmodMat <- matrix(rep(Pmods, times = n_c), ncol = n_pp_plot, byrow = TRUE)
  }
  
  
	ncolors <- 100
	color_fac <- round(mixing_fract[mixed_match, prot], length(ncolors)+1)
	mix_frac_set <- union(mix_frac_set, color_fac)
	match_fact <- c(paste(unique(color_fac), "x", " match", sep = ""), "protein", "non-match")
	all_cols <- c("protein", paste(color_fac, "x", " match", sep = ""), rep("non-match", times = length(div_match)))
	prot_cols <- apply((rep(1, times = n_c) %*% t(c(1:n_pp_plot))), c(1,2), function(col){all_cols[col]})
	
  cond_offset <- data.frame(limitation = unique(conditions$limitation), offset = (0:(length(unique(conditions$limitation))-1))*2)
	
  xpos <- (1:n_c) %*% t(rep(1, n_pp_plot)) + rep(1, n_c) %*% t(((1:n_pp_plot)-1)*n_c) + sapply(rep(conditions$limitation, each = n_pp_plot), function(offset){cond_offset$offset[cond_offset$limitation == offset]})
  all_data <- data.frame(abund = c(t(prot_abundM)), abundMax = c(t(prot_max)), abundMin = c(t(prot_min)), xpos = c(xpos), color_fact = c(t(prot_cols)), stringsAsFactors = FALSE)
	if(label_phospho){
	  all_data <- cbind(all_data, phospho_sites = c(t(PmodMat)))
    }
	
  #switch protein name to be displayed to the name in plotting names if included
  if(!is.null(plotting_names)){
    prot_Display = unname(plotting_names[prot])
    }else{
      prot_Display = colnames(prot_abund)[prot]
      }
  
  all_data_bind <- rbind(all_data_bind, cbind(proteinNumber = prot_Display, all_data))
  abline_df <- rbind(abline_df, data.frame(xpos = xpos[cumsum(rle(rep(conditions$limitation, each = n_pp_plot))$lengths[-length(cond_offset[,1])])]+1.5, proteinNumber = prot_Display, stringsAsFactors = FALSE))
	
	#determine the position of condition labels
	cond_bounds <- c(0, xpos[cumsum(rle(rep(conditions$limitation, each = n_pp_plot))$lengths[-length(cond_offset[,1])])]+1, max(xpos))
	nutlab_df <- rbind(nutlab_df, data.frame(proteinNumber = prot_Display, xpos = sapply(1:length(cond_offset[,1]), function(cond){(cond_bounds[cond] + cond_bounds[cond + 1])/2}), ypos = range(c(all_data$abundMin, all_data$abundMax), na.rm = TRUE)[1] + diff(range(c(all_data$abundMin, all_data$abundMax), na.rm = TRUE))/40, label = toupper(cond_offset$limitation), stringsAsFactors = FALSE))
}

all_data_bind$proteinNumber <- as.factor(all_data_bind$proteinNumber)
abline_df$proteinNumber <- as.factor(abline_df$proteinNumber)
nutlab_df$proteinNumber <- as.factor(nutlab_df$proteinNumber)

 


match_cols <- c(colorpanel(ncolors + 1, low = "gainsboro", mid = "chocolate1", high = "firebrick1")[unique(sort(mix_frac_set))*ncolors + 1], "seagreen1", "gold")
#match_cols <- c(colorpanel(ncolors + 1, low = "gray70", high = "gray0")[unique(sort(mix_frac_set))*ncolors + 1], "red", "gold")
match_fact <- c(paste(sort(mix_frac_set), "x", " match", sep = ""), "protein", "non-match")

all_data_bind$color_fact <- as.factor(all_data_bind$color_fact)

plotter <- ggplot(all_data_bind)
if(label_phospho){plotter <- plotter + scale_linetype_manual(values = c("No P-site known" = 'solid', "Known phospho-site" = 'dotted')) + labs(linetype = "Phosphorylation Status")}
plotter <- plotter + xlab("Sample * Peptide") + ylab("Relative Abundance") + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())
plotter <- plotter + labs(colour = "Color label") + opts(axis.text.y = theme_text(colour = "red"), panel.background=theme_rect(fill="seashell"), strip.background = theme_rect(fill="tan1")) + scale_colour_manual(values = match_cols, limits = match_fact)
plotter <- plotter + facet_wrap( ~ proteinNumber, ncol = num.cols, scales = "free")
if(label_phospho){
  print(plotter + geom_linerange(aes(xpos, abund, ymin = abundMin, ymax = abundMax, colour = color_fact, linetype = factor(phospho_sites))) + geom_text(data = nutlab_df, aes(x = xpos, y = ypos, label = label, proteinNumber = proteinNumber), size = 3) + geom_vline(data = abline_df, aes(xintercept = xpos, proteinNumber = proteinNumber)))  
  }else{
    print(plotter + geom_linerange(aes(xpos, abund, ymin = abundMin, ymax = abundMax, colour = color_fact)) + geom_text(data = nutlab_df, aes(x = xpos, y = ypos, label = label, proteinNumber = proteinNumber), size = 3) + geom_vline(data = abline_df, aes(xintercept = xpos, proteinNumber = proteinNumber)))  
    }


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

knownPtest <- function(pepSummary, phospho_sites_red){
  
  #For a data.frame with peptides and the associated protein, determine which are known phospho-sites (return T)
  pSiteTest <- psite_table(pepSummary, phospho_sites_red)
  pSiteTest[1,] != 0
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



plot_proteinConnections <- function(overlapMat, subsumedIDs = NULL){
  
  require("igraph")
  
  graph_genes <- colnames(overlapMat)
  subsumedIDsN <- c(1:length(graph_genes))[graph_genes %in% subsumedIDs]
  
  graph_genes_common <- unname(orf2common(graph_genes))
  
  overlapMatWeights <- overlapMat/diag(overlapMat) 
  colnames(overlapMatWeights) <- rownames(overlapMatWeights) <- graph_genes_common
  
  prot_overlap_layout2 <- graph.adjacency(overlapMatWeights, weighted=T)
  
  
  # Set vertex attributes
  V(prot_overlap_layout2)$label <- V(prot_overlap_layout2)$name
  V(prot_overlap_layout2)$label.color <- rgb(0,0,.2,.6)
  V(prot_overlap_layout2)$size <- diag(overlapMat)/8
  V(prot_overlap_layout2)$label.cex <- .3
  V(prot_overlap_layout2)$frame.color <- NA
  V(prot_overlap_layout2)$color <- rgb(0,0.5,0.5,.8)
  
  
  V(prot_overlap_layout2)$label.color[subsumedIDsN] <- rgb(1,0,0,1)
  
  # Set edge attributes
  E(prot_overlap_layout2)$arrow.size <- 0
  
  # Set edge gamma according to edge weight
  egam <- (E(prot_overlap_layout2)$weight+.1)/max(E(prot_overlap_layout2)$weight+.1)
  E(prot_overlap_layout2)$color <- rgb(1,0,0,egam)
  
  prot_overlap_layout2 <- simplify(prot_overlap_layout2, remove.multiple=FALSE, remove.loops=TRUE)
  
  pdf("protOverlap.pdf")
  graph_layout = layout.fruchterman.reingold(prot_overlap_layout2,niter=500,area=vcount(prot_overlap_layout2)^2.3,repulserad=vcount(prot_overlap_layout2)^3.2)
  plot(prot_overlap_layout2, main = "Overlap of peptides", layout=graph_layout, vertex.label.dist=0.2, edge.width=2)
  dev.off()
}



##### simulations #######

varianceGP <- function(){
  
  require(msm)
  
  #### generative process for peptide variance ####
  I <- 1000 # number of peptides, each with their own dispersion
  C <- 100 # number of conditions
  max_reps <- 3 # number of possible replicates (fewer will be present to model missing values)
  #CV_fold <- 10
  
  ### knowns ###
  #known_sample_dispersion <- matrix(rlnorm(I*C*max_reps, 0, 0.5), ncol = C*max_reps, nrow = I)
  known_sample_dispersion <- matrix(rlnorm(I*C*max_reps, 0, 0), ncol = C*max_reps, nrow = I)
  
  
  ### unknowns ###
  peptide_overdispersion <- rlnorm(I, 0, 1)  # this is the parameter that needs to be estimated
  total_dispersion <- known_sample_dispersion * (peptide_overdispersion %*% t(rep(1, C*max_reps)))
  #pmissingval <- rbeta(I, 2, 1)
  pmissingval <- rbeta(I, 10, 0.001)
  sample_effects <- matrix(rnorm(I*C, 0, 1), ncol = C, nrow = I)
  
  ### input missing values ###
  for(i in 1:I){
    missingVals <- rbinom(C*max_reps, 1, pmissingval[i]) == 0
    total_dispersion[i,missingVals] <- NA
    }
  

  normality_tests <- data.frame(KS_A = rep(NA, I), KS_B = NA, KS_oracle = NA, KS_oracle_Nadjust = NA)
  peptide_ODdiff <- data.frame(A = rep(NA, I), B = NA, nquant_samples = NA)
  
  for(i in 1:I){
    
    ### generate technical replicates ###
    
    mock_data = data.frame(sample = rep(1:C, each = max_reps), effect = rep(sample_effects[i,], each = max_reps), var = total_dispersion[i,], known_dispersion = known_sample_dispersion[i,])
    mock_data <- mock_data[!is.na(mock_data$var),]
    mock_data$observed <- rnorm(length(mock_data[,1]), mock_data$effect, sqrt(mock_data$var))
    
    if(sum(table(mock_data$sample) >= 2) >= 5){
      
      ### fitted value is equal to the observations weighted by the inverse of their dispersion
      gls_model <- gls(observed ~ factor(sample), weights = ~known_dispersion, data = mock_data)
      #weighted.mean(mock_data$observed[mock_data$sample == 2], 1/mock_data$known_dispersion[mock_data$sample == 2])
      
      cond_est <-  summary(gls_model)$coef + ifelse(names(summary(gls_model)$coef) == "(Intercept)", 0, gls_model$coef[1]); names(cond_est) <- unique(mock_data$sample)
      mock_data$fitted <- sapply(mock_data$sample, function(x){cond_est[names(cond_est) == x]})
  
      repeated_conds <- mock_data[mock_data$sample %in% names(table(mock_data$sample))[table(mock_data$sample) >= 2],]
      repeated_conds$nreps <- sapply(repeated_conds$sample, function(x) sum(repeated_conds$sample == x))
      
      ### MLE of dispersion
      
      ODest1 <- mean(gls_model$resid^2 * repeated_conds$nreps/(repeated_conds$nreps - 1) * (1/mock_data$known_dispersion))
      ODest2 <- mean(gls_model$resid^2 * sapply(repeated_conds$sample, function(x){sample_size_correction(1/repeated_conds$known_dispersion[repeated_conds$sample == x])}) * (1/mock_data$known_dispersion))
      
      
      z1 <- (repeated_conds$observed - repeated_conds$fitted)*sqrt(repeated_conds$nreps/(repeated_conds$nreps - 1))/sqrt(repeated_conds$known_dispersion * ODest1)
      z2 <- (repeated_conds$observed - repeated_conds$fitted)*sqrt(sapply(repeated_conds$sample, function(x){sample_size_correction(1/repeated_conds$known_dispersion[repeated_conds$sample == x])}))/sqrt(repeated_conds$known_dispersion * ODest2)
      z3 <- (repeated_conds$observed - repeated_conds$fitted)*sqrt(repeated_conds$nreps/(repeated_conds$nreps - 1))/sqrt(repeated_conds$known_dispersion * peptide_overdispersion[i])
      ### adjust for weight in unbiased correction
      z4 <- (repeated_conds$observed - repeated_conds$fitted)*sapply(repeated_conds$sample, function(x){sqrt(sample_size_correction(1/repeated_conds$known_dispersion[repeated_conds$sample == x]))})/sqrt(repeated_conds$known_dispersion * peptide_overdispersion[i])
      
      normality_tests[i,] <- c(ks.test(z1[z1 >= 0], ptnorm, lower = 0)$p, ks.test(z2[z2 >= 0], ptnorm, lower = 0)$p, ks.test(z3[z3 >= 0], ptnorm, lower = 0)$p, ks.test(z4[z4 >= 0], ptnorm, lower = 0)$p)
      
      peptide_ODdiff[i,] <- c(c(ODest1, ODest2), length(repeated_conds[,1]))
    }
  }
  cor(cbind(peptide_ODdiff[,1:2], peptide_overdispersion), use = "complete.obs")
  
  plot(log2(peptide_ODdiff[,2])~ log2(peptide_overdispersion), col = green2red(C*max_reps + 1)[peptide_ODdiff$nquant_samples], pch = 16)
  
  
  plot((peptide_overdispersion - peptide_ODdiff[,1])/peptide_ODdiff[,1] ~ log2(peptide_ODdiff[,1]), pch = 16, cex = 0.5, col = green2red(C*max_reps + 1)[peptide_ODdiff$nquant_samples]) #fractional error in dispersion prediction
  plot((peptide_overdispersion - peptide_ODdiff[,1])/peptide_ODdiff[,1] ~ peptide_ODdiff$nquant_samples, pch = 16, cex = 0.5)
  
}





### normality test is uniform when known_sample_dispersion is constant 

https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
weighted.var2 <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
}


weighted.var <- function(x, w, na.rm = FALSE) {
    if (na.rm) {
        w <- w[i <- !is.na(x)]
        x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
na.rm)
}

sample_size_correction <- function(w) {
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    sum.w*(sum.w / (sum.w^2 - sum.w2)) 
}


weighted.var2(repeated_conds$observed[repeated_conds$sample == 1], repeated_conds$known_dispersion[repeated_conds$sample == 1])
weighted.var(repeated_conds$observed[repeated_conds$sample == 1], 1:3)



W = 1/repeated_conds$known_dispersion
sum(W)^2 - sum(W^2)





sum(1/W)





