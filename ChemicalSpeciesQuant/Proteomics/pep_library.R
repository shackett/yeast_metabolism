
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

plot_protein_lattice <- function(prots, possibleMap, prot_abund, prot_prec, uniquePepMean, uniquePepPrecision, mixing_fract, num.cols = 5){



mix_frac_set <- c()
all_data_bind <- NULL
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
	
	xpos <- (1:n_c) %*% t(rep(1, n_pp_plot)) + rep(1, n_c) %*% t(((1:n_pp_plot)-1)*n_c + ((1:n_pp_plot)-1)*2)
	all_data <- data.frame(abund = c(t(prot_abundM)), abundMax = c(t(prot_max)), abundMin = c(t(prot_min)), xpos = c(xpos), color_fact = c(t(prot_cols)), stringsAsFactors = FALSE)
	all_data_bind <- rbind(all_data_bind, cbind(proteinNumber = colnames(prot_abund)[prot], all_data))
}
all_data_bind$proteinNumber <- as.factor(all_data_bind$proteinNumber)

match_cols <- c(rainbow(ncolors + 1, start = 0, end = 0.8)[unique(sort(mix_frac_set))*ncolors + 1], "black", "darkgray")
match_fact <- c(paste(sort(mix_frac_set), "x", " match", sep = ""), "protein", "non-match")

all_data_bind$color_fact <- as.factor(all_data_bind$color_fact)


plotter <- ggplot(all_data_bind, aes(xpos, abund, ymin = abundMin, ymax = abundMax, colour = color_fact))
plotter <- plotter + xlab("sample * peptide") + ylab("relative abundance") 
plotter <- plotter + labs(colour = "Color label")
plotter <- plotter + facet_wrap( ~ proteinNumber, ncol = num.cols, scales = "free")
print(plotter + geom_linerange()
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


