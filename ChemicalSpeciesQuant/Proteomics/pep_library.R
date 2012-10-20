
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

logL <- function(prot_abund, uniquePepMean, mixing_fract, uniquePepPrecision){
	
	#calculate the complete data log-likelihood
	p_est <- prot_abund %*% t(mixing_fract)
	pos_lik <- -0.5*sum(uniquePepPrecision*(uniquePepMean - p_est)^2)
	prior_adj <- sum(mixing_fract*prior_mat_likadj)
	pos_lik + prior_adj
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

plot_protein_lattice <- function(prots, possibleMap, prot_abund, prot_prec, uniquePepMean, uniquePepPrecision, mixing_fract){

prots <- c(26:50)
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
	all_data_bind <- rbind(all_data_bind, cbind(proteinNumber = prot, all_data))
}
all_data_bind$proteinNumber <- as.factor(all_data_bind$proteinNumber)

match_cols <- c(rainbow(ncolors + 1, start = 0, end = 0.8)[unique(sort(mix_frac_set))*ncolors + 1], "black", "darkgray")
match_fact <- c(paste(sort(mix_frac_set), "x", " match", sep = ""), "protein", "non-match")

all_data_bind$color_fact <- as.factor(all_data_bind$color_fact)


plotter <- ggplot(all_data_bind, aes(xpos, abund, ymin = abundMin, ymax = abundMax, colour = color_fact))
plotter <- plotter + xlab("sample * peptide") + ylab("relative abundance") 
plotter <- plotter + labs(colour = "Color label")
plotter <- plotter + facet_wrap( ~ proteinNumber, ncol = 5, scales = "free")
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


#from Cookbook for R
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, by.row=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

