
#library(norm)

# Functions

plotting_fxn <- function(){
	library(hexbin)
	library(gplots)
	library(colorRamps)	
	}

matrix_fxn <- function(){
	library(Matrix)
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
