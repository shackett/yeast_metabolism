
########### Function used in flux prediction and rxn species conversions ###########

write_stoiMat = function(metabolites, reactions, corrFile, rxnFile, internal_names = FALSE){

#create the stoichiometry matrix, indicating the change of metabolites (rows) per chemical reaction (columns)

	stoiMat <- matrix(data = 0, ncol = length(reactions), nrow = length(metabolites))

	metName = rep(NA, times = length(metabolites))
	for (i in 1:length(metabolites)){
		metName[i] = as.character(corrFile$SpeciesName[corrFile$SpeciesID %in% metabolites[i]][1])
		}
	rxnName = rep(NA, times = length(reactions))
	for (i in 1:length(reactions)){
		rxnName[i] = as.character(rxnFile$Reaction[rxnFile$ReactionID %in% reactions[i]][1])
		}
	rownames(stoiMat) <- metName	
	colnames(stoiMat) <- rxnName

	for (i in 1:length(rxnFile[,1])){
		stoiMat[c(1:length(metabolites))[metabolites == rxnStoi$Metabolite[i]], c(1:length(reactions))[reactions == rxnStoi$ReactionID[i]]] <- rxnStoi$StoiCoef[i]		}

	if (internal_names == FALSE){
		save(stoiMat, file = "yeast_stoi.R")} else {
		rownames(stoiMat) <- metabolites	
		colnames(stoiMat) <- reactions
		save(stoiMat, file = "yeast_stoi.R")
			}
	}

perfect.match <- function(source, query, corrFile){
all_char <- "[[:graph:][:space:]]"
	
tmp <- corrFile[grep(source, query)[!(grep(source, query) %in% union(grep(paste(all_char, source, sep = ""), query), grep(paste(source, all_char, sep = ""), query)))],]
if(length(tmp[,1]) == 0){tmp <- corrFile[grep(source, query, fixed = TRUE),]}
tmp
	}


rxn_search = function(stoiMat, search_string, is_rxn = TRUE, index = FALSE){
	#search by metabolite or reactant and return all reactions and nonzero metabolites.
  #stoiMat rows and columns must be named with metabolite and enzyme common names: named_stoi
	if (is_rxn == TRUE){
		colz = grep(search_string, colnames(stoiMat), fixed = TRUE)
		} else {
		met = grep(search_string, rownames(stoiMat), fixed = TRUE)
		if (length(met) == 1){
			colz = c(1:length(stoiMat[1,]))[stoiMat[met,] != 0]
			} else {
			colz = c(1:length(stoiMat[1,]))[apply(stoiMat[met,], 2, is.not.zero)]
		}}
	
	if(length(colz) == 0){
		print("no hits")
		} else {
			if(index == TRUE){
				colz
				} else {
			
			rxns = stoiMat[,colz]
			if(is.vector(rxns)){
				c(colz, rxns[rxns != 0])
				} else {
					output <- rbind(colz, rxns[apply(rxns, 1, is.not.zero),])
					colnames(output) = colnames(stoiMat)[colz]
					output
					}}
		}
	}

flip.rxn <- function(reactions, joint.stoi){
	stoi <- joint.stoi
	stoi[,colnames(joint.stoi) %in% reactions] <- stoi[,colnames(joint.stoi) %in% reactions]*-1
	stoi
	}
	
is.not.zero = function(vec){
	length(vec[vec!=0]) != 0
	}	
	
######### fxns to convert between IDs and species ######

metIDtoSpec <- function(meta){
	sapply(meta, function(x){
		corrFile$SpeciesName[corrFile$SpeciesID == x]
		})}

rxnIDtoEnz <- function(rxn){
	sapply(rxn, function(x){
		rxnFile$Reaction[rxnFile$ReactionID == x][1]
		})}
		
rxnIDtoGene <- function(rxns){
  sapply(rxns, function(rx){
    paste(unique(strsplit(paste(rxnFile[rxnFile$ReactionID == rx,]$MetName[is.na(rxnFile[rxnFile$ReactionID == rx,]$StoiCoef)], collapse = ':'), ':')[[1]]), collapse = '/')
    })}

metToCHEBI <- function(mets){
	#associate species IDs and CHEBI ids where available/applicable
	if(length(grep("chebi", rxnparFile[,3][rxnparFile[,1] == corrFile$SpeciesType[corrFile$SpeciesID %in% mets]])) == 0){
		NA
		}else{
	unlist(strsplit(rxnparFile[,3][rxnparFile[,1] == corrFile$SpeciesType[corrFile$SpeciesID %in% mets]], split = "%3A"))[2]	
	}}

rxnIDtoSGD <- function(rxnIDs){
  #output the compartment where a reaction occurs followed by all of the genes involved in the rxn  
  
  output <- t(sapply(rxnIDs, function(rxnID){
    tmp <- rxnFile[rxnFile$ReactionID == rxnID,]
    c(tmp$Compartment[1], paste(unique(strsplit(paste(tmp$MetName[is.na(tmp$StoiCoef)], collapse = ':'), ':')[[1]]), collapse = ':'))
  }))
output
}


eval_mets <- function(query_met, grep_it = FALSE){
	#find all of the reactions that a greped or exact matched metabolite participates in and then get all of the other metabolites also in those reactions
	if(grep_it == TRUE){
		met_matches <- grep(query_met, metIDtoSpec(rownames(stoiMat)))
		}else{
			met_matches <- c(1:length(stoiMat[,1]))[metIDtoSpec(rownames(stoiMat)) %in% query_met]
			}
	if(length(met_matches) == 0){print("miss")}
	if(length(met_matches) == 1){
		eval_mat <- stoiMat[apply(stoiMat[,stoiMat[met_matches,] != 0] != 0, 1, sum) != 0,stoiMat[met_matches,] != 0]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))
		}else{
		eval_mat <- stoiMat[apply(stoiMat[,apply(stoiMat[met_matches,] != 0, 2, sum) != 0] != 0, 1, sum) != 0, apply(stoiMat[met_matches,] != 0, 2, sum) != 0]; rownames(eval_mat) <- metIDtoSpec(rownames(eval_mat)); colnames(eval_mat) <- rxnIDtoEnz(colnames(eval_mat))
		}
	eval_mat
	}		
  
  
  
reaction_info <- function(rxnName){
    
  ##### write a function to list:
  # reactants -> products
  # Reaction name and designation
  # Thermodynamics
  # KEGG and EC reactions name
    
  rxnStoi <- stoiMat[,colnames(stoiMat) == rxnName][stoiMat[,colnames(stoiMat) == rxnName] != 0]
  speciesNames <- metIDtoSpec(names(rxnStoi))
  rxnDir <- reversibleRx$reversible[reversibleRx$rx == rxnName]
  if(rxnDir == 1){rxnDir <- " -> "}
  if(rxnDir == 0){rxnDir <- " <=> "}
  if(rxnDir == -1){rxnDir <- " <- "}
    
  substrate_prep <- paste(sapply(c(1:length(rxnStoi[rxnStoi < 0])), function(x){
    tmp <- (rxnStoi[rxnStoi < 0] * -1)[x]
    if(tmp == 1){tmp <- ''}
    paste(tmp, speciesNames[rxnStoi < 0][x])
  }), collapse = ' + ')
    
  product_prep <- paste(sapply(c(1:length(rxnStoi[rxnStoi > 0])), function(x){
    tmp <- (rxnStoi[rxnStoi > 0])[x]
    if(tmp == 1){tmp <- ''}
    paste(tmp, speciesNames[rxnStoi > 0][x])
  }), collapse = ' + ')
    
  rxList <- list()
  rxList$reaction = unname(rxnIDtoEnz(rxnName))
  rxList$enzymes = unname(rxnIDtoGene(rxnName))
  rxList$stoichiometry = paste(substrate_prep, rxnDir, product_prep)
  rxList$thermo = reversibleRx[reversibleRx$rx == rxnName,]
  rxList  
  }  



########## Functions used in optimization of fitted flux versus actual flux ############

species_plot <- function(run_rxn, flux_fit, chemostatInfo){
  #generate a list of four plots:
  #1: Flux predicted from FBA versus parametric form
  #2: Flux predicted from FBA - actual
  #3: Relationship between metabolites/enzyme levels and condition
  #4: Relationship between metabolite/enzyme levels and flux carried
  
  output_plots <- list()
  
  scatter_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "azure"), legend.position = "right", 
      panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "pink"), axis.ticks = element_line(colour = "pink"), strip.background = element_rect(fill = "cyan"),
      legend.key.size = unit(3, "line"), legend.text = element_text(size = 40, face = "bold"))

  flux_plot <- data.frame(actual = run_rxn$flux, predicted = flux_fit$fitted_flux, condition = chemostatInfo$limitation, DR = chemostatInfo$actualDR)
  output_plots$flux_plot <- ggplot() + geom_path(data = flux_plot, aes(x = actual, y = predicted, col = condition, size = 2)) + scatter_theme +
    geom_point(data = flux_plot, aes(x = actual, y = predicted, col = condition, size = DR*50)) + scale_size_identity() + 
    scale_color_brewer("Limitation", palette = "Set2") + ggtitle("Flux predicted from FBA versus parametric form") + geom_abline(intercept = 0, slope = 1, size = 2)
  
  output_plots$FBA_flux <- ggplot() + geom_path(data = flux_plot, aes(x = DR, y = actual, col = condition, size = 2)) + scatter_theme +
    geom_point(data = flux_plot, aes(x = DR, y = actual, col = condition, size = DR*50)) + scale_size_identity() + 
    scale_color_brewer("Limitation", palette = "Set2") + ggtitle("Flux predicted from FBA - actual")
    
  all_species <- data.frame(run_rxn$enzymes, run_rxn$metabolites)
  colnames(all_species) <- run_rxn$all_species$commonName
  
  all_species_tab <- run_rxn$all_species[colSums(all_species != 1) != 0,]
  all_species <- all_species[,colSums(all_species != 1) != 0]
  
  species_df <- melt(data.frame(all_species, condition = chemostatInfo$limitation, DR = chemostatInfo$actualDR, flux = run_rxn$flux), id.vars = c("condition", "DR", "flux"))
  
  output_plots$species <- ggplot() + geom_path(data = species_df, aes(x = DR, y = value, col = condition, size = 2)) + facet_wrap(~ variable, scale = "free_y") +
    scatter_theme + scale_size_identity() + scale_color_brewer("Limitation", palette = "Set2") + scale_y_continuous("Relative concentration") +
    ggtitle("Relationship between metabolites/enzyme levels and condition")
  
  output_plots$flux_species <- ggplot() + geom_path(data = species_df, aes(x = value, y = flux, col = condition, size = 2)) + facet_wrap(~ variable, scale = "free") +
    scatter_theme + scale_size_identity() + scale_color_brewer("Limitation", palette = "Set2") + scale_x_continuous("Relative concentration") +
    geom_point(data = species_df, aes(x = value, y = flux, col = condition, size =  DR*30)) + ggtitle("Relationship between metabolite/enzyme levels and flux carried")
  
  output_plots
  }






flux_fitting <- function(run_rxn, par_markov_chain, par_likelihood){
  # predict flux based upon parameter sets to determine how much variance in flux can be accounted for using the prediction
  param_interval <- exp(apply(par_markov_chain, 2, function(x){quantile(x, probs = c(0.025, 0.975))}))
  param_interval <- data.frame(cbind(t(param_interval), median = exp(apply(par_markov_chain, 2, median)), MLE = exp(par_markov_chain[which.max(par_likelihood$likelihood),])))
  
  par_stack <- rep(1, n_c) %*% t(exp(par_markov_chain[which.max(par_likelihood$likelihood),])); colnames(par_stack) <- run_rxn$kineticPars$formulaName
  occupancy_vals <- data.frame(run_rxn$metabolites, par_stack)
  predOcc <- model.matrix(run_rxn$occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
 
  enzyme_activity <- (predOcc %*% t(rep(1, sum(run_rxn$all_species$SpeciesType == "Enzyme"))))*run_rxn$enzymes #occupany of enzymes * relative abundance of enzymes
  flux_fit <- nnls(enzyme_activity, run_rxn$flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  
  fit_summary <- data.frame(residDF = sum(run_rxn$flux != 0) - length(par_stack[1,]), parametricFit = NA, NNLS = NA, LS = NA, LS_met = NA, LS_enzyme = NA, TSS = NA)
  
  ### using flux fitted from the median parameter set, how much variance is explained
  fit_summary$parametricFit = anova(lm(run_rxn$flux ~ flux_fit$fitted))$S[1]
  
  ### using flux fitted using non-negative least squares regression, how much variance is explained ### metabolite abundances are corrected for whether the metabolite is a product (*-1) or reactant (*1)
  NNLSmetab <- run_rxn$metabolites * -1*(rep(1, n_c) %*% t(sapply(colnames(run_rxn$metabolites), function(x){run_rxn$rxnSummary$rxnStoi[names(run_rxn$rxnSummary$rxnStoi) == names(run_rxn$rxnSummary$metsID2tID)[run_rxn$rxnSummary$metsID2tID == x]]})))
    
  if(all(run_rxn$flux < 0)){
    NNLSanova <- anova(lm(-1*run_rxn$flux ~ nnls(as.matrix(data.frame(NNLSmetab, run_rxn$enzymes)), -1*run_rxn$flux)$fitted))
    fit_summary$NNLS <- ifelse(length(NNLSanova$S) != 1, NNLSanova$S[1], NA)
    }else{
      NNLSanova <- anova(lm(run_rxn$flux ~ nnls(as.matrix(data.frame(NNLSmetab, run_rxn$enzymes)), run_rxn$flux)$fitted))
      fit_summary$NNLS = ifelse(length(NNLSanova$S) != 1, NNLSanova$S[1], NA)
    }
  
  
  ### using LS regression, how much variance is explained 
  fit_summary$LS_met = anova(lm(run_rxn$flux ~ run_rxn$metabolites))$S[1]
  fit_summary$LS_enzyme = anova(lm(run_rxn$flux ~ run_rxn$enzymes))$S[1]
  fit_summary$LS = sum(anova(lm(run_rxn$flux ~ run_rxn$metabolites + run_rxn$enzymes))$S[1:2])
  fit_summary$TSS = sum(anova(lm(run_rxn$flux ~ run_rxn$metabolites + run_rxn$enzymes))$S)
  
  output <- list()
    output$fit_summary <- fit_summary
    output$param_interval <- param_interval
    output$fitted_flux <- flux_fit$fitted
  output
  }


param_compare <- function(run_rxn, par_markov_chain, par_likelihood){
  # visualize the joint and marginal distribution of parameter values from the markov chain
  
  par_combinations <- expand.grid(1:length(run_rxn$kineticPars[,1]), 1:length(run_rxn$kineticPars[,1]))
  like_comparison <- ifelse(par_combinations[,1] == par_combinations[,2], TRUE, FALSE)
  
  max_likelihood <- par_markov_chain[which.max(par_likelihood$likelihood),]
  
  par_comp_like <- NULL
  for(i in 1:sum(like_comparison)){
    par_comp_like <- rbind(par_comp_like, data.frame(xval = par_markov_chain[,par_combinations[like_comparison,][i,1]], parameter_1 = colnames(par_markov_chain)[par_combinations[like_comparison,][i,1]],
         parameter_2 = colnames(par_markov_chain)[par_combinations[like_comparison,][i,1]]))
      }
  
  par_comp_dissimilar <- NULL
  for(i in 1:sum(!like_comparison)){
    par_comp_dissimilar <- rbind(par_comp_dissimilar, data.frame(xval = par_markov_chain[,par_combinations[!like_comparison,][i,1]], yval = par_markov_chain[,par_combinations[!like_comparison,][i,2]], 
          parameter_1 = colnames(par_markov_chain)[par_combinations[!like_comparison,][i,1]], parameter_2 = colnames(par_markov_chain)[par_combinations[!like_comparison,][i,2]]))
      }
  
  MLEbarplot <- data.frame(xval = max_likelihood[par_combinations[like_comparison,1]], parameter_1 = colnames(par_markov_chain)[par_combinations[like_comparison,1]],
      parameter_2 = colnames(par_markov_chain)[par_combinations[like_comparison,1]])
  MLEpoints <- data.frame(xval = max_likelihood[par_combinations[!like_comparison,1]], yval = max_likelihood[par_combinations[!like_comparison,2]],
      parameter_1 = colnames(par_markov_chain)[par_combinations[!like_comparison,1]], parameter_2 = colnames(par_markov_chain)[par_combinations[!like_comparison,2]])
  
  
  
  #### determine the maximum bin from the histogram so that values can be scaled to the bivariate histogram values ###

  par_hist_binwidth = 0.2
  
  max_density <- max(apply(par_markov_chain, 2, function(x){max(table(round(x/par_hist_binwidth)))}))
  
  density_trans_inv <- function(x){x*(max_density/20) + max_density/2}
  density_trans <- function(x){(x - max_density/2)/(max_density/20)}
  
  par_comp_dissimilar$yval <- density_trans_inv(par_comp_dissimilar$yval)
  MLEpoints$yval <- density_trans_inv(MLEpoints$yval)
  
  hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "aliceblue"), 
      legend.position = "top", strip.background = element_rect(fill = "cornflowerblue"), strip.text = element_text(color = "cornsilk"), panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line"), axis.title = element_blank()) 

  ggplot() + geom_hex(data = par_comp_dissimilar, aes(x = xval, y = yval)) + geom_bar(data = par_comp_like, aes(x = xval), binwidth = par_hist_binwidth, col = "black") + facet_grid(parameter_2 ~ parameter_1, scales = "fixed") + hex_theme +
    scale_fill_gradientn(name = "Counts", colours = c("white", "darkgoldenrod1", "chocolate1", "firebrick1", "black")) +
    scale_x_continuous(NULL, expand = c(0.02,0.02)) + scale_y_continuous(NULL, expand = c(0.01,0.01), labels = density_trans, breaks = density_trans_inv(seq(-10, 10, by = 5))) +
    geom_vline(data = MLEbarplot, aes(xintercept = xval), col = "cornflowerblue", size = 2) + geom_point(data = MLEpoints, aes(x = xval, y = yval), size = 2, col = "cornflowerblue")

  }

par_draw <- function(updates){
  #### update parameters using their prior (given by kineticParPrior) - update those those parameters whose index is in "updates" ####
  
  draw <- current_pars
  for(par_n in updates){
    if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- runif(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
      } else if(kineticParPrior$distribution[par_n] == "unif"){
      draw[par_n] <- rnorm(1, kineticParPrior$par_1[par_n], kineticParPrior$par_2[par_n])
      }
    }
  draw
  }


lik_calc <- function(proposed_params){
  #### determine the likelihood of predicted flux as a function of metabolite abundance and kinetics parameters relative to actual flux ####
  
  par_stack <- rep(1, n_c) %*% t(proposed_params); colnames(par_stack) <- kineticPars$formulaName
  par_stack <- exp(par_stack)
  occupancy_vals <- data.frame(met_abund, par_stack)
  
  predOcc <- model.matrix(occupancyEq, data = occupancy_vals)[,1] #predict occupancy as a function of metabolites and kinetic constants based upon the occupancy equation
  enzyme_activity <- (predOcc %*% t(rep(1, sum(all_species$SpeciesType == "Enzyme"))))*enzyme_abund #occupany of enzymes * relative abundance of enzymes
  
  flux_fit <- nnls(enzyme_activity, flux) #fit flux ~ enzyme*occupancy using non-negative least squares (all enzymes have activity > 0, though negative flux can occur through occupancy)
  fit_resid_error <- sqrt(mean((flux_fit$resid - mean(flux_fit$resid))^2))
  
  sum(dnorm(flux, flux_fit$fitted, fit_resid_error, log = TRUE))
  
  }

  