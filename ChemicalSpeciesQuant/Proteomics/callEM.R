setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
source("pep_library.R")
matrix_fxn()
load("EMimport.Rdata")
library(limSolve)

#initialize theta (mixing_fraction) and pi (prob non-divergent peptide) and alpha (support for a protein existing)

mixing_fract <- unique_mappingMat/rowSums(unique_mappingMat)
sparse_mapping <- Matrix(unique_mappingMat)

pi_fit <- rep(1, times = n_p)
#record the log-likelihood of a match versus non-match
pi_lik_comp <- data.frame(match = rep(NA, times = n_p), divergent = rep(NA, times = n_p))

#### for alpha, determine which proteins are a subset or subsumable (no non-shared peptides)

sharedPeps <- rowSums(sparse_mapping) != 1
subSumable <- ifelse(colSums(sparse_mapping[!sharedPeps,]) == 0, 1, 0)
alpha_pres <- rep(1, times = n_prot)

#decompose overlapping peptides into sub-bipartite graphs
subGraphs <- list()
relFract <- list()
#identify overlapping proteins and seperate by ORF name
overlapMat <- t(sparse_mapping) %*% sparse_mapping
overlapMat <- overlapMat[rowSums(overlapMat != 0) != 1,rowSums(overlapMat != 0) != 1]

entry_split <- c()
entries_remaining <- c(1:length(overlapMat[,1]))
while(length(entries_remaining) != 0){
  members <- entries_remaining[1]
  nmembersChange <- 1
  
  while(nmembersChange != 0){
    new_members <- c(1:length(overlapMat[,1]))[rowSums(matrix(overlapMat[members,], nrow = length(overlapMat[,1]), byrow = TRUE)) != 0]
    nmembersChange <- length(new_members) - length(members)
    members <- new_members
    }
  
  entries_remaining <- entries_remaining[!(entries_remaining %in% members)]
  entry_split <- c(entry_split, paste(members, collapse = "_"))
}

for(entry in 1:length(entry_split)){
  subGraphs[[entry]] <- rownames(overlapMat)[as.numeric(strsplit(entry_split[entry], split = "_")[[1]])]
  relFract[[entry]] <- rep((1/length(subGraphs[[entry]])), times =  length(subGraphs[[entry]]))  
}


#plot the shared peptides across groups of proteins
plot_proteinConnections(overlapMat)

#save the estimates of mixing fractions generated on a per-peptide basis
mixing_fract_inconsistent <- mixing_fract
ambig_peps <- c(1:n_p)[rowSums(sparse_mapping) != 1]  
  

#override determination - switch the protein abundance with the most consistent single peptide

prot_abund_over <- matrix(NA, ncol = n_prot, nrow = n_c)
for(a_prot in c(1:n_prot)){
  condMed <- apply(matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c), 1, median)
  bestPepMatch <- matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)[,which.max(cor(condMed, matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)))]
  bestPepMatch[bestPepMatch == 0] <- condMed[bestPepMatch == 0]
  prot_abund_over[,a_prot] <- bestPepMatch
}
  
### Iteration ###

whole_data_logL <- NULL
previous_it <- -Inf
continue_it <- TRUE
set_alpha <- FALSE # determines which proteins can be dropped at the time of convergence

t <- 1
while(continue_it){
	
	#update protein abundances: using integrated likelihood
	#update protein means (c x j)
	prot_abund = ((uniquePepMean*uniquePepPrecision) %*% mixing_fract)/(uniquePepPrecision %*% mixing_fract)
	prot_abund <- as.matrix(Matrix(prot_abund, sparse = FALSE))
	prot_abund[is.nan(prot_abund)] <- 0
	
  #if all peptides attributed to a protein are called "divergent" then initialize its trend as the median of matched peptides for each condition
	divOverride <- colSums(mixing_fract * pi_fit) == 0; divOverride[alpha_pres == 0] <- FALSE
	prot_abund[,divOverride] <- prot_abund_over[,divOverride]
  
	#update protein precision
	prot_prec <- uniquePepPrecision %*% mixing_fract
	
	#update pi - relative prob peptide match based on bayes factors
  
  prot_matchLik <- -(1/2)*uniquePepPrecision*(((prot_abund %*% diag(alpha_pres) %*% t(mixing_fract)) - uniquePepMean)^2)
	pi_lik_comp$match <- colSums(prot_matchLik); pi_lik_comp$divergent <- log(prior_p_div)
  #table(apply(pi_lik_comp, 1, diff) < 0)
	pi_fit <- ifelse(apply(pi_lik_comp, 1, diff) < 0, 1, 0)
  
  #update mixing fract based on weighted quadratic programming
	
  for(peptide in ambig_peps){
    n_matches <- sum(sparse_mapping[peptide,]*alpha_pres)
    solution <- lsei(A = prot_abund[,sparse_mapping[peptide,]*alpha_pres == 1]*uniquePepPrecision[,peptide], B = uniquePepMean[,peptide] * uniquePepPrecision[,peptide], E = t(rep(1, times = n_matches)), F = 1, G = diag(n_matches), H = rep(0, times = n_matches))$X
    mixing_fract_inconsistent[peptide,sparse_mapping[peptide,]*alpha_pres == 1] <- solution
    }
  
  #update mixing fractions based on weighted quadratic programming over each isolated bipartite graph
  for(graph_part in 1:length(subGraphs)){
    
    relProt <- colnames(sparse_mapping) %in% subGraphs[[graph_part]]
    #remove invalid proteins from consideration (those with alpha = 0)
    relProt[relProt] <- alpha_pres[relProt] == 1
    
    if(sum(relProt) <= 1){
      mixing_fract[,relProt][sparse_mapping[,relProt] != 0] <- 1
      next
      }
    
    relPep <- rowSums(sparse_mapping[,relProt] != 0) != 0
    
    subProtAbund <- prot_abund[,relProt]
    subPepAbund <- uniquePepMean[,relPep]
    subPepPrec <- uniquePepPrecision[,relPep]
    subMap <- sparse_mapping[relPep, relProt]
    
    #attempt 2
    #find mixing proporitons for observed combination of peptide-protein matches
    
    pepCombos <- apply(subMap, 1, function(x){paste(x, collapse = "")})
    
    uniqueCombos <- unique(pepCombos[apply(subMap, 1, sum) != 1])
    
    combo_mixes <- matrix(0, nrow = length(uniqueCombos), ncol = length(subMap[1,])); rownames(combo_mixes) <- uniqueCombos; colnames(combo_mixes) <- colnames(subMap); consistentMixes <- combo_mixes
    combo_precision <- c()
    for(a_combo in 1:length(uniqueCombos)){
      comboMatches <- pepCombos %in% uniqueCombos[a_combo]
      matchedCombo <- as.numeric(strsplit(uniqueCombos[a_combo], split = "")[[1]])
      
      protStack <- NULL
      pepVec <- NULL
      precTotal <- 0
      
      for(pepEntries in c(1:length(subMap[,1]))[comboMatches]){
        #protein trends weighted by peptide precisions and poriton of total mixing fraction for matched proteins in this peptide
        protStack <- rbind(protStack, prot_abund[,relProt] * t(t(rep(1, n_c))) %*% t(subMap[pepEntries,]) * t(t(subPepPrec[,pepEntries])) %*% rep(1, length(matchedCombo)))
        #peptide trends weighted by peptide precisions
        pepVec <- c(pepVec, subPepAbund[,pepEntries] * subPepPrec[,pepEntries])
        
        precTotal <- precTotal + sum(subPepPrec[,pepEntries])
      }
      
      protStack <- protStack[,matchedCombo == 1]
      
      solution <- lsei(A = protStack, B = pepVec, E = rep(1, sum(matchedCombo)), F = 1, G = diag(rep(1, sum(matchedCombo))), H = rep(0, times = sum(matchedCombo)))$X
      combo_mixes[a_combo,matchedCombo==1] <- solution
      combo_mixes[a_combo,matchedCombo==0] <- NA
      combo_precision <- c(combo_precision, precTotal)
      }
    
    #find consistent mixtures across combinations, weighting by the precision of estimation
      if(length(uniqueCombos) == 1){
        comboMatches <- pepCombos %in% uniqueCombos
        matchedCombo <- as.numeric(strsplit(uniqueCombos, split = "")[[1]])
        mixing_fract[relPep, relProt][comboMatches,matchedCombo == 1] <- t(t(rep(1, sum(comboMatches)))) %*% combo_mixes
      
        }else{
        for(a_combo in 1:length(uniqueCombos)){
          comboMatches <- pepCombos %in% uniqueCombos[a_combo]
          matchedCombo <- as.numeric(strsplit(uniqueCombos[a_combo], split = "")[[1]])
          subMixes <- combo_mixes[,matchedCombo == 1]/rowSums(combo_mixes[,matchedCombo == 1], na.rm = TRUE)
          wsubMix <- subMixes * combo_precision
          
          mix_converge <- FALSE
          
          pseudoRowTotal <- rowSums(wsubMix, na.rm = TRUE)
          pseudoFracts <- wsubMix
          pi_tot_row <- rep(1, times = length(uniqueCombos))
          pi_prop <- rep(1/length(wsubMix[1,]), length(wsubMix[1,]))
          
          while(mix_converge == FALSE){
            for(i in 1:length(wsubMix[,1])){
              for(j in 1:length(wsubMix[1,])){
                if(is.na(wsubMix[i,j])){
                  pseudoFracts[i,j] <- pseudoRowTotal[i] * pi_prop[j]
                }
              }
            }
            pseudoRowTotal <- rowSums(pseudoFracts)
            
            new_pi <- colSums(pseudoFracts, na.rm = TRUE)/sum(colSums(pseudoFracts, na.rm = TRUE))
            if(sum((pi_prop - new_pi)^2) < 10^-12){mix_converge <- TRUE}
            pi_prop <- new_pi
          }
          mixing_fract[relPep, relProt][comboMatches,matchedCombo == 1] <- t(t(rep(1, sum(comboMatches)))) %*% pi_prop
        }
      }
  }  
    
    
  if(set_alpha == TRUE){ 
  
    #update alpha - whether a protein is present, based on bayes factors - do this once upon initial convergence and then continue running
    
    for(protein in c(1:n_prot)[subSumable == 1]){
      #find all peptides that match the possibly absent protein
      
      alpha_neg <- alpha_pres; alpha_neg[protein] <- 0
      alpha_pos <- alpha_pres; alpha_pos[protein] <- 1
      
      lpos <- 0
      lneg <- 0
      
      for(peptide in c(1:n_p)[sparse_mapping[,protein] == 1]){
        
        n_matches_pos <- sum(sparse_mapping[peptide,]*alpha_pos)
        lpos <- lpos + sum((apply(prot_abund[,sparse_mapping[peptide,]*alpha_pos == 1]*matrix(mixing_fract[peptide,sparse_mapping[peptide,]*alpha_pos == 1], ncol = n_matches_pos, nrow = n_c, byrow = TRUE), 1, sum)
        - uniquePepMean[,peptide])^2*uniquePepPrecision[,peptide]*(-1/2))
                            
        n_matches_neg <- sum(sparse_mapping[peptide,]*alpha_neg)
        #rescale the negative mixing fractions to account for the fraction lost due to one proteins invalidation
        neg_mixing <- mixing_fract[peptide,sparse_mapping[peptide,]*alpha_neg == 1]
        if(sum(neg_mixing) == 0){
          neg_mixing <- rep(1/length(neg_mixing), times = length(neg_mixing))
          }
        neg_mixing <- neg_mixing/sum(neg_mixing)
        
        lneg <- lneg + sum((apply(prot_abund[,sparse_mapping[peptide,]*alpha_neg == 1]*matrix(neg_mixing, ncol = n_matches_neg, nrow = n_c, byrow = TRUE), 1, sum)
                            - uniquePepMean[,peptide])^2*uniquePepPrecision[,peptide]*(-1/2))
          }
        if(lneg >= lpos + log(prior_p_div)){
          alpha_pres[protein] <- 0 
          mixing_fract[,protein] <- 0
        }
      }
    set_alpha <- "SET"
    print(paste(sum(alpha_pres == 0), "proteins were removed"))
    }
  
	#update complete data log-likelihood
	new_log_lik <- sum(apply(-(1/2)*uniquePepPrecision*(((prot_abund %*% diag(alpha_pres) %*% t(mixing_fract)) - uniquePepMean)^2), 2, sum)*pi_fit + (1-pi_fit)*log(prior_p_div)) + sum(alpha_pres[subSumable == 1])*log(prior_p_div)
  
	#check for convergence
	
	if(abs(new_log_lik - previous_it) < 0.01 | (t > 50 + 50*(set_alpha == "SET"))){
	  
    if(set_alpha == FALSE){
      set_alpha <- TRUE
      whole_data_logL <- c(whole_data_logL, new_log_lik)
      t <- t + 1
      print("removing unsupported proteins")
    }else{
      continue_it <- FALSE
      whole_data_logL <- c(whole_data_logL, new_log_lik)
      t <- t + 1
      print("done")
    }
    }else{
  	  whole_data_logL <- c(whole_data_logL, new_log_lik)
  	  print(paste("iteration", t, "log likelihood", round(new_log_lik, 3)))
  	  t <- t + 1
  	  previous_it <- new_log_lik
  	}
  #plot(c(prot_abund) ~ c(prot_abund_over), cex = 0.5, pch = 16)
}

max_state <- apply(mixing_fract, 1, which.max)
div_max <- pi_fit == 0

#number of peptides not-conforming to some general protein trend
table(div_max)

# for each model compare the best peptide-protein match with the non-match likelihood*prior

likdiff <- pi_lik_comp$match - pi_lik_comp$divergent
names(likdiff) <- rownames(mixing_fract)


likdiff_df <- data.frame(likelihood = likdiff, matched = ifelse(likdiff >= 0, "Protein-match", "Divergent-trend"))
likdiff_df$likelihood[likdiff_df$likelihood <= -400] <- -400


library(ggplot2)
likdiff_plot <- ggplot(likdiff_df, aes(x = likelihood, fill = matched))
likdiff_plot + geom_histogram()


#### are mixing fraction for proteins with more than 1 overlapping peptide comparable #####
#generate a square matrix of overlap
prot_overlap <- t(sparse_mapping) %*% sparse_mapping
mixed_pairs <- triu(prot_overlap); diag(mixed_pairs) <- 0
mixed_pairs2 <- NULL
for(prot in c(1:n_prot)[rowSums(mixed_pairs > 3) != 0]){
  matches <- c(1:n_prot)[mixed_pairs[prot,] > 3]
  matchN <- mixed_pairs[prot,][mixed_pairs[prot,] > 3]
  mixed_pairs2 <- rbind(mixed_pairs2, data.frame(protA = prot, protB = matches, nMatches = matchN))  
  }
mixed_pairs2 <- mixed_pairs2[order(mixed_pairs2$nMatches, decreasing = TRUE),]

mixSD <- rep(NA, times = length(mixed_pairs2[,1]))
mixHL <- data.frame(mixDom = rep(NA, times = length(mixed_pairs2[,1])), N = NA)
boxplot_list <- list()
boxplot_df <- NULL
consistent_est <- NULL

for(pairs in 1:length(mixed_pairs2[,1])){
  shared_peps <- c(1:n_p)[rowSums(sparse_mapping[,c(mixed_pairs2$protA[pairs], mixed_pairs2$protB[pairs])]) == 2]
  paired_mixes <- mixing_fract_inconsistent[shared_peps,c(mixed_pairs2$protA[pairs], mixed_pairs2$protB[pairs])]
  paired_mixes_consistent <- mixing_fract[shared_peps,c(mixed_pairs2$protA[pairs], mixed_pairs2$protB[pairs])]
  boxplot_df <- rbind(boxplot_df, data.frame(Pair = pairs, FractionB = paired_mixes[,2], SharedPeptides = length(paired_mixes[,1])))
  consistent_est <- rbind(consistent_est, data.frame(Pair = pairs, FractionB = paired_mixes_consistent[1,2]/sum(paired_mixes_consistent[1,c(1:2)])))
  }

mixingFracPlot <- ggplot(boxplot_df)
mixingFracPlot + geom_boxplot(aes(x = factor(Pair), y = FractionB, fill = SharedPeptides, colour = "#3366FF"), position = "dodge") + scale_color_discrete(guide = "none") +
  scale_x_discrete('Protein pairs with more than 3 overlapping peptides') + scale_y_continuous('Mixing proprotion of second protein of pair') +
  geom_point(data = consistent_est, aes(x = factor(Pair), y = FractionB), size = 3, colour = I('chartreuse')) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

library(reshape)

ambigProteinPatterns <- NULL
for(pairs in 31:40){
  
  abunion <- c(1:n_p)[rowSums(sparse_mapping[,c(mixed_pairs2$protA[pairs], mixed_pairs2$protB[pairs])]) != 0]
  
  proteinAtt <- apply(sparse_mapping[abunion,c(mixed_pairs2$protA[pairs], mixed_pairs2$protB[pairs])], 1, function(x){
    if(sum(x) == 2){"SHARED"}else{
      if(x[1] == 1){"PROTEIN A"}else{"PROTEIN B"}
      }
    })
  proteinAtt[div_max[abunion]] <- "DIVERGENT"
  
  peptidepatterns <- melt(uniquePepMean[,abunion])
  colnames(peptidepatterns) <- c("condition", "peptide", "abundance")
  ambigProteinPatterns <- rbind(ambigProteinPatterns, cbind(pair = pairs, peptidepatterns, owner = rep(proteinAtt, each = n_c)))
  }

ambigProteinPlot <- ggplot(ambigProteinPatterns, aes(x = condition, y = abundance, colour = owner, alpha = 0.1, size = 3))
ambigProteinPlot + geom_point() + facet_wrap( ~ pair, ncol = 2, scales = "free") + scale_colour_manual(values = c("SHARED" = "ORANGE","PROTEIN A" = "RED","PROTEIN B" = "YELLOW", "DIVERGENT" = "PURPLE"))


#for each peptide compare a model where peptide patterns are equivalent to protein patterns with one where a peptide doesn't match - the divergent model will fit perfectly vs. the alternative where the protein fits

divergentPep_summary <- data.frame(peptide = names(likdiff)[div_max], prot_matches = NA, likelihood = unname(likdiff[div_max]), nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)

poor_match_redMat <- sparse_mapping[div_max,]

for(pep in 1:length(poor_match_redMat[,1])){
	divergentPep_summary$prot_matches[pep] <- paste(colnames(poor_match_redMat)[poor_match_redMat[pep,] == 1], collapse = "/")
	divergentPep_summary$nS[pep] <- length(grep("S", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	divergentPep_summary$nR[pep] <- length(grep("R", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	divergentPep_summary$nY[pep] <- length(grep("Y", unlist(strsplit(divergentPep_summary$peptide[pep], ""))))
	}
divergentPep_summary$possible_phosphoSite <- divergentPep_summary$nS != 0 | divergentPep_summary$nR != 0 | divergentPep_summary$nY != 0


background_SRYfreq <- data.frame(peptide = names(likdiff)[!div_max], prot_matches = NA, likelihood = unname(likdiff[!div_max]), nS = NA, nR = NA, nY = NA, stringsAsFactors = FALSE)

good_match_redMat <- sparse_mapping[!div_max,]

for(pep in 1:length(good_match_redMat[,1])){
	background_SRYfreq$prot_matches[pep] <- paste(colnames(good_match_redMat)[good_match_redMat[pep,] == 1], collapse = "/")
	background_SRYfreq$nS[pep] <- length(grep("S", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	background_SRYfreq$nR[pep] <- length(grep("R", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	background_SRYfreq$nY[pep] <- length(grep("Y", unlist(strsplit(background_SRYfreq$peptide[pep], ""))))
	}
background_SRYfreq$possible_phosphoSite <- background_SRYfreq$nS != 0 | background_SRYfreq$nR != 0 | background_SRYfreq$nY != 0



#how many peptides are informing a protein trend, with only the largest effect considered
barplot(table(table(max_state)), xlab = "Number of peptides for which a protein dominates pattern", ylab = "Frequency")


if(unq_matches_only == TRUE){
	save(prot_abund, prot_prec, sparse_mapping, mixing_fract, whole_data_logL, max_state, div_max, divergentPep_summary, background_SRYfreq, likdiff_df, pi_fit, alpha_pres, file = "EMoutputUnq.Rdata")
	}else{
		save(prot_abund, prot_prec, sparse_mapping, mixing_fract, whole_data_logL, max_state, div_max, divergentPep_summary, background_SRYfreq, likdiff_df, pi_fit, alpha_pres, file = "EMoutputDeg.Rdata")
	}


