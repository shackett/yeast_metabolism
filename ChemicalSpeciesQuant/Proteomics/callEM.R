setwd('~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics')
source("pep_library.R")
matrix_fxn()
load("EMimport.Rdata")
library(limSolve)
library(data.table)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = F)

options(digits = 15)


#prior_p_div <- exp(-1*qchisq(prior_bound, n_c))  ### this should probably be adjusted to the 
n_pep_samples <- colSums(uniquePepMean != 0) 
prior_bound <- 1 - 10^-9
prior_p_div <- exp(-1*qchisq(prior_bound, n_pep_samples)) 
prior_prot_sub <- exp(-1*qchisq(prior_bound, n_c)) 



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



#save the estimates of mixing fractions generated on a per-peptide basis
mixing_fract_inconsistent <- mixing_fract
ambig_peps <- c(1:n_p)[rowSums(sparse_mapping) != 1]  
  


### Iteration ###

whole_data_logL <- NULL
previous_it <- -Inf
continue_it <- TRUE
set_alpha <- FALSE # determines which proteins can be dropped at the time of convergence
param_track <- list()
param_track$mixing_fract[[1]] <- as.matrix(mixing_fract) #for some reason initializing with a sparse matrix doesn't work, but overwriting the first object with a sparse matrix makes the others behave

t <- 1
while(continue_it){
	
	#update protein abundances: using integrated likelihood
	#update protein means (c x j)
	
  prot_abund = ((uniquePepMean*uniquePepPrecision) %*% Diagonal(n_p, pi_fit) %*% mixing_fract)/(uniquePepPrecision %*% Diagonal(n_p, pi_fit) %*% mixing_fract)
	prot_abund <- as.matrix(Matrix(prot_abund, sparse = FALSE))
	prot_abund[is.nan(prot_abund)] <- 0
	
  
  #if all peptides attributed to a protein are called "divergent" then initialize its trend as the median of matched peptides for each condition
	#divOverride <- colSums(mixing_fract * pi_fit) == 0; divOverride[alpha_pres == 0] <- FALSE
	#prot_abund[,divOverride] <- prot_abund_over[,divOverride]
  
	#update protein precision
	prot_prec <- uniquePepPrecision %*% Diagonal(n_p, pi_fit) %*% mixing_fract
	
	#update pi - relative prob peptide match based on bayes factors
  
  prot_matchLik <- -(1/2)*uniquePepPrecision*(((prot_abund %*% Diagonal(n_prot, alpha_pres) %*% t(mixing_fract)) - uniquePepMean)^2)
	pi_lik_comp$match <- colSums(prot_matchLik); pi_lik_comp$divergent <- log(prior_p_div)
  pi_fit <- ifelse(apply(pi_lik_comp, 1, diff) < 0, 1, 0)
  
  ### update mixing fract based on weighted quadratic programming ### 
  ### This gives a solution to the mixing fraction which is specific to mixing multiple proteins to match a single peptide rather than all peptides shared by the proteins
	### This is saved for comparison to the consistent mixtures but is not directly used
  
  for(peptide in ambig_peps){
    n_matches <- sum(sparse_mapping[peptide,]*alpha_pres)
    solution <- lsei(A = prot_abund[,sparse_mapping[peptide,]*alpha_pres == 1]*uniquePepPrecision[,peptide], B = uniquePepMean[,peptide] * uniquePepPrecision[,peptide], E = t(rep(1, times = n_matches)), F = 1, G = diag(n_matches), H = rep(0, times = n_matches))$X
    mixing_fract_inconsistent[peptide,sparse_mapping[peptide,]*alpha_pres == 1] <- solution
    }
  
  ### update mixing fractions based on weighted quadratic programming over each isolated bipartite graph ####
  for(graph_part in 1:length(subGraphs)){
    
    relProt <- colnames(sparse_mapping) %in% subGraphs[[graph_part]]
    ### remove invalid proteins from consideration (those with alpha = 0)
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
    #mixing_fract[relPep, relProt]
    ### Find mixing proportions for observed combination of peptide-protein matches
    
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
        ### protein trends weighted by peptide precisions and portion of total mixing fraction for matched proteins in this peptide
        protStack <- rbind(protStack, prot_abund[,relProt] * t(t(rep(1, n_c))) %*% t(subMap[pepEntries,]) * t(t(subPepPrec[,pepEntries])) %*% rep(1, length(matchedCombo)))
        ### peptide trends weighted by peptide precisions
        pepVec <- c(pepVec, subPepAbund[,pepEntries] * subPepPrec[,pepEntries])
        
        precTotal <- precTotal + sum(subPepPrec[,pepEntries])
      }
      
      protStack <- protStack[,matchedCombo == 1]
      
      solution <- lsei(A = protStack, B = pepVec, E = rep(1, sum(matchedCombo)), F = 1, G = diag(rep(1, sum(matchedCombo))), H = rep(0, times = sum(matchedCombo)))$X
      combo_mixes[a_combo,matchedCombo==1] <- solution
      combo_mixes[a_combo,matchedCombo==0] <- NA
      combo_precision <- c(combo_precision, precTotal)
      }
    
    ### Find consistent mixtures across combinations, weighting by the precision of estimation
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
    #for peptides which are uncontested set mixing fraction at 1 - this corrects for the point when some peptides become unambiguous after degenerate peptides are removed 
    mixing_fract[relPep, relProt][rowSums(sparse_mapping[relPep, relProt]) == 1,] <- sparse_mapping[relPep, relProt][rowSums(sparse_mapping[relPep, relProt]) == 1,]
  }  
    
    
  if(set_alpha == TRUE){ 
  
    #set parameters as the set with the highest likelihood
    
    param_track$alpha_pres[[which.max(whole_data_logL)]] -> alpha_pres
    param_track$mixing_fract[[which.max(whole_data_logL)]] -> mixing_fract
    param_track$prot_abund[[which.max(whole_data_logL)]] -> prot_abund
    param_track$pi_fit[[which.max(whole_data_logL)]] -> pi_fit
    print(paste(c("overwriting parameters with those from iteration", c(1:length(whole_data_logL))[which.max(whole_data_logL)]), collapse = " "))
    
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
        if(lneg >= lpos + log(prior_prot_sub)){
          alpha_pres[protein] <- 0 
          mixing_fract[,protein] <- 0
        }
      }
    set_alpha <- "SET"
    print(paste(sum(alpha_pres == 0), "proteins were removed"))
    }
  
	### update complete data log-likelihood
	new_log_lik <- sum(apply(-(1/2)*uniquePepPrecision*(((prot_abund %*% Diagonal(n_prot, alpha_pres) %*% t(mixing_fract)) - uniquePepMean)^2), 2, sum)*pi_fit + (1-pi_fit)*log(prior_p_div)) + sum(alpha_pres[subSumable == 1])*log(prior_prot_sub)
  
  ### Store parameter sets corresponding to the likelihood at each iteration
  
  param_track$alpha_pres[[t]] <- alpha_pres
  param_track$mixing_fract[[t]] <- mixing_fract#as.matrix(mixing_fract)
  param_track$prot_abund[[t]] <- prot_abund
  param_track$pi_fit[[t]] <- pi_fit
  
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
      
      param_track$alpha_pres[[which.max(whole_data_logL)]] -> alpha_pres
      param_track$mixing_fract[[which.max(whole_data_logL)]] -> mixing_fract
      param_track$prot_abund[[which.max(whole_data_logL)]] -> prot_abund
      param_track$pi_fit[[which.max(whole_data_logL)]] -> pi_fit
      print(paste(c("overwriting parameters with those from iteration", c(1:length(whole_data_logL))[which.max(whole_data_logL)]), collapse = " "))
      print("done")
    }
    }else{
  	  whole_data_logL <- c(whole_data_logL, new_log_lik)
  	  print(paste("iteration", t, "log likelihood", round(new_log_lik, 3)))
  	  t <- t + 1
  	  previous_it <- new_log_lik
  	}
}
save(list = ls(), file = "tmp2.Rdata")
#load("tmp.Rdata")


### Rescue divergent peptides which have a high correlation with protein, but may be have innappropriately estimated precision or have a different offset between H/L than other peptides (e.g. degradation) ####
 
comp_protein <- t(prot_abund %*% t(mixing_fract[pi_fit == 0,]))
divergent_peptides <- t(uniquePepMean[,pi_fit == 0])
comp_protein[comp_protein == 0] <- NA; divergent_peptides[divergent_peptides == 0] <- NA

ppCor <- sapply(c(1:length(comp_protein[,1])), function(x){cor(comp_protein[x,], divergent_peptides[x,], use = "pairwise.complete.obs")})
pi_fit[pi_fit == 0][ppCor > 0.6 & !is.na(ppCor)] <- 1

print(paste(c("Rescued", sum(ppCor > 0.6 & !is.na(ppCor)), "divergent peptides which were moderately correlated with protein"), collapse = " "))

# Recalculate protein trends and precision with newly included peptides #

prot_abund = ((uniquePepMean*uniquePepPrecision) %*% Diagonal(n_p, pi_fit) %*% mixing_fract)/(uniquePepPrecision %*% Diagonal(n_p, pi_fit) %*% mixing_fract)
prot_abund <- as.matrix(Matrix(prot_abund, sparse = FALSE))
prot_abund[is.nan(prot_abund)] <- 0

prot_prec <- uniquePepPrecision %*% Diagonal(n_p, pi_fit) %*% mixing_fract
	

### for non-subsumable proteins - check to see if zero peptides carry weight - if so take the median abundance and the corresponding precision ####
## overright a protein using the peptide with the highest correlation to the median profile
## for missing values, overwrite by the median and corresponding precision

#index_check <- (1:n_p)[rowSums(mixing_fract) == min(rowSums(mixing_fract))]
#index_check2 <- colnames(sparse_mapping)[colSums(sparse_mapping[index_check,]) != 0]

NinformedPeps <- colSums(mixing_fract[pi_fit == 1, alpha_pres == 1])

print(paste(c(sum(NinformedPeps == 0), "proteins reset to most representative peptide"), collapse = " "))

prot_overwrites <- c(1:n_prot)[alpha_pres == 1][NinformedPeps == 0]
#table(sparse_mapping[,prot_overwrites])
#pi_fit[sparse_mapping[,prot_overwrites] == 1]
#mixing_fract[sparse_mapping[,prot_overwrites] == 1,prot_overwrites]

if(length(prot_overwrites) != 0){

  prot_abund_over <- prot_prec_over <- matrix(NA, ncol = length(prot_overwrites), nrow = n_c)
  index_overwrite <- rep(NA, length(prot_overwrites))
  for(a_prot in prot_overwrites){
    condMed <- apply(matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c), 1, median) #across all matched peptides which has the median relative abundance
    condMetPrec <- sapply(1:n_c, function(acond){
      min_dist <- signif(abs(matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)[acond,] - condMed[acond]), 6)
      1/mean(1/matrix(uniquePepPrecision[,sparse_mapping[,a_prot] == 1], nrow = n_c)[acond,][min_dist == min(min_dist)]) #inverse variance of the mean
    })
    bestPepMatch <- matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)[,which.max(cor(condMed, matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)))] #use correlation between each peptide and median vector to determine which peptide is most representatitive
    bestPepPrec <- matrix(uniquePepPrecision[,sparse_mapping[,a_prot] == 1], nrow = n_c)[,which.max(cor(condMed, matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)))] #the above peptides precisio
    bestPepMatch[bestPepMatch == 0] <- condMed[bestPepMatch == 0]
    bestPepPrec[bestPepMatch == 0] <- condMetPrec[bestPepMatch == 0]
    
    prot_abund_over[,a_prot == prot_overwrites] <- bestPepMatch
    prot_prec_over[,a_prot == prot_overwrites] <- bestPepPrec
    index_overwrite[a_prot == prot_overwrites] <- c(1:n_p)[sparse_mapping[,a_prot] == 1][which.max(cor(condMed, matrix(uniquePepMean[,sparse_mapping[,a_prot] == 1], nrow = n_c)))]
  }
  
  prot_abund[,alpha_pres == 1][,NinformedPeps == 0] <- prot_abund_over
  prot_prec[,alpha_pres == 1][,NinformedPeps == 0] <- prot_prec_over
  pi_fit[index_overwrite] <- 1
  mixing_fract[index_overwrite, prot_overwrites] <- 1
}
###

max_state <- apply(mixing_fract, 1, which.max)
div_max <- pi_fit == 0

missingValPrec <- data.table(divergent = div_max, precision = apply(uniquePepPrecision, 2, mean))
ggplot(missingValPrec, aes(x = log2(precision))) + facet_wrap(~ divergent) + geom_bar()
ggplot(missingValPrec, aes(y = log2(precision), x = factor(divergent))) + geom_violin()

#plot the shared peptides across groups of proteins

plot_proteinConnections(overlapMat, subsumedIDs = colnames(prot_abund)[alpha_pres == 0])



# for each model compare the best peptide-protein match with the non-match likelihood*prior

likdiff <- pi_lik_comp$match - pi_lik_comp$divergent
names(likdiff) <- rownames(mixing_fract)


likdiff_df <- data.frame(likelihood = likdiff, matched = ifelse(likdiff >= 0, "Protein-match", "Divergent-trend"))
likdiff_df$matched[likdiff_df$matched == "Divergent-trend"][pi_fit[likdiff_df$matched == "Divergent-trend"] == 1] <- "Correlation-overwrite"

likdiff_df$likelihood[likdiff_df$likelihood <= -400] <- -400


library(ggplot2)
likdiff_plot <- ggplot(likdiff_df, aes(x = likelihood, fill = matched))
likdiff_plot + geom_histogram(binwidth = 5)


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

hex_theme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "gray90"), legend.position = "top", 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.line = element_blank(), legend.key.width = unit(6, "line")) 


mixingFracPlot <- ggplot(boxplot_df)
mixingFracPlot + geom_boxplot(aes(x = factor(Pair), y = FractionB, fill = SharedPeptides, colour = "#3366FF"), position = "dodge") + scale_color_discrete(guide = "none") +
  scale_x_discrete('Protein pairs with more than 3 overlapping peptides') + scale_y_continuous('Mixing proprotion of second protein of pair', expand = c(0.005,0.005)) +
  geom_point(data = consistent_est, aes(x = factor(Pair), y = FractionB), size = 3, colour = I('chartreuse')) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + hex_theme



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


