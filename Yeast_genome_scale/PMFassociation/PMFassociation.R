#### Assessing the extent to which flux can be simply explained through either :
# A) Trade-offs between metabolite abundances and enzyme levels maintain robustness in flux
# B) Rate-limiting enzymes or flux-sensing metabolites

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/PMFassociation")

# libraries and sources
library(qvalue)
library(ggplot2)
library(grid) # units
library(scales) # percent
library(data.table)
library(gridExtra) # for combining ggplot objects
library(reshape2)
library(igraph)
library(org.Sc.sgd.db)

source("../FBA_lib.R")

options(stringsAsFactors = F)

scatterPlotTheme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(color = "black", fill = "white"), legend.position = "none",
                          axis.title.x = element_text(vjust = -0.3), axis.title.y = element_text(vjust = 0.25),
                          panel.grid.minor = element_blank(), legend.key.width = unit(6, "line"), panel.grid.major = element_line(colour = "black"), axis.ticks = element_line(colour = "black"), strip.background = element_rect(fill = "cyan"),
                          axis.text = element_text(color = "blacK"),
                          panel.margin = unit(2, "lines")) 

n_example_plots <- 6

c2o <- toTable(org.Sc.sgdCOMMON2ORF) # mapping between systematic and common yeast names


##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@
##### Trade-offs between metabolite abundances and enzyme levels maintain robustness in flux ######
##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@

load("../flux_cache/paramCI.Rdata")
load("../paramOptim.Rdata")

JME_plot <- function(J,M,E,comp_summary){
  
  JMEdf <- data.frame(J = J, ME = M + E)
  JMEdf <- JMEdf - rep(1, nrow(JMEdf)) %*% t(apply(JMEdf, 2, mean)) # center flux and M*E
  
  # aesthetics
  
  sampleInfo <- data.frame(condition = rownames(JMEdf), limitation = substr(rownames(JMEdf), 1, 1), DR = substr(rownames(JMEdf), 2, 5))
  sampleInfo$limitation <- factor(sampleInfo$limitation, levels = unique(sampleInfo$limitation))
  sampleInfo$size <- sqrt(as.numeric(sampleInfo$DR)*200)
  
  JMEdf <- cbind(JMEdf, sampleInfo)
  
  ggplot(JMEdf, aes(x = ME, y = J, color = limitation, size = size)) + geom_point() + geom_abline(intercept = 0, slope = 1) +
    scatterPlotTheme + scale_color_brewer(palette = "Set2") + scale_size_identity() +
    scale_y_continuous(nameReformat(comp_summary$rx, 30)) + scale_x_continuous(paste(comp_summary$enzyme_common, comp_summary$metabolite, sep = " & "))
  
}

MEreactionMechs <- reactionInfo$rMech[reactionInfo$modification == ""]

MEreactionData <- rxnList_form[names(rxnList_form) %in% MEreactionMechs]

# all flux > 0 | all flux < 0
# dealing with multiple E or M : all pair-wise

# two comparisons for every paired substrate and enzyme
# 1 ) log(J) ~ log(E) + log(M) - arbitrary dependence via bootstrapped pearson correlation
# 2 ) log(J) propto log(E) + log(M) - as implied by mass-action kinetics
# % variance explained from explained SS / TSS | slope = 1
# p-value via permutation testing 

nbs_perm <- 10^6 # runtime ~6h
#nbs_perm <- 10000
ME_assoc <- NULL
ME_plot <- list()

#i <- which(names(MEreactionData) == "r_0195-rm")

for(i in 1:length(MEreactionData)){
  
  a_reactionData <- MEreactionData[[i]]
  
  relFlux <- (a_reactionData$flux$FVAmax + a_reactionData$flux$FVAmin)/2
  
  if(!(all(rowSums(a_reactionData$flux > 0) == 3) | all(rowSums(a_reactionData$flux < 0) == 3))){
    # if reaction net flux direction depends on condition then skip reaction
    # log is taken to stabilize variance and remove multiplicative constant
    next    
  }
  
  if(all(relFlux > 0)){
    rxnSubstrates <- a_reactionData$rxnMet[,a_reactionData$rxnStoi < 0, drop = F]
  }else{
    rxnSubstrates <- a_reactionData$rxnMet[,a_reactionData$rxnStoi > 0, drop = F]
    relFlux <- relFlux * -1
  }
  # drop non-measured species
  rxnSubstrates <- rxnSubstrates[,colSums(is.na(rxnSubstrates)) != nrow(rxnSubstrates), drop = F]
  if(ncol(rxnSubstrates) == 0){
    next # no substrates
  }
  
  rxnEnzymes <- t(a_reactionData$enzymeAbund)
  
  MEcombos <- expand.grid(1:ncol(rxnSubstrates), 1:ncol(rxnEnzymes))
  for(k in 1:nrow(MEcombos)){
    
    # Correlation
    
    J = log(relFlux)
    M = rxnSubstrates[,MEcombos[k,1]]
    E = rxnEnzymes[,MEcombos[k,2]]
    
    MEcor <- cor(J, M + E)
    
    MEcor_bs <- sapply(1:nbs_perm, function(x){
      rs <- sample(1:length(J), replace = T)
      cor(J[rs], M[rs] + E[rs])
    })
    
    # two-tailed test for correlation != 0 with a pseudo-count to avoid p = 0
    cor_sig <- min(1 - 2*abs(0.5 - sum(MEcor_bs > 0) / nbs_perm) + 1/nbs_perm, 1)
    
    # Proportionality
    
    varEx <- (sum((J - mean(J))^2) - sum((J - c(M+E) - mean(J - c(M+E)))^2))/sum((J - mean(J))^2)  # (TSS - RSS)/TSS
    
    varEx_perm <- sapply(1:nbs_perm, function(x){
      rs <- sample(1:length(J), replace = F)
      MEperm <- c(M+E)[rs]
      (sum((J - mean(J))^2) - sum((J - MEperm - mean(J - MEperm))^2))/sum((J - mean(J))^2)  # (TSS - RSS)/TSS
    })
    
    varEx_sig <- min((sum(varEx < varEx_perm) + 1)/nbs_perm, 1)
    
    comp_summary <- data.frame(
      rxNum = ifelse(is.null(ME_assoc), 1, nrow(ME_assoc) + 1),
      rx = a_reactionData$reaction, 
      metabolite = a_reactionData$metNames[names(a_reactionData$metNames) == colnames(rxnSubstrates)[MEcombos[k,1]]],
      enzyme_systematic = colnames(rxnEnzymes)[MEcombos[k,2]],
      enzyme_common = c2o$gene_name[chmatch(colnames(rxnEnzymes)[MEcombos[k,2]], c2o$systematic_name)],
      corr = MEcor, corr_pvalue = cor_sig,
      varEx = varEx, varex_pvalue = varEx_sig)
    
    ME_assoc <- rbind(ME_assoc, comp_summary)
    
    ME_plot[[length(ME_plot)+1]] <- JME_plot(J,M,E,comp_summary)
    
  }
  
  print(paste(i, "reactions done"))
}

ME_assoc$corr_qvalue <- qvalue(ME_assoc$corr_pvalue)$q
ME_assoc$varex_qvalue <- qvalue(ME_assoc$varex_pvalue)$q

ME_assoc$corr_color <- ifelse(ME_assoc$corr_qvalue > 0.1, "darkgray", NA)
ME_assoc$corr_color[is.na(ME_assoc$corr_color)] <- ifelse(ME_assoc$corr[is.na(ME_assoc$corr_color)] > 0, "darkgoldenrod1", "cornflowerblue")

ME_assoc$varex_color <- ifelse(ME_assoc$varex_qvalue > 0.1, "darkgray", "darkgoldenrod1")

# Correlation between M+E and Flux
ggplot() + geom_point(data = ME_assoc, aes(y = -1*log(corr_pvalue, base = 10), x = corr, color = corr_color, alpha = 0.7), size = 4) + geom_hline(yintercept = c(0, 2, 4, 6)) + geom_vline(xintercept = seq(-1,1,by = 0.5)) +
  geom_hline(y = -1*log(max(ME_assoc$corr_pvalue[ME_assoc$corr_qvalue < 0.1]), base = 10), color = "RED", size = 2) +
  scale_y_continuous(expression(-log[10]~pvalue), expand = c(0,0), breaks = c(0, 2, 4, 6), limits = c(0,6)) + scale_x_continuous("Pearson Correlation", limits = c(-1, 1), expand = c(0,0)) +
  scale_alpha_identity() + scale_color_identity() + scale_size_identity() +
  scatterPlotTheme
ggsave("MEcorr.pdf", height = 8, width = 8)

# Generate summary plots for the N genes with the greatest absolute correlation - only take the greatest pair for each gene
best_corr <- sapply(unique(ME_assoc$rx), function(x){
  ME_assoc$rxNum[ME_assoc$rx == x][which.max(abs(ME_assoc$corr[ME_assoc$rx == x]))]
})

MEcorr_dataTable <- data.table(ME_assoc[best_corr, c('rxNum', 'rx', 'metabolite', 'enzyme_systematic', 'corr', 'corr_pvalue', 'corr_qvalue')])
MEcorr_dataTable <- MEcorr_dataTable[corr_qvalue < 0.1,]
MEcorr_dataTable[, absCorr := abs(corr) ,]

MEcorr_dataTable <- MEcorr_dataTable[order(-absCorr)]

pdf(file = "MEcorr_examples.pdf", height = 20, width = 15)
do.call(grid.arrange,  ME_plot[MEcorr_dataTable$rxNum[1:n_example_plots]])
dev.off()


# Fraction of variance explained given proportionality
#ME_assoc$varEx[ME_assoc$varEx < 0] <- 0 # For associations that predict worse than the mean, floor them to this value

ggplot() + geom_hline(yintercept = seq(0,6, by = 1)) + geom_vline(xintercept = seq(0,1,by = 0.25)) +
  geom_point(data = ME_assoc, aes(y = -1*log(varex_pvalue, base = 10), x = varEx), color = "black", size = 6, alpha = 1) +
  geom_point(data = ME_assoc, aes(y = -1*log(varex_pvalue, base = 10), x = varEx, color = varex_color), size = 5, alpha = 1) +
  geom_hline(y = -1*log(max(ME_assoc$varex_pvalue[ME_assoc$varex_qvalue < 0.1]), base = 10), color = "RED", size = 2) +
  scale_y_continuous(expression(-log[10]~pvalue), expand = c(0.01,0), breaks = c(0, 2, 4, 6), limits = c(0,6)) + scale_x_continuous("Variance Explained", limits = c(0, 1), expand = c(0.01,0), labels = percent) +
  scale_alpha_identity() + scale_color_identity() + scale_size_identity() +
  scatterPlotTheme 
ggsave("MEvarex.pdf", height = 12, width = 12)


best_varex <- sapply(unique(ME_assoc$rx), function(x){
  ME_assoc$rxNum[ME_assoc$rx == x][which.max(ME_assoc$varEx[ME_assoc$rx == x])]
})

MEvarex_dataTable <- data.table(ME_assoc[best_varex, c('rxNum', 'rx', 'metabolite', 'enzyme_systematic', 'varEx', 'varex_pvalue', 'varex_qvalue')])
MEvarex_dataTable <- MEvarex_dataTable[varex_qvalue < 0.1,]

MEvarex_dataTable <- MEvarex_dataTable[order(-varEx)]

pdf(file = "MEvarex_examples.pdf", height = 20, width = 15)
do.call(grid.arrange,  ME_plot[MEvarex_dataTable$rxNum[1:n_example_plots]])
dev.off()

# write a summary table

ME_summary_table <- ME_assoc[,c('rx', 'metabolite', 'enzyme_common', 'enzyme_systematic', 'corr', 'corr_pvalue', 'corr_qvalue', 'varEx', 'varex_pvalue', 'varex_qvalue')]
ME_summary_table$corr_qvalue <- round(ME_summary_table$corr_qvalue, 3)
ME_summary_table$varex_qvalue <- round(ME_summary_table$varex_qvalue, 3)

write.table(ME_summary_table, file = "FluxVersusME.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@###@
##### Determine which metabolite or protein is most correlated with pathway flux for each pathway #####
##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@###@

JMorE_plot <- function(PWname, PWflux, predictorName, predictorConc){
  
  # scale predictor by median and express flux as fraction of maximal activity
  JMorEdf <- data.frame(J = PWflux / max(abs(PWflux)), predictor = predictorConc / median(predictorConc))
  
  # aesthetics
  
  sampleInfo <- data.frame(condition = rownames(JMorEdf), limitation = substr(rownames(JMorEdf), 1, 1), DR = substr(rownames(JMorEdf), 2, 5))
  sampleInfo$limitation <- factor(sampleInfo$limitation, levels = unique(sampleInfo$limitation))
  sampleInfo$size <- sqrt(as.numeric(sampleInfo$DR)*200)
  
  JMorEdf <- cbind(JMorEdf, sampleInfo)
  
  dependence_regression <- lm(data = JMorEdf, J ~ predictor)$coef
  
  ggplot(JMorEdf, aes(x = predictor, y = J, color = limitation, size = size)) + geom_point() + geom_abline(intercept = dependence_regression[1], slope = dependence_regression[2], color = "darkgoldenrod3", size = 2) +
    geom_hline(y = 0, size = 1) + geom_vline(x = 0, size = 1) + 
    scatterPlotTheme + scale_color_brewer(palette = "Set2") + scale_size_identity() + expand_limits(x = 0, y = 0) +
    scale_y_continuous(nameReformat(PWname, 30), expand = c(0.01,0), label = percent) + scale_x_continuous(predictorName, expand = c(0.01,0))
  
}

# topology: metabolism_stoichiometry - set of all directed reactions

load("../flux_cache/yeast_stoi_directed.Rdata") # S: metabolism stoichiometry
load("../flux_cache/reconstructionWithCustom.Rdata") # metabolic reconstruction files

carried_flux <- read.table("../flux_cache/fluxCarriedSimple.tsv", header = T, sep = "\t")
S_carried <- S; colnames(S_carried) <- sub('_[FR]$', '', colnames(S_carried))
S_carried <- S_carried[,colnames(S_carried) %in% rownames(carried_flux)] # reactions which carry flux
S_carried <- S_carried[rowSums(S_carried != 0) != 0,] # metabolites which are utilized
S_carried <- S_carried[!(rownames(S_carried) %in% corrFile$SpeciesID[grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', corrFile$SpeciesName)]),] # remove cofactors and common species
#grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', corrFile$SpeciesName, value = T)


# convert stoichiometric matrix to a directed bipartite graph
S_graph <- melt(t(S_carried), varnames = c("source", "sink"))
S_graph <- S_graph[S_graph$value != 0,] # reactions - metabolite links
S_graph$source <- as.character(S_graph$source); S_graph$sink <- as.character(S_graph$sink) # class coercion

S_graph[S_graph$value < 0,] <- S_graph[S_graph$value < 0,][,c(2,1,3)] # for consumed metabolites, invert direction
S_graph <- S_graph[,-3]

S_igraph <- graph.data.frame(S_graph)
V(S_igraph)$type <- V(S_igraph) %in% colnames(S_carried)

sort(betweenness(S_igraph), decreasing = T)[1:4] # remove water, H+, ATP, ADP

# define a pathway operationally as a set of reactions satisfying:
# connected through other pathway members
# proportional flux through all members (pearson correlation = 1)

# define clusters of highly correlated reactions

fluxCorrelation <- abs(cor(t(carried_flux[rownames(carried_flux) %in% colnames(S_carried),])))
reaction_adjacency <- graph.adjacency(fluxCorrelation > 0.999)
reaction_clusters <- igraph::clusters(reaction_adjacency, mode = "weak") # find clusters of correlaated reacrtions

# identify connected subgraphs using a subset of reactions defined by these correlated reaction clusters (without restricting metabolities)

pathway_set <- list()

for(clust_n in 1:reaction_clusters$no){
  
  if(length(V(reaction_adjacency)$name[reaction_clusters$membership == clust_n]) < 3){next} # skip minute clusters
  
  # reduce stoichiometric matrix to the relevent set of reactions
  pw_subset <- S_carried[,colnames(S_carried) %in% V(reaction_adjacency)$name[reaction_clusters$membership == clust_n]]
  pw_subset <- pw_subset[,!duplicated(colnames(pw_subset))]
  pw_subset <- pw_subset[rowSums(pw_subset != 0) != 0,]
  
  if(nrow(pw_subset) < 4){next}
  
  # determine the weakly connected subgraphs within this set of reactions
  S_graph_subset <- melt(t(pw_subset), varnames = c("source", "sink"))
  S_graph_subset <- S_graph_subset[S_graph_subset$value != 0,] # reactions - metabolite links
  S_graph_subset$source <- as.character(S_graph_subset$source); S_graph_subset$sink <- as.character(S_graph_subset$sink) # class coercion

  S_graph_subset[S_graph_subset$value < 0,] <- S_graph_subset[S_graph_subset$value < 0,][,c(2,1,3)] # for consumed metabolites, invert direction
  S_graph_subset <- S_graph_subset[,-3]
  
  S_igraph_subset <- graph.data.frame(S_graph_subset)
  V(S_igraph_subset)$type <- V(S_igraph_subset) %in% colnames(pw_subset)

  reaction_clusters_subset <- igraph::clusters(S_igraph_subset, mode = "weak")
  # save clusters which involve more than 3 reactions
  cluster_rxns <- lapply(1:reaction_clusters_subset$no, function(x){
    clust_members <- V(S_igraph_subset)$name[reaction_clusters_subset$membership == x]
    clust_members[clust_members %in% colnames(pw_subset)]
    })
  
  for(a_clust in 1:reaction_clusters_subset$no){
    rxns <- cluster_rxns[[a_clust]]
    
    if(length(rxns) >= 3){
      pw_output <- pw_subset[,colnames(pw_subset) %in% rxns]
      pw_output <- pw_output[rowSums(pw_output != 0) != 0,]
      
      pathway_set <- append(pathway_set, list(stoi = pw_output))
      
      }
    }
  }



# load metabolomics and proteomics

rxn_pathways <- read.delim("../flux_cache/reactionPathways.tsv")
metabolite_abundance <- read.delim('../flux_cache/tab_boer_log2rel.txt')

enzyme_abund <- read.delim("../companionFiles/proteinAbundance.tsv")
rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
rxn_enzyme_groups <- read.delim("../flux_cache/rxn_enzyme_groups.tsv")
  
#nbs <- 10000
nbs <- 10^6 # ~24 hours

rxn_pathways_carryFlux <- rxn_pathways[rxn_pathways$reactionID %in% rownames(carried_flux),]
rxn_pathways_class_list <- sapply(rxn_pathways_carryFlux$pathway, function(x){strsplit(x, split = '__')[[1]]})
names(rxn_pathways_class_list) <- rxn_pathways_carryFlux$reactionID
pw_counts <- table(unlist(rxn_pathways_class_list))


for(a_pathway in 1:length(pathway_set)){
  stoi_subset <- pathway_set[[a_pathway]]
  
  pw_flux <- apply(carried_flux[rownames(carried_flux) %in% colnames(stoi_subset),], 2, median)
  pw_flux <- pw_flux * ifelse(median(pw_flux) < 0, -1, 1) # set positive direction as direction carried in the majority of conditions
  
  # determine the name of this pathway using hypergeometric
  
  pw_annotations <- unname(unlist(rxn_pathways_class_list[colnames(stoi_subset)]))
  if(is.null(pw_annotations)){
    pw <- "Unspecified"
  }else{
    pw_annot <- table(pw_annotations)
    pw_p <- rep(NA, length(pw_annot))
    for(i in 1:length(pw_annot)){
      pw_p[i] <- 1 - phyper(q = pw_annot[i], m = pw_counts[names(pw_counts) == names(pw_annot[i])], n = sum(pw_counts) - pw_annot[i], k = sum(pw_annot))
    }
    pw <- names(pw_annot[pw_p == min(pw_p)][which.max(pw_annot[pw_p == min(pw_p)])])
  }
  pw <- paste0(pw, " (", ncol(stoi_subset), ")") 
  
  
  rx_names <- unique(rxnFile$Reaction[rxnFile$ReactionID %in% colnames(stoi_subset)])  
  
  # metabolites
  matching_met_index <- sapply(unique(corrFile$SpeciesType[corrFile$SpeciesID %in% rownames(stoi_subset)]), function(a_met){grep(a_met, metabolite_abundance$tID)})
  matching_met_df <- metabolite_abundance[unique(unlist(matching_met_index)),, drop = F]
  
  metabolite_figs <- list()
  
  if(nrow(matching_met_df) > 0){
    matching_met_matrix <- 2^matching_met_df[,grep('[PCNLU][0-9.]{4}', colnames(matching_met_df))]
    if(!all(names(pw_flux) %in% colnames(matching_met_matrix))){stop("check ordering")}
    matching_met_matrix <- matching_met_matrix[,chmatch(names(pw_flux), colnames(matching_met_matrix)), drop = F]
    
    pathway_corr <- t(sapply(1:nrow(matching_met_matrix), function(a_row){
      # regress metabolite/enzyme on flux 
      # determine a null distribution using a pivotal t-statistic-based bootstrap
      
      pw_reg <- lm(pw_flux ~ matching_met_matrix[a_row,])
      reg_resids <- pw_reg$resid * sqrt(length(pw_flux)/pw_reg$df.residual)
       
      X <- as.matrix(cbind(1, matching_met_matrix[a_row,]))
      
      pivotal_t <- t(sapply(1:nbs, function(z){ 
        pw_null <- I(pw_reg$fitted + sample(reg_resids, replace = T))
        B = solve(t(X)%*%X)%*%t(X)%*%pw_null # solve for regression coefficients
        c(B[2],
        sqrt(sum((pw_null - X%*%B)^2)/(length(pw_flux) - 2) * solve(t(X)%*%X)[2,2])) # solve for se of slope
      }))
      
      pivotal_dist <- pivotal_t[,1]/pivotal_t[,2]
      
      data.frame(pwnum = a_pathway, class = "metabolite", specie = matching_met_df$Metabolite[a_row], corr = cor(pw_flux, matching_met_matrix[a_row,]), pval = (1  - abs(0.5 - sum(pivotal_dist > 0)/nbs)*2) + (1/nbs))
      
    }))
    
    # pairwise-comparisons of M and P
    
    for(a_met in 1:nrow(matching_met_matrix)){
      
      metabolite_figs[[matching_met_df$Metabolite[a_met]]] <- JMorE_plot(PWname = pw, PWflux = pw_flux,
                                                                         predictorName = matching_met_df$Metabolite[a_met],
                                                                         predictorConc = matching_met_matrix[a_met,])
      
    }
    
  }else{
    pathway_corr <- NULL
  }
  
  # enzymes
  enzyme_subset <- 2^enzyme_abund[rownames(enzyme_abund) %in% unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% colnames(stoi_subset)]),,drop = F]
  colnames(enzyme_subset) <- toupper(colnames(enzyme_subset))
  
  enzyme_figs <- list()
  
  if(nrow(enzyme_subset) > 0){
    if(!all(names(pw_flux) %in% colnames(enzyme_subset))){stop("check naming")}
    matching_enzyme_matrix <- enzyme_subset[,chmatch(names(pw_flux), colnames(enzyme_subset)), drop = F]
    
    enzyme_corr <- t(sapply(1:nrow(matching_enzyme_matrix), function(a_row){
      # regress metabolite/enzyme on flux 
      # determine a null distribution using a pivotal t-statistic-based bootstrap
      
      pw_reg <- lm(pw_flux ~ matching_enzyme_matrix[a_row,])
      reg_resids <- pw_reg$resid * sqrt(length(pw_flux)/pw_reg$df.residual)
      
      X <- as.matrix(cbind(1, matching_enzyme_matrix[a_row,]))
      
      pivotal_t <- t(sapply(1:nbs, function(z){ 
        pw_null <- I(pw_reg$fitted + sample(reg_resids, replace = T))
        B = solve(t(X)%*%X)%*%t(X)%*%pw_null # solve for regression coefficients
        c(B[2],
        sqrt(sum((pw_null - X%*%B)^2)/(ncol(matching_enzyme_matrix) - 2) * solve(t(X)%*%X)[2,2])) # solve for se of slope
      }))
      
      pivotal_dist <- pivotal_t[,1]/pivotal_t[,2]
      
      data.frame(pwnum = a_pathway, class = "enzyme", specie = rownames(matching_enzyme_matrix)[a_row], corr = cor(pw_flux, matching_enzyme_matrix[a_row,]), pval = (1  - abs(0.5 - sum(pivotal_dist > 0)/nbs)*2) + (1/nbs))
      
    }))
    
    pathway_corr <- rbind(pathway_corr, enzyme_corr)
    
    # pairwise-comparisons of M and P
    
    for(a_enz in 1:nrow(matching_enzyme_matrix)){
      
      enzyme_figs[[rownames(enzyme_subset)[a_enz]]] <- JMorE_plot(PWname = pw, PWflux = pw_flux,
                                                                  predictorName = c2o$gene_name[chmatch(rownames(matching_enzyme_matrix)[a_enz], c2o$systematic_name)],
                                                                  predictorConc = matching_enzyme_matrix[a_enz,])
      
    }
    
    
  }
  pathway_set[[a_pathway]] <- list(pw_name = pw, rx_names = rx_names, stoi = pathway_set[[a_pathway]], corr = pathway_corr, flux = pw_flux, metabolite_figures = metabolite_figs, enzyme_figures = enzyme_figs)
  
  print(paste(a_pathway, "pathways analyzed"))
  
}


pw_associations <- sapply(1:length(pathway_set), function(x){pathway_set[[x]]$corr})
pw_associations <- do.call(rbind, pw_associations)
pw_associations <- as.data.frame(apply(pw_associations, c(1,2), as.character))
pw_associations$qval <- qvalue(as.numeric(pw_associations$pval))$qvalues

pw_associations$color <- ifelse(pw_associations$qval > 0.1, "darkgray", NA)
pw_associations$color[is.na(pw_associations$color)] <- ifelse(pw_associations$corr[is.na(pw_associations$color)] > 0, "darkgoldenrod1", "cornflowerblue")

# Correlation between pathway flux and pathway metabolites or enzymes
ggplot() + geom_point(data = pw_associations, aes(y = -1*log(as.numeric(pval), base = 10), x = as.numeric(corr), color = color, alpha = 0.7, shape = class), size = 4) + geom_hline(yintercept = c(0, 2, 4, 6)) + geom_vline(xintercept = seq(-1,1,by = 0.5)) +
  geom_hline(y = -1*log(max(as.numeric(pw_associations$pval)[pw_associations$qval < 0.1]), base = 10), color = "RED", size = 2) + facet_grid(class ~ .) + 
  scale_y_continuous(expression(-log[10]~pvalue), expand = c(0.01,0), breaks = c(0, 2, 4, 6), limits = c(0,6)) + scale_x_continuous("Pearson Correlation", limits = c(-1, 1), expand = c(0,0)) +
  scale_alpha_identity() + scale_color_identity() + scale_size_identity() +
  scatterPlotTheme
ggsave("pathwayCorr.pdf", height = 14, width = 9)


# Look at examples of the most correlated (- or +) metabolites or enzymes with pathway flux

pw_associations <- data.table(pw_associations)
pw_associations[,pwnum := as.numeric(pwnum)]; pw_associations[,corr := as.numeric(corr)] ; pw_associations[,pval := as.numeric(pval)] ; pw_associations[,qval := as.numeric(qval)]


pw_associations_strongest <- pw_associations[,list(specie = specie[which.max(abs(corr))],
                                                   corr = corr[which.max(abs(corr))],
                                                   pval = pval[which.max(abs(corr))],
                                                   qval = qval[which.max(abs(corr))]
), by = c("pwnum", "class")]

  
met_subset <- pw_associations_strongest[class == "metabolite",]
met_subset <- met_subset[order(-abs(corr))]

pdf(file = "pathwayCorr_mets.pdf", height = 20, width = 15)
do.call(grid.arrange,
  lapply(1:n_example_plots, function(i){
    pathway_set[[met_subset[i,pwnum]]]$metabolite_figures[[met_subset[i,specie]]]
  })
)
dev.off()
  
enz_subset <- pw_associations_strongest[class == "enzyme",]
enz_subset <- enz_subset[order(-abs(corr))]

pdf(file = "pathwayCorr_enzs.pdf", height = 20, width = 15)
do.call(grid.arrange,
  lapply(1:n_example_plots, function(i){
    pathway_set[[enz_subset[i,pwnum]]]$enzyme_figures[[enz_subset[i,specie]]]
  })
)
dev.off()

# write pathway associations to a file
pw_corr_output <- pw_associations[,c('pwnum', 'class', 'specie', 'corr', 'pval', 'qval'), with = F]

pw_corr_output[, pval := round(pval, 6) ]
pw_corr_output[, qval := round(qval, 6) ]

pw_corr_output[, pwName := pathway_set[[pwnum]]$pw_name, by = "pwnum"]
pw_corr_output[, pwReactions := paste(pathway_set[[pwnum]]$rx_names, collapse = "; "), by = "pwnum"]

write.table(pw_corr_output, file = "PW_corr_table.tsv", sep = "\t", row.names = F, col.names = T, quote = F)










library(colorRamps)
heatmap.2(abs(cor(t(carried_flux))), trace = "none")
heatmap.2((fluxCorrelation > 0.999)*1, trace = "none")
table(abs(cor(t(carried_flux))) > 0.999)



# igraph tcl-tk interface ?
tkplot(S_igraph)
tkplot(S_igraph)
layout.norm # rescale = F
# load
# add species
# save
# subset for relevant species