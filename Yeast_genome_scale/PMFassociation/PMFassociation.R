#### Assessing the extent to which flux can be simply explained through either :
# A) Trade-offs between metabolite abundances and enzyme levels maintain robustness in flux
# B) Rate-limiting enzymes or flux-sensing metabolites

setwd("~/Desktop/Rabinowitz/FBA_SRH/Yeast_genome_scale/PMFassociation")

# libraries and sources
library(qvalue)
library(ggplot2)
library(grid)

scatterPlotTheme <- theme(text = element_text(size = 23, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(color = "black", fill = "white"), legend.position = "none", 
                          panel.grid.minor = element_blank(), legend.key.width = unit(6, "line"), panel.grid.major = element_line(colour = "black"), axis.ticks = element_line(colour = "black"), strip.background = element_rect(fill = "cyan"),
                          axis.text = element_text(color = "blacK")) 

JME_plot <- function(J,M,E){
  
  JMEdf <- data.frame(J = J, ME = M + E)
  JMEdf <- JMEdf - rep(1, nrow(JMEdf)) %*% t(apply(JMEdf, 2, mean)) # center flux and M*E
  
  ggplot(JMEdf, aes(x = ME, y = J)) + geom_point() + geom_abline(intercept = 0, slope = 1)
  
}

##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@
##### Trade-offs between metabolite abundances and enzyme levels maintain robustness in flux ######
##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@##@

load("../flux_cache/paramCI.Rdata")
load("../paramOptim.Rdata")

MMreactionMechs <- reactionInfo$rMech[reactionInfo$modification == ""]

MMreactionData <- rxnList_form[names(rxnList_form) %in% MMreactionMechs]

# all flux > 0 | all flux < 0
# dealing with multiple E or M : all pair-wise
nbs <- 10000
ME_TO <- NULL
i <- which(names(MMreactionData) == "r_0195-rm")
for(i in 1:length(MMreactionData)){
  
  a_reactionData <- MMreactionData[[i]]
  
  relFlux <- (a_reactionData$flux$FVAmax + a_reactionData$flux$FVAmin)/2
  
  if(!(all(relFlux > 0) | all(relFlux < 0))){
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
  
  rxnEnzymes <- t(a_reactionData$enzymeComplexes)
  
  MEcombos <- expand.grid(1:ncol(rxnSubstrates), 1:ncol(rxnEnzymes))
  for(k in 1:nrow(MEcombos)){
    
    MEcor <- cor(log(relFlux), rxnSubstrates[,MEcombos[k,1]] + rxnEnzymes[,MEcombos[k,2]])
    
    J = log(relFlux)
    M = rxnSubstrates[,MEcombos[k,1]]
    E = rxnEnzymes[,MEcombos[k,2]]
    
    MEcor_bs <- sapply(1:nbs, function(x){
      rs <- sample(1:length(J), replace = T)
      cor(J[rs], M[rs] + E[rs])
    })
    
    # two-tailed test for correlation != 0 with a pseudo-count to avoid p = 0
    cor_sig <- min(1 - 2*abs(0.5 - sum(MEcor_bs > 0) / nbs) + 1/nbs, 1)
    
    ME_TO <- rbind(ME_TO, 
                   data.frame(rx = names(MMreactionData)[i], metabolite = colnames(rxnSubstrates)[MEcombos[k,1]], enzyme = colnames(rxnEnzymes)[MEcombos[k,2]], 
                              corr = MEcor, pvalue = cor_sig)
    )
  }
}

ME_TO$qvalue <- qvalue(ME_TO$pvalue)$q
#ME_TO_sig <- ME_TO[ME_TO$qvalue < 0.05,]

ME_TO$color <- ifelse(ME_TO$qvalue > 0.1, "darkgray", NA)
ME_TO$color[is.na(ME_TO$color)] <- ifelse(ME_TO$corr[is.na(ME_TO$color)] > 0, "darkgoldenrod1", "cornflowerblue")
  
ggplot() + geom_point(data = ME_TO, aes(y = -1*log(pvalue, base = 10), x = corr, color = color, alpha = 0.7), size = 4) + geom_hline(yintercept = c(0, 2, 4)) + geom_vline(xintercept = seq(-1,1,by = 0.5)) +
  geom_hline(y = -1*log(max(ME_TO$pvalue[ME_TO$qvalue < 0.1]), base = 10), color = "RED", size = 2) +
  scale_y_continuous(expression(-log[10]~pvalue), expand = c(0,0), breaks = c(0, 2, 4), limits = c(0,4)) + scale_x_continuous("Pearson Correlation", limits = c(-1, 1), expand = c(0,0)) +
  scale_alpha_identity() + scale_color_identity() + scale_size_identity() +
  scatterPlotTheme




##### Determine which metabolite or protein is most correlated with pathway flux for each pathway #####
# topology: metabolism_stoichiometry - set of all directed reactions

library(igraph)

load("flux_cache/yeast_stoi_directed.Rdata") # S: metabolism stoichiometry
load("flux_cache/reconstructionWithCustom.Rdata") # metabolic reconstruction files

carried_flux <- read.table("Flux_analysis/fluxCarriedSimple.tsv", header = T, sep = "\t")
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
  
  if(length(V(reaction_adjacency)$name[reaction_clusters$membership == clust_n]) < 3){next} # catch minute clusters
  
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

metabolite_abundance <- read.delim('flux_cache/tab_boer_log2rel.txt')

enzyme_abund <- read.delim("./companionFiles/proteinAbundance.tsv")
rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
rxn_enzyme_groups <- read.delim("./flux_cache/rxn_enzyme_groups.tsv")
  
nbs <- 10000

 
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
  
  rx_names <- unique(rxnFile$Reaction[rxnFile$ReactionID %in% colnames(stoi_subset)])  
  
# metabolites
  matching_met_index <- sapply(unique(corrFile$SpeciesType[corrFile$SpeciesID %in% rownames(stoi_subset)]), function(a_met){grep(a_met, metabolite_abundance$tID)})
  matching_met_df <- metabolite_abundance[unique(unlist(matching_met_index)),, drop = F]
  
  if(nrow(matching_met_df) > 0){
    matching_met_matrix <- 2^matching_met_df[,grep('[PCNLU][0-9.]{4}', colnames(matching_met_df))]
    if(!all(names(pw_flux) %in% colnames(matching_met_matrix))){stop("check ordering")}
    matching_met_matrix <- matching_met_matrix[,chmatch(colnames(matching_met_matrix), names(pw_flux)), drop = F]
    
    pathway_corr <- t(sapply(1:nrow(matching_met_matrix), function(a_row){
      # regress metabolite/enzyme on flux 
      # determine a null distribution using a pivotal t-statistic-based bootstrap
      
      pw_reg <- lm(pw_flux ~ matching_met_matrix[a_row,])
      reg_resids <- pw_reg$resid * sqrt(length(pw_flux)/pw_reg$df.residual)
      
      pivotal_t <- t(sapply(1:nbs, function(z){ 
        pw_null <- I(pw_reg$fitted + sample(reg_resids, replace = T))
        summary(lm(pw_null ~ matching_met_matrix[a_row,]))$coef[2,1:2]
        }))
      pivotal_dist <- pivotal_t[,1]/pivotal_t[,2]
      
      data.frame(pwnum = a_pathway, class = "metabolite", specie = matching_met_df$Metabolite[a_row], corr = cor(pw_flux, matching_met_matrix[a_row,]), pval = (1  - abs(0.5 - sum(pivotal_dist > 0)/nbs)*2) + (1/nbs))
      
      }))
      
    #pathway_corr <- t(sapply(1:nrow(matching_met_matrix), function(a_row){
    #  null = sapply(1:nperm, function(z){cor(pw_flux, sample(matching_met_matrix[a_row,]))})
    #  data.frame(pwnum = a_pathway, class = "metabolite", specie = matching_met_df$Metabolite[a_row], corr = cor(pw_flux, matching_met_matrix[a_row,]), pval = 1 - sum(cor(pw_flux, matching_met_matrix[a_row,]) > null)/(nperm + 1))
    #  }))
    
    }else{
      pathway_corr <- NULL
    }
  
  # enzymes
  enzyme_subset <- 2^enzyme_abund[rownames(enzyme_abund) %in% unique(rxn_enzyme_groups$enzyme[rxn_enzyme_groups$reaction %in% colnames(stoi_subset)]),,drop = F]
  colnames(enzyme_subset) <- toupper(colnames(enzyme_subset))
  
  if(nrow(enzyme_subset) > 0){
    if(!all(names(pw_flux) %in% colnames(enzyme_subset))){stop("check ordering")}
    matching_enzyme_matrix <- enzyme_subset[,chmatch(colnames(enzyme_subset), names(pw_flux)), drop = F]
    
    enzyme_corr <- t(sapply(1:nrow(matching_enzyme_matrix), function(a_row){
      # regress metabolite/enzyme on flux 
      # determine a null distribution using a pivotal t-statistic-based bootstrap
      
      pw_reg <- lm(pw_flux ~ matching_enzyme_matrix[a_row,])
      reg_resids <- pw_reg$resid * sqrt(length(pw_flux)/pw_reg$df.residual)
      
      pivotal_t <- t(sapply(1:nbs, function(z){ 
        pw_null <- I(pw_reg$fitted + sample(reg_resids, replace = T))
        summary(lm(pw_null ~ matching_enzyme_matrix[a_row,]))$coef[2,1:2]
        }))
      pivotal_dist <- pivotal_t[,1]/pivotal_t[,2]
      
      data.frame(pwnum = a_pathway, class = "enzyme", specie = rownames(matching_enzyme_matrix)[a_row], corr = cor(pw_flux, matching_enzyme_matrix[a_row,]), pval = (1  - abs(0.5 - sum(pivotal_dist > 0)/nbs)*2) + (1/nbs))
      
      }))
    
    #enzyme_corr <- t(sapply(1:nrow(matching_enzyme_matrix), function(a_row){
    #  null = sapply(1:nperm, function(z){cor(pw_flux, sample(matching_enzyme_matrix[a_row,]))})
    #  data.frame(pwnum = a_pathway, class = "enzyme", specie = rownames(matching_enzyme_matrix)[a_row], corr = cor(pw_flux, matching_enzyme_matrix[a_row,]), pval = 1 - sum(cor(pw_flux, matching_enzyme_matrix[a_row,]) > null)/(nperm + 1))
    #  }))
    pathway_corr <- rbind(pathway_corr, enzyme_corr)
    }
  pathway_set[[a_pathway]] <- list(pw_name = pw, rx_names = rx_names, stoi = pathway_set[[a_pathway]], corr = pathway_corr, flux = pw_flux)
}

library(qvalue)

pw_associations <- sapply(1:length(pathway_set), function(x){pathway_set[[x]][[2]]})
pw_associations <- do.call(rbind, pw_associations)
pw_associations <- as.data.frame(apply(pw_associations, c(1,2), as.character))
pw_associations$qval <- qvalue(as.numeric(pw_associations$pval))$qvalues

# some metabolites track flux through pathways but many enzymes and metabolites are anti-correlated with flux-carried

sig_pw_associations <- pw_associations[pw_associations$qval < 0.05,]

sig_met_pw_associations <- sig_pw_associations[sig_pw_associations$class == "metabolite",]

# plot significant examples

pw_flux <- apply(carried_flux[rownames(carried_flux) %in% colnames(stoi_subset),], 2, median)
pw_flux <- pw_flux * ifelse(median(pw_flux) < 0, -1, 1) # set positive direction as direction carried in the majority of conditions










library(colorRamps)
heatmap.2(abs(cor(t(carried_flux))), trace = "none", col = blue2yellow(100))
heatmap.2((fluxCorrelation > 0.999)*1, trace = "none", col = green2red(2))
table(abs(cor(t(carried_flux))) > 0.999)



# igraph tcl-tk interface ?
tkplot(S_igraph)
tkplot(S_igraph)
layout.norm # rescale = F
# load
# add species
# save
# subset for relevant species