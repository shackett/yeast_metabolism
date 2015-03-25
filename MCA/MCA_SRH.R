######## Study control properties of sets of continuous sets of reactions using MCA #########

setwd("~/Desktop/Rabinowitz/FBA_SRH/MCA")
load("../Yeast_genome_scale/flux_cache/paramCI.Rdata") # load reaction summaries 
load("../Yeast_genome_scale/flux_cache/elasticityData.Rdata")

###
library(data.table)
library(reshape2)
library(ggplot2)
library(grid)

options(stringsAsFactors = F)

### For now, manually define pathways to analyze

pathway_models <-  c("Glycolysis with FBP feed-forward", # glycolysis without FBP feed-forward
                     "Glycolysis", # glycolysis with FBP feed-forward
                     "Arginine Biosynthesis", # Arginine biosynthesis
                     "Cysteine Biosynthesis" # Cysteine biosynthesis
)


for(a_model in pathway_models){
  
  if(a_model == "Arginine Biosynthesis"){
    
    modelRxns <- c(OTCase = "r_0816-rm-t_0461-inh-noncomp", # OTCase inhibited by alanine
                   ASS = "r_0208-rm", # arginosuccinate synthase
                   ASL = "r_0207-rm" # arginosuccinate lyase
    )
    
    modelMet <- c(t_0471 = "L-citrulline", t_0010 = "(N(omega)-L\narginino)succinic acid")
    
    modelRxnNames <- c(OTCase = 'Alanine -| OTCase') # optional renaming
    equivalent_mets <- NULL
  }
  
  if(a_model == "Cysteine Biosynthesis"){
    
    modelRxns <- c(CBS = "r_0309-rm", # cystathionine beta synthase
                   CBL = "r_0310-rm-t_0461-inh-noncomp" # cystathionine beta lyase
    )
    modelMet <- c(t_0472 = "L-cystathionine")
    
    modelRxnNames <- c(CBL = 'Alanine -| CBL') # optional renaming
    equivalent_mets <- NULL
  }
  
  if(a_model %in% c("Glycolysis", "Glycolysis with FBP feed-forward")){
    
    modelRxns <- c(PFK = "r_0886-rm-t_0452-inh-uncomp", # PFK
                   ALD = "r_0450-rm-t_0234-inh-noncomp", # ALD
                   GAPDH = "r_0486-rm", # GAPDH
                   PGK = "r_0892-rm", # PGK
                   #PGM = "r_0893-rm-t_0139-inh-noncomp", # PGM
                   PyK = "r_0962-rm-t_0276-inh-noncomp" # PyK
    )
    
    modelRxnNames <- c(PFK = 'Isocitrate -| PFK',
                       ALD = 'AMP -| ALD',
                       PGM = '3PG -| PGM',
                       PyK = 'Citrate -| PyK') # optional renaming
    
    if(a_model == "Glycolysis with FBP feed-forward"){
      modelRxns[names(modelRxns) == "PyK"] <- "r_0962-rm-t_0290-act-mm"
      modelRxnNames['PyK'] <- 'FBP -+ PyK'
      
    }
    modelMet <- c(t_0290 = "D-fructose 1,6\nbisphosphate", t_0386 = "glyceraldehyde 3\nphosphate", t_0039 = "1,3-bisphospho-D\nglycerate",
                  t_0139 = "3-phosphoglycerate")
    
    equivalent_mets <- data.frame(subsumes = c("glyceraldehyde 3\nphosphate", "3-phosphoglycerate"), consumed = c("dihydroxyacetone\nphosphate", "phosphoenolpyruvate"), descrip = c("DHAP -> GAP", "PEP -> 3PG"))
      
  }
  
  # Find distributions of elasticities w.r.t. each pathway metabolite for each reaction
  
  rxn_elasticities <- lapply(unname(modelRxns), function(x){
    rxn_elasticity <- ELdata[[x]]
    rxn_elasticity <- rxn_elasticity[rxn_elasticity$specie %in% union(equivalent_mets$consumed, modelMet),]
    rxn_elasticity$specie <- as.character(rxn_elasticity$specie)
    
    # If both equivalent metabolites are in a reaction then average their elasticities
    # If only a consumed metabolite is present then just change its name
    if(!is.null(equivalent_mets)){
      paired_equivalent_mets <- equivalent_mets$subsumes %in% unique(rxn_elasticity$specie) & equivalent_mets$consumed %in% unique(rxn_elasticity$specie)
      
      if(any(paired_equivalent_mets)){
        for(i in c(1:nrow(equivalent_mets))[paired_equivalent_mets]){
          equiv_met_elasticity <- data.table(rxn_elasticity[rxn_elasticity$specie %in% c(equivalent_mets$subsumes[i], equivalent_mets$consumed[i]),])
          equiv_met_elasticity <- equiv_met_elasticity[,list(specie = equivalent_mets$subsumes[i], Elasticity = mean(Elasticity)), by = c("condition", "markovSample")]
          setcolorder(equiv_met_elasticity, colnames(rxn_elasticity))
          
          rxn_elasticity <- rxn_elasticity[!(rxn_elasticity$specie %in% c(equivalent_mets$subsumes[i], equivalent_mets$consumed[i])),]
          rxn_elasticity <- rbind(rxn_elasticity, as.data.frame(equiv_met_elasticity))
        }
      }
      
      # Rename un-paired consumed metabolites
      if(any(rxn_elasticity$specie %in% equivalent_mets$consumed)){
      rxn_elasticity$specie[rxn_elasticity$specie %in% equivalent_mets$consumed] <- equivalent_mets$subsumes[chmatch(as.character(rxn_elasticity$specie[rxn_elasticity$specie %in% equivalent_mets$consumed]), equivalent_mets$consumed)]
      }
    }
      
    if(nrow(rxn_elasticity) == 0){
      warning("No metabolite match for reaction ", x)
      NULL
    }else{
      data.frame(reaction = x, rxn_elasticity[rxn_elasticity$specie %in% modelMet,])
    }
  })
  rxn_elasticities <- do.call("rbind", rxn_elasticities)
  
  missing_mets <- modelMet[!(modelMet %in% unique(rxn_elasticities$specie))]
  if(length(missing_mets) != 0){
    for(i in missing_mets){
      warning("No matches found for ", i, ": please check name")
    }
  }
  
  # Add additional element to setup the summation rule
  MCAsummation <- data.frame(specie = "summation", unique(rxn_elasticities[,c('condition', 'reaction', 'markovSample')]), Elasticity = 1)
  rxn_elasticities <- rbind(MCAsummation, rxn_elasticities)
  
  MCAmat <- acast(rxn_elasticities, specie ~ reaction ~ markovSample ~ condition, value.var = "Elasticity", fill = 0)
  # MCAmat[,,1,1]
  # For each condition & markov sample calculate control coefficients
  controlCoef <- apply(MCAmat, c(3,4), function(x){solve(x)[,'summation']})
  controlCoef_melt <- melt(controlCoef, varnames =  c("Reaction", "MarkovSample", "Condition"), value.name = "controlCoefficient")
  
  # Beautification - spruce up the naming
  
  # Reactions:
  rxnLabels <- data.frame(RXN = names(modelRxns), RXN_ID = unname(modelRxns), RXN_RENAME = names(modelRxns))
  rxnLabels$RXN_RENAME[rxnLabels$RXN %in% names(modelRxnNames)] <- unname(modelRxnNames)[chmatch(rxnLabels$RXN[rxnLabels$RXN %in% names(modelRxnNames)], names(modelRxnNames))]
  controlCoef_melt$rxnName <- factor(rxnLabels$RXN_RENAME[chmatch(as.character(controlCoef_melt$Reaction), rxnLabels$RXN_ID)], levels = rxnLabels$RXN_RENAME)
  
  # Conditions
  controlCoef_melt$Condition = factor(controlCoef_melt$Condition, levels = paste0(rep(c("P", "C", "N", "L", "U"), each = 5), rep(c("0.05", "0.11", "0.16", "0.22", "0.30"), times = 5)))
  controlCoef_melt$Limitation <- factor(substr(controlCoef_melt$Condition, 1, 1), levels = c("P", "C", "N", "L", "U"))
  
  
  boxplot_theme <- theme(text = element_text(size = 25, face = "bold"), title = element_text(size = 25, face = "bold"), panel.background = element_rect(fill = "mintcream"), axis.title.y = element_text(vjust = 0.35), 
                         legend.position = "top", panel.grid = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(size = 1, color = "chartreuse4"),
                         axis.text = element_text(color = "black", hjust = 1), axis.text.x = element_text(size = 20, angle = 90, hjust = 1), axis.line = element_line(size = 2, color = "chartreuse4"), axis.line.x = element_blank(), 
                         panel.margin = unit(1.5, "lines"), strip.background = element_rect(fill = "darkseagreen2"), strip.text = element_text(size = 25, vjust = 0.6, colour = "darkblue"))
  
  ggplot(controlCoef_melt, aes(x = Condition, y = controlCoefficient, fill = Limitation)) + geom_hline(yintercept = -0.02, size = 4, color = "chartreuse4") +
    #geom_jitter(position = position_jitter(width = .4), size = 0.3, alpha = 0.5) +
    geom_boxplot(outlier.size = 0, alpha = 0.7) + facet_grid(rxnName ~ .) +
    boxplot_theme + scale_fill_brewer(palette = "Set2", guide = F) +
    scale_x_discrete("") + scale_y_continuous("Pathway Control Coefficient", expand = c(0,0))
  
  ggsave(file = paste0("MCA_plots/", a_model, ".pdf"), width = 10, height = nrow(rxnLabels) * 4)
  
}