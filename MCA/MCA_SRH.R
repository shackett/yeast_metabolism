######## Study control properties of sets of continuous sets of reactions using MCA #########

setwd("~/Desktop/Rabinowitz/FBA_SRH/MCA")
load("../Yeast_genome_scale/flux_cache/paramCI.Rdata") # load reaction summaries 
load("../Yeast_genome_scale/flux_cache/elasticityData.Rdata")

###
library(data.table)
library(reshape2)
library(ggplot2)

options(stringsAsFactors = F)

### For now, manually define pathways to analyze

pathway_models <-  c("GlycFBP", # glycolysis without FBP feed-forward
                     "GlycSTD", # glycolysis with FBP feed-forward
                     "ArgBio", # Arginine biosynthesis
                     "CysBio" # Cysteine biosynthesis
)


for(a_model in pathway_models){
  
  if(a_model == "ArgBio"){
    
    modelRxns <- c(OTCase = "r_0816-rm-t_0461-inh-noncomp", # OTCase inhibited by alanine
                   ASS = "r_0208-rm", # arginosuccinate synthase
                   ASL = "r_0207-rm" # arginosuccinate lyase
    )
    
    modelMet <- c(t_0471 = "L-citrulline", t_0010 = "(N(omega)-L\narginino)succinic acid")
    
    modelRxnNames <- c(OTCase = 'Alanine -| OTCase') # optional renaming
    equivalent_mets <- NULL
  }
  
  if(a_model == "CysBio"){
    
    modelRxns <- c(CBS = "r_0309-rm", # cystathionine beta synthase
                   CBL = "r_0310-rm-t_0461-inh-noncomp" # cystathionine beta lyase
    )
    modelMet <- c(t_0472 = "L-cystathionine")
    
    modelRxnNames <- c(CBL = 'Alanine -| CBL') # optional renaming
    equivalent_mets <- NULL
  }
  
  if(a_model %in% c("GlycSTD", "GlycFBP")){
    
    modelRxns <- c(PFK = "r_0886-rm-t_0452-inh-uncomp", # PFK
                   ALD = "r_0450-rm-t_0234-inh-noncomp", # ALD
                   GAPDH = "r_0486-rm", # GAPDH
                   PGK = "r_0892-rm", # PGK
                   PGM = "r_0893-rm-t_0139-inh-noncomp", # PGM
                   PyK = "r_0962-rm-t_0276-inh-noncomp" # PyK
    )
    
    modelRxnNames <- c(PFK = 'Isocitrate -| PFK',
                       ALD = 'AMP -| ALD',
                       PGM = '3PG -| PGM',
                       PyK = 'Citrate -| PyK') # optional renaming
    
    if(a_model == "GlycFBP"){
      modelRxns[names(modelRxns) == "PyK"] <- "r_0962-rm-t_0290-act-mm"
      modelRxnNames['PyK'] <- 'FBP -+ PyK'
      
    }
    modelMet <- c(t_0290 = "D-fructose 1,6-bisphosphate", t_0386 = "glyceraldehyde 3-phosphate", t_0039 = "1,3-bisphospho-D-glycerate",
                  t_0139 = "3-phosphoglycerate", t_0617 = "phosphoenolpyruvate")
    
    equivalent_mets <- data.frame(subsumes = c("t_0386", "t_0617"), consumed = c("t_0331", "t_0096"), descrip = c("DHAP -> GAP", "2PG -> PEP"))
      
  }
  
  # Find distributions of elasticities w.r.t. each pathway metabolite for each reaction
  
  rxn_elasticities <- lapply(unname(modelRxns), function(x){
    rxn_elasticity <- ELdata[[x]]
    rxn_elasticity <- rxn_elasticity[rxn_elasticity$specie %in% modelMet,]
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
  rxn_elasticities <- rbind(rxn_elasticities, MCAsummation)
  
  MCAmat <- acast(rxn_elasticities, specie ~ reaction ~ markovSample ~ condition, value.var = "Elasticity", fill = 0)
  # For each condition & markov sample calculate control coefficients
  controlCoef <- apply(MCAmat, c(3,4), function(x){solve(x)[,'summation']})
  controlCoef_melt <- melt(controlCoef, varnames =  c("Reaction", "MarkovSample", "Condition"), value.name = "controlCoefficient")
  
  # Beautification - spruce up the naming
  
  controlCoef_melt$Condition = factor(controlCoef_melt$Condition, levels = paste0(rep(c("P", "C", "N", "L", "U"), each = 5), rep(c("0.05", "0.11", "0.16", "0.22", "0.30"), times = 5)))
  controlCoef_melt$Limitation <- factor(substr(controlCoef_melt$Condition, 1, 1), levels = c("P", "C", "N", "L", "U"))
  
  cleanName <- names(modelRxns)[chmatch(as.character(controlCoef_melt$Reaction), modelRxns)]
  cleanName[cleanName %in% names(modelRxnNames)] <- unname(modelRxnNames)[chmatch(cleanName[cleanName %in% names(modelRxnNames)], names(modelRxnNames))]
  controlCoef_melt$rxnName <- cleanName
  
  ggplot(controlCoef_melt, aes(x = Condition, y = controlCoefficient, fill = Limitation)) + geom_boxplot(outlier.size = 0) + facet_grid(rxnName ~ .) +
  geom_jitter(position = position_jitter(width = .2), size = 0.3)
  
}