# Match brenda to chebiIDs by using the synonymes from the CHEBI_complete_sdf file
options(stringsAsFactors=FALSE)


if (file.exists('./flux_cache/brendamat.tsv')){
  brendamat <- read.delim('./flux_cache/brendamat.tsv')
} else {
  
  # Import all lists from Brenda, put it together to a large matrix with ligandID - name
  brendamat <- matrix(nrow = 0,ncol=2)
  
  all_activators <- read.delim('companionFiles/all_activators.tsv',sep="\t",header = TRUE)
  brendamat <- rbind(brendamat,cbind(all_activators$ligandStructureId,all_activators$activatingCompound))
  rm(all_activators)
  
  all_ki <- read.delim('companionFiles/all_ki.tsv',sep="\t",header = TRUE)
  brendamat <- rbind(brendamat,cbind(all_ki$ligandStructureId,all_ki$inhibitor))
  rm(all_ki)
  
  all_km <- read.delim('companionFiles/all_km.tsv',sep="\t",header = TRUE)
  brendamat <- rbind(brendamat,cbind(all_km$ligandStructureId,all_km$substrate))
  rm(all_km)
  
  # clean the matrix up
  brendamat <- brendamat[!is.na(brendamat[ ,1]), ]
  brendamat <- unique(brendamat)
  
  # sort the matrix
  brendamat <- data.frame(brendamat[order(as.integer(brendamat[ ,1]),brendamat[ ,2], decreasing = FALSE), ])
  colnames(brendamat) <- c('ligandID','name')
  
  
  # match chebis
  
  brendamat$chebi <- sapply(brendamat$name,MatchName2Chebi) # this could be sped-up
  
  write.table(brendamat,'./flux_cache/brendamat.tsv',sep='\t', quote = F)
}


### Match brenda and tIDs

# match chebi-chebi
# match fuzzy_chebi-fuzzy_chebi (a more degenerate match)

brendamat <- brendamat %>% tbl_df() %>% filter(!is.na(chebi)) %>% rowwise() %>%
  mutate(fuzchebi = Chebi2fuzchebi(chebi))

# join BRENDA compound IDs with model CHEBI IDs
chebi_match <- listTID %>% inner_join(brendamat, by = c("CHEBI" = "chebi"))
fuz_chebi_match <- listTID %>% inner_join(brendamat, by = c("fuzCHEBI" = "fuzchebi"))

matchTIDbrenda <- rbind(
  chebi_match %>% dplyr::select(SpeciesType, ligandID),
  fuz_chebi_match %>% dplyr::select(SpeciesType, ligandID)
) %>% unique() %>% dplyr::rename(TID = SpeciesType)

# Sanity check : name associated with TID is similar to BRENDA association
#namesMatch <- cbind(sapply(matchTIDbrenda$TID,function(x){paste(listTID$SpeciesName[ listTID$SpeciesType %in% x],collapse='; ')}),
#                    sapply(matchTIDbrenda$ligandID,function(x){paste(brendamat$name[ brendamat$ligandID %in% x],collapse='; ')}))


### Make a activator/inhibitor summary table
# make a activator/inhibitor table with columns: rxn,EC,name,ligandID,tID,chebi,kegg,measured (2=abs,1=rel,0=not),
#modtype (act/inh), subtype (all/comp), origin (brenda, textbook), organism,literature,man validated

modTable <- data.frame(rxn=character(),
                       EC=character(),
                       name=character(),
                       ligandID=numeric(),
                       tID =character(),
                       chebi =numeric(),
                       kegg = character(),
                       measured = factor(levels = c("abs","rel","not")),
                       modtype = factor(levels = c("act","inh")),
                       subtype = factor(levels = c("allosteric","competitive","uncompetitive","noncompetitive","mixed")),
                       k = numeric(),
                       kmax = numeric(),
                       origin = character(),
                       organism = character(),
                       literature = numeric(),
                       commentary = character(),
                       hill = numeric(),
                       manvalidated = factor(levels = c("yes","no")),
                       stringsAsFactors=FALSE)


all_activators <- read.delim('companionFiles/all_activators.tsv',sep="\t",header = TRUE)
all_ki <- read.delim('companionFiles/all_ki.tsv',sep="\t",header = TRUE)

all_brenda <- data.frame(EC = character(),
                         name = character(),
                         modtype = character(),
                         ligandID = character(),
                         organism = character(),
                         k = numeric(),
                         kmax = numeric(),
                         literature = character(),
                         commentary = character()
)

ECs <- data.frame(c(all_activators$ecNumber,all_ki$ecNumber))
colnames(ECs) <- 'EC'
all_brenda<- merge(ECs,all_brenda,by.x='EC',by.y='EC',all.x=T,sort=F)
rm(ECs)
all_brenda$name <- c( all_activators$activatingCompound, all_ki$inhibitor )
all_brenda$ligandID <- c( all_activators$ligandStructureId, all_ki$ligandStructureId)
all_brenda$organism <- c( all_activators$organism, all_ki$organism)
all_brenda$literature <- c( all_activators$literature, all_ki$literature)
all_brenda$commentary <- c(all_activators$commentary, all_ki$commentary)
all_brenda$k <- c(rep(NA,nrow(all_activators)), all_ki$kiValue)
all_brenda$kmax <- c(rep(NA,nrow(all_activators)), all_ki$kiValueMaximum)
all_brenda$modtype <- c(rep("act",nrow(all_activators)),rep("inh",nrow(all_ki)))

all_brenda <- all_brenda[!is.na(all_brenda$ligandID),]
all_brenda$tID <- sapply(all_brenda$ligandID, function(x){paste(matchTIDbrenda$TID[matchTIDbrenda$ligandID == x], collapse = "/")})
all_brenda <- all_brenda[all_brenda$tID != "",]
all_brenda <- data.table(unique(all_brenda))
all_brenda[,isYeast := ifelse(organism == "Saccharomyces cerevisiae", T, F),]
all_brenda <- all_brenda[is.na(k) | k > 0,,]


rm(all_activators,all_ki)

#commonOrganisms <- names(sort(table(all_brenda$organism), decreasing = T)[1:10])
#all_brenda <- all_brenda[all_brenda$organism %in% commonOrganisms,]

####### Summarize Ki and Km values for use in future comparisions #####

if(!file.exists("flux_cache/metaboliteAffinities.tsv")){
  
  modifiers_summary <- all_brenda[,list(log10mean = mean(log10(k), na.rm = T),
                                        sdOflog10 = sd(log10(k), na.rm = T),
                                        nQuant = length(k[!is.na(k)]),
                                        nQual = length(k)), by = c("EC", "tID", "modtype", "isYeast")]
  # modifiers_summary <- modifiers_summary[!(modifiers_summary$nQual == 1 & modifiers_summary$nQuant == 0),]
  
  #ggplot(melt(table(modifiers_summary$nQuant)), aes(x = Var1, y = log2(value))) + geom_point()
  
  all_km <- read.delim('companionFiles/all_km.tsv',sep="\t",header = TRUE)
  all_km <- all_km[!is.na(all_km$ligandStructureId),]
  all_km$tID <- sapply(all_km$ligandStructureId, function(x){paste(matchTIDbrenda$TID[matchTIDbrenda$ligandID == x], collapse = "/")})
  all_km <- all_km[all_km$tID != "",]
  all_km <- data.table(unique(all_km))
  all_km[,isYeast := ifelse(organism == "Saccharomyces cerevisiae", T, F),]
  setnames(all_km, old = "ecNumber", new = "EC")
  all_km <- all_km[all_km$kmValue > 0,,] #some Km values set at -999
  
  substrates_summary <- all_km[,list(log10mean = mean(log10(kmValue)), sdOflog10 = sd(log10(kmValue)), nQuant = length(kmValue), nQual = length(kmValue), modtype = NA), by = c("EC", "tID", "isYeast")]
  setcolorder(substrates_summary, names(modifiers_summary))
  
  all_affinities <- list(modifiers_summary, substrates_summary)
  all_affinities <- rbindlist(all_affinities)
  
  all_affinities[,speciesType := ifelse(is.na(modtype), "substrate", "regulator"), by = "modtype"]
  
  
  all_affinities$reactions <- sapply(all_affinities$EC, function(x){
    paste(rxnParYeast$ReactionID[grep(paste(x, '($|,)', sep = ""), rxnParYeast$EC)], collapse = "/")
  })
  
  all_affinities <- all_affinities[all_affinities$reactions != "",]
  
  write.table(all_affinities, "flux_cache/metaboliteAffinities.tsv", row.names = F, sep = "\t", quote = F)
  
}else{
  all_affinities <- read.delim("flux_cache/metaboliteAffinities.tsv")
}

#### Specify specific regulatory interactions ####
# only look at modifiers measured in yeast or with 2+ qualitative measurements in other organisms

all_brenda <- all_affinities %>% filter(isYeast | nQual > 1)

# split up degenerate species and reactions
all_brenda <- lapply(1:nrow(all_brenda), function(x){
  cbind(
    all_brenda %>% dplyr::slice(x) %>% dplyr::select(EC, modtype),
    expand.grid(tID = strsplit(all_brenda$tID[x], split = '/')[[1]],
                rxn = strsplit(all_brenda$reactions[x], split = '/')[[1]])
  )
})

all_brenda <- do.call("rbind", all_brenda)

all_brenda <- all_brenda %>% group_by(tID) %>% left_join(listTID %>% dplyr::select(tID = SpeciesType, chebi = CHEBI, name = SpeciesName), by = "tID") %>%
  mutate(hill = 1, subtype = NA, origin = "Brenda")

all_brenda <- all_brenda %>% dplyr::select(EC, rxn, tID, chebi, name, hill, modtype, subtype, origin) %>% unique()

modTable <- as.data.frame(all_brenda) #coerce back to a data.frame for compatability issues