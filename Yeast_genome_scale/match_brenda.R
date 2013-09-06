# Match brenda to chebiIDs by using the synonymes from the CHEBI_complete_sdf file
options(stringsAsFactors=FALSE)


if (file.exists('./flux_cache/brendamat.tsv')){
  brendamat <- read.table('./flux_cache/brendamat.tsv',sep='\t')
} else {
  #source('./RabinowitzData/Data_files/FBA_lib.R')
  #source('./git/vzcode/find_chebi/find_chebi.R', chdir = T)
  
  
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
  
  brendamat$chebi <- sapply(brendamat$name,MatchName2Chebi)
  
  
  write.table(brendamat,'./flux_cache/brendamat.tsv',sep='\t')
}


### Match brenda and tIDs

pairRows <- MatchLists(listTID$CHEBI,brendamat$chebi)
pairRows <- rbind(MatchLists(listTID$fuzCHEBI,sapply(brendamat$chebi,Chebi2fuzchebi)),pairRows)

pairRows <- unique(pairRows)

matchTIDbrenda <- data.frame(cbind(listTID$SpeciesType[pairRows[,1]],brendamat$ligandID[pairRows[,2]]))

matchTIDbrenda <- unique(matchTIDbrenda)
colnames(matchTIDbrenda) <- c('TID','ligandID')

namesMatch <- cbind(sapply(matchTIDbrenda$TID,function(x){paste(listTID$SpeciesName[ listTID$SpeciesType %in% x],collapse='; ')}),
                    sapply(matchTIDbrenda$ligandID,function(x){paste(brendamat$name[ brendamat$ligandID %in% x],collapse='; ')}))


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

rm(all_activators,all_ki)

#commonOrganisms <- names(sort(table(all_brenda$organism), decreasing = T)[1:10])
#all_brenda <- all_brenda[all_brenda$organism %in% commonOrganisms,]

####### Summarize Ki and Km values for use in future comparisions #####

if(!file.exists("flux_cache/metaboliteAffinities.tsv")){
  
  all_brenda$ChEBI <- brendamat$chebi[chmatch(as.character(all_brenda$ligandID), as.character(brendamat$ligandID))]
  all_brenda <- all_brenda[!is.na(all_brenda$ChEBI),]
  all_brenda <- data.table(unique(all_brenda))
  all_brenda[,isYeast := ifelse(organism == "Saccharomyces cerevisiae", T, F),]
  all_brenda <- all_brenda[is.na(k) | k > 0,,]
  
  
  modifiers_summary <- all_brenda[,list(log10mean = mean(log10(k), na.rm = T), sdOflog10 = sd(log10(k), na.rm = T), nQuant = length(k[!is.na(k)]), nQual = length(k)), by = c("EC", "ChEBI", "modtype", "isYeast")]
  modifiers_summary <- modifiers_summary[!(modifiers_summary$nQual == 1 & modifiers_summary$nQuant == 0),]
  
  #ggplot(melt(table(modifiers_summary$nQuant)), aes(x = Var1, y = log2(value))) + geom_point()
  
  all_km <- read.delim('companionFiles/all_km.tsv',sep="\t",header = TRUE)
  all_km <- all_km[!is.na(all_km$ligandStructureId),]
  all_km$ChEBI <- brendamat$chebi[chmatch(as.character(all_km$ligandStructureId), as.character(brendamat$ligandID))]
  all_km <- all_km[!is.na(all_km$ChEBI),]
  all_km <- data.table(unique(all_km))
  all_km[,isYeast := ifelse(organism == "Saccharomyces cerevisiae", T, F),]
  setnames(all_km, old = "ecNumber", new = "EC")
  all_km <- all_km[all_km$kmValue > 0,,] #some Km values set at -999
  
  substrates_summary <- all_km[,list(log10mean = mean(log10(kmValue)), sdOflog10 = sd(log10(kmValue)), nQuant = length(kmValue), nQual = length(kmValue), modtype = NA), by = c("EC", "ChEBI", "isYeast")]
  setcolorder(substrates_summary, names(modifiers_summary))
  
  all_affinities <- list(modifiers_summary, substrates_summary)
  all_affinities <- rbindlist(all_affinities)
  
  all_affinities[,speciesType := ifelse(is.na(modtype), "substrate", "regulator"), by = "modtype"]
  
  
  all_affinities$reactions <- sapply(all_affinities$EC, function(x){
    paste(rxnParYeast$ReactionID[grep(paste(x, '($|,)', sep = ""), rxnParYeast$EC)], collapse = "/")
  })
  all_affinities$tIDs <- as.character(sapply(all_affinities$ChEBI, function(x){
    z = paste(listTID$SpeciesType[!is.na(listTID$CHEBI)][listTID$CHEBI[!is.na(listTID$CHEBI)] == x], collapse = "/")
    if(length(z) == 0){NA}else{z}
  }))
  
  all_affinities <- all_affinities[tIDs != "",,]
  
  write.table(all_affinities, "flux_cache/metaboliteAffinities.tsv", sep = "\t")
}

###########translate get all genes for all reactions to EC numbers
# load Map between reactions to EC numbers
# (from match_compounds)

row.names(rxnParYeast) <- rxnParYeast$ReactionID
rxnMatch <- list()

for (rxn in rxnParYeast$ReactionID[!is.na(rxnParYeast$EC)]){
    idx <- which(all_brenda$EC %in% strsplit(rxnParYeast[rxn,'EC'],',')[[1]])
    if (length(idx) > 0){
      rxnMatch[[rxn]] <- idx
    }
}

# initialize the matrix
rxnVec <- unlist(unlist(sapply(names(rxnMatch),function(x){
  rep(x,length(rxnMatch[[x]]))
})))

rxnVec <-data.frame(rxnVec)
colnames(rxnVec) <- 'rxn'
modTable <- merge(rxnVec,modTable,all.x =T,sort=F)
rm(rxnVec)

rxnMatchIdx <- unlist(rxnMatch)

modTable[,c('EC','name','ligandID','modtype','organism','literature','commentary','k','kmax')] <-
  all_brenda[rxnMatchIdx,c('EC','name','ligandID','modtype','organism','literature','commentary','k','kmax')]


for (ligID in unique(modTable$ligandID[!is.na(modTable$ligandID)])){
  tIDs <- paste(matchTIDbrenda$TID[ matchTIDbrenda$ligandID == ligID],collapse='/')
  if (tIDs != ''){
    modTable$tID[ modTable$ligandID == ligID] <- tIDs
  }
}

# Split the entries with more then 1 tID, such that each entry only has 1 tID

for(idx in grep('/',modTable$tID)){
  splittID <- strsplit(modTable$tID[idx],'/')[[1]]
  modTable$tID[idx] <- splittID[1]
  tModLines <- modTable[rep(idx,length(splittID)-1),]
  tModLines$tID <- splittID[-1]
  modTable <- rbind(modTable,tModLines)
}

modTable$origin <- 'Brenda'
modTable$hill <- 1


