

options(stringsAsFactors=FALSE)

source('./find_chebi/find_chebi.R', chdir = T)
#load('../../RabinowitzData/Data_files/condition_model_setup.Rdata')

library(stringr)

#### Load SBML files describing metabolites, rxn stoichiometry ####

rxnFile <- read.delim("companionFiles/rxn_yeast.tsv")
rxnparFile = read.delim("companionFiles/rxn_par_yeast.tsv")
corrFile = read.delim("companionFiles/spec_yeast.tsv")
specparFile = read.delim("companionFiles/species_par_yeast.tsv")
fluxDirFile = read.delim("companionFiles/flux_dir_yeast.tsv")

customList <- parse_custom("companionFiles/customRxns.txt")

rxnFile <- rbind(rxnFile, customList$rxnFile)
rxnparFile <- rbind(rxnparFile, customList$rxnparFile)
corrFile <- rbind(corrFile, customList$corrFile)
specparFile <- rbind(specparFile, customList$specparFile)
fluxDirFile <- rbind(fluxDirFile, customList$fluxDirFile)




tab_yifan <- read.table('companionFiles/yeast_absolute_concentration_chemo.txt', sep="\t",header=TRUE)
# load the old boer table in order to have the meta information
tab_boer <- read.table('companionFiles/BoerMetabolites.txt', sep="\t",header=TRUE)
tab_boer <- tab_boer[-1,c('Metabolite','KEGG')]

# load recalculated boer mean, SD and residual correlation 
tab_boer_mean <- read.table('companionFiles/boerMean.tsv', sep="\t",header=TRUE)
tab_boer_sd <- read.table('companionFiles/boerSD.tsv', sep="\t",header=TRUE)
tab_boer_corr <- read.table('companionFiles/boerCorr.tsv', sep="\t",header=TRUE)

# delete the condtition PO4.0.061
tab_boer_mean <- tab_boer_mean[,-chmatch('PO4.0.061',colnames(tab_boer_mean))]
tab_boer_sd <- tab_boer_sd[,-chmatch('PO4.0.061',colnames(tab_boer_sd))]
# in the mean/sd files the metabolites are called slightly different, add these names to tab_boer
namenew <- sort(tab_boer$Metabolite)
nameold <- row.names(tab_boer_mean)
nameold[ chmatch('OMP',nameold)] <- "orotidine-phosphate"
nameold <- sort(nameold)
# manually check for match!
cbind(namenew,nameold)
names(nameold) <- namenew
nameold[ chmatch("orotidine-phosphate",nameold)] <-'OMP'
tab_boer$altnames <- nameold[tab_boer$Metabolite]
rm(nameold,namenew)
### functions

MatchLists <- function(list1,list2){
  #compares two lists, returns a table with two columns
  # 1. row of element in list1, 2. row of element list 2
  pairRows <- matrix(nrow=0,ncol=2)
  tmprow <- matrix(nrow=1,ncol=2)
  
  for (i in 1:length(list1)){
    idx <- which(list2 %in% list1[i] & !is.na(list2))
    if (length(idx) >0){
      for(j in 1:length(idx)){
        tmprow[1] <- i
        tmprow[2] <- idx[j]
        pairRows <- rbind(pairRows,tmprow)
        
      }
    }
  }
  return(pairRows)
}

tIDToCHEBIs <- function(tIDs){
  cfil <- grepl('chebi\\/CHEBI:',specparFile[,3])
  chebis <- sapply(tIDs,function(tID){
    sIDs <- corrFile$SpeciesID[ corrFile$SpeciesType == tID]
    fil <- specparFile[,1] %in% sIDs & cfil
    chebis <- unique(specparFile[fil,3])
    chebis <- sapply(chebis,function(url){
      strsplit(url,'CHEBI:')[[1]][2]
    })
    return(paste(chebis,collapse=';'))
  })
  return(chebis)
}

tIDToKEGGs <- function(tIDs){
  kfil <- grepl('kegg\\.compound\\/',specparFile[,3])
  keggs <- sapply(tIDs,function(tID){
    sIDs <- corrFile$SpeciesID[ corrFile$SpeciesType == tID]
    fil <- specparFile[,1] %in% sIDs & kfil
    keggs <- unique(specparFile[fil,3])
    keggs <- sapply(keggs,function(url){
      strsplit(url,'compound/')[[1]][2]
    })
    return(paste(keggs,collapse=';'))
  })
  return(keggs)
}

## start code

# Phosphate condicitions can later be estimated, so reserve a line

p <- tab_boer[1,]; p[] <- NA; p$Metabolite <- 'phosphate'; p$KEGG='C00009';
tab_boer <- rbind(tab_boer,p)

# add 3PG can be latter added with the values of PEP
p <- tab_boer[1,]; p[] <- NA; p$Metabolite <- 'phosphoenolpyruvate'; p$KEGG='C00074';
tab_boer <- rbind(tab_boer,p)

rm(p)



# according to the FBA_run_full_reco file the rows in the stoiMat file correspond to "metabolites <- unique(rxnStoi$Metabolite)"
# -> the s_#### is used as an identifier
# -> the mapping from s_#### to the name is done in the corrFile
# -> actually the t_#### id is the one unique for the compound, s_#### ids are combination of metabolite and compartments
# -> unique(rxnStoi$Metabolite) -> the s_#### identifier is used => as expected the model also takes account for the compartments

# plan: make a table with BoerName \t YifanName \t ModelName \n with the model names using the t_### identifiers

# 1) get a list of of unique t_ identifiers and their compound name

listTID <- corrFile[!duplicated(corrFile$SpeciesType),c(3,1,2)]
listTID$SpeciesName <- gsub(' \\[.*\\]$','',listTID$SpeciesName)

# only metabolites
listTID <- listTID[substr(listTID$SpeciesID,1,1) == 's',]

listTID$CHEBI <- tIDToCHEBIs(listTID$SpeciesType)
listTID$KEGG <- tIDToKEGGs(listTID$SpeciesType)
# make sure only primary chebis are used
listTID$CHEBI <- sapply(listTID$CHEBI, Chebi2pchebi)

# match the TID without CHEBI by name pefect match only matching

listTID$CHEBI[is.na(listTID$CHEBI)] <- sapply(listTID$SpeciesName[is.na(listTID$CHEBI)], function(x){MatchName2Chebi(x,1)})

listTID$CHEBIname <- sapply(listTID$CHEBI, Chebi2name)

# create the fuzzy chebi

listTID$fuzCHEBI <- sapply(listTID$CHEBI,Chebi2fuzchebi)

# # match a CHEBI for each entry that has already one Chebi measured
# listTID$fuzmCHEBI[!is.na(listTID$CHEBI)] <-sapply(listTID$SpeciesName[!is.na(listTID$CHEBI)],function(x){
# x<-MatchName2Chebi(x)
# Chebi2fuzchebi(x)
# })

#a <- listTID[ listTID$fuzmCHEBI != listTID$fuzCHEBI, ]

# add alist with the rownumber

listTID$row <- 1:length(listTID$SpeciesType)


# 2) make a list (listBoer) for all compounds in tab_boer, with rowBoer (row in tab_boer), name in Boer, keggID, chebiID, chebiname

listBoer <- data.frame(matrix(ncol=3,nrow=0))
tmpline <- data.frame(matrix(ncol=3,nrow=1))

for (i in 1:nrow(tab_boer)){
  tmpnames <- unlist(strsplit(tab_boer$Metabolite[i],'/'))
  tmpkeggs <- tab_boer$KEGG[i]
  
  
  if (length(tmpnames) == 1){
    tmpline[1] <- i
    tmpline[2] <- tmpnames
    tmpline[3] <- tmpkeggs
    listBoer <- rbind(listBoer,tmpline)
  } else {
    tmpkeggs <- unlist(strsplit(tmpkeggs,'/'))
    for (j in 1:length(tmpnames)){
      tmpline[1] <- i
      tmpline[2] <- tmpnames[j]
      tmpline[3] <- tmpkeggs[j]
      listBoer <- rbind(listBoer,tmpline)
    }
  }
}

rm(tmpkeggs)
rm(tmpnames)
# split the yet unsplit KeggIDs in listBoer

tmp <- c()
for (i in 1:length(listBoer$X1)){
  tmpkeggs <- listBoer[i,3]
  tmpkeggs <- unlist(strsplit(tmpkeggs,'/'))
  
  if (length(tmpkeggs) > 1){
    for (j in 1:length(tmpkeggs)){
      tmpline[1] <- listBoer[i,1]
      tmpline[2] <- listBoer[i,2]
      tmpline[3] <- tmpkeggs[j]
      listBoer <- rbind(listBoer,tmpline)
    }
    tmp <- c(tmp,i)
  }
}

listBoer <- listBoer[-tmp, ]

colnames(listBoer) <- c('row','name','KEGG')

# get the chebi IDs
# Try first via KEGG
listBoer$CHEBI <- sapply(listBoer$KEGG,Kegg2chebi)

# match the few yet unmatched by name
listBoer$CHEBI[is.na(listBoer$CHEBI)] <- sapply(listBoer$name[is.na(listBoer$CHEBI)],MatchName2Chebi)

# do a fuzzy chebi on the exicting chebis

listBoer$fuzCHEBI <- sapply(listBoer$CHEBI,Chebi2fuzchebi)

#match all boer names which already have a chebi again by name matching to create an alternative label

listBoer$fuzmCHEBI <- sapply(listBoer$name,MatchName2Chebi)
listBoer$fuzmCHEBI <- sapply(listBoer$fuzmCHEBI,Chebi2fuzchebi)



#### 3) make a list (listYifan) for all compounds in tab_yifan, rowYifan (row in tab_Yifan), name in Yifan, chebiID, chebiName

listYifan <- data.frame(1:length(tab_yifan$Compound),tab_yifan$Compound)
colnames(listYifan) <- c('row','name')

listYifan$CHEBI <- sapply(listYifan$name,MatchName2Chebi)

listYifan$CHEBIname <-sapply(listYifan$CHEBI,Chebi2name)

listYifan$fuzCHEBI <- sapply(listYifan$CHEBI,Chebi2fuzchebi)


#### Match the lists:

# match the Boer and TID
pairRows <- MatchLists(listBoer$KEGG,listTID$KEGG)
matchBoerTID <- cbind(listBoer$row[pairRows[,1]],listTID$row[pairRows[,2]])

pairRows <- MatchLists(listBoer$fuzCHEBI,listTID$fuzCHEBI)
matchBoerTID <- cbind(listBoer$row[pairRows[,1]],listTID$row[pairRows[,2]])

pairRows <- MatchLists(listBoer$fuzmCHEBI,listTID$fuzCHEBI)
matchBoerTID <- rbind(matchBoerTID,cbind(listBoer$row[pairRows[,1]],listTID$row[pairRows[,2]]))

matchBoerTID <- unique(matchBoerTID)


#tab_boer$Metabolite[-matchBoerTID[,1]]

# match the Boer and Yifan

pairRows <- MatchLists(listBoer$fuzCHEBI,listYifan$fuzCHEBI)
matchBoerYifan <- cbind(listBoer$row[pairRows[,1]],listTID$row[pairRows[,2]])

pairRows <- MatchLists(listBoer$fuzmCHEBI,listYifan$fuzCHEBI)
matchBoerYifan <- cbind(listBoer$row[pairRows[,1]],listTID$row[pairRows[,2]])


matchBoerYifan <- unique(matchBoerYifan)

#tab_boer$Metabolite[-matchBoerYifan[,1]]
#test <- unique(matrix(c(tab_boer$Metabolite[matchBoerYifan[ ,1]],listYifan$name[matchBoerYifan[ ,2]]),,2))
tab_yifan$Compound[-matchBoerYifan[,2]]

# match the Yifan and TID
pairRows <- MatchLists(listYifan$CHEBI,listTID$CHEBI)
matchYifanTID <- cbind(listYifan$row[pairRows[,1]],listTID$row[pairRows[,2]])

pairRows <- MatchLists(listYifan$fuzCHEBI,listTID$fuzCHEBI)
matchYifanTID <- cbind(listYifan$row[pairRows[,1]],listTID$row[pairRows[,2]])

matchYifanTID <- unique(matchYifanTID)

#tab_yifan$Compound[-matchYifanTID[,1]]


# look at the matching with primary names

test <- unique(matrix(c(tab_boer$Metabolite[matchBoerTID[ ,1]],listTID$SpeciesName[matchBoerTID[ ,2]]),,2))
#test <- unique(matrix(c(tab_boer$Metabolite[matchBoerYifan[ ,1]],listYifan$name[matchBoerYifan[ ,2]]),,2))

## create a tID to Boer table dictionary
# go through the reactions and try to translate the sIDs to Boer measurements
#make a temporary copy
tmatchBoerTID <- matchBoerTID
# every TID should correspond to only one BoerID
# -> resolve ambiguities by allways taking the match with the most similar name
ambRow <- tmatchBoerTID[duplicated(tmatchBoerTID[,2]),2]
for (rowTID in ambRow){
  fil <- tmatchBoerTID[,2] == rowTID
  rowBoer <- tmatchBoerTID[fil,1]
  tmp <- which.min(levenshteinDist(listTID$SpeciesName[rowTID],tab_boer$Metabolite[rowBoer]))
  tmp <- which(tmatchBoerTID[,1] %in% rowBoer[-tmp] & fil)
  tmatchBoerTID <- tmatchBoerTID[-tmp, ]
}

# make a tID to Boer dictionary
dicttIDBoer <- tmatchBoerTID[,1]
names(dicttIDBoer) <- listTID$SpeciesType[tmatchBoerTID[,2]]

### Make the Boer Table absolute ###
# 1) correct for the different dilution rates in the Boer-experiments and the experiments where the protein concentrations have been measured
# 2) use the Yifan table to translate the metabolites to absolute values, where possible


if (!file.exists('flux_cache/metaboliteTables.RData')){
    
  # This takes a while, so make it cacheable
  
  # the 25 chemostat conditions of interest and their actual growth rates:
  chemostatInfo <- chemostatInfo[!(chemostatInfo$condition %in% c("p0.05H1", "p0.05H2")),]
  
  # temporarly save and remove the info of the additional metabolites
  tMet <- tab_boer[-(1:nrow(tab_boer_mean)),]
  metabolomicsSD <- metabolomicsData <- tab_boer[(1:nrow(tab_boer_mean)),]
  rownames(metabolomicsData) <- metabolomicsData$Metabolite
  metabolomicsData <- cbind(metabolomicsData,tab_boer_mean[match(metabolomicsData$altnames,rownames(tab_boer_mean)),])

  metabolomicsMatrix <-as.matrix(metabolomicsData[,-c(1:3)])
  class(metabolomicsMatrix) <- "numeric"
  rownames(metabolomicsMatrix) <- metabolomicsData$Metabolite
  
  met_cond_match <- data.frame(standard = c("P", "C", "N", "L", "U"), boer = c("PO4", "GLU", "NH4", "LEU", "URA"))
  
  metabolomics_conds <- data.frame(metCond = colnames(metabolomicsMatrix))
  metabolomics_conds$met_cond <- sapply(metabolomics_conds$metCond, function(x){strsplit(x, split = "\\.")[[1]][1]})
  metabolomics_conds$cond_rename <- sapply(metabolomics_conds$met_cond, function(x){met_cond_match$standard[met_cond_match$boer == x]})
  metabolomics_conds$DR <- as.numeric(unname(sapply(metabolomics_conds$metCond, function(x){paste(strsplit(x, split = "\\.")[[1]][-1], collapse = ".")})))
  
  metabolomics_conds$cond_rename <- factor(metabolomics_conds$cond_rename, levels = met_cond_match$standard)
  
  reorder_vec <- order(metabolomics_conds$cond_rename, metabolomics_conds$cond_rename)
  
  metabolomics_conds <- metabolomics_conds[reorder_vec,]
  metabolomics_conds$DRgoal <- chemostatInfo$DRgoal
  metabolomics_conds$actualDRprot <- chemostatInfo$actualDR
  
  metabolomicsMatrix <- metabolomicsMatrix[,reorder_vec]
  colnames(metabolomicsMatrix) <- mapply(function(x,y){paste(c(x,y), collapse = "")}, x = as.character(metabolomics_conds$cond_rename), y = metabolomics_conds$DRgoal)
  colnames(metabolomicsMatrix) <- sapply(colnames(metabolomicsMatrix),function(x){
    while (nchar(x) <5){
      x <- paste(x,'0',sep='')
    }
    return(x)
  })
  
  n_c <- length(metabolomicsMatrix[1,])
  n_m <- length(metabolomicsMatrix[,1])
  
  ### Determine how many significant principal components exist in the metabolomics matrix ###
  
  matrix_svd <- svd(metabolomicsMatrix)
  plot(matrix_svd$d^2/sum(matrix_svd$d^2)) #scree-plot - fraction of variance explained by each PC
  
  library(missMDA)
  #determine how many significant principal components should be included based on repeated random sub-sampling validation
  pcrange <- c(2,18)
  npc.compare <- estim_ncpPCA(metabolomicsMatrix, ncp.min = pcrange[1], ncp.max = pcrange[2], method.cv = 'Kfold', pNA = 0.10, nbsim = 100)
  npc <- (pcrange[1]:pcrange[2])[npc.compare$criterion < (max(npc.compare$criterion) - min(npc.compare$criterion))*0.01 + min(npc.compare$criterion)][1]
  
  plot(npc.compare$criterion ~ c(pcrange[1]:pcrange[2]), pch = 16, ylab = "MS error of prediction", xlab = "number of PCs")
  abline(v = npc, col = "RED", lwd = 2)
  
  metSVD <- svd(metabolomicsMatrix)
  metMatrixProj <- metSVD$u[,1:npc] %*% diag(metSVD$d[1:npc]) %*% t(metSVD$v[,1:npc])
  
  DR_change_mat <- matrix(0, nrow = n_c, ncol = n_c)
  #colnames(DR_change_mat) <- rownames(prot_cond);
  #rownames(DR_change_mat) <- colnames(transcript.condition)
  #colnames(DR_change_mat) <- metabolomics_conds$actualDRprot
  #rownames(DR_change_mat)<- metabolomics_conds$DR
  for(cond in 1:n_c){
    #find the 2 closest DR within the same limitation
    c_match <- c(1:n_c)[metabolomics_conds$cond_rename == metabolomics_conds$cond_rename[cond]]
    flanking_match <- c_match[order(abs(metabolomics_conds[c_match,]$DR - metabolomics_conds$actualDRprot[cond]))[1:2]]
    
    lb_diff <- (metabolomics_conds$actualDRprot[cond] - metabolomics_conds$DR[flanking_match][1])/diff(metabolomics_conds$DR[flanking_match])
    DR_change_mat[flanking_match,cond] <- c((1-lb_diff), lb_diff)
  }
  
  remapped_metabolites <- metMatrixProj %*% DR_change_mat
  colnames(remapped_metabolites) <- unname(colnames(metabolomicsMatrix))
  
  metabolomicsData_remapped <- data.frame(metabolomicsData[,1:3], remapped_metabolites)
  
  tab_boer <- metabolomicsData_remapped
  
  tMet[,colnames(tab_boer)[!colnames(tab_boer) %in% colnames(tMet)]] <- NA
  tab_boer <- rbind(tab_boer,tMet)
  ### convert the relative to quantitative compounds (where available)
  # use the matchboeryifan dictionary created above
  
  nMet <- nrow(tab_boer)
  
  metOrigin <- rep('rel',nMet)
  
  metOrigin[matchBoerYifan[,1]] <- 'abs'
  
  
  refCol = which(colnames(tab_boer) == 'C0.30')
  for(i in 1:nrow(matchBoerYifan)){
    tab_boer[matchBoerYifan[i,1],4:28] <- tab_boer[matchBoerYifan[i,1],4:28] + log2(tab_yifan$c_lim_conc[matchBoerYifan[i,2]]) - tab_boer[matchBoerYifan[i,1],refCol]
  }
  
  # Recalculate the phosphate concentration with the absolute ATP/ADP values
  
  atp <- 2^ as.numeric(tab_boer[ tab_boer$Metabolite == 'ATP' & !is.na(tab_boer$Metabolite),4:28]) #M
  adp <- 2^ as.numeric(tab_boer[ tab_boer$Metabolite == 'ADP'& !is.na(tab_boer$Metabolite),4:28]) #M
  h2o <- 1 # (assumption pure water)
  
  
  dG0 = 37.9 # kJ / mol,eQuilibrator
  dG = 57 # kJ/mol Bionumbers http://bionumbers.hms.harvard.edu//bionumber.aspx?id=100775&ver=0
  
  R = 8.3144621*10^(-3) # kJ/(mol*K)
  Tm= 310.15 # K
  
  p = (atp*h2o)/adp * exp((dG0-dG)/(R*Tm))
  
  names(p) <- colnames(tab_boer[ tab_boer$Metabolite == 'ATP',4:28])
  
  tab_boer[tab_boer[,1] == 'phosphate',4:28] <- log2(p)
  metOrigin[tab_boer[,1] == 'phosphate'] <- 'abs'
  
  # add 3pg as phosphoenolpyruvate
  tab_boer[tab_boer[,1] == 'phosphoenolpyruvate',4:28]<- tab_boer[tab_boer[,1] == '3-phospho-D-glycerate',4:28]
  names(metOrigin) <- tab_boer[,1]
  
 # add phosphate and phosphoenolpyruvate to the sd table
  
  rm(atp,adp,h2o,dG,dG0,R,Tm)
  
  
  ## get the confidence interval for the boer data
  
  #boerComplete = read.table('./Data/matchCompounds/boer_data_2.txt', sep="\t",header=TRUE)
  
  
  # save the absolute boer table
  
  abs_tab_boer <- tab_boer[metOrigin == 'abs',]
  abs_tab_boer[,4:28] <- 2^abs_tab_boer[,4:28]
  
  abs_tab_boer$tID <- sapply(abs_tab_boer$Metabolite,function(x){
    paste(listTID$SpeciesType[matchBoerTID[matchBoerTID[,1] == which(tab_boer$Metabolite== x ),2]],collapse=";")
  })
  write.table(abs_tab_boer,'flux_cache/tab_boer_abs.txt',sep='\t')
  
  save(tab_boer,metOrigin,file='flux_cache/metaboliteTables.RData')
} else load('flux_cache/metaboliteTables.RData')

### Check the stoichiometric model ########################

# create the stoiMat
getStoiMat<- function(metabolites, reactions, corrFile, rxnFile){
  #create the stoichiometry matrix, indicating the change of metabolites (rows) per chemical reaction (columns)
  stoiMat <- matrix(data = 0, ncol = length(reactions), nrow = length(metabolites))
  rxnFile <- rxnFile[rxnFile$Metabolite %in% metabolites,]
  rownames(stoiMat) <- metabolites
  colnames(stoiMat) <- reactions
  for (rxn in reactions){
    rFil <- rxnFile$ReactionID == rxn & !is.na(rxnFile$StoiCoef)
    mets <- rxnFile$Metabolite[rFil]
    for (met in mets){
      if (length(rxnFile$StoiCoef[rxnFile$Metabolite == met &rFil]) > 1){
        print(c(met,rxn))
      }
      stoiMat[met,rxn] <- rxnFile$StoiCoef[rxnFile$Metabolite == met &rFil]
      
    }
  }
  return(stoiMat)
}

rxns <- sort(unique(rxnFile$ReactionID))
mets <- sort(unique(rxnFile$Metabolite[ !is.na(rxnFile$StoiCoef)]))

stoiMat <- getStoiMat(mets, rxns, corrFile, rxnFile)
stoiMat <- apply(stoiMat, c(1,2), as.numeric)
rm(rxns,mets)



## functions:

sIDtoName <-function(sID){
  
  unname(sapply(sID, function(x){
    corrFile$SpeciesName[corrFile$SpeciesID == x]
  }))
  
}

tIDtoName <-function(tID){
  unname(sapply(tID, function(x){
    listTID$SpeciesName[listTID$SpeciesType == x]
  }))
}

NametotID<-function(tID, nr=1){
  unname(sapply(tID, function(x){
    x <- gsub(' \\[.*\\]$','',x)
    x<- listTID$SpeciesType[listTID$SpeciesName == x]
    if (nr ==1){
      x <- as.integer(sub('t_','',x[1]))
    }
    return(x[1])
  }))
}

sIDtotID <- function(sIDs){
  
  unname(sapply(sIDs, function(x){
    corrFile$SpeciesType[corrFile$SpeciesID == x]
  }))
}

Name2TIDchebi <- function(name,match =1){
  # converts the name to chebi from the consensus model
  unname(sapply(name, function(x){
    x <- gsub(' \\[.*\\]$','',x)
    tmp <- unique(listTID$CHEBI[listTID$SpeciesName == x])
    # if the name is not in the consensus model, look for a name
    # by matching
    
    if (length(tmp) == 0){
      if (match == 1){
        tmp<- MatchName2Chebi(x)
      } else {tmp <- NA}
    }
    return(tmp[1])
    
  }))
  
}

name2stoi <- function(names,rxnName){
  rxnStoi <- stoiMat[,colnames(stoiMat) == rxnName][stoiMat[,colnames(stoiMat) == rxnName] != 0]
  if(reversibleRx$reversible[reversibleRx$rx == rxnName] < 0){
    rxnStoi <- rxnStoi * -1
  }
  speciesNames <- metIDtoSpec(names(rxnStoi))
  return(sapply(names,function(x){
    rxnStoi[speciesNames == x]
  }))
}


# Plan:
# 1) make a list of all t_ids which have are measured -> done
# 2) make a filter for the stoichiometric matrix for all measured substrates
# 3) make a filter for things that are not measured anyway (e.g. H+, Proteins)
# 4) make a histogram of how many reactants are missing per reaction
# 1)
boer_t <- unique(listTID$SpeciesType[matchBoerTID[,2]])

boer_s <- corrFile$SpeciesID[corrFile$SpeciesType %in% boer_t]

boer_filt <- row.names(stoiMat) %in% boer_s

# some compounds are small and thus probaly dont matter if their missing
# e.g. H+, H2O
#est_names <- c("H+","water",'ammonium','diphosphate','oxygen','carbon dioxide','phosphate')
smallMol <- c("t_0398","t_0399",'t_0233','t_0332','t_0591','t_0249','t_0608')
est_names <-smallMol

est_s <- corrFile$SpeciesID[ corrFile$SpeciesType %in% est_names]
est_filt <- row.names(stoiMat) %in% est_s

fin_filt <- boer_filt | est_filt #| nch_filt

stoiMat_bin <-stoiMat !=0

stoiMat_sub <- stoiMat
stoiMat_sub <- stoiMat_sub < 0

# 1) Create a dataframe rxn_table with the columns rxn_id, yname,
#factor:0,1,sub,prod (0: no reactand missing, 1: one reactand missing, sub: no substrate missing, prod: no product missing)
rxn_table<-data.frame(row.names = NULL)

#filter for "at least one measured reactand present in the reaction"
filtmin1 <- colSums(stoiMat_bin[fin_filt,]) > 0
filtmin1sub <- colSums(stoiMat_sub[fin_filt,]) > 0

rxn_none<-unname(colnames(stoiMat[,colSums(stoiMat_bin[!fin_filt,]) == 0 & filtmin1]))
rxn_one <- colnames(stoiMat[,colSums(stoiMat_bin[!fin_filt,]) == 1 & filtmin1])
rxn_sub <-unname(colnames(stoiMat[,colSums(stoiMat_sub[!fin_filt,]) ==0 & filtmin1sub]))
rxn_table <- data.frame(c(rxn_none,rxn_one,rxn_sub),row.names = NULL)
colnames(rxn_table) <- "rx"
rxn_table$rx <- as.character(rxn_table$rx)

rxn_table$type <- c(rep(0,length(rxn_none)),rep(1,length(rxn_one)),rep("sub0",length(rxn_sub)))
rxn_table$type <- as.factor(rxn_table$type)
rxn_table <- rxn_table[ !duplicated(rxn_table$rx),]

# # In the cases where only one reactand is missing, which metabolites are mostly missing?
# onefilt <- colSums(stoiMat_bin[!fin_filt,]) == 1
#
# #histogram per sID
# summary(as.factor(rowSums(stoiMat_bin[!fin_filt,onefilt])))
# nraffrxn_s <-rowSums(stoiMat_bin[!fin_filt,onefilt])
# fil <-nraffrxn_s >2
# row_tID_nm<- sIDtotID(row.names(stoiMat_bin[!fin_filt,]))
#
# uni_tID_nm <- unique(row_tID_nm)
#
# nraffrxn_t<- vector(mode = "integer", length = length(uni_tID_nm))
#
# for(i in 1:length(uni_tID_nm)){
# tID <- uni_tID_nm[i]
# filt <- row_tID_nm == tID
# nraffrxn_t[i] <- sum(nraffrxn_s[filt])
# }
# summary(as.factor(nraffrxn_t))
# # which things are missing in more than 4 reactions
# tIDtoName(uni_tID_nm[nraffrxn_t >4])

# now make a filter to exclude transport reactions

row_tID<- sIDtotID(row.names(stoiMat_bin))
uni_tID <- unique(row_tID)
isTransport <- rep(FALSE,dim(stoiMat)[2])

for(i in 1:length(uni_tID)){
  tID <- uni_tID[i]
  filt <- row_tID == tID
  if(sum(filt) > 1){
    idx <- which(abs(colSums(stoiMat[filt,])) < colSums(stoiMat_bin[filt,]))
    isTransport[idx]<- TRUE
  }
}

TransportRxn <- colnames(stoiMat[,isTransport])
## create the table filtered for transport reactions
# + filtered for custom reactions

cfilt <-!rxn_table$rx %in% c('r_0255','r_0256','r_0568','r_0569') # Version 7: first 2 are catalase, second 2 are diphosphatase

tfilt <- !rxn_table$rx %in% TransportRxn
rxn_table <- rxn_table[cfilt & tfilt,]

# SUMMARY
# How many reactions have at least one mesured enyzme when Including transport & redundant reactions:
#table(rxn_table$type, rxn_table$min1gene)
# How many have all measured:
#table(rxn_table$type, rxn_table$allmeas)

# How many reactions have at least one mesured enyzme when excluding transport & redundant reactions:
#table(rxn_table_filt$type, rxn_table_filt$min1gene)
# How many have all measured:
#table(rxn_table_filt$type, rxn_table_filt$allmeas)

# how many reactions are there in total?

#sum(!(uni_rxn %in% TransportRxn))

rm(stoiMat_bin,,stoiMat_sub,est_filt,est_names,isTransport)

########## Map between reactions, Kegg RxnIDs and EC numbers ####

rxnEnzGroup <-read.table('./flux_cache/rxn_enzyme_groups.tsv',header=T,sep='\t')

## For each reaction in the consensus reconstruction, determine which pathways and ECs are associated ######
if (file.exists('./flux_cache/rxnParYeast.tsv')){
  rxnParYeast <- read.table('./flux_cache/rxnParYeast.tsv',header=T,sep='\t')
} else {
  rxnParYeast <- rxnparFile
  
  # get the kegg gene to reaction mappings
  
  gene2keggID = read.delim("http://rest.kegg.jp/link/reaction/sce", header = FALSE)
  colnames(gene2keggID) = c('gene','keggID')
  gene2keggID$keggID <- str_extract(gene2keggID$keggID,'R[0-9][0-9][0-9][0-9][0-9]')
  gene2keggID$gene <- gsub('sce:','',gene2keggID$gene)
  
  
  keggID2pathways = read.delim("http://rest.kegg.jp/link/pathway/reaction", header = FALSE)
  colnames(keggID2pathways) <- c('keggID','pathway')
  keggID2pathways$pathway <- str_extract(keggID2pathways$pathway,'[0-9][0-9][0-9][0-9][0-9]')
  keggID2pathways$keggID <- str_extract(keggID2pathways$keggID,'R[0-9][0-9][0-9][0-9][0-9]')
  
  
  pathway2names = read.delim("http://rest.kegg.jp/list/pathway/sce", header = FALSE)
  colnames(pathway2names) <- c('pathway','pathname')
  pathway2names$pathname <- sapply(pathway2names$pathname, function(x){strsplit(x, split = " - Sacc")[[1]][1]})
  pathway2names$pathname <- sapply(pathway2names$pathname, function(x){strsplit(x, split = " - yeast")[[1]][1]})
  pathway2names$pathway <- str_extract(pathway2names$pathway,'[0-9][0-9][0-9][0-9][0-9]')
  
  
  keggID2reaction = read.delim("http://rest.kegg.jp/link/reaction/ec", header = FALSE)
  colnames(keggID2reaction) <- c('EC','reaction')
  keggID2reaction$reaction <- str_extract(keggID2reaction$reaction,'R[0-9][0-9][0-9][0-9][0-9]')
  keggID2reaction$EC <- gsub('ec:','',keggID2reaction$EC)
  
  
  # get the gene2Ec map from bioconductor
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("org.Sc.sgd.db")
  library("org.Sc.sgd.db")
  # convert Yeast genes to a list
  # Get the ORF identifiers that are mapped to an EC number
  mapped_genes <- mappedkeys(org.Sc.sgdENZYME)
  # Convert to a list
  ygene2EnzMap<- as.list(org.Sc.sgdENZYME[mapped_genes])
  rm(mapped_genes)
  
  # get the kegg rxn ID
  rxnParYeast$keggID <- sapply(rxnParYeast$Annotation,function(x){
    str_extract(x,'R[0-9][0-9][0-9][0-9][0-9]')
  })
  rxnParYeast$keggID[ rxnParYeast$keggID == ''] <- NA
  
  # for reactions without kegg rxn IDs from the model, get them via the genes
  fil = is.na(rxnParYeast$keggID)
  rxnParYeast$keggID[fil] <- sapply(rxnParYeast$ReactionID[fil],function(x){
    genes <- rxnEnzGroup$enzyme[ rxnEnzGroup$reaction == x]
    rxns <- unique(gene2keggID$keggID[ gene2keggID$gene %in% genes])
    if (length(rxns) == 0) {return(NA)}
    paste(rxns,collapse=',')
  })
  
  
  # get the EC numbers by the kegg reaction ID
  rxnParYeast$EC <- sapply(rxnParYeast$keggID,function(kID){
    kID <- strsplit(kID,',')[[1]]
    pID <- keggID2reaction$EC[ keggID2reaction$reaction %in% kID]
    pID <- paste(unique(pID),collapse=',')
    if (length(pID) ==0 || pID =='') {return(NA)}
    return(pID)
  })
  
  # in the cases where the kegg reaction ID did not yield an EC, get it via the gene
  fil <- is.na(rxnParYeast$EC)
  rxnParYeast$EC[fil] <- sapply(rxnParYeast$ReactionID[fil],function(x){
    genes <- rxnEnzGroup$enzyme[rxnEnzGroup$reaction == x]
    if (length(genes) == 0) {return(NA)}
    EC <- ygene2EnzMap[genes]
    EC <- unlist(EC[!sapply(EC,is.null)])
    if (length(EC) == 0) {return(NA)}
    return(paste(unique(EC),collapse=','))
  })
  
  # for the reaction without a kegg reaction ID but with an EC, try to get the keggID via EC
  fil <- is.na(rxnParYeast$keggID) & !is.na(rxnParYeast$EC)
  rxnParYeast$keggID[fil] <- sapply(rxnParYeast$EC[fil],function(EC){
    EC <- strsplit(EC,',')[[1]]
    pID <- keggID2reaction$reaction[ keggID2reaction$EC %in% EC]
    pID <- paste(unique(pID),collapse=',')
    if (length(pID) ==0 || pID =='') {return(NA)}
    return(pID)
  })
  
  # get the pathway IDs from the keggrxnIDs
  rxnParYeast$pathway <- sapply(rxnParYeast$keggID,function(kID){
    kID <- strsplit(kID,',')[[1]]
    pID <- keggID2pathways$pathway[ keggID2pathways$keggID %in% kID]
    pID <- paste(unique(pID),collapse=',')
    if (length(pID) ==0 || pID =='') {return(NA)}
    return(pID)
  })
  
  # get the pathyway names
  rxnParYeast$pathname <- sapply(rxnParYeast$pathway,function(kID){
    kID <- strsplit(kID,',')[[1]]
    pID <- pathway2names$pathname[ pathway2names$pathway %in% kID]
    pID <- paste(unique(pID),collapse='__')
    if (length(pID) ==0 || pID =='') {return(NA)}
    return(pID)
  })
  
  write.table(rxnParYeast,file='./flux_cache/rxnParYeast.tsv',col.names=T,row.names=F,sep='\t')
  rm(keggID2pathways,keggID2reaction,pathway2names,gene2keggID)
}


# parse all EC numbers in one file -> such that they can be looked up in brenda
tmpEC <- unique(unlist(strsplit(rxnParYeast$EC,",")))
tmpEC <- tmpEC[!is.na(tmpEC)]
tmpEC <- tmpEC[grep('-', tmpEC, invert = T)] # remove generic EC designations

write.table(tmpEC,file= "./flux_cache/ECnr.txt", col.names=F,row.names=F, quote = F)

rm(tmpEC)


### Structural similarity #####

# Match substrates and products that have a structural similarity and thus
# could occupy the same active sites

# a) classify reactions into 1. a+b -> c or a -> b+c 2. a->b 3. a+b -> c+d
# -> the 4. case needs further investigation as it has to be determined which product
# competes with which substrate for the binding site

# 1) make a function to get rxn - substrate - product - reversibility - num sub - num prod
# without H+ and water

rxn_reactand <- function(rxnName){
  
  rxnStoi <- stoiMat[,colnames(stoiMat) == rxnName][stoiMat[,colnames(stoiMat) == rxnName] != 0]
  speciesNames <- metIDtoSpec(names(rxnStoi))
  
  # remove small molecules
  #c("H+","water",'ammonium','diphosphate','oxygen','carbon dioxide')
  fil <- !(speciesNames %in% corrFile$SpeciesName[ corrFile$SpeciesType %in% smallMol] )
  rxnStoi <- rxnStoi[fil]
  speciesNames <- speciesNames[fil]
  
  rxnDir <- reversibleRx$reversible[reversibleRx$rx == rxnName]
  
  rctlst <- list()
  rctlst$reaction <- rxnName
  rctlst$substrate <- paste(speciesNames[rxnStoi < 0], collapse = ' / ')
  rctlst$product <- paste(speciesNames[rxnStoi > 0], collapse = ' / ')
  rctlst$rev <- rxnDir
  rctlst$numsub <- abs(sum(rxnStoi < 0))
  rctlst$numprod <- abs(sum(rxnStoi > 0))
  
  rctlst
  
}


rxn_small_mol <- function(rxnName){
  
  rxnStoi <- stoiMat[,colnames(stoiMat) == rxnName][stoiMat[,colnames(stoiMat) == rxnName] != 0]
  speciesNames <- metIDtoSpec(names(rxnStoi))
  
  # small molecules only
  fil <- (speciesNames %in% corrFile$SpeciesName[ corrFile$SpeciesType %in% smallMol] )
  rxnStoi <- rxnStoi[fil]
  speciesNames <- speciesNames[fil]
  
  rxnDir <- reversibleRx$reversible[reversibleRx$rx == rxnName]
  
  
  rctlst <- list()
  rctlst$reaction <- rxnName
  rctlst$substrate <- paste(speciesNames[rxnStoi < 0], collapse = ' / ')
  rctlst$product <- paste(speciesNames[rxnStoi > 0], collapse = ' / ')
  rctlst$rev <- rxnDir
  rctlst$numsub <- abs(sum(rxnStoi < 0))
  rctlst$numprod <- abs(sum(rxnStoi > 0))
  rctlst
  
}

# load a tool to calculate chemical similarities

#source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script.
#biocLite("ChemmineR") # Installs the package.
library("ChemmineR")

#biocLite("fmcsR")
library("fmcsR")

# functions:

cmpCpds <- function(input,inp='names',err = 0){ #input must be a vector of two names
  if (inp == 'names'){
    input <- gsub(' \\[.*\\]$','',input)
    chebi <- unlist(as.character(Name2TIDchebi(input)))
  } else if( inp == 'tid'){
    chebi <- sapply(input,function(x){
      return(as.character(listTID$CHEBI[ listTID$SpeciesType == x]))
    })
  } else if( inp != 'chebi'){
    warning('cmpCpds: wrong input type')
  }
  
  if(!chebi[1] %in% cid(apset)){
    chebi[1] <- as.character(Chebi2fuzchebi(chebi[1]))
  }
  
  if(!chebi[2] %in% cid(apset)){
    chebi[2] <- as.character(Chebi2fuzchebi(chebi[2]))
  }
  
  chebi <- chebi[chebi %in% cid(apset)]
  
  if(length(chebi) != 2){
    #print('chebi not found')
    if (err == 0 ) {
      return(0) # if chebi is a wrong input or one of the chebi is not found, 0 is returned
    } else {return(-1)}
  }else{
    return(cmp.similarity(apset[chebi[1]],apset[chebi[2]]))
    #return(fmcs(sdfset[chebi[1]],sdfset[chebi[2]],fast=TRUE)[["Tanimoto_Coefficient"]])
  }
  
}

# Read Chebi compound library

if(file.exists("./flux_cache/sdf_set_Mod.RData")){
  load("./flux_cache/sdf_set_Mod.RData")
} else {
  sdfset <- read.SDFset("../Yeast_reconstruction/Sequences/ChEBI_lite_3star.sdf")
  valid <- validSDF(sdfset)
  sdfset <- sdfset[valid]
  
  
  # make a list of all CHEBI vs IDX that are present in the GSM, Boer or Yifan
  gsmCHEBI <- unname(listTID$fuzCHEBI) #gsm
  boerCHEBI <- unname(c(listBoer$fuzCHEBI,listBoer$fuzmCHEBI))
  yifanCHEBI <- unname(listYifan$fuzCHEBI)
  
  CHEBIlist <- unique(c(gsmCHEBI,boerCHEBI,yifanCHEBI))
  
  # get all chebi that can result to the fuzzy chebi in the list
  
  CHEBIlist <- unique(c(CHEBIlist,chebifuzsyndict$secfCHEBI[ chebifuzsyndict$pfCHEBI %in% CHEBIlist]))
  
  # get all chebi that could be secondary chebi from the affected
  CHEBIlist <- unique(c(CHEBIlist,chebisecdict$secCHEBI[ chebisecdict$secCHEBI %in% CHEBIlist]))
  
  # get CHEBI from sdf
  sdfCHEBI = datablocktag(sdfset, tag="ChEBI ID")
  sdfCHEBI = unname(sapply(sdfCHEBI,function(x){
    sub('CHEBI:','',x)
  }))
  
  # reduce the dataset to only entries present in the gsm/boer/yifan + their fuzzy and secondary synonyms
  
  idx <- which(sdfCHEBI %in% CHEBIlist)
  
  sdfset <- sdfset[idx]
  
  sdfCHEBI = datablocktag(sdfset, tag="ChEBI ID")
  sdfCHEBI = unname(sapply(sdfCHEBI,function(x){
    sub('CHEBI:','',x)
  }))
  
  cid(sdfset) <- sdfCHEBI
  
  # convert these entries to ap, that can be used for structure comparison
  
  apset <- sdf2ap(sdfset)
  
  
  # clear entries that failed to return APs
  
  fil <- !sapply(as(apset, "list"), length)==1
  apset <- apset[fil]
  sdfset <- sdfset[fil]
  
  # the goal would be to have for each fuzzy chebi an entry with structure
  # -> by this it is guaranteed that each fuz chebi that has any parent with a structure, also returns a structure
  tmpCHEBI<- sapply(cid(apset),Chebi2fuzchebi)
  # get the fuzzy chebis that have no entry
  tmpCHEBI <- tmpCHEBI[!tmpCHEBI %in% cid(apset)]
  
  # get the entries of the parents of the fuzzy chebis
  tmpCHEBI <- chebifuzsyndict$secfCHEBI[chebifuzsyndict$pfCHEBI %in% tmpCHEBI]
  
  # take only the ones that are in the sdf/apset
  tmpCHEBI <- tmpCHEBI[tmpCHEBI %in% cid(apset)]
  tmpapset <- apset[as.character(tmpCHEBI)]
  tmpsdfset <- sdfset[as.character(tmpCHEBI)]
  
  # convert the tmpCHEBI again to fuzzy chebi
  tmpCHEBI <- sapply(tmpCHEBI,Chebi2fuzchebi)
  
  # make a filter for the unique
  fil <- !duplicated(tmpCHEBI)
  tmpCHEBI <- tmpCHEBI[fil]
  tmpapset <- tmpapset[fil]
  tmpsdfset <- tmpsdfset[fil]
  tmpCHEBI <- unname(as.character(tmpCHEBI))
  cid(tmpapset) <- tmpCHEBI
  cid(tmpsdfset) <-tmpCHEBI
  
  
  apset[tmpCHEBI] <-tmpapset
  sdfset <- c(sdfset,tmpsdfset)
  
  save(sdfset,apset,file="./flux_cache/sdf_set_Mod.RData")
}




#### Reaction substrates and products are matched based on their chemical similarity ####

# do a function for matching:

matchs2p <- function(rctlstinp){
  rct_s2p <- vector()
  for(i in 1:length(rctlstinp$reaction)){
    subs <- unlist(strsplit(as.character(rctlstinp$substrate[i]),' / ',fixed = TRUE))
    prod <- unlist(strsplit(as.character(rctlstinp$product[i]),' / ',fixed = TRUE))
    
    s= length(subs)
    p=length(prod)
    s2pvec <- list()
    
    if(s==p & s > 1){# reactions with equal number of substrates and products != 1
      
      # match all substrates to exactly one product
      idiff <- 1:(s*p)
      
      diffvec <- abs(apply(expand.grid(subs, prod), 1, cmpCpds))
      
      for(j in 1:s){
        
        ind = which(diffvec == max(diffvec))[1]
        sidx = ((ind-1) %% s) + 1
        pidx =(floor((ind-1) / s) +1)
        
        # as we want strictly pairwise matches, find the idx of substrates
        # and products affected and them set to -1
        diffvec[(((idiff-1) %% s) + 1) == sidx] <- -1
        diffvec[(floor((idiff-1) / s) +1) == pidx] <- -1
        
        s2pvec[j] <- paste(c(sidx,pidx), collapse = ':')
        
        
      }
      
      rct_s2p[i] <-paste(s2pvec, collapse = ', ')
      
    }else if(s==p){rct_s2p[i] <-paste(c(s,p), collapse = ':') # reactions with 1:1
                   
    }else if(s == 1){ # reactions with only one substrate
      for(j in 1:p){
        s2pvec[j] <- paste(c(s,j), collapse = ':')
      }
      rct_s2p[i] <-paste(s2pvec, collapse = ', ')
      
      
    }else if(p == 1){ # reactions with only one product
      for(j in 1:s){
        s2pvec[j] <- paste(c(j,p), collapse = ':')
      }
      rct_s2p[i] <-paste(s2pvec, collapse = ', ')
      
    }else{ # reactions with different nr of substrates then products (s != p)
      # -> match all from the larger of both (e.g. s) the other side
      
      diffvec <- abs(apply(expand.grid(subs, prod), 1, cmpCpds))
      idiff <- 1:(s*p)
      
      for(j in 1:max(s,p)){
        
        ind = which(diffvec == max(diffvec))[1]
        sidx = ((ind-1) %% s) + 1
        pidx =(floor((ind-1) / s) +1)
        
        # we want one side to be completly matched
        if(s>p){
          diffvec[(((idiff-1) %% s) + 1) == sidx] <- -1
          
        }else{
          diffvec[(floor((idiff-1) / s) +1) == pidx] <- -1
        }
        
        
        s2pvec[j] <- paste(c(sidx,pidx), collapse = ':')
      }
      
      rct_s2p[i] <-paste(s2pvec, collapse = ', ')
      
      
    }
    
  }
  names(rct_s2p) <- rctlstinp$reaction
  return(rct_s2p)
}

s2ptable <- function(rct_s2p,rctlst){
  # Arrange the information in a nicer format:
  # a matrix with rxnID, substrate, binding site, stoichiometry
  
  outmat <- matrix(ncol = 4)
  for(j in 1:length(rct_s2p)){
    x <- rct_s2p[j]
    
    subs <- unlist(strsplit(as.character(rctlst$substrate[names(x)]),' / ',fixed = TRUE))
    prod <- unlist(strsplit(as.character(rctlst$product[names(x)]),' / ',fixed = TRUE))
    
    stoisubs <- unname(name2stoi(subs,names(x)))
    stoiprod <- unname(name2stoi(prod,names(x)))
    
    pairs <- unlist(strsplit(x,', '))
    
    mat <-matrix(nrow=2*length(pairs),ncol=4)
    
    bidx = 0 # binding site id
    rowid =1
    for(i in 1:length(pairs)){
      sidx <- as.integer(unlist(strsplit(pairs[i],':')))
      pidx <- sidx[2]
      sidx <- sidx[1]
      
      
      if(sum(mat[,2] %in% c(subs[sidx],prod[pidx]))==0){
        bidx <- bidx+1
        tb <-bidx
        
      }else{
        tb <- unique(mat[mat[,2] %in% c(subs[sidx],prod[pidx]),3])
      }
      
      if(sum(mat[,2] %in% subs[sidx])==0){
        
        mat[rowid,2] <- subs[sidx]
        mat[rowid,3] <- tb
        mat[rowid,4] <- stoisubs[sidx]
        rowid = rowid+1
      }
      
      if(sum(mat[,2] %in% prod[pidx])==0){
        
        mat[rowid,2] <- prod[pidx]
        mat[rowid,3] <- tb
        mat[rowid,4] <- stoiprod[pidx]
        rowid = rowid+1
      }
    }
    
    # remove empty rows
    mat <- mat[rowSums(is.na(mat)) != ncol(mat),]
    
    mat[,1] = names(x)
    
    outmat <- rbind(outmat,mat)
  }
  outmat <- outmat[rowSums(is.na(outmat)) != ncol(outmat),]
  
  
  colnames(outmat) <- c("ReactionID","Substrate","BindingSite","Stoi")
  
  return(outmat)
  
}

addsm2rct <-function(rctlstinp,rct_inp_s2p){
  for(i in 1:length(rctlstinp$reaction)){
    x <- unname(as.character(rctlstinp$reaction[i]))
    
    subs <- unlist(strsplit(as.character(rctlstinp$substrate[x]),' / ',fixed = TRUE))
    prod <- unlist(strsplit(as.character(rctlstinp$product[x]),' / ',fixed = TRUE))
    
    rct <- c(subs,prod)
    
    templ <- rct_inp_s2p[1,]
    for(j in 1:length(rct)){
      if(sum(rct_inp_s2p$Substrate[rct_inp_s2p$ReactionID == x] == rct[j]) ==0){
        templ$ReactionID <- x
        templ$Substrate <- rct[j]
        templ$BindingSite <- max(as.integer(rct_inp_s2p$BindingSite[rct_inp_s2p$ReactionID == x]))+1
        templ$Stoi <- unname(name2stoi(rct[j],x))
        templ$Validated <- "not validated"
        
        rct_inp_s2p <- rbind(rct_inp_s2p,templ)
      }
      
    }
  }
  return(rct_inp_s2p)
}

discrCondRxn <- function(rct_inp,cutoff = 0){
  ## after some discussions with Jun, we realized that reactions
  # A+B -> C have to be further discriminated:
  # a) A+B -> C+ C (C also inhibits A & B)
  # b) 'A'+B -> 'C' (C inhibits only the active site of A)
  # c) 'A' + 'B' -> 'C' (C inhibits both active site of A and B)
  
  
  rct_inp$Stoi <- as.integer(rct_inp$Stoi)
  unirxn <- unique(rct_inp$ReactionID)
  
  for(i in 1:length(unirxn)){
    rxnfil <- rct_inp$ReactionID == unirxn[i]
    uniSite <- rct_inp$BindingSite[rxnfil]
    
    for (j in 1:length(uniSite)){
      sitefil <- rct_inp$BindingSite == uniSite[j] & rxnfil
      
      # if there are more than two reactands connected, we have a 2:1 or situation,
      # that has to be resolved
      if(sum(sitefil) > 2){
        #if the stoichiometry sums up to 0, that means we have the case a) which needs
        #no further corrections
        # in the other case it has to be discriminated between the 2 possible forms according
        # to b) or c)
        if(sum(rct_inp$Stoi[sitefil]) != 0){
          
          dir <- sign(sum(rct_inp$Stoi[sitefil]))
          
          # get the condensation products name
          cond <- rct_inp$Substrate[rct_inp$Stoi*dir <0 & sitefil]
          
          # get the educts name
          educ <- rct_inp$Substrate[rct_inp$Stoi*dir >0 & sitefil]
          
          if(length(educ) > 2){ print(c(i,j,'Attention more then 2 educts, code will fail'))}
          
          diffvec <- abs(apply(expand.grid(educ, cond), 1, cmpCpds))
          
          # if both educts are equally similar (diff < cutoff*max(diffvec)), we keep them together
          
          if(abs(diff(diffvec)) > max(diffvec)* cutoff ){
            # change less site of the less similar compound, such that it is not longer inhibited by
            # the educt
            rct_inp$BindingSite[rct_inp$Substrate ==educ[diffvec == min(diffvec)] & sitefil] <-
              as.double(rct_inp$BindingSite[rct_inp$Substrate ==educ[diffvec == min(diffvec)] & sitefil])+0.1
            rct_inp$Validated[rct_inp$Substrate ==educ[diffvec == min(diffvec)] & sitefil] <-
              "not validated"
            
          }
          
        }
        
      }
      
    }
    
  }
  return(rct_inp)
  
}

rct2Unisites <- function(rctInput){
  #makes the substrate and product sites unique (1 unique site per substrate and 1 per product)
  rctInput$BindingSite <- as.numeric(rctInput$BindingSite)
  nUni <- duplicated(rctInput[,c('ReactionID','BindingSite','Stoi')])
  while(any(nUni)){
    rctInput$BindingSite[nUni] <- rctInput$BindingSite[nUni]+0.1
    nUni <- duplicated(rctInput[,c('ReactionID','BindingSite','Stoi')])
  }
  return(rctInput)
  
}

tab2ReactionForms <- function(rctInput,mode){
  # converts the reaction table to reaction forms.
  # Possible modi are currently:
  # rm = reversible menten with feedback
  # cc = convinience kinetics (Liebermeister)
  
  if (!mode %in% c('cc','rm')){
    error("no suitable mode selected.")
  }
  rctInput
  pythinput = vector(mode="character",length=nrow(rctInput))
  for (i in 1:nrow(rctInput)){
    pythinput[i] = paste(rctInput[i,c("ReactionID","Reversible","Type","SubstrateID","BindingSite","Stoi","Hill","Subtype")],collapse='\t')
  }
  
  mode = as.character(mode)
  pythinput = append(pythinput,mode,after=0)
  pythout <- system('python genReactionForms.py',intern=T,input=pythinput)
  
  rxnf <- list()
  for (i in 1:length(pythout)){
    eval(parse(text=pythout[i]))
  }
  for (entry in names(rxnf)){
    environment(rxnf[[entry]]$rxnForm) <- .GlobalEnv;
    rxnf[[entry]]$rxnFormData <- rctInput[ rctInput$ReactionID == substr(entry,1,6),c("ReactionID","Reversible","Type","SubstrateID","BindingSite","Stoi","Hill","Subtype")]
    rxnf[[entry]]$rxnFormData$EqType <- mode
  }
  return(rxnf)
}



##### Start the matching:


# temporary create reversibleRx self, instead of loading it

reversibleRx <- fluxDirFile

colnames(reversibleRx) <- c('rx','reversible','FBAbounds')

reversibleRx$reversible[reversibleRx$reversible == 'true'] <- 0
reversibleRx$reversible[reversibleRx$reversible == 'false'] <- 1
#### match everything but small molecules

rctlst = data.frame(t(sapply(rxn_table$rx,rxn_reactand)))

rctlst <- rctlst[rctlst$numsub!= 0 &rctlst$numprod != 0,]

rct_s2p <- matchs2p(rctlst)

rct_s2p <- data.frame(s2ptable(rct_s2p,rctlst),stringsAsFactors = FALSE)

rct_s2p$Validated <- "not validated"


write.csv(rct_s2p,"./flux_cache/rxn_bindingsites.csv")

## load the document with the manual validations

rct_val_s2p <- read.table("./flux_cache/rxn_bindingsites_validated.csv",sep="\t",header=TRUE,stringsAsFactors = FALSE)
rct_val_s2p <- rct_val_s2p[,colnames(rct_val_s2p) != 'X']
# make a list with all reactions that have been validated and all reactions that have been corrected

valrxn <- unique(rct_val_s2p$ReactionID[rct_val_s2p$Validated == "manually validated"])
corrxn <- unique(rct_val_s2p$ReactionID[rct_val_s2p$Validated == "manually corrected"])

rct_s2p$Validated[rct_s2p$ReactionID %in% valrxn] <- "manually validated"



#apply the corrections
for(i in 1: length(corrxn)){
  x <- corrxn[i]
  if(sum(rct_s2p$ReactionID %in% x)){
    rct_s2p <- rct_s2p[!rct_s2p$ReactionID %in% x,]
    rct_s2p <- rbind(rct_s2p,rct_val_s2p[rct_val_s2p$ReactionID %in% x,])
  }
}



## now add the small molecules
rctlst_sm = data.frame(t(sapply(unique( rct_s2p$ReactionID),rxn_small_mol)))
rctlst_sm <- rctlst_sm[rowSums(rctlst_sm[,c(2,3)] =='') != 2,]

# add each small molecule as in a separate active sites
rct_s2p <-addsm2rct(rctlst_sm,rct_s2p)

# sort the reactions

rct_s2p <- rct_s2p[order(rct_s2p$ReactionID,rct_s2p$BindingSite,rct_s2p$Stoi,rct_s2p$Substrate,decreasing= FALSE),]

## after some discussions with Jun, we realized that reactions
# A+B -> C have to be further discriminated:
# a) A+B -> C+ C (C also inhibits A & B)
# b) 'A'+B -> 'C' (C inhibits only the active site of A)
# c) 'A' + 'B' -> 'C' (C inhibits both active site of A and B)
# -> discriminate the condensation reactions
#rct_s2p <- discrCondRxn(rct_s2p)
# Does not work

# sort the reactions -> in order to get easy understandable tables
rct_s2p <- rct_s2p[order(rct_s2p$ReactionID,rct_s2p$BindingSite,rct_s2p$Stoi,rct_s2p$Substrate,decreasing= FALSE),]

# make the sites unique
rct_s2p <- rct2Unisites(rct_s2p)

# add colums for hill coefficients (only used for activators/inhibitors)

rct_s2p$Hill = 1

# add type: reactand, activator, inhibitor

rct_s2p$Type = 'rct'

# add substrate id

rct_s2p$SubstrateID = NametotID(rct_s2p$Substrate,0)
rct_s2p$SubstrateType = gsub(' \\[.*\\]$','',rct_s2p$Substrate)

rct_s2p$Reversible = sapply(rct_s2p$ReactionID,function(x){
  reversibleRx$reversible[ reversibleRx$rx == x]
})

rct_s2p$BindingSite <- as.numeric(rct_s2p$BindingSite)
rct_s2p$Stoi <- as.numeric(rct_s2p$Stoi)
rct_s2p$Subtype <- NA
rct_s2p$Subtype[ rct_s2p$Stoi <0] <- 'substrate'
rct_s2p$Subtype[ rct_s2p$Stoi >0] <- 'product'

#write.table(rct_s2p[,c("ReactionID","Reversible","Type","SubstrateID","BindingSite","Stoi","Hill")],file='./Python/testinp.tsv',row.names=F,col.names=F,sep='\t')


#### Add modifiers ####

if (file.exists('./flux_cache/modTable.tsv') & file.exists('./flux_cache/metaboliteAffinities.tsv')){
  modTable <- read.table('./flux_cache/modTable.tsv',header=T,sep='\t')
} else {
  
  source('./match_brenda.R')
  
  library(ggplot2)
  
  splilistTID <- function(input,split){
    unname(sapply(input,function(x){strsplit(as.character(x),split,fixed = TRUE)}))
  }
  
  Chebi2ModelName <- function(chebis){
    # only works for model Chebis
    unname(sapply(chebis,function(x){
      listTID$SpeciesName[listTID$CHEBI == x & !is.na(listTID$CHEBI)][1]
    }))
  }
  
  ### check if the modulators are in are in boer or yifan ###
  isModMeas <- function(modTable){
    boerTID <- listTID$SpeciesType[matchBoerTID[!matchBoerTID[,1] %in% matchBoerYifan[,1],2]]
    yifanTID <- listTID$SpeciesType[matchBoerTID[matchBoerTID[,1] %in% matchBoerYifan[,1],2]]
    
    rfil <- !is.na(modTable$tID)
    for (tIDs in unique(modTable$tID[rfil])){
      tfil <- rfil & modTable$tID == tIDs
      tIDs <- splilistTID(tIDs,'/')[[1]]
      if (any(tIDs %in% yifanTID)){ meas <- 'abs'
      } else if (any(tIDs %in% boerTID)){ meas <- 'rel'
      } else {meas <- 'not'}
      modTable$measured[tfil] <- meas
    }
    return(modTable)
  }
  
  Mod2reactionEq <- function(modTable,formMode,allInhMods=F){
    # formMode:
    # cc for convinience kinetics (only supports uncompetitive inhibition)
    # rm for reversible menten with substrate inhibition
    
    # if allInhMods = T:
    # write equations for comp, uncom and noncomp inhibition for all inhibitors that have a maxSim target
    
    if (!allInhMods){
      modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype','subtype')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
    } else {
      modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
      modTab[modTab$modtype == 'inh','subtype'] <- 'uncompetitive'
      tmodTab <-modTab[modTab$modtype == 'inh' & !is.na(modTab$SimMatch),]
      tmodTab$subtype <- 'noncompetitive'
      modTab <- rbind(modTab,tmodTab)
      tmodTab$subtype <- 'competitive'
      modTab <- rbind(modTab,tmodTab)
    }
    
    
    rxnForms <- list()
    for (i in 1:nrow(modTab)){
      rctLine <- rct_s2p[1,]
      rctLine[] <-NA
      rctTab <- rct_s2p[rct_s2p$ReactionID == modTab$rxn[i],]
      if (nrow(rctTab) > 0){
        rctLine$ReactionID <- modTab$rxn[i]
        rctLine$SubstrateID <- modTab$tID[i]
        rctLine$Substrate <- tIDtoName(strsplit(modTab$tID[i],'/')[[1]][1])
        rctLine$Stoi <- 0
        rctLine$Hill <- modTab$hill[i]
        rctLine$Type <- modTab$modtype[i]
        rctLine$Reversible <- rct_s2p$Reversible[ rct_s2p$ReactionID == modTab$rxn[i]][1]
        rctLine$Subtype <- modTab$subtype[i]
        
        subMod <- 'allo'
        if (modTab$modtype[i] == 'inh' & !is.na(modTab$subtype[i]) & modTab$subtype[i] == 'competitive'){
          rctLine$BindingSite <- rct_s2p$BindingSite[ rct_s2p$ReactionID == modTab$rxn[i] & rct_s2p$SubstrateType == modTab$SimMatch[i]]
          subMod <- 'comp'
        } else if (modTab$modtype[i] == 'inh' & !is.na(modTab$subtype[i]) & modTab$subtype[i] == 'noncompetitive'){
          rctLine$BindingSite <- rct_s2p$BindingSite[ rct_s2p$ReactionID == modTab$rxn[i] & rct_s2p$SubstrateType == modTab$SimMatch[i]]
          subMod <- 'noncomp'
        } else if (modTab$modtype[i] == 'inh'){
          rctLine$BindingSite <- 0
          subMod <- 'uncomp'
          rctLine$Subtype <- 'uncompetitive'
        } else {
          rctLine$BindingSite <- 0
          rctLine$Subtype <- 'allosteric'
        }
        rctTab <- rbind(rctTab,rctLine)
        print(i)
        entry <- paste(modTab$rxn[i],'-',as.character(formMode),'-',strsplit(modTab$tID[i],'/')[[1]][1],'-',modTab$modtype[i],'-',subMod,sep='')
        rxnForms[[entry]]<- tab2ReactionForms(rctTab,formMode)[[1]]
        #rxnForms[[entry]]$modTable <- modTable[modTable$rxn == modTab$rxn[i] & modTable$tID == modTab$tID[i] & modTable$modtype == modTab$modtype[i] &
        # !is.na(modTable$rxn) & !is.na(modTable$tID),]
      }
    }
    return(rxnForms)
  }
  
  
  ### calculate the chemical similarity of the inhibitors & activators ###
  
  # make a cosimilarity matrix between all compounds in the model
  
  ## 1) make a distmatAPt of all compounds: ##
  # use atom pair tanimoto similarity (APt)
  #dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmatAPt.rda")
  
  if(file.exists("flux_cache/distmatAPt.rda")){
    load("flux_cache/distmatAPt.rda")
    }else{
      dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="flux_cache/distmatAPt.rda")
      load("flux_cache/distmatAPt.rda")
      rm(dummy)
      }
  
  
  distmatAPt <- 1-distmat
  rm(distmat)
  
  row.names(distmatAPt) <- cid(apset)
  colnames(distmatAPt) <- cid(apset)
  trashF <- !cid(apset) %in% Name2TIDchebi(c('H2O','hydrogen sulfide','hydroxide','H+'))
  gsmF <- cid(apset) %in% listTID$CHEBI & trashF
  measF <- cid(apset) %in% listTID$CHEBI[matchBoerTID[,2]]
  
  distmatAPt <- distmatAPt[gsmF,gsmF]
  

  
  ### get the APt similarities of the brenda compounds ###
  # (APt = atom pair tanimoto similarity)
  distfil <- listTID$CHEBI %in% row.names(distmatAPt)
  infil <- modTable$rxn %in% rct_s2p$ReactionID & !is.na(modTable$tID)
  
  for (rxn in unique(modTable$rxn[infil])){
    rfil <- infil & modTable$rxn == rxn
    for (tIDs in unique(modTable$tID[rfil])){
      
      tfil <- rfil & modTable$tID == tIDs
      tIDs <- splilistTID(tIDs,'/')[[1]]
      #rxntIDs <- NametotID(rct_s2p$Substrate[rct_s2p$ReactionID == rxn],0)
      #sim <- max(abs(apply(expand.grid(tIDs,rxntIDs), 1, function(x){cmpCpds(x,'tid')})),na.rm=T)
      rxnChebis <- as.character(listTID$CHEBI[ listTID$SpeciesType %in% rct_s2p$SubstrateID[rct_s2p$ReactionID %in% rxn] & distfil])
      chebis <- as.character(listTID$CHEBI[ listTID$SpeciesType %in% tIDs & distfil ])
      
      if (length(chebis) != 0 & length(rxnChebis) != 0){
        c
        tmp <- abs(apply(expand.grid(chebis,rxnChebis), 1, function(x){distmatAPt[x[1],x[2]]}))
        sim <- max(tmp,na.rm=T)
        matchChebi <- rxnChebis[floor((which.max(tmp)-1) / length(chebis))+1]
      } else sim <- NA
      
      modTable$maxSimAPt[tfil] <- sim
      if (sim > 0 & !is.na(sim)){
        modTable$SimMatch[tfil] <- Chebi2ModelName(matchChebi)
      }
      print(sim)
    }
  }
  
  
  
# ## make some interesting checks if inhibitors or activators are rather similar to the reactands
#
# #fil <- !is.na(modTable$maxSimAPt) & !duplicated(modTable[, c('rxn','ligandID')])
# #sum(modTable$maxSimAPt > 0.1 & modTable$maxSimAPt < 1 & modTable$modtype == 'inh' & fil)
# #sum(modTable$maxSimAPt > 0.1 & modTable$maxSimAPt < 1 & modTable$modtype == 'act' & fil)
# #hist(modTable$maxSimAPt[modTable$modtype == 'inh' & fil],1000,freq=F,ylim=c(0,20))
# #hist(modTable$maxSimAPt[modTable$modtype == 'act' & fil],1000,freq=F,ylim=c(0,20))
#
# # make a nice ggplot
# data = 'maxSimAPt'
# tplotMat <- modTable[!duplicated(modTable[,c('rxn','tID')]) & !is.na(modTable[,data]) & modTable$origin =='Brenda',]
#
# fil = tplotMat$modtype =='act'
# tHist <- hist(tplotMat[fil,data],seq(0,1,0.05),plot=F)
#
# plotMat <- data.frame(tHist$counts/sum(fil),rep('activator',length(tHist$counts)),tHist$mids)
# colnames(plotMat) <- c('Density','Type','Similarity')
# tHist <- hist(tplotMat[!fil,data],seq(0,1,0.05),plot=F)
# tmp <- data.frame(tHist$counts/sum(!fil),rep('inhibitor',length(tHist$counts)),tHist$mids)
# colnames(tmp) <- c('Density','Type','Similarity')
# plotMat <- rbind(plotMat,tmp)
#
# histPlot <- ggplot(plotMat,aes(x = Similarity,y=Density,fill = Type))
# histPlot <- histPlot + geom_bar(position = 'dodge',stat="identity",) +
# #histPlot <- histPlot + stat_bin(binwidth = 0.05,aes(y=..density..,), position = 'dodge') +
# #scale_fill_manual(name = 'Type', values = c('activator' = 'orangered2','inhibitor'='royalblue2'))+
# xlab('maximal atom pair tanimoto similarity with reactands')
#
# plot(histPlot)
# rm(tmp,tplotMat)
  #median(modTable$maxSimAPt[fil & modTable$modtype == 'inh'& modTable$maxSimAPt < 1])
  #median(modTable$maxSimAPt[fil & modTable$modtype == 'act'& modTable$maxSimAPt < 1])
  
  
  ## Textmine the Brenda commentaries:
  # Look which inhibitory commentaries contain the keywords: allo,compe
  
  fil <- modTable$modtype == 'inh' & modTable$origin == 'Brenda' & !is.na(modTable$commentary)
  idx <- which(fil)
  idx[grepl('(un-comp)|(uncomp)',modTable$commentary[fil])]
  sum(grepl('(uncomp)',modTable$commentary[fil]))
  ### Get some putative competitive inhibitors
  #get through all the reactions in the reactionfile and check if there are compounds with a similarity of > x
  # that do not belong to the the reactands and also are not the same compound as one of the reactands
  
  
  ## 1) Use APt as a similarity measurementd
  # use Apt > 0.5 for compounds to get included in the list
  
  
  simList <- list()
  for (rxn in unique(rct_s2p$ReactionID)){
    rfil <- rct_s2p$ReactionID %in% rxn
    chebi<- sapply(rct_s2p$Substrate[rfil],Name2TIDchebi)
    chebi <- as.character(chebi[chebi %in% row.names(distmatAPt)])
    simIdx <- which(distmatAPt[,chebi] > 0.5)
    simChebi <- row.names(distmatAPt)[(simIdx-1) %% nrow(distmatAPt)+1]
    simIdx <- simIdx[!simChebi %in% chebi ]
    simChebi <- simChebi[!simChebi %in% chebi ]
    
    simList[[rxn]] <- distmatAPt[,chebi][simIdx]
    names(simList[[rxn]]) <- simChebi
    
    # for duplicated chebis take only the one with the highest similarity
    while(any(duplicated(names(simList[[rxn]])))){
      dupChebi <- names(simList[[rxn]])[duplicated(names(simList[[rxn]]))][1]
      fil <- names(simList[[rxn]] ) %in% dupChebi
      simList[[rxn]][fil] <- max(simList[[rxn]][fil])
      simList[[rxn]] <- simList[[rxn]][!(duplicated(names(simList[[rxn]])) & names(simList[[rxn]]) %in% dupChebi)]
    }
  }
  
  tmpTable <- modTable[1:length(unlist(simList)),]
  tmpTable[] <- NA
  
  tmpTable$maxSimAPt <- unlist(simList)
  tmpTable$rxn <- sapply(names(unlist(simList)),function(x){strsplit(x,'\\.')[[1]][1]})
  tmpTable$chebi <- sapply(names(unlist(simList)),function(x){strsplit(x,'\\.')[[1]][2]})
  tmpTable$origin <- 'APtSimilarity'
  tmpTable$modtype <- 'inh'
  tmpTable$subtype <- 'competitive'
  tmpTable$name <- Chebi2ModelName(tmpTable$chebi)
  tmpTable$tID <- NametotID(tmpTable$name,0)
  tmpTable$hill <- 1
  
  tmpTable$SimMatch<- apply(tmpTable[,c('rxn','chebi')],1,function(x){
    rfil <- rct_s2p$ReactionID %in% x[1]
    chebi<- sapply(rct_s2p$Substrate[rfil],Name2TIDchebi)
    chebi <- as.character(chebi[chebi %in% row.names(distmatAPt)])
    idx <- which.max(distmatAPt[chebi,x[[2]]])
    matchChebi <- chebi[floor((idx-1) / length(chebi)) + 1]
    name <- Chebi2ModelName(matchChebi)
    return(name)
  })
  
  
  modTable <- rbind(modTable,tmpTable)
  
  
  
  # check if measured
  
  modTable <- isModMeas(modTable)
  modTable <- modTable[!duplicated(modTable),]
  rm(simList,tmpTable)
  
  write.table(modTable,file='./flux_cache/modTable.tsv',col.names=T,row.names=F,sep='\t')
}

#### Write the reaction laws and the rxnf structure with metabolite and protein measurements ####
# Functions

addInfo2rxnf <- function(rxnf){
  
  ## This function adds all the metabolite/enzyme/reaction infos to the rxnf list
  temprown <- colnames(tab_boer[,4:28])
  nCond <- length(temprown)
  
  # prepare the rxnParYeast
  rownames(rxnParYeast) <- rxnParYeast$ReactionID
  # load the proteomics data if necessary
  if (!'enzyme_abund' %in% ls()){
  enzyme_abund <- read.delim("./Data/matchCompounds/proteinAbundance.tsv")
  rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
  }
  # assert that the proteomics data is ordered the same as the metabolite data
  colnames(enzyme_abund) <- toupper(colnames(enzyme_abund))
  idx <- sapply(colnames(tab_boer[,4:28]),function(x){
    chmatch(x,colnames(enzyme_abund))
  })
  enzyme_abund <- enzyme_abund[,idx]
  
  # Use the information from the rxnFormData to pull the rest of the data

  for (entry in names(rxnf)){
    nMet <- length(unique(rxnf[[entry]]$rxnFormData$SubstrateID))
    mettIDs <- unique(rxnf[[entry]]$rxnFormData$SubstrateID)
    # Add the metabolite table
    rxnf[[entry]]$rxnMet <- matrix(NA,nCond,nMet)
    colnames(rxnf[[entry]]$rxnMet) <- mettIDs
    row.names(rxnf[[entry]]$rxnMet) <- temprown
    rxnf[[entry]]$rxnMet <- data.frame(rxnf[[entry]]$rxnMet)
    rxnf[[entry]]$originMet <- character(length = nMet)
    names(rxnf[[entry]]$originMet) <- mettIDs
    
    for (tID in colnames(rxnf[[entry]]$rxnMet)){
      rowBoer <- dicttIDBoer[tID]
      if (is.na(rowBoer)){
        rxnf[[entry]]$rxnMet[,tID] <- NA
        rxnf[[entry]]$originMet[tID] <- 'nm'
      } else {
        rxnf[[entry]]$rxnMet[,tID] <- as.numeric(tab_boer[rowBoer,4:28])
        rxnf[[entry]]$originMet[tID] <- metOrigin[rowBoer]
      }
    }
    # Add the rxnID
    rxnf[[entry]]$rxnID <- substr(entry,1,6)
    
    # Add the stoichiometries & binding sites
    rxnf[[entry]]$rxnStoi <- numeric(length=nMet)
    names(rxnf[[entry]]$rxnStoi) <- mettIDs
    for (tID in rxnf[[entry]]$rxnFormData$SubstrateID[rxnf[[entry]]$rxnFormData$Type == 'rct']){
      rxnf[[entry]]$rxnStoi[tID] <- rxnf[[entry]]$rxnFormData$Stoi[rxnf[[entry]]$rxnFormData$Type == 'rct' & rxnf[[entry]]$rxnFormData$SubstrateID == tID]
    }
    
    # Add the names
    rxnf[[entry]]$metNames <- tIDtoName(mettIDs)
    names(rxnf[[entry]]$metNames) <- mettIDs
    rxnf[[entry]]$listEntry <- entry
    
    # add the pathway
    fil <- rxnParYeast$ReactionID == substr(entry,1,6)
    rxnf[[entry]]$pathway <- rxnParYeast$pathname[fil]
    # add the EC numbers
    rxnf[[entry]]$EC <- rxnParYeast$EC[fil]
    
    # genes
    fil <- rxnEnzGroup$reaction == substr(entry,1,6)
    rxnf[[entry]]$genes <- paste(rxnEnzGroup$enzyme[fil] ,collapse='/')
    if(is.null(rxnf[[entry]]$genes)){rxnf[[entry]]$genes <- ''}
    # enzyme groups
    rxnf[[entry]]$enzGroup <- rxnEnzGroup$group[fil]
    names(rxnf[[entry]]$enzGroup) <- rxnEnzGroup$enzyme[fil]
    # Protein measurements
    rxnf[[entry]]$enzymeAbund = enzyme_abund[unlist(sapply(rxnEnzGroup$enzyme[fil], function(x){grep(x, rownames(enzyme_abund))})),]
  }
  return(rxnf)
}

# get the reaction forms

rxnForms <- tab2ReactionForms(rct_s2p,'rm') # with product inhibition
rxnf <- list()
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}
rxnForms <- tab2ReactionForms(rct_s2p,'cc') # convinience kinetics
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}
## Write all the rate laws with all Brenda modifications and all possible inhibitor modes

rxnForms <- Mod2reactionEq(modTable[ modTable$origin == 'Brenda', ],'rm',T)
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}


## Write out allosteric activator and inhibitor reaction equations for a generic regulator

deNovoRegulators <- data.frame(rxn = unique(rct_s2p$ReactionID), name = "Hypothetical Regulator", tID = "t_X", modtype = NA, subtype = NA, measured = "rel", origin = "novelMetSearch", hill = 1)
deNovoRegulators <- rbind(deNovoRegulators, deNovoRegulators)    
deNovoRegulators$modtype <- rep(c("act", "inh"), each = length(unique(rct_s2p$ReactionID)))

rxnForms <- Mod2reactionEq(deNovoRegulators,'rm',F)

    
Mod2reactionEq <- function(modTable,formMode,allInhMods=F){
    # formMode:
    # cc for convinience kinetics (only supports uncompetitive inhibition)
    # rm for reversible menten with substrate inhibition
    
    # if allInhMods = T:
    # write equations for comp, uncom and noncomp inhibition for all inhibitors that have a maxSim target
    
    if (!allInhMods){
      modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype','subtype')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
    } else {
      modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
      modTab[modTab$modtype == 'inh','subtype'] <- 'uncompetitive'
      tmodTab <-modTab[modTab$modtype == 'inh' & !is.na(modTab$SimMatch),]
      tmodTab$subtype <- 'noncompetitive'
      modTab <- rbind(modTab,tmodTab)
      tmodTab$subtype <- 'competitive'
      modTab <- rbind(modTab,tmodTab)
    }
    
    
    rxnForms <- list()
    for (i in 1:nrow(modTab)){
      rctLine <- rct_s2p[1,]
      rctLine[] <-NA
      rctTab <- rct_s2p[rct_s2p$ReactionID == modTab$rxn[i],]
      if (nrow(rctTab) > 0){
        rctLine$ReactionID <- modTab$rxn[i]
        rctLine$SubstrateID <- modTab$tID[i]
        rctLine$Substrate <- tIDtoName(strsplit(modTab$tID[i],'/')[[1]][1])
        rctLine$Stoi <- 0
        rctLine$Hill <- modTab$hill[i]
        rctLine$Type <- modTab$modtype[i]
        rctLine$Reversible <- rct_s2p$Reversible[ rct_s2p$ReactionID == modTab$rxn[i]][1]
        rctLine$Subtype <- modTab$subtype[i]
        
        subMod <- 'allo'
        if (modTab$modtype[i] == 'inh' & !is.na(modTab$subtype[i]) & modTab$subtype[i] == 'competitive'){
          rctLine$BindingSite <- rct_s2p$BindingSite[ rct_s2p$ReactionID == modTab$rxn[i] & rct_s2p$SubstrateType == modTab$SimMatch[i]]
          subMod <- 'comp'
        } else if (modTab$modtype[i] == 'inh' & !is.na(modTab$subtype[i]) & modTab$subtype[i] == 'noncompetitive'){
          rctLine$BindingSite <- rct_s2p$BindingSite[ rct_s2p$ReactionID == modTab$rxn[i] & rct_s2p$SubstrateType == modTab$SimMatch[i]]
          subMod <- 'noncomp'
        } else if (modTab$modtype[i] == 'inh'){
          rctLine$BindingSite <- 0
          subMod <- 'uncomp'
          rctLine$Subtype <- 'uncompetitive'
        } else {
          rctLine$BindingSite <- 0
          rctLine$Subtype <- 'allosteric'
        }
        rctTab <- rbind(rctTab,rctLine)
        print(i)
        entry <- paste(modTab$rxn[i],'-',as.character(formMode),'-',strsplit(modTab$tID[i],'/')[[1]][1],'-',modTab$modtype[i],'-',subMod,sep='')
        rxnForms[[entry]]<- tab2ReactionForms(rctTab,formMode)[[1]]
        #rxnForms[[entry]]$modTable <- modTable[modTable$rxn == modTab$rxn[i] & modTable$tID == modTab$tID[i] & modTable$modtype == modTab$modtype[i] &
        # !is.na(modTable$rxn) & !is.na(modTable$tID),]
      }
    }
    return(rxnForms)
  }
      
    
    
      
    
rxnForms <- Mod2reactionEq(modTable[ modTable$origin == 'Brenda', ],'rm',T)
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}

    
rm(rxnForms)
    
    
# add the infos to the rxn structure
rxnf <- addInfo2rxnf(rxnf)

save(rxnf,file='./flux_cache/rxnf_formulametab.rdata')




### Varia ###
## sanity check for absolute concentrations: adenylate energy charge

atp <- 2^as.numeric(tab_boer[ tab_boer$Metabolite == 'ATP',4:28])
adp <- 2^as.numeric(tab_boer[ tab_boer$Metabolite == 'ADP',4:28])
amp <- 2^as.numeric(tab_boer[ tab_boer$Metabolite == 'AMP',4:28])

energyC = (atp + 0.5 *adp)/(atp+adp+amp)

nad <- 2^tab_boer[ tab_boer$Metabolite == 'NAD+',4:28]
nadh <- 2^tab_boer[ tab_boer$Metabolite == 'NADH',4:28]

nad/(nad + nadh)
#### varia

### make an awesome heatmap with all gsm vs boer metabolites

# make a distmatAPt of all compounds:
#dummy <- cmp.cluster(db=apset, cutoff=0, save.distances="distmatAPt.rda")
#load("distmatAPt.rda")

#row.names(distmatAPt) <- cid(apset)
#colnames(distmatAPt) <- cid(apset)

## test various stuff
#plot(sdfset[name2chebi('FMN')])

#cmpCpds(c("ATP","GTP"))

# MCS
#a<-fmcs(sdfset[Name2TIDchebi("ATP")],sdfset[Name2TIDchebi("GTP")]) # gives back an MCS object

#b <- fmcs(sdfset[Name2TIDchebi("ATP")],sdfset[Name2TIDchebi("GTP")],fast=TRUE) # gives back a vector

#plotMCS(a)

# library(gplots)
#
#
# gsmF <- row.names(distmatAPt) %in% gsmCHEBI
# boerF <- row.names(distmatAPt) %in% boerCHEBI
#
# distmatAPt_names <- distmatAPt
# row.names(distmatAPt_names) <- Chebi2name(cid(apset))
# colnames(distmatAPt_names) <- row.names(distmatAPt_names)
#
# hcgsm <- hclust(as.dist(distmatAPt_names[gsmF,gsmF]), method="single")
# hcgsm[["labels"]] <- Chebi2name(cid(apset[gsmF])) # Assign correct item labels
# hcboer <- hclust(as.dist(distmatAPt_names[boerF,boerF]), method="single")
# hcboer[["labels"]] <- Chebi2name(cid(apset[boerF])) # Assign correct item labels
#
# heatmap.2(distmatAPt_names[gsmF,boerF], Rowv=as.dendrogram(hcgsm), Colv=as.dendrogram(hcboer), col=colorpanel(40, "darkblue", "yellow", "white"), density.info="none", trace="none")
#

# fun functions:

SimByName <- function(name1,name2){
  chebi1 <- as.character(MatchName2Chebi(name1))
  chebi2 <- as.character(MatchName2Chebi(name2))
  
  if (all(c(chebi1,chebi2) %in% cid(sdfset))){
    mcs <- fmcs(sdfset[chebi1],sdfset[chebi2])
    print(mcs)
    print(cmp.similarity(apset[chebi1],apset[chebi2],2))
    plotMCS(mcs)
  } else print('cpd not found!')
}


# <- sapply(a,Chebi2name)

row.names(distmatAPt_names) <- a
colnames(distmatAPt_names) <- row.names(distmatAPt_names)

hcgsm <- hclust(as.dist(distmatAPt_names[gsmF,gsmF]), method="single")
# hcgsm[["labels"]] <- Chebi2name(cid(apset[gsmF])) # Assign correct item labels
# hcboer <- hclust(as.dist(distmatAPt_names[boerF,boerF]), method="single")
# hcboer[["labels"]] <- Chebi2name(cid(apset[boerF])) # Assign correct item labels
#
# heatmap.2(distmatAPt_names[gsmF,boerF], Rowv=as.dendrogram(hcgsm), Colv=as.dendrogram(hcboer), col=colorpanel(40, "darkblue", "yellow", "white"), density.info="none", trace="none")
#

outtab <- character(length(rxnf))
for (i in 1:length(rxnf)){
  outtab[i] <- paste('$',rxnf[[i]]$rxnFormTex,'$\\newline',sep='')
}
write.table(outtab,file='./latexformula.txt',sep='\t',row.names=F,col.names=F)

