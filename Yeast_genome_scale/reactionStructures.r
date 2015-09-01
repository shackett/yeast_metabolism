#### Intended to be sourced from FBGA.R ####
# Generate reaction forms - MM and with additional regulation

options(stringsAsFactors=FALSE)

source("cleanNalign.R")
load("flux_cache/metaboliteTables.RData")
tab_boer <- read.delim("flux_cache/tab_boer.txt")
load("flux_cache/yeast_stoi.Rdata")

# Functions

sIDtoName <- function(sID){
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

if(!(all(table(tab_boer$SpeciesType) == 1))){
  stop("Non-unique mapping between metabolomics data and tID - check metaboliteSummaries.R")
  }

boer_t <- tab_boer$SpeciesType
boer_s <- corrFile$SpeciesID[corrFile$SpeciesType %in% boer_t]
boer_filt <- row.names(stoiMat) %in% boer_s

# some compounds are small and thus probaly dont matter if they're missing
# e.g. H+, H2O
#est_names <- c("H+","water",'ammonium','diphosphate','oxygen','carbon dioxide','phosphate')

smallMol <- c("t_0398","t_0399",'t_0233','t_0332','t_0591','t_0249','t_0608')
est_s <- corrFile$SpeciesID[ corrFile$SpeciesType %in% smallMol]
est_filt <- row.names(stoiMat) %in% est_s

fin_filt <- boer_filt | est_filt #| nch_filt

stoiMat_bin <- stoiMat !=0

stoiMat_sub <- stoiMat
stoiMat_sub <- stoiMat_sub < 0

# 1) Create a dataframe rxn_table with the columns rxn_id, yname,
#factor:0,1,sub,prod (0: no reactand missing, 1: one reactand missing, sub: no substrate missing, prod: no product missing)
rxn_table<-data.frame(row.names = NULL)

#filter for "at least one measured reactand present in the reaction"
filtmin1 <- colSums(stoiMat_bin[fin_filt,]) > 0
filtmin1sub <- colSums(stoiMat_sub[fin_filt,]) > 0

rxn_none <-unname(colnames(stoiMat[,colSums(stoiMat_bin[!fin_filt,]) == 0 & filtmin1]))
rxn_one <- colnames(stoiMat[,colSums(stoiMat_bin[!fin_filt,]) == 1 & filtmin1])
rxn_sub <-unname(colnames(stoiMat[,colSums(stoiMat_sub[!fin_filt,]) ==0 & filtmin1sub]))
rxn_table <- data.frame(c(rxn_none,rxn_one,rxn_sub),row.names = NULL)
colnames(rxn_table) <- "rx"
rxn_table$rx <- as.character(rxn_table$rx)

rxn_table$type <- c(rep(0,length(rxn_none)),rep(1,length(rxn_one)),rep("sub0",length(rxn_sub)))
rxn_table$type <- as.factor(rxn_table$type)
rxn_table <- rxn_table[ !duplicated(rxn_table$rx),]

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

rm(stoiMat_bin,stoiMat_sub,est_filt,isTransport)


#### Map between reactions, Kegg RxnIDs and EC numbers ####

rxnEnzGroup <-read.table('./flux_cache/rxn_enzyme_groups.tsv',header=T,sep='\t')


#### For each reaction in the consensus reconstruction, determine which pathways and ECs are associated ######
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


#### Structural similarity #####

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
      
    }else{ # reactions with different nr of substrates than products (s != p)
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
  
  if(!mode %in% c('cc','rm')){
    error("no suitable mode selected.")
  }
  if(!("enzymeInvolved" %in% colnames(rctInput))){rctInput$enzymeInvolved <- NA }
  
  pythinput = vector(mode="character",length=nrow(rctInput))
  for (i in 1:nrow(rctInput)){
    pythinput[i] = paste(rctInput[i,c("ReactionID","Reversible","Type","SubstrateID","BindingSite","Stoi","Hill","Subtype","enzymeInvolved")],collapse='\t')
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
    
    rxnFormSave <- rctInput[ rctInput$ReactionID == substr(entry,1,6),c("ReactionID","Reversible","Type","SubstrateID","BindingSite","Stoi","Hill","Subtype","enzymeInvolved")]
    rxnFormSave$EqType <- mode
    rxnFormSave$affinityParameter <- mapply(function(a,b,c){paste0("K_", rxnFormSave$ReactionID[1], "_", a,b,c)},
      a = rxnFormSave$SubstrateID,
      b = sapply(rxnFormSave$Type, function(x){ifelse(x != "rct", paste0("_", x), "")}),
      c = sapply(rxnFormSave$enzymeInvolved, function(x){ifelse(!is.na(x), paste0("_", x), "")})
      )
    
    rxnf[[entry]]$rxnFormData <- rxnFormSave
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

# Verfiy that these reactions are from the current version of the model
rct_val_test <- rct_val_s2p %>% dplyr::select(ReactionID, SpeciesName = Substrate) %>% left_join(corrFile %>% dplyr::select(SpeciesID, SpeciesName))
rct_val_test_match <- data.frame(rxn = unique(rct_val_test$ReactionID), rxn_match = sapply(unique(rct_val_test$ReactionID), function(aRxn){
  nmetMatches <- rct_val_test %>% filter(ReactionID == aRxn) %>% inner_join(rxnFile %>% dplyr::select(SpeciesID = Metabolite, ReactionID), by = 'SpeciesID') %>%
    group_by(ReactionID.y) %>% dplyr::summarize(nmatch = n())
  bestMatches <- nmetMatches %>% filter(nmatch == max(nmatch)) %>% dplyr::select(ReactionID = ReactionID.y) %>% unlist() %>% unname()
  if(aRxn %in% bestMatches){aRxn}else{bestMatches[1]}
})
)
if(!all(rct_val_test_match$rxn == rct_val_test_match$rxn_match)){
  stop("use the results of rct_val_test_match to translate the manually updated reaction to their new names")
  }

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
rct_s2p$Subtype[ rct_s2p$Stoi < 0] <- 'substrate'
rct_s2p$Subtype[ rct_s2p$Stoi > 0] <- 'product'

#write.table(rct_s2p[,c("ReactionID","Reversible","Type","SubstrateID","BindingSite","Stoi","Hill")],file='./Python/testinp.tsv',row.names=F,col.names=F,sep='\t')

  splilistTID <- function(input,split){
    unname(sapply(input,function(x){strsplit(as.character(x),split,fixed = TRUE)}))
  }
  
  Chebi2ModelName <- function(chebis){
    # only works for model Chebis
    unname(sapply(chebis,function(x){
      listTID$SpeciesName[listTID$CHEBI == x & !is.na(listTID$CHEBI)][1]
    }))
  }

isModMeas <- function(modTable){
  
  ### check if the modulators were measured in either the relative or absolute abundance ###
  
  relMeas <- tab_boer$SpeciesType[is.na(tab_boer$concConv)]
  absMeas <- tab_boer$SpeciesType[!is.na(tab_boer$concConv)]
  
  rfil <- !is.na(modTable$tID)
  for (tIDs in unique(modTable$tID[rfil])){
    tfil <- rfil & modTable$tID == tIDs
    tIDs <- splilistTID(tIDs,'/')[[1]]
    if (any(tIDs %in% absMeas)){ meas <- 'abs'
    } else if (any(tIDs %in% relMeas)){ meas <- 'rel'
    } else {meas <- 'not'}
    modTable$measured[tfil] <- meas
  }
  return(modTable)
}

#### Add modifiers ####

if (file.exists('./flux_cache/modTable.tsv') & file.exists('./flux_cache/metaboliteAffinities.tsv')){
  modTable <- read.delim('./flux_cache/modTable.tsv',header=T,sep='\t')
} else {
  
  # 1) match rID-EC-BRENDA annotated modifiers - inhibitors and activators
  # 2) generate a table with summary of kinetic parameters - Km and Ki
  # 3) return tID modifiers -> reaction matches
  
  source('./match_brenda.R')
  
  ### Determine  the chemical similarity of BRENDA compunds as well others compunds to identify possible competitive inhibitors ###
  # make a cosimilarity matrix between all compounds in the model
  
  ## make a distmatAPt of all compounds: ##
  # use atom pair tanimoto similarity (APt)
  
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
  measF <- cid(apset) %in% listTID$CHEBI[listTID$SpeciesType %in% boer_t]
  
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
  
  ### Generate a list of metabolites which could be potential competitive inhibitors of a substrate based upon chemical similarity
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
  
  #save(list = ls(), file = "tmp.Rdata")
  
  modTable <- rbind(modTable,tmpTable)
  
  # check if measured
  
  modTable <- isModMeas(modTable)
  modTable <- modTable[!duplicated(modTable),]
  rm(simList,tmpTable)
  
  write.table(modTable,file='./flux_cache/modTable.tsv',col.names=T,row.names=F,sep='\t', quote = F)
}


#### Write the reaction laws and the rxnf structure with metabolite and protein measurements ####
# Functions

addInfo2rxnf <- function(rxnf){
  
  ## This function adds all the metabolite/enzyme/reaction info to the rxnf list
  temprown <- grep('[PCNLU][.0-9]{4}', colnames(tab_boer), value = T)
  n_c <- length(temprown)
  
  # prepare the rxnParYeast
  rownames(rxnParYeast) <- rxnParYeast$ReactionID
  
  # load the protein relative abundance, standard deviation and association to reactions
  enzyme_abund <- read.delim("./companionFiles/proteinAbundance.tsv")
  rownames(enzyme_abund) <- enzyme_abund$Gene; enzyme_abund <- enzyme_abund[,-1]
  enzyme_precision <- read.delim("./companionFiles/proteinPrecision.tsv", sep = " ")
  
  rxn_enzyme_groups <- read.delim("./flux_cache/rxn_enzyme_groups.tsv")
  
  # assert that the proteomics data is ordered the same as the metabolite data
  colnames(enzyme_abund) <- toupper(colnames(enzyme_abund))
  enzyme_abund <- enzyme_abund[,chmatch(temprown, colnames(enzyme_abund))]
  
  colnames(enzyme_precision) <- toupper(colnames(enzyme_precision))
  enzyme_precision <- enzyme_precision[,chmatch(temprown, colnames(enzyme_precision))]
  
  if(!all(colnames(enzyme_abund) == colnames(enzyme_precision))){stop("proteomics point estimates and precision are misaligned")}
  if(!all(temprown == colnames(enzyme_abund))){stop("proteomics and metabolomics are misaligned")}
  
  # Use the information from the rxnFormData to pull the rest of the data

  for (entry in names(rxnf)){
    nMet <- length(unique(rxnf[[entry]]$rxnFormData$SubstrateID))
    mettIDs <- unique(rxnf[[entry]]$rxnFormData$SubstrateID)
    # Add the metabolite table
    rxnf[[entry]]$rxnMet <- matrix(NA,n_c,nMet)
    colnames(rxnf[[entry]]$rxnMet) <- mettIDs
    row.names(rxnf[[entry]]$rxnMet) <- temprown
    rxnf[[entry]]$rxnMet <- data.frame(rxnf[[entry]]$rxnMet)
    rxnf[[entry]]$originMet <- character(length = nMet)
    names(rxnf[[entry]]$originMet) <- mettIDs
    
    for (tID in colnames(rxnf[[entry]]$rxnMet)){
      rowBoer <- which(tab_boer$SpeciesType == tID)
      if (length(rowBoer) == 0){
        rxnf[[entry]]$rxnMet[,tID] <- NA
        rxnf[[entry]]$originMet[tID] <- 'nm'
      } else {
        rxnf[[entry]]$rxnMet[,tID] <- as.numeric(tab_boer[rowBoer,colnames(tab_boer) %in% temprown]) + ifelse(is.na(tab_boer$concConv[rowBoer]), 0, log2(tab_boer[rowBoer,'concConv']))
        rxnf[[entry]]$originMet[tID] <- ifelse(is.na(tab_boer$concConv[rowBoer]), "rel", "abs")
      }
    }
    
    # Add the standard deviation of metabolites and their residual correlations
    
    measuredTid <- rxnf[[entry]]$originMet[rxnf[[entry]]$originMet != "nm"]
    
    # Remove reaction form if there arent any measured metabolites (only applicable for irreversible reactions, others should have products or substrates)
    if(length(measuredTid) == 0){
      rxnf[[entry]] <- NULL
      next
    }
    
    rxnf[[entry]]$met_SD <- matrix(unlist(remapped_SD[chmatch(names(measuredTid), tab_boer$SpeciesType),]), nrow = n_c, ncol = length(measuredTid), byrow = T)
    rownames(rxnf[[entry]]$met_SD) <- temprown; colnames(rxnf[[entry]]$met_SD) <- names(measuredTid)
    
    rxnf[[entry]]$met_Corr <- expanded_met_correlations[chmatch(names(measuredTid), tab_boer$SpeciesType),chmatch(names(measuredTid), tab_boer$SpeciesType)]
    
    # Add the rxnID
    rxnf[[entry]]$rxnID <- substr(entry,1,6)
    
    # Add the stoichiometries & binding sites
    rxnf[[entry]]$rxnStoi <- numeric(length=nMet)
    names(rxnf[[entry]]$rxnStoi) <- mettIDs
    for (tID in unique(rxnf[[entry]]$rxnFormData$SubstrateID[rxnf[[entry]]$rxnFormData$Type == 'rct'])){
      rxnf[[entry]]$rxnStoi[tID] <- unique(rxnf[[entry]]$rxnFormData$Stoi[rxnf[[entry]]$rxnFormData$Type == 'rct' & rxnf[[entry]]$rxnFormData$SubstrateID == tID])
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
    rxnf[[entry]]$enzGroup <- LETTERS[rxnEnzGroup$group[fil]] #recode enzyme groups with letters to avoid problems associated with leading numbers when interpreted in equation
    names(rxnf[[entry]]$enzGroup) <- rxnEnzGroup$enzyme[fil]
    # Protein measurements - estimate of mean and precision (inverse varianc) - precision more naturally allows calculation of weighted mean
    rxnf[[entry]]$enzymeAbund = enzyme_abund[unlist(sapply(unique(rxnEnzGroup$enzyme[fil]), function(x){grep(x, rownames(enzyme_abund))})),]
    rxnf[[entry]]$enzymePrec = enzyme_precision[unlist(sapply(unique(rxnEnzGroup$enzyme[fil]), function(x){grep(x, rownames(enzyme_abund))})),]
    
    # Aggregate protein relative abundances to form complexes weighted by measurement precision
    
    measuredEnzymeGroups <- rxnf[[entry]]$enzGroup[names(rxnf[[entry]]$enzGroup) %in% rownames(rxnf[[entry]]$enzymeAbund)]
    
    if(length(unique(measuredEnzymeGroups)) > 1){
      
      # filter subsumable protein groups to generate a parsimonious set of protein groups
      meg_matrix <- matrix(0, ncol = length(unique(names(measuredEnzymeGroups))), nrow = length(unique(measuredEnzymeGroups)))
      rownames(meg_matrix) <- unique(measuredEnzymeGroups); colnames(meg_matrix) <- unique(names(measuredEnzymeGroups))
      for(i in unique(measuredEnzymeGroups)){
        meg_matrix[rownames(meg_matrix) == i, colnames(meg_matrix) %in% names(measuredEnzymeGroups)[measuredEnzymeGroups == i]] <- 1
        }
      
      # collapse equivalent protein sets and rename set as joint
      meg_matrix_unique <- unique(meg_matrix)
      for(pset in rownames(meg_matrix)[duplicated(meg_matrix)]){
        row_match <- apply(meg_matrix_unique, 1, function(x){all(x ==  meg_matrix[rownames(meg_matrix) == pset,])})
        rownames(meg_matrix_unique)[row_match] <- paste(rownames(meg_matrix_unique)[row_match], pset, sep = "/")
        }
      
      # for each row, see if it is a subset of proteins in any other group, and has 1+ non-ascertained components -> subsume into another complex 
      prot_group_subset <- rep(NA, nrow(meg_matrix_unique))
      for(i in 1:nrow(meg_matrix_unique)){
        complement_groups <- matrix(meg_matrix_unique[-i,meg_matrix_unique[i,] == 1], ncol = sum(meg_matrix_unique[i,] == 1))
        prot_group_subset[i] <- any(apply(complement_groups, 1, function(x){all(x == 1)}))
        }
      incomplete_groups <- unique(rxnf[[entry]]$enzGroup[!(names(rxnf[[entry]]$enzGroup) %in% rownames(rxnf[[entry]]$enzymeAbund))]) #identify groups with 1+ unascertained proteins
      incomplete_unique_groups <- sapply(rownames(meg_matrix_unique), function(x){all(strsplit(x, split = '/')[[1]] %in% incomplete_groups)}) # after combining equivalent groups, are there any which are fully measured
      
      measuredEnzymeGroups <- measuredEnzymeGroups[!(measuredEnzymeGroups %in% rownames(meg_matrix_unique)[(prot_group_subset & incomplete_unique_groups)])] # remove protein groups which are both subsumable and have missing components
      
    }
    
    # take a weighted average of protein RA to determine complex RA - weighting with precision
    enzymeComplexes <- as.data.frame(matrix(NA, ncol = n_c, nrow = length(unique(measuredEnzymeGroups))))
    rownames(enzymeComplexes) <- sapply(unique(measuredEnzymeGroups), function(x){paste(x, paste(names(measuredEnzymeGroups)[measuredEnzymeGroups == x], collapse = "."), sep = "_")})
    rownames(enzymeComplexes) <- gsub('-', '_', rownames(enzymeComplexes)) #hyphens cause problems later when interpretted mathematically as subtraction
    
    colnames(enzymeComplexes) <- colnames(enzyme_abund)
    enzymeComplex_SD <- enzymeComplexes
    
    for(i in unique(measuredEnzymeGroups)){
      
      E <- rxnf[[entry]]$enzymeAbund[rownames(rxnf[[entry]]$enzymeAbund) %in% names(measuredEnzymeGroups)[measuredEnzymeGroups == i],]
      W <- rxnf[[entry]]$enzymePrec[rownames(rxnf[[entry]]$enzymeAbund) %in% names(measuredEnzymeGroups)[measuredEnzymeGroups == i],]
      
      if(nrow(E) == 1){
        
        # unambiguous case if only a single specie
        enzymeComplexes[c(1:length(unique(measuredEnzymeGroups)))[unique(measuredEnzymeGroups) == i],] <- unlist(E)
        enzymeComplex_SD[c(1:length(unique(measuredEnzymeGroups)))[unique(measuredEnzymeGroups) == i],] <- unlist(1/sqrt(W))
        
      }else{
        
        # if multiple proteins are combined, combine them according to their precision
        enzymeComplexes[c(1:length(unique(measuredEnzymeGroups)))[unique(measuredEnzymeGroups) == i],] <- apply(E * W, 2, sum)/apply(W, 2, sum)
        enzymeComplex_SD[c(1:length(unique(measuredEnzymeGroups)))[unique(measuredEnzymeGroups) == i],] <- 1/sqrt(apply(W, 2, sum))
        
      }
    }
    
    rxnf[[entry]]$enzymeComplexes <- enzymeComplexes
    
    rxnf[[entry]]$all_species_SD <- cbind(rxnf[[entry]]$met_SD, t(enzymeComplex_SD))
    
    all_species_corr <- diag(ncol(rxnf[[entry]]$all_species_SD))
    colnames(all_species_corr) <- rownames(all_species_corr) <- colnames(rxnf[[entry]]$all_species_SD)
    
    all_species_corr[1:ncol(rxnf[[entry]]$met_SD),1:ncol(rxnf[[entry]]$met_SD)] <- rxnf[[entry]]$met_Corr
    
    rxnf[[entry]]$all_species_corr <- all_species_corr
  
  }
  return(rxnf)
}

Mod2reactionEq <- function(modTable,formMode, allInhMods=F, allActMods=F){
  # formMode:
  # cc for convinience kinetics (only supports uncompetitive inhibition)
  # rm for reversible menten with substrate inhibition
  
  # if allInhMods = T:
  # write equations for comp, uncom and noncomp inhibition for all inhibitors that have a maxSim target
  
  # if allActMods = T:
  # evaluate kinetics with both convenience-kinetics style allosteric activation and conventional MM
  
  if (!allInhMods){
    modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype','subtype','hill')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
  } else {
    modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype','hill')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
    modTab$subtype[modTab$modtype == 'inh'] <- 'uncompetitive'
    tmodTab <- modTab[modTab$modtype == 'inh' & !is.na(modTab$SimMatch),]
    tmodTab$subtype <- 'noncompetitive'
    modTab <- rbind(modTab,tmodTab)
    tmodTab$subtype <- 'competitive'
    modTab <- rbind(modTab,tmodTab)
  }
  
  if (!allActMods){
    modTab[modTab$modtype == "act",'subtype'] <- "mm"
  }else{
    modTab[modTab$modtype == "act",'subtype'] <- "mm"
    tModTab <- modTab[modTab$modtype == "act",]
    tModTab[tModTab$modtype == "act",'subtype'] <- "cc"
    modTab <- rbind(modTab, tModTab)
  }
  
  rxnForms <- list()
  for (i in 1:nrow(modTab)){
    rctLine <- rct_s2p[1,]
    rctLine[] <- NA
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
      } else if (modTab$modtype[i] == "act" & modTab$subtype[i] == "mm"){
        rctLine$BindingSite <- 0
        subMod <- 'mm'
        rctLine$Subtype <- 'mm'
      } else if (modTab$modtype[i] == "act" & modTab$subtype[i] == "cc"){
        rctLine$BindingSite <- 0
        subMod <- 'cc'
        rctLine$Subtype <- 'cc'
      } else {
        stop("Invalid regulator modtype/subtype")
        }
      
      rctTab <- rbind(rctTab,rctLine)
      print(i)
      entry <- paste(modTab$rxn[i],'-',as.character(formMode),'-',strsplit(modTab$tID[i],'/')[[1]][1],'-',modTab$modtype[i],'-',subMod,ifelse(is.na(modTab$hill[i]) | modTab$hill[i] == 1, '', '_ultra'), sep='')
      rxnForms[[entry]]<- tab2ReactionForms(rctTab,formMode)[[1]]
      
    }
  }
  return(rxnForms)
}

Mod2reactionEqComplex <- function(modTable){
  
  # Allow for the creation of more complex reaction forms - e.g.
  # multiple regulators
  # isoenzymes which respond differently to regulators
  # isoenzymes which respond differently to substrates
  
  rxnForms <- list()
  regulators <- unique(modTable$TechnicalName)
  
  for(a_reg in regulators){
    
    rctTab <- rct_s2p[rct_s2p$ReactionID == modTable$rxn[modTable$TechnicalName == a_reg][1],]
    if (nrow(rctTab) == 0){
      stop(paste(regulators, "reaction not recognized"))
      next
    }
    
    formMode <- modTable$eqType[modTable$TechnicalName == a_reg][1] # rm or cc
    
    modRegulators <- modTable[modTable$TechnicalName == a_reg & modTable$modtype %in% c("act", "inh"),]
    modReg_indecies <- (1:nrow(modTable))[modTable$TechnicalName == a_reg & modTable$modtype %in% c("act", "inh")]
    if (nrow(modRegulators) != 0){
      
      infoRegulators <- rct_s2p[1:sum(modTable$TechnicalName == a_reg & modTable$modtype %in% c("act", "inh")),]
      infoRegulators[] <- NA
      infoRegulators$ReactionID <- modRegulators$rxn
      infoRegulators$SubstrateID <- modRegulators$tID
      infoRegulators$Substrate <- tIDtoName(modRegulators$tID)
      infoRegulators$Stoi <- 0
      infoRegulators$Hill <- modRegulators$Hill
      infoRegulators$Type <- modRegulators$modtype
      infoRegulators$Reversible <- rct_s2p$Reversible[ rct_s2p$ReactionID == infoRegulators$ReactionID[1]][1]
      infoRegulators$Subtype <- modRegulators$subtype
      infoRegulators$BindingSite <- modRegulators$BindingSite
      
    }
    
    # Are isoenzymes differentially impacted by kinetics or regulation ?
    
    isoenzymes <- modTable$enzymeInvolved[modTable$TechnicalName == a_reg]
    if(any(!is.na(isoenzymes))){
      if(length(isoenzymes) == 1){
        # if only one isoenzyme is provided, assume that there are others and differ kinetic w.r.t. the identified isoenzyme
        isoenzymes <- c(isoenzymes, "other")
        }
      # generate isozyme-specific reaction forms
      isoenzymeForms <- list()
      for(isoenzyme in isoenzymes){
        
        # isoenzyme specific changes in primary reaction species
        isoenzymeMod <- modTable[modTable$TechnicalName == a_reg & modTable$enzymeInvolved == isoenzyme,]
        isoenzymeMod <- isoenzymeMod[!(isoenzymeMod$modtype %in% c("act", "inh")),]
        # isoenzyme specific regulation
        isoenzymeReg <- infoRegulators[c(1:length(modReg_indecies))[modTable[modReg_indecies,]$enzymeInvolved == isoenzyme],]
         if(nrow(isoenzymeReg) != 0){isoenzymeReg$enzymeInvolved <- isoenzyme}
        
        rctTab_enzyme <- rctTab
        rctTab_enzyme$enzymeInvolved <- NA
        
        if(nrow(isoenzymeMod) != 0){
          for(i in 1:nrow(isoenzymeMod)){
            rctTab_enzyme$enzymeInvolved[rctTab_enzyme$SubstrateID == isoenzymeMod$tID[i]] <- paste0(isoenzyme)
            rctTab_enzyme$Hill[rctTab_enzyme$SubstrateID == isoenzymeMod$tID[i]] <- isoenzymeMod$Hill[i]
          }}
        
        rctTab_enzyme <- rbind(rctTab_enzyme, isoenzymeReg)
        
        isoenzymeForms[[isoenzyme]] <- tab2ReactionForms(rctTab_enzyme,formMode)[[1]]
        
        }
      
      isoenzymeList <- list()
      isoenzymeList$rxnForm <- lapply(isoenzymeForms, function(x){x$rxnForm})
      isoenzymeList$rxnFormTex <- lapply(isoenzymeForms, function(x){x$rxnFormTex})
      
      rxnFormData <- melt(lapply(isoenzymeForms, function(x){x$rxnFormData}), level = 1, id.vars = colnames( isoenzymeForms[[1]]$rxnFormData))
      rxnFormData <- unique(rxnFormData[,colnames(rxnFormData) != "L1"])
      isoenzymeList$rxnFormData <- rxnFormData
        
      rxnForms[[modTable$TechnicalName[modTable$TechnicalName == a_reg][1]]] <- isoenzymeList
      
      }else{
        ####### Fix tab2reactionForms ########
        rctTab <- rbind(rctTab, infoRegulators)
        rxnForms[[modTable$TechnicalName[modTable$TechnicalName == a_reg][1]]] <- tab2ReactionForms(rctTab,formMode)[[1]]
        
      }
    }
    return(rxnForms)
  }

#Mod2reactionEq(significant_regulation, modTable[ modTable$origin == 'Brenda', ],'rm',allInhMods=T, allActMods=T)

Mod2reactionEq_pairwise <- function(significant_regulation, modTable,formMode, allInhMods=F, allActMods=F){
  # formMode:
  # cc for convinience kinetics (only supports uncompetitive inhibition)
  # rm for reversible menten with substrate inhibition
  
  # if allInhMods = T:
  # write equations for comp, uncom and noncomp inhibition for all inhibitors that have a maxSim target
  
  # if allActMods = T:
  # evaluate kinetics with both convenience-kinetics style allosteric activation and conventional MM
  
  if (!allInhMods){
    modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype','subtype','hill')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
  } else {
    modTab <- modTable[!is.na(modTable$measured) & !duplicated(modTable[,c('rxn','tID','modtype','hill')]) & !(modTable$measured == 'not') & modTable$rxn %in% rct_s2p$ReactionID,]
    modTab[modTab$modtype == 'inh','subtype'] <- 'uncompetitive'
    tmodTab <- modTab[modTab$modtype == 'inh' & !is.na(modTab$SimMatch),]
    tmodTab$subtype <- 'noncompetitive'
    modTab <- rbind(modTab,tmodTab)
    tmodTab$subtype <- 'competitive'
    modTab <- rbind(modTab,tmodTab)
  }
  
  if (!allActMods){
    modTab[modTab$modtype == "act",'subtype'] <- "mm"
  }else{
    modTab[modTab$modtype == "act",'subtype'] <- "mm"
    tModTab <- modTab[modTab$modtype == "act",]
    tModTab[tModTab$modtype == "act",'subtype'] <- "cc"
    modTab <- rbind(modTab, tModTab)
  }
  
  # proposed secondary regulation
  modTab <- modTab %>% filter(rxn %in% significant_regulation$rxn)
  
  # proposing pairwise regulation
  
  all_pairwise_regulation <- lapply(unique(significant_regulation$rxn), function(a_rxn){
    rxn_possible <- modTab %>% filter(rxn == a_rxn)
    rxn_sig <- significant_regulation %>% filter(rxn == a_rxn)
    rxn_sig <- rxn_sig %>% left_join(rxn_possible, by = c("rxn", "tID", "modtype", "subtype"))
    
    rxn_pairs <- expand.grid(sig_n = 1:nrow(rxn_sig), propose_n = 1:nrow(rxn_possible))
    
    invalid_pairs <- lapply(1:nrow(rxn_pairs), function(x){
      if(rxn_sig[rxn_pairs[x,1],'tID'] == rxn_possible[rxn_pairs[x,2],'tID']){
        c("") # only regulation by two different metabolites is considered
      }else{
        paste(sort(c(paste(rxn_sig[rxn_pairs[x,1], c('tID', 'subtype')], collapse = "_"), 
                     paste(rxn_possible[rxn_pairs[x,2], c('tID', 'subtype')], collapse = "_"))), collapse = "_")
      }
    })
    invalid_pairs <- unlist(invalid_pairs)
    rxn_pairs <- rxn_pairs[!duplicated(invalid_pairs) & invalid_pairs != "",]
    if(nrow(rxn_pairs) == 0){
      return(NULL) 
    }else{
      lapply(1:nrow(rxn_pairs), function(x){
        rbind(rxn_sig[rxn_pairs[x,1],], rxn_possible[rxn_pairs[x,2],])
      })
    }
  })
  
  pairwise_forms <- all_pairwise_regulation[[1]]
  for(i in 2:length(all_pairwise_regulation)){
    pairwise_forms <- append(pairwise_forms, all_pairwise_regulation[[i]])
  }
  
  rxnForms <- list()
  for (i in 1:length(pairwise_forms)){
    rctLine <- rct_s2p[1:nrow(pairwise_forms[[i]]),]
    rctLine[] <- NA
    rctTab <- rct_s2p[rct_s2p$ReactionID == pairwise_forms[[i]]$rxn[1],]
    if (nrow(rctTab) > 0){
      rctLine$ReactionID <- pairwise_forms[[i]]$rxn
      rctLine$SubstrateID <- pairwise_forms[[i]]$tID
      rctLine$Substrate <- tIDtoName(unlist(strsplit(pairwise_forms[[i]]$tID,'/')))
      rctLine$Stoi <- 0
      rctLine$Hill <- pairwise_forms[[i]]$hill
      rctLine$Type <- pairwise_forms[[i]]$modtype
      rctLine$Reversible <- rct_s2p$Reversible[ rct_s2p$ReactionID == rctLine$ReactionID[1]][1]
      rctLine$Subtype <- pairwise_forms[[i]]$subtype
      
      subMod <- c()
      for(a_reg in 1:nrow(pairwise_forms[[i]])){
        
        if (pairwise_forms[[i]]$modtype[a_reg] == 'inh' & !is.na(pairwise_forms[[i]]$subtype[a_reg]) & pairwise_forms[[i]]$subtype[a_reg] == 'competitive'){
          rctLine$BindingSite[a_reg] <- rct_s2p$BindingSite[ rct_s2p$ReactionID == pairwise_forms[[i]]$rxn[a_reg] & rct_s2p$SubstrateType == pairwise_forms[[i]]$SimMatch[a_reg]]
          subMod <- c(subMod, 'comp')
        } else if (pairwise_forms[[i]]$modtype[a_reg] == 'inh' & !is.na(pairwise_forms[[i]]$subtype[a_reg]) & pairwise_forms[[i]]$subtype[a_reg] == 'noncompetitive'){
          rctLine$BindingSite[a_reg] <- rct_s2p$BindingSite[ rct_s2p$ReactionID == pairwise_forms[[i]]$rxn[a_reg] & rct_s2p$SubstrateType == pairwise_forms[[i]]$SimMatch[a_reg]]
          subMod <- c(subMod, 'noncomp')
        } else if (pairwise_forms[[i]]$modtype[a_reg] == 'inh'){
          rctLine$BindingSite[a_reg] <- 0
          subMod <- c(subMod, 'uncomp')
        } else if (pairwise_forms[[i]]$modtype[a_reg] == "act" & pairwise_forms[[i]]$subtype[a_reg] == "mm"){
          rctLine$BindingSite[a_reg] <- 0
          subMod <- c(subMod, 'mm')
        } else if (pairwise_forms[[i]]$modtype[a_reg] == "act" & pairwise_forms[[i]]$subtype[a_reg] == "cc"){
          rctLine$BindingSite[a_reg] <- 0
          subMod <- c(subMod, 'cc')
        } else {
          stop("Invalid regulator modtype/subtype")
        } 
      }
      
      rctTab <- rbind(rctTab,rctLine)
      print(i)
      entry <- paste(rctTab$ReactionID[1],'-',as.character(formMode),'-pairwise-',
                     paste(paste(unlist(strsplit(rctLine$SubstrateID, split = '/')), '-', rctLine$Type, '-', subMod,
                                 ifelse(is.na(rctLine$Hill) | rctLine$Hill == 1, '', '_ultra'), sep = ''), collapse = "+"), sep = '')
      rxnForms[[entry]]<- tab2ReactionForms(rctTab,formMode)[[1]]
      
    }
  }
  return(rxnForms)
}

##### Generate all reaction forms #####

rxnf <- list() # a list

# Generate reaaction forms just using irreversible michaelis-menten kinetics

# For reactions flowing in the canonical forward direction : evaluate when V > 0
rct_s <- rct_s2p[rct_s2p$Subtype == "substrate",]
rct_s$Reversible <- 1
rxnForms <- tab2ReactionForms(rct_s,'rm')
names(rxnForms) <- unname(sapply(names(rxnForms), function(x){sub('rm$', 'im-forward', x)}))
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}

# For reactions flowing in reverse : evaluate when V < 0
rct_p <- rct_s2p[rct_s2p$Subtype == "product",]
rct_p$Reversible <- -1
rxnForms <- tab2ReactionForms(rct_p,'rm')
names(rxnForms) <- unname(sapply(names(rxnForms), function(x){sub('rm$', 'im-reverse', x)}))
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}

## Generate reaction forms for reversible michaelis-menten and convenience kinetics
rxnForms <- tab2ReactionForms(rct_s2p,'rm') # with product inhibition
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}

rxnForms <- tab2ReactionForms(rct_s2p,'cc') # convenience kinetics
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}

## Write all the rate laws with all Brenda modifications and all possible inhibitor modes

# all BRENDA regulators (tID), measured or not for all reactions in reconstruction (rxn)
allo_table <- modTable %>% dplyr::select(rxn, name, tID, modtype, subtype, measured, origin, hill, SimMatch) %>%
  filter(origin == "Brenda")

## Supplement with gold standard regulation and manually added regulation if not already present

gs_regulation <- read.table('companionFiles/gold_standard_regulation.txt', sep = "\t", header = T)

gs_regulation <- gs_regulation %>% isModMeas() %>% mutate(modtype = ifelse(type == "inhibitor", "inh", "act"),
                                                          origin = "gold standard", hill = 1, subtype = NA) %>%
  dplyr::select(rxn = reaction, name = regulator, tID, modtype, subtype, measured, origin, hill)

gs_missing <- gs_regulation %>% dplyr::select(rxn, tID, modtype) %>% anti_join(allo_table, by = c("rxn", "tID", "modtype")) %>% left_join(gs_regulation, by = c("rxn", "tID", "modtype")) %>%
  mutate(SimMatch = NA)

allo_table <- rbind(allo_table, gs_missing)

rxnForms <- Mod2reactionEq(allo_table,'rm',T,F)
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}
Brenda_rxn_forms <- names(rxnForms) # save 

## for allosteric activation and uncompetitive inhibition, also look at a varaible hill coefficient
alloModTable <- allo_table
alloModTable$hill <- 0
rxnForms <- Mod2reactionEq(alloModTable,'rm',F,F)
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}

## Write out allosteric activator and inhibitor reaction equations for a hypothetical regulator whose patterns of variation are governed by metabolomic PCs.
deNovoRegulators <- data.frame(rxn = unique(rct_s2p$ReactionID), name = "Hypothetical Regulator", tID = "t_metX", modtype = NA, subtype = NA, measured = "rel", origin = "novelMetSearch", hill = 1)
deNovoRegulators <- rbind(deNovoRegulators, deNovoRegulators)
deNovoRegulators$modtype <- rep(c("act", "inh"), each = length(unique(rct_s2p$ReactionID)))
deNovoRegulators_bind <- deNovoRegulators; deNovoRegulators_bind$hill <- 0
deNovoRegulators <- rbind(deNovoRegulators, deNovoRegulators_bind)
rxnForms <- Mod2reactionEq(deNovoRegulators,'rm',F,F)
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}
## add information for the hypothetical metabolite, t_metX, to information list
listTID <- rbind(listTID, data.frame(SpeciesType = "t_metX", SpeciesID = "s_metX", SpeciesName = "Hypothetical metabolite X", CHEBI = "", KEGG = "", CHEBIname = "", fuzCHEBI = ""))

## Add additional regulators of interest from manual_modTable
#manualRegulators <- read.table('companionFiles/manual_modTable.txt', sep = "\t", header = T)
#rxnForms <- Mod2reactionEq(manualRegulators,'rm',F)
#for (x in names(rxnForms)){
#  rxnf[[x]] <- rxnForms[[x]]
#}
#rm(rxnForms)

## Supervised verification of multiple regulators or isoenzyme specific kinetics/regulation
manualRegulators <- read.table('companionFiles/manual_ComplexRegulation.txt', sep = "\t", header = T)
rxnForms <- Mod2reactionEqComplex(manualRegulators)
for (x in names(rxnForms)){
  rxnf[[x]] <- rxnForms[[x]]
}
rm(rxnForms)

## Look at pairwise interactions involving regulators that improve fit x all BRENDA

if(add_pairwise_regulation & file.exists("flux_cache/paramCI.Rdata")){

  # pull down al ready evaluated regulation
  load("flux_cache/paramCI.Rdata")

  library(dplyr)
  library(tidyr)
  
  #reactionInfo %>% filter(!is.na(Qvalue) & Qvalue < 0.1) %>% filter(modelType == "regulator") %>% filter(!(rMech %in% Brenda_rxn_forms))
  
  # pull out single/highly-suggestive single regulators
  significant_regulation <- reactionInfo %>% filter(!is.na(Qvalue) & Qvalue < 0.1) %>%
    filter(rMech %in% Brenda_rxn_forms) %>% dplyr::select(rxn = reaction, modification) %>%
    separate(modification, into = c("tID", "modtype", "subtype"), sep = "-")  %>%
    mutate(subtype = gsub('comp', 'competitive', subtype))
  
  rxnForms <- Mod2reactionEq_pairwise(significant_regulation, allo_table,'rm',allInhMods=T, allActMods=F)
  
  for (x in names(rxnForms)){
    rxnf[[x]] <- rxnForms[[x]]
  }
  rm(rxnForms)
}

# add the infos to the rxn structure
rxnf <- addInfo2rxnf(rxnf)
save(rxnf,file='./flux_cache/rxnf_formulametab.rdata')

# output latex reaction forms

#outtab <- character(length(rxnf))
#for (i in 1:length(rxnf)){
#  outtab[i] <- paste(names(rxnf)[i],' ---- $',rxnf[[i]]$rxnFormTex,'$\\newline',sep='')
#}
#write.table(outtab,file='./latexformula.txt',sep='\t',row.names=F,col.names=F)

