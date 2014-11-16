# This is a collection of functions used to find CHEBI IDs by name, KEGG ID etc.

# The script must  with "source('script location', chdir = T)


# It loads several tables acquired from CHEBI, KEGG or manual anotations and creates 4 tables:

# chebisecdict -> a table with the primary CHEBI and the secondary CHEBI IDs according to CHEBI
# chebisyndict -> a table primary CHEBI and a list of synonym names
# chebinamedict -> a table with primary CHEBI and their 'official' name according to CHEBI
# chebifuzsyndict -> a table with a primary 'fuzzy' CHEBI and a list of so called 'fuzzy' primary CHEBI synonymes
#                   that have been found by grouping CHEBIs that have at least one synonym in common and the same
#                   chemical formula
# chebikeggdict -> a table with primary CHEBI and their associated KEGG id that have been found by
#                  a previously match with the KEGG rest api (which is rather slow)


# The relevant functions are:

# MatchName2Chebi(name, cutoff = 0.8) -> match a chebi to a given name using the synonym list. 
#                                        A cutoff of 1 means only perfect matches allowed

# Name2chebi(name) -> matches a recommended CHEBI name to a CHEBI. ONLY works for recommended Chebi names,
#                      for other cases consider MatchName2Chebi

# Chebi2name(chebiID) -> looks up the recommended CHEBI name for a primary CHEBI

# Chebi2pchebi(chebiID) -> looks up the recommended primary CHEBI for a CHEBI

# Chebi2fuzchebi(chebiID) -> converts a CHEBI to its fuzzy primary synonym 
#                   => these not so specific CHEBI can be use to match lists with more specific CHEBIs
#                   => stereochemistry and charges are neglected to a large extend

# Chebi2kegg(chebiID) -> Converts a CHEBIID to a KEGG ID by looking it up in the previouly created table.
#                       => if this does not work (return NA) consider  to use Querychebi2kegg

# Kegg2chebi(keggID) -> the reverse of Chebi2kegg

# Querychebi2kegg -> querys the kegg from the CHEBI via the Kegg rest api. This is slow thus consider trying
#                     Chebi2kegg first.


# Querykegg2chebi - > the reverse of Querychebi2kegg

# Chebi2shortname(chebiID) -> looks up the shortest Synonym for a CHEBI using the fuzzy chebi





##
options(stringsAsFactors=FALSE)

path = '.'

library(RecordLinkage)
library(Matrix)


# 2 will overwrite previously saved data -> make sure to start with an empty workspace

if (file.exists(paste(path,'/cache/find_chebiData.Rdata',sep=''))){
  load(paste(path,'/cache/find_chebiData.Rdata',sep=''))
} else {
  
  ## Functions
  
  GetClosestMatch = function(string, stringVector){
    
    out = list()
    idx <- match(string,stringVector)
    
    if (length(idx) >0){  
      out[[1]] = idx
      out[[2]] = 1
    } else {
      
      distance = levenshteinSim(string, stringVector)
      md <- max(distance,na.rm = TRUE)
      out = list()
      out[[1]] = which(distance == md)
      out[[2]] = md
    }
    return(out)
  }
  
  MatchName2Chebi <- function(name,cutoff = 0.8){
    # MatchName2Chebi(name, cutoff = 0.8) -> match a chebi to a given name using the synonym list. 
    #                                        A cutoff of 1 means only perfect matches allowed
    name <- tolower(as.character(name))
    name <- gsub("[^a-zA-Z0-9+]","",name)
    name <- gsub(" ","",name)
    
    idx <- which(chebisyndict$name == name)
    
    if (length(idx) >0){ 
      print(1)
      return(chebisyndict$CHEBI[unlist(idx[1])])
      
    } else if (cutoff < 1 & !is.na(name)){
      
      distance = levenshteinSim(name, chebisyndict$name)
      md <- max(distance,na.rm = TRUE)
      print(md)
      if(md > cutoff){
        return(chebisyndict$CHEBI[distance == md][1])
      } else return(NA)
    } else {return(NA)}
    
  }
  

  
  Chebi2pchebi <- function(chebi){
    # Chebi2pchebi(chebiID) -> looks up the recommended primary CHEBI for a CHEBI
    chebi <- as.integer(chebi)
    
    # converts the chebi id to a primary chebi ID
    pchebi <- chebisecdict$CHEBI[chebisecdict$secCHEBI==chebi]
    if(length(pchebi) >0){
      return(pchebi[1])
    } else { 
      print('Warning no primary chebi found!')
      print(chebi)
      return(chebi)
    }
  }
  
  Chebi2fuzchebi <- function(chebi){
    # Chebi2fuzchebi(chebiID) -> converts a CHEBI to its fuzzy primary synonym 
    #                   => these not so specific CHEBI can be use to match lists with more specific CHEBIs
    #                   => stereochemistry and charges are neglected to a large extend
    chebi <- as.integer(chebi)
    chebi <- Chebi2pchebi(chebi)
    # converts the chebi id to a primary fuzzy chebi ID 
    # (=all chebi which share a high quality annotated synonym are united to one chebi )
    # (e.g. to make atp2- and atp have the same chebi)
    pchebi <- chebifuzsyndict$pfCHEBI[chebifuzsyndict$secfCHEBI == chebi]
    if(length(pchebi) >0){
      return(pchebi[1])
    } else { 
      return(chebi)
    }
  }
  
  Chebi2kegg <- function(chebi){
    # Chebi2kegg(chebiID) -> Converts a CHEBIID to a KEGG ID by looking it up in the previouly created table.
    #                       => if this does not work (return NA) consider  to use Querychebi2kegg
    
    chebi <- as.integer(chebi)
    chebi <- Chebi2pchebi(chebi)
    kegg <- chebikeggdict$KEGG[chebikeggdict$CHEBI == chebi]
    if(length(kegg) >0){
      return(kegg[1])
    } else { 
      return(NA)
    }
  }
  
  Kegg2chebi <- function(kegg){
    # Kegg2chebi(keggID) -> the reverse of Chebi2kegg
    chebi <- chebikeggdict$CHEBI[chebikeggdict$KEGG == kegg]
    if(length(chebi)>0){
      return(chebi[1])
    } else { 
      return(NA)
    }
  }
  
  Chebi2name <- function(chebi){
    # Chebi2name(chebiID) -> looks up the recommended CHEBI name for a primary CHEBI
    chebi <- Chebi2pchebi(chebi)
    name <- chebinamedict$name[chebinamedict$CHEBI == chebi]
    if(length(name) >0){
      return(name[1])
    } else { 
      print('Warning, no primary chebi name')
      print(chebi)
      return(NA)
    }
  }
  
  Name2chebi <- function(name){
    # Name2chebi(name) -> matches a recommended CHEBI name to a CHEBI. ONLY works for recommended Chebi names,
    #                      for other cases consider MatchName2Chebi
    
    chebi <- chebinamedict$CHEBI[chebinamedict$name == name]
    if(length(chebi) >0){
      return(chebi[1])
    } else { 
      print('Only works for primary CHEBI names, for others use MatchName2Chebi')
      return(NA)
    }
  }
  
  Chebi2shortname <- function(chebiID){
    # Chebi2shortname(chebiID) -> looks up the shortest Synonym for a CHEBI using the fuzzy Chebi
    # look up the fuzzy chebi id
    chebiID <- Chebi2fuzchebi(chebiID)
    # look up the chebiIDs that belong to this fuzzy CHEBI
    
    if (chebiID %in%chebifuzsyndict$pfCHEBI ){
    chebiID <- chebifuzsyndict$secfCHEBI[chebifuzsyndict$pfCHEBI == chebiID]
    }
    # find the shortest synonym
    
    shorsyn <- chebisyndict$name[chebisyndict$CHEBI %in% chebiID & (chebisyndict$origin %% 1000) > 49]
    return(shorsyn[which.min(nchar(shorsyn))])
  }
  
  Querychebi2kegg <- function(chebiID){
    # Querychebi2kegg -> querys the kegg from the CHEBI via the Kegg rest api. This is slow thus consider trying
    #                     Chebi2kegg first.
    keggID <- NA
    if(!is.na(chebiID)){
      keggTab <- try(read.table(paste('http://rest.kegg.jp/conv/compound/chebi:',chebiID,sep="")),silent = TRUE)
      if(class(keggTab) != "try-error"){
        
        keggID <-as.character(keggTab[1,2])
      }
      
    }
    return(sub('cpd:','',keggID)) 
  }
  
  
  Querykegg2chebi <- function(keggID){
    # Querykegg2chebi - > the reverse of Querychebi2kegg
    chebiID <- NA
    if(!is.na(keggID)){
      chebTab <- try(read.table(paste('http://rest.kegg.jp/conv/chebi/compound:',keggID,sep="")),silent = TRUE)
      if(class(chebTab) != "try-error"){
        
        chebiID <-as.character(chebTab[1,2])
      }
      
    }
    return(sub('chebi:','',chebiID)) 
  }
  
  
Chebi2synkegg <- function(chebi){
  # this function looks whether any of the fuzzy synonyms have a kegg and
  # assigns the kegg with the compound of the most similar name
  fuzChebi <- Chebi2fuzchebi(chebi)
  synChebis <-chebifuzsyndict$secfCHEBI[ chebifuzsyndict$pfCHEBI %in% fuzChebi]
  synKegg <-  chebikeggdict$KEGG[ chebikeggdict$CHEBI %in% synChebis]
  synKegg <- synKegg[!is.na(synKegg)]
  
  if (length(synKegg) > 0){
    return(synKegg[1])
  } else {return(NA)}
}

  
  # Temporary functions:
  
  Adddict <- function(inpFile, dict = matrix(nrow=0,ncol=3)){
    inp <- read.table(inpFile, sep="\t",header = TRUE)
    colnames(inp) <- c('CHEBI','name','origin')
    dict <- rbind(dict,inp)
    
    return(dict)
  }
  
  RateName <- function(namesinp){
    
    # rate the names according to:
    # If it is a NAME -> 200
    # if it is a synonyme -> look at the source: KEGG, ChEBI, MetaCyc, IUPAC, UniProt = 50
    #                                            SUBMITTER = 20
    #                                            all others = 30
    rating <- numeric(length(namesinp$ID))
    q1 <- c("KEGG COMPOUND","ChEBI","KEGG DRUG","KEGG GLYCAN","MetaCyc","IUPAC","UniProt")
    q2 <- c("SUBMITTER")
    
    for(i in 1:length(namesinp$ID)){
      if(namesinp$TYPE[i] == 'NAME'){
        rating[i] <- 200
        
      } else {
        if (namesinp$SOURCE[i] %in% q1){
          rating[i] <- 200
        } else if(namesinp$SOURCE[i] %in% q2){
          rating[i] <- 20
        } else {
          rating[i] <- 30
        }
        
      }
      
      
    }
    return(rating)
  }
  
  
  
  
  ### Create the primary to secondary chebi dictionary
  # we have the CHEBIcomplete_secondaryCHEBI table
  
  #chebisecdict <- read.delim(paste(path,"CHEBIcomplete_secondaryCHEBI.tsv",sep=''), sep="\t",header = TRUE)
  # use the information from ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv
  inp <- read.table(paste(path,"/data_chebi_website/CHEBI_compounds_web.tsv",sep=''),sep ="\t",header = TRUE)
  
  chebisecdict <- data.frame(as.integer(inp$PARENT_ID),as.integer(gsub("CHEBI:",'',inp$ID)))
  # warning NAs introduced by coercion is normal here
  colnames(chebisecdict) <- c('CHEBI','secCHEBI')
  
  # finish the secondary to primary dictionary
  chebisecdict$CHEBI[is.na(chebisecdict$CHEBI)] <- chebisecdict$secCHEBI[is.na(chebisecdict$CHEBI)]
  
  chebisecdict <- unique(chebisecdict)
  
  
  ### Create the synonyms dictionary:
  
  ## Create the synonym table by reading in the CHEBI_names_official table, from ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/
  
  
  namesinp <- read.table(paste(path,'/data_chebi_website/CHEBI_names_web.tsv',sep=''),sep='\t',header = TRUE)
  
  
  chebisyndict <- data.frame(as.integer(namesinp$COMPOUND_ID),namesinp$NAME)
  colnames(chebisyndict) <- c('CHEBI','name')
  
  # rate the names according to:
  # If it is a NAME -> 200
  # if it is a synonyme -> look at the source: KEGG, ChEBI, MetaCyc, IUPAC = 50
  #                                            SUBMITTER = 20
  #                                            all others = 30
  
  
  chebisyndict$origin <- RateName(namesinp)
  
  chebisyndict = as.data.frame(chebisyndict)
  chebisyndict$origin = as.integer(chebisyndict$origin)
  chebisyndict$CHEBI = as.integer(chebisyndict$CHEBI)
  
  # Add the chebi names
  # Things that have no parent ID are main CHEBIs and have a main name
  tmp <- data.frame(as.integer(gsub("CHEBI:",'',inp$ID)),inp$NAME) # inp is from above
  tmp$origin <- 300 # 300= it is a main name
  # for compounds having less than 3 stars, the name should have a rather low origin
  # -> we want to minimize the amount of non 3 star compounds
  tmp$origin[as.integer(inp$STAR) < 3] <- 40
  
  tmp <- tmp[is.na(as.integer(inp$PARENT_ID)),]
  # warning NA introduced by coercion is normal here
  colnames(tmp)<- c('CHEBI','name','origin')
  
  
  
  chebisyndict <- rbind(chebisyndict,tmp)
  rm(tmp)
  rm(inp)
  # Add the additional synonymes
  
  chebisyndict <- Adddict(paste(path,'/data_manual/CHEBI_manualnames.tsv',sep=''),chebisyndict) # manual annotations
  chebisyndict <- Adddict(paste(path,'/data_manual/CHEBI_rxnparnames.tsv',sep=''),chebisyndict) # chebi name matches from the stoich model
  
  chebisyndict$origin = as.integer(chebisyndict$origin)
  chebisyndict$CHEBI = as.integer(chebisyndict$CHEBI)
  
  # sort the chbi and with respect to the origin -> the higher originss first
  chebisyndict = chebisyndict[order(chebisyndict$CHEBI,chebisyndict$origin,decreasing=TRUE), ]
  
  
  
  
  # Filter remove lines with have a negative origin and their duplicates
  
  for(i in which(chebisyndict$origin < 0)){
    tmprow <- chebisyndict[i,c(1,2)]
    rows <- which(apply(mapply(chebisyndict[,c(1,2)],tmprow, FUN="=="), MARGIN=1, FUN=all));
    chebisyndict <- chebisyndict[-rows,]
  }
  
  
  # create CHEBI main name -> take the name with the highest quality origin
  #sort according to the origin
  chebisyndict = chebisyndict[order(chebisyndict$CHEBI,chebisyndict$origin,decreasing=TRUE), ]
  
  #chebinamedict <- read.delim(paste(path,"CHEBIcomplete_CHEBIname.tsv",sep=''), sep="\t",header = TRUE)
  chebinamedict <- chebisyndict
  chebinamedict$CHEBI <- as.integer(chebinamedict$CHEBI)
  # -> finish it later
  
  
  
  
  # create a dict to unite similar CHEBIs eg. ATP, ATP(2-) etc.
  # -> use higher quality synonyms (origin >= 50) and unite all CHEBIs which share them
  #    under the name of the one synonym with the shortest CHEBIname
  
  chebisyndict$name <- tolower(gsub(' ','',chebisyndict$name))
  
  hqchebisyndict <- chebisyndict[chebisyndict$origin >= 40, ]
  
  # Fuzzy search should be stereo unspecific -> remove indicators
  # remove (-nr) at the end of words -> 
  hqchebisyndict$name <- gsub('\\b\\([0-9][-,+]\\)$','',hqchebisyndict$name)
  # remove (R)-/(S)- in front of words
  hqchebisyndict$name <- gsub('^\\([a-z]\\)\\-','',hqchebisyndict$name)
  # remove L-/R- in front of words
  hqchebisyndict$name <- gsub('^[l,d]-','',hqchebisyndict$name)
  #remove  zwitterion end
  hqchebisyndict$name <- gsub('zwitterion$','',hqchebisyndict$name)
  #remove acid at end
  hqchebisyndict$name <- gsub('\\Bicacid$','ate',hqchebisyndict$name)
  # Remove alpha-D-
  hqchebisyndict$name <- gsub('^(beta|alpha)-[l,d]-','',hqchebisyndict$name)
  hqchebisyndict$name <- gsub('-(beta|alpha)-[l,d]-','',hqchebisyndict$name)
  # Remove N(0-9)-acetyl -> N-acetyl
  hqchebisyndict$name <- gsub('^[n,o]\\(.\\)-','n-',hqchebisyndict$name)
  # Remove trans/cis
  # from start
  hqchebisyndict$name <- gsub('^(trans|cis)-','',hqchebisyndict$name)
  # from the middle 
  hqchebisyndict$name <- gsub('-(trans|cis)-','',hqchebisyndict$name)
  # Remove -L-/-D-
  hqchebisyndict$name <- gsub('-[l,d]-','',hqchebisyndict$name)
  # Remove -(S/R)-
  hqchebisyndict$name <- gsub('-\\([l,d]\\)\\-','',hqchebisyndict$name)
  
  fil <- rowSums(chebisyndict[chebisyndict$origin >= 40, ] == hqchebisyndict) < 3
  fil[is.na(fil)] <-1
  hqchebisyndict$origin[fil] <- hqchebisyndict$origin[fil]-1
  
  chebisyndict <- unique(rbind(chebisyndict,hqchebisyndict))
  # -> I use this later
  
  
  ### load the chemical formulas
  
  cheminp <- read.table(paste(path,'/data_chebi_website/CHEBI_chemical_data_web.tsv',sep=''),sep='\t',header = TRUE)
  cheminp <- cheminp[cheminp$TYPE == 'FORMULA', ]
  
  cheminp <- cheminp[order(cheminp$COMPOUND_ID,cheminp$SOURCE), ]
  
  cheminp <- cheminp[!duplicated(cheminp$COMPOUND_ID), ]
  
  cheminp <- data.frame(cheminp$COMPOUND_ID,cheminp$CHEMICAL_DATA)
  colnames(cheminp) <- c('CHEBI','formula')
  
  cheminp <- cheminp[!cheminp$formula == '.',]
  cheminp <- unique(cheminp)
  
  
  # add 1000 to the origin of all compounds with any formula
  # -> if the names are ambiguous, we want to prefer stuff with formulas, as these are higher quality compounds
  
  fil <- chebisyndict$CHEBI %in% cheminp$CHEBI
  chebisyndict$origin[fil] <- chebisyndict$origin[fil] + 1000
  
  # remove everything not alphanumeric (expt +) and spaces from the names
  
  chebisyndict$name <- tolower(as.character(chebisyndict$name))
  chebisyndict$name <- gsub("[^a-zA-Z0-9+]","",chebisyndict$name)

  
  
  # sort the names alphabetically and with respect to the origin -> the higher origins first
  
  chebisyndict = chebisyndict[order(chebisyndict$name,chebisyndict$origin,decreasing=TRUE), ]
  
  
  
  
  ##### Proceed with the syndict remove all duplicates, keep the one with the highest origin
  
  chebisyndict = chebisyndict[!duplicated(chebisyndict$name), ]
  
  # the creation of the dictionary is completed
  #clean up:
  rm(Adddict)
  rm(RateName)
  rm(namesinp)
  
  ### Add all not synonym CHEBIs that have no entry as secondary chebi as primary chebis 
  
  tmp <- unique(chebisyndict$CHEBI[!chebisyndict$CHEBI %in% chebisecdict$secCHEBI])
  tmp <- cbind(tmp,tmp)
  colnames(tmp) <- c('CHEBI','secCHEBI')
  chebisecdict <- rbind(chebisecdict,tmp)
  chebisecdict <- unique(chebisecdict)
  
  rm(tmp)
  
  
  ### Create the chebikeggdict
  chebikeggdict <- read.table(paste(path,"/data_from_sdf/CHEBI3star_KeggID.tsv",sep=''), sep="\t",header = TRUE)
  chebikeggdict <- rbind(chebikeggdict,read.table(paste(path,"/data_from_sdf/CHEBI3star_KeggID_chebi.tsv",sep=''), sep="\t",header = TRUE))
  chebikeggdict <- rbind(chebikeggdict,read.table(paste(path,"/data_from_sdf/CHEBIcomplete_KeggID_chebi.tsv",sep=''), sep="\t",header = TRUE))
  chebikeggdict <- rbind(chebikeggdict,read.table(paste(path,"/data_from_sdf/CHEBIcomplete_KeggID.tsv",sep=''), sep="\t",header = TRUE))
  chebikeggdict <- rbind(chebikeggdict,read.table(paste(path,"/data_manual/CHEBImanual_KeggID.tsv",sep=''), sep="\t",header = TRUE))
  
  #chebikeggdict <- rbind(chebikeggdict,read.table(paste(path,"/data_from_sdf/CHEBIcomplete_KeggID.tsv",sep=''), sep="\t",header = TRUE))
  
  
  
  colnames(chebikeggdict) <- c('CHEBI','KEGG')
  chebikeggdict$KEGG <- gsub('cpd:','',chebikeggdict$KEGG)
  chebikeggdict$CHEBI <- as.integer(chebikeggdict$CHEBI)
  
  ### Create the CHEBI fuzzy synonyme dictionaire -> combine all chebi that have at least
  # one high quality synonyme in common
  
  
  
  hqchebisyndict <- hqchebisyndict[,-3]
  hqchebisyndict <- unique(hqchebisyndict)
    
  
  # remove everything non alphanummerical (probably it could be a good idea to keep the '+')
  hqchebisyndict$name <- tolower(gsub("[^a-zA-Z0-9]","",hqchebisyndict$name))
  
  hqchebisyndict <- unique(hqchebisyndict)
  
  
  
  
  
  fil <- duplicated(hqchebisyndict$name)
  
  uninames <- unique(hqchebisyndict$name[fil])
  
  # make the dictionary smaller
  hqchebisyndict <- hqchebisyndict[hqchebisyndict$name %in% uninames, ]
  # translate all to primary chebi
  hqchebisyndict$CHEBI[!hqchebisyndict$CHEBI %in% chebisecdict$CHEBI] <- sapply(hqchebisyndict$CHEBI[!hqchebisyndict$CHEBI %in% chebisecdict$CHEBI],Chebi2pchebi)
  hqchebisyndict <- unique(hqchebisyndict)
  
  # Add the chemical CHNO formula
  
  cheminp$CHEBI <- sapply(as.integer(cheminp$CHEBI),Chebi2pchebi)
  # sort it by chebi and formula length -> we prefer the shortest formula
  cheminp <- cheminp[order(cheminp$CHEBI,sapply(cheminp$formula,nchar)),]
  cheminp <- cheminp[!duplicated(cheminp$CHEBI), ]
  # remove everything but CNO from the formulas
  # remove all two letter compounds eg Na, Cl etc
  cheminp$formula <- gsub('[A-Z][a-z][0-9]?[0-9]?[0-9]?','',cheminp$formula)
  # remove everything but N C 
  cheminp$formula <- gsub('[^NCO0-9][0-9]?[0-9]?[0-9]?','',cheminp$formula)
  
  
  #add the formula to all names -> only compounds with same name and fromula are considered equal
  fil <- hqchebisyndict$CHEBI %in% cheminp$CHEBI
  
  for (i in which(fil)){
    hqchebisyndict$name[i] <- paste(hqchebisyndict$name[i],cheminp$formula[cheminp$CHEBI == hqchebisyndict$CHEBI[i]],sep='')
  }
  
  rm(cheminp,fil)
  
  # remove 1 letter words (otherwise they are problems with AA one letter codes)
  hqchebisyndict <- hqchebisyndict[nchar(hqchebisyndict$name) > 1, ]
  
  fil <- duplicated(hqchebisyndict$name)
  
  uninames <- unique(hqchebisyndict$name[fil])
  nameidx <- setNames(1:length(uninames),uninames)
  
  # make the dictionary smaller again
  hqchebisyndict <- hqchebisyndict[hqchebisyndict$name %in% uninames, ]
  
  
  
  unichebi <- unique(hqchebisyndict$CHEBI[hqchebisyndict$name %in% uninames])
  chebiidx <- setNames(1:length(unichebi),unichebi)
  
  # create a matrix uniname_idx vs unichebi_idx 
  
  tmp <- nameidx[hqchebisyndict$name[hqchebisyndict$name %in% uninames]]
  pairs <- matrix(c(tmp,chebiidx[as.character(hqchebisyndict$CHEBI[hqchebisyndict$name %in% uninames])]),ncol=2)
  
  #pairs <- sparseMatrix(pairmat[ ,1],pairmat[ ,2]) # rows=name, col= chebi
  
  pairmat <- Matrix(FALSE,nrow = length(uninames),ncol=length(unichebi),sparse = TRUE)
  
  pairs <- unique(pairs)
  pairmat[pairs] <- TRUE
  
  rm(pairs)
  
  namemat <- pairmat %*%  t(pairmat)  # gives the squared chebi chebi matrix, 1 whenever 2 chebis are in one group 
  
  rm(pairmat)
  
  #i=TRUE
  #while(i){
  for (i in 1:20){
    namemat <- namemat %*% namemat
  }
  #oldnamemat <- namemat != 0
  #namemat <- namemat %*% namemat
  namemat <- namemat != 0
  
  # if the mats do not change anymore we have reached the point where each group
  # will have no element in common with another group, when merging all groups that have a 1 in a row
  # i= !all(namemat==oldnamemat)
  #} -> I dont do this check as it uses a lot of memory and all groups should be merged after 20 iterations
  
  # We have no the namemat, with the dimension group x group.
  # Each groups are now merged rowwise: -> in the end
  
  vec <- 1:length(uninames)
  chebifuzsyndict <- matrix(nrow=0,ncol=2) # this will have the fuzzy parent vs the daughter chebis
  tmprow <-  matrix(nrow=1,ncol=2)
  while (length(vec) > 0){
    i=vec[1]
    idx <- which(namemat[,i])
    syns <- names(nameidx[idx])
    chebi <- hqchebisyndict$CHEBI[hqchebisyndict$name %in% syns]
    chebi <- unique(chebi)
    pname <- sapply(chebi,Chebi2name)
    pchebi <- chebi[which.min(nchar(pname))]
    tmprow[1] <- pchebi
    for (k in 1:length(chebi)){
      tmprow[2] <- chebi[k]
      chebifuzsyndict <- rbind(chebifuzsyndict,tmprow)
    }
    vec <- setdiff(vec,idx)
  }
  
  rm(hqchebisyndict)
  
  chebifuzsyndict <- data.frame(as.integer(chebifuzsyndict[,1]),as.integer(chebifuzsyndict[,2]))
  colnames(chebifuzsyndict) <- c('pfCHEBI','secfCHEBI')
  
  ## Make sure that allways primary chebis are used as identifiers in the list
  chebikeggdict$CHEBI[!chebikeggdict$CHEBI %in% chebisecdict$CHEBI] <- sapply(chebikeggdict$CHEBI[!chebikeggdict$CHEBI %in% chebisecdict$CHEBI],Chebi2pchebi)
  chebisyndict$CHEBI[!chebisyndict$CHEBI %in% chebisecdict$CHEBI] <- sapply(chebisyndict$CHEBI[!chebisyndict$CHEBI %in% chebisecdict$CHEBI],Chebi2pchebi)
  chebinamedict$CHEBI[!chebinamedict$CHEBI %in% chebisecdict$CHEBI] <- sapply(chebinamedict$CHEBI[!chebinamedict$CHEBI %in% chebisecdict$CHEBI],Chebi2pchebi)
  
  # make again sure everything is unique and rows with nan only are removed
  chebikeggdict <- unique(chebikeggdict[rowSums(is.na(chebikeggdict)) != ncol(chebikeggdict),])
  chebisyndict <- unique(chebisyndict[rowSums(is.na(chebisyndict)) != ncol(chebisyndict),])
  chebinamedict <-unique(chebinamedict[rowSums(is.na(chebinamedict)) != ncol(chebinamedict),])
  
  ## finish the chebinamedict
  # remove all synonyms
  #sort according to the origin
  chebinamedict = chebinamedict[order(chebinamedict$CHEBI,chebinamedict$origin,decreasing=TRUE), ]
  # remove duplicated
  chebinamedict <- chebinamedict[!duplicated(chebinamedict$CHEBI), ]
  
  
 
 # clean up
 rm(vec,pchebi,tmprow,chebi,chebiidx,fil,i,idx,k,nameidx,namemat,pname,syns,tmp,unichebi,uninames)
 files <- c("Chebi2fuzchebi", "Chebi2kegg" ,  "Chebi2name" ,  "Chebi2pchebi" , "Chebi2shortname",
 "Chebi2synkegg" , "chebifuzsyndict", "chebikeggdict",  "chebinamedict" , "chebisecdict",  
 "chebisyndict" , "GetClosestMatch" ,"Kegg2chebi"   ,"MatchName2Chebi" ,"Name2chebi",   
  "Querychebi2kegg" ,"Querykegg2chebi" )
 save(list=files, file=paste(path,'/cache/find_chebiData.Rdata',sep=''))
}




# check in the synonym list, wheter there are CHEBI with no main name and annotate that

##### varia: read chebis names from the reactionparfile and save them as manual annotated names

# rxnparchebi = matrix(nrow = 0, ncol=2)
# colnames(rxnparchebi) <- c('CHEBI','name')
# tmp = matrix(nrow = 1, ncol=2)
# for(i in 1:length(rxnparFile$V1)){
#   if(length(grep('chebi',rxnparFile[i, 3])) > 0){
#     tmp[1,2] <- rxnparFile[i,2]
#     tmp[1,1] <-strsplit(rxnparFile[i,3], split = "%3A")[[1]][2]
#     rxnparchebi <- rbind(rxnparchebi,tmp)
#   }
# }
# 
# write.table(rxnparchebi,file = paste(path,'CHEBI_rxnparnames.tsv',sep=''),
#             row.names = FALSE, quote=FALSE,sep = "\t")

