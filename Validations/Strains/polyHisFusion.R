setwd("~/Desktop/Rabinowitz/FBA_SRH/Validations/Strains")

library(seqinr) #read.fasta / translate


YFGs <- c("YJL130C", "YJR109C")

# load each ORF with 1000 bp of 5' and 3' flanking

ORFs_extended <- read.fasta("../../Yeast_reconstruction/Sequences/orf_genomic_1000.tsv")
ORFs <- read.fasta("../../Yeast_reconstruction/Sequences/orf_genomic.fasta")

if(!all(YFGs %in% names(ORFs) & YFGs %in% names(ORFs_extended))){
  stop("Some genes not recognized")
}

# ordered restriction sites that are available - the first two RE without a site in the ORF will be chosen
restrictionSites <- data.frame(RE = c("SpeI", "XmaI", "XhoI"), site = c("actagt", "cccggg", "ctcgag"))

primerOutput <- c("6xHis extension primers")
primerWarnings <- NULL
YFG_RE <- data.frame(ORF = YFGs, upstreamRE = NA, downstreamRE = NA)

for(YFG in YFGs){
  
  # check ORF for restriction sites
  # only check coding strand because REs of interest recognize reverse pallindrome
  concatenatedORF <- paste(ORFs[[YFG]], collapse = "")
  
  for(RE in restrictionSites$RE){
    if(gregexpr(restrictionSites$site[restrictionSites$RE == RE], concatenatedORF)[[1]][1] != -1){
      warning(paste(length(gregexpr(restrictionSites$site[restrictionSites$RE == RE], concatenatedORF)[[1]]), RE, "site(s) in", YFG))
    }else{
      YFG_RE[YFG_RE$ORF == YFG, c(2:3)][is.na(YFG_RE[YFG_RE$ORF == YFG, c(2:3)])][1] <- RE
    }
    if(all(!is.na(YFG_RE[YFG_RE$ORF == YFG, c(2:3)]))){
      next
    }
  }
  
  if(any(is.na(YFG_RE[YFG_RE$ORF == YFG, c(2:3)]))){
    warning(paste(sum(is.na(YFG_RE[YFG_RE$ORF == YFG, c(2:3)])), "more valid restriction enzyme(s) needed for", YFG))
    primerWarnings <- c(primerWarnings, paste(sum(is.na(YFG_RE[YFG_RE$ORF == YFG, c(2:3)])), "more valid restriction enzyme(s) needed for", YFG))
    next
  }
  
  # In frame 3' fusion
  # -Stop
  # 6xHis (CAU, CAC)
  # + Stop (TGA)
  # downstream restriction site
  # an extra set of trinucleotides to aid in RE docking
  
  terminalTag <- paste(c("catcaccatcaccatcattga", restrictionSites$site[restrictionSites$RE == YFG_RE$downstreamRE[YFG_RE$ORF == YFG]], "gtc"), collapse = "")
  seqinr::translate(s2c(terminalTag)) # verify
  
  translated_ORF <- seqinr::translate(ORFs[[YFG]]) # for fun
  
  # generate forward primer
  # five prime GAC to aid in RE docking
  # upstream restriction site
  # stretch of homology is 20+ bp ending in a C or G
  
  ORFfivePrime <- c(s2c("gac"), s2c(restrictionSites$site[restrictionSites$RE == YFG_RE$upstreamRE[YFG_RE$ORF == YFG]]), ORFs[[YFG]][c(1:(19 + which(ORFs[[YFG]][-c(1:19)] %in% c("c", "g"))[1]))])
  primerOutput <- c(primerOutput, "", YFG, paste(YFG_RE$upstreamRE[YFG_RE$ORF == YFG], "- Forward:"), paste("5\' ", toupper(paste(ORFfivePrime, collapse = "")), " 3\'", sep = ""))
  
  
  # remove stop codon
  # find the 3' reverse complement
  
  ORFthreePrime <- rev(rev(ORFs[[YFG]])[-c(1:3)][1:(19 + which(rev(ORFs[[YFG]])[-c(1:3)][-c(1:19)] %in% c("c", "g"))[1])])
  ORFfusion <- c(ORFthreePrime, s2c(terminalTag))
  if(length(ORFfusion) %% 3 == 0){
    seqinr::translate(ORFfusion)
  }else{
    seqinr::translate(ORFfusion[-c(1:(length(ORFfusion) %% 3))])
  }
  
  # generate reverse primer
  primerOutput <- c(primerOutput, YFG, paste(YFG_RE$downstreamRE[YFG_RE$ORF == YFG], "- Reverse:"), paste("5\' ", toupper(paste(rev(comp(ORFfusion)), collapse = "")), " 3\'", sep = ""))
  
}

if(length(primerWarnings) != 0){
  primerOutput <- c(primerOutput, "", "Warings:", primerWarnings)
}


write.table(primerOutput, "FusionPrimers.tsv", quote = F, row.names = F, col.names = F)



