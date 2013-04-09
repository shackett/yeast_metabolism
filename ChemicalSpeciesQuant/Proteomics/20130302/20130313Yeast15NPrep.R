#setwd("/home/jgoya/Dropbox/ChemostatYeastProteomics/20130302/")

setwd("~/Desktop/Rabinowitz/FBA_SRH/ChemicalSpeciesQuant/Proteomics/20130302")

rbind(read.delim("20130313Yeast15NRound1.qtsv"),
      read.delim("20130313Yeast15NRound2.qtsv"))->QTSV
QTSV$Iter<-as.factor(c(rep(0, 1026710), rep(1, 590122), rep(2, 400606), rep(3, 307210), 
                       rep(0, 190324), rep(1, 586134), rep(2, 233278), rep(3, 104250)))
QTSV$Round<-as.factor(c(rep(1, 2324648), rep(2, 1113986)))


PreparePeptides<-function(Peptides) {
  
  Peptides$File<-gsub("(.+)\\.\\d{5}\\.\\d{5}\\.\\d", "\\1", Peptides$SpectrumID, perl=T)
  as.factor(Peptides$File)->Peptides$File
  
  Peptides$Charge<-gsub("(.+)\\.\\d{5}\\.\\d{5}\\.(\\d)", "\\2", Peptides$SpectrumID, perl=T)
  
  
  Peptides$PepID<-paste(Peptides$Sequence, Peptides$Charge, sep=".")
  as.factor(Peptides$PepID)->Peptides$PepID
  
  return(Peptides)
}

FilterPeptides<-function(Peptides, minCount) {
  ## get unique File:PepID pairs
  unique(Peptides[c("File", "PepID")])->PepFileClean
  
  ## make boolean matrix of PepID and File paris
  PepFileCount<-colSums(xtabs( ~ File + PepID, PepFileClean))
  
  return(droplevels(Peptides[ Peptides$PepID %in% names(PepFileCount[ PepFileCount>=minCount]),]))
}

QTSV<-PreparePeptides(QTSV)
QTSV$QuantID<-paste(QTSV$Iter, QTSV$File, QTSV$PepIndex, sep=".")
QTSV$SN<-QTSV$Quant/QTSV$SD
#QTSV<-droplevels(QTSV[ QTSV$QuantID %in% unique(QTSV[QTSV$SN>10,c("QuantID")]),])
#QTSV<-FilterPeptides(QTSV, 2)

library(reshape2)
QTSV_SD<-dcast(QTSV, QuantID ~ QuantClass, value.var="SD")
QTSV_SD$SD<-sqrt(QTSV_SD$N14^2 + QTSV_SD$N15^2)
QTSV_SD<-QTSV_SD[c("QuantID", "SD")]

QTSV_SN<-dcast(QTSV, QuantID ~ QuantClass, value.var="SN")
QTSV_SN$SN<-sqrt(QTSV_SD$N14^2 + QTSV_SD$N15^2)
QTSV_SN<-QTSV_SD[c("QuantID", "SD")]

QTSV_Ratios<-dcast(QTSV, QuantID ~ QuantClass, value.var="Quant")
QTSV_Ratios$Quant<-QTSV_Ratios$N14 + QTSV_Ratios$N15
QTSV_Ratios<-merge(QTSV_Ratios, QTSV_SD)
QTSV_Ratios$Log2LvH<-log2(QTSV_Ratios$N14/QTSV_Ratios$N15)
QTSV_Ratios<-merge(QTSV_Ratios, unique(QTSV[c("QuantID", "SpectrumID", "Sequence", "ProtAcc", "Iter")]))

  
QTSV_Ratios$File<-as.factor(gsub("(.+)\\.\\d+\\.\\d+\\.\\d", "\\1", QTSV_Ratios$SpectrumID, perl=T))
QTSV_Ratios$Set<-gsub(".+_([CLNPU]0\\.\\d{2})_.+", "\\1", QTSV_Ratios$File)
QTSV_Ratios$Set<-gsub("P030_\\d{2}", "P0.30", QTSV_Ratios$Set)
QTSV_Ratios<-QTSV_Ratios[grep("N14N15", QTSV_Ratios$Set, invert=T),]
as.factor(QTSV_Ratios$Set)->QTSV_Ratios$Set
  
QTSV_Ratios$Charge<-gsub("(.+)\\.\\d+\\.\\d+\\.(\\d)", "\\2", QTSV_Ratios$SpectrumID, perl=T)
  
QTSV_Ratios$PepID<-paste(QTSV_Ratios$Sequence, QTSV_Ratios$Charge, sep=".")
as.factor(QTSV_Ratios$PepID)->QTSV_Ratios$PepID

QTSV_Ratios$PepID.File<-paste(QTSV_Ratios$PepID, QTSV_Ratios$File, sep=".")

unique(QTSV_Ratios[c("PepID", "File")])->PepID.File
library(data.table)

QTSV_Ratios<-data.table(QTSV_Ratios)
setkeyv(QTSV_Ratios, c("PepID.File"))
system.time({QTSV_Ratios.Uni<-QTSV_Ratios[,data.table(QuantID=.SD[which.max(.SD$N14+.SD$N15)]$QuantID),by="PepID.File"]
            setkeyv(QTSV_Ratios.Uni, c("PepID.File", "QuantID"))
            setkeyv(QTSV_Ratios, c("PepID.File", "QuantID"))
            QTSV_Ratios.Uni<-QTSV_Ratios[QTSV_Ratios.Uni]})

QTSV_Ratios.Uni$Rep<-ifelse(grepl("OG-", QTSV_Ratios.Uni$File), "QTOF6538", 
                            ifelse(grepl("OGfr", QTSV_Ratios.Uni$File), "QTOF6550-1", 
                                   ifelse(grepl("2Factor", QTSV_Ratios.Uni$File), "QTOF6550-2", 
                                          ifelse(grepl("P030", QTSV_Ratios.Uni$File), "ORBI", "ERROR"))))

QTSV_Ratios.Uni<-data.table(QTSV_Ratios.Uni)
setkeyv(QTSV_Ratios.Uni, c("PepID", "Set", "Rep"))
system.time({QTSV_Ratios.Uni.Set<-QTSV_Ratios.Uni[,data.frame(LQuant=sum(.SD$N14), HQuant=sum(.SD$N15)),by=c("PepID", "Set", "Rep")]})

save(QTSV_Ratios.Uni.Set, file="20130313QTSV_RatiosUniSet.Rdata")

## Remove methionine oxidations?
#QTSV_Ratios.Uni.Set<-QTSV_Ratios.Uni.Set[!grepl("m", QTSV_Ratios.Uni.Set$PepID),]

load("20130313QTSV_RatiosUniSet.Rdata")
QTSV_Ratios.Uni.Set$Log2LvH<-log2(QTSV_Ratios.Uni.Set$LQuant/QTSV_Ratios.Uni.Set$HQuant)
library(reshape2)
acast(QTSV_Ratios.Uni.Set, PepID~Set + Rep, value.var="Log2LvH")->PepMatrix
acast(QTSV_Ratios.Uni.Set, PepID~Set + Rep, value.var="HQuant")->heavyIC
acast(QTSV_Ratios.Uni.Set, PepID~Set + Rep, value.var="LQuant")->lightIC

#PepMatrix[rowSums(is.finite(PepMatrix))==16,]->test2
#test2[sample(rownames(test2)),]->test2

#pdf("20120806QTOF_HQuant.pdf", height=8.5, width=11)
#l_ply(1:228, function (i) {
#print(ggplot(melt(test2[((i-1)*10+1):(i*10),], value.name="HQuant", varnames=c("PeptideID", "Sample"))) + aes(x=Sample, y=log10(HQuant), group=PeptideID) + geom_line(aes(col=PeptideID)))
#})
#dev.off()

library(ggplot2)
ggplot(QTSV_Ratios.Uni.Set) + facet_grid(Rep~Set) + geom_histogram(aes(x=Log2LvH)) + xlim(-5,5)


rbind(read.delim("SpecRep20130106YeastMerge_Condensed.xls", 
                 skip=(grep("Experiment name", readLines("SpecRep20130106YeastMerge_Condensed.xls", n=1000))-1)),
      read.delim("SpecRep20130131Yeast2FactorRound2_10mDa.xls", 
                 skip=(grep("Experiment name", readLines("SpecRep20130131Yeast2FactorRound2_10mDa.xls", n=1000))-1)))->SpecRep

SpecRep$PepID<-as.factor(paste(SpecRep$Peptide.sequence, SpecRep$Spectrum.charge, sep="."))
unique(SpecRep[c("PepID", "Protein.accession.numbers", "Other.Proteins")])->ProtPep
ProtPep<-droplevels(ProtPep[ ProtPep$PepID %in% levels(droplevels(QTSV_Ratios.Uni.Set$PepID)),])

library(plyr)
ddply(ProtPep, .(PepID), function(df) {
  unique(unlist(c(strsplit(as.character(df$Protein.accession.numbers), split=","), strsplit(as.character(df$Other.Proteins), split=","))))->ProtAcc
  data.frame(ProtAcc=ProtAcc)
})->ProtPep
ProtPep$Val<-1
acast(ProtPep, PepID~ProtAcc, value.var="Val", fill=0)->ProtPepMatrix

save(PepMatrix, ProtPepMatrix, heavyIC, lightIC, file="20130313ProtPepMatrices.Rdata")
