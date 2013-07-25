import sys
import itertools
import re

KEGG = open('keggDict.tsv','r')
organism = 'SCE:'


"""
for each KEGG reaction, match the KEGG ID to the EC number and match the each yeast gene to the KEGG ID
output yeast systematic name \t KEGG ID \t EC number
"""

KEGGtoEC = dict()
KEGGtoRxn = dict()
SYSTtoKEGG = dict()
SYSTtoSGD = dict()

RID = re.compile('RN:')
 

#count = 0
for ln in KEGG:
	ln = re.sub('GENES', '', ln)
	ln = re.sub(' +', ' ', ln)
	
	ln_split = ln.split(' ')
	
	if len(ln_split) == 1:
		continue
	
	if ln_split[0] == "ENTRY":
		keggID = ln_split[1]
		KEGGtoEC[keggID] = 'NA'
		KEGGtoRxn[keggID] = 'NA'
		
	if RID.search(ln):
	
		#print re.sub('\n', '', ln.split('RN:')[1])[1]
		#print re.sub('\n', '', "".join(ln.split('RN:')[1]))
		KEGGtoRxn[keggID] = ''.join(re.sub('\n|', '', ln.split('RN:')[1]))
		
	
	if ln_split[0] == "DEFINITION":
		KEGGtoEC[keggID] = re.sub('\n', '', " ".join(ln_split[1:]))
	
	if ln_split[1] == organism:
		if keggID == "K07019":
			print ln_split
		for orgID in ln_split[2:]:
			geneID = re.sub('\n', '', orgID)
			geneIDsplit = geneID.split('(')
			
			if len(geneIDsplit) == 1:
				SYSTtoSGD[geneIDsplit[0]] = 'NA'
			else:
				SYSTtoSGD[geneIDsplit[0]] = re.sub('\)', '', geneIDsplit[1])
			SYSTtoKEGG[geneIDsplit[0]] = keggID
		
		
	#count += 1
	#if count > 10000:
	#	break

outf = open('yeastNameDict.tsv', 'w')
outf.write("SYST\tSGD\tKEGG\tNAME\tRID\n")
	
for gene in SYSTtoKEGG.keys():
	outf.write("%s\t%s\t%s\t%s\t%s\n" % (
	gene, SYSTtoSGD[gene], SYSTtoKEGG[gene], KEGGtoEC[SYSTtoKEGG[gene]], KEGGtoRxn[SYSTtoKEGG[gene]]
	))
outf.close()



