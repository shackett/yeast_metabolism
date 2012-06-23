"""
Take a list of yeast ChEBI IDs - determine corresponding possible KEGG ids
Compare the list of n sets of ChEBI-associated KEGG ids with a list of ecoli KEGG IDs
and return which KEGG ids match
"""

#qsub -l 1day -cwd -sync n python ecoliYeastMatch.py

yeast_names = 'yeast_abs_chebi_id.csv'
out_file = 'yeast_abs_ids.dict.csv'

import ../EcoliYeastMatch/CTSpy
import re


"""
Import yeast names, associate name with synonyms and then synonyms with chebi 
and kegg ids
"""

for l in open(ecoli_kegg):
	line = l.strip().split(",")
	
	
quit()
sys.exit()	

synNametoKEGG = dict()
synNametoChEBI = dict()

cts = CTSpy.CTS()

for id in chebi_ids:
	import_pass = False
	while import_pass == False:
	 	try:
			directCHtoKEGG[id] = cts.convert("chebi", "kegg", id)
			import_pass = True
		except:
			import_pass = False
	
	import_pass = False
	syn_set = set()
	while import_pass == False:
		try:
			ch_syn = cts.convert("chebi", "names", id)
			import_pass = True
		except:
			import_pass = False
			
	for syn in ch_syn:
		import_pass = False
		while import_pass == False:
			try:
				keggids = cts.convert("name", "kegg", syn.split(" - ")[0])
				import_pass = True
			except:
				import_pass = False
		for aKegg in keggids:
			syn_set.add(aKegg)
	synCHtoKEGG[id] = syn_set


"""
Associate yeast metabolite IDs with ecoli metabolite IDs via id-chebi-kegg-kegg-id
"""

yeast_ecoli_id_dict = dict()

for yeastID in yeast_id_dict.keys():
	corr_ids = set()
	yeast_ecoli_id_dict[yeastID] = corr_ids
	
	for eID in synCHtoKEGG[yeast_id_dict[yeastID]]:
		if eID in unique_ecoli_kegg.keys():
			for el in unique_ecoli_kegg[eID]:
				yeast_ecoli_id_dict[yeastID].add(el)
		

outf = open(out_file, "w")		
outf.write("yeastID\tecoliID" + "\n")
			
for l in yeast_ecoli_id_dict.keys():
	outf.write(l + "\t" + ",".join(yeast_ecoli_id_dict[l]) + "\n")
outf.close()	 
	
		
