"""
Take a list of yeast ChEBI IDs - determine corresponding possible KEGG ids
Compare the list of n sets of ChEBI-associated KEGG ids with a list of ecoli KEGG IDs
and return which KEGG ids match
"""


#qsub -l 1day -cwd -sync n python ecoliYeastMatch.py

yeast_names = 'yifan_abs_names.tsv'
out_file = 'yeast_abs_ids.dict.tsv'

import CTSpy
import re


"""
Import yeast names, associate name with synonyms and then synonyms with chebi 
and kegg ids
"""
	
synNametoKEGG = dict()
synNametoChEBI = dict()

cts = CTSpy.CTS()

l_count = 1
for id in open(yeast_names):
    id = id[:-1].split("\t")[0]
    
    ch_syn = id
    
    syn_set_kegg = set()
    syn_set_chebi = set()
    
    for syn in ch_syn:
        
        import_pass = False
        while import_pass == False:
            try:
                keggids = cts.convert("name", "kegg", syn.split(" - ")[0])
                chebiids = cts.convert("name", "chebi", syn.split(" - ")[0])
                import_pass = True
            except:
                import_pass = False
        
        for aKegg in keggids:
            syn_set_kegg.add(aKegg)
        for aChebi in chebiids:
            if re.match("CHEBI:", aChebi):
                ch_rem = re.search("CHEBI:", aChebi)
                aChebi = aChebi[:ch_rem.start()] + aChebi[ch_rem.end():]
            syn_set_chebi.add(aChebi)
            
    synNametoKEGG[id] = syn_set_kegg
    synNametoChEBI[id] = syn_set_chebi
    
    #l_count += 1
    #if l_count > 2:
     #   break
				

outf = open(out_file, "w")		
outf.write("Name\tChEBI ID\tKEGG ID" + "\n")

for id in open(yeast_names):
	id = id[:-1].split("\t")[0]
	print id
	outf.write(id + "\t" + ",".join(synNametoChEBI[id]) + "\t" + ",".join(synNametoKEGG[id]) + "\n")	
outf.close()	 
	
		
