"""
Take a list of yeast ChEBI IDs - determine corresponding possible KEGG ids
Compare the list of n sets of ChEBI-associated KEGG ids with a list of ecoli KEGG IDs
and return which KEGG ids match
"""

ecoli_kegg = 'ecoli_kegg_id.csv'
yeast_chebi = 'yeast_chebi_id.csv'

import CTSpy
import re

""" import ecoli KEGG ids, associate internal name with KEGG id """

""" 
Invert the ecoli dictionary so KEGG ids point to a set of the 
corresponding species names 
"""

ecoli_id_dict = dict()
unique_ecoli_kegg = dict()
lcount = 1
for l in open(ecoli_kegg):
	line = l.strip().split(",")
	ecoli_id_dict[line[0]] = line[1]
	
	if line[1] in unique_ecoli_kegg.keys():
		print 'woot'
	else:
		new_el = set()
		print line[0]
		unique_ecoli_kegg[line[1]] = new_el.add(line[0])
		print unique_ecoli_kegg[line[1]]
		print 'wee'
	
	
	
	lcount += 1
	if lcount > 50:
		break



print ecoli_id_dict.values()


quit()

""" import yeast ChEBI ids, associate internal name with ChEBI id """

yeast_id_dict = dict()

chebi_ids = set()
lcount = 1
for l in open(yeast_chebi):
	line = l[:-1].split(",")
	chebi_ids.add(line[1])
	
	lcount += 1
	if lcount > 5:
		break

directCHtoKEGG = dict()
synCHtoKEGG = dict()

cts = CTSpy.CTS()

for id in chebi_ids:
	directCHtoKEGG[id] = cts.convert("chebi", "kegg", id)
	
	syn_set = set()
	ch_syn = cts.convert("chebi", "names", id)
	for syn in ch_syn:
		keggids = cts.convert("name", "kegg", syn.split(" - ")[0])
		for aKegg in keggids:
			syn_set.add(aKegg)
	synCHtoKEGG[id] = syn_set

print synCHtoKEGG.keys()
print synCHtoKEGG.values()
	