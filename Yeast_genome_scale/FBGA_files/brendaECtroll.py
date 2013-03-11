""" Using a newline delimited list of E.C. numbers and a specified organism: pull down kinetic information  """

""" pip install -e "git+http://github.com/pelletier/SOAPpy.git@develop#egg=SOAPpy """

import string
import re
from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file

def removeNonAscii(s): return "".join(filter(lambda x: ord(x)<128, s)) #remove non Ascii charcters - the angstrom symbol

ec_file = 'ECnr.txt' #a text file with E.C. numbers to be queried against BRENDA
#ec_file = 'ECnr_troublemakerRemoved.txt' #a text file with E.C. numbers to be queried against BRENDA
#ec_file = 'ECnr_red.txt' #use a more reduced set for 
organism = 'Saccharomyces cerevisiae'
use_organism = 'FALSE' # either look for a specific organism (listed in the organism field) or find constants for all organims
params_to_generate = set(["km", "ki", "activators"]) #which kinetic constants should be searched for

if use_organism == 'TRUE':
	outfile = "_".join((organism.split(' ')[0].strip()[0], organism.split(' ')[1]))
else:
	outfile = "all"
 
endpointURL = "http://www.brenda-enzymes.info/soap2/brenda_server.php"
client = SOAPProxy(endpointURL)

observedLigands = set()

if "km" in params_to_generate:

	""" Determine the Km values for all metabolites associated with an E.C. number """

	km_out = open("".join((outfile, "_km.tsv")), "w")
	km_out.write("ecNumber\tkmValue\tkmValueMaximum\tsubstrate\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

	for l in open(ec_file):
		if use_organism == 'TRUE':
			resultString = client.getKmValue("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
		else:
			resultString = client.getKmValue("".join(("ecNumber*", re.split('\n',  l)[0]))).split('!')
	 
		for entry in resultString: #split each Brenda Km entry into a separate line
		
		
			all_fields = entry.split('#')
			all_fields.remove('')
		
			if len(all_fields) != 8:
				continue
			observedLigands.add(all_fields[7].split('*')[1]) #add ligandID to set of ligand IDs
		
			entry_write = []
			for field in all_fields:
				entry_write.append(field.split('*')[1])
			entry_write[4] = removeNonAscii(entry_write[4]) #some description include non ascii characters
		
			km_out.write("\t".join(entry_write) + "\n")
	
	km_out.close()

if "ki" in params_to_generate:	

	""" Determine the Ki values for all metabolites associated with an E.C. number """

	ki_out = open("".join((outfile, "_ki.tsv")), "w")
	ki_out.write("ecNumber\tkiValue\tkiValueMaximum\tinhibitor\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

	for l in open(ec_file):
		if use_organism == 'TRUE':
			resultString = client.getKiValue("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
		else:
			resultString = client.getKiValue("".join(("ecNumber*", re.split('\n',  l)[0]))).split('!')
		 
		for entry in resultString: #split each Brenda Km entry into a separate line
		
			all_fields = entry.split('#')
			all_fields.remove('')
		
			if len(all_fields) != 8:
				continue
			observedLigands.add(all_fields[7].split('*')[1]) #add ligandID to set of ligand IDs
		
			entry_write = []
			for field in all_fields:
				entry_write.append(field.split('*')[1])
			entry_write[4] = removeNonAscii(entry_write[4]) #some description include non ascii characters
		
		
			ki_out.write("\t".join(entry_write) + "\n")
	
	ki_out.close()

if "activators" in params_to_generate:

	""" Determine the activity of activators values for all metabolites associated with an E.C. number """

	act_out = open("".join((outfile, "_activators.tsv")), "w")
	act_out.write("ecNumber\tactivatingCompound\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

	for l in open(ec_file):
		if use_organism == 'TRUE':
			resultString = client.getActivatingCompound("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
		else:
			resultString = client.getActivatingCompound("".join(("ecNumber*", re.split('\n',  l)[0]))).split('!')
	 	
		for entry in resultString: #split each Brenda Km entry into a separate line
			print entry
			
			entry_write = []
		
			all_fields = entry.split('#')
			all_fields.remove('')
		
			if len(all_fields) != 6:
				continue
			observedLigands.add(all_fields[5].split('*')[1]) #add ligandID to set of ligand IDs
		
			for field in all_fields:
				entry_write.append(field.split('*')[1])
			entry_write[2] = removeNonAscii(entry_write[2]) #some description include non ascii characters
		
			act_out.write("\t".join(entry_write) + "\n")
	
	act_out.close()



