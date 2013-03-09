""" Using a newline delimited list of E.C. numbers and a specified organism: pull down kinetic information  """

""" pip install -e "git+http://github.com/pelletier/SOAPpy.git@develop#egg=SOAPpy """

import string
import re
from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file


ec_file = 'ECnr.txt'
organism = 'Saccharomyces cerevisiae'

endpointURL = "http://www.brenda-enzymes.info/soap2/brenda_server.php"
client = SOAPProxy(endpointURL)


observedLigands = set()

""" Determine the Km values for all metabolites associated with an E.C. number """

km_out = open("km_file.tsv", "w")
km_out.write("ecNumber\tkmValue\tkmValueMaximum\tsubstrate\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

for l in open(ec_file):
	resultString = client.getKmValue("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
	 
	for entry in resultString: #split each Brenda Km entry into a separate line
		
		
		all_fields = entry.split('#')
		all_fields.remove('')
		
		if len(all_fields) == 8:
			observedLigands.add(all_fields[7].split('*')[1]) #add ligandID to set of ligand IDs
		
		entry_write = []
		for field in all_fields:
			entry_write.append(field.split('*')[1])
		
		km_out.write("\t".join(entry_write) + "\n")
	
km_out.close()
	

""" Determine the Ki values for all metabolites associated with an E.C. number """

ki_out = open("ki_file.tsv", "w")
ki_out.write("ecNumber\tkiValue\tkiValueMaximum\tinhibitor\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

for l in open(ec_file):
	resultString = client.getKiValue("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
	 
	for entry in resultString: #split each Brenda Km entry into a separate line
		
		all_fields = entry.split('#')
		all_fields.remove('')
		
		if len(all_fields) == 8:
			observedLigands.add(all_fields[7].split('*')[1]) #add ligandID to set of ligand IDs
		
		entry_write = []
		for field in all_fields:
			entry_write.append(field.split('*')[1])
		
		ki_out.write("\t".join(entry_write) + "\n")
	
ki_out.close()


""" Determine the activity of activators values for all metabolites associated with an E.C. number """

act_out = open("activator_file.tsv", "w")
act_out.write("ecNumber\tactivatingCompound\torganism\tligandStructureId\tliterature" + "\n")

for l in open(ec_file):
	resultString = client.getActivatingCompound("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
	 
	for entry in resultString: #split each Brenda Km entry into a separate line
		
		entry_write = []
		
		all_fields = entry.split('#')
		all_fields.remove('')
		
		if len(all_fields) == 6:
			observedLigands.add(all_fields[5].split('*')[1]) #add ligandID to set of ligand IDs
		
		for field in all_fields:
			entry_write.append(field.split('*')[1])
		
		act_out.write("\t".join(entry_write) + "\n")
	
act_out.close()



