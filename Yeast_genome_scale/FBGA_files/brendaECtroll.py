""" Using a newline delimited list of E.C. numbers and a specified organism: pull down kinetic information  """

""" pip install -e "git+http://github.com/pelletier/SOAPpy.git@develop#egg=SOAPpy """

import string
import re
from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file

def removeNonAscii(s): return "".join(filter(lambda x: ord(x)<128, s)) #remove non Ascii charcters - the angstrom symbol

ec_file = 'ECnr.txt' #a text file with E.C. numbers to be queried against BRENDA
organism = 'Saccharomyces cerevisiae'
use_organism = 'FALSE' # either look for a specific organism (listed in the organism field) or find constants for all organims
#params_to_generate = set(["km", "ki", "activators"]) #which kinetic constants should be searched for
params_to_generate = set(["ki", "activators"]) #which kinetic constants should be searched for
catchup = 'FALSE'

if use_organism == 'TRUE':
	outfile = "_".join((organism.split(' ')[0].strip()[0], organism.split(' ')[1]))
else:
	outfile = "all"
 
endpointURL = "http://www.brenda-enzymes.org/soap/brenda_server.php"
client = SOAPProxy(endpointURL)

observedLigands = set()

if "km" in params_to_generate:
	
	""" Determine the Km values for all metabolites associated with an E.C. number """

	km_out = open("".join((outfile, "_km.tsv")), "w")
	km_out.write("ecNumber\tkmValue\tkmValueMaximum\tsubstrate\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

	for l in open(ec_file):
		
		print l
		
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
			entry_write[3] = removeNonAscii(entry_write[3]) #some description include non ascii characters
			entry_write[4] = removeNonAscii(entry_write[4]) #some description include non ascii characters
		
			km_out.write("\t".join(entry_write) + "\n")
	
	km_out.close()

if "ki" in params_to_generate:	

	""" Determine the Ki values for all metabolites associated with an E.C. number """
	
	if catchup == 'TRUE': # Because the API may putz out part way through run, readin the old output and start from final EC number
		past_run = open("".join((outfile, "_ki.tsv")), "r")
		ec_finished = set()
		for l in past_run:
			ec_finished.add(l.split("\t")[0])
			last = l.split("\t")[0]
		ec_finished.remove(last)
		print ec_finished
	
		ki_out = open("".join((outfile, "_ki.tsv")), "a")
	
	else:
		ki_out = open("".join((outfile, "_ki.tsv")), "w")
		ki_out.write("ecNumber\tkiValue\tkiValueMaximum\tinhibitor\tcommentary\torganism\tligandStructureId\tliterature" + "\n")
	
	
	for l in open(ec_file):
		
		ec_current = l.split("\t")[0].strip()
		print ec_current
		
		if catchup == 'TRUE':
			if ec_current in ec_finished:
				print 'previously completed'
				continue
		
		if use_organism == 'TRUE':
			resultString = client.getKiValue("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
			resultString2 = client.getInhibitors("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
		
		else:
			
			resultString = client.getKiValue("".join(("ecNumber*", re.split('\n',  l)[0]))).split('!')
		 	
		 	if ec_current in ['3.5.3.1', '3.6.3.44', '2.7.1.32', '6.3.2.2', '2.7.1.91', '3.1.4.4','2.5.1.47','3.1.3.25']: # some EC numbers have a faulty inhibitor database which cannot be queried
		 		print "skip"
		 		resultString2 = ''
		 	else:	
		 		resultString2 = client.getInhibitors("".join(("ecNumber*", re.split('\n',  l)[0]))).split('!')
		
		# for entries with a km value
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
			#print entry_write
			ki_out.write("\t".join(entry_write) + "\n")
			
		# for entries with only a species annotation
		for entry in resultString2:
			
			all_fields = entry.split('#')
			all_fields.remove('')
			
			if len(all_fields) != 6:
				continue
			
			entry_write = []
			for field in all_fields:
				entry_write.append(field.split('*')[1])
			entry_write[1] = removeNonAscii(entry_write[1]) #some names include non ascii characters
			entry_write[2] = removeNonAscii(entry_write[2]) #some description include non ascii characters
			
			ki_out.write("%s\t%s\t\t%s\t%s\t%s\t%s\t%s\n" % (
				entry_write[0], 'NA', entry_write[1], entry_write[2], entry_write[3], entry_write[4], entry_write[5]
				))

	ki_out.close()

if "activators" in params_to_generate:

	""" Determine the activity of activators values for all metabolites associated with an E.C. number """

	act_out = open("".join((outfile, "_activators.tsv")), "w")
	act_out.write("ecNumber\tactivatingCompound\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

	for l in open(ec_file):
		
		ec_current = l.split("\t")[0].strip()
		print ec_current
		
		if ec_current in ['3.6.3.14']:
			print "skip"
			continue
		
		if use_organism == 'TRUE':
			resultString = client.getActivatingCompound("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
		else:
			resultString = client.getActivatingCompound("".join(("ecNumber*", re.split('\n',  l)[0]))).split('!')
	 	
		for entry in resultString: #split each Brenda Km entry into a separate line
			
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



