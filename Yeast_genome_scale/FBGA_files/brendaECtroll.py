""" Using a newline delimited list of E.C. numbers and a specified organism: pull down kinetic information  """

""" pip install -e "git+http://github.com/pelletier/SOAPpy.git@develop#egg=SOAPpy """

import string
import re
from SOAPpy import WSDL ## for extracting the URL of the endpoint (server script) from the WSDL file
from SOAPpy import SOAPProxy ## for usage without WSDL file


ec_file = 'ecnums.txt'
organism = 'Saccharomyces cerevisiae'

endpointURL = "http://www.brenda-enzymes.info/soap2/brenda_server.php"
client = SOAPProxy(endpointURL)


ki_values = dict()

km_out = open("km_file.tsv", "w")
km_out.write("ecNumber\tkmValue\tkmValueMaximum\tsubstrate\tcommentary\torganism\tligandStructureId\tliterature" + "\n")

for l in open(ec_file):
	resultString = client.getKmValue("".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))).split('!')
	 
	#print "".join(("ecNumber*", re.split('\n',  l)[0], "#organism*", organism))
	#print resultString[0]
	
	for entry in resultString: #split each Brenda Km entry into a separate line
		entry_write = []
		
		all_fields = entry.split('#')
		all_fields.remove('')
		
		for field in all_fields:
			entry_write.append(field.split('*')[1])
		
		km_out.write("\t".join(entry_write) + "\n")
	
km_out.close()
	





#print resultString.split('!')



resultString = client.getKmValue("ecNumber*1.1.1.1#organism*Saccharomyces cerevisiae")
#print resultString.split('!')
#print woot[1]

getActivatingCompound(String)


resultString = client.getKiValue("ecNumber*1.1.1.1#organism*Saccharomyces cerevisiae")
#print resultString.split('!')[0]


