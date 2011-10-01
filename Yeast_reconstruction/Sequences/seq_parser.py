import re
import sys

def RNAfile_Parse(mRNA, RNAoutf):

	""" Loop Through a sequence of fasta expression files and write the base
	composition, gene name and gene IDs """
	
	riboNucs = ["A", "T", "G", "C"]
	RNAcomp = dict()

	for line in mRNA:
		if line[0] == ">": 
			l_split = line.split(" ")
			gene_Name = l_split[1]
			gene_ID = re.split(':|,',  l_split[2])[1]
			l_array = [0,0,0,0]
		else:
			for i in range(len(riboNucs)):
				l_array[i] = l_array[i] + line.count(riboNucs[i])
				RNAcomp[gene_Name] = (gene_ID, l_array)
	
	outf = open(RNAoutf, "w")
	outf.write("Name\tID\tA\tU\tG\tC\n")
	
	for gene in RNAcomp.keys():
		outf.write("%s\t%s\t%d\t%d\t%d\t%d\n" % (
		gene, RNAcomp[gene][0], RNAcomp[gene][1][0], RNAcomp[gene][1][1], 
		RNAcomp[gene][1][2], RNAcomp[gene][1][3])
		)
	outf.close()
	
###############

def ORFfile_Parse(protein, PROToutf):

	""" Loop Through a sequence of sgd orf files and write the amino acid
	composition, gene name and gene IDs """


	aminoAcids = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", 
	"A", "V", "I", "L", "M", "F", "Y", "W"]
	
	PROTcomp = dict()
	for line in protein:
		if line[0] == ">": 
			l_split = line.split(" ")
			gene_Name = l_split[1]
			gene_ID = re.split(':|,',  l_split[2])[1]
			l_array = [0 for number in xrange(len(aminoAcids))]
		else:
			for i in range(len(aminoAcids)):
				l_array[i] = l_array[i] + line.count(aminoAcids[i])
				PROTcomp[gene_Name] = (gene_ID, l_array)
	
	outf = open(PROToutf, "w")
	outf.write("Name\tID\t%s\n" % ("\t".join(aminoAcids)))
	
	for gene in PROTcomp.keys():
		outf.write("%s\t%s" % (gene, PROTcomp[gene][0]))
		outf.write("\t%d\t%d\t%d\t%d\t%d" % (PROTcomp[gene][1][0], PROTcomp[gene][1][1], PROTcomp[gene][1][2], PROTcomp[gene][1][3], PROTcomp[gene][1][4]))
		outf.write("\t%d\t%d\t%d\t%d\t%d" % (PROTcomp[gene][1][5], PROTcomp[gene][1][6], PROTcomp[gene][1][7], PROTcomp[gene][1][8], PROTcomp[gene][1][9]))
		outf.write("\t%d\t%d\t%d\t%d\t%d" % (PROTcomp[gene][1][10], PROTcomp[gene][1][11], PROTcomp[gene][1][12], PROTcomp[gene][1][13], PROTcomp[gene][1][14]))
		outf.write("\t%d\t%d\t%d\t%d\t%d" % (PROTcomp[gene][1][15], PROTcomp[gene][1][16], PROTcomp[gene][1][17], PROTcomp[gene][1][18], PROTcomp[gene][1][19]))
		outf.write("\t%d\n" % (PROTcomp[gene][1][20]))
	
	outf.close()

###############
	
def CHEBIparser(metabolites, METAoutf, write_ele_comp):

	""" Loop through a CHEBI file and write a file indicating the
	correspondence between a CHEBI ID the metabolite ID"""


	ID_match = re.compile("<ChEBI ID>")
	Name_match = re.compile("<ChEBI Name>")
	Formula_match = re.compile("<Formulae>")

	last_ID_match = "FALSE"
	last_Name_match = "FALSE"
	last_Formula_match = "FALSE"
	
	meta_dict = dict()
	formula_dict = dict()
	for line in metabolites:
		if last_ID_match == "TRUE":
			cur_ID = line.strip().split(":")[1]
			last_ID_match = "FALSE"
		elif last_Name_match == "TRUE":
			meta_dict[cur_ID] = line.strip() 
			last_Name_match = "FALSE"
		elif last_Formula_match == "TRUE":
			formula_dict[cur_ID] = line.strip()
			last_Formula_match = "FALSE"
		elif ID_match.search(line):
			last_ID_match = "TRUE"
		elif Name_match.search(line):
			last_Name_match = "TRUE"
		elif Formula_match.search(line):
			last_Formula_match = "TRUE"

	for line in meta_dict.keys():
		if not (formula_dict.has_key(line)):
			formula_dict[line] = "$$NA$$"
	
	all_elements = []
	for line in formula_dict.keys():
		all_elements.extend(elements_formula(formula_dict[line], determine_comp = "TRUE"))
		print set(all_elements)


	outf = open(METAoutf, "w")
	outf.write("ID\tName\n")	
			
	for meta in meta_dict.keys():
		outf.write("%s\t%s\t%s\n" % (meta, meta_dict[meta], formula_dict[meta]))

	outf.close()

def elements_formula(formula, determine_comp = "FALSE", *elements):

	#remove formulas that are either missing ("$"), have adducts (".") or are polymers (")n")	
	if formula[0] == '$':
		sys.exit(0)
	add = re.compile('[.]')
	poly = re.compile('[)n]')
	print add
	print poly
	print add == "none"
	print poly == "none"
	
	if (add.search(formula)) | (poly.search(formula)):
		sys.exit(0)


	# seperate elements and counts
	
	lc = re.compile('[a-z]')
	uc = re.compile('[A-Z]')
	alpha = re.compile('[A-Za-z]')
	num = re.compile('[0-9]+')
	char_type = []
	character = []
	for char in formula:
		character.extend(char)
		if uc.search(char):
			char_type.extend("u")			
		elif num.search(char):
			char_type.extend("n")
		elif lc.search(char):
			char_type.extend("l")
		
	formula = []
	last_char = "NA"
	for i in range(len(char_type)):
		if (char_type[i] == "u") & (last_char != "u"):
			formula.append(character[i])
			last_char = "u"
		elif (char_type[i] == "u") & (last_char == "u"):
			formula.append("1")
			formula.append(character[i])
			last_char = "u"
		elif char_type[i] == "l":
			formula.append("".join((formula.pop(), character[i])))
			last_char = "l"
		elif (char_type[i] == "n") & (last_char != "n"):
			formula.append(character[i])
			last_char = "n"
		elif (char_type[i] == "n") & (last_char == "n"):
			formula.append(str(int(formula.pop())*10 + int(character[i])))			
			last_char = "n"
	
	# create a correspondence dictionary between elements and counts
	
	el_comp = dict()
	for i in formula:
		if alpha.search(i):
			element = i
		elif num.search(i):
			el_comp[element] = i	
		else:
			print i
	if not el_comp.has_key(element):
		el_comp[element] = 1
	
	# conditionally write out the elements found in the compound. Otherwise use a list of all
	# elements found in the dataset and write out counts of each as tsv.
	
	if determine_comp == "TRUE":
		return el_comp.keys()
		
		
		


#form = "C10O12ZHe3P"
#ele = elements_formula(form, determine_comp = "TRUE")




####### Main ##########

#mRNA = open('orf_genomic.fasta','r')
#RNAoutf = 'RNAcomposition.tsv'
#RNAfile_Parse(mRNA, RNAoutf)

#protein = open('orf_trans.fasta','r')
#PROToutf = 'PROTcomposition.tsv'
#ORFfile_Parse(protein, PROToutf)

metabolites = open('ChEBI_complete.sdf','r')
METAoutf = 'METAdict.tsv'
CHEBIparser(metabolites, METAoutf, write_ele_comp = "TRUE")


