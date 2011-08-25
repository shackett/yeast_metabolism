import sys
import itertools
from libsbml import *


### FUNCTIONS

def concat(lst):
	return itertools.chain(*lst)


def get_children(node, name=None):
	"""get all children of an XMLNode, optionally specifying the name"""
	children = [node.getChild(i) for i in range(node.getNumChildren())]
	if name != None:
		return [c for c in children if c.getName() == name]
	return children


def get_chain_down(node, *names):
	"""
	given an XMLNode and a list of names, return all children fitting
	a chain of children with those name
	"""
	if len(names) == 0:
		return [node]
	children = get_children(node, names[0])
	if len(children) == 0:
		# there were no children with that name
		return []
	return concat([get_chain_down(c, *names[1:]) for c in children])


def write_resource_file(model, outfile):
	"""
	write a tab delimited file where each line is one resource, in the form
	SpeciesInternalID\SpeciesName\tResource
	"""
	#Species types are a specific metabolite/enzyme... ignoring localization.  
	#This framework also includes the references tying the internal names to
	#external standardized metabolite and yeast protein names
	outf = open(outfile, "w")
	speciesT = model.getListOfSpeciesTypes()
	for ele in speciesT:
		#print dir(ele.getAnnotation())
		#.getAttrURI()
		#print "\n".join(dir(ele.getAnnotation().getChild()))
		annotation = ele.getAnnotation()
		if annotation == None:
			print "No annotation:", ele.getName()
			continue
		#print dir(annotation)
		children = get_chain_down(annotation, "RDF", "Description", None, "Bag", "li")
		for c in children:
			ret = "\t".join((ele.getId(), ele.getName(),
							 c.getAttributes().getValue("rdf:resource")))
			outf.write(ret + "\n")


def IDdict(model):
	""" Create a corresponence between species type designation and ID """
	name_corr = dict()
	listSpecies = model.getListOfSpecies()
	for lst in listSpecies:
		name_corr[lst.getId()] = lst.getSpeciesType()
	return name_corr

### MAIN ###

#if __name__ == "__main__":
#	xmlfile = sys.argv[1]
#	outfile = sys.argv[2]
#	reader = SBMLReader()
#	doc = reader.readSBML(xmlfile)
#	model = doc.getModel()
#	write_resource_file(model, outfile)



xmlfile = sys.argv[1]
outfile = sys.argv[2]
reader = SBMLReader()
doc = reader.readSBML(xmlfile)
model = doc.getModel()

#reaction_species_write(model, spec_outfile, rxn_outfile)

listrxn = model.getListOfReactions() 

rxn_outfile = "".join(("rxn_", outfile, ".txt"))
rxn_par_outfile = "".join(("rxn_par_", outfile, ".txt"))


rxn_outf = open(rxn_outfile, 'w')
rxn_outf.close()	
par_outf = open(rxn_par_outfile, 'w')
par_outf.close()

for rxn in listrxn:
	
	### Get the reactants and products and their stoichiometries, along with any rxn modifiers
	
	rxn_outf = open(rxn_outfile, 'a')
	for lst in rxn.getListOfReactants():
		rxn_outf.write("%s\t%s\t%s\n" % (rxn.getName(), lst.getSpecies(), 
		lst.getStoichiometry()*-1))
	for lst in rxn.getListOfProducts():
		rxn_outf.write("%s\t%s\t%s\n" % (rxn.getName(), lst.getSpecies(), 
		lst.getStoichiometry()))
	rxn_outf.close()
	
	par_outf = open(rxn_par_outfile, 'a')
	parameters = list(rxn.getKineticLaw().getListOfParameters())
	if len(parameters) > 0:
		for p in parameters:
			par_outf.write("%s\t%s\t%s\t%s\n" % (rxn.getName(), p.getId(), p.getValue(), p.getUnits()))
	par_outf.close()
	
	




#listParams = model.getListOfParameters()
#for par in listParams:
#	print "\t", [par.getId(), par.getValue(), par.getUnits()]
	
	

	
		
	
	