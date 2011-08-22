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
reader = SBMLReader()
doc = reader.readSBML(xmlfile)
model = doc.getModel()

listrxn = model.getListOfReactions() 
for rxn in listrxn:
	print rxn.getName()

	### Get the reactants and products and their stoichiometries, along with any rxn modifiers
	
	for lst in rxn.getListOfReactants():
		print [lst.getSpecies(), lst.getStoichiometry()*-1]
	for lst in rxn.getListOfProducts():
		print [lst.getSpecies(), lst.getStoichiometry()]
	print [m.getSpecies() for m in rxn.getListOfModifiers()]
	
	### Get The parameters of the rxn ###
	
	parameters = list(rxn.getKineticLaw().getListOfParameters())
	if len(parameters) > 0:
		print [(p.getId(), p.getValue(), p.getUnits()) for p in parameters]
	print rxn.getKineticLaw().getFormula()
	break




listParams = model.getListOfParameters()
for par in listParams:
	print "\t", [par.getId(), par.getValue(), par.getUnits()]
	
	

	
		
	
	