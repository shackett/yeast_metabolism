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
	""" Create a correspondence between species type designation and ID """
	name_corr = dict()
	listSpecies = model.getListOfSpecies()
	for lst in listSpecies:
		name_corr[lst.getId()] = lst.getSpeciesType()
	return name_corr

def write_rxn(model, outfile, rxn_par_found = "TRUE", rxn_kineticForm = "TRUE"
, use_modifiers = "TRUE"):

	"""writes a rxn-file describing the stoichiometry and (optionally) modifiers of a
	chemical reaction.  writes a spec-file describing the association between a metabolite/
	protein internal name in a rxn/compartment and a single chemical-specific identifier; 
	also notes the compartment.  optionally writes reaction parameters and reaction kinetics."""
	
	listrxn = model.getListOfReactions()
	listspecT = model.getListOfSpeciesTypes()
	listspec = model.getListOfSpecies()
	
	rxn_outfile = "".join(("rxn_", outfile, ".tsv"))
	rxn_outf = open(rxn_outfile, 'w')
	rxn_outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("Reaction", "ReactionID", "Compartment", "Metabolite", "MetName", 
	"StoiCoef"))
	rxn_outf.close()	
	
	spec_outfile = "".join(("spec_", outfile, ".tsv"))
	spec_outf = open(spec_outfile, 'w')
	spec_outf.write("%s\t%s\t%s\t%s\n" % ("SpeciesID", "SpeciesName", "SpeciesType"
	, "Compartment"))
	spec_outf.close()
	
	if rxn_par_found == "TRUE":
		rxn_par_outfile = "".join(("rxn_par_", outfile, ".tsv"))
		par_outf = open(rxn_par_outfile, 'w')
		par_outf.write("%s\t%s\t%s\t%s\n" % ("ReactionID", "ParID", "Value", "Units"))
		par_outf.close()
	
		rxn_kinetic_outfile = "".join(("rxn_kinetic_", outfile, ".tsv"))
		kin_outf = open(rxn_kinetic_outfile, 'w')
		kin_outf.write("%s\t%s\t%s\n" % ("ReactionID", "ParID", "Formula"))
		kin_outf.close()
	
	#correspondence between internal name and a common name
	met_namedict = dict()
	for spect in listspecT:
		met_namedict[spect.getId()] = spect.getName()
	
	#correspondence between internal name, unique chemical name and compartment
	comp_dict = dict()
	name_trans = dict()
	spec_outf = open(spec_outfile, 'a')	
	for spec in listspec:
		comp_dict[spec.getId()] = spec.getCompartment()
		name_trans[spec.getId()] = spec.getSpeciesType()
		spec_outf.write("%s\t%s\t%s\t%s\n" % (spec.getId(), 
		met_namedict[spec.getSpeciesType()], spec.getSpeciesType(), spec.getCompartment()))
	spec_outf.close()	
	
	for rxn in listrxn:
		
	### Get the reactants and products and their stoichiometries, along with any rxn modifiers
		rxn_location = reaction_class(rxn, comp_dict)
		
		rxn_outf = open(rxn_outfile, 'a')
		for lst in rxn.getListOfReactants():
			rxn_outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (rxn.getName(), rxn.getId(), rxn_location, lst.getSpecies(), 
			met_namedict[name_trans[lst.getSpecies()]], lst.getStoichiometry()*-1))
			
		for lst in rxn.getListOfProducts():
			rxn_outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (rxn.getName(), rxn.getId(), rxn_location, lst.getSpecies(), 
			met_namedict[name_trans[lst.getSpecies()]], lst.getStoichiometry()))
		
		if use_modifiers == "TRUE":
			for lst in rxn.getListOfModifiers():
				rxn_outf.write("%s\t%s\t%s\t%s\t%s\n" % (rxn.getName(), rxn.getId(), rxn_location, lst.getSpecies(), met_namedict[name_trans[lst.getSpecies()]]))
		rxn_outf.close()	
				
		if rxn_par_found == "TRUE":
		
			parameters = list(rxn.getKineticLaw().getListOfParameters())
			if len(parameters) > 0:
				par_outf = open(rxn_par_outfile, 'a')
				for p in parameters:
					par_outf.write("%s\t%s\t%s\t%s\n" % (rxn.getId(), p.getId(), 
					p.getValue(), p.getUnits()))
				par_outf.close()
				
			kin_outf = open(rxn_kinetic_outfile, 'a')
			kin_outf.write("%s\t%s\n" % (rxn.getId(), rxn.getKineticLaw().getFormula()))
			kin_outf.close()
			

def reaction_class(rxn, comp_dict):
	
	spec_list = [el.getSpecies() for el in rxn.getListOfReactants()]
	spec_list.extend([el.getSpecies() for el in rxn.getListOfProducts()])
	comp = [comp_dict[el] for el in spec_list]
	
	if len(unique_list(comp)) == 1:
		return unique_list(comp)[0]
	else:
		return "exchange"
	
	
	
def unique_list(seq):
    keys = {}
    for e in seq:
        keys[e] = 1
    return keys.keys()	
	
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

species_par_outfile = "".join(("species_par_", outfile, ".tsv"))
write_resource_file(model, species_par_outfile)

write_rxn(model, outfile, rxn_par_found = "TRUE", rxn_kineticForm = "TRUE", 
use_modifiers = "TRUE")


 

	
	
			

	
		
	
	