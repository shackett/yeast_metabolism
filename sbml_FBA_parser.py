
import sys
sys.path.append("/Users/Sean/Documents/MiscProgramsPackages/build/src/bindings/python")

import itertools
import re
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
	speciesT = model.getListOfSpecies()
	
	
	count = 0
	for ele in speciesT:
		
		count += 1
		
		annotation = ele.getAnnotation()
		
		if annotation == None:
			print "No annotation:", ele.getName()
			continue
		
		children = get_chain_down(annotation, "RDF", "Description", None, "Bag", "li")
		
		for c in children:			
			ret = "\t".join((ele.getId(), ele.getName(),
							 c.getAttributes().getValue("resource")))
			outf.write(ret + "\n")
		
		#if count > 10:
	#		break
	outf.close()
	
def write_compartment_file(model, outfile):
	"""
	write a tab delimited file where each line is one resource, in the form
	CompartmentID\tComparmentName
	"""

	comp_outfile = "".join(("comp_", outfile, ".tsv"))
	comp_outf = open(comp_outfile, 'w')
	comp_outf.write("%s\t%s\n" % ("compID", "compName"))
	comp_outf.close()

	compT = model.getListOfCompartments()
	comp_outf = open(comp_outfile, 'a')
	for ele in compT:
		comp_outf.write("%s\t%s\n" % (ele.getId(), ele.getName()))
	comp_outf.close()




def IDdict(model):
	""" Create a correspondence between species type designation and ID """
	name_corr = dict()
	listSpecies = model.getListOfSpecies()
	for lst in listSpecies:
		name_corr[lst.getId()] = lst.getSpeciesType()
	return name_corr



def write_rxn(model, outfile, TSdict, rIDtoEnz, use_modifiers = "TRUE"):

	"""writes a rxn-file describing the stoichiometry and (optionally) modifiers of a
	chemical reaction.  writes a spec-file describing the association between a metabolite/
	protein internal name in a rxn/compartment and a single chemical-specific identifier; 
	also notes the compartment.  optionally writes reaction parameters and reaction kinetics."""
	
	listrxn = model.getListOfReactions()
	listspec = model.getListOfSpecies()
	listcomp = model.getListOfCompartments()
	
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
	
	rxn_par_outfile = "".join(("rxn_par_", outfile, ".tsv"))
	rxn_par_outf = open(rxn_par_outfile, 'w')
	rxn_par_outf.write("%s\t%s\t%s\n" % ("ReactionID", "Enzymes", "Annotation"))
	rxn_par_outf.close()
		
	#correspondence between internal name (sID) and a common name
	met_namedict = dict()
	for spect in listspec:
		met_namedict[spect.getId()] = spect.getName()
		
	
	#correspondence between internal name, unique chemical name and compartment
	comp_dict = dict()
	name_trans = dict()
	spec_outf = open(spec_outfile, 'a')	
	for spec in listspec:
		comp_dict[spec.getId()] = spec.getCompartment()
		name_trans[spec.getId()] = TSdict[spec.getId()]
		
		spec_outf.write("%s\t%s\t%s\t%s\n" % (spec.getId(), 
		met_namedict[spec.getId()], TSdict[spec.getId()], spec.getCompartment()))
	spec_outf.close()	
	
	rxn_outf = open(rxn_outfile, 'a')
	for rxn in listrxn:
		
		### Get the reactants and products and their stoichiometries, along with any rxn modifiers
		rxn_location = reaction_class(rxn, comp_dict)
		rID = rxn.getId()
		
		for lst in rxn.getListOfReactants():
			
			rxn_outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (rxn.getName(), rID, rxn_location, lst.getSpecies(), 
			met_namedict[lst.getSpecies()], lst.getStoichiometry()*-1))
			
		for lst in rxn.getListOfProducts():
			rxn_outf.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (rxn.getName(), rID, rxn_location, lst.getSpecies(), 
			met_namedict[lst.getSpecies()], lst.getStoichiometry()))
		
		if use_modifiers == "TRUE":
			for lst in rxn.getListOfModifiers():
				rxn_outf.write("%s\t%s\t%s\t%s\t%s\n" % (rxn.getName(), rID, rxn_location, lst.getSpecies(), met_namedict[lst.getSpecies()]))
	
	rxn_outf.close()
	
	
	
	
	rxn_par_outf = open(rxn_par_outfile, 'a')
	for rxn in listrxn:
	
		### write reaction annotations - chiefly kegg rxn and enzymes
	
		rID = rxn.getId()
		
		annotation = rxn.getAnnotation()
		
		if annotation == None:
			print "No annotation:", rxn.getName()
			rxn_ref = 'NA'
		else:
			
			children = get_chain_down(annotation, "RDF", "Description", None, "Bag", "li")
			rxn_ref = []
			for c in children:
				rxn_ref.append(c.getAttributes().getValue("resource"))
			rxn_ref = ",".join(rxn_ref)
						
		if rID not in rIDtoEnz.keys():
			Enz = 'NA'
		else:
			Enz = rIDtoEnz[rID]
		
		rxn_par_outf.write("%s\t%s\t%s\n" % (rID, Enz, rxn_ref))
		
		
	rxn_par_outf.close()	
		
		
				
		
			

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


def define_groups(textModel):

	""" Read through and define groups (since they are not defined in libsbml and
	create a dictionary associating sID with tID """

	tID = re.compile('groups:group groups')
	sID = re.compile('groups:member')
	TSdict = dict()

	for line in textModel:
	
		if tID.search(line):
			tIDstore = re.search('groups:id="(.*?)" groups:name=', line).group(1)
		
		elif sID.search(line):
			TSdict[re.search('groups:symbol="(.*?)"/>', line).group(1)] = tIDstore	
		
	return TSdict




def write_directionality(textModel, outfile):

	""" Read through and scrape FluxBounds recording rID, bound and constraint direction """
	
	flux_bound_outfile = "".join(("flux_dir_", outfile, ".tsv"))
	flux_outf = open(flux_bound_outfile, 'w')
	flux_outf.write("%s\t%s\t%s\n" % ("ReactionID", "Reversible", "FluxBound"))
	
	rIDrev = dict()
	rxLine = re.compile('<reaction')
	
	rIDbound = dict()
	fluxBoundLine = re.compile('fbc:fluxBound fbc:reaction')
	
	for line in textModel:
		
		if rxLine.search(line):
			rIDrev[re.search('" id="(.*?)" name="', line).group(1)] = re.search('reversible="(.*?)" fast=', line).group(1)
		if fluxBoundLine.search(line):
			rIDbound[re.search('fbc:reaction="(.*?)" fbc:', line).group(1)] = re.search('fbc:operation="(.*?)" fbc:value', line).group(1)
	
	for ID in rIDrev.keys():
		if ID not in rIDbound.keys():
			rIDbound[ID] = "reversible"
			
	for ID in rIDrev.keys():
		flux_outf.write("%s\t%s\t%s\n" % (ID, rIDrev[ID], rIDbound[ID]))




def reaction_enzymes(textModel):
	
	rIDfind = re.compile('reaction metaid')
	Enzfind = re.compile('<p>GENE_ASSOCIATION:')
	
	rIDtoEnz = dict()
	
	for line in textModel:
		
		if rIDfind.search(line):
			rIDstore = re.search('" id="(.*?)" name="', line).group(1)
		
		if Enzfind.search(line):
			rIDtoEnz[rIDstore] = re.search('<p>GENE_ASSOCIATION:(.*?)</p>', line).group(1)
	
	return rIDtoEnz
	
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
textModel = list(open(xmlfile,'r'))






species_par_outfile = "".join(("species_par_", outfile, ".tsv"))
write_resource_file(model, species_par_outfile)

write_compartment_file(model, outfile)

write_directionality(textModel, outfile)


TSdict = define_groups(textModel)
rIDtoEnz = reaction_enzymes(textModel)
write_rxn(model, outfile, TSdict, rIDtoEnz, use_modifiers = "TRUE")


 

	
	
			

	
		
	
	