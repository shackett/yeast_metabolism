"""
CTS.py

A thin interface for the Chemical Translation Service. Requires the suds
library.

NOTE:
from/to options are
cas, chebi, cid, formula, hmdb, inchi, inchikey, kegg, lipidmap, mass,
name, sid, skeleton, smiles
"""

import sys
import re
import itertools
import unittest

from suds.client import Client


def concat(lst):
	"""concatenate a list of lists"""
	return list(itertools.chain(*lst))


class CTS(object):
	def __init__(self):
		self.client = Client("http://uranus.fiehnlab.ucdavis.edu:8080/" +
							 "cts/services/transformWeb?wsdl")

	def convert(self, from_type, to_type, term):
		"""
		given a type to convert from and to, and a search term, return a
		list of all synonyms in that type
		"""
		# this step does the logic:
		ret = self.client.service.convertFromTo(from_type, to_type, term)
		print ret
		# at this point, will be in a form that looks like
		# [[KEGG: C00186], [KEGG: C01432, KEGG: D00111]]
		# extract only the IDs (the parts after the :s)
		return concat([[m.groups()[0]
							for m in re.finditer(": ([^\,\]]+)", str(e))]
								for e in ret])

	def junk(self):
		"""random stuff we're trying"""
		#print self.client.service.transformToDeepJSONString("name", "mass", "lactic acid")


### TESTS

class CTSTest(unittest.TestCase):
	def test_convert(self):
		"""some manual checks of conversions"""
		cts = CTS()
		
		# test as a set rather than as a list, in case the order changes
		self.assertEqual(set(cts.convert("name", "cas", "lactic acid")),
						 set(['50-21-5', '26811-96-1', '79-33-4',
						 '79-33-4ChEBI', '79-33-4', '1715-99-7', '598-82-3',
						 '50-21-5', '50-21-5DrugBank', '50-21-5ChEBI']))

		self.assertEqual(set(cts.convert("name", "kegg", "lactic acid")),
						 set(["C00186", "C01432", "D00111"]))


if __name__ == "__main__":
	unittest.main()

