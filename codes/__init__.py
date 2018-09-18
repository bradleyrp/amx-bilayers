#!/usr/bin/env python

"""
AUTOMACS BILAYERS

Codes for making protein-membrane systems using automacs.
"""

from .bilayer import *
from .bilayer_bookkeeping import *
from .adhere_protein_bilayer import *

if False:
	import glob,os,importlib,re
	import ortho
	from ortho.imports import importer
	for i in glob.glob(os.path.join(os.path.dirname(__file__),'*.py')):
		print('status','importing %s'%i)
		#mod = importlib.import_module(re.match(r'^(.*?)\.py',i).group(1))
		mod = ortho.importer(re.match(r'^(.*?)\.py',i).group(1))
		globals().update(**dict([(i,j) 
			for i,j in mod.__dict__.items() if not i.startswith('_')]))
	import pdb;pdb.set_trace()
