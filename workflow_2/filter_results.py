#!/usr/bin/env python3
DESCRIPTION = '''
Filter results from MAGI/MAGI2 using compound_score, e_score_r2g (reaction-to-gene), e_score_g2r (gene-to-reaction), and reciprocal_score.
'''
import sys
import os
import argparse
import logging
import gzip
import pandas as pd
from compound_to_reaction import prepare_smiles


## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-i', '--input', metavar='magi_results.csv', 
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Input [gzip] MAGI results file (default: stdin)'
	)
	parser.add_argument('-o', '--output', metavar='magi_results.filtered.tsv', 
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'), 
		help='Output [gzip] filtered MAGI results file (default: stdout)'
	)
	parser.add_argument('--compound_score', 
		required=False, default=1, type=float,
		help='Filter compound_score == x (default: %(default)s)'
	)
	parser.add_argument('--e_score_r2g',
		required=False, default=5, type=float,
		help='Filter e_score_r2g > x (default: %(default)s)'
	)
	parser.add_argument('--e_score_g2r',
		required=False, default=5, type=float,
		help='Filter e_score_g2r > x (default: %(default)s)'
	)
	parser.add_argument('--reciprocal_score',
		required=False, default=2, type=float,
		help='Filter reciprocal_score == x (default: %(default)s)'
	)
	parser.add_argument('--debug', 
		required=False, action='store_true', 
		help='Print DEBUG info (default: %(default)s)'
	)
	args = parser.parse_args()
	
	## Set up basic debugger
	logFormat = "[%(levelname)s]: %(message)s"
	logging.basicConfig(format=logFormat, stream=sys.stderr, level=logging.INFO)
	if args.debug:
		logging.getLogger().setLevel(logging.DEBUG)
	
	logging.debug('%s', args) ## DEBUG
	
	os.environ['NUMEXPR_MAX_THREADS'] = str(1)
	
	with args.input as input_fh, args.output as output_fh:
		filter_results(input_fh, output_fh, args.compound_score, args.e_score_r2g, args.e_score_g2r, args.reciprocal_score)



def filter_results(input_fh, output_fh, compound_score_cutoff, e_score_r2g_cutoff, e_score_g2r_cutoff, reciprocal_score_cutoff):
	df = pd.read_csv(input_fh)
	logging.info("Loaded results files with %s rows.", len(df.index)) ## INFO
	df_filtered = df[
			(df['compound_score']   == compound_score_cutoff  ) & 
			(df['e_score_r2g']       > e_score_r2g_cutoff     ) & 
			(df['e_score_g2r']       > e_score_g2r_cutoff     ) & 
			(df['reciprocal_score'] == reciprocal_score_cutoff)
		]
	logging.info(" - %s rows remain after filtering.", len(df_filtered.index)) ## INFO
	df_filtered.to_csv(output_fh, sep='\t', index=False)



def safe_InchiToSmiles(InChI):
	# Converting InChI to SMILES (None if any error occurs)
	try:
		SMILES = MolToSmiles(MolFromInchi(InChI, sanitize=True))
	except Exception as e:
		logging.info("Problem converting %s into SMILES - setting as None.", InChI)
		SMILES = None
	# Try to convert SMILES to mol
	#   - Unfortinatly the SMILES produced can be invalud (for some reason) and kills MAGI2. 
	#     Best to check now that we can actually do the conversion now, save having to handle
	#     problem later. If the conversion fails, return None.
	try:
		prepare_smiles(SMILES)
	except Exception as e:
		logging.info("SMILES: %s failed prepare_smiles conversion - setting as None.", SMILES)
		print(e)
		SMILES = None
	return(SMILES)



def load_mz_values(input_fh):
	mz_list = []
	with input_fh as infile:
		for line in infile:
			line = line.strip()
			if not line or line.startswith('#'):
				continue
			if str(line) == "original_compound":
				continue
			mz_list.append(float(line))
	return mz_list



class File(object):
	'''
	Context Manager class for opening stdin/stdout/normal/gzip files.

	 - Will check that file exists if mode='r'
	 - Will open using either normal open() or gzip.open() if *.gz extension detected.
	 - Designed to be handled by a 'with' statement (other wise __enter__() method wont 
	    be run and the file handle wont be returned)
	
	NOTE:
		- Can't use .close() directly on this class unless you uncomment the close() method
		- Can't use this class with a 'for' loop unless you uncomment the __iter__() method
			- In this case you should also uncomment the close() method as a 'for'
			   loop does not automatically cloase files, so you will have to do this 
			   manually.
		- __iter__() and close() are commented out by default as it is better to use a 'with' 
		   statement instead as it will automatically close files when finished/an exception 
		   occures. 
		- Without __iter__() and close() this object will return an error when directly closed 
		   or you attempt to use it with a 'for' loop. This is to force the use of a 'with' 
		   statement instead. 
	
	Code based off of context manager tutorial from: https://book.pythontips.com/en/latest/context_managers.html
	'''
	def __init__(self, file_name, mode):
		## Upon initializing class open file (using gzip if needed)
		self.file_name = file_name
		self.mode = mode
		
		## Check file exists if mode='r'
		if not os.path.exists(self.file_name) and mode == 'r':
			raise argparse.ArgumentTypeError("The file %s does not exist!" % self.file_name)
	
		## Open with gzip if it has the *.gz extension, else open normally (including stdin)
		try:
			if self.file_name.endswith(".gz"):
				#print "Opening gzip compressed file (mode: %s): %s" % (self.mode, self.file_name) ## DEBUG
				self.file_obj = gzip.open(self.file_name, self.mode+'b')
			else:
				#print "Opening normal file (mode: %s): %s" % (self.mode, self.file_name) ## DEBUG
				self.file_obj = open(self.file_name, self.mode, encoding='utf-8-sig')
		except IOError as e:
			raise argparse.ArgumentTypeError('%s' % e)
	def __enter__(self):
		## Run When 'with' statement uses this class.
		#print "__enter__: %s" % (self.file_name) ## DEBUG
		return self.file_obj
	def __exit__(self, type, value, traceback):
		## Run when 'with' statement is done with object. Either because file has been exhausted, we are done writing, or an error has been encountered.
		#print "__exit__: %s" % (self.file_name) ## DEBUG
		self.file_obj.close()


SCRIPT_DIR=os.path.dirname(os.path.realpath(__file__))
UNIQUE_COMPOUNDS_FILE=SCRIPT_DIR+"/unique_compounds.csv.gz"

ADDUCT_INFO = {'[2M+H]': {'charge': '1',
              'color': '#fabebe',
              'common': True,
              'comp_num': '2',
              'mass': '1.0073'},
             '[2M-H]': {'charge': '-1',
              'color': '#008080',
              'common': True,
              'comp_num': '2',
              'mass': '-1.0073'},
             '[M+2H]': {'charge': '2',
              'color': '#ffe119',
              'common': True,
              'comp_num': '1',
              'mass': '2.0146'},
             '[M+2Na]': {'charge': '2',
              'color': '#fffac8',
              'common': False,
              'comp_num': '1',
              'mass': '45.9784'},
             '[M+Cl]': {'charge': '-1',
              'color': '#d2f53c',
              'common': True,
              'comp_num': '1',
              'mass': '34.9694'},
             '[M+H-H2O]': {'charge': '1',
              'color': '#911eb4',
              'common': True,
              'comp_num': '1',
              'mass': '-17.0033'},
             '[M+H]': {'charge': '1',
              'color': '#3cb44b',
              'common': True,
              'comp_num': '1',
              'mass': '1.0073'},
             '[M+K]': {'charge': '1',
              'color': '#aa6e28',
              'common': False,
              'comp_num': '1',
              'mass': '38.963158'},
             '[M+NH4]': {'charge': '1',
              'color': '#0082c8',
              'common': True,
              'comp_num': '1',
              'mass': '18.0338'},
             '[M+Na]': {'charge': '1',
              'color': '#f58231',
              'common': True,
              'comp_num': '1',
              'mass': '22.9892'},
             '[M+acetate]': {'charge': '-1',
              'color': '#808000',
              'common': False,
              'comp_num': '1',
              'mass': '59.0139'},
             '[M-2H]': {'charge': '-2',
              'color': '#f032e6',
              'common': True,
              'comp_num': '1',
              'mass': '-2.014552904'},
             '[M-H+2Na]': {'charge': '1',
              'color': '#000080',
              'common': False,
              'comp_num': '1',
              'mass': '44.9711'},
             '[M-H+Cl]': {'charge': '-2',
              'color': '#ffd8b1',
              'common': False,
              'comp_num': '1',
              'mass': '33.9621'},
             '[M-H+Na]': {'charge': '0',
              'color': '#e6beff',
              'common': False,
              'comp_num': '1',
              'mass': '21.98194425'},
             '[M-H]': {'charge': '-1',
              'color': '#46f0f0',
              'common': True,
              'comp_num': '1',
              'mass': '-1.0073'},
             '[M-e]': {'charge': '1',
              'color': '#aaffc3',
              'common': False,
              'comp_num': '1',
              'mass': '-0.0005'},
             '[M]': {'charge': '0',
              'color': '#e6194b',
              'common': True,
              'comp_num': '1',
              'mass': '0'}}



if __name__ == '__main__':
        main()
