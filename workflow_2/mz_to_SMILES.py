#!/usr/bin/env python3
DESCRIPTION = '''
Convert m/z values to all possible structures given a set of standard adducts.
'''
import sys
import os
import argparse
import logging
import gzip
import pandas as pd
from rdkit.Chem import MolFromInchi, MolToSmiles


## Pass arguments.
def main():
	## Pass command line arguments. 
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
	parser.add_argument('-i', '--input', metavar='mz_values.txt', 
		required=False, default=sys.stdin, type=lambda x: File(x, 'r'), 
		help='Input [gzip] file listing m/z values to convert to SMILES (default: stdin)'
	)
	parser.add_argument('-o', '--out', metavar='smiles.txt', 
		required=False, default=sys.stdout, type=lambda x: File(x, 'w'), 
		help='Output [gzip] file with SMILES (default: stdout)'
	)
	parser.add_argument('-n', '--ncpus', 
		required=False, default=4, type=int, 
		help='Number of threads to use (default: %(default)s)'
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
	
	os.environ['NUMEXPR_MAX_THREADS'] = str(args.ncpus)
	
	mz_list = load_mz_values(args.input)
	mz_to_SMILES(mz_list, args.out)



def mz_to_SMILES(mz_list, output_fh):
	cols = ['formula','inchi', 'inchi_key','mono_isotopic_molecular_weight', 'name',]
	cpds = pd.read_csv(UNIQUE_COMPOUNDS_FILE, usecols=cols)
	cpds['merge_key'] = 1
	MY_CHARGE = '1'# or '-1' for NEGATIVE MODE
	PPM_CUTOFF = 5 # five parts per million mass accuracy is fairly standard for Orbitrap mass spectrometers
	
	## Process adduct info
	adducts = {k:v for k,v in ADDUCT_INFO.items() if (v['charge']=='1') & (v['comp_num']=='1') & (v['common']==True)}
	adducts = pd.DataFrame(adducts).T
	adducts.index.name = 'adduct'
	adducts.reset_index(inplace=True,drop=False)
	adducts['mass'] = adducts['mass'].astype(float)
	adducts['merge_key'] = 1
	if MY_CHARGE == '1':
		adducts.loc[-1] = ['[M+]',1,'',True,1,0,1]
	
	## Convert target m/z values
	targets = pd.DataFrame(data=mz_list,columns=['mz'])
	targets['merge_key'] = 1
	targets = pd.merge(adducts,targets,on='merge_key',how='outer')
	targets['theoretical_mass'] = targets['mz'] + targets['mass']
	
	## Convert
	n = 1000  #chunk row size
	out = []
	for i in range(0,cpds.shape[0],n):
		temp = pd.merge(targets,cpds[i:i+n],how='outer',on='merge_key')
		temp['error'] = abs(temp['theoretical_mass'] - temp['mono_isotopic_molecular_weight']) / temp['mono_isotopic_molecular_weight'] * 1e6
		temp = temp[temp['error']<PPM_CUTOFF]
		if temp.shape[0]>0:
			out.append(temp)
	
	## Convert output table ready for MAGI2
	out = pd.concat(out)
	out['original_compound'] = out['inchi'].apply(lambda x: safe_InchiToSmiles(x))
	with output_fh as outfile:
		out.to_csv(outfile, sep='\t', index=False)



def safe_InchiToSmiles(InChI):
	try:
		SMILES = MolToSmiles(MolFromInchi(InChI))
	except Exception as e:
		logging.info("Problem converting " + InChI + " into SMILES - setting as None.")
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


UNIQUE_COMPOUNDS_FILE="unique_compounds.csv.gz"
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
