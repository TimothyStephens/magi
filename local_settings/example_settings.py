"""
Two files need to be created.  Neither of these will be in the repo.

You must create a file called "local_settings.py" should be put in this directory with a single line in it.
This line will define a settings.py file that has all the variables declared to interact with the database

For example:
	SETTINGS_FILE = 'user_settings_bpb_laptop'

In this directory the 
For local setup, in the SETTINGS_FILE 'user_settings_bpb_laptop.py' you might find content like this:
	db_password_file = None
	if db_password_file:
	     with open(db_passwd_file) as fid:
	         pw = fid.read().strip()

	db_username = None

	db_name = 'bpb_workspace.db'

	blastbin = '/Users/bpb/repos/metatlas_reactions/bin/mac_ncbi_blast_2.6.0/bin' # path to blast binary
	fasta_path = '/Users/bpb/repos/metatlas_reactions/docs/example_notebooks/magi_refseq.fa' #path to fasta file
	magi_blast_path = '/Users/bpb/repos/metatlas_reactions/data' #where you want the db stored
	refseq_db_path = '/Users/bpb/repos/metatlas_reactions/data/BLAST_dbs/magi_refseq.db'

	sql_program = 'sqlite'

For production setup as database administrator, in the SETTINGS_FILE 'user_settings_admin_NERSC.py' you might find content like this:
	import os
	metatlas_rxn_path = '/Users/bpb/repos/metatlas_reactions/'
	db_password_file = os.path.join(metatlas_rxn_path,'db_admin_auth.txt')
	if db_password_file:
	     with open(db_passwd_file) as fid:
	         pw = fid.read().strip()
	db_username = 'admin'
	db_name = 'metatlas_rxn.db'
	sql_program = 'mysql+pymysql'
	#db_path = '%s://%s:%s@scidb1.nersc.gov/%s' % (sql_program,db_username,pw, db_name)

Finally, for production setup as database user, in the SETTINGS_FILE 'user_settings_user_NERSC.py' you might find content like this:
	import os
	metatlas_rxn_path = '/Users/bpb/repos/metatlas_reactions/'
	db_password_file = os.path.join(metatlas_rxn_path,'db_user_auth.txt')
	if db_password_file:
	     with open(db_passwd_file) as fid:
	         pw = fid.read().strip()
	db_username = 'metatlas_rxn_user'
	db_name = 'metatlas_rxn.db'
	sql_program = 'mysql+pymysql'
	#db_path = '%s://%s:%s@scidb1.nersc.gov/%s' % (sql_program,db_username,pw, db_name)

To use these settings in code in this repo, simply do this:
	from local_settings import local_settings as settings_loc
	my_settings = getattr(__import__('local_settings', fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)

	blastbin = my_settings['blastbin'] # path to blast binary
	fasta_path = my_settings['fasta_path'] #path to fasta file
	magi_blast_path = my_settings['magi_blast_path'] #where you want the db stored
	refseq_db_path = my_settings['refseq_db_path']

""" 