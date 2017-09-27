import os
import subprocess

repo_path = os.getcwd()

fname = raw_input('USER INPUT: Settings Name (leave blank for default): ')

# step one: extract the db tarball
print 'Extracting database files...'
db_dir = os.path.join(repo_path, 'workflow', 'database')
# wipe out existing db directory
if os.path.isdir(db_dir):
    subprocess.call(['rm', '-r', db_dir])
os.makedirs(db_dir)
for zipfile in ['blastdb.zip', 'pickle_files.zip']:
    cmd = ['unzip', '../%s' %(zipfile)]
    subprocess.call(
        cmd,
        cwd = os.path.join(repo_path, 'workflow/database')
        )
print 'Done'

# step two: make local_settings.py file

if fname == '':
    fname = 'user_settings'

lines = ["import os\n\n"]
lines.append("repo_location = '%s'\n\n" % (repo_path))
lines.append("blastbin =          os.path.join(repo_location, 'workflow/blastbin')\n")
lines.append("\n")
lines.append("refseq_path =       os.path.join(repo_location, 'workflow/database/mrs_reaction_filtered_refseq_db_newrxns_actinofixed.pkl')\n")
lines.append("refseq_db =         os.path.join(repo_location, 'workflow/database/mrs_reaction_filtered_refseq_db_fasta_blastdb_actinofix')\n")
lines.append("mrs_reaction_path = os.path.join(repo_location, 'workflow/database/mrs_reaction_newrxns_added_actinofix.pkl')\n")
lines.append("compounds_df =      os.path.join(repo_location, 'workflow/database/unique_compounds_groups_magi.pkl')\n")
lines.append("mst_path =          os.path.join(repo_location, 'workflow/database/graph.pkl')\n")
lines.append("chemnet_pickle =    os.path.join(repo_location, 'workflow/database/compound_groups.pkl')\n")
lines.append("c2r =               os.path.join(repo_location, 'workflow/database/c2r.pkl')\n")
lines.append("\n")
lines.append("magi_results_storage = os.path.join(repo_location, 'outputs')\n")

with open('local_settings/%s.py' %(fname), 'w') as f:
    for line in lines:
        f.write(line)
with open('local_settings/local_settings.py', 'a') as f:
    f.write("SETTINGS_FILE= '%s'\n" % (fname))
print 'Successfully wrote local settings files'

# steps 3 and 4: change localsettings paths
def change_localsettings_path(fname):
    """
    Requires the sys.path.insert to be all on one line
    Requires "from local_settings import local_settings as settings_loc"
        to immediately follow sys.path.insert
    overwrites existing file!
    """
    
    with open(fname, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('sys.path.insert'):
            if 'from local_settings import local_settings as settings_loc' in lines[i+1]:
                break
    if i+1 == len(lines):
        raise RuntimeError('Could not find the line to replace in %s!' % (fname))
    lines[i] = "sys.path.insert(0, '%s')\n" % (repo_path)
    with open(fname, 'w') as f:
        f.write(''.join(lines))
    print 'Successfully adjusted %s' % (fname)
files = [
    'workflow/workflow_helpers.py',
    'workflow/magi_workflow_20170519.py'
]
for fname in files:
    change_localsettings_path(fname)

print 'Setup Successful!'
