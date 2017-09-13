"""
This is a set of functions used for the main MAGI workflow.
Written by Onur Erbilgin and Ben Bowen at Lawrence Berkeley National Lab

Datafiles used by these functions are loaded into memory, and use the
paths provded in the local_settings file (see documentation). Specific
variables in the local_settings files used are:
- blastbin: path to directory houseing the binary of blastp and
            makeblastdb programs
- refseq_path: path to reaction reference sequence dataframe
- refseq_db: path to reaction reference sequence BLAST database
- mrs_reaction_path: path to reaction dataframe
- compounds_df: path to compounds dataframe
- chemnet_pickle: path to pickle file describing compound groups for
                  chemical network
- mst_path: path to minimum spanning tree describing the chemical
            network
"""

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All the setup junk
import sys
# local settings path
sys.path.insert(
    0,
    '/global/homes/e/erbilgin/repos/magi/')
from local_settings import local_settings as settings_loc

# needed for rdkit
sys.path.insert(
    0,
    '/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages'
    )
from rdkit import Chem
# turn off rdkit warnings
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

from molvs.standardize import enumerate_tautomers_smiles
import pandas as pd
import numpy as np
import subprocess
import time
import multiprocessing as mp
import os
import pickle
import re
import networkx as nx
from shutil import rmtree

my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)

def load_dataframe(fname, filetype=None, key=None):
    """
    Uses the appropriate pandas function to load a file based on the
    file extension:
    .pkl or .pickle: pickle file
    .xls or .xlsx: excel file
    .csv: comma separated file
    .txt, .tab, or .tsv: tab separated file
    .h5 or .hdf5: HDF5 formatted files. To load these, you must also
                  pass a key argument!

    Inputs
    ------
    fname: path to the file to be loaded as a pandas dataframe
    filetype: Used to override the autodetection of a filetype based on
              its file extension. Accepts file extensions listed above,
              but without the preceding "." (e.g. pass it "pkl")
    key: They key for the table to be loaded in the HDF5 file. Only
         required/used when loading HDF5 tables. 

    Outputs
    -------
    df: the loaded pandas dataframe

    WARNING: for binary file types (pickle and hdf5), be wary of
    saving and loading from different versions of pandas, as this will
    very likely break the loader.
    """

    if filetype is None:
        file_ext = fname.split('.')[-1]
    else:
        file_ext = filetype
    if file_ext in ['pkl', 'pickle']:
        df = pd.read_pickle(fname)
    elif file_ext in ['xls', 'xlsx']:
        df = pd.read_excel(fname)
    elif file_ext in ['csv']:
        df = pd.read_csv(fname)
    elif file_ext in ['txt', 'tab', 'tsv']:
        df = pd.read_csv(fname, sep='\t')
    elif file_ext in ['h5', 'hdf5']:
        if key is None:
            raise IOError('"key" argument must be used when loading\
                an HDF5 file')
        else:
            df = pd.read_hdf(fname, key)
    else:
        raise IOError('could not infer what type of file %s is... please \
            make it a csv, tab, or pickle file')

    # remove rows that are empty
    df = df[~pd.isnull(df).all(axis=1)]
    return df

# JGI custom compiled blast binaries
blastbin = my_settings.blastbin

print '!!! loading refseq and reaction tables...'
# this table is only refseqs that are found in mrs-reaction
refseq_path = my_settings.refseq_path
print '!!! Reference sequences in this file:', refseq_path
refseq = load_dataframe(refseq_path)
refseq.dropna(inplace=True)
print '!!!', len(refseq), 'reference sequences'

# path to the refseq database (for reverse blasting)
refseq_dbpath = my_settings.refseq_db

# third version of mrs-reaction
# with additional metacyc reactions manually added
# (those that had compounds with R-groups)
mrs_reaction_path = my_settings.mrs_reaction_path

print '!!! MRS-Reaction:', mrs_reaction_path
mrs_reaction = load_dataframe(mrs_reaction_path)
print '!!!', len(mrs_reaction), 'reactions'
print '!!!', len(mrs_reaction[mrs_reaction['refseq_id'] != '']), \
        'reactions with a refseq'
print '!!!', len(mrs_reaction[mrs_reaction['ECs'] != '']), 'reactions with an EC'

print '!!! loading compound table'
compounds = load_dataframe(my_settings.compounds_df)

#chemnet files
print '!!! loading chemnet files'
with open(my_settings.chemnet_pickle, 'r') as fid:
    cpd_group_data = pickle.load(fid)
cpd_df = pd.DataFrame(cpd_group_data, columns=['inchikey'])
cpd_group_lookup = cpd_group_data

# load the MST chemical network
with open(my_settings.mst_path, 'r') as f:
    net = pickle.load(f)

# regex expression to test inchikey input
# Character 9 of the second block must be "S" (standard inchi)
# Character 10 of the second block must be "A" (version 1)
inchikey_test = '[A-Z]{14}-[A-Z]{8}SA-[A-Z]{1}'

# set up joiner tables
a = mrs_reaction.index.values
b = mrs_reaction['refseq_id'].values
to_blowup = zip(a,b)
data = []
for pair in to_blowup:
    rxn_id = pair[0]
    refseqs = pair[1]
    for r in refseqs.split('|'):
        if r != '':
            data.append([rxn_id, r])

rxn_refseq_join = pd.DataFrame(data, columns=['reaction_id', 'refseq_id'])

print '!!! All databases loaded into memory'
# setup junk has been loaded
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def reformat_pactolus(df, original_compound=None, compound_score=None):
    """
    Reformats pactolus output table to be direcly portable into MAGI.
    1.  changes column "inchi_key_y" or "inchi_key" to
        "original_compound"
    2.  changes column "score" to "compound_score"

    Inputs
    ------
    df: Pactolus output table as Pandas DataFrame
    original_compound: column name corresponding to the original
                       compound inchi key. If None, will look for
                       matches to both "inchi_key_y" or "inchi_key".
    compound_score: column name corresponding to the pactolus score. If
                    None, will look for match to "score"

    Outputs:
    df: reformatted Pactolus table
    """

    cols = df.columns.values
    if original_compound is None:
        old_cols = pd.Series(['inchi_key_y', 'inchi_key'])
        tmp = old_cols[old_cols.isin(cols)].values
        if len(tmp) == 0:
            raise RuntimeError('no columns named "inchi_key_y" or "inchi_key" \
                in Pactolus table to rename')
        elif len(tmp) > 1:
            raise RuntimeError('both "inchi_key_y" and "inchi_key" found in \
                Pactolus table columns')
        else:
            original_compound = tmp[0]

    i = np.argwhere(cols == original_compound)
    if len(i) == 0:
        raise RuntimeError('There is no column named "%s", please double \
            check your input.' % (original_compound))
    else:
        i = i[0][0]
    cols[i] = 'original_compound'
    
    if compound_score is None:
        compound_score = 'score'
    
    i = np.argwhere(cols == compound_score)
    if len(i) == 0:
        raise RuntimeError('There is no column named "%s", please double \
            check your input.' % (original_compound))
    else:
        i = i[0][0]
    cols[i] = 'compound_score'
    df.columns = cols
    return df

def ec_parse(x):
    """
    Cleans up IMG gene table's Enzyme column by pulling out and
    separating EC numbers. Can be used as a pandas.apply() function.

    Inputs
    ------
    x: a string

    Outputs
    -------
    out: a string of EC numbers separated by a "|" character, or only the
         "|" character if the input is null or empty
    """

    if pd.isnull(x):
        return '|'

    out = '|'
    if x != '':
        if '<<>>' in x:
            one = x.split('<<>>')
            for two in one:
                try:
                    ec = two.split('EC:')[1].split('=')[0]
                except:
                    ec = ''
                out += ec+'|'
        else:
            try:
                ec = x.split('EC:')[1].split('=')[0]
            except:
                ec = ''
            out += ec+'|'
        return out
    else:
        return out

def load_genome(fasta, MAGI_PATH, annotation_file=None):
    """
    This function will take some standard fasta a input and convert it
    to an appropriate genome dataframe.
    This function will also construct a blast database for the genome.

    Inputs
    ------
    fasta: path to standard fasta file of all the genes.
           FASTA sequence header format should be:
           >UNIQUE_GENE_IDENTIFIER OTHER_INFORMATION
           Note the space between the unique identifier and other info

    MAGI_PATH: path to the root directory of MAGI Data. If None, the
               BLAST database is not made, and the gene table is not
               stored.

    annotation_file: path to an optional .TAB file that contains
                     annotations for the genes in the fasta file. One
                     column must be named "Gene_ID" and correspond to
                     the UNIQUE_GENE_IDENTIFIER in the fasta headers.

    Outputs
    -------
    genome: gene sequence table that is merged with the annotation table
            if one is provided. This table is saved as a pickle file to
            MAGI_PATH/gene_fastas/filename.pkl
    db_path: path to the genome's BLAST database
    """
    # TODO: find a way to handle windows text files (\n\r for new lines)
    # TODO: find a way to handle additional '>' characters in the header

    # process the fasta file
    with open(fasta, 'r') as f:
        genes = f.read()
    if genes.strip()[0] != '>':
        raise RuntimeError('%s does not appear to be a FASTA file' % (fasta))

    data = []
    for entry in genes.split('>')[1:]:
        h = '>' + entry.split('\n')[0]
        img = h.split(' ')[0][1:]
        s = ''.join(entry.split('\n')[1:])
        data.append([img, h, s])
    genome = pd.DataFrame(data, columns=['Gene_ID', 'header', 'sequence'])
    if genome['Gene_ID'].duplicated().any():
        first_dup = genome[genome['Gene_ID'].duplicated()].head(1)
        first_header = first_dup.iloc[0, 1]
        first_identifier = first_dup.iloc[0, 0]
        raise RuntimeError('There are duplicated Gene_ID fields! please check\
            your unique gene identifiers in the fasta headers. For the first\
            FASTA sequence of %s, I parsed the unique gene identifier to be \
            %s' % (repr(first_header), repr(first_identifier)))

    if annotation_file is not None:
        # Make gene info table
        annotation_table = load_dataframe(annotation_file)
        # turn gene IDs into strings to keep things consistent
        annotation_table['Gene_ID'] = annotation_table['Gene_ID'].apply(str)

        # clean up the nans
        annotation_table.fillna(value='', inplace=True)

        # Really only if the annotation file was an IMG-derived table
        try:
            annotation_table['EC'] = annotation_table['Enzyme'].apply(ec_parse)

            # move the new EC column to front for ease of visualization
            cols = annotation_table.columns.tolist()
            newcols = cols[-1:] + cols[:-1]
            annotation_table = annotation_table[newcols]
        except KeyError:
            print 'Could not find a column corresponding to EC annotations, \
                    skipping EC parsing'

        genome = pd.merge(genome, annotation_table, on='Gene_ID', how='left')

    print '!@# FASTA file loaded with', len(genome), 'genes'

    genome.set_index('Gene_ID', inplace=True, drop=True)
    genome_name = fasta.split('/')[-1].split('.')[0]
    if MAGI_PATH is not None:
        gene_fasta_path = os.path.join(MAGI_PATH, 'gene_fastas')
        if not os.path.isdir(gene_fasta_path):
            os.makedirs(gene_fasta_path)
        gene_seq_path = os.path.join(gene_fasta_path,
            '%s_sequences.pkl' % (genome_name))
        # gene_seq_path = '%s/gene_fastas/%s_sequences.pkl' \
        #                 % (MAGI_PATH, genome_name)
        genome.to_pickle(gene_seq_path)
        print '!!! saved gene_sequence table here:', gene_seq_path

        # Make the blast database of the fasta if it doesn't already exist
        # this command makes the blast database for a given fasta file.
        makeblastdb_path = '%s/makeblastdb' % (blastbin)
        fasta_path = fasta
        db_path = MAGI_PATH + '/BLAST_dbs/' + fasta_path.split('/')[-1].split(
            '.')[0] + '.db'
        print '!!! blast database stored here:', db_path
        make_db_command = '%s -in %s -out %s -dbtype prot' \
            % (makeblastdb_path, fasta_path, db_path)
        if not os.path.isfile(db_path + '.pin'):
            blastp = subprocess.Popen(
                make_db_command,
                shell=True, stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if blastp.stdout is not None:
                print blastp.stdout.read()
            if blastp.stderr is not None:
                print blastp.stderr.read()
    else:
        db_path = None
    return genome, db_path

def partition_indexes(totalsize, numberofpartitions, offset=0):
    """
    Used to split up an iterable into multiple chunks
    This function returns the indices to create chunks of roughly the
    same size

    From:
    http://stackoverflow.com/questions/33878939/
        get-indices-of-roughly-equal-sized-chunks

    Inputs
    ------
    totalsize: length of the iterable to partition
    numberofpartitions: number of partitions to split the iterable into
    offset: number to add to the partition indices. Used when you're
            partitioning a slice of something, for example if you want
            to parititon array[100:200] into 4 partitions, this function
            would return [[0, 25], [25, 50], [50, 75], [75, 100]] if
            offset were kept default. So you would set offset to 100,
            which would return
            [[100, 125], [125, 150], [150, 175], [175, 200]]
    """

    chunk_indices = []
    # Compute the chunk size (integer division; i.e. assuming Python 2.7)
    chunksize = totalsize / numberofpartitions
    # How many chunks need an extra 1 added to the size?
    remainder = totalsize - chunksize * numberofpartitions
    a = 0
    for i in xrange(numberofpartitions):
        b = a + chunksize + (i < remainder)
        # Yield the inclusive-inclusive range
#         chunk_indices.append([a, b - 1])
        chunk_indices.append([a, b])
        a = b

    offset_chunk_indices = []
    for indices in chunk_indices:
        i = indices[0] + offset
        j = indices[1] + offset
        offset_chunk_indices.append([i, j])

    return offset_chunk_indices

def tabulate_blast(results):
    """
    Convert a blast results text output into a dataframe

    Inputs
    ------
    results: BLAST results output when run with the following parameter:
             -outfmt \'10 qacc sacc qcovs length ppos evalue bitscore\'

    Outputs
    -------
    df: Pandas dataframe summarizing the BLAST results
    """

    cols = ['query acc.', 'subject acc.', '% query coverage per subject',
            'alignment length', '% positives', 'evalue', 'bit score']
    result_list = results.split('\n')[:-1]
    data = []
    for r in result_list:
        elements = r.split(',')
        data.append(elements)
    df = pd.DataFrame(data, columns=cols)
    return df

def multi_blast(query_list, query_full_table, database_path, result_path,
                cpu=2, raise_blast_error=True):
    """
    Large-scale BLASTP with multiple processes.
    Partitions the query list into N lists, then uses a shell script to
    open several BLASTP processes, saves them to disk and then collects
    them.

    Inputs
    ------
    query_list: list of indices to query_full_table that you want blasted
    query_full_table:   pandas DataFrame representing protein sequences
                        Must have a column named "sequence"
    database: path to blast database to be used
    result_path: where you want results saved
    cpu: number of processes to open
    scriptpath: path to shell script that opens multiple blast processes
    raise_blast_error: boolean, False will still print BLAST errors/warnings,
                        but will not raise an error

    Outputs
    -------
    blast_results: pandas DataFrame that stores all blast results.
    """
    db_name = database_path.split('/')[-1]

    # set up the blast search strings
    idxes = partition_indexes(len(query_list), cpu)
    mplist = []
    for pair in idxes:
        i, j = pair
        seq = ''
        for gid in query_list[i:j]:
            if gid in query_full_table.index:
                seq += '>' + gid + '\n'
                seq += query_full_table.loc[gid, 'sequence']
                seq += '\n'
            else:
                print '%s not in sequence database!' % (gid)
        mplist.append(seq)

    # call a shell script that opens n blasts,
    # instead of opening pool of workers
    # first need to save the input seq as a temporary file
    cwd = result_path + '/' + 'multi_blast_files'
    if not os.path.isdir(cwd):
        os.makedirs(cwd)
    for i, seq in enumerate(mplist):
        with open('%s/tmp_seq_%s.faa' % (cwd, i), 'w') as f:
            f.write(seq)

    # then call the job
    scriptpath = '%s/recip_blaster.sh' % (blastbin)
    print '!!! blast script:', scriptpath
    print '!!! # processes to open:', cpu
    print '!!! database path:', database_path
    print '!!! results stored:', cwd
    sys.stdout.flush()
    subprocess.call(
        '%s %s %s %s %s 2> %s/blasterr__%s.txt'
        % (scriptpath, cpu, database_path, cwd, my_settings.repo_location,
            result_path, db_name),
        shell=True)
    with open('%s/blasterr__%s.txt' % (result_path, db_name), 'r') as f:
        msg = f.read()
    if msg != '':
        print '!@# WARNING: BLAST script error message:'
        print msg
        print '-'*80
        if raise_blast_error:
            raise RuntimeError('blast had an error message')
    sys.stdout.flush()

    # collect the results
    results = ''
    for i, seq in enumerate(mplist):
        with open('%s/tmp_out_blasted_%s.txt' % (cwd, i), 'r') as f:
            results += f.read()
    results_table = tabulate_blast(results)

    # convert e-value to blast score
    results_table['e_score'] = results_table['evalue'].apply(float).apply(
        np.log10) * -1
    # convert perfect blast scores to 200
    idx = results_table[results_table['e_score'] == np.inf].index
    results_table.loc[idx, 'e_score'] = 200.

    # clean up temp folders
    rmtree(cwd)

    return results_table

def keep_top_blast(escore_series, filt=0.85):
    """
    Function to keep only the top scores of a blast result.
    Meant to be an "apply" function for a pandas groupby object.
    Can also be used as a regular "apply" function.

    Ex. Usage:
    gene_groups = df.groupby('query acc.')
    multidx = gene_groups['e_score'].apply(keep_top_blast).index
    idx = multidx.levels[1]
    df.loc[idx] # this is now the filtered df

    Inputs
    ------
    escore_series: Pandas Series
    filt: lower cutoff of max value to keep

    Outputs
    -------
    final: Pandas Series
    """
    max_score = escore_series.max()
    final = escore_series[escore_series >= (max_score * filt)]
    return final

def tautomer_finder(compound_mol, result='split', raise_errors=False):
    """
    enumerates tautomers of a given compound

    MolVS tautomer enumerator flattens compounds, so the default only
    returns the first block of the inchikey

    Inputs
    ------
    compound_mol: RdKit Mol object of the compound to tautomerize
    result: "split" returns the first block of an inchikey
            "full" returns the full inchikey
            "smiles" returns the smiles string
            "inchi" returns the InChI string 
            "mol" returns RdKit Mol object
    reaise_errors: when False, it still prints the error as a warning

    Outputs
    -------
    tautomer_list: list of unique tautomers, in the output format
                   specified by result argument
    """

    if not isinstance(compound_mol, type(Chem.Mol())):
        raise RuntimeError('The input is not an rdkit Mol object!')

    if not isinstance(result, str):
        raise RuntimeError('"result" arg must be a string!')
    if result.lower() not in ['split', 'full', 'smiles', 'inchi', 'mol']:
        raise RuntimeError('"result" arg must be either "split", "full", \
            "smiles", "inchi", or "mol"')

    compound_smiles = Chem.MolToSmiles(compound_mol, isomericSmiles=True)

    # some compounds break the tautomerizer
    try:
        enumerated_tautomers = enumerate_tautomers_smiles(compound_smiles)
    except TypeError, e:
        if e.message == 'tuple indices must be integers, not NoneType':
            enumerated_tautomers = []
    except Exception as e:
        if raise_errors is False:
            print 'WARNING: %s could not be tautomerized; %s' \
                    % (Chem.InchiToInchiKey(Chem.MolToInchi(compound_mol)),
                        e.args)
            enumerated_tautomers = [compound_smiles]
        else:
            raise

    # convert smiles into inchikeys
    tautomer_list = []
    for ts in enumerated_tautomers:
        if result.lower() in ['split', 'full']:
            i = Chem.InchiToInchiKey(Chem.MolToInchi(Chem.MolFromSmiles(ts)))
            if result.lower() == 'split':
                tautomer_list.append(i.split('-')[0])
            elif result.lower() == 'full':
                tautomer_list.append(i)
            else:
                raise RuntimeError(
                    'could not determine what to do with %s' % (result))
        elif result.lower() == 'smiles':
            tautomer_list.append(ts)
        elif result.lower() == 'inchi':
            tautomer_list.append(Chem.MolToInchi(Chem.MolFromSmiles(ts)))
        elif result.lower() == 'mol':
            tautomer_list.append(Chem.MolFromSmiles(ts))
        else:
            raise RuntimeError(
                'could not determine what to do with %s' % (result))

    return list(set(tautomer_list))

def neighbor_finder(inchikey, chemical_network=net, cpd_group=None, level=2):
    """
    finds neighbors of a given compound to the given level in the
    chemical network

    Inputs
    ------
    inchikey:           Full or partial inchikey for the compound to be
                        searched
    chemical_network:   MST chemical network
    cpd_group:          if the compound group of inchikey is already
                        known, this will override the cpd_group
                        neighbor_finder
    level:              to what level to search the network

    Outputs
    -------
    returns a list of tuples describing the network level and inchikeys:
    (level, [list of inchikeys at that level])
    """

    if cpd_group is None:
        # find the compound group/node number
        # flatten the inchikey
        inchikey_query = inchikey.split('-')[0]

        # find the compound group of the inchikey
        # there should only be one!!! If there are multiple, this will
        # only take the first group found
        for idx, istring in enumerate(cpd_group_lookup):
            if inchikey_query in istring:
                cpd_group = idx
                break
    neighbor_groups = []
    if cpd_group is not None:
        # compound group number is also the node ID; find the compound
        # groups that are neighbors in the MST
        neighbor_node_dict = nx.single_source_shortest_path_length(
            net, cpd_group, level)
        transformed_neighbor_node_dict = {}
        for k in list(set(neighbor_node_dict.values())):
            transformed_neighbor_node_dict[k] = []
        for k, v in neighbor_node_dict.iteritems():
            if v != 0:
                transformed_neighbor_node_dict[v].append(k)
        for level, node_list in transformed_neighbor_node_dict.iteritems():
            neighbor_inchikey_list = []
            for neighbor_cpd_group in node_list:
                neighbor_inchikey_str = cpd_group_lookup[neighbor_cpd_group]
                for neighbor_inchikey in neighbor_inchikey_str.split('///'):
                    neighbor_inchikey_list.append(neighbor_inchikey)
            neighbor_groups.append((level, neighbor_inchikey_list))

    else:
        print '695 WARNING: Could not find "%s" in the chemical network' \
              % (inchikey)

    return neighbor_groups

def find_reactions_of_compound(inchikey, rxn_db=mrs_reaction, 
                               compound_col='allcpd_ikeys'):
    """
    Given an inchikey input, find reactions that contain the compound

    Inputs
    ------
    inchikey: can be one, two, or three block inchi
    rxn_db: reaction database table
    compound_col: column in rxn_db that contains a list of inchikeys of
                  all compounds involved in that reaction. Elements in
                  this column must be a string for now (e.g. not a list)

    Outputs
    -------
    A list of reaction indices, or None if the inchikey was not found in
    any reactions.
    """
    # TODO: use join tables instead of .str.contains() method

    reactions = rxn_db[rxn_db[compound_col].str.contains(
        inchikey)]
    if reactions.shape[0] != 0:
        reaction_idx_list = reactions.index.tolist()
        return reaction_idx_list
    else:
        return None

def enumerate_compound_results(original_compound, compound_results,
                               reaction_idx_list, level=0,
                               neighbor='', note=''):
    """
    Inputs
    ------
    original_compound: initial compound that began the search
    compound_results:   dictionary with the following keys:
                        'original_compound', 'level', 'neighbor',
                        'reaction_id', 'note' where each key-value is a
                        list
    reaction_idx_list:  the ouptut of find_reaction_of_compound()
    level:              The chemical network search level
    neighbor:           The chemical network neighbor of
                        original_compound

    Outputs:
    --------
    An updated compound_results dictionary
    """

    if reaction_idx_list is not None:
        for rxn_idx in reaction_idx_list:
            compound_results['original_compound'].append(original_compound)
            compound_results['level'].append(level)
            compound_results['neighbor'].append(neighbor)
            compound_results['reaction_id'].append(rxn_idx)
            compound_results['note'].append(note)
    else:
        compound_results['original_compound'].append(original_compound)
        compound_results['level'].append(level)
        compound_results['neighbor'].append(neighbor)
        compound_results['reaction_id'].append('')
        compound_results['note'].append(note)

    return compound_results

def mol_from_inchikey(inchikey):
    """
    Returns an RdKit Mol object from an inchikey input via the compound
    DataFrame

    Inputs
    ------
    inchikey: a standard InChI key

    Outputs
    -------
    Rdkit Mol object or None, if the inchikey input was not found in the
    compound DataFrame
    """

    matched_inchis = compounds[compounds['inchi_key'].str.contains(
        inchikey, regex=False)]['inchi'].values
    if matched_inchis.any():
        # in case there are duplicates, just take the first one
        inchi = str(matched_inchis[0])
        compound_mol = Chem.MolFromInchi(inchi)
        return compound_mol
    else:
        return None

def connect_compound_to_reaction(inchikey, tautomer=True, neighbor_level=2):
    """
    Connects a compound, its neighbors, and optionally its and its
    neighbors' tautomers to a reaction

    replaces part of old searcher_thorough() function and incorporates
    neighbor_searcher

    Inputs
    ------
    inchikey:       compound's inchikey
    mrs_reaction:   reaction database dataframe
    tautomer:       boolean; whether or not to find tautomers of input
                    compound (and it's neighbors)
    neighbor_level: to what level to search the chemical network
                    0: does not search the network
                    1: finds the immediate neighbors
                    2: finds the neighbor's neighbors
                    etc.

    Outputs
    -------
    a dataframe describing the compound and what reactions it connected
    through, through which neighbors if applicable. The final DataFrame
    is sorted by direct/flat tautomer, and then
    original_compound-reaction duplicates are dropped.

    """

    # test validity of input
    if not re.match(inchikey_test, inchikey):
        raise RuntimeError('"%s" is not a valid inchikey or is not the correct\
         version of inchi : xxxxxxxxxxxxxx-xxxxxxxxSA-X' % (inchikey))

    # convert inchikey into two-block
    search_inchikey = '-'.join(inchikey.split('-')[:2])

    # find any direct matches
    direct_reaction_idx_list = find_reactions_of_compound(search_inchikey)

    # initialize the results dict
    compound_results = {
        'original_compound': [],
        'level': [],
        'neighbor': [],
        'reaction_id': [],
        'note': []
        }

    compound_results = enumerate_compound_results(
        inchikey, compound_results, direct_reaction_idx_list,
        level=0, neighbor='', note='direct')

    if tautomer:
        # make an rdkit mol of the compound
        compound_mol = mol_from_inchikey(inchikey)
        if compound_mol is None:
            print 'WARNING: Could not find "%s" in the compound database; \
                skipping tautomer and neighbor searching' % (inchikey)
            return pd.DataFrame(compound_results)

        # get tautomers of the compound
        tautomer_list = tautomer_finder(compound_mol)

        # find reactions for the tautomers
        tautomer_search_pattern = '|'.join(tautomer_list)
        tautomer_reaction_idx_list = find_reactions_of_compound(
            tautomer_search_pattern)
        compound_results = enumerate_compound_results(
            inchikey, compound_results, tautomer_reaction_idx_list, level=0,
            neighbor='', note='flat tautomer')

    # find neighbors
    if neighbor_level != 0:
        neighbor_groups = neighbor_finder(inchikey, level=neighbor_level)

        # next look for reaction matches to neighbors
        for level, neighbor_compound_list in neighbor_groups:
            for neighbor_inchikey in neighbor_compound_list:
                # first get the direct matches
                reaction_idx_list = find_reactions_of_compound(
                    neighbor_inchikey)
                compound_results = enumerate_compound_results(
                    inchikey, compound_results, reaction_idx_list,
                    level=level, neighbor=neighbor_inchikey, note='direct')
                if tautomer:
                    # then the flat tautomer searches
                    neighbor_mol = mol_from_inchikey(neighbor_inchikey)
                    if neighbor_mol is not None:
                        tautomer_list = tautomer_finder(neighbor_mol)
                        tautomer_search_pattern = '|'.join(tautomer_list)
                        tautomer_reaction_idx_list = \
                            find_reactions_of_compound(tautomer_search_pattern)
                        compound_results = enumerate_compound_results(
                            inchikey, compound_results,
                            tautomer_reaction_idx_list, level=level,
                            neighbor=neighbor_inchikey, note='flat tautomer')

    compound_reaction_result_df = pd.DataFrame(compound_results)
    # now perform cleanup on the df
    # first sort the df so that all the direct hits come first
    compound_reaction_result_df.sort_values(
        ['note', 'level'], ascending=[True, True], inplace=True)
    # then drop any duplicated compound-reaction associations, keeping
    # only the direct hits
    compound_reaction_result_df.drop_duplicates(
        ['original_compound', 'reaction_id'], inplace=True)
    return compound_reaction_result_df

def refseq_to_reactions(blast_results, refseq_col):
    """
    enumerates the reaction(s) that a refseq is associated with

    Inputs
    ------
    blast_results: output of multi_blast()
    refseq_col: column name that has the reaction reference sequences
                either "query acc." or "subject acc."
                "query acc." for reaction_to_gene
                "subject acc." for gene_to_reaction

    Outputs
    -------
    result_table: pandas dataframe that connects BLAST results to a
                  reaction ID in the reaction database
    """

    if refseq_col not in blast_results.columns:
        raise RuntimeError('%s is not a column in the dataframe!'
                           % (refseq_col))

    result_table = pd.merge(
        blast_results, rxn_refseq_join, 
        left_on=refseq_col, right_on='refseq_id', how='left')
    result_table.drop('refseq_id', axis=1, inplace=True)

    # keep only the best hit to each reaction, for each query
    result_table.sort_values(
        ['query acc.', 'e_score'],
        ascending=[True, False]).drop_duplicates(['query acc.', 'reaction_id'])

    return result_table

def reciprocal_agreement(df, forward_name='reaction_id_r2g',
                         reverse_name='reaction_id_g2r',
                         closeness_threshold=0.75):
    """
    Finds and scores reciprocal agreement converging on a reaction
    between forward and reverse homology searches.
    
    Inputs
    ------
    df: dataframe containing forward and reverse homology searches
    forward_name: column name in df corresponding to the reaction ID for
                  the forward search
    reverse_name: column name in df corresponding to the reaction ID for
                  the reverse search
    closeness_threshold: Determines cutoff to call a reciprocal search
                         as "close" if the reaction IDs between forward
                         and reverse searches do not agree.
                         For example, if there is no agreement but the
                         forward score is 100 and reverse score is 90,
                         this would be scored as "close" if the
                         closeness_threshold was less than or equal to
                         0.9

    Outputs
    -------
    df: input df but with a new column named "reciprocal_score"
        corresponding to scored reciprocal agreement:
        2.0: agreement between forward and reverse searches
        1.0: forward and reverse searches disagree on reaction, but the
             homology scores are close as judged by closeness_threshold
        0.1: either a forward or reverse search does not exist (homology
             score below threshold)
        0.01: forward and reverse searches disagree on reaction and are
              not close
    """
    # find agreement
    agree_idx = df[df['reaction_id_r2g'] == df['reaction_id_g2r']].index
    df.loc[agree_idx, 'reciprocal_score'] = 2.
    # disagreement
    disagree = df[df['reaction_id_r2g'] != df['reaction_id_g2r']].index
    slc = df.loc[disagree]
    # close disagreements get a medium score
    close = (
        slc[['e_score_r2g', 'e_score_g2r']].min(axis=1)
        >= (slc[['e_score_r2g', 'e_score_g2r']].max(axis=1) * closeness_threshold))
    close_idx = slc.loc[close].index
    df.loc[close_idx, 'reciprocal_score'] = 1
    # very different disagreements get a low score
    wrong_idx = df[pd.isnull(df['reciprocal_score'])].index
    df.loc[wrong_idx, 'reciprocal_score'] = 0.01
    # if one direction did not get a blast score, 
    # change reciprocal score to 0.1 - not wrong, but not close either.
    incomparable_idx = df.loc[pd.isnull(df[['e_score_r2g',
        'e_score_g2r']]).any(axis=1)].index
    df.loc[incomparable_idx, 'reciprocal_score'] = 0.1

    return df

def homology_score(df, forward_name='e_score_r2g', reverse_name='e_score_g2r'):
    """
    Calculates homology score for a compound-gene association:
    Takes the sum of forward and reverse scores, subtracts the difference
    between them: (forward + reverse) - abs(forward-reverse)

    deprecated:
    average(forward_score, reverse_score) - abs(forward_score - reverse_score)

    Inputs
    ------
    df: dataframe containing forward and reverse homology scores
    forward_name: column name in df corresponding to the forward score
    reverse_name: column name in df corresponding to the reverse score

    Output
    ------
    array of combined homology scores
    """

    forward = df[forward_name].values.astype(float)
    reverse = df[reverse_name].values.astype(float)

    score = (forward + reverse) - abs(forward-reverse)

    return score

def magi_score(s, weights=None):
    """
    Calculates the geometric mean of an array of numbers.
    Used to calculate one score even if sub-scores are scaled differently.

    Inputs
    ------
    s: 1D array-like list of scores
    weights: 1D array-like list of weights

    Outputs
    -------
    The (weighted) geometric mean
    """

    if not isinstance(s, np.ndarray):
        s = np.asarray(s)
    s = s.astype(float)
    if weights:
        if not isinstance(weights, np.ndarray):
            weights = np.asarray(weights)
        weights = weights.astype(float)
    else:
        weights = np.ones(s.shape)

    return np.exp(np.sum(np.multiply(weights, np.log(s))) / weights.sum())

def mass_from_inchikey(inchikey_list, compound_db=compounds,
                       inchi_key_col='inchi_key',
                       mass_col='mono_isotopic_molecular_weight'):
    """
    Given a list of inchikeys, returns a dataframe with their
    monoisotopic molecular weights

    Inputs
    ------
    inchikey_list: array-like list of inchikeys
    compound_db: compound database
    inchi_key_col: name of column corresponding to inchi_key in
                   compound_db
    mass_col: name of column corresponding to monoisotopic molecular
              weight in compound_db

    Outputs
    -------
    df: dataframe with inchi_key and monoisotopic_mass as columns
    """
    df = compounds.set_index(inchi_key_col).loc[inchikey_list]
    df = df.reset_index()[[inchi_key_col, mass_col]]
    df.drop_duplicates(inchi_key_col, inplace=True)
    return df

def ppm_window(mass, ppm=5, result='bounds'):
    """
    Given a mass and a ppm error, returns lower and upper bounds
    corresponding to that ppm window.

    Inputs
    ------
    mass:   monoisotopic mass
    ppm:    the ppm error to construct the window
    result: "bounds" returns a list of lower and upper bounds
            "error" returns the amu value of the error, if you wanted to
            add and/or subtract the error as you wish

    Outputs
    -------
    Either a list of lower and upper mass bounds, or a single value
    corresponding to the amu error
    """
    error = ppm/1e6 * mass
    lower_bound = mass - error
    upper_bound = mass + error
    if result.lower() == 'bounds':
        return [lower_bound, upper_bound]
    elif result.lower() == 'error':
        return error
    else:
        raise RuntimeError(
            '%s is not a valid result argument' %(result)
            )

def ppm_error(mass, theoretical_mass):
    """
    Returns the ppm error of a given observed mass, and theoretical mass
    """
    if not isinstance(theoretical_mass, float):
        theoretical_mass = float(theoretical_mass)
    ppm = (mass - theoretical_mass) / theoretical_mass * 1e6
    return abs(ppm)

def accurate_mass_search(mass_list, ppm=10, compound_db=compounds):
    """
    Given an input list of monoisotopic masses, finds all compounds in
    the database matching within the given ppm error.
    
    Inputs
    ------
    mass_list: array like list of monoisotopic masses
    ppm: ppm error cuttof for accurate mass search
    
    Outputs
    -------
    df: dataframe with columns query_mass, target_mass, ppm,
        original_compound, compound_score corresponding to the mass_list
        masses, found masses, ppm difference, inchikeys for those
        compounds, and a function of the ppm error as a compound score

    compound_score = ppm_cutoff + 1 - ppm_difference
    """

    data = []
    for m in mass_list:
        slc = compound_db[abs(compound_db['mono_isotopic_molecular_weight']\
            - m) <= ppm_window(m, ppm=ppm, result='error')]
        vals = slc[['mono_isotopic_molecular_weight', 'inchi_key']].values
        if len(vals) > 0:
            for v in vals:
                m2 = v[0]
                inchikey = v[1]
                p = ppm_error(m, m2)
                s = ppm + 1. - p
                data.append([m, m2,p, inchikey, s])
        else:
            data.append([m, pd.np.nan, pd.np.nan, 'None', pd.np.nan])
    df = pd.DataFrame(
        data,
        columns=['query_mass', 'target_mass', 'ppm',
                 'original_compound', 'compound_score'])
    return df