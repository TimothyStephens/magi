import os
import pandas as pd
import numpy as np
import sys
import subprocess
import warnings
import datetime
import argparse
import time
import json
from multiprocessing import cpu_count as counting_cpus
import csv

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from local_settings import local_settings_magi2 as settings_loc

def is_existing_file(filepath):
    """
    Checks if a file exists and return absolute path if it exists

    Inputs
    ------
    filepath: a relative or absolute path to a file.

    Outputs
    -------
    An error is raised if the path does not exist and the absolute path is returned if the file exists.
    """
    if not os.path.exists(filepath):
        msg = "{0} does not exist".format(filepath)
        raise argparse.ArgumentTypeError(msg)
    else:
        return os.path.abspath(filepath)
    
def is_database(db_path):
    """
    Checks if the BLAST genome database exists.

    Inputs
    ------
    db_path: a relative or absolute path to a BLASTP database. The path should end with .db
            In this folder, a .db.phr file, a .db.pin file and a .db.psq file should be present.

    Outputs
    -------
    An error is raised if the database does not exist and the absolute path is returned if the database exists.
    """

    for file_extension in [".phr", ".pin", ".psq"]:
        db_file = db_path + file_extension
        if not os.path.exists(db_file):
            msg = "{0} does not exist".format(db_file)
            raise argparse.ArgumentTypeError(msg)
    return os.path.abspath(db_path)

def percentage_values_to_decimal(percentage):
    """Turns the blast filter and reciprocal closeness percentages 
    into decimal numbers"""
    try:
        percentage = int(percentage)
    except:
        msg = "Please enter an integer value"
        raise argparse.ArgumentTypeError(msg)        
    if percentage > 100:
        msg = "Max value is 100"
        raise argparse.ArgumentTypeError(msg)
    elif percentage < 0:
        msg = "Value cannot be negative"
        raise argparse.ArgumentTypeError(msg)
    else:
        decimal = percentage/100.
    return decimal

def positive_number(number):
    """Checks if a number is positive"""
    try:
        number = float(number)
    except:
        msg = "Please enter a numeric value"
        raise argparse.ArgumentTypeError(msg)        
    if number < 0:
        msg = "Value cannot be negative"
        raise argparse.ArgumentTypeError(msg)
    else:
        return number

def str2bool(value):
    """
    source: adapted from https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    """
    if isinstance(value, bool):
       return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_arguments():
    """
    This is the MAGI argument parser that is used in all workflows. 
    It checks if all required arguments are passed, if numbers fall within MAGI-approved ranges and if files exist.

    Outputs
    -------
    An argparse.args object containing all arguments. 
    """
    def set_cpu_count(cpu_count):
        max_cpu = counting_cpus()  
        if cpu_count == 0:
            cpu_count = max_cpu    
        if cpu_count > max_cpu:
            msg = "ERROR: You have exceeded the cpus on this machine ({})".format(max_cpu)
            raise argparse.ArgumentTypeError(msg)
        return cpu_count

    """parse arguments"""
    parser = argparse.ArgumentParser()
    # required arguments
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-f', '--fasta', type=is_existing_file,
        help='path to fasta file of genes in sample')
    required_args.add_argument('-c', '--compounds', type=is_existing_file,
        help='path to observed compounds file')
    required_args.add_argument('--jsonfile', type=is_existing_file,
        help='path to .json file with magi input parameters. Use instead of any arguments here. MAGI is extremely sensitive to the keys and values in this file.')
    
    # jump-start the script after certain computations
    start_halfway_args = parser.add_argument_group('Arguments to jump-start the script after certain computations')
    start_halfway_args.add_argument('--gene_to_reaction', type=is_existing_file,
        help='path to gene_to_reaction file, must be in pickle format')
    start_halfway_args.add_argument('--compound_to_reaction', type=is_existing_file,
        help='path to compound_to_reaction file, must be in pickle format')
    start_halfway_args.add_argument('--reaction_to_gene', type=is_existing_file,
        help='path to reaction_to_gene file, must be in pickle format')
    start_halfway_args.add_argument('--merged_before_score', type=is_existing_file,
        help='path to merged_before_score table, must be in hdf5 format,\
        with the key "merged_before_score"')
    start_halfway_args.add_argument('--genome_db', help = "path to genome .db files", type=is_database)

    # Use this if only a part of the workflow should be run
    stop_halfway_args = parser.add_argument_group('Arguments to run a part of the script')
    stop_halfway_args.add_argument('--gene_to_reaction_only',
        help="Use this parameter if you are only interested in the gene to reaction search", 
        action='store_true', default=False)
    stop_halfway_args.add_argument('--compound_to_reaction_only',
        help="Use this parameter if you are only interested in the compound to reaction search", 
        action='store_true', default=False)
    
    # Arguments to run multiple parts of the workflow sequentially
    pipeline_args = parser.add_argument_group('Arguments to control running MAGI scripts sequentially')
    pipeline_args.add_argument('--not_first_script', default = False, action='store_true', help="set this to true if this is not the first script in a sequential magi run.")
    pipeline_args.add_argument('--intermediate_files_dir', help="path to intermediate files directory")
    
    # optional runtime variables
    optional_args = parser.add_argument_group("Optional runtime variables")
    optional_args.add_argument('-a', '--annotations', type=is_existing_file,
        help='path to annotation file for genes in sample', 
        default=None)
    optional_args.add_argument('-n', '--cpu_count', 
        help='number of cpus to use for multiprocessing. Default is to use max!', 
        type=int, default=0)
    optional_args.add_argument('-o', '--output', 
        help='path to a custom output', 
        type=str)
    optional_args.add_argument('-l', '--level', 
        help='how many levels deep to search the chemical network', 
        type=int, choices=[0,1,2,3], default=2)
    optional_args.add_argument('--legacy', dest='legacy', action='store_true',
        help='use legacy tautomer searching; default is no')
    optional_args.add_argument('--no-legacy', dest='legacy', action='store_false',
        help='use precomputed compound-to-reaction; default is yes')
    optional_args.set_defaults(legacy=False)
    optional_args.add_argument('--mute', 
        help='mutes pandas warnings', 
        action='store_true')
    optional_args.add_argument('--pactolus', 
        help='Flag to tell MAGI that the compounds input is a pactolus file', 
        action='store_true')
    optional_args.add_argument('--test', 
        help='TBD: run MAGI only on the first # of pactolus compounds', 
        type=int)
    optional_args.add_argument('--debug', 
        help='TBD: prints a lot of info', 
        action='store_true')
    optional_args.add_argument('--blast_filter', 
        help='How stringent to filter the top BLAST results, as percent;\
        default is 85 meaning that only BLAST results within 85%% of the top\
        result will be taken.', 
        type=percentage_values_to_decimal, default=0.85)
    optional_args.add_argument('--reciprocal_closeness', 
        help='Cutoff to call a reciprocal disagreement as "close", as percent;\
        default is 75 meaning that a reciprocal disagreement will be classified\
        as "close" if the lower blast score (e score) is within 75%% of the higher\
        score', 
        type=percentage_values_to_decimal, default=0.75)
    optional_args.add_argument('--final_weights', 
        help='Defined weights to weight the final scoring for the scores:\
        compound_score, similarity, diameter, reciprocal_score homology_score reaction_connection', 
        type=positive_number, nargs=6, default=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    optional_args.add_argument('--chemnet_penalty', 
        help='Base factor in the chemical network search level penalty', 
        type=positive_number, default=4)
    optional_args.add_argument('--intermediate_files',
        help='What directory within --output to store intermediate files',
        type=str, default='intermediate_files')
    
    # Parameters for accurate mass search
    mass_search_args = parser.add_argument_group("Arguments for the optional accurate mass search")
    mass_search_args.add_argument('--is_mass_search',
                        action='store_true', default=False,
                        help = "Set this parameter if m/z values need to be transformed into candidate compounds."
                        )
    mass_search_args.add_argument('--polarity',
                        type=str, choices=['pos','neg','neut'], default=None,
                        help = "Specify if the masses were measured in negative mode, positive mode or if they have been transformed to neutral masses."
                        )
    mass_search_args.add_argument('--adducts_pos',
                        type=str, default=None,
                        help="Specify which positive adducts to investigate if the polarity is positive."
                        )
    mass_search_args.add_argument('--adducts_neg',
                        type=str, default=None,
                        help="Specify which positive adducts to investigate if the polarity is negative."
                        )
    mass_search_args.add_argument('--ppm', 
        help='The ppm cutoff for the accurate mass search. Default is 10 ppm.', 
        type=int, default=10)

    # MAGI c2r 2.0 parameters
    c2r_args = parser.add_argument_group("Arguments for the optional accurate mass search")
    c2r_args.add_argument('--diameter', 
        help="Minimum diameter to use for retro rules reactions", type = int, default = 12) # TODO: add check to be in range and even
    c2r_args.add_argument('--fingerprint', 
        help="fingerprint radius for Morgan molecular fingerprint", type = int, default = 3)
    c2r_args.add_argument('--similarity_cutoff', 
        help="Minimum similarity cutoff", type = float, default = 0.6)
    c2r_args.add_argument('--use_precomputed_reactions', 
        help="Use of previously computed reactions. Default is True. Very slow if set to False", type = str2bool, default = True)
    
    args = parser.parse_args()
    
    # Check parameters and set number of required CPUs
    if args.fasta is None and args.compounds is None and args.jsonfile is None and args.not_first_script is False:
        parser.error('ERROR: either FASTA or metabolites file is required, or a json file with input parameters.')
    args.cpu_count = set_cpu_count(args.cpu_count)
    return args

def make_output_dirs(output_dir=None, fasta_file=None, compounds_file=None, intermediate_files='intermediate_files'):
    """
    Set up where MAGI results will be stored. This creates an output directory, intermediate files directory and starts
    the overall MAGI program timer.

    Inputs
    ------
    output_dir:     a relative or absolute path to the output directory. 
                    This directory will be created if it does not yet exist.
    fasta_file:     if no output directory is specified, the output will be stored at the location of the fasta file, 
                    in a folder with the date and the fasta file name as its name.
    compounds_file: if no fasta file and output directory are specified, the output will be stored at the location of the fasta file, 
                    in a folder with the date and the compounds file name as its name.
    intermediate_files: What directory within the output directory to store intermediate files.

    Outputs
    -------
    output_dir:     absolute path to the output directory
    intermediate_files_dir: absolute path to the intermediate files directory
    """
    if output_dir is None:
        # autoname the directory based on fasta, or compound file
        # this will change eventually
        if fasta_file is not None:
            experiment_name = os.path.splitext(os.path.basename(fasta_file))[0]
            experiment_dir = os.path.abspath(os.path.dirname(fasta_file))
        else:
            experiment_name = os.path.splitext(os.path.basename(compounds_file))[0]
            experiment_dir = os.path.abspath(os.path.dirname(compounds_file))
        today = datetime.datetime.now()
        experiment_name += today.strftime('_%Y%m%d')
        experiment_path = os.path.join(experiment_dir, experiment_name)
    else:
        experiment_path = output_dir
    
    output_dir = os.path.abspath(experiment_path)
    intermediate_files_dir = os.path.join(output_dir, intermediate_files)
    
    print( '!!! Saving all results here: {}'.format(output_dir))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(intermediate_files_dir):
        os.makedirs(intermediate_files_dir)
    main_start = time.time() # overall program timer
    with open(os.path.join(intermediate_files_dir, "timer.txt"),"w") as timerfile:
        timerfile.write(str(main_start))
    return output_dir, intermediate_files_dir

def use_json_as_magi_input(jsonfilename, magi_parameters):
    """
    Set parameters from a MAGI job file to the current MAGI run.
    
    Inputs
    ------
    jsonfilename:       path to json file
    magi_parameters:    dictionary with magi parameters

    Outputs
    ------
    magi_parameters:    dictionary with updated parameters from the json file. No parameter checking is performed here.
    """
    with open(jsonfilename, 'r') as jsonfile:
        json_input = json.load(jsonfile)
    for key, value in json_input.items():
        if key == "blast_cutoff":
            magi_parameters["blast_filter"] = percentage_values_to_decimal(value)
        elif key == "reciprocal_cutoff":
            magi_parameters["reciprocal_closeness"] = percentage_values_to_decimal(value)
        elif key == "metabolite_file":
            magi_parameters["compounds"] = is_existing_file(value)
        elif key == "fasta_file":
            magi_parameters["fasta"] = is_existing_file(value)
        elif key == "score_weight_compound":
            magi_parameters["final_weights"][0] = positive_number(value)
        elif key == "score_weight_reciprocal":
            magi_parameters["final_weights"][1] = positive_number(value)
        elif key == "score_weight_homology":
            magi_parameters["final_weights"][2] = positive_number(value)
        elif key == "score_weight_rxnconnect":
            magi_parameters["final_weights"][3] = positive_number(value)
        elif key == "network_level":
            magi_parameters["level"] = value
        elif key == "output_directory":
            magi_parameters["output"] = value
        elif key == "ppm":
            magi_parameters["ppm_cutoff"] = value
        else:
            magi_parameters[key] = value
    return magi_parameters

def print_version_info():
    """
    Print versions of modules that may be troublesome for magi.
    """
    print('!!! Python version:'+ sys.version)
    print('!!! numpy version: '+ np.__version__)
    print('!!! pandas version:'+ pd.__version__)
    #print('!!! pickle version:'+ pickle.__version__)
    print('#'*80)

def print_parameters(args):
    """
    Print parameters used in the MAGI run.
    """    
    print('~~~~~PARAMETERS~~~~~~')

    # print your paths in stdout for logging
    print('@@@ FASTA file input: %s' %(args.fasta))
    print('@@@ Compound input: %s' %(args.compounds))
    if args.annotations is not None:
        print( '@@@ annotations input: %s' %(args.annotations))
    # print parameters
    if args.fasta is not None: 
        print( '@@@ BLAST filter: %s' % (args.blast_filter))
    if args.compounds is not None:
        print( '@@@ Using precomputed compound results: %s' % (not args.use_precomputed_reactions))
        print( '@@@ Similarity cutoff %s' % (args.similarity_cutoff))
        print( '@@@ Minimum diameter: %s' % (args.diameter))
        print( '@@@ Fingerprint radius: %s' % (args.fingerprint))
    print( '@@@ MAGI score weights: %s' % (args.final_weights))
    print( '@@@ Using %s CPUs' % (args.cpu_count))
    
    if args.gene_to_reaction is not None:
        print( '@@@ gene_to_reaction input: %s' %(args.gene_to_reaction))
    
    if args.compound_to_reaction is not None:
        print( '@@@ compound_to_reaction input: %s' %(args.compound_to_reaction))
    
    if args.reaction_to_gene is not None:
        print( '@@@ reaction_to_gene input: %s' %(args.reaction_to_gene))
    if args.mute:
        print( '!!! Warnings are muted')
        warnings.filterwarnings('ignore')
        
def general_magi_preparation():
    """
    This function prepares for a MAGI run. It:
        - parses arguments
        - makes an output file directory if the argument --not_first_script is not used.
        - prints versions of possibly troublesome modules if the argument --not_first_script is not used.
        - prints input parameters if the argument --not_first_script is not used.

    Outputs
    -------
    magi_parameters: A dictionary with parameters for the MAGI run.
    """
    args = parse_arguments()
    magi_parameters = vars(args)
    # Read parameters from json file if needed
    if magi_parameters["jsonfile"] is not None:
        magi_parameters = use_json_as_magi_input(magi_parameters["jsonfile"], magi_parameters)
    if magi_parameters["not_first_script"]:
        if magi_parameters["output"] is None:
            raise RuntimeError("Enter output file directory")
        else:
            output_dir_to_use = os.path.abspath(magi_parameters["output"])
            # Read used parameters from previous run
            magi_parameters = use_json_as_magi_input(os.path.join(output_dir_to_use, "used_parameters.json"), magi_parameters)
            # Overwrite output directory
            magi_parameters["output_dir"] = output_dir_to_use
            magi_parameters["intermediate_files_dir"] = os.path.join(output_dir_to_use, magi_parameters["intermediate_files"])
    else:
        print_version_info()
        print_parameters(args)
        output_dir, intermediate_files_dir = make_output_dirs(output_dir=magi_parameters["output"], 
                                                          fasta_file=magi_parameters["fasta"], 
                                                          compounds_file=magi_parameters["compounds"], 
                                                          intermediate_files=magi_parameters["intermediate_files"])
        magi_parameters["output_dir"] = output_dir
        magi_parameters["intermediate_files_dir"] = intermediate_files_dir
    
    # Check if only metabolite or fasta file is given
    if magi_parameters["fasta"] is None and magi_parameters["compounds"] is not None:
        magi_parameters["compound_to_reaction_only"] = True
    if magi_parameters["compounds"] is None and magi_parameters["fasta"] is not None:
        magi_parameters["gene_to_reaction_only"] = True
    # Write used parameters to file
    with open(os.path.join(magi_parameters["output_dir"], "used_parameters.json"), "w") as jsonfile:
        json.dump(magi_parameters, jsonfile)
    return magi_parameters
    
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
        file_ext = os.path.splitext(fname)[1][1:]
    else:
        file_ext = filetype
    try:
        if file_ext in ['pkl', 'pickle']:
            df = pd.read_pickle(fname)
        elif file_ext in ['xls', 'xlsx']:
            df = pd.read_excel(fname)
        elif file_ext in ['csv']:
            df = pd.read_csv(fname, engine='python', quoting=csv.QUOTE_NONE)
        elif file_ext in ['txt', 'tab', 'tsv']:
            df = pd.read_csv(fname, sep='\t', engine='python', quoting=csv.QUOTE_NONE)
        elif file_ext in ['h5', 'hdf5']:
            if key is None:
                raise IOError('"key" argument must be used when loading\
                    an HDF5 file')
            else:
                df = pd.read_hdf(fname, key)
        else:
            raise IOError('could not infer what type of file %s is... please \
                make it a csv, tab, or pickle file'%fname)
    except pd.errors.EmptyDataError:
        # Exit when file is empty
        msg = "The file {} seems to be empty. Specify no compounds file if you do not want to run the compound to reaction search.".format(fname)
        sys.exit(msg)

    # remove rows that are empty
    df = df[~pd.isnull(df).all(axis=1)]
    return df

def load_compound_results(compounds_file): 
    """ 
    load compound results.

    Inputs
    ------
    compounds_file: path to the file with compounds.
    """
    print( '\n!!! LOADING COMPOUNDS')
    compounds = load_dataframe(compounds_file)

    # Check if column with original compounds is present.
    if 'original_compound' not in compounds.columns:
        msg = 'Could not find "original_compound" as a column, please rename the column corresponding to SMILES for the compounds.'
        sys.exit(msg)
    # TODO: Check if column contains SMILES?

    # remove any missing compounds
    compounds = compounds[~pd.isnull(compounds['original_compound'])]
    compounds.fillna('', inplace=True)
    # remove duplicates
    compounds = compounds[compounds['original_compound'].isin(compounds['original_compound'].unique())]

    used_compounds_count = compounds['original_compound'].nunique()
    print( '!@# {} total input compounds to search\n'.format(used_compounds_count))

    if 'compound_score' not in compounds.columns:
        print( 'WARNING: "compound_score" not found as a column; assuming that there is no score for compounds, and setting the compound scores to 1.0')
        compounds['compound_score'] = 1.0
    else:
        compounds['compound_score'] = compounds['compound_score'].apply(float)
    return compounds

def get_settings():
    """
    Function to get local settings.

    Outputs
    -------
    my_settings: a module object with local settings.
    """
    my_settings = getattr(
    __import__(
        'local_settings',
        fromlist=[settings_loc.SETTINGS_FILE]), settings_loc.SETTINGS_FILE)
    return my_settings

def get_intermediate_file_path(output_dir, variable_of_interest):
    """
    This function will return the path that matches a path name in the used_parameters.json file.
    This is used to find paths to intermediate files that were made in previous magi runs.

    Inputs
    ------
    output_dir: The path to the output directory with a used_parameters.json file.
    variable_of_interest: The name of the variable for which the path needs to be found.

    Outputs
    -------
    variable_path:  The path stored in the used_parameters.json file.
    """
    with open(os.path.join(output_dir, "used_parameters.json"), "r") as jsonfile:
        used_parameters = json.load(jsonfile)
    # Try to find variable of interest
    try:
        variable_path = used_parameters[variable_of_interest]
    except:
        raise RuntimeError("Could not find object {} in used_parameters.json".format(variable_of_interest))
    # Return path if the file exists.
    if os.path.exists(variable_path) or is_database(variable_path):
        return variable_path
    else:
        raise OSError("File not found: {}".format(variable_path))

def write_intermediate_file_path(output_dir, variable_of_interest, variable_path_of_interest):
    """
    This function reads a file called used_parameters.json with paths stored and adds the new variable and its file path.

    Inputs
    ------
    output_dir: The path to the output files directory
    variable_of_interest: The name of the variable for which the path needs to be stored.
    path_of_interest: The file path that contains the variable_of_interest file.
    """
    # Read the json dictionary with file path object names and paths if it exists
    with open(os.path.join(output_dir, "used_parameters.json"), "r") as jsonfile:
        used_parameters = json.load(jsonfile)

    # Write the new path. This overwrites the old path.
    used_parameters[variable_of_interest] = variable_path_of_interest
    with open(os.path.join(output_dir, "used_parameters.json"), "w") as jsonfile:
         json.dump(used_parameters, jsonfile)

## TODO: Check if this is necessary
def load_mrs_reaction():
    """
    This function loads the mrs reaction database. 

    Outputs
    ------
    mrs_reaction: a pandas dataframe with information on reactions. Each row represents a reaction.
    """
    my_settings = get_settings()
    mrs_reaction_path = my_settings.mrs_reaction_path
    print( '!!! MRS-Reaction: {}'.format(mrs_reaction_path))
    mrs_reaction = load_dataframe(mrs_reaction_path)
    return mrs_reaction


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
