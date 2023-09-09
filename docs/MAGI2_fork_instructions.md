# Using this fork of MAGI2
---
## Instructions:
1. Clone repository and set up environment
```bash
git clone git@github.com:TimothyStephens/magi.git
cd magi
conda env create -f magi_2_env.yml
source activate magi_2
python setup_magi2.py
```

2. Install NCBI BLAST
Two [NCBI BLAST](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/)
binaries are required to run MAGI.
These tools should have already been installed by conda, so you can use the following commands to copy the binaries into the magi bin directory.

```bash
P=$(which blastp); cp $P workflow/blastbin/
P=$(which makeblastdb); cp $P workflow/blastbin/
```

3. Run MAGI2
The MAGI2 scripts are in `workflow_2/`. Within this directory their is a script `magi2.sh` which is the main MAGI2 working script.
To get a detailed description of all options for this script, simply run `./magi2.sh -h`.

This script requires at least a fasta file (`-f`), and a feature file with either m/z (`-m`) values or SMILES (`-s`).

Example run using m/z values. The m/z file should be single column and not have a column header.
```bash
./magi2.sh -f tests/test_pep.fa -m tests/test_mz.txt
```

Example run using SMILES values, which need to be in a column with a header called 'original_compound'.
```bash
./magi2.sh -f tests/test_pep.fa -m tests/test_SMILES.txt
```

The final (filtered) MAGI2 results will be in.
```bash
output_magi2.filtered_magi_results.csv
output_magi2.filtered_magi_gene_results.csv
output_magi2.filtered_magi_compound_results.csv
```
This selects for just results with compound_score = 1, e_score_r2g (reaction-to-gene) > 5, e_score_g2r (gene-to-reaction) > 5, and reciprocal_score = 2.

If you want to use different cutoffs, you can rerun filtering with the `filter_results.py` script.

