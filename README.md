# clumping

The `clumping` package reutilizes the TCR clumping and lit db matching tools from `CoNGA` and adapts them from use on paired `TIRTLseq` outputs.

`--lit_match` matches to TCRs of known antigen specifities.

`--mc_match` matches to metaCoNGA convergent clumps annotates clonotypes by the associated metadata.

The package requires `pandas`, `numpy`,`scipy`, and `scikit-learn`. Actually, I don't know if the last two are required to run but need to clean this up more.

# Installation and running

### Clone the repository

`git clone https://github.com/sschattgen/clumping.git`

### Install the package into your environment. We recommend using a conda environment.

`cd <path>/clumping`
###
`pip install -e .`

### Compile the C++ TCRdist and neighbor tools

`cd <path>/clumping/tcrdist_cpp && make`


### The tools could be run interactively or using the script `run_pipeline.py`:

`python scripts/run_pipeline.py --paired_data_tsv_file <path-to-tsv> --organism human --outfile_prefix <path-and-prefix-for-outputs>` and setting a flag for `--all`, `--clumping`, `--lit_match`, or `--mc_match`.

### Input format

The pipeline currently takes TIRTLseq output directly, but any paired TCR sequences in tabular format could be used by setting the column names to: the following:

`va` - alpha variable gene
`ja` - alpha joining gene
`vb` - beta variable gene
`jb` - beta joining gene
`cdr3a` - alpha CDR3 amino acid sequence
`cdr3b` - beta CDR3 amino acid sequence

For the test run:

`cd <path>/clumping`

`python scripts/run_pipeline.py --paired_data_tsv_file test/test.tsv --organism human --outfile_prefix test/tmp_TIRTL --all`

Currently, three files are output: `_metaCoNGA_matches.tsv` containing the metaCoNGA matches, `_lit_matches.tsv` containing the lit matches, `_clumping.tsv` containing the clumping groups, an `_clean_clones.tsv` as the cleaned version of the input used for the clumping and matching.

### Organism

Literature and metaCoNGA matching can only be performed using human abTCR data. Clumping will work for ab, gd, and ig data across all supported species.

### Contact 

Please email `stefan.schattgen` at `stjude` dot `org` for help. 
