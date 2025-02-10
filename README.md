# clumping

The `clumping` package reutilizes the TCR clumping and lit db matching tools from `CoNGA` and adapts them from use on paired `TIRTLseq` outputs.

The package requires `pandas`, `numpy`,`scipy`, and `scikit-learn`. Actually, I don't know if the last two are required to run but need to clean this up more.

# Installation and running

### Clone the repository

`git clone https://github.com/sschattgen/clumping.git`

### Install the package into your environment (I installed it in my conga conda env)

`cd <path>/clumping`
###
`pip install -e .`

### Compile the C++ TCRdist and neighbor tools

`cd <path>/clumping/tcrdist_cpp && make`


### The tools could be run interactively or using the script `run_pipeline.py`:

`python scripts/run_pipeline.py --paired_data_tsv_file <path-to-tsv> --organism human --outfile_prefix <path-and-prefix-for-outputs>` and setting a flag for `--all`, `--clumping`, or `--lit_match`.



For the test run:

`cd <path>/clumping`
###
`python scripts/run_pipeline.py --paired_data_tsv_file test/test.tsv --organism human --outfile_prefix test/tmp_TIRTL --all`

Currently three files are output: `_lit_matches.tsv` containing the lit matches, `_clumping.tsv` containing the clumping groups, an `_clean_clones.tsv` as the cleaned version of the input used for the clumping and matching.

