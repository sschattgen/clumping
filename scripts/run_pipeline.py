######################################################################################88
import argparse
from os.path import exists
import sys
import os

parser = argparse.ArgumentParser(
    description='Run the CoNGA TCR clumping and lit match on TIRTLseq output.',
    epilog = f'''
==============================================================================

see https://github.com/sschattgen/clumping for more information.

Requires that a tsv formatted paired clones file output from Step 2 of TIRTL processing

Minimal command line arguments would be:

    --paired_data_tsv_file
    --organism
    --outfile_prefix (string that will be prepended to all outputs)

Examples:

    python {sys.argv[0]} --paired_data_tsv_file test.tsv --organism human --outfile_prefix tmp_TIRTL

    ''',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

parser.add_argument('--paired_data_tsv_file',
                    help='Input file with the paired TIRTLseq data')

parser.add_argument('--organism',
                    choices=['mouse', 'human', 'mouse_gd', 'human_gd',
                             'human_ig', 'rhesus', 'rhesus_gd'])

parser.add_argument('--outfile_prefix',
                    help='string that will be prepended to all output files'
                    ' and images')

args = parser.parse_args()

import pandas as pd
import clumping
import os

# read in the data and preprocess it
df = pd.read_table(args.paired_data_tsv_file)
df = clumping.preprocess.cleanup_paired_chains(df, args.organism)
tcrs = clumping.preprocess.tcrs_from_dataframe_helper(df, add_j_and_nucseq=True)
path_clean_clones = os.path.join(args.outfile_prefix + '_clean_clones.tsv')
df.to_csv(path_clean_clones, index = False, sep = '\t')

# clumping
path_clumping = os.path.join(args.outfile_prefix + '_clumping.tsv')
clump_results = clumping.tcr_clumping.find_tcr_clumping(tcrs, args.organism)
clump_results.to_csv(path_clumping, index = False, sep = '\t')

# db lit matches
path_lit_matches = os.path.join(args.outfile_prefix + '_lit_matches.tsv')
lit_matches = clumping.tcr_clumping.match_tcrs_to_db_tcrs(df, args.organism)
lit_matches.to_csv(path_lit_matches, index = False, sep = '\t')
