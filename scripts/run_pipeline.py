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
    --all (run all available modes)

Examples:

    python {sys.argv[0]} --paired_data_tsv_file test.tsv --organism human --outfile_prefix tmp_TIRTL --all

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

parser.add_argument('--all',help='run clumping, lit match, and metaCoNGA match',
    action='store_true'
    )

parser.add_argument('--clumping',help='run TCR clumping only',
    action='store_true'
    )

parser.add_argument('--lit_match',help='run lit match only',
    action='store_true'
    )

parser.add_argument('--mc_match',help='run matching to metaCoNGA clumps only',
    action='store_true'
    )

args = parser.parse_args()

import pandas as pd
import clumping
import os

if not (args.all or args.clumping or args.lit_match or args.mc_match):
    raise ValueError("Must specify at least one of: all, clumping, lit_match, or mc_match")

if args.all:
    print("Processing clumping, lit matching, and metaCoNGA matching...")
    # Perform actions for all data
elif args.clumping:
    print("Performing clumping...")
    # Perform clumping-specific actions
elif args.lit_match:
    print("Performing lit matching...")
    # Perform literal matching-specific actions
elif args.mc_match:
    print("Performing metaCoNGA matching...")

if args.organism != 'human':
    if args.all:
        print('Literature and metaCoNGA matching only available for human. Only clumping will be run.')
        args.all = False
    if (args.lit_match or args.mc_match):
        print('Literature and metaCoNGA matching only available for human. Exiting.')
        sys.exit()

# read in the data and preprocess it
df = pd.read_table(args.paired_data_tsv_file)
df = clumping.preprocess.cleanup_paired_chains(df, args.organism)
tcrs = clumping.preprocess.tcrs_from_dataframe_helper(df, add_j_and_nucseq=True)
path_clean_clones = os.path.join(args.outfile_prefix + '_clean_clones.tsv')
df.to_csv(path_clean_clones, index = False, sep = '\t')

# clumping
if args.all or args.clumping:
    path_clumping = os.path.join(args.outfile_prefix + '_clumping.tsv')
    clump_results = clumping.tcr_clumping.find_tcr_clumping(tcrs, args.organism)
    clump_results.to_csv(path_clumping, index = False, sep = '\t')

# db lit matches
if args.all or args.lit_match:
    path_lit_matches = os.path.join(args.outfile_prefix + '_lit_matches.tsv')
    lit_matches = clumping.tcr_clumping.match_tcrs_to_db_tcrs(df, args.organism)
    lit_matches.to_csv(path_lit_matches, index = False, sep = '\t')

# metaconga
if args.all or args.mc_match:
    path_mc_matches = os.path.join(args.outfile_prefix + '_metaCoNGA_matches.tsv')
    mc_matches = clumping.tcr_clumping.match_tcrs_to_metaconga_clumps(df, args.organism)
    mc_matches.to_csv(path_mc_matches, index = False, sep = '\t')


