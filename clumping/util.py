import numpy as np
import sys
from os import system
import os.path
from pathlib import Path
import os
from scipy.sparse import issparse
from scipy.stats import mannwhitneyu
from collections import Counter, OrderedDict
import subprocess

# try not to have any conga imports here
#

# convenience paths
path_to_conga = Path(__file__).parent
assert os.path.isdir( path_to_conga )

path_to_data = Path.joinpath( path_to_conga, 'data')
assert os.path.isdir( path_to_data )

path_to_tcrdist_cpp = Path.joinpath( path_to_conga.parents[0] ,'tcrdist_cpp')
path_to_tcrdist_cpp_bin = Path.joinpath( path_to_tcrdist_cpp ,'bin')
path_to_tcrdist_cpp_db = Path.joinpath( path_to_tcrdist_cpp ,'db')
assert os.path.isdir( path_to_tcrdist_cpp_bin ) and os.path.isdir( path_to_tcrdist_cpp_db )


def tcrdist_cpp_available():
    if os.name == 'posix':
        return os.path.exists(Path.joinpath( path_to_tcrdist_cpp_bin ,'find_neighbors'))
    else:
        return os.path.exists(Path.joinpath( path_to_tcrdist_cpp_bin ,'find_neighbors.exe'))

# this is the (OPTIONAL) obs key used to store a subject-specific identifier
# right now this is only used to prevent condensing of clonotypes that span subjects,
# which is pretty unlikely but could happen with certain populations (e.g., MAIT cells)
#
SUBJECT_ID_OBS_KEY = 'subject_id'

# not a big deal, but if we have protein ie antibody data we use these to mask it out
GENE_EXPRESSION_FEATURE_TYPE = 'Gene Expression'
ANTIBODY_CAPTURE_FEATURE_TYPE = 'Antibody Capture'

EXPECTED_FEATURE_TYPES = [GENE_EXPRESSION_FEATURE_TYPE, ANTIBODY_CAPTURE_FEATURE_TYPE]

FUNNY_MOUSE_TRBV_GENE = '5830405F06Rik' # actually seems to be a tcr v gene transcript or correlate with one
FUNNY_HUMAN_IG_GENES = [
    'AC233755.1', 'AC233755.2', # seem to be associated with one or more IGHV genes
    'CH17-224D4.2', # chr14 bac, suspiciously high correlation with tcr features??
    'IGLL5', # correlated with IGLJ1
]
FUNNY_HUMAN_TR_GENES = [
    'TRD-AS1', # overlaps TRA
    'ENSG00000251002', # overlaps TRA (alt name for TRD-AS1)
    'ENSG00000288882', # overlaps TRB
    'ENSG00000289938', # overlaps TRB
]

def run_command( cmd, verbose=False ):

    if verbose:
        print('util.run_command: cmd=', cmd)

    if os.name == 'posix':
        system(cmd)
    else:
        cmd_run = 'cmd /c' + cmd
        subprocess.check_call(list(cmd_run.split(' ')))


# different types of repertoire data we might have
TCR_AB_VDJ_TYPE = 'TCR_AB_VDJ_TYPE'
TCR_GD_VDJ_TYPE = 'TCR_GD_VDJ_TYPE'
IG_VDJ_TYPE = 'IG_VDJ_TYPE'

organism2vdj_type = {
    'human':TCR_AB_VDJ_TYPE,
    'mouse':TCR_AB_VDJ_TYPE,
    'human_gd':TCR_GD_VDJ_TYPE,
    'mouse_gd':TCR_GD_VDJ_TYPE,
    'human_ig':IG_VDJ_TYPE,
    'mouse_ig':IG_VDJ_TYPE,
    'rhesus':TCR_AB_VDJ_TYPE,
    'rhesus_gd':TCR_GD_VDJ_TYPE,
}

def is_vdj_gene( gene_upper, organism, include_constant_regions=False ):
    # for filtering out TR or IG gene names from GEX prior to processing
    # or for skipping such genes in the graph_vs_features analysis
    vdj_type = organism2vdj_type[organism]

    gene = gene_upper.lower()
    if vdj_type == TCR_AB_VDJ_TYPE:
        return ( gene.startswith('trav') or gene.startswith('trbv') or
                 gene.startswith('traj') or gene.startswith('trbj') or
                 gene.startswith('trbd') or gene_upper == FUNNY_MOUSE_TRBV_GENE or
                 gene_upper in FUNNY_HUMAN_TR_GENES or
                 ( include_constant_regions and (gene.startswith('trac') or gene.startswith('trbc'))))

    elif vdj_type == TCR_GD_VDJ_TYPE:
        return ( gene.startswith('trav') or gene.startswith('trdv') or
                 gene.startswith('traj') or gene.startswith('trdj') or
                 gene.startswith('trgv') or gene.startswith('trgj') or
                 gene.startswith('tcrg-') or gene.startswith('trdd') or
                 ( include_constant_regions and (gene.startswith('trdc') or gene.startswith('trgc'))))
    elif vdj_type == IG_VDJ_TYPE:
        return ( gene.startswith('ighv') or gene.startswith('iglv') or gene.startswith('igkv') or
                 gene.startswith('ighj') or gene.startswith('iglj') or gene.startswith('igkj') or
                 gene.startswith('ighd') or gene_upper in FUNNY_HUMAN_IG_GENES or
                 ( include_constant_regions and
                   ( gene.startswith('ighc') or gene.startswith('iglc') or gene.startswith('igkc'))))
    else:
        print('unrecognized vdj_type:', vdj_type)
        exit()
    return None

def make_clones_file( tcrs, outfilename, subject = 'UNK', epitope = 'UNK_E' ):
    ''' This may not have all the standard fields
    Right now just adding the fields we need in order for make_tcr_logo.py to work...
    '''
    gene_fields = ['{}{}_{}'.format(x,y,z) for x in 'vj' for y in 'ab' for z in ['gene', 'genes']]
    outfields = 'clone_id subject epitope cdr3a cdr3a_nucseq cdr3b cdr3b_nucseq'.split() + gene_fields

    out = open(outfilename, 'w')
    out.write('\t'.join(outfields)+'\n')
    for ii,(atcr, btcr) in enumerate(tcrs):
        outl = { 'clone_id': 'clone_{}'.format(ii+1),
                 'subject': subject,
                 'epitope': epitope,
                 'va_gene': atcr[0],
                 'va_genes': atcr[0],
                 'ja_gene': atcr[1],
                 'ja_genes': atcr[1],
                 'cdr3a': atcr[2],
                 'cdr3a_nucseq': atcr[3],
                 'vb_gene': btcr[0],
                 'vb_genes': btcr[0],
                 'jb_gene': btcr[1],
                 'jb_genes': btcr[1],
                 'cdr3b': btcr[2],
                 'cdr3b_nucseq': btcr[3],
                 }
        out.write('\t'.join( outl[x] for x in outfields)+'\n')
    out.close()


