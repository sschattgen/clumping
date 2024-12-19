######################################################################################88
import random
import pandas as pd
from os.path import exists
from pathlib import Path
from collections import Counter
from sklearn.metrics import pairwise_distances
from sklearn.utils import sparsefuncs
from sklearn.decomposition import KernelPCA
import numpy as np
import scipy
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform, cdist
from scipy.sparse import issparse, csr_matrix
import sys
import os
from sys import exit
from . import util
from .tcrdist.tcr_distances import TcrDistCalculator
from .util import tcrdist_cpp_available
from .tcrdist.all_genes import all_genes


tcr_keys = 'va ja cdr3a cdr3a_nucseq vb jb cdr3b cdr3b_nucseq'.split()

MIN_CDR3_LEN = 6

def add_allele(df):
    for i in ['va','ja','vb','jb']:
        df[i] = df[i].apply(lambda x: x + "*01")
    return df

def cleanup_paired_chains(df, organism):

    try:
        df = df.rename(columns={ 'alpha_nuc':'cdr3a_nucseq','beta_nuc':'cdr3b_nucseq'})
    except:
        df = df.rename(columns={ 'alpha_nuc_seq':'cdr3a_nucseq','beta_nuc_seq':'cdr3b_nucseq'})

    df = df[(~df['cdr3a'].str.contains('_|\*')) & (~df['cdr3b'].str.contains('_|\*'))].copy()
    df = df.drop_duplicates(subset=['va','ja','vb','jb','cdr3a','cdr3b'])

    #sanity check that the len of aa seq goes into len of nuc seq 3 times
    df['cdr3a_nucseq_len'] = df['cdr3a_nucseq'].apply(lambda x: len(x))
    df['cdr3a_len'] = df['cdr3a'].apply(lambda x: len(x))
    df['len_match_alpha'] = df['cdr3a_nucseq_len'] == df['cdr3a_len']*3
    assert df['len_match_alpha'].all()
    
    df['cdr3b_nucseq_len'] = df['cdr3b_nucseq'].apply(lambda x: len(x))
    df['cdr3b_len'] = df['cdr3b'].apply(lambda x: len(x))
    df['len_match_beta'] = df['cdr3b_nucseq_len'] == df['cdr3b_len']*3
    assert df['len_match_beta'].all()

    df = add_allele(df)

    print(f'{df.shape[0]} clones prior to filtering for {organism} genes')

    df = df[df['va'].isin(all_genes[organism].keys())]
    df = df[df['ja'].isin(all_genes[organism].keys())]
    df = df[df['vb'].isin(all_genes[organism].keys())]
    df = df[df['jb'].isin(all_genes[organism].keys())]

    print(f'{df.shape[0]} clones after filtering for {organism} genes')

    df = df[(df['cdr3a_len'] >= MIN_CDR3_LEN) & (df['cdr3b_len'] >= MIN_CDR3_LEN)]

    print(f'{df.shape[0]} clones after filtering for CDR3 length =< {MIN_CDR3_LEN}')


    return df

def retrieve_tcrs_from_df(df):
    tcrs = []
    arrays = [ df[x] for x in tcr_keys ] # the keys set the columns, but doesn't matter how unpacked below
    for va,ja,cdr3a,cdr3a_nucseq,vb,jb,cdr3b,cdr3b_nucseq in zip(*arrays):
        tcrs.append(((va, ja, cdr3a, cdr3a_nucseq.lower()),
                        (vb, jb, cdr3b, cdr3b_nucseq.lower()) ) )
    return tcrs

def tcrs_from_dataframe_helper(df, add_j_and_nucseq=False):
    ''' Helper function that creates a tcrs list from a dataframe
    the dataframe should have columns va, cdr3a, vb, cdr3b

    if add_j_and_nucseq is True, then the dataframe should also have the
    columns ja, cdr3a_nucseq, jb, cdr3b_nucseq
    '''
    if add_j_and_nucseq:
        return [ ( (x.va, x.ja, x.cdr3a, x.cdr3a_nucseq.lower()),
                   (x.vb, x.jb, x.cdr3b, x.cdr3b_nucseq.lower()) )
                 for x in df.itertuples() ]
    else:
        return [ ( (x.va, None, x.cdr3a), (x.vb, None, x.cdr3b) )
                 for x in df.itertuples() ]


def setup_tcr_groups_for_tcrs( tcrs ):
    atcrs = sorted( set( x[0] for x in tcrs ) )
    btcrs = sorted( set( x[1] for x in tcrs ) )

    atcr2agroup = dict( (y,x) for x,y in enumerate(atcrs))
    btcr2bgroup = dict( (y,x) for x,y in enumerate(btcrs))

    agroups = np.array( [ atcr2agroup[x[0]] for x in tcrs] )
    bgroups = np.array( [ btcr2bgroup[x[1]] for x in tcrs] )

    return agroups, bgroups

