################################################################################
##
## this file is a place to put code that is experimental
##  or otherwise under development / not ready for prime time
##
##

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
from scipy.stats import hypergeom, mannwhitneyu, linregress, norm, ttest_ind, poisson
import scipy.sparse as sps
from scipy.sparse import issparse, csr_matrix
from statsmodels.stats.multitest import multipletests
from collections import Counter, OrderedDict
from . import preprocess
from . import util
from . import tcr_clumping
from .tcrdist.all_genes import all_genes
import sys
import pandas as pd
from sys import exit
import time #debugging
import random
from os.path import exists
import os
from pathlib import Path

def _get_pvalue_from_rvalue(r, n):
    ''' hacky helper stolen from scipy.stats.linregress
    alternative='two-sided'
    '''
    from scipy import stats
    TINY = 1.0e-20
    df = n - 2  # Number of degrees of freedom
    # n-2 degrees of freedom because 2 has been used up
    # to estimate the mean and standard deviation
    t = r * np.sqrt(df / ((1.0 - r + TINY)*(1.0 + r + TINY)))
    #t, prob = _ttest_finish(df, t, alternative)
    prob = 2*stats.t.cdf(-abs(t), df)
    return prob


def find_tcr_clumping_single_chain(
        tcrs,
        organism,
        tmpfile_prefix = 'tmp',
        radii = [6, 12, 24, 48],
        num_random_samples_multiplier = 100, # ie, 100 * num_clones
        pvalue_threshold = 1.0,
        verbose=True,
        preserve_vj_pairings = False,
        bg_tcrs = None, # usually better to leave this None
):
    ''' Returns a pandas dataframe with the following columns:
    - chain (A or B)
    - clone_index
    - nbr_radius
    - pvalue_adj
    - num_nbrs
    - expected_num_nbrs
    - raw_count
    - va, ja, cdr3a, vb, jb, cdr3b (ie, the 6 tcr cols for clone_index clone)
    - clumping_group: clonotypes within each other's significant nbr_radii are linked
    - clump_type: string, either 'global' or 'intra_gex_cluster' (latter only if also_find_clumps_within_gex_clusters=T)

    '''

    outprefix = f'{tmpfile_prefix}_{random.random()}_tcr_clumping'

    num_clones = len(tcrs)
    num_random_samples = num_random_samples_multiplier * num_clones

    if bg_tcrs is None:
        bg_tcrs = tcrs

    agroups, bgroups = preprocess.setup_tcr_groups_for_tcrs(tcrs)

    tmpfiles = [] # for cleanup

    ## compute background neighbor counts at the specified radii
    tcr_clumping.estimate_background_tcrdist_distributions(
        organism, tcrs, max(radii), num_random_samples,
        tmpfile_prefix=outprefix,
        preserve_vj_pairings=preserve_vj_pairings,
        save_unpaired_dists=True,
        tcrs_for_background_generation=bg_tcrs,
        #nocleanup=True,
    )

    afile = f'{outprefix}_dists.txt_A.txt'
    bfile = f'{outprefix}_dists.txt_B.txt'
    if not (exists(afile) and exists(bfile)):
        print('ERROR find_tcr_clumping_single_chain::',
              'estimate_background_tcrdist_distributions failed')
        return pd.DataFrame()
    tmpfiles.extend([afile, bfile])

    acounts = np.cumsum(np.loadtxt(afile, dtype=int), axis=1)[:,radii]
    bcounts = np.cumsum(np.loadtxt(bfile, dtype=int), axis=1)[:,radii]

    # newcol = np.full((len(tcrs),1), num_random_samples)
    # acounts = np.hstack([acounts[:,radii], newcol])
    # bcounts = np.hstack([bcounts[:,radii], newcol])


    big_dfl = [] # will hold results for both chains as dataframes

    tcrs_file = outprefix +'_tcrs.tsv'
    pd.DataFrame({
        'va'   : [x[0][0] for x in tcrs],
        'cdr3a': [x[0][2] for x in tcrs],
        'vb'   : [x[1][0] for x in tcrs],
        'cdr3b': [x[1][2] for x in tcrs],
    }).to_csv(tcrs_file, sep='\t', index=False)
    tmpfiles.append(tcrs_file)

    for chain, chain_groups, bg_counts in [['A', agroups, acounts],
                                           ['B', bgroups, bcounts]]:

        # find neighbors in fg tcrs up to max(radii) ############################

        if os.name == 'posix':
            exe = Path.joinpath(
                Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors_single_chain')
        else:
            exe = Path.joinpath(
                Path(util.path_to_tcrdist_cpp_bin) , 'find_neighbors_single_chain.exe')


        db_filename = Path.joinpath(
            Path(util.path_to_tcrdist_cpp_db), f'tcrdist_info_{organism}.txt')

        tcrdist_threshold = max(radii)

        cmd = (f'{exe} -f {tcrs_file} -t {tcrdist_threshold} -d {db_filename}'
               f' -c {chain} -o {outprefix}')

        util.run_command(cmd, verbose=True)

        nbr_indices_filename = f'{outprefix}_nbr{tcrdist_threshold}_indices.txt'
        nbr_distances_filename = f'{outprefix}_nbr{tcrdist_threshold}_distances.txt'
        tmpfiles.extend([nbr_indices_filename, nbr_distances_filename])

        if not exists(nbr_indices_filename) or not exists(nbr_distances_filename):
            print('find_neighbors failed:', exists(nbr_indices_filename),
                  exists(nbr_distances_filename))
            exit(1)

        all_nbrs = []
        all_distances = []
        for group, line1, line2 in zip(chain_groups,
                                       open(nbr_indices_filename,'r'),
                                       open(nbr_distances_filename,'r')):
            nbrs = [int(x) for x in line1.split()]
            dists = [int(x) for x in line2.split()]
            assert len(nbrs) == len(dists)
            mask = [chain_groups[x] != group for x in nbrs]
            all_nbrs.append([x for x,m in zip(nbrs,mask) if m])
            all_distances.append([x for x,m in zip(dists,mask) if m])
        assert len(all_nbrs) == num_clones

        # use poisson to find nbrhoods with more tcrs than expected;
        #  have to handle agroups/bgroups
        dfl = []

        is_clumped = np.full((num_clones,), False)

        all_raw_pvalues = np.full((num_clones, len(radii)), 1.0)

        for ii in range(num_clones):
            ii_bg_counts = bg_counts[ii]
            ii_dists = all_distances[ii]
            for irad, radius in enumerate(radii):
                num_nbrs = sum(x<=radius for x in ii_dists)
                if num_nbrs<1:
                    continue
                max_nbrs = np.sum(chain_groups != chain_groups[ii])
                pval = hypergeom.sf(
                    num_nbrs-1, # total fg nbr tcrs (-1)
                    max_nbrs+num_random_samples, # total fg+bg tcrs
                    num_nbrs+ii_bg_counts[irad], # total nbr tcrs
                    max_nbrs) # total fg tcrs
                mu = max_nbrs * ii_bg_counts[irad]/num_random_samples
                all_raw_pvalues[ii, irad] = pval
                pval *= len(radii) * num_clones # simple multiple test correction
                #hgeom_pval *= len(radii) * num_clones # simple multiple test correction
                if pval <= pvalue_threshold:
                    is_clumped[ii] = True
                    raw_count = ii_bg_counts[irad]
                    if verbose:
                        atcr_str = ' '.join(tcrs[ii][0][:3])
                        btcr_str = ' '.join(tcrs[ii][1][:3])
                        print(f'tcr_nbrs_global: {num_nbrs:2d} {mu:9.6f}',
                              f'radius: {radius:2d} pval: {pval:9.1e}',
                              f'{raw_count:9.1f} tcr: {atcr_str} {btcr_str}')

                    dfl.append( OrderedDict(
                        chain=chain,
                        clone_index=ii,
                        nbr_radius=radius,
                        pvalue_adj=pval,
                        #hgeom_pvalue_adj=hgeom_pval,
                        num_nbrs=num_nbrs,
                        expected_num_nbrs=mu,
                        raw_count=raw_count,
                        va   =tcrs[ii][0][0],
                        ja   =tcrs[ii][0][1],
                        cdr3a=tcrs[ii][0][2],
                        vb   =tcrs[ii][1][0],
                        jb   =tcrs[ii][1][1],
                        cdr3b=tcrs[ii][1][2],
                    ))

        results_df = pd.DataFrame(dfl)
        if results_df.shape[0] == 0:
            big_dfl.append(results_df)
            continue # to next chain


        # compute FDR values in addition to the simple adjusted pvalues
        _, fdr_values, _, _ = multipletests(
            all_raw_pvalues.reshape(-1), alpha=0.05, method='fdr_bh')
        fdr_values = fdr_values.reshape((num_clones, len(radii))).min(axis=1)
        # right now we don't assign fdr values for intra-gex cluster clumping
        results_df['clonotype_fdr_value'] = [
            fdr_values[x] for x in results_df.clone_index]

        # identify groups of related hits?
        all_clumped_nbrs = {}
        for l in results_df.itertuples():
            ii = l.clone_index
            radius = l.nbr_radius
            clumped_nbrs = set(x for x,y in zip(all_nbrs[ii], all_distances[ii])
                               if y<= radius and is_clumped[x])
            clumped_nbrs.add(ii)
            if ii in all_clumped_nbrs:
                all_clumped_nbrs[ii] = all_clumped_nbrs[ii] | clumped_nbrs
            else:
                all_clumped_nbrs[ii] = clumped_nbrs


        clumped_inds = sorted(all_clumped_nbrs.keys())
        assert len(clumped_inds) == np.sum(is_clumped)

        # make nbrs symmetric
        for ii in clumped_inds:
            for nbr in all_clumped_nbrs[ii]:
                assert nbr in all_clumped_nbrs
                all_clumped_nbrs[nbr].add(ii)

        all_smallest_nbr = {}
        for ii in clumped_inds:
            all_smallest_nbr[ii] = min(all_clumped_nbrs[ii])

        while True:
            updated = False
            for ii in clumped_inds:
                nbr = all_smallest_nbr[ii]
                new_nbr = min(nbr, np.min([all_smallest_nbr[x]
                                           for x in all_clumped_nbrs[ii]]))
                if nbr != new_nbr:
                    all_smallest_nbr[ii] = new_nbr
                    updated = True
            if not updated:
                break

        # define clusters, choose cluster centers
        clusters = np.array([0]*num_clones) # 0 if not clumped

        cluster_number=0
        cluster_sizes = Counter()
        for ii in clumped_inds:
            nbr = all_smallest_nbr[ii]
            if ii==nbr:
                cluster_number += 1
                members = [ x for x,y in all_smallest_nbr.items() if y==nbr]
                clusters[members] = cluster_number
                cluster_sizes[cluster_number] = len(members)

        for ii, nbrs in all_clumped_nbrs.items():
            for nbr in nbrs:
                # confirm single-linkage clusters
                assert clusters[ii] == clusters[nbr]

        assert not np.any(clusters[is_clumped]==0)
        assert np.all(clusters[~is_clumped]==0)

        # reorder the clumping groups by size
        remap = {x[0]:i+1 for i,x in enumerate(cluster_sizes.most_common())}
        remap[0] = 0

        results_df['clumping_group'] = [ remap[clusters[x.clone_index]]
                                         for x in results_df.itertuples()]

        results_df.sort_values('pvalue_adj', inplace=True)

        big_dfl.append(results_df)

    # cleanup the tmpfiles
    for tmpfile in tmpfiles:
        if exists(tmpfile):
            os.remove(tmpfile)

    results = pd.concat(big_dfl)
    return results


#def make_single_chain_clumping_logos(
#        results_df, # generated by find_tcr_clumping_single_chain
#        adata,
#        nbrs_gex,
#        nbrs_tcr,
#        outfile_prefix,
#        min_cluster_size_for_logos=3,
#        pvalue_threshold_for_logos=1.0, #pvalues are crude bonferroni corrected
#        max_color_pvalue = 1e-16, # in the UMAP figure; no brighter beyond there
#        **logo_plot_args,
#):
#    ''' Make tcr clumping logos for each chain individually, based
#    on previously calculated single-chain clumping results dataframe
#    '''
#    num_clones = adata.shape[0]#

#    for chain in 'AB':
#        fake_clusters_gex = np.zeros((num_clones,)).astype(int)
#        fake_clusters_tcr = np.zeros((num_clones,)).astype(int)
#        clumping_pvals = np.full( (num_clones,), num_clones).astype(float)#
#

#        for l in results_df.itertuples():
#            if l.chain == chain:
#                clumping_pvals[ l.clone_index] = min(l.pvalue_adj,
#                                                     clumping_pvals[l.clone_index])
#                fake_clusters_tcr[l.clone_index] = l.clumping_group#

#        pngfile = f'{outfile_prefix}_{chain}_chain_clumping_logos.png'#

#        if 'rank_genes_uns_tag' not in logo_plot_args:
#            logo_plot_args['rank_genes_uns_tag'] = f'rg_tcr_clumping_{chain}_biclusters'#

#        if 'conga_scores_name' not in logo_plot_args:
#            logo_plot_args['conga_scores_name'] = f'TCR{chain} clumping'#

#        if 'show_real_clusters_gex' not in logo_plot_args:
#            logo_plot_args['show_real_clusters_gex'] = True#

#        plotting.make_cluster_logo_plots_figure(
#            adata, clumping_pvals, pvalue_threshold_for_logos,
#            fake_clusters_gex, fake_clusters_tcr, nbrs_gex, nbrs_tcr,
#            min_cluster_size_for_logos, pngfile,
#            **logo_plot_args)

