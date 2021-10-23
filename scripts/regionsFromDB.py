import pyBigWig, re, gffutils, os, sys, collections, pandas, glob, random
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats

"""
This script includes functions meant to take a gffutils database file
or object, and given some parameters produce lists of tuples of genomic intervals.

Specifically, what is expected is to call a load_cov function:

def load_cov(bw_fnames, db, featuretype='transcript', subset=False, target_size=10, params={}):
    sa_txpt = bigwigToNormArray.signalToNormArray(bw_fnames, bw_uses_chr_in_chrom_names=True)
    feat_select = regionsFromDB.get_feat_select(db, featuretype, params, subset=subset)
    n_feat = len(feat_select)
    cov_txpt = sa_txpt.read_coverages_in_regions(
        feat_select, n_features=n_feat, target_size=target_size)    
    return cov_txpt
    
load_cov is passed a list of bigwig filenames with an expected format, 
the passed ones being the positive strand files, and a database object,
and the feature type to subset from the database object, transcript, exon 
or intron or three prime UTR or five prime UTR etc., and a dictionary of 
parameters. 

This function might better be copied and pasted into
where it it is called more directly but it is included here for storage.

It creates a signalToNormArray object, then uses the functions in the script to
define a list of regions, then loads in the length normalized Signal from the bigwig
files and returns a dictionary object of bigwigs filename to two-dimensional
numpy array.

params = {
    # Either all ('all'), pass a set ('set'), or take the top ('top').
    'which_RNAs': 'all',
    
    # Either no norm (just sum RPMs, 'none') or each RNA to total ('sum').
    'how_norm': 'sum',
    
    # If which_RNAs='set', need this param (ignored otherwise):
    'gene_names': set(),
    
    # If which_RNAs='top', need this param (ignored otherwise):
    'top': 100,
}

# Tidy version to copy/paste:
#params = {'which_RNAs': 'all', 'how_norm': 'sum', 'gene_names': set(), 'top': 100}

The scripts exonTools, regionDef, and signalMatrix were created specifically to
test signals around differential exons. These scripts regionsFromDB and bigwigToNormArray
reflect the more general case of normalizing coverage across genomic elements.
"""

def an_iv_list(feat, strand=False):
    """Given A feature object from gffutils return an interval list."""
    if not strand:
        return [feat.chrom, feat.start, feat.end]
    return [feat.chrom, feat.start, feat.end, feat.strand]


def iv_iterable(iterable_of_feats, strand=False):
    """Generator for intervals from gffutils iterator of features."""
    for feat in iterable_of_feats:
        yield an_iv_list(feat, strand=strand)


def get_feat_select(db, featuretype, params, subset=False):
    """Given a gffutils database, featuretype and params dictionary,
    return a list of intervals."""
    
    feat_select = db.features_of_type(featuretype)

    if 'which_RNAs' not in params or params['which_RNAs']=='all':
        feat_select = list(iv_iterable(feat_select, strand=True)) 
        
    elif params['which_RNAs']=='set':
        feat_select = [x for x in feat_select if x.attributes['gene_name'][0] in params['gene_names']]
        feat_select = list(iv_iterable(feat_select, strand=True)) 

    elif params['which_RNAs']=='protein_coding':
        feat_select = [x for x in feat_select if x.attributes['transcript_type'][0]=='protein_coding']
        feat_select = list(iv_iterable(feat_select, strand=True)) 
        
    n_feat = len(feat_select)
    if subset and subset < n_feat:
        feat_select = random.sample(feat_select, k=subset)

    print(f"{len(feat_select)} features to load.")
    
    return feat_select




