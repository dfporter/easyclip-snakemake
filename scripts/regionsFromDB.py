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
and the feature type to subset from the database object (transcript, exon 
or intron or three prime UTR or five prime UTR etc.), and a dictionary of 
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
    """Given A feature object from gffutils return an interval list.
    Subtracting 1 from the start because gffs are 1-based closed.
    This makes 0-based half-open.
    """
    if not strand:
        return [feat.chrom, feat.start-1, feat.end]
    return [feat.chrom, feat.start-1, feat.end, feat.strand]


def iv_iterable(iterable_of_feats, strand=False):
    """Generator for intervals from gffutils iterator of features."""
    for feat in iterable_of_feats:
        yield an_iv_list(feat, strand=strand)

def one_txpt_per_gene(db, feat_select, how='longest'):
    gene_to_txpt = collections.defaultdict(set)

    if how == 'random':
        for feat in feat_select:
            gene_to_txpt[feat.attributes._d['gene_name'][0]].add(feat.attributes._d['transcript_id'][0])
        for gene in gene_to_txpt:
            gene_to_txpt[gene] = random.choice(list(gene_to_txpt[gene]))
        
        # Now subset:
        feat_select = [
            feat for feat in feat_select \
            if gene_to_txpt[feat.attributes._d['gene_name'][0]]==feat.attributes._d['transcript_id'][0]]
        
    elif how == 'longest':
        txpt_to_len = {}
        for feat in feat_select:
            gene = feat.attributes._d['gene_name'][0]
            txpt = feat.attributes._d['transcript_id'][0]
            txpt_to_len[txpt] = feat.end - feat.start
            
            if gene in gene_to_txpt:
                if txpt_to_len[txpt] > txpt_to_len[gene_to_txpt[gene]]:
                    gene_to_txpt[gene] = txpt
            else:
                gene_to_txpt[gene] = txpt
            
        # Now subset:
        print(f"Subsetting to one txpt per gene, starting with {len(feat_select)} features...")
        feat_select = [
            feat for feat in feat_select \
            if gene_to_txpt[feat.attributes._d['gene_name'][0]]==feat.attributes._d['transcript_id'][0]]
        print(f"... Ended up with {len(feat_select)} features...")
        
    return feat_select

def as_mature_rna(db, feat_select, strand=True):
    subregions = []
    for feat in feat_select:
        exons = [an_iv_list(x, strand=strand) for x in \
                 db.children(feat.id, featuretype='exon', order_by='start')]
        if len(exons) > 0 and len(exons[0])>0:
            subregions.append(exons)
    return subregions

def organize_regions_from_gene_together(db, feat_select, featuretype, strand=True):
    subregions = []
    for feat in feat_select:
        exons = [an_iv_list(x, strand=strand) for x in \
                     db.children(feat.id, featuretype=featuretype, order_by='start')]
        if len(exons) > 0 and len(exons[0])>0:
            subregions.append(exons)
    return subregions

def filtering(db, feat_select, params, _one_txpt_per_gene=False):
    
    if 'which_RNAs' not in params or params['which_RNAs']=='all':
        pass
    elif params['which_RNAs']=='set':
        feat_select = [x for x in feat_select if x.attributes['gene_name'][0] in params['gene_names']]
    elif params['which_RNAs']=='protein_coding':
        feat_select = [x for x in feat_select if x.attributes['transcript_type'][0]=='protein_coding']
    
    if 'one_txpt_per_gene' in params and params['one_txpt_per_gene']:
        
        feat_select = one_txpt_per_gene(db, feat_select, how=params['one_txpt_per_gene'])
    
    return feat_select

def get_feat_select(
    db, featuretype, params, subset=False, _one_txpt_per_gene=False, mature_rna=False):
    """Given a gffutils database, featuretype and params dictionary,
    return a list of intervals."""
    
    mature_rna = params['mature_rna'] if ('mature_rna' in params) else False
            
    if featuretype == 'transcript':
        feat_select = db.features_of_type(featuretype)
        feat_select = filtering(db, feat_select, params)
        
        if mature_rna:
            feat_select = as_mature_rna(db, feat_select, strand=True)
        else:
            feat_select = list(iv_iterable(feat_select, strand=True)) 

    else:
        feat_select = db.features_of_type('transcript')
        feat_select = filtering(db, feat_select, params)
        feat_select = organize_regions_from_gene_together(db, feat_select, featuretype, strand=True)

    n_feat = len(feat_select)
    if subset and subset < n_feat:
        feat_select = random.sample(feat_select, k=subset)

    print(f"{len(feat_select)} features to load.")
    if len(feat_select):
        print(f"The first is {feat_select[0]}")
    
    return feat_select




