import os, sys, re, glob, pandas, importlib, shutil, math, pysam, pickle, collections, dill, gffutils
import numpy as np
import matplotlib.pyplot as plt
import scripts
import scripts.exp
import scripts.init_with_config
import scripts.rnaDataFileMaker
import scripts.make_repeats_chrom
import scripts.split_bam
import scripts.random_sequence

# Load config file.
configfile: "config.yaml"


#########################################################
# Set global path values and initialize exp object.
#########################################################

"""
# If loading the config outside snakemake is ever desired, use this code block:
configfile = "config.yaml"
config = {}
with open(configfile) as f:
    for li in f:
        li = li.split('#')[0].rstrip('\n')
        if ':' not in li:
            continue
        s = [x.strip(' ') for x in li.split(':')]
        config[s[0]] = s[1]
"""

config = scripts.init_with_config.init_with_config(config)

# Create exp object for organizing some tasks.
ex = scripts.exp.exp(name=config['name'], file_paths=config)

# Read the samples.txt file. This defines barcodes, the studied protein, ect.
df = pandas.read_csv(ex.file_paths['samples'], sep='\t')

df['Gene'] = [re.sub(' ', '-', x) for x in df['Gene']]  # Get rid of spaces.

# Define the samples list used throughout the workflow.
samples = [f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}" for exp,protein,l5_bc,l3_bc,rep in zip(
    df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'])]

proteins = list(set(df['Gene'])) 

# Three levels of peak cutoffs.
peak_cutoffs = [1, 2, 3]

# Constrain sample wildcards to not contain '/'.
wildcard_constraints:
    sample = "[^\/\.]+",
    protein= "[^\/\.]+",
    peak_cutoff="[123]+",
    strand="[\+-]",

# Read the sample sheet into the exp object.
ex.read_scheme(config['samples'])

# Define these global variables used numerous times in the workflow.
TOP_DIR = config['top_dir'].rstrip('/')
SAMS_DIR = config['sams'].rstrip('/')  # Holds both bam and sam, as well as STAR output logs.
FASTQ_DIR = config['fastq'].rstrip('/')  # The directory to put fastq intermediates in.
BIGWIG = config['bigwig'].rstrip('/')  # Holds bigwig outputs.

# Path to the clipper excecutable.
CLIPPER_PATH = "~/anaconda3/envs/clipper3/bin/clipper"
if 'clipper' in config:
    CLIPPER_PATH = config['clipper']
    
def protein_in_fname(fname):
    if '/' in fname:
        base = fname.split('/')[-1]
        return '_'.join(base.split('_')[1:-3])
    return '_'.join(fname.split('_')[1:-3])

def control_bigwigs():
    _df = pandas.read_csv("random_controls/samples.txt", sep='\t')
    return [f"random_controls/bigwig/{exp}_{gene}_{l5}_{l3}.bigwig" for \
            exp,gene,l5,l3 in zip(_df.Experiment, _df.Gene, _df.L5_BC, _df.L3_BC)]

def control_bigwigs_3prime():
    _df = pandas.read_csv("random_controls/samples.txt", sep='\t')
    plus = [f"random_controls/bigwig/3prime/{exp}_{gene}_{l5}_{l3}.+.bigwig" for \
            exp,gene,l5,l3 in zip(_df.Experiment, _df.Gene, _df.L5_BC, _df.L3_BC)]
    minus = [f"random_controls/bigwig/3prime/{exp}_{gene}_{l5}_{l3}.-.bigwig" for \
            exp,gene,l5,l3 in zip(_df.Experiment, _df.Gene, _df.L5_BC, _df.L3_BC)]
    return plus + minus


#########################################################
# Load genome annotations.
#########################################################
# The version in data/processed/features.db has introns added.
db_fname = "data/processed/features.db"
db = gffutils.FeatureDB(db_fname, keep_order=True)


#########################################################
# Begin rules.
#########################################################

rule all:
    input:
        txpt = TOP_DIR + "/data/processed/normalized_signal/txpt/len_norm.dill",
        intron = TOP_DIR + "/data/processed/normalized_signal/intron/len_norm.dill",
        exon = TOP_DIR + "/data/processed/normalized_signal/exon/len_norm.dill",
        utr5 = TOP_DIR + "/data/processed/normalized_signal/utr5/len_norm.dill",
        utr3 = TOP_DIR + "/data/processed/normalized_signal/utr3/len_norm.dill",
        
        txptr = TOP_DIR + "/data/processed/normalized_signal/txpt/len_norm.rows_subset.norm_by_row.dill",
        intronr = TOP_DIR + "/data/processed/normalized_signal/intron/len_norm.rows_subset.norm_by_row.dill",
        exonr = TOP_DIR + "/data/processed/normalized_signal/exon/len_norm.rows_subset.norm_by_row.dill",
        utr5r = TOP_DIR + "/data/processed/normalized_signal/utr5/len_norm.rows_subset.norm_by_row.dill",
        utr3r = TOP_DIR + "/data/processed/normalized_signal/utr3/len_norm.rows_subset.norm_by_row.dill",
    run:
        shell("echo Completed!")

rule make_raw_signal_objects:
    input:
        plus = BIGWIG + "/3prime/{sample}.+.bigwig",
        minus = BIGWIG + "/3prime/{sample}.-.bigwig",
    output:
        txpt = TOP_DIR + "/data/processed/normalized_signal/txpt/raw_reads.{sample}.dill",
        intron = TOP_DIR + "/data/processed/normalized_signal/intron/raw_reads.{sample}.dill",
        exon = TOP_DIR + "/data/processed/normalized_signal/exon/raw_reads.{sample}.dill",
        utr5 = TOP_DIR + "/data/processed/normalized_signal/utr5/raw_reads.{sample}.dill",
        utr3 = TOP_DIR + "/data/processed/normalized_signal/utr3/raw_reads.{sample}.dill",
    run:
        #params = {'which_RNAs': 'set', 'how_norm': 'sum', 'one_txpt_per_gene': 'longest',
        #          'gene_names': set(["NOC2L", "DDX11L1", "SAMD11"]), 'top': 100, 'mature_rna': True}        
        params = {
            # Either all ('all'), pass a set ('set'), or take the top ('top').
            'which_RNAs': 'top',

            # Either no norm (just sum RPMs, 'none') or each RNA to total ('sum').
            'how_norm': 'sum',

            # If which_RNAs='set', need this param (ignored otherwise):
            'gene_names': set(),

            # If which_RNAs='top', need this param (ignored otherwise):
            'top': 100,
            
            'one_txpt_per_gene': 'longest',
            'mature_rna': True,
        }

        # Tidy version to copy/paste:
        #params = {'which_RNAs': 'all', 'how_norm': 'sum', 'gene_names': set(), 'top': 100}
        import pyBigWig, re, gffutils, os, sys, collections, pandas, glob, random
        import numpy as np
        import matplotlib.pyplot as plt
        import scipy
        import scipy.stats

        import scripts.bigwigToNormArray as bigwigToNormArray
        import scripts.regionsFromDB as regionsFromDB
        
        # The version in data/processed/features.db has introns added.
        db_fname = "data/processed/features.db"
        db_fname = "data/processed/features.utrs_fixed.db"
        db = gffutils.FeatureDB(db_fname, keep_order=True)
        
        #importlib.reload(bigwigToNormArray); importlib.reload(regionsFromDB)
        #from bigwigToNormArray import bigwigToNormArray

        def load_cov(bw_fnames, db, featuretype='transcript', subset=False, target_size=10, params={}):
            sa_txpt = bigwigToNormArray.signalToNormArray(bw_fnames, bw_uses_chr_in_chrom_names=True)
            feat_select = regionsFromDB.get_feat_select(db, featuretype, params, subset=subset)
            n_feat = len(feat_select)
            
            # coverages: dict of {bw_name: 2D numpy array}
            cov_txpt = sa_txpt.read_coverages_in_regions(
                feat_select, n_features=n_feat, target_size=target_size, regions_are_split=bool(
                'mature_rna' in params and params['mature_rna']==True))
            
            to_fname = {
                'five_prime_utr': str(output.utr5), 'three_prime_utr': str(output.utr3),
                'transcript': str(output.txpt), 'exon': str(output.exon), 'intron': str(output.intron)}
            
            with open(to_fname[featuretype], 'wb') as fh:
                dill.dump(cov_txpt, fh)
                
            return cov_txpt

        bw_fnames = [str(input.plus), str(input.minus)]
        params = {'which_RNAs': 'protein_coding', 'how_norm': 'sum', 'gene_names': set(), 'top': 100}

        utr5_cov = load_cov(bw_fnames, db, featuretype='five_prime_utr', subset=10000, target_size=100, params=params)
        utr3_cov = load_cov(bw_fnames, db, featuretype='three_prime_utr', subset=10000, target_size=100, params=params)
        txpt_cov = load_cov(bw_fnames, db, featuretype='transcript', subset=False, target_size=100, params=params)

        intron_cov = load_cov(bw_fnames, db, featuretype='intron', subset=10000, target_size=100, params=params)
        exon_cov = load_cov(bw_fnames, db, featuretype='exon', subset=10000, target_size=100, params=params)

rule combine_raw_signal_objects:
    input:
        txpt = expand(TOP_DIR + "/data/processed/normalized_signal/txpt/raw_reads.{sample}.dill", sample=samples),
        intron = expand(TOP_DIR + "/data/processed/normalized_signal/intron/raw_reads.{sample}.dill", sample=samples),
        exon = expand(TOP_DIR + "/data/processed/normalized_signal/exon/raw_reads.{sample}.dill", sample=samples),
        utr5 = expand(TOP_DIR + "/data/processed/normalized_signal/utr5/raw_reads.{sample}.dill", sample=samples),
        utr3 = expand(TOP_DIR + "/data/processed/normalized_signal/utr3/raw_reads.{sample}.dill", sample=samples),
    output:
        txpt = TOP_DIR + "/data/processed/normalized_signal/txpt/len_norm.dill",
        intron = TOP_DIR + "/data/processed/normalized_signal/intron/len_norm.dill",
        exon = TOP_DIR + "/data/processed/normalized_signal/exon/len_norm.dill",
        utr5 = TOP_DIR + "/data/processed/normalized_signal/utr5/len_norm.dill",
        utr3 = TOP_DIR + "/data/processed/normalized_signal/utr3/len_norm.dill",
    run:
        regions = ['txpt', 'intron', 'exon', 'utr5', 'utr3']
        for region in regions:
            cov = {}
            input_fnames = [
                TOP_DIR + f"/data/processed/normalized_signal/{region}/raw_reads.{sample}.dill" \
                for sample in samples]
            for _covs in [dill.load(open(x, 'rb')) for x in input_fnames]:
                cov.update(_covs)
                
            with open(TOP_DIR + f"/data/processed/normalized_signal/{region}/len_norm.dill", 'wb') as fh:
                dill.dump(cov, fh)

rule normalize_signal_objects:
    input:
        txpt = TOP_DIR + "/data/processed/normalized_signal/txpt/len_norm.dill",
        intron = TOP_DIR + "/data/processed/normalized_signal/intron/len_norm.dill",
        exon = TOP_DIR + "/data/processed/normalized_signal/exon/len_norm.dill",
        utr5 = TOP_DIR + "/data/processed/normalized_signal/utr5/len_norm.dill",
        utr3 = TOP_DIR + "/data/processed/normalized_signal/utr3/len_norm.dill",
    output:
        txpt = TOP_DIR + "/data/processed/normalized_signal/txpt/len_norm.rows_subset.norm_by_row.dill",
        intron = TOP_DIR + "/data/processed/normalized_signal/intron/len_norm.rows_subset.norm_by_row.dill",
        exon = TOP_DIR + "/data/processed/normalized_signal/exon/len_norm.rows_subset.norm_by_row.dill",
        utr5 = TOP_DIR + "/data/processed/normalized_signal/utr5/len_norm.rows_subset.norm_by_row.dill",
        utr3 = TOP_DIR + "/data/processed/normalized_signal/utr3/len_norm.rows_subset.norm_by_row.dill",
        #intron = expand(TOP_DIR + "/data/processed/normalized_signal/intron/len_norm.{protein}.dill", protein=proteins),
        #exon = expand(TOP_DIR + "/data/processed/normalized_signal/exon/len_norm.{protein}.dill", protein=proteins),
        #utr5 = expand(TOP_DIR + "/data/processed/normalized_signal/utr5/len_norm.{protein}.dill", protein=proteins),
        #utr3 = expand(TOP_DIR + "/data/processed/normalized_signal/utr3/len_norm.{protein}.dill", protein=proteins),
    run:
        # For smoothing.
        from scipy.ndimage.filters import uniform_filter1d

        def sort_2d_array_by_sum(arrays):
            b = np.sum(arrays, axis = 1)
            idx = b.argsort()[::-1]
            return arrays[idx,]

        def norm_by_row(arr):
            sum_of_rows = arr.sum(axis=1)
            return arr / sum_of_rows[:, np.newaxis]

        def subset_rows_by_cutoff(arr, c=1):
            sum_of_rows = arr.sum(axis=1)
            print(arr[sum_of_rows>c].shape)
            return arr[sum_of_rows>c]

        def combine_groups(cov: dict, groups: dict):
            """groups: {group_name: [bw_fname1, ...]}.
            cov: {bw_fname1: 2D array, ...}
            Returns: {group_name: 2D array of average, ...}
            """
            means = {}
            #to_ave = {group: [x for x in cov if x in group] for group in groups}
            for group, bw_fnames in groups.items():
                bw_fname = bw_fnames[0]
                if any([bool(x not in cov) for x in bw_fnames]):
                    print(f"Tried to average {bw_fnames} in coverage object but at least one was missing.")
                    print(f"Coverage: {cov}.\n\n keys in cov object: {cov.keys()}")
                    continue
                means[group] = np.zeros(cov[bw_fname].shape)
                for row_n in range(0, len(means[group])):
                    means[group][row_n] = np.nanmean([cov[bw][row_n,:] for bw in bw_fnames], axis=0)
            return means

        def norm_and_write_object(cov, groups, output_dill, out_figname):

            mean_cov_txpt = combine_groups(cov, groups)
            # Plot smoothed coverage:
            #fig, ax = plt.subplots(figsize=(20,10))
            top_n = 500
            processed = {}
            for col in mean_cov_txpt:
                #plt.plot(sampleMeans[col], label=col, alpha=0.5)
                a = mean_cov_txpt[col]
                a = sort_2d_array_by_sum(a)[:top_n]
                print(a.shape)
            #    a = subset_rows_by_cutoff(a)
                a = norm_by_row(mean_cov_txpt[col])
                x = np.nansum(a, axis=0)
                x = x/np.nansum(x)
                
                processed[col] = x
                
            with open(output_dill, 'wb') as fh:
                dill.dump(processed, fh)

            for col, x in processed.items():
                N = 5
                plt.plot(uniform_filter1d(x, size=N), label=col, alpha=0.5, linewidth=4)
            plt.legend()
            plt.savefig(out_figname.split('.pdf')[0] + f".top_{top_n}.{N}_smooth.pdf")
            plt.show(); plt.clf()
            
        groups = {protein_in_fname(s): [] for s in samples}
        for label in groups:
            groups[label] = [x for x in samples if protein_in_fname(x)==label]
        regions = ['txpt', 'intron', 'exon', 'utr5', 'utr3']
        for region in regions:
            cov = dill.load(open(TOP_DIR + f"/data/processed/normalized_signal/{region}/len_norm.dill", 'rb'))
            #for group, samples in groups.items():
                #input_fnames = [
                #    TOP_DIR + f"/data/processed/normalized_signal/txpt/raw_reads.{sample}.dill" \
                #    for sample in samples]
            out_figname = TOP_DIR + f"/outs/figs/norm_signal/{region}.pdf"
            os.makedirs(os.path.dirname(out_figname), exist_ok=True)
            output_dill = TOP_DIR + f"/data/processed/normalized_signal/{region}/len_norm.rows_subset.norm_by_row.dill"
            norm_and_write_object(cov, groups, output_dill, out_figname)
#'''
