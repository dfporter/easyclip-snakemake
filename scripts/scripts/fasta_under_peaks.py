import pandas, os, re, HTSeq, argparse, random, glob
#sys, typing, matplotlib,  , pysam, collections, subprocess, pyBigWig, re, scipy, importlib
#import seaborn as sns
#import numpy as np


def load_peaks(macs2_dir):

    dfs = {}
    for fname in glob.glob(f'{macs2_dir}/*narrow*'):
        base = os.path.basename(fname)

        try:
            df = pandas.read_csv(fname, sep='\t', header=None, skiprows=1)
        except:
            print(f"Parsing error for {fname}. Skipped.")
            continue
            
        df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'enrich', 'p', 'q', 'peak_offset']
        df['iv'] = list(zip(df['chrom'], df['start'], df['end']))
        # Drop weird chromosomes.
        df = df.loc[[ bool(re.match('chr[0-9YX]+', c, re.IGNORECASE) is not None) for c in df['chrom'] ], :]
        #df = df.sample(n=min([df.shape[0], 10000]), replace=False)
        #assign_strand_and_coverage(df, base)
        
        dfs[os.path.basename(fname)] = df
        
    return dfs

########################################################################
# Add the genomic sequence under peaks and write:
# 1. The annotated peaks files.
# 2. Fasta files for motif enrichment calculations with homer.
########################################################################


def add_seq(iv):
    if iv[0] in genomic_fasta:
        _seq = genomic_fasta[iv[0]][int(iv[1]):int(iv[2])]
        if len(iv) > 3 and iv[3] == '-':
            _seq = rc(_seq)
        return _seq
    else:
        return ''

    
def write_fastas(dfs, genomic_fasta, out_folder):
    os.makedirs(out_folder, exist_ok=True)

    for k, df in dfs.items():
        df['iv'] = list(zip(df['chrom'], df['start'], df['end'], df['strand']))
        df['seq'] = [add_seq(iv) for iv in df.iv]

        df['seq_len'] = [len(x) for x in df.seq]
        dfs[k] = df.loc[df['seq_len']>0, :]

        outf = f"{out_folder}/{k.split('_peaks')[0]}.fa"
        with open(outf, 'w') as f:
            for iv, seq in zip(df.iv, df.seq):
                seq_name = "_".join([str(x) for x in iv])
                f.write(f'>{seq_name}\n{seq}\n')

    ann_peaks_dir = f'{out_folder}/annPeaks/'
    os.makedirs(ann_peaks_dir, exist_ok=True)
    for fname in dfs:
        dfs[fname].to_excel(f"{ann_peaks_dir}/{fname}.xlsx")

        
def write_random_seqs(seqs, outfname, minimum=20000):

    alphabet = ['A', 'T', 'C', 'G']

    n_seqs = len(seqs) if len(seqs) > minimum else minimum

    if len(seqs) == 0:
        return ''
    
    with open(outfname, 'w') as outfh:
        for n, seq in enumerate(seqs):
            outlines = ">{}\n".format(n)
            outlines += "".join(random.choices(alphabet, k=len(seq))) + '\n'
            outfh.write(outlines)
    
    return outlines


def load_fasta(genomic_fasta_fname):
    return dict( (s[1], s[0]) for s in HTSeq.FastaReader(genomic_fasta_fname, raw_iterator=True) )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract DNA sequences for a folder of MACS2 peaks files.')

    parser.add_argument('--fasta', help='Genomic fasta file', type=str)

    parser.add_argument('--peaks_dir', help='Directory containing peaks files.', type=str)

    parser.add_argument('--out_folder', help = 'Directory to write fasta files.', type=str)

    args = parser.parse_args()
    
    dfs = load_peaks(args.peaks_dir)

    genomic_fasta = load_fasta(args.fasta)
    
    write_fastas(dfs, genomic_fasta, args.out_folder)
    
    random_seqs_dir = f"{args.out_folder}/randomControls/"
    os.makedirs(random_seqs_dir, exist_ok=True)
    for name, df in dfs.items():
        write_random_seqs(list(df['seq']), f"{random_seqs_dir}/{name.split('_peaks')[0]}.fa")    
