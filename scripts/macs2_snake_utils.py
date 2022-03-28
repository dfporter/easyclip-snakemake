import os, sys, re, glob, pandas, importlib, shutil, math, pysam, pickle, collections, HTSeq, random
import numpy as np

def load_narrowPeaks(fname):

    column_names = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'enrich', 'p', 'q', 'point_source']
    
    try:
        df = pandas.read_csv(fname, sep='\t')
    except:
        df = pandas.DataFrame([{
            'chrom': '1', 'start': 10, 'end': 110, 'name': 'no_peaks', 'score': 0, 'strand': '+', 
            'enrich': 0, 'p': 0, 'q': 0, 'point_source': 70}])
        print(f"Couldn't load {fname}. Making an dummy dataframe instead.")

    # Raw MACS2 output has no header.
    df.columns = column_names[:len(df.columns)]
    
    # To match with the sequence in the fasta file of sequences under peaks.
    df['fasta_seq_title'] = [f"{chrom}_{start}_{end}" for chrom, start, end in zip(df.chrom, df.start, df.end)]
    
    # To paste into IGV:
    df['IGV_location'] = [f"{chrom}:{start}-{end}" for chrom, start, end in zip(df.chrom, df.start, df.end)]
    
    return df

def peaks_as_bdg(fname, bdg_fname):
    base = os.path.basename(fname).split(".narrowPeak")[0]
    df = load_narrowPeaks(fname)
    ga = HTSeq.GenomicArray('auto', stranded=False)
    for chrom, start, end, q in zip(df.chrom, df.start, df.end, df.q):
        ga[HTSeq.GenomicInterval(chrom, start, end)] += q
    ga.write_bedgraph_file(bdg_fname)
        
def load_fasta(fasta_fname):
    return dict( (s.name, s) for s in HTSeq.FastaReader(fasta_fname) )

def write_random_seqs(seqs, outfname, minimum=20000):

    alphabet = ['A', 'T', 'C', 'G']

    n_seqs = len(seqs) if len(seqs) > minimum else minimum
    
    outlines = []
    
    print(f'write_random_seqs: len(seqs)={len(seqs)}. n_seqs={n_seqs}')
    
    if len(seqs) > minimum:
        for n, input_seq in enumerate(seqs):
            seq = "".join(random.choices(alphabet, k=len(input_seq)))
            outlines.append(f">{n}\n{seq}\n")

    elif len(seqs) > 9:
        ave_len = np.median([len(x) for x in seqs])
        for n in range(n_seqs):
            seq = "".join(random.choices(alphabet, k=int(ave_len)))
            outlines.append(f">{n}\n{seq}\n")
            
    else:
        for n in range(10):
            seq = "".join(random.choices(alphabet, k=200))
            outlines.append(f">{n}\n{seq}\n")        
            
    with open(outfname, 'w') as outfh:
        outfh.write(''.join(outlines))
        
    return outlines

def find_chrom(chrom, skipped_chroms, genome):
    chrom = str(chrom)
    if chrom not in genome:
        _chrom = 'chr' + chrom
        if _chrom not in genome:
            skipped_chroms.add(chrom)
            return False
        return _chrom
    return chrom

def add_seq(df, genome):
    skipped_chroms = set()
    seqs = []
    for chrom, start, end, strand in zip(df.chrom, df.start, df.end, df.strand):
        chrom = find_chrom(chrom, skipped_chroms, genome)
        if not chrom:
            seqs.append('N')
            continue
        if strand == '-':
            seq = genome[chrom][start:end].get_reverse_complement().seq.decode()
        else:
            seq = genome[chrom][start:end].seq.decode()
        seqs.append(seq)
    #print(seqs)
    df['seq'] = seqs
    print(f"skipped {skipped_chroms}")
    return df

def add_seq_transcriptome(df, genome):
    skipped_chroms = set()
    seqs = []
    for chrom, start, end in zip(df.chrom, df.start, df.end):
        chrom = find_chrom(chrom, skipped_chroms, genome)
        if not chrom:
            seqs.append('N')
            continue
        seq = genome[chrom][start:end].seq.decode()
        seqs.append(seq)
    #print(seqs)
    df['seq'] = seqs
    print(f"skipped {skipped_chroms}")
    return df