used_sequencing_runs = {

    'l1': {'Instrument': 'Illumina Hiseq 2500', 'read length': 100,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane1/fastq/raw/R2.fastq.gz',
    },

    'l2': {'Instrument': 'Illumina Hiseq 2500', 'read length': 100,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/hiseq/200522_medgenome/lane2/fastq/raw/R2.fastq.gz',
    },

    'm0': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200326/fastq/raw/R2.fastq.gz',
    },

    'm1': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200420/fastq/raw/R2.fastq.gz',
    },

    'm2': {'Instrument': 'Illumina Miseq', 'read length': 76,
    'location': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/',
    'R1_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/fastq/raw/R1.fastq.gz',
    'R2_fastq': '/oak/stanford/groups/khavari/users/dfporter/seq/miseq/Runs/200425/fastq/raw/R2.fastq.gz',
    },
}

import os, sys, re, glob, pandas

import exp
import metaExp
from pprint import pprint as pp
import numpy as np

def to_R2(fname: str) -> str:
    """Return the second mate pair filename given the R1 mate filename"""
    b = os.path.basename(fname) if '/' in fname else fname
    b = b.split('_')
    b = '_'.join([b[0], b[1], b[2], 'R2', b[3]])
    if '/' in fname:
        b = os.path.dirname(fname) + '/' + b
    return b

def read_a_read(fh):
    
    try:
        name = next(fh).rstrip('\n')
    except:
        return False, False, False
    
    seq = next(fh).rstrip('\n')
    next(fh)  # The ever-useful '+' line in fastq files.
    qual = next(fh).rstrip('\n')

    return name, seq, qual

def concat(fastq_file_r1, fastq_file_r2, out_read1_fname, out_read2_fname):
    
    fastq_fh_r2 = open(fastq_file_r2, 'r')
    fastq_fh_r1 = open(fastq_file_r1, 'r')
        
    out_read1 = open(out_read1_fname, 'a')
    out_read2 = open(out_read2_fname, 'a')
    
    print(f"Writing to \n{out_read1}\n{out_read2}")
    
    while True:
        name1, seq1, qual1 = read_a_read(fastq_fh_r1)
        name2, seq2, qual2 = read_a_read(fastq_fh_r2)
        
        if not name1:  # End of file.
            break
        
        #prefix = re.search('(@[\w-]+rand=[ATCGN]+)-', name1.split(' ')[0]).groups()[0]
        #name1 = prefix + ' ' + name1.split(' ')[1]
        #name2 = prefix + ' ' + name2.split(' ')[1]
        
        out_read1.write(f'{name1}\n{seq1}\n+\n{qual1}\n')
        out_read2.write(f'{name2}\n{seq2}\n+\n{qual2}\n')


def proclaim(cmd):
    print(cmd)
    os.system(cmd)


def extract_medgenome_sequencing_and_miseq(meta, out_top_dir='combined/'):
    """lane1, lane2, miseq0, miseq1, miseq2"""
    
    to_combine = {}  # {protein: {basename -> exp_name}}
    filenames_to_concat = {}  # {basename -> (R1 file path, R2 file path)}
    
    for exp_name in ['l1', 'l2', 'm0', 'm1', 'm2']:
        df = meta.exps[exp_name].scheme.scheme_df
        
        for protein, long_fname_R1, long_fname_R2 in zip(
            df.Gene, df['long_fname_R1'], df['long_fname_R2']):
            
            to_combine.setdefault(protein, {})
            to_combine[protein].setdefault(long_fname_R1, [])
            filenames_to_concat.setdefault(long_fname_R1, [])
            
            to_combine[protein][long_fname_R1].append(exp_name)
            filenames_to_concat[long_fname_R1].append(
                (meta.exps[exp_name].file_paths['cutadapt'] + '/split/' + long_fname_R1,
                meta.exps[exp_name].file_paths['cutadapt'] + '/split/' + long_fname_R2)
                )
            
    pp(to_combine)

    os.makedirs(f"{out_top_dir}/prezip/", exist_ok=True)
    os.makedirs(f"{out_top_dir}/zip/", exist_ok=True)

    for protein in to_combine:
        for long_fname_R1 in to_combine[protein]:
            miseqs = [x for x in to_combine[protein][long_fname_R1] if x[0]=='m']
            hiseqs = [x for x in to_combine[protein][long_fname_R1] if x[0]=='l']
            b = long_fname_R1.split('_')
            
            b[2] = 'mi' + b[2]
            miseqs_r1 = f"{out_top_dir}/prezip/{'_'.join(b)}"
            miseqs_r2 = to_R2(miseqs_r1)
            b[2] = 'hi' + b[2][2:]
            hiseqs_r1 = f"{out_top_dir}/prezip/{'_'.join(b)}"
            hiseqs_r2 = to_R2(hiseqs_r1)

            if len(miseqs):
                [os.system(f'touch {fname}') for fname in [miseqs_r1, miseqs_r2]]
            if len(hiseqs):
                [os.system(f'touch {fname}') for fname in [hiseqs_r1, hiseqs_r2]]

            for exp_name in miseqs:
                ca = meta.exps[exp_name].file_paths['cutadapt'] + '/split/'
                concat(
                    f'{ca}/{long_fname_R1}', f"{ca}/{to_R2(long_fname_R1)}",
                    miseqs_r1, miseqs_r2)
            for exp_name in hiseqs:
                ca = meta.exps[exp_name].file_paths['cutadapt'] + '/split/'
                concat(
                    f'{ca}/{long_fname_R1}', f"{ca}/{to_R2(long_fname_R1)}",
                    hiseqs_r1, hiseqs_r2)

            if len(miseqs):
                cmd = f'gzip < {miseqs_r1} > '
                cmd += f'{out_top_dir}/zip/{os.path.basename(miseqs_r1)}.gz'
                proclaim(cmd)
                cmd = f'gzip < {miseqs_r2} > '
                cmd += f'{out_top_dir}/zip/{os.path.basename(miseqs_r2)}.gz'
                proclaim(cmd)

            if len(hiseqs):
                cmd = f'gzip < {hiseqs_r1} > '
                cmd += f'{out_top_dir}/zip/{os.path.basename(hiseqs_r1)}.gz'
                proclaim(cmd)
                cmd = f'gzip < {hiseqs_r2} > '
                cmd += f'{out_top_dir}/zip/{os.path.basename(hiseqs_r2)}.gz'
                proclaim(cmd)

def make_samples_table_for_geo_upload(dirname):
    fnames = [os.path.basename(x) for x in glob.glob(dirname + '/*fastq*')]
    geo_upload_columns = [
    'Sample name', 'title', 'source name', 'organism',
    'characteristics: Expression',
    'characteristics: Exp',
    'characteristics: L3-BC', 
    'characteristics: L5-BC',
    'molecule', 'description', 'processed data file', 'raw file']
    r1_fnames = [x for x in fnames if '_R2' not in x]
    r2_fnames = [to_R2(x) for x in r1_fnames]
    print(fnames)
    print(r1_fnames)
    print(r2_fnames)
    geo = pandas.DataFrame(columns=geo_upload_columns, index=range(len(r1_fnames)))

def write_table_of_insert_lengths(dirname):
    fnames = [x for x in glob.glob(dirname + '/*fastq*')]
    r1_fnames = [x for x in fnames if ('_R2' not in os.path.basename(x))]
    df = pandas.DataFrame(index=range(len(r1_fnames)), columns=[
        'file name 1', 'file name 2', 'average insert size', 'standard deviation'])
    df['file name 1'] = r1_fnames
    df['file name 2'] = [to_R2(x) for x in r1_fnames]

    def get_insert_lengths(filename):
        lengths = []
        print(f"Getting insert lengths from {filename}...")

        with open(filename, 'r') as fh:
            lines = []
            for line in fh:
                lines.append(line.rstrip())
                if len(lines) == 4:
                    lengths.append(len(lines[1]))
                    lines = []

        arr = np.array(lengths)
        mean, stddev = arr.mean(), np.var(arr)**0.5
        mean, stddev = f"{mean:.2f}", f"{stddev:.2f}" 
        return mean, stddev 
    
    tups = [get_insert_lengths(x) for x in r1_fnames]
    df['average insert size'] = [x[0] for x in tups]
    df['standard deviation'] = [x[1] for x in tups]

    df.to_csv('insert_sizes.csv')

def run():
    seq_dir = "/oak/stanford/groups/khavari/users/dfporter/seq/"
    if not os.path.exists(seq_dir):
        raise IOError(f"input seq directory does not exist: {seq_dir}")
        
    meta = metaExp.metaExp(
        file_paths={'top_dir': f'{seq_dir}/2020_seq_for_geo/'})

    for name, paths in used_sequencing_runs.items():
        meta.add_exp(
            exp.exp(
                name=name,
                file_paths=exp.exp.autogenerate_paths(paths['location'])
            )
       )
        meta.exps[name].read_scheme()

    extract_medgenome_sequencing_and_miseq(
        meta, out_top_dir=f'{seq_dir}/2020_seq_for_geo/')
    write_table_of_insert_lengths(f'{seq_dir}/2020_seq_for_geo/prezip/')
    
    #meta.combine_scheme_files()

if __name__ == '__main__':
    run()


                                                                                                                                                                                                                                            
    
