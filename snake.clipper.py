'''
# Clipper installation:
git clone https://github.com/YeoLab/clipper
cd clipper
conda env create -f environment3.yml
conda activate clipper3
python setup.py install

# Replace the samtools import version in the environment3.yml file with 1.6
# (or some other version that works) to prevernt an error message whem samtools called.

# Prevernts error.
conda install -c bioconda samtools=1.6 --force-reinstall
'''

import os, sys, re, glob, pandas, importlib, shutil
#from RepEnrich2 import RepEnrich2
import scripts
import scripts.exp
import scripts.rnaDataFileMaker
import scripts.make_repeats_chrom
import scripts.split_bam
import scripts.random_sequence

# Load config file.
configfile: "config.rb.yaml"


#########################################################
# Set global path values and initialize exp object.
#########################################################

# The config.yaml must contain a name definition.
if ('top_dir' not in config) and ('name' not in config):
    raise ValueError(
        f"Need name value in {configfile}. Prefer name & top_dir both set. top_dir=name when only name set.")

# Set top_dir from name if not given.    
if ('top_dir' not in config) and ('name' in config):
    config['top_dir'] = config['name']
    os.makedirs(config['top_dir'], exist_ok=True)

# A star index is required.   
if ('STAR_index' in config) and ('STAR index' not in config):
    config['STAR index'] = config['STAR_index']

# Set some default path values as relative to the top_dir specified in config.yaml.
_to = lambda x: config['top_dir'].rstrip('/') + f'/{x}'
defaults = {
    'samples': _to('samples.txt'), 'beds': _to('beds'), 'fastq': _to('fastq'), 'sams': _to('sams'),
    'bigwig': _to('bigwig'), 'bedgraph': _to('bedgraph'), 'counts': _to('counts.txt'), 'STAR': 'STAR',
    'outs': _to('outs/'), 'counts': _to('outs/counts/'),
    'scheme': _to('samples.txt'),
    'run_clipper': 'false',  # Lower case false because this is used as a string in a bash command.
}

for key in [_ for _ in defaults.keys() if (_ not in config)]:
    config[key] = defaults[key]
    
# Create exp object for organizing some tasks.
ex = scripts.exp.exp(name=config['name'], file_paths=config)

# If a feature_gtf is not specified in the config.yaml, we give it a default
# of assets/reference/gencode.v29.basic.annotation.gtf. The workflow will attempt
# to download this file if it is not present. If the star index was made with a 
# different annotation than gencode v29, some weird results might happen in
# assigning reads to genes.
if not os.path.exists(config['feature_gtf']):
    print("feature gtf not found: ", config['feature_gtf'])
    print("setting config['feature_gtf'] to match this default gtf: assets/reference/gencode.v29.basic.annotation.gtf")
    config['feature_gtf'] = "assets/reference/gencode.v29.basic.annotation.gtf"

# Read the samples.txt file. This defines barcodes, the studied protein, ect.
df = pandas.read_csv(ex.file_paths['samples'], sep='\t')

df['Gene'] = [re.sub(' ', '-', x) for x in df['Gene']]  # Get rid of spaces.

# Define the samples list used throughout the workflow.
samples = [f"{exp}_{protein}_{l5_bc}_{l3_bc}" for exp,protein,l5_bc,l3_bc in zip(
    df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'])]

proteins = list(set(df['Gene'])) 

# Constrain sample wildcards to not contain '/'.
wildcard_constraints:
    sample = "[^\/\.]+",
    protein= "[^\/\.]+",

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
    
#########################################################
# Begin rules.
#########################################################

rule all:
    input:
        "assets/reference/Dfam_curatedonly.embl",
        "assets/repeats_star_index/Genome",
        config['feature_gtf'],
        "assets/reference/featureCounts_formatted.gtf",
        "data/processed/features.db",
        "data/processed/repeats_and_longest_txpt_per_gene.bed",
        ex.file_paths['R1_fastq'],
        SAMS_DIR + '/all_reads.bam',
        FASTQ_DIR + '/umis_moved/R1.fastq.gz',
        FASTQ_DIR + '/ready_to_map/R1.fastq.gz',
        expand(SAMS_DIR + "/split/{sample}.bam", sample=samples),
        expand(SAMS_DIR + "/dedup/{sample}.bam", sample=samples),
        expand(SAMS_DIR + '/genome_only/{sample}.bam', sample=samples),
        bed = expand(TOP_DIR + '/peaks/clipper/{sample}.bed', sample=samples),
    run:
        shell("echo Completed!")

include: 'rules/references.smk'

#########################################################
# Run clipper.
#########################################################s
 
rule call_clipper:
    input:
        bam = SAMS_DIR + '/genome_only/{sample}.bam',
        bed12 = 'data/processed/repeats_and_longest_txpt_per_gene.bed',
    output:
        bed = TOP_DIR + '/peaks/clipper/{sample}.bed'
    conda:
        'clipper/environment3.yml'
    threads: 12
    shell:
        CLIPPER_PATH + " -b {input.bam} -o {output.bed} -s GRCh38_v29;"

        