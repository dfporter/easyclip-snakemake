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
import scripts.init_with_config
import scripts.exp
import scripts.rnaDataFileMaker
import scripts.make_repeats_chrom
import scripts.split_bam
import scripts.random_sequence

# Load config file.
configfile: "config.yaml"


#########################################################
# Set global path values and initialize exp object.
#########################################################

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
        SAMS_DIR + '/all_reads.bam',
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

        