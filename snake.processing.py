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

# Constrain sample wildcards to not contain '/' or '.'.
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
RAW_FASTQ_R1_INPUTS = [config['Fastq_folder'].rstrip('/') + f"/{R1}" for R1 in df['R1_fastq']]
RAW_FASTQ_R2_INPUTS = [config['Fastq_folder'].rstrip('/') + f"/{R2}" for R2 in df['R2_fastq']]
PCR_INDEX_SET = list(set([os.path.basename(x).split("R1.fastq.gz")[0] for x in  RAW_FASTQ_R1_INPUTS]))
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
        #expand(BIGWIG + "/{sample}.bigwig", sample=samples),
        expand(SAMS_DIR + "/3end/{sample}.bam", sample=samples),
        expand(BIGWIG + "/3prime/{sample}.+.bigwig", sample=samples),
        expand(BIGWIG + "/3prime/{sample}.-.bigwig", sample=samples),
        expand(SAMS_DIR + '/genome_only/{sample}.bam', sample=samples),
    run:
        shell("echo Completed!")

include: 'rules/references.smk'
include: 'rules/read_preprocessing_and_mapping.smk'


#########################################################
# Mapped read format conversions.
#########################################################            

rule bedGraphToBigWig:
    input:
        "data/processed/chrom.sizes",
        plus = ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.+.wig",  # Bedgraph. wig extension so IGV loads.
        minus = ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.-.wig",  # Bedgraph. wig extension so IGV loads.
    output:
        plus = BIGWIG + "/3prime/{sample}.+.bigwig",
        minus = BIGWIG + "/3prime/{sample}.-.bigwig",
    run:
        # Assuming input is sorted.
        shell("bedSort {input.plus} {input.plus}")
        shell("bedSort {input.minus} {input.minus}")
        shell("bedGraphToBigWig {input.plus} data/processed/chrom.sizes {output.plus}")
        shell("bedGraphToBigWig {input.minus} data/processed/chrom.sizes {output.minus}")

rule convert_to_3prime_bedgraph:
    """5' end of read 1 -> 3' end of cDNA. This outputs the 3' end of the cDNA molecule."""
    input:
        SAMS_DIR + "/dedup/{sample}.bam"
    output:
        plus = ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.+.wig",  # Bedgraph. wig extension so IGV loads.
        minus = ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.-.wig",  # Bedgraph. wig extension so IGV loads.
    run:
        shell("bedtools genomecov -bg -strand + -5 -ibam {input} > {output.plus}")
        shell("bedSort {output.plus} {output.plus}")
        shell("bedtools genomecov -bg -strand - -5 -ibam {input} > {output.minus}")
        shell("bedSort {output.minus} {output.minus}")
        
rule convert_bam_to_3prime_end_only:
    """Not using. Not sure if this is correct."""
    input:
        SAMS_DIR + "/dedup/{sample}.bam"
    output:
        bam = SAMS_DIR + "/3end/{sample}.bam"
    run:
        output_bam = str(output.bam)
        sam_output = os.path.splitext(output_bam)[0] + '.sam'
        shell("python scripts/convert_bam_to_RT_stop_only.py {input} {sam_output}")
        shell("samtools view -Sb {sam_output} > {output_bam}")
        shell("rm {sam_output}")
        shell("samtools sort -o {output_bam}.sort {output_bam}")
        shell("mv {output_bam}.sort {output_bam}")
        shell("samtools index {output_bam}")
        
rule remove_non_chromosomal_reads:
    input:
        bam = SAMS_DIR + '/dedup/{sample}.bam',
    output:
        bam = SAMS_DIR + '/genome_only/{sample}.bam',
    run:
        shell("samtools view -b {input.bam} chr{{1..22}} > {output.bam}")
        
       
rule bamToBigwig:
    input:
        bams = SAMS_DIR + "/dedup/{sample}.bam"
    output:
        bigwig = BIGWIG + "/{sample}.bigwig"
    shell:
        #"bamCoverage --binSize 10 -b {input} -o {output} -of bigwig"
        "bamCoverage --binSize 1 -b {input} -o {output} -of bigwig"
