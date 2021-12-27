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

import os, sys, re, glob, pandas, importlib, shutil, math, pysam, pickle, collections
import numpy as np
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

#PROTEINS = list(set(df['Gene'])) 

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


MERGED_BAMS_DIR = SAMS_DIR + "/merged"
MERGED_BIGWIG = BIGWIG + '/merged'
BEDGRAPH = ex.file_paths['bedgraph'].rstrip('/') + '/merged'

# Path to the clipper excecutable.
CLIPPER_PATH = "~/anaconda3/envs/clipper3/bin/clipper"
if 'clipper' in config:
    CLIPPER_PATH = config['clipper']

#########################################################
# Custom merge settings.
#########################################################

custom_merge = pandas.read_csv(config['custom_merge'], sep='\t')

PROTEINS = list(set(custom_merge['Name']))
sample_to_set_of_custom_merges = dict(custom_merge.groupby("sam_path")['Name'].apply(set))

#sample_to_set_of_custom_merges.update({
#    os.path.basename(k).split('.bam')[0]: v for k,v in sample_to_set_of_custom_merges.items()})

def bams_for_set_in_combination_sheet(wildcards):
    sam_paths = [sam for sam in custom_merge['sam_path'] if (
        wildcards.protein in sample_to_set_of_custom_merges[sam])]
    return sam_paths

def protein_in_fname(fname):
    if '/' in fname:
        base = fname.split('/')[-1]
        return '_'.join(base.split('_')[1:-3])
    return '_'.join(fname.split('_')[1:-3])

def control_bigwigs():
    _df = pandas.read_csv("random_controls/samples.txt", sep='\t')
    return [f"random_controls/bigwig/{exp}_{gene}_{rep}_{l5}_{l3}.bigwig" for \
            exp,gene,rep,l5,l3 in zip(_df.Experiment, _df.Gene, _df.Replicate, _df.L5_BC, _df.L3_BC)]

def control_bigwigs_3prime():
    _df = pandas.read_csv("random_controls/samples.txt", sep='\t')
    plus = [f"random_controls/bigwig/3prime/{exp}_{gene}_{rep}_{l5}_{l3}.+.bigwig" for \
            exp,gene,rep,l5,l3 in zip(_df.Experiment, _df.Gene, _df.Replicate, _df.L5_BC, _df.L3_BC)]
    minus = [f"random_controls/bigwig/3prime/{exp}_{gene}_{rep}_{l5}_{l3}.-.bigwig" for \
            exp,gene,rep,l5,l3 in zip(_df.Experiment, _df.Gene, _df.Replicate, _df.L5_BC, _df.L3_BC)]
    return plus + minus

def bams_for_protein(wildcards):
    samples_of_protein = [sample for sample in samples if protein_in_fname(sample)==wildcards.protein]
    bam_fnames = [SAMS_DIR + f"/dedup/{sample}.bam" for sample in samples_of_protein]
    return bam_fnames

    


#########################################################
# Begin rules.
#########################################################

rule all:
    input:
        config['data']+"/custom_merge_total_read_numbers.txt",
        expand(BEDGRAPH + "/{protein}.+.wig", protein=PROTEINS),
        expand(BEDGRAPH + "/{protein}.-.wig", protein=PROTEINS),
        expand(MERGED_BIGWIG + "/3prime/{protein}.+.bigwig", protein=PROTEINS),
        expand(MERGED_BIGWIG + "/3prime/{protein}.-.bigwig", protein=PROTEINS),
        #expand(MERGED_BIGWIG + "/{protein}.bigwig", protein=PROTEINS),
        expand(MERGED_BIGWIG + "/rpm/{protein}.bigwig", protein=PROTEINS),
        merged_bam = expand(MERGED_BAMS_DIR + "/{protein}.bam", protein=PROTEINS),
    shell:
        "echo Completed!"
        
rule merge_bams:
    input:
        bams_to_merge = bams_for_set_in_combination_sheet
    output:
        merged_bam = MERGED_BAMS_DIR + "/{protein}.bam"
    run:
        #"echo samtools merge -o {output.merged_bam} " + " ".join(input.bams_to_merge)
        shell("samtools merge {output.merged_bam} " + " ".join(input.bams_to_merge))
        shell("samtools sort -o {output.merged_bam}.sort {output.merged_bam}")
        shell("mv {output.merged_bam}.sort {output.merged_bam}")
        shell("samtools index {output.merged_bam}")
        
#########################################################
# Mapped read format conversions.
#########################################################            

rule bedGraphToBigWig:
    input:
        "data/processed/chrom.sizes",
        plus = BEDGRAPH + "/{protein}.+.wig",  # Bedgraph. wig extension so IGV loads.
        minus = BEDGRAPH + "/{protein}.-.wig",  # Bedgraph. wig extension so IGV loads.
    output:
        plus = MERGED_BIGWIG + "/3prime/{protein}.+.bigwig",
        minus = MERGED_BIGWIG + "/3prime/{protein}.-.bigwig",
    run:
        # Assuming input is sorted.
        shell("bedSort {input.plus} {input.plus}")
        shell("bedSort {input.minus} {input.minus}")
        shell("bedGraphToBigWig {input.plus} data/processed/chrom.sizes {output.plus}")
        shell("bedGraphToBigWig {input.minus} data/processed/chrom.sizes {output.minus}")
"""
subprocess.check_output(['samtools', 'idxstats', bamfile])
for li in c.decode().split('\n'):
    s = li.split('\t')
    n_reads += int(s[3])
    n_reads += int(s[4])
"""
rule total_read_numbers:
    input:
        plus = expand(BEDGRAPH + "/{protein}.+.wig",  protein=PROTEINS),
        minus = expand(BEDGRAPH + "/{protein}.-.wig", protein=PROTEINS),
    output:
        total_reads = config['data']+"/custom_merge_total_read_numbers.txt",
    run:
        import scripts.total_read_numbers
        scripts.total_read_numbers.total_read_numbers(
            folder=BEDGRAPH, outfile=str(output.total_reads))
        
rule convert_to_3prime_bedgraph:
    """5' end of read 1 -> 3' end of cDNA. This outputs the 3' end of the cDNA molecule."""
    input:
        MERGED_BAMS_DIR + "/{protein}.bam"
    output:
        plus =  BEDGRAPH + "/{protein}.+.wig",  # Bedgraph. wig extension so IGV loads.
        minus = BEDGRAPH + "/{protein}.-.wig",  # Bedgraph. wig extension so IGV loads.
    run:
        shell("bedtools genomecov -bg -strand + -5 -ibam {input} > {output.plus}")
        shell("bedSort {output.plus} {output.plus}")
        shell("bedtools genomecov -bg -strand - -5 -ibam {input} > {output.minus}")
        shell("bedSort {output.minus} {output.minus}")
        
def scale_factor(wildcards):
    if os.path.exists(config['data']+"/custom_merge_total_read_numbers.txt"):
        df = pandas.read_csv(config['data']+"/custom_merge_total_read_numbers.txt", sep='\t')
    else:
        print(f"scale_factor() called before " + config['data']+"/custom_merge_total_read_numbers.txt" + \
              " created. Shouldn't happen after initial function compilation.")
        return "not_defined"
    df["basename"] = [x.split('/')[-1].split('.')[0] for x in df.Dataset]
    return str(1000000/df.loc[df['basename']==wildcards.protein, 'Total read number'].iloc[0])
            
rule bigWigToScaledWig:
    input:
        bigwig = MERGED_BIGWIG + "/{protein}.bigwig",
        total_reads = config['data'] + "/merged_total_read_numbers.txt",
    output:
        wig = TOP_DIR + "/scaledWig/{protein}.bigwig"
    params:
        scaling_factor = scale_factor
    conda:
        "envs/wiggletools.yml"
    shell:
        "wiggletools write {output.wig} scale {params.scaling_factor} {input.bigwig}"

rule scaledWigToBigWig:
    input:
        wig = TOP_DIR + "/scaledWig/{protein}.bigwig",
        chrom_sizes = "data/processed/chrom.sizes"
    output:
        bigwig = MERGED_BIGWIG + "/rpm/{protein}.bigwig"
    shell:
        "wigToBigWig {input.wig} {input.chrom_sizes} {output.bigwig}"
        
rule bamToBigwig:
    input:
        bams = MERGED_BAMS_DIR + "/{protein}.bam"
    output:
        bigwig = MERGED_BIGWIG + "/{protein}.bigwig"
    shell:
        #"bamCoverage --binSize 10 -b {input} -o {output} -of bigwig"
        "bamCoverage --binSize 1 -b {input} -o {output} -of bigwig"