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

import os, sys, re, glob, pandas, importlib, shutil, math, pysam, pickle
import numpy as np
import scripts
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
config = init_with_config.init_with_config(config)
   
# Create exp object for organizing some tasks.
ex = scripts.exp.exp(name=config['name'], file_paths=config)

# Read the samples.txt file. This defines barcodes, the studied protein, ect.
df = pandas.read_csv(ex.file_paths['samples'], sep='\t')

df['Gene'] = [re.sub(' ', '-', x) for x in df['Gene']]  # Get rid of spaces.

# Define the samples list used throughout the workflow.
samples = [f"{exp}_{protein}_{l5_bc}_{l3_bc}" for exp,protein,l5_bc,l3_bc in zip(
    df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'])]

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
        return '_'.join(base.split('_')[1:-2])
    return '_'.join(fname.split('_')[1:-2])

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
        #expand(BIGWIG + "/{sample}.bigwig", sample=samples),
        #expand(SAMS_DIR + "/3end/{sample}.bam", sample=samples),
        config['counts'].rstrip('/') + "/bigwig_3prime_counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/bam_3prime_counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/featureCounts_on_bams.txt",
        expand(SAMS_DIR + '/genome_only/{sample}.bam', sample=samples),

        expand(TOP_DIR + "/outs/homer/{sample}.2/{sample}.2", sample=samples),
        expand(TOP_DIR + "/peaks/clipper/filtered/fastas/{sample}.2.fa", sample=samples), 
        expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed', sample=samples),
        expand(TOP_DIR + "/outs/homer/merged.{protein}.1/{protein}.1", protein=proteins),
        expand(TOP_DIR + "/outs/homer/merged.{protein}.2/{protein}.2", protein=proteins),
        expand(TOP_DIR + "/outs/homer/merged.{protein}.3/{protein}.3", protein=proteins),
        
        expand(TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/merged.{protein}.1.bed", protein=proteins),
        expand(TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/merged.{protein}.2.bed", protein=proteins),
        expand(TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/merged.{protein}.3.bed", protein=proteins),
        
        expand(TOP_DIR + "/outs/homer/no_ncRNA_filtered.merged.{protein}.1/{protein}.1/homerMotifs.motifs12", protein=proteins),
        expand(TOP_DIR + "/outs/homer/no_ncRNA_filtered.merged.{protein}.2/{protein}.2/homerMotifs.motifs12", protein=proteins),
        expand(TOP_DIR + "/outs/homer/no_ncRNA_filtered.merged.{protein}.3/{protein}.3/homerMotifs.motifs12", protein=proteins),
        
        expand(TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/annotated/merged.{protein}.3.bed', protein=proteins),
        expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.3.bed', sample=samples, strand=["+", "-"]),
        
        rna_data = TOP_DIR + '/data/rna.data',
        counts = TOP_DIR + '/outs/counts/ann_counts.bedgraphs.txt',
        reads_per_million = TOP_DIR + '/outs/counts/reads_per_million.bedgraphs.txt',
        pvals_per_read = TOP_DIR + "/tables/pvals_per_read.xlsx",
        
        per_read_xlsx = TOP_DIR + "/tables/RNA targets per read vs random non-RBPs.xlsx",
        #per_protein_xlsx = TOP_DIR + "/tables/RNA targets per protein vs random non-RBPs.xlsx",
    run:
        shell("echo Completed!")

include: 'rules/references.smk'
    

rule make_ncRNA_gtf:
    input:
        gtf = "assets/reference/only_tsl_1_and_NA.gtf",
    output:
        ncRNA = "data/processed/ncRNA.only_tsl1_and_NA.bed",
        coding = "data/processed/mRNA.only_tsl1_and_NA.bed",
    run:
        df = pandas.read_csv(input.gtf, comment='#', sep='\t', header=None)
        df['biotype'] = [re.search('transcript_type "([^"]+)"', s) for s in df[8]]
        df['biotype'] = [str(s) if s==None else s.group(1) for s in df['biotype']]
        df['coding'] = [True if x=='protein_coding' else False for x in df.biotype]

        nc_df = df[~df['coding']]
        coding = df[df['coding']]

        nc_df.loc[:, [0, 3, 4, 8, 5, 6]].to_csv(output.ncRNA, sep='\t', header=None, index=None)
        coding.loc[:, [0, 3, 4, 8, 5, 6]].to_csv(output.coding, sep='\t', header=None, index=None)


def raw_peak_files_of_protein(wildcards):
    _sample_names = [x for x in samples if protein_in_fname(x).upper()==wildcards.protein.upper()]
    return [TOP_DIR + f'/peaks/clipper/{x}.bed' for x in _sample_names]

def assign_to_gene_and_add_height(bed_fname, bamfnames):        
    df = pandas.read_csv(bed_fname, sep='\t', header=None)
    df['Gene'] = [m.group(1) if (m := re.search('gene_name "([^"]+)"', x)) is not None else '.' for x in df[9]]
    df['len'] = [b-a for a,b in zip(df[1], df[2])]

    df = df[df['Gene']!='.']
    df = df[df['len']>0]

    bamfiles = {bam_fname:pysam.AlignmentFile(bam_fname, "rb" ) for bam_fname in bamfnames}

    for bam_fname, bam_fh in bamfiles.items():
        df[f'max_depth|{bam_fname}'] = [
            np.max(np.sum(bam_fh.count_coverage(chrm, start, end), axis=0)) for chrm,start,end in zip(df[0], df[1], df[2])
        ]

    [bamfh.close() for bamfh in bamfiles.values()]

    if len(bamfnames) > 1:
        cols = [x for x in df.columns if 'max_depth|' in str(x)]
        df['max_depth_all'] = df.loc[:, cols].apply(np.nanmean, axis=1)#[  np.nanmean(arr) for arr in df.loc[:,cols].values  ]
    else:
        df['max_depth_all'] = df[f'max_depth|{bamfnames[0]}']
        
    df = df.sort_values(by='max_depth_all', ascending=False)

    df = df.drop_duplicates(subset=[0, 1, 2, 'Gene'], keep='first')

    return df
     
def bam_files_of_protein_wildcards(wildcards):
    _sample_names = [x for x in samples if protein_in_fname(x).upper()==wildcards.protein.upper()]
    return [TOP_DIR + f'/sams/dedup/{x}.bam' for x in _sample_names]

def bam_files_of_protein(protein):
    _sample_names = [x for x in samples if protein_in_fname(x).upper()==protein.upper()]
    return [TOP_DIR + f'/sams/dedup/{x}.bam' for x in _sample_names]  


rule mRNA_peaks:
    input:
        bamfnames = bam_files_of_protein_wildcards,
        one = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.1.bed',
        two = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.2.bed', 
        three = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.3.bed', 
    output:
        one = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/annotated/merged.{protein}.1.bed',
        two = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/annotated/merged.{protein}.2.bed', 
        three = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/annotated/merged.{protein}.3.bed',    
    run:
        for input_bed, output_bed in zip(
            [input.one, input.two, input.three], [output.one, output.two, output.three]):
            
            input_bed, output_bed = str(input_bed), str(output_bed)
            
            # Produces multiple lines per peak when overlapping multiple genes.
            protein = output_bed.split('/')[-1].split('.')[1]
            basename = output_bed.split('/')[-1]
            os.makedirs('/tmp', exist_ok=True)
            shell(f"bedtools intersect -a {input_bed} -b data/processed/mRNA.only_tsl1_and_NA.bed -loj > tmp/{basename}")
            df = assign_to_gene_and_add_height(f"tmp/{basename}", bam_files_of_protein(protein))
            df.to_csv(output_bed, sep='\t', index=False)
        
rule remove_peaks_with_ncRNA:
    input:
        ncRNA = "data/processed/ncRNA.only_tsl1_and_NA.bed",
        one = TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.1.bed',
        two = TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.2.bed', 
        three = TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.3.bed', 
    output:
        one = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.1.bed',
        two = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.2.bed', 
        three = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.3.bed', 
    run:
        #os.makedirs(TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/', exist_ok=True)
        
        # -A: Remove entire feature in A if there is overlap with B.
        # -f: Fraction of A (per kb) overlapping with B required to call an overlap. 1E-9 = 1 bp.
        # -s: Stranded.
        shell("bedtools subtract -a {input.one} -b {input.ncRNA} -A -f 1E-8  > {output.one}")
        shell("bedtools subtract -a {input.two} -b {input.ncRNA} -A -f 1E-8 > {output.two}")
        shell("bedtools subtract -a {input.three} -b {input.ncRNA} -A -f 1E-8 > {output.three}")
        
rule write_fastas_merged_no_ncRNA:
    input:
        one = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.1.bed',
        two = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.2.bed', 
        three = TOP_DIR + '/peaks/clipper/no_ncRNA_filtered/merged.{protein}.3.bed', 
    output:
        peak_fa_1 = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/merged.{protein}.1.fa",
        peak_fa_2 = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/merged.{protein}.2.fa", 
        peak_fa_3 = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/merged.{protein}.3.fa", 
        randoms_1 = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/randomControls/merged.{protein}.1.fa",
        randoms_2 = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/randomControls/merged.{protein}.2.fa", 
        randoms_3 = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/randomControls/merged.{protein}.3.fa", 
    run:
        # -s: force strandedness.
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.one} -fo {output.peak_fa_1} -s")
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.two} -fo {output.peak_fa_2} -s")
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.three} -fo {output.peak_fa_3} -s")
        
        for peak_fasta, random_fa in [
            (output.peak_fa_1, output.randoms_1),
            (output.peak_fa_2, output.randoms_2),
            (output.peak_fa_3, output.randoms_3)]:
            
            if os.path.exists(peak_fasta):
                with open(peak_fasta) as f:
                    seqs = [x for x in f.readlines() if x[0]!='>']
                    scripts.random_sequence.write_random_seqs(seqs, random_fa)     

rule call_homer_merged_no_ncRNA:
    input:
        peak_fasta = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/merged.{protein}.{peak_cutoff}.fa",
        randoms = TOP_DIR + "/peaks/clipper/no_ncRNA_filtered/fastas/randomControls/merged.{protein}.{peak_cutoff}.fa",
    output:
        homer_results = TOP_DIR + "/outs/homer/no_ncRNA_filtered.merged.{protein}.{peak_cutoff}/{protein}.{peak_cutoff}/homerMotifs.motifs12",
    log:
        TOP_DIR + "/outs/homer/logs/no_ncRNA_filtered.{protein}.{peak_cutoff}.log"
    conda:
        "envs/homer.yml"
    params:
        folder = TOP_DIR + "/outs/homer/no_ncRNA_filtered.merged.{protein}.{peak_cutoff}/{protein}.{peak_cutoff}",
    shell:
        #"findMotifs.pl {input.peak_fasta} fasta <output_folder> -fasta {input.randoms} -rna -homer1"
        "findMotifs.pl {input.peak_fasta} fasta {params.folder} -fasta {input.randoms} -rna -homer1"
        " -mcheck data/processed/filtered_homer_known.rna.motifs"

        
#########################################################
# Run clipper.
#########################################################
    
rule call_homer:
    input:
        peak_fasta = TOP_DIR + "/peaks/clipper/filtered/fastas/{sample}.{peak_cutoff}.fa",
        randoms = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/{sample}.{peak_cutoff}.fa",
    output:
        homer_results = TOP_DIR + "/outs/homer/{sample}.{peak_cutoff}/{sample}.{peak_cutoff}",
    log: TOP_DIR + "/outs/homer/logs/{sample}.{peak_cutoff}.log"
    conda:
        "envs/homer.yml"
    shell:
        "homer2 denovo -i {input.peak_fasta} -b {input.randoms} -len 6,8,9,12 -S 10 -strand + -o {output.homer_results} &> {log}"

rule call_homer_merged:
    input:
        peak_fasta = TOP_DIR + "/peaks/clipper/filtered/fastas/merged.{protein}.{peak_cutoff}.fa",
        randoms = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/merged.{protein}.{peak_cutoff}.fa",
    output:
        homer_results = TOP_DIR + "/outs/homer/merged.{protein}.{peak_cutoff}/{protein}.{peak_cutoff}",
    log:
        TOP_DIR + "/outs/homer/logs/{protein}.{peak_cutoff}.log"
    conda:
        "envs/homer.yml"
    shell:
        "homer2 denovo -i {input.peak_fasta} -b {input.randoms} -len 6 -S 10 -strand + -o {output.homer_results} &> {log}"

rule split_peak_files_by_strand:
    input:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.1.bed', sample=samples),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed', sample=samples),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.3.bed', sample=samples),
    output:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.1.bed', sample=samples, strand=["+", "-"]),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.2.bed', sample=samples, strand=["+", "-"]),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.3.bed', sample=samples, strand=["+", "-"]),
    run:
        for sample in samples:
            for level in ['1', '2', '3']:  # For each level of stringency in peak calling.
                
                for strand in ['+', '-']:
                    in_file = TOP_DIR + f'/peaks/clipper/filtered/{sample}.{level}.bed'
                    
                    out_file = TOP_DIR + f'/peaks/clipper/filtered/{sample}.+.{level}.bed'
                    cmd = f"awk \'$6==\"+\"\' " + in_file + f" > {out_file}"
                    print(cmd)
                    os.system(cmd)

                    out_file = TOP_DIR + f'/peaks/clipper/filtered/{sample}.-.{level}.bed'
                    cmd = f"awk \'$6==\"-\"\' " + in_file + f" > {out_file}"
                    print(cmd)
                    os.system(cmd)
                    
rule intersect_peaks:  
    input:
        #one = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.1.bed', sample=samples),
        #two = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed', sample=samples),
        #three = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.3.bed', sample=samples),
        one = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.1.bed', sample=samples, strand=["+", "-"]),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.2.bed', sample=samples, strand=["+", "-"]),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.{strand}.3.bed', sample=samples, strand=["+", "-"]),
    output:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.{strand}.1.bed', protein=proteins, strand=["+", "-"]),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.{strand}.2.bed', protein=proteins, strand=["+", "-"]),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.{strand}.3.bed', protein=proteins, strand=["+", "-"]),
    run:
        for protein in proteins:  # For each CLIP'd protein.
            for level in ['1', '2', '3']:  # For each level of stringency in peak calling.
                
                for strand in ['+', '-']:
                    output_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.{strand}.{level}.bed'
                    _samples = [x for x in samples if protein_in_fname(x).upper()==protein.upper()]  # Replicate names.
                    print('---' * 14)
                    print(f'samples matching {protein}: {_samples}')

                    # Input replicate peak filenames.
                    _files = [TOP_DIR + f'/peaks/clipper/filtered/{_sample}.{strand}.{level}.bed' for _sample in _samples]
                    print(protein, _samples)

                    if len(_files) > 1:  # If multiple replicates for this protein...
                        thresh = math.ceil(0.75 * len(_files))  # Minimum replicate number.
                        filestring = " ".join(_files)
                        output_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.{strand}.{level}.bed'

                        # {{ is converted to { by snakemake (avoids snakemake interpreting {} itself).
                        cmdstring = f"bedtools multiinter -i {filestring} | awk \'{{{{ if($4 >= {thresh}) {{{{ print }}}} }}}}\' | bedtools merge > {output_file}"
                        shell(cmdstring)  # Produce merged file.

                    else:  # If only one replicate, just copy that as the merged file.
                        shell(f"cp {_files[0]} {output_file}")

rule merge_intersected_peaks_from_both_strands:
    input:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.{strand}.1.bed', protein=proteins, strand=["+", "-"]),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.{strand}.2.bed', protein=proteins, strand=["+", "-"]),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.{strand}.3.bed', protein=proteins, strand=["+", "-"]),
    output:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.1.bed', protein=proteins),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.2.bed', protein=proteins),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.3.bed', protein=proteins),
    run:
        for protein in proteins:  # For each CLIP'd protein.
            for level in ['1', '2', '3']:  # For each level of stringency in peak calling.
                                
                output_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.{level}.bed'
                
                input_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.+.{level}.bed'
                # cmd = f"awk \'$6==\"-\"\' " + in_file + f" > {out_file}"
                cmd = f"awk \'BEGIN{{FS=OFS=\"\\t\"}}{{print $0 OFS \".\" OFS \".\" OFS \"+\"}}\' " + \
                    f' {input_file} > {output_file}.+'
                print(cmd)
                os.system(cmd)

                input_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.-.{level}.bed'
                cmd = f"awk \'BEGIN{{FS=OFS=\"\\t\"}}{{print $0 OFS \".\" OFS \".\" OFS \"-\"}}\' " + \
                    f' {input_file} > {output_file}.-'
                print(cmd)
                os.system(cmd)
                
                cmd = f"cat {output_file}.+ {output_file}.- > {output_file}"
                print(cmd)
                os.system(cmd)
                
rule write_fastas:
    input:
        bed = TOP_DIR + '/peaks/clipper/filtered/{sample}.{peak_cutoff}.bed'
    output:
        peak_fastas = TOP_DIR + "/peaks/clipper/filtered/fastas/{sample}.{peak_cutoff}.fa", 
        randoms = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/{sample}.{peak_cutoff}.fa",
    run:
        #os.makedirs(os.path.dirname(output.peak_fastas), exist_ok=True)
        # -s: force strandedness.
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.bed} -fo {output.peak_fastas} -s")
        
        if os.path.exists(output.peak_fastas):
            with open(output.peak_fastas) as f:
                seqs = [x for x in f.readlines() if x[0]!='>']
                scripts.random_sequence.write_random_seqs(seqs, output.randoms)
                
rule write_fastas_merged:
    input:
        one = TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.1.bed',
        two = TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.2.bed', 
        three = TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.3.bed', 
    output:
        peak_fa_1 = TOP_DIR + "/peaks/clipper/filtered/fastas/merged.{protein}.1.fa",
        peak_fa_2 = TOP_DIR + "/peaks/clipper/filtered/fastas/merged.{protein}.2.fa", 
        peak_fa_3 = TOP_DIR + "/peaks/clipper/filtered/fastas/merged.{protein}.3.fa", 
        randoms_1 = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/merged.{protein}.1.fa",
        randoms_2 = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/merged.{protein}.2.fa", 
        randoms_3 = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/merged.{protein}.3.fa", 
    run:
        # -s: force strandedness.
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.one} -fo {output.peak_fa_1} -s")
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.two} -fo {output.peak_fa_2} -s")
        shell("bedtools getfasta -fi " + config['genomic_fasta'] + " -bed {input.three} -fo {output.peak_fa_3} -s")
        
        for peak_fasta, random_fa in [
            (output.peak_fa_1, output.randoms_1),
            (output.peak_fa_2, output.randoms_2),
            (output.peak_fa_3, output.randoms_3)]:
            
            if os.path.exists(peak_fasta):
                with open(peak_fasta) as f:
                    seqs = [x for x in f.readlines() if x[0]!='>']
                    scripts.random_sequence.write_random_seqs(seqs, random_fa)                
        
"""
rule call_clipper:
    input:
        bam = SAMS_DIR + '/genome_only/{sample}.bam',
        bed12 = 'data/processed/repeats_and_longest_txpt_per_gene.bed',
    output:
        bed = TOP_DIR + '/peaks/clipper/{sample}.bed'
    conda:
        'clipper/environment3.yml'
    threads: 8
    shell:
        CLIPPER_PATH + " -b {input.bam} -o {output.bed} -s GRCh38_v29"
"""

rule filter_clipper:
    input:
        bed = TOP_DIR + '/peaks/clipper/{sample}.bed'
    output:
        one = TOP_DIR + '/peaks/clipper/filtered/{sample}.1.bed',
        two = TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed',
        three = TOP_DIR + '/peaks/clipper/filtered/{sample}.3.bed',
    run:
        import scripts.filter_clipper
        scripts.filter_clipper.filter_clipper(str(input.bed))
        
        
rule remove_non_chromosomal_reads:
    input:
        bam = SAMS_DIR + '/dedup/{sample}.bam',
    output:
        bam = SAMS_DIR + '/genome_only/{sample}.bam',
    run:
        shell("samtools view -b {input.bam} chr{{1..22}} > {output.bam}")