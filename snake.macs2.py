
import os, sys, re, glob, pandas, importlib, shutil, math, pysam, pickle, collections, HTSeq, random, gffutils
import numpy as np
import scripts
import scripts.exp
import scripts.init_with_config
import scripts.rnaDataFileMaker
import scripts.make_repeats_chrom
import scripts.split_bam
import scripts.random_sequence

from scripts.macs2_snake_utils import load_narrowPeaks, peaks_as_bdg, load_fasta, write_random_seqs 
from scripts.macs2_snake_utils import find_chrom, add_seq, add_seq_transcriptome

# Load config file.
configfile: "config.yaml"

# Note: to get a graph of the workflow, run something like:
# snakemake -s snake.macs2.py --configfile configs/config.first_ddx21.yaml -j 1 --rulegraph -np | dot -Tpdf > dag.pdf

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

PROTEINS = list(set(df['Gene'])) 

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
MACS_PEAKS = config['top_dir'].rstrip('/') + '/outs/peaks/macs2'
BEDGRAPH = ex.file_paths['bedgraph'].rstrip('/') + '/merged'

#########################################################
# Custom merge settings.
#########################################################

# If run with: --config custom_merge="custom_merge.txt" 

if 'custom_merge' in config:
    custom_merge = pandas.read_csv(config['custom_merge'], sep='\t')
    
    PROTEINS = list(set(custom_merge['Name']))
    sample_to_set_of_custom_merges = dict(custom_merge.groupby("sam_path")['Name'].apply(set))

    #MERGED_BIGWIG = BIGWIG + '/merged'
    #BEDGRAPH = ex.file_paths['bedgraph'].rstrip('/') + '/merged'
    
    
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

#########################################################
# Begin rules.
#########################################################

rule all:
    input:
        "assets/reference/gencode.v39.transcripts.basic.fa",
        expand(MACS_PEAKS + "/{protein}_peaks.narrowPeak", protein=PROTEINS),
        expand(MACS_PEAKS + "/filtered/{protein}_peaks.narrowPeak", protein=PROTEINS),
        expand(MACS_PEAKS + "/bedgraph/{protein}.wig", protein=PROTEINS),
        expand(MACS_PEAKS + "/filtered/with_strand/no_ncrna/{protein}_peaks.narrowPeak", protein=PROTEINS),
        expand(MACS_PEAKS + "/fastas/genomic_no_ncrna/{protein}.fa", protein=PROTEINS),
        expand(MACS_PEAKS + "/fastas/randomControls/{protein}.fa", protein=PROTEINS),
        
        expand(MACS_PEAKS + "/homer/combined_short/{protein}/", protein=PROTEINS),
        expand(MACS_PEAKS + "/dreme/combined_short/{protein}/", protein=PROTEINS),
        expand(MACS_PEAKS + "/homer/nonexonic/{protein}/", protein=PROTEINS),
        expand(MACS_PEAKS + "/dreme/nonexonic/{protein}/", protein=PROTEINS),
        expand(MACS_PEAKS + "/homer/genomic_no_ncrna/{protein}/", protein=PROTEINS),
        expand(MACS_PEAKS + "/dreme/genomic_no_ncrna/{protein}/", protein=PROTEINS),
        
        #expand(BEDGRAPH + "/transcriptome/{protein}.wig", protein=PROTEINS),
        expand(SAMS_DIR + "/merged/no_repeats/{protein}.bam", protein=PROTEINS),
        expand(FASTQ_DIR + "/remapping_to_transcriptome/{protein}.fastq", protein=PROTEINS),
        expand(TOP_DIR + "/outs/salmon_quant/{protein}_quant", protein=PROTEINS),
        expand(SAMS_DIR + "/merged/transcriptome/{protein}.Aligned.out.sam", protein=PROTEINS),
        expand(MACS_PEAKS + "/transcriptome/filtered/{protein}_peaks.narrowPeak", protein=PROTEINS),
        expand(MACS_PEAKS + "/fastas/transcriptome/{protein}.fa", protein=PROTEINS),
        expand(MACS_PEAKS + "/filtered/nonexonic_stranded_no_ncrna/{protein}_peaks.narrowPeak", protein=PROTEINS),
        expand(MACS_PEAKS + "/fastas/combined/{protein}.fa", protein=PROTEINS),
        
        expand(MACS_PEAKS + "/filtered/nonexonic_stranded_no_ncrna/{protein}_peaks.narrowPeak", protein=PROTEINS),
        
        #expand(BEDGRAPH + "/transcriptome/{protein}.wig", protein=PROTEINS),
        
        # Peak bedgraphs.
        expand(BEDGRAPH + "/peaks/filtered/no_ncrna/{protein}.wig", protein=PROTEINS),
        expand(BEDGRAPH + "/peaks/filtered/nonexonic_stranded_no_ncrna/{protein}.wig", protein=PROTEINS),
        #expand(BEDGRAPH + "/peaks/unfiltered_genome/{protein}.wig", protein=PROTEINS),
        #expand(BEDGRAPH + "/transcriptome_peaks/unfiltered/{protein}.wig", protein=PROTEINS),
        expand(BEDGRAPH + "/transcriptome_peaks/filtered/{protein}.wig", protein=PROTEINS),
    run:
        shell("echo Completed!")

#include: 'rules/macs2.smk'


"""
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz
gunzip gencode.v39.transcripts.fa.gz
conda activate salmon
salmon index -t gencode.v39.transcripts.fa -i gencode_v39_salmon -k 31 --gencode

for fq in raw_data/*/*_1.fq.gz;
do

echo ~~~
echo ${fq}
base=$(basename ${fq} _1.fq.gz)
echo $base

echo "Processing sample ${base}"
salmon quant -i gencode_v39_salmon -l A \
         -1 fastq/${base}/${base}_1.fq.gz \
         -2 fastq/${base}/${base}_2.fq.gz \
         -p 8 --validateMappings -o quants/${base}_quant
done

"""
#########################################################
# Create transcriptome reference.
#########################################################

rule prepare_genePred_11_col_format:
    output:
        reformatted_genePred = 'assets/reference/gencode.v39.basic.annotation.11_col_format.genePred'
    run:
        shell('wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.basic.annotation.gff3.gz')
        shell('gunzip gencode.v39.basic.annotation.gff3.gz')
        shell('mv gencode.v39.basic.annotation.gff3 assets/reference/')
        shell('gff3ToGenePred assets/reference/gencode.v39.basic.annotation.gff3 assets/reference/gencode.v39.basic.annotation.genePred')
        

        def reformatted_genePred(fname):
            df = pandas.read_csv(fname, sep='\t', header=None)
            col_order = ['geneName', 'transcript_name', 'chrom',
                        'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd',
                        'exonCount', 'exonStarts', 'exonEnds']
            df['geneName'] = df[11] #[x.split(':')[1] for x in df[11]]
            df['transcript_name'] = df[0] #[x.split(':')[1] for x in df[0]]
            df['chrom'] = df[1]
            df['strand'] = df[2]
            df['txStart'] = df[3]
            df['txEnd'] = df[4]
            df['cdsStart'] = df[5]
            df['cdsEnd'] = df[6]
            df['exonCount'] = df[7]
            df['exonStarts'] = df[8]
            df['exonEnds'] = df[9]
            df = df[col_order]
            return df
        
        df = reformatted_genePred('assets/reference/gencode.v39.basic.annotation.genePred')

        df.to_csv(str(output.reformatted_genePred), sep='\t', header=None, index=False)

rule prepare_transcriptome_fasta:
    output:
        fasta = "assets/reference/gencode.v39.transcripts.basic.fa",
        gtf = "assets/reference/gencode.v39.basic.annotation.gtf",
    run:
        shell("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.basic.annotation.gtf.gz")
        shell("gunzip gencode.v39.basic.annotation.gtf.gz")
        shell("mv gencode.v39.basic.annotation.gtf assets/reference/")
        shell("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz")
        shell("gunzip gencode.v39.transcripts.fa.gz")
        shell("mv gencode.v39.transcripts.fa assets/reference/")
        
        basic_enst = set()
        with open("assets/reference/gencode.v39.basic.annotation.gtf") as f:
            for li in f:
                if li[0] == '#':
                    continue
                if 'tag "basic"' in li:
                    transcript_id = re.search('transcript_id "([^"]+)"', li)
                    basic_enst.add(transcript_id.group(1))
        outf = open(f"{output.fasta}", 'w')
        
        print(f"Saving only the {len(basic_enst)} ENSTs that are tagged as basic.")
        print(f"Writing to {output.fasta}...")
        wrote_enst = set()
        with open("assets/reference/gencode.v39.transcripts.fa") as f:
            writing_this_enst = False
            for li in f:
                if li[0] == '>':
                    enst = li.split('>')[1].split('|')[0]
                    if enst in basic_enst:
                        wrote_enst.add(enst)
                        writing_this_enst = True
                    else:
                        writing_this_enst = False
                    if writing_this_enst:
                        outf.write(li)
                elif writing_this_enst:
                    outf.write(li)
        outf.close()
        print(f"Wrote {len(wrote_enst)} ENST to fasta.")
        
rule star_transcriptome_index:
    input:
        fa = "assets/reference/gencode.v39.transcripts.basic.fa",
        #_dir = "assets/rDNA_star_index/stops_error"
    output:
        star_index = "assets/transcriptome_gencode39_star/Genome",
    threads:
        8
    shell:
        #shell("mkdir -p assets/transcriptome_gencode39_star")
        config['STAR'] + " --limitGenomeGenerateRAM 100000000000 --runThreadN 8 --runMode genomeGenerate --genomeDir assets/transcriptome_gencode39_star --genomeFastaFiles {input.fa}"

#########################################################
# Produce transcriptome bams.
#########################################################
rule get_genome_only:
    input:
        bam = SAMS_DIR + "/merged/{protein}.bam",
    output:
        bam = SAMS_DIR + "/merged/no_repeats/{protein}.bam",
    run:
        shell("samtools view -b {input.bam} chr{{1..22}} chrX chrY > {output.bam}")
        shell("samtools index {output.bam}")
        
rule fastq_of_genome_bams:
    input:
        bam = SAMS_DIR + "/merged/no_repeats/{protein}.bam",
    output:
        fq = FASTQ_DIR + "/remapping_to_transcriptome/{protein}.fastq",
    run:
        shell("samtools fastq {input.bam} > {output.fq}")
        
rule star_map_to_transcriptome:
    input:
        star_index = "assets/transcriptome_gencode39_star/Genome",
        fq = FASTQ_DIR + "/remapping_to_transcriptome/{protein}.fastq",
    output:
        sam = SAMS_DIR + "/merged/transcriptome/{protein}.Aligned.out.sam",
    threads:
        8
    run:
        cmd = config['STAR']
        cmd += " --genomeDir assets/transcriptome_gencode39_star"
        cmd += ' --runThreadN ' + str(threads)
        cmd += ' --readFilesIn {input.fq}'
        cmd += ' --alignIntronMax 1'  # Max intron size = 1. Setting to 0 causes the default behavior.
        cmd += ' --alignEndsType EndToEnd'  # Not completely sure this is right.
        #cmd += ' --outReadsUnmapped Fastx'
        #if os.path.splitext(input_r1)[-1] == '.gz':
        #    cmd += ' --readFilesCommand zcat'
        cmd += f" --outFileNamePrefix {SAMS_DIR}/merged/transcriptome/{wildcards.protein}."
        print(cmd)
        shell(cmd)

#########################################################
# Transcriptome bams to filtered peaks.
#########################################################
rule transcriptome_sam_to_bam:
    input:
        sam = SAMS_DIR + "/merged/transcriptome/{protein}.Aligned.out.sam",
    output:
        bam = SAMS_DIR + "/merged/transcriptome/{protein}.bam",
    run:
        shell("samtools view -Sb {input.sam} > {output.bam}.unsorted.bam")
        shell("samtools sort -o {output.bam} {output.bam}.unsorted.bam")
        shell("rm {output.bam}.unsorted.bam")
        shell("samtools index {output.bam}")
        
rule macs2_peak_calling_transcriptome:
    input:
        bam = SAMS_DIR + "/merged/transcriptome/{protein}.bam",
    output:
        peaks = MACS_PEAKS + "/transcriptome/{protein}_peaks.narrowPeak",
        #peaks_dir = "outs/macs2_peaks"
    conda:
        "envs/macs2.yml"
    shell:
        # --buffer-size default 100000 (100,000)
        # Minimum memory requested  for reading an alignment file is about:
        # Number of CHROMOSOME * BUFFER_SIZE * 2 Bytes.: 
        #shell(f"macs2 callpeak -t {input}  -f BAM -n {sample} --outdir outs/macs2_peaks/ --nomodel -q 0.1 --bw 100;")
        "macs2 callpeak -t {input.bam} -f BAM -n {wildcards.protein} "#-c assets/inputs/Input_RBFOX2_HepG2_hg38.bam"
              " --outdir {MACS_PEAKS}/transcriptome/ --nomodel -q 0.01 --bw 100 --buffer-size 100;"    
            
rule filter_peaks_transcriptome:
    input:
        peaks = MACS_PEAKS + "/transcriptome/{protein}_peaks.narrowPeak",
    output:
        subset_peaks = MACS_PEAKS + "/transcriptome/filtered/{protein}_peaks.narrowPeak",
    run:
        msgs = [f"{wildcards.protein}: {input.peaks}"]
        df = load_narrowPeaks(input.peaks)
        
        # For DDX21_high_glucose:
        # FDR 10^-12 ->      GCNGCNGC (second best)
        # FDR 10^-15 -> best=GCNGCNGC
        # FDR 10^-20 -> best=GCGGCGGC
        
        low_q = len(df[df['q']>=15].index)
        if low_q > 1000:
            df = df[df['q']>=15]
        else:
            df = df.sort_values(by='q', ascending=False).head(  min([1000, len(df.index)])  )
            
            if len(df.index) == 0:
                msg = f"{wildcards.protein}: No significant peaks at all."
            else:
                msg = f"{wildcards.protein}: Only {low_q} peaks with FDR<10^-15, using the best FDR" + \
                      f" {len(df.index)} peaks instead. Lowest q={min(df['q'])}."
            print(msg); msgs.append(msg)
        
        print(df)
        df.to_csv(output.subset_peaks, sep='\t', index=False)
            


#########################################################
# Genomic bams to filtered peaks.
#########################################################

rule macs2_peak_calling:
    input:
        bam = SAMS_DIR + "/merged/{protein}.bam"
    output:
        peaks = MACS_PEAKS + "/{protein}_peaks.narrowPeak",
        #peaks_dir = "outs/macs2_peaks"
    conda:
        "envs/macs2.yml"
    shell:
        #shell(f"macs2 callpeak -t {input}  -f BAM -n {sample} --outdir outs/macs2_peaks/ --nomodel -q 0.1 --bw 100;")
        "macs2 callpeak -t {input.bam} -f BAM -n {wildcards.protein} "#-c assets/inputs/Input_RBFOX2_HepG2_hg38.bam"
              " --outdir {MACS_PEAKS} --nomodel -q 0.01 --bw 100;"

rule peak_bedgraphs:
    input:
        peaks = MACS_PEAKS + "/{protein}_peaks.narrowPeak",
    output:
        bdg = MACS_PEAKS + "/bedgraph/{protein}.wig",
    run:
        peaks_as_bdg(input.peaks, output.bdg)

rule filter_peaks:
    input:
        peaks = MACS_PEAKS + "/{protein}_peaks.narrowPeak",
    output:
        subset_peaks = MACS_PEAKS + "/filtered/{protein}_peaks.narrowPeak",
    run:
        msgs = [f"{wildcards.protein}: {input.peaks}"]
        df = load_narrowPeaks(input.peaks)
        df = df.loc[[len(x)<6 for x in df.chrom], :]  # Lazy removal of weird chromosomes.
        
        # For DDX21_high_glucose:
        # FDR 10^-12 ->      GCNGCNGC (second best)
        # FDR 10^-15 -> best=GCNGCNGC
        # FDR 10^-20 -> best=GCGGCGGC
        
        low_q = len(df[df['q']>=15].index)
        if low_q > 1000:
            df = df[df['q']>=15]
        elif df.shape[0] > 0:
            df = df.sort_values(by='q', ascending=False).head(  min([1000, len(df.index)])  )
            
            msg = f"{wildcards.protein}: Only {low_q} peaks with FDR<10^-15, using the best FDR" + \
                  f" {len(df.index)} peaks instead. Lowest q={min(df['q'])}."
            print(msg); msgs.append(msg)
        else:
            print(f"{wildcards.protein}: No peaks.")
            
        print(df)
        df.to_csv(output.subset_peaks, sep='\t', index=False)

#########################################################
# Genomic filtered peaks to stranded peaks.
#########################################################
rule determine_strand:
    input:
        peaks = MACS_PEAKS + "/filtered/{protein}_peaks.narrowPeak",
        bam = SAMS_DIR + "/merged/{protein}.bam",
    output:
        peaks = MACS_PEAKS + "/filtered/with_strand/{protein}_peaks.narrowPeak",
    run:
        import pysam
        df = pandas.read_csv(input.peaks, sep='\t', index_col=False)
        can_determine, too_few_reads, too_even = 0, 0, 0
        calls = []
        bam = pysam.AlignmentFile(input.bam, 'rb')
        
        for chrom, start, end in zip(df.chrom, df.start, df.end):
            forward, reverse = 0, 0
            for n, r in enumerate(
                bam.fetch(contig=chrom, start=start, stop=end)):
                if r.is_reverse:
                    reverse += 1
                else:
                    forward += 1
                if n>1000:
                    break
            total = forward + reverse
            if total < 5:
                too_few_reads += 1
                calls.append('.')
            else:
                if forward >= 5 * reverse:
                    can_determine += 1
                    calls.append('+')
                elif reverse >= 5 * forward:
                    can_determine += 1
                    calls.append('-')
                else:
                    too_even += 1
                    calls.append('.')
        df['strand'] = calls
        print(f"Could determine strand: {can_determine}, <5 reads: {too_few_reads}, too even: {too_even}")
        df.to_csv(output.peaks, sep='\t', index=False)

################################################################
# Remove exonic peaks from the chromosomal coordinate peak file.
################################################################

def calc_degree_of_overlap(start, end, exons):
    """Replace with pybedtools or some faster method."""
    peak_len = end - start
    overlap = np.zeros(peak_len)
    for exon in exons:
        for n, pos in enumerate(range(start, end)):
            if exon.start <= pos <= exon.end:
                overlap[n] = 1
    frac_overlap = np.mean(overlap)
    return frac_overlap

def frac_overlap_exon(seqid, start, end, strand, db):
    """Cutoff for overlap is 70% of peak."""
    a = list(db.features_of_type(['exon'], limit=('chr' + seqid, start, end), strand=strand))
    if len(a):
        frac_overlap = calc_degree_of_overlap(start, end, a)
        if frac_overlap >= 0.40:
            return True
        return False
    a = list(db.features_of_type(['exon'], limit=(seqid, start, end), strand=strand))
    if len(a):
        frac_overlap = calc_degree_of_overlap(start, end, a)
        if frac_overlap >= 0.40:
            return True
        return False
    return False    

rule remove_exonic_peaks:
    input:
        peaks = MACS_PEAKS + "/filtered/with_strand/{protein}_peaks.narrowPeak",
    output:
        peaks = MACS_PEAKS + "/filtered/nonexonic_stranded/{protein}_peaks.narrowPeak",
    run:
        print(f"Removing exonic peaks from {input.peaks} and saving to {output.peaks}.")
        # Note that this db object needs to be changed to correspond with the transcriptome.
        # It should be generated from the gencode v39 gtf and saved under TOP_DIR.
        db_fname = "data/processed/features.db"
        db = gffutils.FeatureDB(db_fname, keep_order=True)
        df = pandas.read_csv(str(input.peaks), sep='\t', index_col=False)
        print(len(df.index))
        df['exonic'] = [frac_overlap_exon(*tup, db) for tup in zip(df.chrom, df.start, df.end, df.strand)]
        print(df)
        print(df[df['exonic']])
        non_exonic = df[~df['exonic']] 
        non_exonic.to_csv(str(output.peaks), sep='\t', index=False)

def overlap_ncrna(seqid, start, end, strand, db):
    exons = list(db.features_of_type(['exon'], limit=('chr' + seqid, start, end), strand=strand))
    if len(exons):
        for exon in exons:
            if exon.attributes['gene_type'] != 'protein_coding':
                return True
    return False

def assign_to_gene(seqid, start, end, strand, db):
    # Not used in this workflow currently.
    exons = list(db.features_of_type(['exon'], limit=('chr' + seqid, start, end), strand=strand))
    names = set()
    if len(exons):
        for exon in exons:
            names.add(exon.attributes['gene_name'])
    return names

def peak_loader(fname):
    with open(fname) as f:
        for li in f:
            if 'chrom' in li and ('start' in li) and ('end' in li):
                header = 0
            else:
                header = None
            break
    try:
        df = pandas.read_csv(fname, sep='\t', index_col=False, header=header)
        empty = False
    except:
        df = pandas.DataFrame([{
            'chrom': '1', 'start': 10, 'end': 110, 'strand': '+', 'name': 'no_peaks', 'q': 0}])
        df['fasta_seq_title'] = [f"{chrom}_{start}_{end}" for chrom, start, end in zip(df.chrom, df.start, df.end)]
        empty = True
    if (not empty) and (header is None):
            df = df.rename(
                columns={0: 'chrom', 1: 'start', 2: 'end', 3: 'name', 4: 'score',
                        5: 'strand', 6: 'signalValue', 6: 'p', 7: 'q', 8: 'point_source'})
    return df
    
rule remove_overlap_ncrna_from_non_exonic:
    # In progress.
    input:
        peaks = MACS_PEAKS + "/filtered/nonexonic_stranded/{protein}_peaks.narrowPeak",
        gtf = "assets/reference/gencode.v39.basic.annotation.gtf",
    output:
        peaks = MACS_PEAKS + "/filtered/nonexonic_stranded_no_ncrna/{protein}_peaks.narrowPeak",
    run:
        db_fname = "data/processed/features.db"
        db = gffutils.FeatureDB(db_fname, keep_order=True)
        df = peak_loader(str(input.peaks))
        df['overlaps_ncrna'] = [overlap_ncrna(*tup, db) for tup in zip(df.chrom, df.start, df.end, df.strand)]
        non_exonic = df[~df['overlaps_ncrna']] 
        non_exonic.to_csv(str(output.peaks), sep='\t', index=False)    

rule remove_overlap_ncrna_from_all_genomic:
    # In progress.
    input:
        peaks = MACS_PEAKS + "/filtered/with_strand/{protein}_peaks.narrowPeak",
        gtf = "assets/reference/gencode.v39.basic.annotation.gtf",
    output:
        peaks = MACS_PEAKS + "/filtered/with_strand/no_ncrna/{protein}_peaks.narrowPeak",
    run:
        db_fname = "data/processed/features.db"
        db = gffutils.FeatureDB(db_fname, keep_order=True)
        df = peak_loader(str(input.peaks))
        df['overlaps_ncrna'] = [overlap_ncrna(*tup, db) for tup in zip(df.chrom, df.start, df.end, df.strand)]
        non_exonic = df[~df['overlaps_ncrna']] 
        non_exonic.to_csv(str(output.peaks), sep='\t', index=False)    
        
################################################################
# Convert transcriptomic coordinates to genomic coordinates.
# Also, write their bedgraphs.
################################################################
def genomic_coord(iv, c):
    txpt_id = iv[0]
    gene = c.txpt_to_gene[txpt_id]
    chrom = c.transcript_dict[gene][txpt_id].chrom
    strand = c.transcript_dict[gene][txpt_id].strand
    a = c.transcript_dict[gene][txpt_id].genomic_position(iv[1])
    b = c.transcript_dict[gene][txpt_id].genomic_position(iv[2])
    return [chrom, min([a, b]), max([a, b]), strand]

def find_enst(_str):
    if _str.startswith('ENST'):
        return _str.split('|')[0]
    return re.search('(ENST[\d\.]+)', _str).group(1)

def transcriptome_peak_coordinates_to_genomic(fname):
    import scripts.wckdouglas_sequencing_tools_modified
    importlib.reload(scripts.wckdouglas_sequencing_tools_modified)
    from scripts.wckdouglas_sequencing_tools_modified import Transcript, Exon, converter

    c = converter()
    c.fname = "assets/reference/gencode.v39.basic.annotation.11_col_format.genePred"
    c.MakeTranscriptomeFromRefFlat()
    
    df = peak_loader(fname)
    as_genomic = df.copy()
    col_order = list(as_genomic.columns)
    coords = [
        genomic_coord([find_enst(txpt), start, end], c) for txpt,start,end in zip(df.chrom, df.start, df.end)]
    as_genomic['chrom'] = [x[0] for x in coords]
    as_genomic['start'] = [x[1] for x in coords]
    as_genomic['end'] = [x[2] for x in coords]
    as_genomic['strand'] = [x[3] for x in coords]
    as_genomic['IGV_location'] = [f"{iv[0]}:{iv[1]}-{iv[2]}" for iv in coords]
    as_genomic = as_genomic[col_order]
    return as_genomic

def peaks_to_bedgraph(fname, outfname):
    df = peak_loader(fname)
    ga = HTSeq.GenomicArray('auto', stranded=False)
    print(f"{fname} df -----> ", df)
    for chrom, start, end, q in zip(df.chrom, df.start, df.end, df.q):
        try:
            ga[HTSeq.GenomicInterval(chrom, start, end)] += q
        except:
            print(f"WARNING: couldn't make a genomic interval from: {(chrom, start, end)}")            
    ga.write_bedgraph_file(outfname)

rule genomic_coordinates_of_transcriptomic_peaks_filtered:
    input:
        peaks = MACS_PEAKS + "/transcriptome/filtered/{protein}_peaks.narrowPeak",
    output:
        peaks = MACS_PEAKS + "/transcriptome/filtered/genomic_coord/{protein}_peaks.narrowPeak",
        bdg = BEDGRAPH + "/transcriptome_peaks/filtered/{protein}.wig",
    run:
        df = transcriptome_peak_coordinates_to_genomic(str(input.peaks))
        df.to_csv(str(output.peaks), sep='\t', index=False)
        peaks_to_bedgraph(str(output.peaks), str(output.bdg))
        
#rule genomic_coordinates_of_transcriptomic_peaks_unfiltered:
#    input:
#        peaks = MACS_PEAKS + "/transcriptome/{protein}_peaks.narrowPeak",
#    output:
#        peaks = MACS_PEAKS + "/transcriptome/genomic_coord/{protein}_peaks.narrowPeak",
#        bdg = BEDGRAPH + "/transcriptome_peaks/unfiltered/{protein}.wig",
#    run:
#        df = transcriptome_peak_coordinates_to_genomic(str(input.peaks))
#        df.to_csv(str(output.peaks), sep='\t', index=False)  
#        peaks_to_bedgraph(str(output.peaks), str(output.bdg))
        
#########################################################
# Write bedgraphs for the genomic peaks.
#########################################################        
rule peak_bedgraph_no_ncrna:
    input:
        peaks =MACS_PEAKS + "/filtered/with_strand/no_ncrna/{protein}_peaks.narrowPeak"
    output:
        bdg = BEDGRAPH + "/peaks/filtered/no_ncrna/{protein}.wig"
    run:
        peaks_to_bedgraph(str(input.peaks), str(output.bdg))
        
rule peak_bedgraph_nonexonic_stranded_no_ncrna:
    input:
        peaks = MACS_PEAKS + "/filtered/nonexonic_stranded_no_ncrna/{protein}_peaks.narrowPeak",
    output:
        bdg = BEDGRAPH + "/peaks/filtered/nonexonic_stranded_no_ncrna/{protein}.wig", 
    run:
        peaks_to_bedgraph(str(input.peaks), str(output.bdg))
        
#rule peak_bedgraph_unfiltered_genome:
#    input:
#            peaks = MACS_PEAKS + "/{protein}_peaks.narrowPeak", 
#    output:
#        bdg = BEDGRAPH + "/peaks/unfiltered_genome/{protein}.wig", 
#    run:
#        peaks_to_bedgraph(str(input.peaks), str(output.bdg))
#########################################################
# Write fastas.
#########################################################


def write_fa(df, out_fa):
    if df.shape[0] > 0:
        with open(out_fa, 'w') as f:
            f.write(''.join(
                [f">{iv}\n{seq}\n" for iv, seq in zip(df.fasta_seq_title, df.seq)]
            ))
    else:
        with open(out_fa, 'w') as f:
            f.write('>No_peaks\nNNNNNNNNNNNN\n>No_peaks2\nNNNNNN\n')
    
rule write_fastas:
    input:
        peaks = expand(MACS_PEAKS + "/filtered/with_strand/no_ncrna/{protein}_peaks.narrowPeak", protein=PROTEINS),
    output:
        peak_fastas = expand(MACS_PEAKS + "/fastas/genomic_no_ncrna/{protein}.fa", protein=PROTEINS),
        randoms = expand(MACS_PEAKS + "/fastas/randomControls/{protein}.fa", protein=PROTEINS),
    run:
        genome = load_fasta(config['genomic_fasta'])
        peaks_files = [str(x) for x in input.peaks]
        fa_files = [str(x) for x in output.peak_fastas]
        random_files = [str(x) for x in output.randoms]
        
        for fname, out_fa, out_random in zip(peaks_files, fa_files, random_files):
            df = peak_loader(fname)
            df = df[df['strand']!='.']
            df = add_seq(df, genome)
            df = df.loc[[len(seq)>9 for seq in df.seq], :]

            write_fa(df, out_fa)
            write_random_seqs(df['seq'], out_random)

rule write_fastas_genomic_nonexonic:
    input:
        peaks = expand(MACS_PEAKS + "/filtered/nonexonic_stranded_no_ncrna/{protein}_peaks.narrowPeak", protein=PROTEINS),
    output:
        peak_fastas = expand(MACS_PEAKS + "/fastas/nonexonic_stranded_no_ncrna/{protein}.fa", protein=PROTEINS),
    run:
        genome = load_fasta(config['genomic_fasta'])
        peaks_files = [str(x) for x in input.peaks]
        fa_files = [str(x) for x in output.peak_fastas]
        
        for fname, out_fa in zip(peaks_files, fa_files):
            print(fname)
            df = peak_loader(fname)
            df = df[df['strand']!='.']
            df = add_seq(df, genome)
            df = df.loc[[len(seq)>9 for seq in df.seq], :]
            write_fa(df, out_fa)

rule write_combined_transcriptome_and_nonexonic_fasta:
    input:
        fa1 = expand(MACS_PEAKS + "/fastas/transcriptome/{protein}.fa", protein=PROTEINS),
        fa2 = expand(MACS_PEAKS + "/fastas/nonexonic_stranded_no_ncrna/{protein}.fa", protein=PROTEINS),
    output:
        fa = expand(MACS_PEAKS + "/fastas/combined/{protein}.fa", protein=PROTEINS),
    run:
        txpt = [str(x) for x in input.fa1]
        nonex = [str(x) for x in input.fa2]
        out = [str(x) for x in output.fa]
        for _fa1, _fa2, out_fa in zip(nonex, txpt, out):
            shell(f"cat {_fa1} {_fa2} > {out_fa}")

rule size_filter_fasta:
    input:
        fa = expand(MACS_PEAKS + "/fastas/combined/{protein}.fa", protein=PROTEINS),
    output:
        fa = expand(MACS_PEAKS + "/fastas/combined/short/{protein}.fa", protein=PROTEINS),
    run:
        in_fa = [str(x) for x in input.fa]
        out_fa = [str(x) for x in output.fa]
        for _fa1, _fa2 in zip(in_fa, out_fa):
            outf = open(_fa2, 'w')
            name = ''
            with open(_fa1) as f:
                for li in f:
                    if li[0] == '>':
                        name = li
                    else:  # If the sequence len is <500, write the peak
                        if len(li) < 500:
                            outf.write(f"{name}{li}")
            outf.close()
            
rule write_fastas_transcriptome:
    input:
        peaks = expand(MACS_PEAKS + "/transcriptome/filtered/{protein}_peaks.narrowPeak", protein=PROTEINS),
        fasta = "assets/reference/gencode.v39.transcripts.basic.fa",
    output:
        peak_fastas = expand(MACS_PEAKS + "/fastas/transcriptome/{protein}.fa", protein=PROTEINS),
        randoms = expand(MACS_PEAKS + "/fastas/transcriptome/randomControls/{protein}.fa", protein=PROTEINS),
    run:
        genome = load_fasta(str(input.fasta))
        peaks_files = [str(x) for x in input.peaks]
        fa_files = [str(x) for x in output.peak_fastas]
        random_files = [str(x) for x in output.randoms]
        
        for fname, out_fa, out_random in zip(peaks_files, fa_files, random_files):
            print(fname)
            print(f"Writing to {out_fa} and {out_random}.")
            
            df = peak_loader(fname)
            df = add_seq_transcriptome(df, genome)
            df['fasta_seq_title'] = [x.split('|')[0] + x.split('|')[-1] for x in df.fasta_seq_title]
            #df['fasta_seq_title'] = [re.sub('.', '_', x) for x in df.fasta_seq_title]
            df = df.loc[[len(seq)>9 for seq in df.seq], :]
            
            print(f"len df.index = {len(df.index)}.")
            write_fa(df, out_fa)
            write_random_seqs(df['seq'], out_random)


        
#########################################################
# Motif calling.
#########################################################

#########################################################
# Motif calling - combined, short.
rule call_homer_comb:
    input:
        peak_fasta = MACS_PEAKS + "/fastas/combined/short/{protein}.fa",
        randoms = MACS_PEAKS + "/fastas/randomControls/{protein}.fa",
    output:
        homer_results = directory(MACS_PEAKS +  "/homer/combined_short/{protein}/"),
    conda:
        "envs/homer.yml"
    shell:
        "findMotifs.pl {input.peak_fasta} fasta {output.homer_results} -fasta {input.randoms} -rna -homer1 -len 6,7,8"

rule call_dreme_comb:
    input:
        fasta = MACS_PEAKS + "/fastas/combined/short/{protein}.fa",
    output:
        folder = directory(MACS_PEAKS +  "/dreme/combined_short/{protein}/"),
    conda:
        "envs/meme.yml"        
    shell:
        "dreme -p {input.fasta} -rna -norc -oc {output.folder} -e 0.000001 -t 300"
        
#########################################################
# Motif calling - genomic.
rule call_homer:
    input:
        peak_fasta = MACS_PEAKS + "/fastas/genomic_no_ncrna/{protein}.fa",
        randoms = MACS_PEAKS + "/fastas/randomControls/{protein}.fa",
    output:
        homer_results = directory(MACS_PEAKS +  "/homer/genomic_no_ncrna/{protein}/"),
    conda:
        "envs/homer.yml"
    shell:
        "findMotifs.pl {input.peak_fasta} fasta {output.homer_results} -fasta {input.randoms} -rna -homer1 -len 6,7,8"
        
rule call_dreme:
    input:
        fasta = MACS_PEAKS + "/fastas/genomic_no_ncrna/{protein}.fa",
    output:
        folder = directory(MACS_PEAKS +  "/dreme/genomic_no_ncrna/{protein}/"),
    conda:
        "envs/meme.yml"        
    shell:
        "dreme -p {input.fasta} -rna -norc -oc {output.folder} -e 0.000001 -t 300"
        
#########################################################
# Motif calling - nonexonic.
rule call_homer_nonexonic:
    input:
        peak_fasta = MACS_PEAKS + "/fastas/nonexonic_stranded_no_ncrna/{protein}.fa",
        randoms = MACS_PEAKS + "/fastas/randomControls/{protein}.fa",
    output:
        homer_results = directory(MACS_PEAKS +  "/homer/nonexonic_stranded_no_ncrna/{protein}/"),
    conda:
        "envs/homer.yml"
    shell:
        "findMotifs.pl {input.peak_fasta} fasta {output.homer_results} -fasta {input.randoms} -rna -homer1 -len 6,7,8"
        
rule call_dreme_nonexonic:
    input:
        fasta = MACS_PEAKS + "/fastas/nonexonic_stranded_no_ncrna/{protein}.fa",
    output:
        folder = directory(MACS_PEAKS +  "/dreme/nonexonic_stranded_no_ncrna/{protein}/"),
    conda:
        "envs/meme.yml"        
    shell:
        "dreme -p {input.fasta} -rna -norc -oc {output.folder} -e 0.000001 -t 300"
        
rule transcriptome_bedgraph:
    input:
        db_fname = "data/processed/features.db",
        bdg_plus = BEDGRAPH + "/{protein}.+.wig",
        bdg_minus = BEDGRAPH + "/{protein}.-.wig",
    output:
        bdg = BEDGRAPH + "/transcriptome/{protein}.wig",
    run:
        import pysam, gffutils, random
        import scripts.bedgraphs
        # The version in data/processed/features.db has introns added.
        db = gffutils.FeatureDB(str(input.db_fname), keep_order=True)
        
        # feats defines the order of the artificial chromosome.
        feats = list(set([x.id for x in db.features_of_type('transcript')]))

        # We'll assume it's ok to have everything on one giant chromosome?

        print(f"{len(feats)} id values.")

        # These are used for reverse mapping later back to chromosomal coordinates.
        #txptGa = HTSeq.GenomicArrayOfSets('auto', stranded=False)
        exonGa = HTSeq.GenomicArrayOfSets('auto', stranded=False)

        locs = {}
        txpt_locs = {}
        txCoordStart = 0
        chrm_len = 0
        spacer = 1000

        #txptome_ga = HTSeq.GenomicArray('auto', stranded=True)

        ga_input = {f"{wildcards.protein}": scripts.bedgraphs.read_bedgraphs(f"{input.bdg_plus}")}
        ga_outs = {f"{wildcards.protein}": HTSeq.GenomicArray('auto', stranded=False)}

        #for n, txpt in enumerate(random.sample(feats, k=1000)):
        for n, txpt in enumerate(feats):
            if not(n % 1000):
                print(f'{n}/{len(feats)}', txpt)

            # This next line takes a long time.
            exons = sorted(list(db.children(txpt, featuretype='exon')), key=lambda x: x.start)

            txpt_len = sum([ex.end - ex.start for ex in exons])

            ht_exons = [HTSeq.GenomicInterval(ex.chrom.split('chr')[-1], ex.start, ex.end, ex.strand) for ex in exons]

            txCoordStart += spacer
            txpt_locs[txpt] = (txCoordStart, txCoordStart+txpt_len)
            
            # Read coverage. 
            for protein, ga in ga_input.items():
                txCoordExonStart = txCoordStart
                for iv in ht_exons:  # In order.

                    #locs[(txpt, (iv.start, iv.end))] = (txCoordExonStart, txCoordExonStart + iv.end-iv.start)
                    genomeExonStart = iv.start

                    # This exonGa object is holding info to convert txptCoords -> genomicCoords.
                    exonGa[HTSeq.GenomicInterval('txpt', txCoordExonStart, txCoordExonStart + iv.end-iv.start)
                          ] = set([f'{txpt}//{txCoordExonStart}//{iv.chrom}:{iv.start}'])  

                    for _iv, value in ga[iv].steps():

                        if value > 0:
                            txCoordValStart = _iv.start - genomeExonStart + txCoordExonStart
                            txCoordValEnd = _iv.end - genomeExonStart + txCoordExonStart
                            
                            ga_outs[protein][HTSeq.GenomicInterval('txpt', txCoordValStart, txCoordValEnd)] += value

                    txCoordExonStart += iv.end-iv.start

            txCoordStart += txpt_len

        for k, fh in ga_outs.items():
            ga_outs[k].write_bedgraph_file(str(output.bdg))


rule salmon_quant:
    input:
        fq = FASTQ_DIR + "/remapping_to_transcriptome/{protein}.fastq",
    output:
        #bam = SAMS_DIR + "/merged/remapped_to_transcriptome/{protein}.bam",
        quant = directory(TOP_DIR + "/outs/salmon_quant/{protein}_quant/"),
    conda:
        "envs/salmon.yml"
    shell:
        "salmon quant -i ~/oak/seq/genomes/gencode_v39_salmon -l A -r {input.fq} -p 8" + \
        " -o {output.quant}"
"""
head -n 22  .snakemake/conda/472b3364dee5bb1b995b1b3e567d0f0e/share/homer/motifs/groups/all.motifs > ./all.motifs

G4 RNA:     5’-GUU GGG GCG GGC GUU GGG UUU GGG GGG ACG-3’; 
CLIP-RNA:   5’-CAG UGG CUG CUG CUG UGG CCA CGU G-3’; 
Random-RNA: 5’-UUG GUA GUC ACC CCA AAU UGU UAU U-3’. 

CLIP RNA: CAGTG GC TGC TGC TGTGGCCACGTG
G4 RNA:   GTTGGGGCGGGCGTTGGGTTTGGGGGGACG
"""
