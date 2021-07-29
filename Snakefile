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
configfile: "config.yaml"


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
        #expand(BIGWIG + "/{sample}.bigwig", sample=samples),
        #expand(SAMS_DIR + "/3end/{sample}.bam", sample=samples),
        config['counts'].rstrip('/') + "/bigwig_3prime_counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/bam_3prime_counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/featureCounts_on_bams.txt",
        expand(SAMS_DIR + '/genome_only/{sample}.bam', sample=samples),

        # Uncomment this to enable the use of clipper.
        #bed = expand(TOP_DIR + '/peaks/clipper/{sample}.bed', sample=samples),
        expand(TOP_DIR + "/outs/homer/{sample}.2/{sample}.2", sample=samples),
        expand(TOP_DIR + "/peaks/clipper/filtered/fastas/{sample}.2.fa", sample=samples), 
        expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed', sample=samples),
        expand(TOP_DIR + "/outs/homer/merged.{protein}.2/{protein}.2", protein=proteins),
    run:
        shell("echo Completed!")

include: 'rules/references.smk'
include: 'rules/read_preprocessing_and_mapping.smk'


#########################################################
# Run clipper.
#########################################################s
    
rule call_homer:
    input:
        peak_fasta = TOP_DIR + "/peaks/clipper/filtered/fastas/{sample}.2.fa",
        randoms = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/{sample}.2.fa",
    output:
        homer_results = TOP_DIR + "/outs/homer/{sample}.2/{sample}.2",
    log: TOP_DIR + "/outs/homer/logs/{sample}.2.log"
    conda:
        "envs/homer.yml"
    shell:
        "homer2 denovo -i {input.peak_fasta} -b {input.randoms} -len 6 -S 10 -strand + -o {output.homer_results} &> {log}"

rule call_homer_merged:
    input:
        peak_fasta = TOP_DIR + "/peaks/clipper/filtered/fastas/merged.{protein}.2.fa",
        randoms = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/merged.{protein}.2.fa",
    output:
        homer_results = TOP_DIR + "/outs/homer/merged.{protein}.2/{protein}.2",
    log:
        TOP_DIR + "/outs/homer/logs/{protein}.2.log"
    conda:
        "envs/homer.yml"
    shell:
        "homer2 denovo -i {input.peak_fasta} -b {input.randoms} -len 6 -S 10 -strand + -o {output.homer_results} &> {log}"
    
rule intersect_peaks:  
    input:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.1.bed', sample=samples),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed', sample=samples),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/{sample}.3.bed', sample=samples),
    output:
        one = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.1.bed', protein=proteins),
        two = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.2.bed', protein=proteins),
        three = expand(TOP_DIR + '/peaks/clipper/filtered/merged.{protein}.3.bed', protein=proteins),
    run:
        for protein in proteins:  # For each CLIP'd protein.
            for level in ['1', '2', '3']:  # For each level of stringency in peak calling.
                output_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.{level}.bed'
                _samples = [x for x in samples if x.split('_')[1]==protein]  # Replicate names.
                print('---')
                # Input replicate peak filenames.
                _files = [TOP_DIR + f'/peaks/clipper/filtered/{_sample}.{level}.bed' for _sample in _samples]
                print(protein, _samples)
                
                if len(_files) > 1:  # If multiple replicates for this protein...
                    thresh = math.ceil(0.75 * len(_files))  # Minimum replicate number.
                    filestring = " ".join(_files)
                    output_file = TOP_DIR + f'/peaks/clipper/filtered/merged.{protein}.{level}.bed'
                    
                    # {{ is converted to { by snakemake (avoids snakemake interpreting {} itself).
                    cmdstring = f"bedtools multiinter -i {filestring} | awk \'{{{{ if($4 >= {thresh}) {{{{ print }}}} }}}}\' | bedtools merge > {output_file}"
                    shell(cmdstring)  # Produce merged file.
                    
                else:  # If only one replicate, just copy that as the merged file.
                    shell(f"cp {_files[0]} {output_file}")

rule write_fastas:
    input:
        bed = TOP_DIR + '/peaks/clipper/filtered/{sample}.2.bed'
    output:
        peak_fastas = TOP_DIR + "/peaks/clipper/filtered/fastas/{sample}.2.fa", 
        randoms = TOP_DIR + "/peaks/clipper/filtered/fastas/randomControls/{sample}.2.fa",
    run:
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
        
#########################################################
# Read to gene assignment.
#########################################################
        
rule reads_per_gene_using_gffutils_bw:
    input:
        gtf = "data/processed/repeats_and_longest_txpt_per_gene.db",
        plus = expand(BIGWIG + "/3prime/{sample}.+.bigwig", sample=samples),
        minus = expand(BIGWIG + "/3prime/{sample}.-.bigwig", sample=samples)        
    output:
        config['counts'].rstrip('/') + "/bigwig_3prime_counts_transcripts.txt",
    run:
        input_dir = BIGWIG + "/3prime/"     
        shell("python scripts/reads_to_genes.py {input_dir} {input.gtf} {output}")
        
rule reads_per_gene_using_gffutils_bam:
    input:
        gtf = "data/processed/repeats_and_longest_txpt_per_gene.db",
        bam = expand(SAMS_DIR + "/3end/{sample}.bam", sample=samples)
    output:
        config['counts'].rstrip('/') + "/bam_3prime_counts_transcripts.txt",
    run:
        input_dir = SAMS_DIR + "/3end/"    
        shell("python scripts/reads_to_genes.py {input_dir} {input.gtf} {output}")
        
rule featureCounts_on_bams:
    input:
        bams = expand(SAMS_DIR + "/genome_only/{sample}.bam", sample=samples),        
        gtf = "assets/reference/featureCounts_formatted.gtf",
    output:
        counts = config['counts'].rstrip('/') + "/featureCounts_on_bams.txt",
    run:
        shell("featureCounts -t exon -g gene_name -F SAF -a {input.gtf} -o {output.counts} {input.bams}") 
        
rule create_featureCounts_formatted_gtf_from_regular_gtf:
    input:
        "assets/reference/only_tsl_1_and_NA.gtf"
    output:
        "assets/reference/featureCounts_formatted.gtf"
    shell:
        "python scripts/gtf_to_featureCounts_formatted_gtf.py {input} {output}"

rule subset_gtf_to_only_tsl1_and_NA:      
    input:
        config['feature_gtf']
    output:
        "assets/reference/only_tsl_1_and_NA.gtf"
    shell:
        "python scripts/subset_gtf_to_only_tsl1_and_NA.py {input} {output}"

#########################################################
# Old processing pipeline. Not run in this version.
#########################################################  
rule make_RNAs_object:
    input:
        gtf = "repeats_and_longest_txpt_per_gene.gff",
    output:
        TOP_DIR.rstrip('/') + '/data/rna.data'
    run:
        # Make data/rna.data file to assign reads to genes.
        import scripts.rnaDataFileMaker

        maker = scripts.rnaDataFileMaker.rnaDataFileMaker()
        RNAs = maker.make_from_gtf_file(gtf_filename=str(input.gtf))
        # -> outputs data/rna.data. 
        # This holds gtf object information - just regions, really.

        # Make counts file.
        # Outputs a data/bed_x.data file that holds signal, w/o RNA information:
        ex.make_signal_data_file()
        # Assigns signal to genes and writes counts.txt:
        ex.make_scheme_signal_RNA_data_files(rna_data_object=RNAs)

        # Make ann_counts.txt file. This has simplified
        # column names and biotypes.
        ex.annotate_counts_file()
    
rule reads_per_gene_statistics_vs_controls:
    run:
        import scripts.negativeCounts
        import scripts.positiveCounts
        import scripts.statsForCountsNB
        # If never run before:
        negatives = scripts.negativeCounts.negativeCounts(
            negative_metadata)#, xl_rate_fname='/Users/dp/pma/percentCrosslinked.xlsx')

        # Optional: write_txt=True to write some txt's of the data.
        negatives.save(write_object=True, write_txt=True)

        # If loading:
        negatives = scripts.negativeCounts.negativeCounts.load()

        # If never run before:
        positives = scripts.positiveCounts.positiveCounts(
            positive_metadata)#, xl_rate_fname='/Users/dp/pma/percentCrosslinked.xlsx')
        positives.save(write_object=True, write_txt=True)

        # If loading:
        positives = scripts.positiveCounts.positiveCounts.load()

        #nb = scripts.statsForCountsNB.statsForCountsNB.load()
        nb = scripts.statsForCountsNB.statsForCountsNB(
            negatives=negatives, positives=positives, data_dir=positive_metadata.top_dir + '/data/')
        nb.calculate_pvalues(which='per_read')
        nb.write_targets(which='per_read', outfname='default')
        nb.calculate_pvalues(which='per_protein')
        nb.write_targets(which='per_protein', outfname='default')
        nb.save()
        
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
"""        
rule bamToBigwig:
    input:
        bams = SAMS_DIR + "/dedup/{sample}.bam"
    output:
        bigwig = BIGWIG + "/{sample}.bigwig"
    shell:
        #"bamCoverage --binSize 10 -b {input} -o {output} -of bigwig"
        "bamCoverage --binSize 1 -b {input} -o {output} -of bigwig"
        """