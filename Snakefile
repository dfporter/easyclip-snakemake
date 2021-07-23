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

# Constrain sample wildcards to not contain '/'.
wildcard_constraints:
    sample="[^\/]+",

# Read the sample sheet into the exp object.
ex.read_scheme(config['samples'])

# Define these global variables used numerous times in the workflow.
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
        expand(BIGWIG + "/{sample}.bigwig", sample=samples),
        expand(SAMS_DIR + "/3end/{sample}.bam", sample=samples),
        "done/read_scheme.done",
        config['counts'].rstrip('/') + "/bigwig_3prime_counts_transcripts.txt",
        config['counts'].rstrip('/') + "/counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/featureCounts_on_bams.txt",
        chr_only = expand(SAMS_DIR + '/genome_only/{sample}.bam', sample=samples),

        # Uncomment this to enable the use of clipper.
        #bed = expand(config['top_dir'] + '/peaks/clipper/{sample}.bed', sample=samples)
    run:
        shell("echo Completed!")

include: 'rules/references.smk'
include: 'rules/read_preprocessing_and_mapping.smk'
    
#########################################################
# Run clipper.
#########################################################

rule call_clipper:
    input:
        bam = SAMS_DIR + '/genome_only/{sample}.bam',
        bed12 = 'data/processed/repeats_and_longest_txpt_per_gene.bed',
    output:
        bed = config['top_dir'] + '/peaks/clipper/{sample}.bed'
    conda:
        'clipper/environment3.yml'
    shell:
        CLIPPER_PATH + " -b {input.bam} -o {output.bed} -s GRCh38_v29"

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
        config['counts'].rstrip('/') + "/counts_transcripts.txt",
    run:
        input_dir = SAMS_DIR + "/3end/"    
        shell("python scripts/reads_to_genes.py {input_dir} {input.gtf} {output}")
        
rule featureCounts_on_bams:
    input:
        bams = expand(SAMS_DIR + "/dedup/{sample}.bam", sample=samples),        
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
        config['top_dir'].rstrip('/') + '/data/rna.data'
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
    input:
        SAMS_DIR + "/dedup/{sample}.bam"
    output:
        plus = ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.+.wig",  # Bedgraph. wig extension so IGV loads.
        minus = ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.-.wig",  # Bedgraph. wig extension so IGV loads.
    run:
        shell("bedtools genomecov -bg -strand + -3 -ibam {input} > {output.plus}")
        shell("bedSort {output.plus} {output.plus}")
        shell("bedtools genomecov -bg -strand - -3 -ibam {input} > {output.minus}")
        shell("bedSort {output.minus} {output.minus}")
        
rule convert_bam_to_3prime_end_only:
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
        
rule bamToBigwig:
    input:
        bams = SAMS_DIR + "/dedup/{sample}.bam"
    output:
        bigwig = BIGWIG + "/{sample}.bigwig"
    shell:
        #"bamCoverage --binSize 10 -b {input} -o {output} -of bigwig"
        "bamCoverage --binSize 1 -b {input} -o {output} -of bigwig"
        