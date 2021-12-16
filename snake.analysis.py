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
        "assets/reference/Dfam_curatedonly.embl",
        "assets/repeats_star_index/Genome",
        config['feature_gtf'],
        "assets/reference/featureCounts_formatted.gtf",
        "data/processed/features.db",
        "data/processed/repeats_and_longest_txpt_per_gene.bed",
        expand(BIGWIG + "/{sample}.bigwig", sample=samples),
        expand(BIGWIG + "/rpm/{sample}.bigwig", sample=samples),
        #expand(SAMS_DIR + "/3end/{sample}.bam", sample=samples),
        #config['counts'].rstrip('/') + "/bigwig_3prime_counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/bam_3prime_counts_transcripts.txt",
        #config['counts'].rstrip('/') + "/featureCounts_on_bams.txt",
        #expand(SAMS_DIR + '/genome_only/{sample}.bam', sample=samples),
        control_bigwigs(),
        control_bigwigs_3prime(),
        
        exon_ht = expand(config['counts'].rstrip('/') + '/htseq_count_raw.exon.{sample}.txt', sample=samples),
        intr_ht = expand(config['counts'].rstrip('/') + '/htseq_count_raw.intron.{sample}.txt', sample=samples),
        #htseq_counts = expand(config['counts'].rstrip('/') + '/htseq_count_raw.{sample}.txt', sample=samples),
        htseq_counts = config['counts'].rstrip('/') + '/htseq_count_raw.all.txt',
        htseq_txt = config['counts'].rstrip('/') + '/htseq_count_rpm.all.txt', 
        htseq_xlsx = config['counts'].rstrip('/') + '/htseq_count_rpm.all.xlsx', 
        
        bdg_plus = expand(ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.+.wig",  sample=samples),
        bdg_minus = expand(ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.-.wig", sample=samples),
        
        control_counts = "random_controls/counts/htseq_count_raw.all.txt",
        rna_data = TOP_DIR + '/data/rna.data',
        counts = TOP_DIR + '/outs/counts/ann_counts.bedgraphs.txt',
        reads_per_million = TOP_DIR + '/outs/counts/reads_per_million.bedgraphs.txt',
        pvals_per_read = TOP_DIR + "/tables/pvals_per_read.xlsx",
        
        per_read_xlsx = TOP_DIR + "/tables/RNA targets per read vs random non-RBPs.xlsx",
        #per_protein_xlsx = TOP_DIR + "/tables/RNA targets per protein vs random non-RBPs.xlsx",
    run:
        shell("echo Completed!")

include: 'rules/references.smk'

#########################################################
# Old processing pipeline. 
#########################################################  
rule make_RNAs_object:
    input:
        gtf = "data/processed/repeats_and_longest_txpt_per_gene.gff",

    output:
        data = TOP_DIR + '/data/rna.data',
    run:
        # Make data/rna.data file to assign reads to genes.
        import scripts.rnaDataFileMaker

        maker = scripts.rnaDataFileMaker.rnaDataFileMaker()
        RNAs = maker.make_from_gtf_file(
            gtf_filename=str(input.gtf), output_data_filename=str(output.data))
        # -> outputs data/rna.data. 
        # This holds gtf object information - just regions, really.

rule raw_counts_file_from_bedgraph:
    input:
        data = TOP_DIR + '/data/rna.data',
        # Bedgraphs. wig extension so IGV loads.
        plus = expand(ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.+.wig",  sample=samples),
        minus = expand(ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.-.wig", sample=samples),
    output:
        counts = TOP_DIR + '/outs/counts/counts.bedgraphs.txt'
    run:
        RNAs = pickle.load(open(str(input.data), 'rb'))

        # Outputs a data/bed_x.data file that holds signal, w/o RNA information:
        # Assigns signal to genes and writes raw counts.txt:
        ex.make_scheme_signal_RNA_data_files(
            rna_data_object=RNAs, no_clobber=True, counts_fname=str(output.counts))
        
rule total_read_numbers:
    input:
        plus = expand(ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.+.wig",  sample=samples),
        minus = expand(ex.file_paths['bedgraph'].rstrip('/') + "/{sample}.-.wig", sample=samples),
    output:
        total_reads = config['data']+"/total_read_numbers.txt",
    run:
        import scripts.total_read_numbers
        scripts.total_read_numbers.total_read_numbers(
            folder=ex.file_paths['bedgraph'], outfile=str(output.total_reads))
        
rule counts_files:
    input:
        counts = TOP_DIR + '/outs/counts/counts.bedgraphs.txt',
        gtf = "data/processed/repeats_and_longest_txpt_per_gene.gff",
        total_reads = config['data']+"/total_read_numbers.txt",
    output:
        counts = TOP_DIR + '/outs/counts/ann_counts.bedgraphs.txt',
        reads_per_million = TOP_DIR + '/outs/counts/reads_per_million.bedgraphs.txt',
    run:
        import scripts.readsPerGene
        from scripts.readsPerGene import rawReadsPerGene, readsPerMillion
        
        rpg = rawReadsPerGene(str(input.counts), scheme_filename=config['samples'])
        rpg.add_biotypes_column(gtf_filename=str(input.gtf))
        
        # Write the biotypes-added raw counts.
        rpg.df.to_csv(str(output.counts), sep='\t')

        rpm = readsPerMillion(rpg, load_total_read_numbers=str(input.total_reads))
        
        # Write the reads per million file.
        rpm.df.to_csv(str(output.reads_per_million), sep='\t')

rule combine_pvals_and_other_gene_information:
    input:
        pvals_per_read = TOP_DIR + "/tables/pvals_per_read.xlsx",
        counts = TOP_DIR + '/outs/counts/ann_counts.bedgraphs.txt',
        reads_per_million = TOP_DIR + '/outs/counts/reads_per_million.bedgraphs.txt',
        total_reads = config['data']+"/total_read_numbers.txt",
    output:
        per_read_xlsx = TOP_DIR + "/tables/RNA targets per read vs random non-RBPs.xlsx",
    run:
        import scripts.readsPerGene
        from scripts.readsPerGene import rawReadsPerGene, readsPerMillion
        raw_reads = rawReadsPerGene(str(input.counts), scheme_filename=config['samples'])
        per_million = readsPerMillion(raw_reads, load_total_read_numbers=str(input.total_reads))
        
        raw_means_df = raw_reads.average_by_protein()
        
        col_order = ["Gene/exon-or-intron", "Gene", "Exon or intron"]
        
        raw_means_df.columns = [x + " (# reads)" for x in raw_means_df.columns]
        raw_means_df['Gene'] = raw_means_df.index
        
        
        per_million_means_df = per_million.average_by_protein()        
        per_million_means_df.columns = [x + " (reads/million)" for x in per_million_means_df.columns]
        per_million_means_df['Gene'] = per_million_means_df.index
        
        pvals_df = pandas.read_excel(str(input.pvals_per_read), index_col=0)
        pvals_df['Gene'] = pvals_df.index
        pvals_df.columns = [x + " (FDR)" if x!='Gene' else x for x in pvals_df.columns]
        
        col_order += [x for x in pvals_df.columns if x!='Gene']
        col_order += [x for x in raw_means_df.columns if x!='Gene']
        col_order += [x for x in per_million_means_df.columns if x!='Gene']
        
        pvals_df = pvals_df.merge(raw_means_df, left_on='Gene', right_on='Gene', how='inner')
        print(pvals_df)
        print('-\n', per_million_means_df)
        pvals_df = pvals_df.merge(per_million_means_df, left_on='Gene', right_on='Gene', how='inner')                            
        
        
        pvals_df["Gene/exon-or-intron"] = pvals_df["Gene"]
        
        pvals_df["Gene"] = [x.split("::")[0] for x in pvals_df["Gene/exon-or-intron"]]
        pvals_df["Exon or intron"] = [x.split("::")[-1] for x in pvals_df["Gene/exon-or-intron"]]
        pvals_df = pvals_df.loc[:, col_order]
        print(pvals_df)
        
        pvals_df.to_excel(str(output.per_read_xlsx), index=False)
        
        
rule reads_per_gene_statistics_vs_controls:
    input:
        #counts = TOP_DIR + '/outs/counts/counts.bedgraphs.txt',
        counts = TOP_DIR + '/outs/counts/htseq_count_raw.all.txt',
        #random_counts = 'random_controls/counts/counts.bedgraphs.txt',
        random_counts = 'random_controls/counts/htseq_count_raw.all.txt',
        cntrl_total_reads = 'random_controls/data/total_read_numbers.txt',
        total_reads = TOP_DIR + '/data/total_read_numbers.txt',
    output:
        pvals_per_read = TOP_DIR + "/tables/pvals_per_read.xlsx",
    run:
        import scripts.negativeCounts
        import scripts.positiveCounts
        import scripts.statsForCountsNB
        
        from types import SimpleNamespace
        
        negative_metadata = SimpleNamespace(**{
            'scheme_file_with_random_proteins': 'random_controls/samples.txt',
            'top_dir': 'random_controls',
            'ann_counts_file': str(input.random_counts),
            'random_proteins': list(set(pandas.read_csv('random_controls/samples.txt', sep='\t')['Gene'])),
        })

        positive_metadata = SimpleNamespace(**{
            'scheme_file': ex.file_paths["samples"],
            'top_dir': TOP_DIR,
            'bed_file_dir': ex.file_paths['bedgraph'].rstrip('/'),
            'ann_counts_file': str(input.counts),
            'positive_proteins': proteins,
        })
 
        print(positive_metadata)
        print(negative_metadata)
        # If never run before:
        negatives = scripts.negativeCounts.negativeCounts(
            negative_metadata)#, xl_rate_fname='/Users/dp/pma/percentCrosslinked.xlsx')

        # Optional: write_txt=True to write some txt's of the data.
        negatives.save(write_object=True, write_txt=True)

        # If loading:
        #negatives = scripts.negativeCounts.negativeCounts.load()

        # If never run before:
        positives = scripts.positiveCounts.positiveCounts(positive_metadata)

        positives.save(write_object=True, write_txt=True)

        # If loading:
        #positives = scripts.positiveCounts.positiveCounts.load()

        #nb = scripts.statsForCountsNB.statsForCountsNB.load()
        nb = scripts.statsForCountsNB.statsForCountsNB(
            negatives=negatives, positives=positives, data_dir=TOP_DIR + '/data/')
        nb.calculate_pvalues(which='per_read')
        
        # Default: outfname = f'{self.positives.metadata.top_dir}/tables/pvals_{which}.xlsx'
        nb.write_pvals_single_file(which='per_read', outfname=str(output.pvals_per_read))
        #nb.write_targets(which='per_read', outfname=str(output.pvals_per_read))
        #nb.calculate_pvalues(which='per_protein')
        #nb.write_targets(which='per_protein', outfname='default')
        
        # Writes to f"{self.positives.metadata.data_folder}/stats.dill"
        nb.save(TOP_DIR + '/data/stats.dill')
        
        
        
#########################################################
# Filter peaks overlapping ncRNA.
#########################################################


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

        
rule remove_non_chromosomal_reads:
    input:
        bam = SAMS_DIR + '/dedup/{sample}.bam',
    output:
        bam = SAMS_DIR + '/genome_only/{sample}.bam',
    run:
        shell("samtools view -b {input.bam} chr{{1..22}} > {output.bam}")

rule biotype_lookup_file:
    input:
        "data/processed/U13369_repeats_and_genome.gtf",
    output:
        "data/processed/biotype_lookup.txt"
    run:
        df = pandas.read_csv(
            "data/processed/U13369_repeats_and_genome.gtf",
            sep='\t', comment='#', header=None)

        def parse(li):
            de = lambda x: x.group(1) if x is not None else 'Unknown'
            biotype = re.search('gene_type "([^"]+)"', li)
            gene_name = re.search('gene_name "([^"]+)"', li)
            gene_id = re.search('gene_id "([^"]+)"', li)
            return (de(gene_id), de(gene_name), de(biotype))

        lookup = [parse(x) for x in df[8]]

        def corrections(name, biotype):
            if 'tRNA' in name:
                return 'tRNA'
            if name=='LSU-rRNA_Hsa' or name=='SSU-rRNA_Hsa':
                return 'rRNA'
            if name in [f"U{n}" for n in range(1,18)]:
                return 'snRNA'
            if 'transcribed_spacer' in name:
                return 'pre-rRNA'
            if name == '7SK':
                return '7SK'
            if 'rRNA_Cel' in name:
                return 'Unknown'
            return biotype
        
        lookup = pandas.DataFrame(lookup, columns=['gene_id', 'gene_name', 'gene_type'])
        lookup = lookup.drop_duplicates()
        lookup['gene_type'] = [corrections(x,y) for x,y in zip(lookup.gene_name, lookup.gene_type)]
        lookup.to_csv("data/processed/biotype_lookup.txt", sep='\t', index=False)
        
#########################################################
# Read to gene assignment.
#########################################################
rule reads_per_gene_from_htseq_exonic:
    input:
        gtf = "data/processed/U13369_repeats_and_genome.gtf",
        bam = SAMS_DIR + "/dedup/{sample}.bam",
        #bams = expand(SAMS_DIR + "/dedup/{sample}.bam", sample=samples)
    output:
        #txt = config['counts'].rstrip('/') + '/htseq_count_raw.txt'
        txt = config['counts'].rstrip('/') + '/htseq_count_raw.exon.{sample}.txt',
        sam = TOP_DIR + '/data/processed/htseq.{sample}.nonexonic.sam',
    params:
        _tmp = TOP_DIR + "/data/processed/htseq.{sample}.sam",
    run:
        # -s yes: stranded. Uses gene_id by default (ENSG for us).
        shell("htseq-count -t exon --additional-attr=gene_name -s yes --format=bam" + \
        " -o {params._tmp} {input.bam} {input.gtf} > {output.txt}")
        shell("grep XF:Z:__no_feature {params._tmp} > {output.sam}")

rule reheader_htseq_nonexonic_sam:
    input:
        bam = SAMS_DIR + "/dedup/{sample}.bam",
        sam = TOP_DIR + '/data/processed/htseq.{sample}.nonexonic.sam',
    output:
        header = TOP_DIR + '/data/processed/bam_headers/{sample}.header.sam',
        sam = TOP_DIR + '/data/processed/htseq.{sample}.nonexonic.header.sam',
    run:
        shell("samtools view -H {input.bam} > {output.header}")
        shell("cat {output.header} {input.sam} > {output.sam}")

        
rule reads_per_gene_from_htseq_intronic:
    input:
        sam = TOP_DIR + '/data/processed/htseq.{sample}.nonexonic.header.sam',
        gtf = "data/processed/U13369_repeats_and_genome.gtf",
    output:
        txt = config['counts'].rstrip('/') + '/htseq_count_raw.intron.{sample}.txt',
    run:
        shell("htseq-count -t intron --additional-attr=gene_name -s yes --format=sam" + \
              " {input.sam} {input.gtf} > {output.txt}")
        
def list_duplicates(seq):
    tally = collections.defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>1)

rule reads_per_gene_from_htseq_combine_samples:
    input:
        exon = expand(config['counts'].rstrip('/') + '/htseq_count_raw.exon.{sample}.txt', sample=samples),
        intron = expand(config['counts'].rstrip('/') + '/htseq_count_raw.intron.{sample}.txt', sample=samples),
    output:
        txt = config['counts'].rstrip('/') + '/htseq_count_raw.all.txt', 
        xlsx = config['counts'].rstrip('/') + '/htseq_count_raw.all.xlsx',
    run:
        
        def merge_counts_files(file_list, left_strip='htseq_count_raw.exon.', right_strip='.txt'):
            dfs = [pandas.read_csv(
                f, index_col=[0], sep='\t', header=None, names=['gene_id', 'gene_name', 
                    os.path.basename(f).split(left_strip)[-1].split(right_strip)[0]]) \
                    for f in file_list]
        
            for df_ex in dfs:
                
                df_ex.gene_name = [
                        _name if (not pandas.isna(_name) and len(_name)) else _id \
                            for _id, _name in zip(df_ex.index, df_ex.gene_name)]
                
                # Check for duplicate gene ids. Cannot merge with non-unique index.
                for dup, locs in list_duplicates(df_ex.gene_name):
                    nameA, nameB = df_ex.iloc[locs[0]].name, df_ex.iloc[locs[0]].name
                    print(f"Duplicate names for {dup} at {locs}. Using gene_ids {nameA} and {nameB}")
                    df_ex.loc[nameA, 'gene_name'] = nameA
                    df_ex.loc[nameB, 'gene_name'] = nameB
                    
                dups = list_duplicates(df_ex.gene_name)
                df_ex.index = df_ex['gene_name']
                del df_ex['gene_name']
                            
            finaldf = pandas.concat(dfs, axis=1, join='inner').sort_index()

            return finaldf
        
        exons = merge_counts_files(input.exon, left_strip='htseq_count_raw.exon.', right_strip='.txt')
        introns = merge_counts_files(input.intron, left_strip='htseq_count_raw.intron.', right_strip='.txt')

        exons.index = [x + '::exon' for x in exons.index]
        introns.index = [x + '::intron' for x in introns.index]

        both = pandas.concat([exons, introns])

        both.to_csv(str(output.txt), sep='\t')
        both.to_excel(str(output.xlsx), engine='openpyxl')

rule reads_per_gene_from_htseq_normalize:
    input:
        total_reads = config['data'] + "/total_read_numbers.txt",
        txt = config['counts'].rstrip('/') + '/htseq_count_raw.all.txt', 
    output:
        txt = config['counts'].rstrip('/') + '/htseq_count_rpm.all.txt', 
        xlsx = config['counts'].rstrip('/') + '/htseq_count_rpm.all.xlsx', 
    run:
        df = pandas.read_csv(str(input.txt), sep='\t', index_col=0)
        per_mill = pandas.read_csv(str(input.total_reads), sep='\t', index_col=0)
        
        for col in df.columns:
            df[col] = df[col] * 1E6/per_mill.loc[col, 'Total read number']
            
        df.to_csv(str(output.txt), sep='\t')
        df.to_excel(str(output.xlsx), engine='openpyxl')
        
rule reads_per_gene_from_htseq_combine_exon_and_intron:
    input:
        exon = config['counts'].rstrip('/') + '/htseq_count_raw.exon.{sample}.txt',
        intron = config['counts'].rstrip('/') + '/htseq_count_raw.intron.{sample}.txt',
    output:
        txt = config['counts'].rstrip('/') + '/htseq_count_raw.{sample}.txt'
    run:
        exon = pandas.read_csv("{input.exon}", header=None, index_col=None, sep='\t')
        intron = pandas.read_csv("{input.intron}", header=None, index_col=None, sep='\t')
        
        for df in [exon, intron]:
            df.columns = ["gene_id", "gene_name", "{sample}"]
            df['gene_name'] = [_name if not pandas.isna(_name) else _id for _id, _name in zip(df.gene_id, df.gene_name)]
        exon['gene_name'] = [x + '::exon' for x in exon.gene_name]
        intron['gene_name'] = [x + '::intron' for x in intron.gene_name]
        both = pandas.concat([exon, intron], ignore_index=True)
        print(both)
        
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

rule download_control_data:
    output:
        #per_million = "random_controls/counts/reads_per_million.bedgraphs.txt",
        counts = "random_controls/counts/htseq_count_raw.all.txt",
        bigwigs = control_bigwigs(),
        bigwig_3prime = control_bigwigs_3prime(),
    run:

        os.makedirs("random_controls/bigwig/3prime", exist_ok=True)
        os.makedirs("random_controls/counts/", exist_ok=True)

        shell(f"wget http://khavarilab.stanford.edu/s/htseq_count_rawall.txt")
        shell("mv htseq_count_rawall.txt random_controls/counts/htseq_count_raw.all.txt")
        #shell("tar -xf random_non_RBP_control_counts.tar.gz -C random_controls/counts/")
        #shell("mv random_controls/counts/counts/* random_controls/counts/")
        #shell("rmdir random_controls/counts/counts/")
        
        shell("sleep 10")
        shell(f"wget http://khavarilab.stanford.edu/s/random_non_RBP_control_bigwig_cDNA_end_1tar.gz")
        shell("mv random_non_RBP_control_bigwig_cDNA_end_1tar.gz random_non_RBP_control_bigwig_cDNA_end_1.tar.gz")
        shell("tar -xf random_non_RBP_control_bigwig_cDNA_end_1.tar.gz -C random_controls/bigwig/3prime/")
        shell("mv random_controls/bigwig/3prime/zips1/* random_controls/bigwig/3prime/")
        shell("rm random_controls/bigwig/3prime/zips1/._*")
        shell("rmdir random_controls/bigwig/3prime/zips1/")
        
        shell("sleep 10")
        shell(f"wget http://khavarilab.stanford.edu/s/random_non_RBP_control_bigwig_cDNA_end_2tar.gz")
        shell("mv random_non_RBP_control_bigwig_cDNA_end_2tar.gz random_non_RBP_control_bigwig_cDNA_end_2.tar.gz")
        shell("tar -xf random_non_RBP_control_bigwig_cDNA_end_2.tar.gz -C random_controls/bigwig/3prime/")
        shell("mv random_controls/bigwig/3prime/zips2/* random_controls/bigwig/3prime/")
        shell("rm random_controls/bigwig/3prime/zips2/._*")
        shell("rmdir random_controls/bigwig/3prime/zips2/")
        
        shell("sleep 10")
        shell(f"wget http://khavarilab.stanford.edu/s/random_non_RBP_control_bigwigs_1tar.gz")
        shell("mv random_non_RBP_control_bigwigs_1tar.gz random_non_RBP_control_bigwigs_1.tar.gz")
        shell("tar -xf random_non_RBP_control_bigwigs_1.tar.gz -C random_controls/bigwig/")
        shell("mv random_controls/bigwig/zips1/* random_controls/bigwig/")
        shell("rm random_controls/bigwig/zips1/._*")
        shell("rmdir random_controls/bigwig/zips1/")
        
        shell("sleep 10")
        shell(f"wget http://khavarilab.stanford.edu/s/random_non_RBP_control_bigwigs_2tar.gz")
        shell("mv random_non_RBP_control_bigwigs_2tar.gz random_non_RBP_control_bigwigs_2.tar.gz")
        shell("tar -xf random_non_RBP_control_bigwigs_2.tar.gz -C random_controls/bigwig/")
        shell("mv random_controls/bigwig/zips2/* random_controls/bigwig/")
        shell("rm random_controls/bigwig/zips2/._*")
        shell("rmdir random_controls/bigwig/zips2/")

        shell("sleep 10")
        shell(f"wget http://khavarilab.stanford.edu/s/random_non_RBP_control_bigwigs_3tar.gz")
        shell("mv random_non_RBP_control_bigwigs_3tar.gz random_non_RBP_control_bigwigs_3.tar.gz")
        shell("tar -xf random_non_RBP_control_bigwigs_3.tar.gz -C random_controls/bigwig/")
        shell("mv random_controls/bigwig/zips3/* random_controls/bigwig/")
        shell("rm random_controls/bigwig/zips3/._*")
        shell("rmdir random_controls/bigwig/zips3/")   
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
        
def scale_factor(wildcards):
    df = pandas.read_csv(config['data']+"/total_read_numbers.txt", sep='\t')
    df["basename"] = [x.split('/')[-1].split('.')[0] for x in df.Dataset]
    return str(1000000/df.loc[df['basename']==wildcards.sample, 'Total read number'].iloc[0])
            
rule bigWigToScaledWig:
    input:
        bigwig = BIGWIG + "/{sample}.bigwig",
        total_reads = config['data']+"/total_read_numbers.txt",
    output:
        wig = TOP_DIR + "/scaledWig/{sample}.bigwig"
    params:
        scaling_factor = scale_factor
    conda:
        "envs/wiggletools.yml"
    shell:
        "wiggletools write {output.wig} scale {params.scaling_factor} {input.bigwig}"

rule scaledWigToBigWig:
    input:
        wig = TOP_DIR + "/scaledWig/{sample}.bigwig",
        chrom_sizes = "data/processed/chrom.sizes"
    output:
        bigwig = BIGWIG + "/rpm/{sample}.bigwig"
    shell:
        "wigToBigWig {input.wig} {input.chrom_sizes} {output.bigwig}"
        
rule bamToBigwig:
    input:
        bams = SAMS_DIR + "/dedup/{sample}.bam"
    output:
        bigwig = BIGWIG + "/{sample}.bigwig"
    shell:
        #"bamCoverage --binSize 10 -b {input} -o {output} -of bigwig"
        "bamCoverage --binSize 1 -b {input} -o {output} -of bigwig"
