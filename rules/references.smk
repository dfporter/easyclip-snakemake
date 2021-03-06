import os, sys, re, glob, pandas, importlib, shutil
import scripts
import scripts.build_genome_db
import scripts.make_repeats_chrom

rule download_genome_fasta:
    output:
        "assets/reference/GRCh38.primary_assembly.genome.fa"
    shell:
        "wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz;"
        "gunzip GRCh38.primary_assembly.genome.fa.gz;"
        "mv GRCh38.primary_assembly.genome.fa assets/reference/"
        
rule files_for_igv_genome:
    input:
        combined_gff = 'data/processed/repeats_and_longest_txpt_per_gene.gff',
        rDNA_gtf = 'data/processed/U13369_rRNA.gtf',
        repeats_fa = 'assets/reference/repeats_chrom.fa',
        rDNA_fa = "assets/reference/U13369_rRNA.fasta",
        genome_fa = config['genomic_fasta'],
    output:
        fa = "data/processed/U13369_repeats_and_genome.fa",
        gtf = "data/processed/U13369_repeats_and_genome.gtf",
    run:
        shell("cat {input.combined_gff} {input.rDNA_gtf} > {output.gtf}")
        shell("cat {input.rDNA_fa} {input.genome_fa} {input.repeats_fa}  > {output.fa}")
        
#########################################################
# References for mapped read->format conversion/gene assignment.
#########################################################   

rule combined_genome_and_repeats_bed12:
    input:
        db_longest = 'data/processed/repeats_and_longest_txpt_per_gene.db',
        db_all = "data/processed/features.db"
    output:
        bed_longest = "data/processed/repeats_and_longest_txpt_per_gene.bed",
        bed_all = "data/processed/features.bed"
    run:
        scripts.build_genome_db.export_bed12(input.db_longest, output.bed_longest)
        scripts.build_genome_db.export_bed12(input.db_all, output.bed_all)
        
rule combine_genome_and_repeats_gtfs:
    input:
        repeats = 'assets/reference/repeats.gtf',
        genome = 'data/processed/longest_txpt_per_gene.gff',
    output:
        combined = 'data/processed/repeats_and_longest_txpt_per_gene.gff',
        combined_db = 'data/processed/repeats_and_longest_txpt_per_gene.db',
    run:
        shell(f'cat {input.genome} {input.repeats} > {output.combined}')
        
        scripts.build_genome_db.create_a_db(str(output.combined), str(output.combined_db))

rule create_features_db:
    input:
        "assets/reference/only_tsl_1_and_NA.gtf"
    output:
        "data/processed/features.db",
        'data/processed/longest_txpt_per_gene.db',
        'data/processed/longest_txpt_per_gene.gff'
    shell:
        "python scripts/build_genome_db.py {input}"
        
rule get_chrom_sizes:
    input:
        bam = SAMS_DIR + '/all_reads.bam',
        bai = SAMS_DIR + '/all_reads.bam.bai'
    output:
        fname = "data/processed/chrom.sizes"
    run:
        shell("samtools idxstats {input.bam} > {output}.raw")
        outf = open(output.fname, 'w')
        with open(f"{output.fname}.raw") as f:
            for li in f:
                if type(li) == type(''):
                    if li[0] != '*':
                        outf.write('\t'.join([str(x) for x in li.split('\t')[:2]]) + '\n')
        outf.close()
        
#########################################################
# References for read mapping.
#########################################################   

rule subset_gtf_to_only_tsl1_and_NA:      
    input:
        config['feature_gtf']
    output:
        "assets/reference/only_tsl_1_and_NA.gtf"
    shell:
        "python scripts/subset_gtf_to_only_tsl1_and_NA.py {input} {output}"
        
rule create_featureCounts_formatted_gtf_from_regular_gtf:
    input:
        "assets/reference/only_tsl_1_and_NA.gtf"
    output:
        "assets/reference/featureCounts_formatted.gtf"
    shell:
        "python scripts/gtf_to_featureCounts_formatted_gtf.py {input} {output}"

rule rDNA_gtf:
    input:
        rRNA_fa = "assets/reference/U13369_rRNA.fasta",
        gff = "assets/reference/U13369_rRNA.gff",
    output:
        gtf = 'data/processed/U13369_rRNA.gtf',
    shell:
        'python scripts/make_rDNA_gtf.py {input.gff} {output.gtf}'

rule mk_dir:
    output:
        "assets/rDNA_star_index/stops_error"
    shell:
        "mkdir -p assets/rDNA_star_index; touch assets/rDNA_star_index/stops_error"
        
rule rDNA_genome:
    input:
        fa = "assets/reference/U13369_rRNA.fasta",
        gtf = 'data/processed/U13369_rRNA.gtf',
        _dir = "assets/rDNA_star_index/stops_error"
    output:
        star_index = "assets/rDNA_star_index/Genome",
    threads:
        16
    shell:
        #"mkdir -p assets/rDNA_star_index/; " + \
        config['STAR'] + " --genomeSAindexNbases 5 --limitGenomeGenerateRAM 100000000000 --runThreadN 16 --runMode genomeGenerate --genomeDir assets/rDNA_star_index --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 75"
        
rule repeats_chromosome:
    input:
        embl = 'assets/reference/Dfam_curatedonly.embl',
    output:
        fa = 'assets/reference/repeats_chrom.fa',
        gtf = 'assets/reference/repeats.gtf',
        repeats = 'assets/repeats_star_index/Genome',
    threads:
        16
    run:
        ep = scripts.make_repeats_chrom.emblParser(input.embl)
        entries_df = ep.parse()
        ep.write_as_chromosome()
        os.makedirs('assets/repeats_star_index/', exist_ok=True)
        shell(config['STAR'] + f" --genomeSAindexNbases 5 --limitGenomeGenerateRAM 100000000000 --runThreadN 16 --runMode genomeGenerate --genomeDir assets/repeats_star_index --genomeFastaFiles {output.fa} --sjdbGTFfile {output.gtf} --sjdbOverhang 75")
# --genomeSAindexNbases 5 is well below the default 14. It speeds up the mapping to small genomes.

rule download_dfam_annotation:
    output:
        "assets/reference/Dfam_curatedonly.embl"
    run:
        #shell("cd assets/reference/")
        shell("wget https://www.dfam.org/releases/Dfam_3.3/families/Dfam_curatedonly.embl.gz")
        shell("gunzip Dfam_curatedonly.embl.gz")
        shell("mv Dfam_curatedonly.embl assets/reference/")
        #shell("cd ../..")

rule download_reference_gtf:
    output:
        config['feature_gtf']
    run:
        shell("echo feature gtf not found:", config['feature_gtf'])
        shell("echo downloading a default: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz")

        shell("wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz")
        shell("gunzip gencode.v29.basic.annotation.gtf.gz")
        shell("mv gencode.v29.basic.annotation.gtf assets/reference/")
        shell("echo setting config['feature_gtf'] to match this default gtf: assets/reference/gencode.v29.basic.annotation.gtf")
        config['feature_gtf'] = "assets/reference/gencode.v29.basic.annotation.gtf"
