import os, sys, re, glob, pandas, importlib, shutil, subprocess
import scripts
import scripts.split_bam

#########################################################
# Read mapping, splitting and duplicate removal. 
#########################################################    
rule map_to_rDNA:
    input:
        fq1 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        fq2 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
        star_repeats_genome = "assets/rDNA_star_index/Genome",
    output:
        by_index = expand(SAMS_DIR + '/{pcr_index}rDNA.bam', pcr_index=PCR_INDEX_SET),
        filt = expand(SAMS_DIR + '/{pcr_index}rDNA.filtered.bam', pcr_index=PCR_INDEX_SET),
        filt_bai = expand(SAMS_DIR + '/{pcr_index}rDNA.filtered.bam.bai', pcr_index=PCR_INDEX_SET),
    threads:
        20
    run:
        for input_r1, input_r2 in zip(input.fq1, input.fq2):
            
            pcr_index = os.path.basename(input_r1).split('R1.fastq.gz')[0]           

            cmd = config['STAR']
            cmd += f" --genomeDir assets/rDNA_star_index"
            cmd += ' --runThreadN ' + str(threads)
            cmd += f' --readFilesIn {input_r1} {input_r2}'
            cmd += ' --alignIntronMax 1'  # Max intron size = 1. Setting to 0 causes the default behavior.
            cmd += ' --alignEndsType EndToEnd'  # Not completely sure this is right for repeats.
            cmd += ' --outReadsUnmapped Fastx'
            if os.path.splitext(input_r1)[-1] == '.gz':
                cmd += ' --readFilesCommand zcat'
            cmd += f" --outFileNamePrefix {SAMS_DIR}/{pcr_index}rDNA."
        
            # Map to repeats with star.
            print(cmd, '\n')
            subprocess.check_output(cmd.split(' '))
            print("Finished mapping.")
            
            subprocess.check_output(
                f"mv {SAMS_DIR}/{pcr_index}rDNA.Aligned.out.sam {SAMS_DIR}/{pcr_index}rDNA.sam".split(' '))
            
            mapper = scripts.mapping.mappingMethods()
            mapper.file_paths = {k:v for k,v in config.items()}
            mapper.filter_sort_and_index_bam(f"{SAMS_DIR}/{pcr_index}rDNA.sam")
            
rule split_rDNA:
    input:
        fq1 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        fq2 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
        bam = SAMS_DIR + '/all_reads.bam',    
        samples = config['samples'],  # Filename of sample info.
        by_index = expand(SAMS_DIR + '/{pcr_index}rDNA.filtered.bam', pcr_index=PCR_INDEX_SET),
        bais = expand(SAMS_DIR + '/{pcr_index}rDNA.filtered.bam.bai', pcr_index=PCR_INDEX_SET),
    params:
        # E.g., BBBBBBNNNNNNN. B=barcode base, N=UMI base.
        l5_inline = config['L5_inline'],
        # E.g., BBBBBNNNNN
        l3_inline = config['L3_inline'],
    output:
        bams = expand(SAMS_DIR + '/rDNA_split/{sample}.bam', sample=samples)
    run:
        df = pandas.read_csv(config['samples'], sep='\t')
        df['Gene'] = [re.sub(' ', '-', x) for x in df['Gene']]  # Get rid of spaces.
        df = df.loc[[type(x)==type('') for x in df['Gene']], :]

        # Define the samples list used throughout the workflow.
        samples = [f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}" for exp,protein,l5_bc,l3_bc,rep in zip(
            df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'])]
        
        # Read name suffix format example: 3773_AGCTAGAAAATCG_AGT
        # AGCTAGAAAATCG_AGT -> AGCTAG=L5 BC. AGT=L3 BC. AAAATCG=UMI.

        print('=' * 140)
        print(df.head())
        
        barcodes = []
        barcode_to_fname = {}
        for exp,protein,l5_bc,l3_bc,rep,r1_fastq in zip(
            df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'], df['R1_fastq']):

            pcr_prefix = r1_fastq.split('R1.fastq.gz')[0]
            barcode = f"{l5_bc}__{l3_bc}|{pcr_prefix}"
            fname = f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}"
            barcode_to_fname[barcode] = fname
            barcodes.append(barcode)
        
        print('-' * 140)
        print(barcode_to_fname)

        for input_r1, input_r2 in zip(input.fq1, input.fq2):
            
            pcr_index = os.path.basename(input_r1).split('R1.fastq.gz')[0]  
            
            _in_bam = SAMS_DIR + f'/{pcr_index}rDNA.filtered.bam'
            scripts.split_bam.split_bam(
                _in_bam, barcodes, barcode_to_fname, SAMS_DIR + "/rDNA_split/",
                l5_inline_pattern=str(params.l5_inline), l3_inline_pattern=str(params.l3_inline))
            
            
rule dedup_rDNA:
    input:
        bam = SAMS_DIR + '/rDNA_split/{sample}.bam',
    output:
        bam = SAMS_DIR + '/rDNA_dedup/{sample}.bam',  # Input to clipper.
        log = SAMS_DIR + '/rDNA_dedup/log/{sample}.log',
        bai = SAMS_DIR + '/rDNA_dedup/{sample}.bam.bai',
    conda:
        '../envs/umitools.yml'
    shell:
        "umi_tools dedup --stdin={input.bam} --stdout={output.bam}"
        ' --log={output.log} --extract-umi-method=read_id;sleep 2;'
        "samtools index {output.bam}"
        
rule mapping:
    input:
        fq1 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        fq2 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
        star_repeats_genome = "assets/repeats_star_index/Genome",
    output:
        by_index = expand(SAMS_DIR + '/{pcr_index}all_reads.bam', pcr_index=PCR_INDEX_SET),
        all_reads = SAMS_DIR + "/all_reads.bam",
    threads:
        20
    run:
        for input_r1, input_r2 in zip(input.fq1, input.fq2):
            
            pcr_index = os.path.basename(input_r1).split('R1.fastq.gz')[0]           

            cmd = config['STAR']
            cmd += f" --genomeDir assets/repeats_star_index"
            cmd += ' --runThreadN ' + str(threads)
            cmd += f' --readFilesIn {input_r1} {input_r2}'
            cmd += ' --alignIntronMax 1'  # Max intron size = 1. Setting to 0 causes the default behavior.
            cmd += ' --alignEndsType EndToEnd'  # Not completely sure this is right for repeats.
            cmd += ' --outReadsUnmapped Fastx'
            if os.path.splitext(input_r1)[-1] == '.gz':
                cmd += ' --readFilesCommand zcat'
            cmd += f" --outFileNamePrefix {SAMS_DIR}/{pcr_index}repeats."
        
            # Map to repeats with star.
            print(cmd, '\n')
            subprocess.check_output(cmd.split(' '))
            print("Finished mapping.")
            
            subprocess.check_output(
                f"mv {SAMS_DIR}/{pcr_index}repeats.Aligned.out.sam {SAMS_DIR}/{pcr_index}repeats.sam".split(' '))
            
            mapper = scripts.mapping.mappingMethods()
            mapper.file_paths = {k:v for k,v in config.items()}
            mapper.filter_sort_and_index_bam(f"{SAMS_DIR}/{pcr_index}repeats.sam")
        
            unmapped_fastq_fname1 = f'{SAMS_DIR}/{pcr_index}repeats.Unmapped.out.mate1'
            unmapped_fastq_fname2 = f'{SAMS_DIR}/{pcr_index}repeats.Unmapped.out.mate2'
            
            # Get the string of the shell command to map to genome with star.
            cmd = config['STAR']
            cmd += f" --genomeDir {config['STAR_index']}"
            cmd += ' --runThreadN ' + str(threads)
            cmd += ' --limitOutSJcollapsed 2000000'
            cmd += f' --readFilesIn {unmapped_fastq_fname1} {unmapped_fastq_fname2}'
            cmd += ' --outReadsUnmapped Fastx'
            if os.path.splitext(unmapped_fastq_fname1)[-1] == '.gz':
                cmd += ' --readFilesCommand zcat'
            # --outReadsUnmapped Fastx: output of unmapped and partially mapped (i.e. mapped only one mate
            # of a paired end read) reads in separate file(s).
            # output in separate fasta/fastq files, Unmapped.out.mate1/2
            cmd += f" --outFileNamePrefix {SAMS_DIR}/{pcr_index}genome."

            # Map to genome with star.
            print(cmd, '\n')
            subprocess.check_output(cmd.split(' '))
            print("Finished mapping.")
            
            subprocess.check_output(
                f"mv {SAMS_DIR}/{pcr_index}genome.Aligned.out.sam {SAMS_DIR}/{pcr_index}genome.sam".split(' '))
            mapper.filter_sort_and_index_bam(f"{SAMS_DIR}/{pcr_index}genome.sam")
            
            # This writes {SAMS_DIR}/{pcr_index}all_reads.bam".
            mapper.merge_repeats_and_genome(SAMS_DIR, pcr_index)
            
        # Merge the bam files representing the different PCR indexes.
        cmd = f'samtools merge -f {SAMS_DIR}/all_reads.unsorted.bam ' + ' '.join([
            f"{SAMS_DIR}/{pcr_index}all_reads.bam" for pcr_index in PCR_INDEX_SET])
        print(cmd)
        subprocess.check_output(cmd.split(' '))
        
        # Sort and index the bam file.
        os.system(f"samtools sort {SAMS_DIR}/all_reads.unsorted.bam > {SAMS_DIR}/all_reads.bam")
        os.system(f"samtools index {SAMS_DIR}/all_reads.bam")
        os.system(f"rm {SAMS_DIR}/all_reads.unsorted.bam")
        
rule dedup:
    input:
        bam = SAMS_DIR + '/split/{sample}.bam'
    output:
        bam = SAMS_DIR + '/dedup/{sample}.bam',  # Input to clipper.
        log = SAMS_DIR + '/dedup/log/{sample}.log',
        bai = SAMS_DIR + '/dedup/{sample}.bam.bai',
    conda:
        '../envs/umitools.yml'
    shell:
        "umi_tools dedup --stdin={input.bam} --stdout={output.bam}"
        ' --log={output.log} --extract-umi-method=read_id;sleep 1;'
        "samtools index {output.bam}"
        
rule split_bam:
    input:
        bam = SAMS_DIR + '/all_reads.bam',    
        samples = config['samples'],  # Filename of sample info.
    params:
        # E.g., BBBBBBNNNNNNN. B=barcode base, N=UMI base.
        l5_inline = config['L5_inline'],
        # E.g., BBBBBNNNNN
        l3_inline = config['L3_inline'],
    output:
        bams = expand(SAMS_DIR + '/split/{sample}.bam', sample=samples)
    run:
        df = pandas.read_csv(config['samples'], sep='\t')
        df['Gene'] = [re.sub(' ', '-', x) for x in df['Gene']]  # Get rid of spaces.
        df = df.loc[[type(x)==type('') for x in df['Gene']], :]

        # Define the samples list used throughout the workflow.
        samples = [f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}" for exp,protein,l5_bc,l3_bc,rep in zip(
            df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'])]
        
        # Read name suffix format example: 3773_AGCTAGAAAATCG_AGT
        # AGCTAGAAAATCG_AGT -> AGCTAG=L5 BC. AGT=L3 BC. AAAATCG=UMI.

        print('=' * 140)
        print(df.head())
        
        barcodes = []
        barcode_to_fname = {}
        for exp,protein,l5_bc,l3_bc,rep,r1_fastq in zip(
            df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'], df['R1_fastq']):

            pcr_prefix = r1_fastq.split('R1.fastq.gz')[0]
            barcode = f"{l5_bc}__{l3_bc}|{pcr_prefix}"
            fname = f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}"
            barcode_to_fname[barcode] = fname
            barcodes.append(barcode)
        
        print('-' * 140)
        print(barcode_to_fname)
        scripts.split_bam.split_bam(
            input.bam, barcodes, barcode_to_fname, SAMS_DIR + "/split/",
            l5_inline_pattern=str(params.l5_inline), l3_inline_pattern=str(params.l3_inline))


#########################################################
# Unmapped read pre-processing.
#########################################################   

rule second_adapter_trim:
    input:
        in1 = expand(FASTQ_DIR + '/umis_moved/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        in2 = expand(FASTQ_DIR + '/umis_moved/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
    output:
        out1 = expand(FASTQ_DIR + '/second_cut/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        out2 = expand(FASTQ_DIR + '/second_cut/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
    threads:
        20
    run:
        # Parameters:
        # E.g., AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT
        adapter1_short_len = min([12, len(config['cutadapt_r1_adapter_seq_to_trim'])])
        adapter1 = config['cutadapt_r1_adapter_seq_to_trim'][:adapter1_short_len]
        # E.g., TACCCTTCGCTTCACACACAAGGGGAAAGAGTGTAGATCTCGGTGGTCGC
        adapter2_short_len = min([12, len(config['cutadapt_r2_adapter_seq_to_trim'])])
        adapter2 = config['cutadapt_r2_adapter_seq_to_trim'][:adapter2_short_len]
        # E.g., BBBBBBNNNNNNN. B=barcode base, N=UMI base.
        l5_inline = config['L5_inline']
        # E.g., BBBBBNNNNN
        l3_inline = config['L3_inline']        
        
        fq = FASTQ_DIR
        os.makedirs(f'{fq}/second_cut/too_short/', exist_ok=True)
        min_length = 14
        
        bases_l5_to_move = len(l5_inline)
        bases_l3_to_move = len(l3_inline)
        
        for input_r1, input_r2, pcr_index, output_r1, output_r2 in zip(
            input.in1, input.in2, PCR_INDEX_SET, output.out1, output.out2):
            print(input_r1, "<-r1")
            basename1 = os.path.splitext(os.path.basename(input_r1))[0]
            basename2 = os.path.splitext(os.path.basename(input_r2))[0]

            # Construct the cutadapt command.
            # Say, -u 13 -U 3 for original adapters. -j 8 -> use 8 threads.
            cmd = f'cutadapt -A {adapter2} -a {adapter1}' + \
            f' --pair-filter=any -j {threads}' + \
            f' --too-short-output {fq}/second_cut/too_short/{basename1}.gz' + \
            f' --too-short-paired-output {fq}/second_cut/too_short/{basename2}.gz'
            cmd += f' --minimum-length {min_length} -o {output_r1} -p {output_r2}'
            cmd += f' {input_r1} {input_r2}'

            print('%' * 100)
            print(cmd)
            res = subprocess.check_output(cmd.split(' '))
            res = res.decode()
            
rule move_umis:
    input:
        in1 = [config['Fastq_folder'].rstrip('/') + f"/{pcr_index}R1.fastq.gz" for pcr_index in PCR_INDEX_SET],
        in2 = [config['Fastq_folder'].rstrip('/') + f"/{pcr_index}R2.fastq.gz" for pcr_index in PCR_INDEX_SET],
    output:
        out1 = expand(FASTQ_DIR + '/umis_moved/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        out2 = expand(FASTQ_DIR + '/umis_moved/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
    threads:
        20
    run:
        print(f"PCR_INDEX_SET={PCR_INDEX_SET}")
        print('threads:', threads)
        # Parameters:
        # E.g., AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT
        adapter1 = config['cutadapt_r1_adapter_seq_to_trim']   
        # E.g., TACCCTTCGCTTCACACACAAGGGGAAAGAGTGTAGATCTCGGTGGTCGC
        adapter2 = config['cutadapt_r2_adapter_seq_to_trim']
        # E.g., BBBBBBNNNNNNN. B=barcode base, N=UMI base.
        l5_inline = config['L5_inline']
        # E.g., BBBBBNNNNN
        l3_inline = config['L3_inline']        
        
        fq = FASTQ_DIR
        os.makedirs(f'{fq}/umis_moved/too_short/', exist_ok=True)
        min_length = 14
        
        bases_l5_to_move = len(l5_inline)
        bases_l3_to_move = len(l3_inline)
        
        for input_r1, input_r2, pcr_index, output_r1, output_r2 in zip(
            input.in1, input.in2, PCR_INDEX_SET, output.out1, output.out2):
            print(input_r1, "<-r1")
            basename1 = os.path.splitext(os.path.basename(input_r1))[0]
            basename2 = os.path.splitext(os.path.basename(input_r2))[0]

            # Construct the cutadapt command.
            # Say, -u 13 -U 3 for original adapters. -j 8 -> use 8 threads.
            cmd = f'cutadapt -A {adapter2} -a {adapter1}' + \
            f' --pair-filter=any -u {bases_l5_to_move} -U {bases_l3_to_move}'
            cmd += " -j " + str(threads)


            print('%' * 100)
            cmd += r" --rename={id}__{r1.cut_prefix}-{r2.cut_prefix}|" + pcr_index + \
            f' --too-short-output {fq}/umis_moved/too_short/{basename1}.gz'
        
            cmd += f' --too-short-paired-output {fq}/umis_moved/too_short/{basename2}.gz'
            cmd += f' --minimum-length {min_length} -o {output_r1} -p {output_r2}'
            cmd += f' {input_r1} {input_r2}'

            print(cmd)
            res = subprocess.check_output(cmd.split(' '))
            res = res.decode()
        

rule cut:
    input:
        #in1 = FASTQ_DIR + '/umis_moved/R1.fastq.gz',
        #in2 = FASTQ_DIR + '/umis_moved/R2.fastq.gz',
        #in1 = expand(FASTQ_DIR + '/second_cut/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        #in2 = expand(FASTQ_DIR + '/second_cut/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
        in1 = expand(FASTQ_DIR + '/umis_moved/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        in2 = expand(FASTQ_DIR + '/umis_moved/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
    params:
        # E.g., BBBBBBNNNNNNN. B=barcode base, N=UMI base.
        l5_inline = config['L5_inline'],
        # E.g., BBBBBNNNNN
        l3_inline = config['L3_inline'],
    output:
        out1 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R1.fastq.gz', pcr_index=PCR_INDEX_SET),
        out2 = expand(FASTQ_DIR + '/ready_to_map/{pcr_index}R2.fastq.gz', pcr_index=PCR_INDEX_SET),
    threads:
        20
    run:

        bases_l5_to_move = len(str(params.l5_inline))
        bases_l3_to_move = len(str(params.l3_inline))        
        
        #fq = ex.file_paths['fastq']
        os.makedirs(f'{FASTQ_DIR}/ready_to_map/too_short/', exist_ok=True)
        min_length = 14
        
        for input_r1, input_r2, output_r1, output_r2 in zip(input.in1, input.in2, output.out1, output.out2):
            basename1 = os.path.splitext(os.path.basename(input_r1))[0]
            basename2 = os.path.splitext(os.path.basename(input_r2))[0]
            
            cmd = f'cutadapt --pair-filter=any -U -{bases_l5_to_move} -u -{bases_l3_to_move}' + \
                " -j " + str(threads) + \
                f' --too-short-output {FASTQ_DIR}/ready_to_map/too_short/{basename1}.gz'
            cmd += f' --too-short-paired-output {FASTQ_DIR}/ready_to_map/too_short/{basename2}.gz'
            cmd += f' --minimum-length {min_length} -o {output_r1} -p {output_r2}'
            cmd += f' {input_r1} {input_r2}'

            print(cmd)
            res = subprocess.check_output(cmd.split(' '))
            res = res.decode()
        
        