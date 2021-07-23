import os, sys, re, glob, pandas, importlib, shutil
import scripts
import scripts.split_bam

#########################################################
# Read mapping, splitting and duplicate removal.
#########################################################    

rule mapping:
    input:
        'done/preprocess_reads.done'
    output:
        done = 'done/mapping.done',
        all_reads = SAMS_DIR + '/all_reads.bam',
    run:
        ex.mapping(clobber=True)
        shell("touch {output.done}")
        
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
        "umi_tools dedup --stdin={input.bam} --stdout={output.bam}" + \
        ' --log={output.log} --extract-umi-method=read_id;sleep 1;'
        "samtools index {output.bam}"
        
rule split_bam:
    input:
        bam = SAMS_DIR + '/all_reads.bam',
    output:
        bams = expand(SAMS_DIR + '/split/{sample}.bam', sample=samples)
    run:
        # Read name suffix format: 3773_AGCTAGAAAATCG_AGT_ACT_CGATTTTCTAGCT
        # AGCTAGAAAATCG_AGT -> AGCTAG=L5 BC. AGT=L3 BC. AAAATCG=UMI.
        scheme = ex.read_scheme()
        print('=' * 140)
        print(scheme.scheme_df)
        barcodes = ['_'.join([str(x),str(y)]) for x,y in \
                    zip(scheme.scheme_df['L5_BC'], scheme.scheme_df['L3_BC'])]
        print('~' * 140)
        print(barcodes)
        print('zip:', list(zip(scheme.scheme_df['L5_BC'], scheme.scheme_df['L3_BC'])))
        
        barcode_to_fname = {
            '_'.join(l5_l3): os.path.splitext(scheme.p6p3_to_long_filename_r1[l5_l3])[0] \
            for l5_l3 in zip(scheme.scheme_df['L5_BC'], scheme.scheme_df['L3_BC'])}
        
        print('-' * 140)
        print(barcode_to_fname)
        scripts.split_bam.split_bam(
            input.bam, barcodes, barcode_to_fname, SAMS_DIR + "/split/")


#########################################################
# Unmapped read pre-processing.
#########################################################   

rule read_scheme:
    output: "done/read_scheme.done"
    run: 
        ex.read_scheme()
        #shell("mkdir done")
        shell("touch {output}")

rule move_umis:
    input:
        in1 = config['R1_fastq'],
        in2 = config['R2_fastq'],
    output:
        out1 = FASTQ_DIR + '/umis_moved/R1.fastq.gz',
        out2 = FASTQ_DIR + '/umis_moved/R2.fastq.gz',
    run:
        ex.read_scheme()
        
        fq = FASTQ_DIR
        os.makedirs(f'{fq}/umis_moved/too_short/', exist_ok=True)
        min_length = 14
        
        basename1 = os.path.splitext(os.path.basename(input.in1))[0]
        basename2 = os.path.splitext(os.path.basename(input.in2))[0]
        cmd = 'cutadapt -A TACCCTTCGCTTCACACACAAGGGGAAAGAGTGTAGATCTCGGTGGTCGC' +\
        ' -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT' +\
        ' --pair-filter=any -u 13 -U 3'
        cmd += r" --rename='{id}_{r1.cut_prefix}-{r2.cut_prefix}" + \
        f' --too-short-output {fq}/umis_moved/too_short/{basename1}.gz'
        cmd += f' --too-short-paired-output {fq}/umis_moved/too_short/{basename2}.gz'
        cmd += f' --minimum-length {min_length} -o {output.out1} -p {output.out2}'
        cmd += f' {input.in1} {input.in2}'

        print(cmd)
        res = subprocess.check_output(cmd.split(' '))
        res = res.decode()
        #"umi_tools extract --stdin={input.in1} --read2-in={input.in2} --extract-method=regex" + \
        #r''' --bc-pattern="(?P<cell_1>.{{6}})(?P<umi_1>.{{7}})" --bc-pattern2="(?P<cell_2>.{{3}})"''' + \
        #" --stdout={output.out1} --read2-out={output.out2}"
        "umi_tools extract --stdin={input.in1} --read2-in={input.in2} --extract-method=string" + \
        r''' --bc-pattern=CCCCCCNNNNNNN --bc-pattern2=NNN''' + \
        " --stdout={output.out1} --read2-out={output.out2}"
        #N = UMI position (required)
        #C = cell barcode position (optional)
        #X = sample position (optional)
        

rule cut:
    input:
        in1 = FASTQ_DIR + '/umis_moved/R1.fastq.gz',
        in2 = FASTQ_DIR + '/umis_moved/R2.fastq.gz',
    output:
        out1 = FASTQ_DIR + '/ready_to_map/R1.fastq.gz',
        out2 = FASTQ_DIR + '/ready_to_map/R2.fastq.gz',
        done = 'done/preprocess_reads.done',
    run:
        ex.read_scheme()
        
        #fq = ex.file_paths['fastq']
        os.makedirs(f'{FASTQ_DIR}/ready_to_map/too_short/', exist_ok=True)
        min_length = 14
        
        basename1 = os.path.splitext(os.path.basename(input.in1))[0]
        basename2 = os.path.splitext(os.path.basename(input.in2))[0]
        cmd = 'cutadapt --pair-filter=any -U -13 -u -3' +\
        f' --too-short-output {FASTQ_DIR}/ready_to_map/too_short/{basename1}.gz'
        cmd += f' --too-short-paired-output {FASTQ_DIR}/ready_to_map/too_short/{basename2}.gz'
        cmd += f' --minimum-length {min_length} -o {output.out1} -p {output.out2}'
        cmd += f' {input.in1} {input.in2}'

        print(cmd)
        res = subprocess.check_output(cmd.split(' '))
        res = res.decode()
        
        shell("touch {output.done}")
        