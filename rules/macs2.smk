rule merge_replicates_peak_calls:
    input:
        peaks = expand("outs/macs2_peaks/{sample}_peaks.narrowPeak", sample = samples["sample"]),
    output:
        merged = expand("outs/merged_macs2/{condition}_peaks.narrowPeak", condition = samples["condition"]),
    run:
        for condition, _samples in condition_to_samples.items():
            print('---')
            _files = [f"outs/macs2_peaks/{sample}_peaks.narrowPeak" for sample in _samples]
            print(condition, _samples)
            if len(_files) == 4:
                temp_file1 = f"outs/{os.path.basename(_files[0])}_{os.path.basename(_files[1])}.temp"
                temp_file2 = f"outs/{os.path.basename(_files[2])}_{os.path.basename(_files[3])}.temp"
                shell(f"bedtools intersect -a {_files[0]} -b {_files[1]} > {temp_file1}")
                shell(f"bedtools intersect -a {_files[2]} -b {_files[3]} > {temp_file2}")
                shell(f"bedtools intersect -a {temp_file1} -b {temp_file2} > outs/merged_macs2/{condition}_peaks.narrowPeak")
                shell(f"rm {temp_file1}")
                shell(f"rm {temp_file2}")
            if len(_files) == 1:
                shell(f"cp {_files[0]} outs/merged_macs2/{condition}_peaks.narrowPeak")
            if len(_files) == 2:
                shell(f"bedtools intersect -a {_files[0]} -b {_files[1]} "
                      "> outs/merged_macs2/{condition}_peaks.narrowPeak")
                #sample = os.path.basename(str(input)).split('.fa')[0]
        
rule call_homer:
    input:
        peak_fasta = "outs/filtered_macs2_peaks/fastas/{sample}.fa",
        peak_fasta_point = "outs/filtered_macs2_peaks/point/fastas/{sample}.fa",
        randoms = "outs/filtered_macs2_peaks/fastas/randomControls/{sample}.fa",
        point_randoms = "outs/filtered_macs2_peaks/point/fastas/randomControls/{sample}.fa",
    output:
        homer_results = "outs/filtered_macs2_peaks/homer/{sample}/{sample}",
        homer_results_point = "outs/filtered_macs2_peaks/point/homer/{sample}/{sample}"
    run:
        sample = os.path.basename(str(input.peak_fasta)).split('.fa')[0]
        cmd = f"homer2 denovo -i {str(input.peak_fasta)} -b {str(input.randoms)}"
        cmd += f" -len 6 -S 10 -strand + -o {str(output.homer_results)}"
        shell(cmd)
        cmd = f"homer2 denovo -i {str(input.peak_fasta_point)} -b {str(input.point_randoms)}"
        cmd += f" -len 6 -S 10 -strand + -o {str(output.homer_results_point)}"
        shell(cmd)
        
rule call_homer_merged:
    input:
        peak_fasta = "outs/filtered_merged_macs2/fastas/{condition}.fa",
        peak_fasta_point = "outs/filtered_merged_macs2/point/fastas/{condition}.fa",
        randoms = "outs/filtered_merged_macs2/fastas/randomControls/{condition}.fa",
        point_randoms = "outs/filtered_merged_macs2/point/fastas/randomControls/{condition}.fa",
    output:
        homer_results = "outs/filtered_merged_macs2/homer/{condition}/{condition}",
        homer_results_point = "outs/filtered_merged_macs2/point/homer/{condition}/{condition}"
    run:
        #condition = os.path.basename(str(input.peak_fasta)).split('.fa')[0]
        cmd = f"homer2 denovo -i {str(input.peak_fasta)} -b {str(input.randoms)}"
        cmd += f" -len 6 -S 10 -strand + -o {str(output.homer_results)}"
        shell(cmd)
        cmd = f"homer2 denovo -i {str(input.peak_fasta_point)} -b {str(input.point_randoms)}"
        cmd += f" -len 6 -S 10 -strand + -o {str(output.homer_results_point)}"
        shell(cmd)    
        
rule filter_merged_peaks:
    input:
        peaks = expand("outs/merged_macs2/{condition}_peaks.narrowPeak", condition = samples["condition"]),
        peak_fastas = expand("outs/merged_macs2/fastas/{condition}.fa", condition = samples["condition"])
    output:
        subset_peaks = expand(
            "outs/filtered_merged_macs2/{condition}_peaks.narrowPeak", condition = samples["condition"]),
        subset_peak_fastas = expand(
            "outs/filtered_merged_macs2/fastas/{condition}.fa", condition = samples["condition"]),
        subset_peak_fastas_randoms = expand(
            "outs/filtered_merged_macs2/fastas/randomControls/{condition}.fa", condition = samples["condition"]),
        point_source = expand(
            "outs/filtered_merged_macs2/point/{condition}_peaks.narrowPeak", condition = samples["condition"]),
        point_source_fastas = expand(
            "outs/filtered_merged_macs2/point/fastas/{condition}.fa", condition = samples["condition"]),
        point_source_fastas_randoms = expand(
            "outs/filtered_merged_macs2/point/fastas/randomControls/{condition}.fa", condition = samples["condition"]),
    run:
        shell(f"python scripts/filter_peaks.py outs/merged_macs2/ outs/filtered_merged_macs2/ 9")

rule filter_peaks:
    input:
        peaks = expand("outs/macs2_peaks/{sample}_peaks.narrowPeak", sample = samples["sample"]),
        peak_fastas = expand("outs/macs2_peaks/fastas/{sample}.fa", sample = samples["sample"])
    output:
        subset_peaks = expand(
            "outs/filtered_macs2_peaks/{sample}_peaks.narrowPeak", sample = samples["sample"]),
        subset_peak_fastas = expand(
            "outs/filtered_macs2_peaks/fastas/{sample}.fa", sample = samples["sample"]),
        subset_peak_fastas_randoms = expand(
            "outs/filtered_macs2_peaks/fastas/randomControls/{sample}.fa", sample = samples["sample"]),
        point_source = expand(
            "outs/filtered_macs2_peaks/point/{sample}_peaks.narrowPeak", sample = samples["sample"]),
        point_source_fastas = expand(
            "outs/filtered_macs2_peaks/point/fastas/{sample}.fa", sample = samples["sample"]),
        point_source_fastas_randoms = expand(
            "outs/filtered_macs2_peaks/point/fastas/randomControls/{sample}.fa", sample = samples["sample"]),
    run:
        shell(f"python scripts/filter_peaks.py outs/macs2_peaks/ outs/filtered_macs2_peaks/ 9")

rule write_merged_peak_fastas:
    input:
        peaks = expand("outs/merged_macs2/{condition}_peaks.narrowPeak", condition = samples["condition"]),
    output:
        peak_fastas = expand("outs/merged_macs2/fastas/{condition}.fa", condition = samples["condition"]),
        randoms = expand("outs/merged_macs2/fastas/randomControls/{condition}.fa", condition = samples["condition"])
    run:
        shell(f"python scripts/fasta_under_peaks.py --fasta {config['genomic_fasta']}"
              f" --peaks_dir outs/merged_macs2/ --out_folder outs/merged_macs2/fastas/")
        
rule write_fastas:
    input:
        peaks = expand("outs/macs2_peaks/{sample}_peaks.narrowPeak", sample = samples["sample"])
    output:
        peak_fastas = expand("outs/macs2_peaks/fastas/{sample}.fa", sample = samples["sample"]),
        randoms = expand("outs/macs2_peaks/fastas/randomControls/{sample}.fa", sample = samples["sample"])
    run:
        shell(f"python scripts/fasta_under_peaks.py --fasta {config['genomic_fasta']}"
              f" --peaks_dir outs/macs2_peaks/ --out_folder outs/macs2_peaks/fastas/")

rule macs2_peak_calling:
    input:
        bws = "outs/align-merge/{sample}.sorted.bam"
    output:
        peaks = "outs/macs2_peaks/{sample}_peaks.narrowPeak",
        #peaks_dir = "outs/macs2_peaks"
    run:
        sample = os.path.basename(str(input)).split('.sorted')[0]
        #shell(f"macs2 callpeak -t {input}  -f BAM -n {sample} --outdir outs/macs2_peaks/ --nomodel -q 0.1 --bw 100;")
        shell(f"macs2 callpeak -t {input}  -f BAM -n {sample} -c assets/inputs/Input_RBFOX2_HepG2_hg38.bam"
              f" --outdir outs/macs2_peaks/ --nomodel -q 0.01 --bw 100;")
