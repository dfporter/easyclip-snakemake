import collections, pandas, os, re, glob, sys, importlib, pickle, subprocess, time, json, pysam
from typing import List, Tuple, Union, Mapping
from pprint import pprint
from pathlib import Path
from argparse import Namespace

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# This RepEnrich2 has to be the modified python3.8+ one:
# https://github.com/dfporter/RepEnrich2/tree/py3

#FIX ENVIRONMENT LATER
#from RepEnrich2 import RepEnrich2
        
import scripts
import scripts.scheme
import scripts.makeBedDataFile
import scripts.scheme_signal_RNAs
import scripts.bedAndFastaStats
import scripts.bedgraphs
import scripts.clip_adapters
import scripts.combine_fastqs_for_mapping
import scripts.sam_to_bed_and_wig
import scripts.collapse_duplicates


importlib.reload(scripts.scheme_signal_RNAs)
importlib.reload(scripts.makeBedDataFile)
importlib.reload(scripts.bedAndFastaStats)
importlib.reload(scripts.scheme)
importlib.reload(scripts.bedgraphs)
importlib.reload(scripts.clip_adapters)
importlib.reload(scripts.combine_fastqs_for_mapping)
importlib.reload(scripts.sam_to_bed_and_wig)
importlib.reload(scripts.collapse_duplicates)

# Special repeats are those in the categories of ['rRNA', 'scRNA', 'snRNA', 'tRNA'].
special_repeats = [
    'tRNA-Arg-CGA', 'tRNA-Ser-TCG', 'HY4', 'tRNA-Met-i', 'tRNA-Gln-CAG', 'U1', 'tRNA-Val-GTA', 'tRNA-Ser-TCY', 'tRNA-Val-GTY', 'tRNA-Gln-CAA_', 'U13_', 'tRNA-Ser-TCA', 'tRNA-Pro-CCA', 'U5', 'tRNA-Arg-CGG',
    'tRNA-Arg-AGA', 'tRNA-Arg-CGA_', 
    'tRNA-Met_', 'tRNA-Leu-CTY', 'SSU-rRNA_Hsa', 'tRNA-Ile-ATT', 'tRNA-Ile-ATC', 'BC200', 'tRNA-Leu-TTG', 'LSU-rRNA_Hsa', 'tRNA-Gly-GGY', 'tRNA-Leu-CTA', 
    'tRNA-Trp-TGG', 'tRNA-Thr-ACG_', 'tRNA-His-CAY', 'tRNA-Ala-GCY_', 'CRP1', 'tRNA-Arg-AGG', 'tRNA-Gln-CAA', '5S', 'tRNA-Glu-GAG', 'LFSINE_Vert', 'U4', 'tRNA-Asn-AAT', 'tRNA-Thr-ACG', 'tRNA-Lys-AAA', 'U6', 
    'tRNA-Pro-CCG', 'tRNA-Ser-AGY', 'U13', 'tRNA-Arg-CGY_', 'tRNA-Ala-GCG', 'tRNA-Ala-GCA', 'tRNA-Thr-ACY_', 'U7', 'tRNA-Phe-TTY', 'tRNA-Leu-CTA_', 'tRNA-Asn-AAC', 'tRNA-Thr-ACA', 'tRNA-Lys-AAG', 'U14', 
    'tRNA-Pro-CCY', 'tRNA-Gly-GGA', 'tRNA-Ser-TCA_', 'tRNA-Met', 'tRNA-His-CAY_', 'tRNA-Glu-GAA', 'tRNA-Leu-CTG', 'HY1', 'tRNA-Glu-GAG_', 'U3', 'tRNA-Tyr-TAC', 'HY3', 'tRNA-Thr-ACY', 'tRNA-Ala-GCY', 'tRNA-Gly-GGG', 'U2', 'tRNA-Cys-TGY',
    'tRNA-Asp-GAY', 'U17', 'tRNA-Ile-ATA', 'tRNA-Val-GTG', 'tRNA-Tyr-TAT', 'tRNA-Leu-TTA', 'HY5', 'U8'
    #'tRNA-SeC(e)-TGA', 'tRNA-Leu-TTA(m)', 'tRNA-Ser-TCA(m)',
]
snRNA = ["U6", "U2", "U13_", "U7", "U3", "U1", "U5", "U13", "U4", "U8", "U17", "U14",]

class starCaller():

    def __init__(self):
        pass

    def star_cmd_to_repeats(
        self, paths=None, read1_fname: str='', read2_fname: str='', threads=10):

        if paths is None:
            paths = self.file_paths

        cmd = paths['STAR']
        cmd += f" --genomeDir assets/repeats_star_index"
        cmd += f' --runThreadN {threads}'
        cmd += f' --readFilesIn {read1_fname} {read2_fname}'
        cmd += ' --alignIntronMax 1'  # Max intron size = 1. Setting to 0 causes the default behavior.
        cmd += ' --alignEndsType EndToEnd'  # Not completely sure this is right for repeats.
        cmd += ' --outReadsUnmapped Fastx'
        if os.path.splitext(read1_fname)[-1] == '.gz':
            cmd += ' --readFilesCommand zcat'
        # --outReadsUnmapped Fastx: output of unmapped and partially mapped (i.e. mapped only one mate
        # of a paired end read) reads in separate file(s).
        # output in separate fasta/fastq files, Unmapped.out.mate1/2
        cmd += f" --outFileNamePrefix {paths['sams']}/repeats."
        return cmd

    def star_cmd_to_genome(
        self, paths=None, read1_fname: str='', read2_fname: str='', threads=10):

        if paths is None:
            paths = self.file_paths

        cmd = paths['STAR']
        cmd += f" --genomeDir {paths['STAR_index']}"
        cmd += f' --runThreadN {threads} --limitOutSJcollapsed 2000000'
        cmd += f' --readFilesIn {read1_fname} {read2_fname}'
        cmd += ' --outReadsUnmapped Fastx'
        if os.path.splitext(read1_fname)[-1] == '.gz':
            cmd += ' --readFilesCommand zcat'
        # --outReadsUnmapped Fastx: output of unmapped and partially mapped (i.e. mapped only one mate
        # of a paired end read) reads in separate file(s).
        # output in separate fasta/fastq files, Unmapped.out.mate1/2
        cmd += f" --outFileNamePrefix {paths['sams']}/genome."
        return cmd
    

class mappingMethods():
    """General mapping methods for repeats and the genome.
    
    mapping(): Entry point method.
        Parameters:
            fastqfile1: str = '', fastqfile2: str = ''
        Outputs:
            Calls the selected mapping method.
    """

    def mapping(
        self, #fastq_directory: str = '',
        fastqfile1: str = '', fastqfile2: str = '',
        cpus: int=10, clobber: Union[str, bool] = True):
        """Map the given fastq files in the manner indicated by which_first.
        """

        print(self.file_paths)
        # Use default file paths if not given.
        if fastqfile1 == '':
            fastqfile1 = self.file_paths['cutadapt'] + '/R1.fastq.gz'
        if fastqfile2 == '':
            fastqfile2 = self.file_paths['cutadapt'] + '/R2.fastq.gz'

        # Initialize a log dict if needed.
        if not(hasattr(self, 'log')):
            self.log = {}

        # Define a sam folder if not already set.
        if 'sams' not in self.file_paths:
            self.file_paths['sams'] = self.file_paths['fastq'] + '/sam/'

        # Make any sam/bed directories that don't exist.
        for _dir in [self.file_paths['beds'], self.file_paths['sams'] + '/split/']:
            os.makedirs(_dir, exist_ok=True)

        # Mapping.
        self.map_all_reads_to_genome_and_repeats(
            fastqfile1=fastqfile1, fastqfile2=fastqfile2, cpus=cpus, clobber=clobber)


    def map_all_reads_to_genome_and_repeats(
        self, fastqfile1: str, fastqfile2: str, cpus: int=10,
        clobber: Union[str, bool] = False) -> None:

        sf = self.file_paths['sams']  # To make this more concise.

        # If either clobber is True for mapping to repeats, or the output file 
        # (combined) doesn't exist, map.
        if (clobber is True or clobber=='repeats') or (
            not os.path.exists(Path(sf, 'repeats.bam'))):

            # Get the string of the shell command to map to repeats with star.
            _starCaller = starCaller()
            cmd = _starCaller.star_cmd_to_repeats(
                paths=self.file_paths,
                read1_fname=fastqfile1, read2_fname=fastqfile2)

            # Map to repeats with star.
            print(cmd, '\n')
            subprocess.check_output(cmd.split(' '))
            print("Finished mapping.")
            
            subprocess.check_output(f"mv {sf}/repeats.Aligned.out.sam {sf}/repeats.sam".split(' '))
            self.filter_sort_and_index_bam(f"{sf}/repeats.sam")

        # If either clobber is True for the genome, or the output file (filtered)
        # from STAR mapping to the genome doesn't exist, map.
        if (clobber is True or clobber=='genome') or (
            not os.path.exists(Path(sf, 'genome.bam'))):
            
            print(f"clobber: {clobber} os.path.exists({sf}/genome.bam", os.path.exists(Path(sf, 'genome.bam')))
                                           
            unmapped_fastq_fname1 = f'{sf}/repeats.Unmapped.out.mate1'
            unmapped_fastq_fname2 = f'{sf}/repeats.Unmapped.out.mate2'
            
            # Get the string of the shell command to map to genome with star.
            _starCaller = starCaller()
            cmd = _starCaller.star_cmd_to_genome(
                paths=self.file_paths,
                read1_fname=unmapped_fastq_fname1, read2_fname=unmapped_fastq_fname2)

            # Map to repeats with star.
            print(cmd, '\n')
            subprocess.check_output(cmd.split(' '))
            print("Finished mapping.")
            
            subprocess.check_output(f"mv {sf}/genome.Aligned.out.sam {sf}/genome.sam".split(' '))
            self.filter_sort_and_index_bam(f"{sf}/genome.sam")

        # If clobber is not False, or the final output bam
        # doesn't exist, combine the outputs into all_reads.bam.
        # clobber=any_string will evaluate in the conditional below to True.
        if (clobber is not False) or (
            not os.path.exists(Path(sf, 'all_reads.bam'))):

            # Compare the outputs of the two independent mappings.
            # self.sam_repeats_filename: all repeats.
            # mapped_filtered_sam: genomic mappings (genome.filtered.bam).
            #self.combine_genomic_and_repeats_bams(
            #    genomic_bam=Path(sf, 'genome.filtered.bam'), repeats_bam=Path(sf, 'repeats.filtered.bam'),
            #    merged_filename=Path(sf, 'all_reads.bam'))

            #self.split_collapse_and_make_beds_and_bedgraphs_from_sam(   
            #    input_sam_file=sf + '/all_reads.bam')

            # Combine the genomic reads and the reads mapping to repeats.
            # -f forces the merge even if merged_filename exists.
            cmd = f'samtools merge -f {sf}/all_reads.bam {sf}/genome.filtered.bam {sf}/repeats.filtered.bam'
            self.proclaim(cmd)

            repeats, genome = f"{sf}/repeats.filtered.bam", f"{sf}/genome.filtered.bam"
            self.proclaim(f"samtools view -o {sf}/repeats.header -H {repeats}") 
            self.proclaim(f"samtools view -o {sf}/genome.header -H {genome}")
            self.proclaim(f"samtools view -o {sf}/repeats.filtered.sam {repeats}")
            self.proclaim(f"samtools view -o {sf}/genome.filtered.sam {genome}")

            header = self.insert_sq_line(f"{sf}/genome.header", self.grab_sq_line(f"{sf}/repeats.header"))
            with open(f"{sf}/combined.header.sam", 'w') as f:
                f.write(header)
            self.proclaim(f"cat {sf}/combined.header.sam {sf}/genome.filtered.sam {sf}/repeats.filtered.sam > {sf}/both.filtered.sam")
            self.proclaim(f"samtools view -o {sf}/all_reads.bam {sf}/both.filtered.sam")
            self.proclaim(f"samtools sort -o {sf}/both.bam {sf}/all_reads.bam")
            self.proclaim(f"mv {sf}/both.bam {sf}/all_reads.bam")
            self.proclaim(f"samtools index {sf}/all_reads.bam")

    def merge_repeats_and_genome(self, sf, prefix):
            # Combine the genomic reads and the reads mapping to repeats.
            # -f forces the merge even if merged_filename exists.
            cmd = f'samtools merge -f {sf}/{prefix}all_reads.bam {sf}/{prefix}genome.filtered.bam {sf}/{prefix}repeats.filtered.bam'
            self.proclaim(cmd)

            repeats, genome = f"{sf}/{prefix}repeats.filtered.bam", f"{sf}/{prefix}genome.filtered.bam"
            self.proclaim(f"samtools view -o {sf}/{prefix}repeats.header -H {repeats}") 
            self.proclaim(f"samtools view -o {sf}/{prefix}genome.header -H {genome}")
            self.proclaim(f"samtools view -o {sf}/{prefix}repeats.filtered.sam {repeats}")
            self.proclaim(f"samtools view -o {sf}/{prefix}genome.filtered.sam {genome}")

            header = self.insert_sq_line(f"{sf}/{prefix}genome.header", self.grab_sq_line(f"{sf}/{prefix}repeats.header"))
            with open(f"{sf}/{prefix}combined.header.sam", 'w') as f:
                f.write(header)
            self.proclaim(f"cat {sf}/{prefix}combined.header.sam {sf}/{prefix}genome.filtered.sam {sf}/{prefix}repeats.filtered.sam > {sf}/{prefix}both.filtered.sam")
            self.proclaim(f"samtools view -o {sf}/{prefix}all_reads.bam {sf}/{prefix}both.filtered.sam")
            self.proclaim(f"samtools sort -o {sf}/{prefix}both.bam {sf}/{prefix}all_reads.bam")
            self.proclaim(f"mv {sf}/{prefix}both.bam {sf}/{prefix}all_reads.bam")
            self.proclaim(f"samtools index {sf}/{prefix}all_reads.bam")

    def filter_sort_and_index_bam(self, samfname):
        base = os.path.splitext(samfname)[0]
        
        cmd = f'samtools view -h -o {base}.bam {base}.sam'
        print(cmd)
        subprocess.check_output(cmd.split(' '))

        # The cmd below excludes secondaries (-F 256) and unmapped reads (-F 4). 
        # It requires first mates (-f 64) and excludes low MAPQ (-q 10).
        # It prints the header (-h) and outputs (-o) to a file rather than STDOUT.
        cmd = 'samtools view -h -F 256 -f 64 -F 4 -q 10'
        cmd += f" -o {base}.filtered.bam" # Output filename
        cmd += f" {base}.bam"  # Input.

        print(cmd)
        subprocess.check_output(cmd.split(' '))

        # Sort and index the bam file.
        cmd = f"samtools sort {base}.filtered.bam > {base}.filtered.bam.sorted"
        self.proclaim(cmd)
        cmd = f"mv {base}.filtered.bam.sorted {base}.filtered.bam"
        self.proclaim(cmd)
        cmd = f"samtools index {base}.filtered.bam"
        self.proclaim(cmd)
    
    # Print and run a system command.
    def proclaim(self, cmd):
        print(cmd)
        os.system(cmd)
        
    def insert_sq_line(self, header_file, sq_line):
        line_n = 0
        sq_start = -1
        sq_end = -1
        lines = []
        with open(header_file) as f:
            for li in f:
                if li[:3] == '@SQ':
                    if sq_start == 0:
                        sq_start = line_n
                if li[:3] == '@PG' and sq_end == -1:
                    sq_end = line_n 
                line_n += 1
                lines.append(li)
                if li[0] != '@':
                    break
        lines = lines[:sq_end] + [sq_line] + lines[sq_end:]
        return ''.join(lines)
    
    def grab_sq_line(self, samfname):
        with open(samfname) as f:
            for li in f:
                if li[:3] == '@SQ':
                    break
        return li