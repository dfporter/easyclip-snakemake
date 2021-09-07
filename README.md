# Readme

This software is a simplification, formalization, and streamlining of easyCLIP read processing, mapping, and initial analysis.
It is intended to be easily usable.

The code used for the paper https://pubmed.ncbi.nlm.nih.gov/33692367/, warts and all, is placed in https://github.com/dfporter/easyclip.

#### Install workflow

Clone this repositiory into a directory

```bash
git clone git@github.com:dfporter/easyclip-snakemake.git
```

Create and activate the conda environment:

```bash
cd easyclip-snakemake
conda env create -f=envs/conda.yml -n easyclip-env
conda activate easyclip-env
```

To prevent errors running snakemake with --use-conda, install mamba:
```bash
conda install mamba -n base -c conda-forge
```


To install CLIPPER:
```bash

# In the easyCLIP directory:
git clone https://github.com/YeoLab/clipper
cd clipper
conda env create -f environment3.yml
conda activate clipper3
python setup.py install

# Replace the samtools import version in the environment3.yml file with 1.6
# (or some other version that works) to prevernt an error message whem samtools called.
# Alternatively, run this line:
conda install -c bioconda samtools=1.6 --force-reinstall

# Return to the easyclip-env conda environment before running the workflow.
cd ..
conda deactivate
conda activate easyclip-snakemake
```

#### Config files and UMIs/in-line barcodes

This is the example configfile:

```yaml
name: testdata
top_dir: testdata
scheme: testdata/samples.txt
Fastq_folder: testdata/fastq/raw
cutadapt_r1_adapter_seq_to_trim: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT
cutadapt_r2_adapter_seq_to_trim: TACCCTTCGCTTCACACACAAGGGGAAAGAGTGTAGATCTCGGTGGTCGC
L5_inline: BBBBBBNNNNNNN
L3_inline: BBB
STAR_index: /oak/stanford/groups/khavari/users/dfporter/pre_2021_projects/genome/GRCh38.gencode.29/star_index/
STAR: STAR
clipper: /storage/khavari/data/dfporter/anaconda3/envs/clipper3/bin/clipper
feature_gtf: assets/reference/gencode.v29.basic.annotation.gtf
genomic_fasta: assets/reference/GRCh38.primary_assembly.genome.fa
```

The two rules used for read preprocessing are in rules/read_preprocessing_and_mapping.smk - "move_umis" and "cut".
These rules use cutadapt to place in-line barcodes/UMIs in the read name and to trim off barcodes/UMIs/adapter sequences from the reads.
The published easyCLIP L5 and L3 are:

```txt
# Example L5, with barcode CGATGT, 5' to 3':
CTTGTGTGTGTGAAGCGAAGGGTA CGATGT NNNNNNN
# Example L3, with barcode TCA (reverse complement of TGA), 5' to 3':
TGA AGATCGGAAGAGCGGTTCAGAAAAAAAAAAAAAAAAAAAAAAAA
```

The 6 bp barcode 7 bp UMI in L5 sum to 13 bp and the L3 adapter has a 3 bp barcode.
This is represented in the config.yaml file with the entries:

```yaml
cutadapt_r1_adapter_seq_to_trim: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT
cutadapt_r2_adapter_seq_to_trim: TACCCTTCGCTTCACACACAAGGGGAAAGAGTGTAGATCTCGGTGGTCGC
L5_inline: BBBBBBNNNNNNN
L3_inline: BBB
```

The L5_inline entry represents the barcode (B) and UMI (N) in the order it would be sequenced in read 1 (that is, 6 barcode bases and then 7 UMI bases).
The L3_inline entry is the same, but in the order it would be sequenced from read 2.

The result of the above settings is that the parameters `-u 13 -U 3` would be passed to cutadapt will remove the 13 bp barcode/UMI from L5, and the 3 bp barcode from L3.
Then, the parameter `--rename={id}__{r1.cut_prefix}-{r2.cut_prefix}` passed to cutadapt will move the 13 bp from read 1 (r1.cut_prefix) and the 3 bp from read 2 (r2.cut_prefix) to the read name.

#### Sample sheets and fastq inputs

A sample sheet (for example, testdata/samples.txt) has the format:
```txt
Experiment	Gene	Replicate	L3_BC	L5_BC	R1_fastq	R2_fastq
Exp91	A1CF	1	TCA	CTGATC	All_R1.fastq.gz	All_R2.fastq.gz
Exp91	A1CF-E34K	1	AGT	CTGATC	All_R1.fastq.gz	All_R2.fastq.gz
```

All of these columns are required.
The L3_BC and L5_BC are the inline barcodes found in read 2 and read 1, respectively.
They should be written in the samples.txt file as they will appear in the raw read sequences.

In the case where PCR indexing has been used, samples must be demultiplexed from the index reads (not the inline barcodes) before running this workflow. The fastq files for each separate index set must be given in the R1_fastq/R2_fastq columns in the samples file.

The R1_fastq/R2_fastq columns in the samples.txt file must be just the file basename and the fastqs must be found in the folder set in the config file:
```yaml
Fastq_folder: testdata/fastq/raw
```
Samples will be assumed to exist in {Fastq_folder}/{R1_fastq} and {Fastq_folder}/{R2_fastq}.

The paired-end reads for each PCR index set (R1_fastq/R2_fastq pair) must have a common prefix followed by R1.fastq.gz or R2.fastq.gz.
That is, the must match the pattern {prefix}R1.fastq.gz and {prefix}R2.fastq.gz, with a common prefix.

#### Run workflow

Running the snakemake workflow with no arguments will run it on the testdata (placed in testdata/, with the config file config.yaml and the samples file testdata/samples.txt):

```
snakemake -j <cpus> --use-conda 
```

The workflow is split into three parts, in order:

```bash
snakemake -j <cpus> --use-conda -s snake.processing.py  # Raw fastq.gz to bam/bw/ect.
snakemake -j <cpus> --use-conda -s snake.clipper.py  # Run clipper to call peaks.
snakemake -j <cpus> --use-conda -s snake.analysis.py  # Analysis.
```
The required inputs are a star genome index to map to the config.yaml file, the genomic fasta for that build, the samples.txt file, and raw fastq.gz files for read1 and read2.

The paths to input files are set by editing config.yaml.
The snakefile can also be set to run a different config file by calling snakemake with a --configfile parameter.


#### Authors
dfporter, rmgarg

