# Readme

This software is a simplification, formalization, and streamlining of easyCLIP read processing, mapping, and initial analysis.
It is intended to be easily usable.

The code used for the paper https://pubmed.ncbi.nlm.nih.gov/33692367/, warts and all, is placed in https://github.com/dfporter/easyclip.

#### Install workflow

Clone this repositiory into a directory

```
$ git clone git@github.com:dfporter/easyclip-v2.git
```

Create and activate the conda environment:

```
$ conda env create -f=envs/conda.yaml -n easyclip-env
$ conda activate easyclip-env
```

#### Run workflow

Running the snakemake workflow with no arguments will run it on the testdata (placed in testdata/, with the config file config.yaml and the samples file testdata/samples.txt):

```
$ snakemake -j <cpus> --use-conda 
```

The workflow is split into three parts, in order:

```bash
$ snakemake -j <cpus> --use-conda -s snake.processing.py  # Raw fastq.gz to bam/bw/ect.
$ snakemake -j <cpus> --use-conda -s snake.clipper.py  # Run clipper to call peaks.
$ snakemake -j <cpus> --use-conda -s snake.analysis.py  # Analysis.
```
The required inputs are a star genome index to map to the config.yaml file, the genomic fasta for that build, the samples.txt file, and raw fastq.gz files for read1 and read2.

The paths to input files are set by editing config.yaml.
The snakefile can also be set to run a different config file by changing the variable configfile at the top of the workflow.

The format of the samples.txt file is given my demonstration in testdata/samples.txt. It must include an Experiment, Gene, Replicate, L3_BC and L5_BC column (and that is all).

#### Authors
dfporter, rmgarg

