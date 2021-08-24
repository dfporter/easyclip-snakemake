rule correlation_heatmap:
    input:
        "outs/counts/rt_stop_bams/counts.txt"
    output:
        pearson = "outs/figs/heatmap_reads_per_transcript_spearman_all_transcripts.pdf",
    run:
        figs_dir = os.path.dirname(str(output.pearson))
        shell(f"python scripts/reads_per_gene_correlations_heatmap.py {str(input)} {figs_dir}")