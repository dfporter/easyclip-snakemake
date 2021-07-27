import pandas, os, sys, re, glob

"""From the wiki on the git:
Clipper outputs a bed8 file:

chromosome, genomic_start, genomic_stop, cluster_name, min_pval, strand, thick_start, thick_stop

cluster_name contains the gene name, the number of reads in the cluster and a counter to give the 
cluster a unique identifier.

The genomic start and genomic start are the positions above the FDR cutoff that we used to call the 
peak. Any reads with their "middle" in this region are counted towards the cluster's read count.
Thick start and thick stop specify the apex of the fitted curve in the peak.
"""

# For empty bed file inputs, write out a single line 'peak' to prevent
# crashes from empty peaks files and to give snakemake an output file.
columns = ['chrom', 'start', 'end', 'name', 'p', 'strand', 'apex_start', 'apex_end']
_d = {'chrom': 'chr1', 'start': 1000, 'end': 1020, 'name': 'placeholder.0_0', 'p': 1,
     'strand': '-', 'apex_start': 1005, 'apex_end': 1010}
placeholder_df = pandas.DataFrame.from_dict({0:_d}, orient='index')
placeholder_df = placeholder_df.loc[:, columns]


def filter_clipper(input_bed_fname):
    
    output_dir = os.path.dirname(input_bed_fname)
    os.makedirs(output_dir + '/filtered/', exist_ok=True)
    cutoffs = {}
    
    # Bed file might not exist, or might be empty, so we just replace with a placeholder
    # df in order to not have empty files. Not a great solution.
    try:
        df = pandas.read_csv(input_bed_fname, sep='\t', header=None)
        df.columns = ['chrom', 'start', 'end', 'name', 'p', 'strand', 'apex_start', 'apex_end']
    except:
        df = placeholder_df.copy()
        
    # name column example: ENSG00000227232.5_0_4 
    df['# reads'] = [int(x.split('_')[-1]) for x in df['name']]
    
    print(f"Inital peak number: {df.shape[0]}.")
    cutoffs = {
        '1': {'p': 1E-2, '# reads': 2}, '2': {'p': 1E-9, '# reads': 5}, '3': {'p': 1E-20, '# reads': 10}}
    
    for cutoff_name, cutoff in cutoffs.items():
        
        # Filter df.
        sub = df.loc[[x<cutoff['p'] for x in df['p']],:]
        sub = sub.loc[[x>=cutoff['# reads'] for x in sub['# reads']],:]
        
        print(f"{cutoff} -> peak number: {sub.shape[0]}")
        
        if sub.shape[0] == 0:  # No peaks.
            placeholder_df.to_csv(f"{output_dir}/filtered/" + \
                os.path.splitext(os.path.basename(input_bed_fname))[0] + f".{cutoff_name}.bed",
                sep='\t', index=False, header=False)
        
        else:  # There is at least one peak.
            # Write df.
            sub.to_csv(
                f"{output_dir}/filtered/" + \
                os.path.splitext(os.path.basename(input_bed_fname))[0] + f".{cutoff_name}.bed",
                sep='\t', index=False, header=False)
        
if __name__ == '__main__':
    filter_clipper(sys.argv[1])