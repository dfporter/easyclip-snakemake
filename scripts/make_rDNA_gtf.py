import pandas, re, csv, sys

# Convert the gff downloaded from ncbi that annotates U13369 to the required
# format and subset/expand the annotations.
    
def edit_gff(input_fname, output_fname):
    
    # Read the gff downloaded from ncbi that annotates U13369.
    df = pandas.read_csv(input_fname, sep='\t', comment='#', header=None)
    
    # Filter out annotations we don't want.
    df['gbkey'] = [re.search('gbkey=([^;]+)', x).group(1) for x in df[8]]
    df = df[df['gbkey']!='variation']
    df = df[df['gbkey']!='repeat_region']
    df = df[df['gbkey']!='Src']
    df = df[df[2]!='sequence_feature']
    df = df[df[2]!='sequence_conflict']
    df = df[df[2]!='protein_binding_site']
    df = df[['similar to' not in x for x in df[8]]]
    df = df[['site of initiation' not in x for x in df[8]]]
    del df['gbkey']

    # Function to convert the info line to the expected format.
    def parse_info(info: str) -> str:
        print(info)
        name = re.search('product=([^;]+)', info)
        if name is not None:
            name = name.group(1)
            name = re.sub(' ', '_', name)
            name = re.sub("'", '', name)
            name = name.strip('"')
            info = f'gene_id "{name}"; transcript_id "{name}"; gene_type "rRNA"; gene_name "{name}"; exon_number 1;"'
            return info
        name = re.search('Note=([^;]+)', info)
        if name is not None:
            name = name.group(1)
            name = re.sub(' ', '_', name).strip('"')
            info = f'gene_id "{name}"; transcript_id "{name}"; gene_type "rRNA"; gene_name "{name}"; exon_number 1;"'
            return info

        return info
    
    # Convert the info line to the expected format.
    df[8] = [parse_info(x) for x in df[8]]
    
    # Take the transcript and rRNA annotations only and make gene/transcript/exon annotations for all.
    genes = df.loc[[x=='transcript' or x=='rRNA' for x in df[2]],:].copy()

    # Gene annotations - splitting the locus into separate genes.
    genes[2] = 'gene'
    name = 'rDNA_downstream_region'
    
    # Add a line to cover the non-transcribed, downstream region of the rDNA locus.
    genes = genes.append(
        {0: 'U13369.1', 1: 'Genbank', 2: 'gene', 3: 13351, 4: 42999, 5: '.', 6: '+', 7: '.',
         8:f'''gene_id "{name}"; transcript_id "{name}"; gene_type "rRNA"; gene_name "{name}; exon_number 1;''',
        "gbkey": "misc_RNA"}, ignore_index=True)
    
    # Add identical lines for transcripts and exons.
    trxpt = genes.copy()
    trxpt[2] = 'transcript'
    exon = genes.copy()
    exon[2] = 'exon'
    
    # Reassemble the dataframe.
    df = pandas.concat([genes, trxpt, exon])
    
    # Write in gtf format.
    df.to_csv(output_fname, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
    
    
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("USAGE: python make_rDNA_genome input.gff output.gtf")
    edit_gff(sys.argv[1], sys.argv[2])
