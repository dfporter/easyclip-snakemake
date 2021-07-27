"""
Building a gffutils database object.

Original document header:
    https://www.biostars.org/p/152517/
    Example of how to work with Ensembl release 81 GTF files, which:
        1) already have genes and transcripts included
        2) have unique IDs for genes, transcripts, and exons in the corresponding
           "<featuretype>_id" attribute
        3) do not have unique IDs for CDS, stop_codon, start_codon, UTR.
    See background info at on database IDs at:
        https://pythonhosted.org/gffutils/database-ids.html
    GTF file from
    ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz

For humans, the GTF file was obtained from:
http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz   

"""

import gffutils, sys, os, subprocess

def first_n_features(data, n=5000):
    """
    Useful for testing: only use the first `n` features of source data
    """
    for i, feature in enumerate(gffutils.iterators.DataIterator(data)):
        if i > n:
            break
        yield feature


# Note: this function is optional; if you don't want these IDs then comment out 
# the lines at [1] below
def subfeature_handler(f):
    """
    Given a gffutils.Feature object (which does not yet have its ID assigned),
    figure out what its ID should be.
    This is intended to be used for CDS, UTR, start_codon, and stop_codon
    features in the Ensembl release 81 GTF files.  I figured a reasonable
    unique ID would consist of the parent transcript and the feature type,
    followed by an autoincrementing number.
    See https://pythonhosted.org/gffutils/database-ids.html#id-spec for
    details and other options.
    """
    return ''.join(
        ['autoincrement:',
         f.attributes['transcript_id'][0],
         '_',
         f.featuretype])


# gffutils can spend a lot of time trying to decide on a unique ID for each
# feature. So we have to give it hints of where to look in the attributes.
#
# We also tell it to use our subfeature_handler function for featuretypes with
# no unique IDs.
id_spec = {
    'exon': 'exon_id',
    'gene': 'gene_id',
    'transcript': 'transcript_id',

    # [1] These aren't needed for speed, but they do give nicer IDs.
    'CDS': [subfeature_handler],
    'stop_codon': [subfeature_handler],
    'start_codon': [subfeature_handler],
    'UTR':  [subfeature_handler],
}


def longest_transcript(db, gff='data/processed/longest_txpt_per_gene.gff'):
    """Write the gff file, given a gffutils database, using only the longest transcript per gene.
    """
    gene_to_name_to_len = {}
    for transc in db.features_of_type('transcript'):
        #if transc.seqid != 'chr22':
        #    continue        
        gene_to_name_to_len.setdefault(transc.attributes.get('gene_name')[0], {})
        gene_to_name_to_len[transc.attributes.get('gene_name')[0]][transc.id] = transc.stop - transc.start
        
    writer = gffutils.gffwriter.GFFWriter(gff)
    for gene in gene_to_name_to_len:
        longest_txpt = sorted(gene_to_name_to_len[gene], key=lambda x: gene_to_name_to_len[gene][x])[-1]
        writer.write_mRNA_children(db, longest_txpt)
    writer.close()

def create_a_db(gtf_filename, output_filename='data/processed/features.db'):
    # Build the features.db object from the gtf.
    os.makedirs('assets/reference/', exist_ok=True)
    os.makedirs('data/processed/', exist_ok=True)
    db = gffutils.create_db(
        gtf_filename,
        output_filename,

        # Since Ensembl GTF files now come with genes and transcripts already in
        # the file, we don't want to spend the time to infer them (which we would
        # need to do in an on-spec GTF file)
        disable_infer_genes=True,
        disable_infer_transcripts=True,

        # Here's where we provide our custom id spec
        id_spec=id_spec,

        # "create_unique" runs a lot faster than "merge"
        # See https://pythonhosted.org/gffutils/database-ids.html#merge-strategy
        # for details.
        merge_strategy='create_unique',
        verbose=True,
        force=True,
    )

def create_db(gtf_filename, output_filename='data/processed/features.db'):
    # 1. Download a gtf of genome annotations.

    # Download the annotation GTF and unzip:
    # wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
    # gunzip Homo_sapiens.GRCh38.104.chr.gtf.gz

    # 2. Build the features.db object from the gtf.
   # gtf_filename = config['feature_gtf'] #"Homo_sapiens.GRCh38.104.chr.gtf"
    os.makedirs('assets/reference/', exist_ok=True)
    os.makedirs('data/processed/', exist_ok=True)
    db = gffutils.create_db(
        gtf_filename,
        'data/processed/features.no_introns.db',

        # Since Ensembl GTF files now come with genes and transcripts already in
        # the file, we don't want to spend the time to infer them (which we would
        # need to do in an on-spec GTF file)
        disable_infer_genes=True,
        disable_infer_transcripts=True,

        # Here's where we provide our custom id spec
        id_spec=id_spec,

        # "create_unique" runs a lot faster than "merge"
        # See https://pythonhosted.org/gffutils/database-ids.html#merge-strategy
        # for details.
        merge_strategy='create_unique',
        verbose=True,
        force=True,
    )

    # 3. Create introns.
    # This takes a long time.

    intronDb = db.create_introns(exon_featuretype='exon', grandparent_featuretype='transcript')

    # Using the recommended strategy of db.update() fails in my hands, so instead we use a work-around.
    # We write a new gtf on just the introns.
    with open('data/processed/introns.gtf', 'w') as f:
        for li in intronDb:
            f.write(str(li) + '\n')

    # Then concatenate the original gtf and the introns:
    cmd = f"cat {gtf_filename} data/processed/introns.gtf > data/processed/annotation_with_introns.gtf"
    print(cmd)
    os.system(cmd)
    
    # Create a new db, this time with introns.
    gtf_filename = "data/processed/annotation_with_introns.gtf"
    db = gffutils.create_db(
        gtf_filename,
        output_filename,

        # Since Ensembl GTF files now come with genes and transcripts already in
        # the file, we don't want to spend the time to infer them (which we would
        # need to do in an on-spec GTF file)
        disable_infer_genes=True,
        disable_infer_transcripts=True,

        # Here's where we provide our custom id spec
        id_spec=id_spec,

        # "create_unique" runs a lot faster than "merge"
        # See https://pythonhosted.org/gffutils/database-ids.html#merge-strategy
        # for details.
        merge_strategy='create_unique',
        verbose=True,
        force=True,
    )

    # Finally, we write a gff of the longest transcript per gene (so one transcript
    # per gene), then use it to make a db file of just the longest txpt per gene.
    # This will also have introns.
    longest_transcript(db, gff='data/processed/longest_txpt_per_gene.gff')
    
    db = gffutils.create_db(
        'data/processed/longest_txpt_per_gene.gff',
        'data/processed/longest_txpt_per_gene.db',

        # Since Ensembl GTF files now come with genes and transcripts already in
        # the file, we don't want to spend the time to infer them (which we would
        # need to do in an on-spec GTF file)
        disable_infer_genes=True,
        disable_infer_transcripts=True,

        # Here's where we provide our custom id spec
        id_spec=id_spec,

        # "create_unique" runs a lot faster than "merge"
        # See https://pythonhosted.org/gffutils/database-ids.html#merge-strategy
        # for details.
        merge_strategy='create_unique',
        verbose=True,
        force=True,
    )

def export_bed12(
    db_fname='data/processed/longest_txpt_per_gene.db',
    bed12_fname=None):
    
    db = gffutils.FeatureDB(db_fname)
    
    if bed12_fname is None:
        bed12_fname = db_fname + '.bed'
    
    with open(bed12_fname, 'w') as outf:
        for transc in db.features_of_type('transcript'):
            outf.write(db.bed12(transc) + '\n')


if __name__ == '__main__':
    
    final_output = 'data/processed/features.db'
    if len(sys.argv) > 2:
        final_output = sys.argv[2]
    
    create_db(sys.argv[1], final_output)

    
