import pyBigWig, re, gffutils, os, sys, collections, pandas, glob, HTSeq, random
import numpy as np
import pysam


def set_bounds(txpt_id, db):
    CDS_bounds = []
    for feat in db.children(txpt_id, featuretype='CDS', order_by='start'):
        if len(CDS_bounds) == 0:
            CDS_bounds = [feat.start, feat.end]
        else:
            CDS_bounds = [ min([feat.start, CDS_bounds[0]]), max([feat.end, CDS_bounds[1]]) ]    
    return CDS_bounds


def find_exons_for_genename(genename, db):
    five_regions = []
    cds_regions = []
    three_regions = []
    cds = set_bounds(genename, db)
    if len(cds) < 2:
        print(f"No CDS found for {genename}.")
        return {'5UTR': five_regions, 'CDS': cds_regions, '3UTR': three_regions}
    for feat in db.children(genename, featuretype='exon'):
        iv = ['', feat.start, feat.end]
        if feat.strand == '+' and iv[1] >= cds[1]: 
            three_regions.append(feat)
        elif feat.strand == '+' and iv[2] <= cds[0]:
            five_regions.append(feat)
        elif feat.strand == '-' and iv[1] <= cds[0]: 
            three_regions.append(feat)
        elif feat.strand == '-' and iv[2] >= cds[1]:
            five_regions.append(feat)
        elif (cds[0] <= iv[1] <= cds[1]) and (cds[0] <= iv[2] <= cds[1]):
            cds_regions.append(feat) 
    return {'5UTR': five_regions, 'CDS': cds_regions, '3UTR': three_regions}


def get_coverage_in_an_interval_for_bigwig_file(bw, interval):
    
    #print(bw_fname, '>>>', *interval[:3])
    try:
        return np.nan_to_num(bw.values(*interval[:3]))
    except:
        return 0.

    
def get_coverage_in_an_interval_for_bam_file(bam_fname, interval):
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    # Get np.array of coverage depth across interval.
    b = np.sum(samfile.count_coverage(*interval[:-1]), axis=0)
    samfile.close()
    return b


def for_every_transcript(db):
    for transc in db.features_of_type('transcript'):
        #if transc.seqid != 'chr22':
        #    continue
        yield [transc.seqid, transc.start, transc.end, transc.strand, 
               transc.attributes.get('transcript_id')[0], transc.attributes.get('gene_name')[0]]


def coverage(fname, interval):
    if interval[3] == '-':
        fname = re.sub('\.\+\.', '.-.', fname)
        fname = re.sub('plus', 'minus', fname)
    else:
        fname = re.sub('\.\-\.', '.+.', fname)
        fname = re.sub('minus', 'plus', fname)
        
    if fname.split('.')[-1] == 'bam':
        return get_coverage_in_an_interval_for_bam_file(fname, interval)
    
    return get_coverage_in_an_interval_for_bigwig_file(fname, interval)


def get_regions(txpt, db):
    # Broken.
    regions = []
    left = txpt.start
    right = txpt.end
    for feat in db.region(txpt):
        if feat.featuretype != 'exon' and feat.featuretype != 'intron':
            continue
        if feat.attributes['transcript_id'][0] == txpt.id:
            continue
        if feat.start > left:
            regions.append([left, feat.start])
            left = min([txpt.end, feat.end])
    if left < right:
        regions.append([left, right])
    return regions


def coverage_across_every_transcript(db, bam_list):
    results = collections.defaultdict(dict)
    #print(f"Iterating through {db.count_features_of_type('transcript')} transcripts...")
          
    for n, info in enumerate(for_every_transcript(db)):
          
        for bw in bam_list:
            results[info[4]][bw] = np.sum(coverage(bw, info[:4]))
          
        not (n % 1000) and print(n, results[info[4]])
        results[info[4]]['gene_name'] = info[5]
          
    results = pandas.DataFrame(results).T
    
    return results


def coverage_across_every_transcript_bw(db, gaos, bw_list):
    results = collections.defaultdict(dict)
    results['_ambiguous'] = {os.path.basename(bw_name): 0 for bw_name in bw_list}
    peak_loc = collections.defaultdict(dict)
    
    print(f"Iterating through " + str(db.count_features_of_type('transcript')) + " transcripts...")
    bws = {
        '+': {os.path.basename(bw_name):pyBigWig.open(bw_name) for bw_name in bw_list},
        '-': {os.path.basename(bw_name):pyBigWig.open(re.sub('\.\+\.', '.-.', bw_name)) for bw_name in bw_list},
    }
    
    for n, txpt in enumerate(db.features_of_type('transcript')):
        txpt_id = txpt.attributes.get('transcript_id')[0]
        
        txpt_iv = HTSeq.GenomicInterval(txpt.seqid, txpt.start, txpt.end, txpt.strand)
        
        for region_iv, genes in gaos[txpt_iv].steps():
            
            if region_iv.end <= region_iv.start:
                continue  # No 0 length intervals.
            
            # genes = txpt_id::intron or txpt_id::exon.
            if len(genes) == 1:
                
                gene = list(genes)[0]
                
                for bw_name, bw in bws[txpt.strand].items():
                    
                    results[gene].setdefault(bw_name, 0)  # Initialize to zero if needed.
                    
                    try:
                        results[gene][bw_name] += np.nansum(bw.values(txpt.seqid, region_iv.start, region_iv.end), axis=0)
                    except:
                        results[gene][bw_name] += 0.  # Chromosome not found, we assume.
                        
            if len(genes) > 1:
                # If one gene is an intron and the other an exon, take the exon. Otherwise it's ambiguous.
                
                _type_counts = collections.defaultdict(set)
                
                for gene in genes:
                    
                    _type_name = gene.split('::')
                    
                    if len(_type_name) < 2:
                        continue
                    
                    _type_counts[_type_name[1]].add(gene)
                    
                if ('exon' in _type_counts) and (len(_type_counts['exon']) == 1):
                    
                    results[list(_type_counts['exon'])[0]].setdefault(bw_name, 0)  # Initialize to zero if needed.
                    
                    try:
                        results[list(_type_counts['exon'])[0]][bw_name] += np.nansum(
                            bw.values(txpt.seqid, region_iv.start, region_iv.end), axis=0)
                    except:
                        print(f"Error obtaining values for {(txpt.seqid, region_iv.start, region_iv.end)}" + \
                              f" in {bw_name}")
                    
                else:
                    results['_ambiguous'].setdefault(bw_name, 0)  # Initialize to zero if needed.
                    
                    try:
                        results['_ambiguous'][bw_name] += np.nansum(
                            bw.values(txpt.seqid, region_iv.start, region_iv.end), axis=0)
                    except:
                        print(f"Error obtaining values for {(txpt.seqid, region_iv.start, region_iv.end)}" + \
                              f" in {bw_name}")
                        
        # End for-each-region-in-txpt block.
        
        not (n % 100) and print(n, gene, results[gene])  # Once in a while print progress.
        
        # Always write the gene name.        
        results[txpt_id + '::exon']['gene_name'] = txpt.attributes.get('gene_name')[0]
        if txpt_id + '::intron' in results:
            results[txpt_id + '::intron']['gene_name'] = txpt.attributes.get('gene_name')[0]
    
    # End for-each-transcript block.
    
    # Close all the bigwig files.
    for strand in bws:
        for bw in bws[strand].values():
            bw.close()
            
    # Turn the {txpt::exon -> {bigwig -> coverage }} dict into a dataframe to return it.
    results = pandas.DataFrame(results).T
    
    # Sort by the sum of numeric columns.
    results['sum'] = results.loc[:, [x for x in results.columns if 'bigwig' in x]].apply(np.nansum, axis=1)
    results = results.sort_values('sum', ascending=False)
    del results['sum']
    
    return results

def load_dbs(db_file):
    print("Loading FeatureDB...")
    #db_file = 'data/processed/features.db'
    db = gffutils.FeatureDB(db_file, keep_order=True)
    
    print("Making HTSeq.GenomicArrayOfSets...")
    # Create a genomic array of sets.
    gaos = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    
    n = 0
    for feat in db.features_of_type(['exon', 'intron']):
        
        if feat.end <= feat.start:
            continue  # No weird stuff.
            
        gaos[HTSeq.GenomicInterval(
            feat.seqid, feat.start, feat.end, feat.strand)] += '::'.join([
            feat.attributes.get('transcript_id')[0], feat.featuretype])
        
        n += 1
        if not n % 10000:
            print(n, feat)
    return db, gaos


if __name__ == '__main__':
    """Debugging:
input_dir = "rb/bigwig/3prime/"
db_file = "data/processed/repeats_and_longest_txpt_per_gene.db"
output_fname = "rb/outs/counts/bigwig_3prime_counts_transcripts.txt"
    
bam_list = glob.glob(input_dir + '/*.bam')

db, gaos = load_dbs(db_file)

import scripts, importlib
import scripts.reads_to_genes
importlib.reload(scripts.reads_to_genes)
from scripts.reads_to_genes import coverage_across_every_transcript_bw

bw_list = glob.glob(input_dir + '*.+.bigwig')
df = coverage_across_every_transcript_bw(db, gaos, bw_list)

    """
    input_dir, db_file, output_fname = sys.argv[1:]
    
    os.makedirs(os.path.dirname(output_fname), exist_ok=True)
    
    bam_list = glob.glob(input_dir + '/*.bam')
    
    db, gaos = load_dbs(db_file)
    
    if len(bam_list):
        df = coverage_across_every_transcript(db, bam_list)
        
    else:
        bw_list = glob.glob(input_dir + '*.+.bigwig')
        df = coverage_across_every_transcript_bw(db, gaos, bw_list)

        if len(bw_list) == 0:
            bam_list = glob.glob(input_dir + '*.+.bw')
            df = coverage_across_every_transcript_bw(db, gaos, bw_list)

        if len(bw_list) == 0:
            print(f"Could not find any files with .bam/.bigwig/.bw extensions in {input_dir}. Quitting.")
            sys.exit()
    
    df.to_csv(output_fname, sep='\t')
