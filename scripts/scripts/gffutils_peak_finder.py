
db = gffutils.FeatureDB('data/processed/longest_txpt_per_gene.db')

genes = list(db.features_of_type('transcript'))[:1000]
#gene = [x for x in db.children('ENST00000509176.5') if x.featuretype=='transcript']
def get_regions(txpt):
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


#to_iv = lambda txpt,r: (txpt.seqid, r[0], r[1], txpt.strand)

bamfnames = ['m0/sams/split/Exp91_FUBP1_GCTCAT_TCA.bam']
bamfiles = {bam_fname:pysam.AlignmentFile(bam_fname, "rb" ) for bam_fname in bamfnames}
signals, peak_loc = {}, {}
for txpt in db.features_of_type('transcript'):
    regions = get_regions(txpt)
    signals[txpt.id] = {}
    peak_loc[txpt.id] = {}
    if len(regions) == 0:
        continue
    for bam, bam_fh in bamfiles.items():
        signals[txpt.id][bam] = np.concatenate([
            np.sum(bam_fh.count_coverage(txpt.seqid, r[0], r[1]), axis=0) for r in regions])
        
        _loc = np.argmax(signals[txpt.id][bam])
        
        if signals[txpt.id][bam][_loc] > 2:
            print(txpt.id)
            peak_loc[txpt.id][bam] = np.argmax(signals[txpt.id][bam])
        else:
            peak_loc[txpt.id][bam] = np.nan

for bamfh in bamfiles.values():
    bamfh.close()
