import numpy as np

def set_bounds(txpt_id, db):
    
    CDS_bounds = []
    
    for feat in db.children(txpt_id, featuretype='CDS', order_by='start'):
        if len(CDS_bounds) == 0:
            CDS_bounds = [feat.start, feat.end]
        else:
            CDS_bounds = [ min([feat.start, CDS_bounds[0]]), max([feat.end, CDS_bounds[1]]) ]
                
    return CDS_bounds

def assign_location(iv, db, gene_name=False):

    if type(iv) != type([]) and np.isnan(iv):
        return np.nan

    cds_bounds = {}
    by_txpt_id = {}
    for feat in list(db.region(region=(iv), completely_within=False)):
        att = feat.attributes

        if gene_name and (feat.attributes['gene_name'][0] != gene_name):
            continue
        if feat.attributes['gene_type'][0] in ['pseudogene', 'nonsense_mediated_decay']: continue
        if "transcript_id" not in feat.attributes: continue
        if feat.attributes['transcript_support_level'][0] != '1': continue

        txpt_id = att['transcript_id'][0]
        if txpt_id not in cds_bounds:
            cds_bounds[txpt_id] = set_bounds(txpt_id, db)
        
        by_txpt_id.setdefault(att['transcript_id'][0], {})
        by_txpt_id[feat.attributes['transcript_id'][0]][feat.featuretype] = (
            feat.start, feat.end, feat.strand, feat.attributes['gene_id'][0])

    votes = {'3UTR only': 0, '5UTR only': 0, '3UTR': 0, '5UTR': 0, 'CDS': 0, 'Unclear': 0}

    for txpt_id, vals in by_txpt_id.items():
        if txpt_id not in cds_bounds or len(cds_bounds[txpt_id]) < 2:
            continue

        cds = cds_bounds[txpt_id]

        if db[txpt_id].strand == '+' and iv[1] >= cds[1]: 
            votes['3UTR'] += 1
        elif db[txpt_id].strand == '+' and iv[2] <= cds[0]:
            votes['5UTR'] += 1
        elif db[txpt_id].strand == '-' and iv[1] <= cds[0]: 
            votes['3UTR'] += 1
        elif db[txpt_id].strand == '-' and iv[2] >= cds[1]:
            votes['5UTR'] += 1
        elif (cds[0] <= iv[1] <= cds[1]) and (cds[0] <= iv[2] <= cds[1]):
            votes['CDS'] += 1
        else:
            votes['Unclear'] += 1
    
    if sum(votes.values()) == 0:
        return 'Unclear'

    return sorted(votes.keys(), key=lambda x: -votes[x])[0]
