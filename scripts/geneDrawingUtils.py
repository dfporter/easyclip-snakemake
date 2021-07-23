import random, importlib
import HTSeq
import numpy as np
import scripts.eclipComparisonUtils
importlib.reload(scripts.eclipComparisonUtils)
ecu = scripts.eclipComparisonUtils
#coverage = get_coverage_in_an_interval_for_bam_file
smooth = ecu.smooth

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



def plot_gene(iv, _ax, db, RNAs):

    for feat in db.region(region=(str(iv[0]), iv[1], iv[2], str(iv[3])), featuretype=['CDS']):
        left = max([iv[1], feat.start]) - iv[1]
        right = min([iv[2], feat.end]) - iv[1]
        lw = 3
        _ax.hlines(xmin=left, xmax=right, y=-1,  lw=lw, colors='k')
        xarr = np.linspace(left+50, right-50, num=min([5, int(right-left)/100]))
        marker = '>' if iv[3] == '+' else '<'
        _ax.plot(xarr, [-1]*len(xarr), lw=0, marker=marker, markeredgewidth=0, markersize=3, color='w')
        
    for feat in db.region(region=(str(iv[0]), iv[1], iv[2], str(iv[3])), featuretype=["UTR"]):
        left = max([iv[1], feat.start]) - iv[1]
        right = min([iv[2], feat.end]) - iv[1]
        lw = 1
        _ax.hlines(xmin=left, xmax=right, y=-1,  lw=lw, colors='k')
        
    wrote_gene = False
    for step, gene in RNAs.gaos[HTSeq.GenomicInterval(  str(iv[0]), iv[1], iv[2], str(iv[3])  )].steps():
        print(step, gene)
        left = step.start - iv[1]
        right = step.end - iv[1]

        if len(gene):
            
            if 'intron' not in list(gene)[0]:
                continue
            _gene = list(gene)[0].split('::')[0]

            lw = 2 if 'exon' in list(gene)[0] else 0.5
            #lw = 6 if exon_cat == 'CDS' else lw
            yoffset = 0 if 'exon' in list(gene)[0] else 0
            style = 'solid' if 'exon' in list(gene)[0] else 'dashed'
            
            color = 'k'
            
            _ax.hlines(xmin=left, xmax=right, y=-1+yoffset,  lw=lw, linestyle=style, colors=color)
            xarr = np.linspace(left+50, right-50, num=int(right-left)/200)
            marker = '>' if iv[3] == '+' else '<'
            #_ax.plot(xarr, [-1]*len(xarr), lw=0, marker=marker, markeredgewidth=0, markersize=3, color='k')            

            (not wrote_gene) and _ax.annotate(list(gene)[0].split('::')[0], xy=(0.05, 0.55), xycoords='axes fraction')
            wrote_gene = True

                
        else:
            _ax.hlines(xmin=left, xmax=right, y=-1, lw=0, colors='r')
    _ax.axis('off')

    if wrote_gene:
        return _gene
    return ''
