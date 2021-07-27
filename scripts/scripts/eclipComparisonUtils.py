import pandas, os, matplotlib, HTSeq, pysam, pyBigWig, re, scipy, importlib
import random
import numpy as np

import sklearn
from sklearn.neighbors import KernelDensity

import scipy.stats
import scipy.signal
from pprint import pprint as pp
from typing import List, Mapping, Union
import matplotlib.pyplot as plt
from liftover import get_lifter

input_dir = '/Users/dp/pma/rbfox2_eCLIP_comparison/data/'

unusable_intervals = []
with open(f'{input_dir}/unusable_intervals.txt', 'r') as f:
    for li in f:
        s = li.split('\t')
        iv = [s[0], int(s[1]), int(s[2])]
        unusable_intervals.append(iv)

# Get the total read numbers in the easyCLIP datasets.
total_read_numbers_easyCLIP_fname = '/Users/dp/pma/dataAndScripts/clip/meta/data/total_read_numbers.txt'
total_read_numbers = pandas.read_csv(total_read_numbers_easyCLIP_fname, sep='\t')
rbfox_read_numbers = total_read_numbers.loc[[bool('fox2' in x) for x in total_read_numbers.Dataset], :]
_rbfox_read_numbers = dict(zip(rbfox_read_numbers.Dataset, rbfox_read_numbers['Total read number']))

rbfox_read_numbers = {}
for k, v in _rbfox_read_numbers.items():
    rbfox_read_numbers[k] = v
    rbfox_read_numbers[k + '_+.bam'] = v
    rbfox_read_numbers[k + '_-.bam'] = v


# Read easyCLIP peaks.
class easyCLIPpeaks():

    def __init__(self, apply_p_cutoff=0.01):

        self.peaksEasy = self.read_easyCLIP_peaks(apply_p_cutoff=apply_p_cutoff)

    def read_easyCLIP_peaks(self, apply_p_cutoff=0.01):
        def is_ok_chrom(_str):
            if type(_str) != type(''):
                return False
            c = _str.split(':')[0]
            if re.match('\d+', c) is not None:
                return True
            if c in ['Y', 'M', 'chrY', 'chrM']:
                return True
            return False


        peaksEasy = pandas.read_excel('/Users/dp/pma/dataAndScripts/clip/meta/tables/File 5 Peak locations.xlsx')
        peaksEasy = peaksEasy.loc[[not pandas.isna(x) for x in peaksEasy.Rbfox2],['Gene name', 'Rbfox2']]
        peaksEasy = peaksEasy.loc[[is_ok_chrom(_) for _ in peaksEasy.Rbfox2], :]

        pvals = pandas.read_excel('/Users/dp/pma/dataAndScripts/clip/meta/tables/File 4 Significance values.xlsx', sheet_name='Per million reads')

        pvals = pvals.loc[:,['Gene name', 'Rbfox2']]
        pvalsd = dict(zip(pvals['Gene name'], pvals['Rbfox2']))
        genes = set([x.split('::')[0] for x in pvals['Gene name']])
        _p = {}
        for gene in genes:
            try:
                intron = pvalsd[gene + '::intron']
                exon = pvalsd[gene + '::exon']
            except:
                print(gene, "?", end=', ')
                _p[gene] = 1
                continue
            _p[gene] = min([intron, exon])
        peaksEasy['Pval'] = [_p[gene] for gene in peaksEasy['Gene name']]

        if apply_p_cutoff:
            #allEasyPeaks = peaksEasy.copy()
            peaksEasy = peaksEasy.loc[peaksEasy['Pval']<apply_p_cutoff, :]
            return peaksEasy
        else:
            return peaksEasy

    def random_easy_peak_iv(self, expand=1000, N_peaks=1):

        peaksEasy = self.peaksEasy

        def from_row(loc):
            chrom, pos = loc.split(':')
            pos, strand = pos.split('/')
            pos = int(pos)

            if pos > expand:
                loc = [chrom, pos-expand, pos+expand, strand]
            else:
                loc = [chrom, 0, pos+expand, strand]

            if 'chr' not in loc[0]:
                loc[0] = 'chr' + loc[0]
                
            return loc
        
        if N_peaks > 1:
            s = list(range(peaksEasy.shape[0]))
            random.shuffle(s)
            locs = [from_row(peaksEasy.iloc[n]['Rbfox2']) for n in s[:N_peaks]]
            return locs
        else:
            loc = peaksEasy.iloc[scipy.stats.randint.rvs(0, peaksEasy.shape[0])]['Rbfox2']
            return from_row(loc)


def get_chrom_lengths(bam_fname: str) -> Mapping[str, int]:
    samfile = pysam.AlignmentFile(bam_fname)
    header = str(samfile.header)
    chrom_lengths = {}
    for li in header.split('\n'):
        if li[:3] == '@SQ':
            s = li.split('\t')
            chrom_lengths[s[1].split(':')[1]] = int(s[2].split(':')[1])
    samfile.close()
    return chrom_lengths


        
def get_random_genomic_intervals(bam_fname: str, interval_len=1E4, N_intervals=1):
    # Get chromosome lengths.
    chrom_lengths = get_chrom_lengths(bam_fname)
    normal_chroms = [f'chr{n}' for n in range(1, 23)] + ['chrM', 'chrY', 'chrX']

    chrom_lengths = {x:y for x,y in chrom_lengths.items() if x in normal_chroms}
    total_genome_size = sum([x for x in chrom_lengths.values()])
    fractions = {chr_name:_len/total_genome_size for chr_name,_len in chrom_lengths.items()}

    intervals = []
    while len(intervals) < N_intervals:
        c = np.random.choice(list(fractions.keys()), size=1000, replace=True, p=list(fractions.values()))[0]
        start = np.random.randint(0, high=chrom_lengths[c]-interval_len-1)
        end = start + interval_len
        strand = np.random.choice(['-', '+'])
        
        # Check if this is allowed.
        usable = True
        for iv in unusable_intervals:
            if iv[0] == c:
                if (iv[1] <= start <= iv[2]) or (iv[1] <= end <= iv[2]):
                    usable = False
        if usable:
            intervals.append([c, int(start), int(end), strand])
            
    return intervals


def get_coverage_in_an_interval_for_bam_file(bam_fname, interval, norm_to_per_million=True):
    #count_coverage(self, contig, start=None, stop=None, region=None, quality_threshold=15, read_callback='all', reference=None, end=None)
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    # Get np.array of coverage depth across interval.
    b = np.sum(samfile.count_coverage(*interval[:-1]), axis=0)
    if norm_to_per_million:
        b = b * 1000000/rbfox_read_numbers[os.path.basename(bam_fname)]
    samfile.close()
    return b


def get_coverage_in_an_interval_for_bigwig_file(bw_fname, interval):
    bw = pyBigWig.open(bw_fname)
    return np.nan_to_num(bw.values(*interval[:-1]))


def coverage(fname, interval, norm_to_per_million=True):
    
    if interval[3] == '-':
        fname = re.sub('_\+\.', '_-.', fname)
    else:
        fname = re.sub('_\-\.', '_+.', fname)

    if fname.split('.')[-1] == 'bam':
        return get_coverage_in_an_interval_for_bam_file(fname, interval, norm_to_per_million=norm_to_per_million)
    return get_coverage_in_an_interval_for_bigwig_file(fname, interval)

def smooth(y, box_pts):
    box_pts = int(box_pts)
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth