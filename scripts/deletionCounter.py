import pandas, os, sys, typing, matplotlib, HTSeq, glob, pysam, collections, subprocess, pyBigWig, re, scipy, importlib
import random, argparse
import seaborn as sns
import numpy as np

import scipy.stats
from pprint import pprint as pp
from typing import List, Mapping, Union

import matplotlib.pyplot as plt


class deletionCounter():
    
    def __init__(self, genomic_fasta=None):
        self.genomic_fasta = genomic_fasta
        
    def write_deletion_positions_to_bedgraph(self, bamfilenames, interval=None, outdir='/Users/dp/pma/rbfox2_eCLIP_comparison/cims/'):
        # NDEL1: chr17:8,431,758-8,472,580
        os.makedirs(outdir, exist_ok=True)
        
        iv = interval[:3]
        
        for bamfilename in bamfilenames:
            
            if interval[-1] == '-':
                samfile = pysam.AlignmentFile(re.sub('_\+\.bam', '_-.bam', bamfilename), "rb" )
            else:
                samfile = pysam.AlignmentFile(bamfilename, "rb" )
            
            ga = HTSeq.GenomicArray('auto', stranded=True)
            
            for r in samfile.fetch(*iv):  # For every read in bam.
                deletions, deleted_bases, insertions, inserted_bases = self.deletions_and_insertions_in_read(r)
                
                for pos, bases in zip(deletions, deleted_bases):
                    ga[HTSeq.GenomicInterval(iv[0], pos, pos+len(deleted_bases), interval[-1])] += 1
                    
            samfile.close()
            
            ga.write_bedgraph_file(f"{outdir}/" + os.path.basename(bamfilename).rstrip('.bam') + ".wig", strand=interval[-1])
            
            #arrs[bamfilename] = self._deletion_locations_for_intervals(bamfilename, intervals, outdir=outdir)

    def deletion_locations_for_intervals(self, bamfilenames, intervals, outdir='/Users/dp/pma/rbfox2_eCLIP_comparison/cims/'):

        if type(bamfilenames) == type([]):
            
            arrs = {}
            for bamfilename in bamfilenames:
                arrs[bamfilename] = self._deletion_locations_for_intervals(bamfilename, intervals, outdir=outdir)
            
            # Switch to be [iv][replicate] instead of [replicate][iv].
            _arrs = {}
            for iv in intervals:
                if len(iv) > 3:
                    iv = iv[:3]
                _arrs[tuple(iv)]  = {bamfilename:arrs[bamfilename][tuple(iv)] for bamfilename in arrs.keys()}
            
            return _arrs
        
        # Not a list of samfiles, just one.
        return self._deletion_locations_for_intervals(bamfilenames, intervals, outdir=outdir)
    
    def _deletion_locations_for_intervals(self, bamfilename, intervals, outdir='/Users/dp/pma/rbfox2_eCLIP_comparison/cims/'):
        
        deletion_locs = collections.defaultdict(list)
        deletion_arrs = {}
        for iv in intervals:
            arr = np.zeros(iv[2] - iv[1])
            
            if iv[-1] == '-':
                samfile = pysam.AlignmentFile(re.sub('_\+\.bam', '_-.bam', bamfilename), "rb" )
            else:
                samfile = pysam.AlignmentFile(bamfilename, "rb" )
            
            if len(iv) > 3:
                iv = iv[:3]
            
            for n, r in enumerate(samfile.fetch(*iv)):
                deletions, deleted_bases, insertions, inserted_bases = self.deletions_and_insertions_in_read(r)
                #print(deletions, deleted_bases)
                if len(deletions):
                    deletions = [x-iv[1] for x in deletions]
                    
    
                    deletion_locs[tuple(iv)].append(deletions)
                    for pos in deletions:
                        if pos < 0:
                        #    print(f"iv={iv}, pos={pos}, orig={pos+iv[1]}")
                            continue
                        if pos < len(arr):
                            arr[pos] += 1
                        #else:
                        #    print(f"Deletion position at {pos}, outside of len {len(arr)} iv {iv}")
            deletion_arrs[tuple(iv)] = arr
            
            samfile.close()

            #plt.plot(arr)
            #plt.show(); plt.clf(); plt.close()
            
        return deletion_arrs
            
        #print(deletion_locs)
        
        
    def deletions_and_insertions_for_intervals(
        self, samfile, bamfilename, intervals, outdir='/Users/dp/pma/rbfox2_eCLIP_comparison/cims/'):

        print(f'# reads to obtain: {len(intervals)}')
        results = []
        n_reads_done = 0
        no_del_or_ins = 0
    
        for iv in intervals:

            if len(iv) == 4:
                iv = iv[:-1]

            for n, r in enumerate(samfile.fetch(*iv)):
                deletions, deleted_bases, insertions, inserted_bases = \
                    self.deletions_and_insertions_in_read(r)

                if deletions == [] and insertions == []:
                    no_del_or_ins += 1
                else:
                    results.append([len(deletions), deleted_bases, len(insertions), inserted_bases])
                n_reads_done += 1

                not (n_reads_done % 100000) and print(f'{n_reads_done}/{len(intervals)}', end=', ')
      
                
        os.makedirs(f'{outdir}/deletions/', exist_ok=True)
        with open(f"{outdir}/deletions/{os.path.basename(bamfilename)}", 'w') as f:
            f.write('Insertions\tInserted_bases\tDeletions\Deleted_bases\tcount\n')
            f.write(f'0\t[]\t0\t[]\t{no_del_or_ins}\n')
            for row in results:
                f.write('\t'.join([str(x) for x in row]) + '\t1\n')
        
    def deletions_and_insertions_for_random_intervals(self, samfile, bamfilename, N_reads=1000, outdir='.'):
        #ivs = get_random_genomic_intervals(bam_fname=bamfilename, interval_len=1E4, N_intervals=1000)
        
        print(f'# reads to obtain: {N_reads}')
        results = []
        n_reads_done = 0
        no_del_or_ins = 0
        while n_reads_done < N_reads:
            
            ivs = get_random_genomic_intervals(bam_fname=bamfilename, interval_len=1E5, N_intervals=1000)

            for iv in ivs:
                #print(iv)

                if len(iv) == 4:
                    iv = iv[:-1]

                for n, r in enumerate(samfile.fetch(*iv)):
                    #print(n, '\n', r)
                    deletions, deleted_bases, insertions, inserted_bases = \
                        self.deletions_and_insertions_in_read(r)
                    
                    if deletions == [] and insertions == []:
                        no_del_or_ins += 1
                    else:
                        results.append([len(deletions), deleted_bases, len(insertions), inserted_bases])
                    n_reads_done += 1
                    
                    not (n_reads_done % 100000) and print(f'{n_reads_done}/{N_reads}', end=', ')
                
                    if n_reads_done >= N_reads:
                        break
                        
                if n_reads_done >= N_reads:
                    print('Finished.')
                    break
                
        os.makedirs(f'{outdir}/deletions/', exist_ok=True)
        with open(f"{outdir}/deletions/{os.path.basename(bamfilename)}", 'w') as f:
            f.write('Insertions\tInserted_bases\tDeletions\Deleted_bases\tcount\n')
            f.write(f'0\t[]\t0\t[]\t{no_del_or_ins}\n')
            for row in results:
                f.write('\t'.join([str(x) for x in row]) + '\t1\n')
                
    def deletions_and_insertions_in_interval(self, interval, samfile):
        
        if len(interval) == 3:
            iv = interval
        elif len(interval) == 4:
            iv = interval[:-1]
        else:
            print(f"Invalid interval len: {interval}")
            
        #print(f"{iv[0]}:{iv[1]}-{iv[2]}")

        for r in samfile.fetch(*iv):
            return self.deletions_and_insertions_in_read(r)
        
    def deletions_and_insertions_in_read(self, r):
        """Return: 
            List[genomic coordinates of deletions], List[bases deleted],
            List[genomic coordinates of insertions], List[bases inserted]
        """
        
        pos_in_query, pos_in_ref = 0, 0
        deletions, insertions = [], []  # Left positions of D/I interval.
        deleted_bases, inserted_bases =  [], []
        for (flag, n_bases) in r.cigar:
            if flag == 2:  # Deletion.
                deletions.append(pos_in_ref + r.pos)  # Genomic coord.
                if self.genomic_fasta is not None:
                    fastaN = pos_in_ref+r.pos
                    # Against the laws of God and man, bam 0-indexes the chrom names so 'chr1' -> int(0).
                    #print(f"chr{r.rname+1}:{fastaN}-{fastaN+n_bases}")
                    if f'chr{r.rname+1}' in self.genomic_fasta:
                        deleted_bases.append(self.genomic_fasta[f'chr{r.rname+1}'][fastaN:fastaN+n_bases])
                    else:
                        deleted_bases.append('X')
                else:
                    deleted_bases.append('X' * n_bases)
                        
            if flag == 1:  # Insertion.
                insertions.append(pos_in_query + r.pos)  # Genomic cood.
                # SEQ in our bams is the sequenced read, from genomic L>R, even if -.
                inserted_bases.append(r.seq[pos_in_query:pos_in_query+n_bases])  
                
            if flag not in [2,3,5,6]:
                pos_in_query += n_bases
            
            if flag not in [1,4,5,6]:
                pos_in_ref += n_bases
        #print(deletions, deleted_bases, insertions, inserted_bases)
        #print([r.seq[x-5:x+2] for x in deletions])
        return deletions, deleted_bases, insertions, inserted_bases
                                   
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Edit read names in a bam file to remove fixed sequences.')
    parser.add_argument('-i', help='Folder of bam files.')
    parser.add_argument('-f', help='Genomic fasta file.', default='None')
    parser.add_argument('-o', help='Folder to write output bedgraphs to.')
    parser.add_argument('--iv', help='Genomic interval in format 3:300-400/+.')
    
    args = parser.parse_args()
    
    if args.f == 'None':
        args.f = None
        
    chrom = args.iv.split(':')[0]
    strand = args.iv.split('/')[-1]
    start, end = args.iv.split(':')[1].split('/')[0].split('-')
    start, end = int(start.replace(',','')), int(end.replace(',',''))
    args.iv = [chrom, start, end, strand]
    
    dc = deletionCounter(genomic_fasta=args.f)
    dc.write_deletion_positions_to_bedgraph(
        list(glob.glob(f"{args.i}/*bam")), interval=[chrom, start, end, strand], outdir=args.o)