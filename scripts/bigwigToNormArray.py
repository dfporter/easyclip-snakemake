import pyBigWig, re, gffutils, os, sys, collections, pandas, glob, random
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats
"""

The scripts exonTools, regionDef, and signalMatrix were created specifically to
test signals around differential exons. These scripts regionsFromDB and bigwigToNormArray
reflect the more general case of normalizing coverage across genomic elements.

"""
class signalToNormArray:
    
    def __init__(self, bw_fnames, bw_uses_chr_in_chrom_names=True):
        self.bw_fnames = bw_fnames
        self.bw_uses_chr_in_chrom_names = bw_uses_chr_in_chrom_names
        
    @staticmethod
    def basename(bw_fname):
        return os.path.splitext(os.path.basename(bw_fname.replace('.-.', '.').replace('.+.', '.')))[0]
    
    def read_coverages_in_regions(self, regions, n_features: int, regions_are_split=False, target_size=50):
        """regions: iterable of [tuple(chr, start, end, strand), ...]
        target_size: what length to normalize each region to.
        returns:
            coverages: dict of {bw_name: 2D numpy array}
            
        If regions_are_split:
            regions: [ [tuple(chr, start, end, strand), tuple(chr, start, end, strand), ...], ...]
        """
        
        coverages = {}
        print(f"Reading coverages in {len(regions)} regions. Regions are split? {regions_are_split}...")
        
        # Initialize.
        for bw in self.bw_fnames:
            coverages[self.basename(bw)] = np.zeros((n_features, target_size))
            
        read_bw_fnames = set()
        # For each bigwig file, read normalized signal.
        for bw_fname in self.bw_fnames:
            
            if os.path.exists(bw_fname):
                print(f"Reading {bw_fname}...")
            else:
                print(f"{bw_fname} does not exist!")
                
            if bw_fname in read_bw_fnames:
                continue  # Already loaded (with -/+ strand file, probably).
            
            base = self.basename(bw_fname)
            
            if '.+.' not in bw_fname and '.-.' not in bw_fname:
                #stranded_bw = False
                
                bw = pyBigWig.open(bw_fname)
                for row_n, iv_tuple in enumerate(regions):
                    coverages[base][row_n] = self.len_norm_signal(
                        bw, iv_tuple, target_size=target_size, regions_are_split=regions_are_split)
                if np.nansum(coverages[base]) > 0:
                    print('->', np.nansum(coverages[base]))
                bw.close()
                
                read_bw_fnames.add(bw_fname)
                
            else:  # Stranded bigwigs.
                bw_pos_fname = re.sub('\.-\.', '.+.', bw_fname)
                bw_neg_fname = re.sub('\.\+\.', '.-.', bw_fname)
                #print(f"For {bw_fname}, loading {bw_pos_fname} and {bw_neg_fname}")
            
                bw_pos, bw_neg = pyBigWig.open(bw_pos_fname), pyBigWig.open(bw_neg_fname)
                
                print(f"Iterating through regions...")
                for row_n, iv_tuple in enumerate(regions):
                    
                    #if not row_n % 10:
                    #print(row_n, iv_tuple, end=',')
                    if regions_are_split:
                        strand = iv_tuple[0][-1]
                    else:
                        strand = iv_tuple[-1]
                    
                    #print(f'strand={strand}')
                    if strand == '+':
                        coverages[base][row_n] = self.len_norm_signal(
                            bw_pos, iv_tuple, target_size=target_size, regions_are_split=regions_are_split)
                    elif strand == '-':
                        coverages[base][row_n] = self.len_norm_signal(
                            bw_neg, iv_tuple, target_size=target_size, regions_are_split=regions_are_split)
                    else:
                        print("Error: had apparently stranded bigwig files, (names had .+. or .-.), but regions" + \
                             f" not stranded. Interval={iv_tuple}. Expected last element + or -. fname: {bw_fname}")
                        
                if np.nansum(coverages[base]) > 0:
                    print('->', np.nansum(coverages[base]))
                bw_pos.close(); bw_neg.close()
                
                read_bw_fnames |= set([bw_pos_fname, bw_neg_fname])
                
        print(f"Finished reading coverages.")
        return coverages

    def norm_len(self, a, target_size=50) -> list:
        
        if len(a) == 0:
            return [np.nan for x in target_size]
        
        if len(a) >= target_size:
            y = np.array_split(a, target_size)

        else:
            f = int(np.ceil(target_size/len(a)))  # Round up.
            new = np.zeros(f * len(a))
            for n, v in enumerate(a):
                new[f*n:(f*n)+f] = v
            y = np.array_split(new, target_size)     
            
        b = [np.mean(x) for x in y]
        #print('b=', b)
#        b = [np.mean(a[int(k*binsize):int(k*binsize)+int(binsize)]) for k in range(0,int(len(a)/binsize))]
        return b

    def is_iv_ok(self, interval, bw):
        #print(f'is interval {interval} ok?')
        if interval[2] == interval[1]:
            print('0 len', interval)
            return False
        if interval[1] < 0:
            print('<0 start', interval)
            return False
        try:
        #print('sig_arr -->')
            sig_arr = np.nan_to_num(bw.values(*interval[:-1]))
        except:
            # Invalid interval bounds.
            print("Invalid interval:", interval)
            return False
        #print("got sig array")
        return True
        
    def len_norm_signal_one_interval(self, bw, interval, target_size=20) -> list:
        #if type(interval) != type(tuple()):
        #    return [np.nan] * target_size
        
        if self.bw_uses_chr_in_chrom_names:
            interval = ('chr' + interval[0].split('chr')[-1], interval[1], interval[2], interval[3])
            
        #print(f"len_norm_signal_one_interval: interval={interval}")
        if not is_iv_ok(interval, bw):
            #print("interval bad")
            return [np.nan] * target_size
        
        if not(np.any(np.isnan(sig_arr))):
            
            sig_arr = self.norm_len(sig_arr, target_size=target_size)
            
            #print("Found signal")
            # Flip if negative strand.
            if interval[-1] == '-':
                return np.flip(sig_arr)
            elif interval[-1] == '+':
                return sig_arr
            else:
                print(f"Error, expected strand info in interval {interval}")
                
        return [np.nan] * target_size
    
    def signal_one_interval(self, bw, interval, target_size=20) -> list:
        if self.bw_uses_chr_in_chrom_names:
            interval = ('chr' + interval[0].split('chr')[-1], interval[1], interval[2], interval[3])
        
        #print(f"signal_one_interval {interval}")
        if not self.is_iv_ok(interval, bw):
            print("bad interval")
            return [np.nan] * target_size
        
        sig_arr = np.nan_to_num(bw.values(*interval[:-1]))
        
        # Flip if negative strand.
        if interval[-1] == '-':
            return np.flip(sig_arr)
        elif interval[-1] == '+':
            return sig_arr
        else:
            print(f"Error, expected strand info in interval {interval}")
        #print("returning nan from signal_one_interval")
        return [np.nan] * target_size
    
    def len_norm_signal(self, bw, interval, target_size=20, regions_are_split=False) -> list:
        #print(f"len_norm_signal inerval={interval}")
        if not regions_are_split:
            return self.len_norm_signal_one_interval(bw, interval, target_size=target_size)
        else:
            arrs = []
            if interval[0][-1] == '-':
                for iv in interval[::-1]:
                    arrs.append(self.signal_one_interval(bw, iv, target_size=target_size))
            else:
                for iv in interval:
                    arrs.append(self.signal_one_interval(bw, iv, target_size=target_size))
            sig_arr = np.concatenate(arrs)
            
            if not(np.any(np.isnan(sig_arr))):
                return self.norm_len(sig_arr, target_size=target_size)
            
            return [np.nan] * target_size
        
class signalToUnnormArray:
    
    def __init__(self, bw_fnames, bw_uses_chr_in_chrom_names=True):
        self.bw_fnames = bw_fnames
        self.bw_uses_chr_in_chrom_names = bw_uses_chr_in_chrom_names
        
    @staticmethod
    def basename(bw_fname):
        return os.path.splitext(os.path.basename(bw_fname.replace('.-.', '.').replace('.+.', '.')))[0]
    
    def read_coverages_in_regions(self, regions, n_features: int, regions_are_split=False, target_size=50):
        """regions: iterable of [tuple(chr, start, end, strand), ...]
        target_size: what length to normalize each region to.
        returns:
            coverages: dict of {bw_name: 2D numpy array}
            
        If regions_are_split:
            regions: [ [tuple(chr, start, end, strand), tuple(chr, start, end, strand), ...], ...]
        """
        
        coverages = {}
        print(f"Reading coverages in {len(regions)} regions. Regions are split? {regions_are_split}...")
        
        # Initialize.
        for bw in self.bw_fnames:
            coverages[self.basename(bw)] = []#np.zeros((n_features, target_size))
            
        read_bw_fnames = set()
        # For each bigwig file, read normalized signal.
        for bw_fname in self.bw_fnames:
            
            if os.path.exists(bw_fname):
                print(f"Reading {bw_fname}...")
            else:
                print(f"{bw_fname} does not exist!")
                
            if bw_fname in read_bw_fnames:
                continue  # Already loaded (with -/+ strand file, probably).
            
            base = self.basename(bw)
            
            if '.+.' not in bw_fname and '.-.' not in bw_fname:
                #stranded_bw = False
                
                bw = pyBigWig.open(bw_fname)
                for row_n, iv_tuple in enumerate(regions):
                    coverages[base].append(
                        self.signal(
                            bw, iv_tuple, target_size=target_size, regions_are_split=regions_are_split))
                bw.close()
                
                read_bw_fnames.add(bw_fname)
                
            else:  # Stranded bigwigs.
                bw_pos_fname = re.sub('\.-\.', '.+.', bw_fname)
                bw_neg_fname = re.sub('\.\+\.', '.-.', bw_fname)
                #print(f"For {bw_fname}, loading {bw_pos_fname} and {bw_neg_fname}")
            
                bw_pos, bw_neg = pyBigWig.open(bw_pos_fname), pyBigWig.open(bw_neg_fname)
                
                for row_n, iv_tuple in enumerate(regions):
                    if regions_are_split:
                        strand = iv_tuple[0][-1]
                    else:
                        strand = iv_tuple[-1]                        

                    if strand == '+':
                        coverages[base].append(
                            self.signal(
                                bw_pos, iv_tuple, target_size=target_size, regions_are_split=regions_are_split))
                    elif strand == '-':
                        coverages[base].append(
                            self.signal(
                                bw_neg, iv_tuple, target_size=target_size, regions_are_split=regions_are_split))
                    else:
                        print("Error: had apparently stranded bigwig files, (names had .+. or .-.), but regions" + \
                             f" not stranded. Interval={iv_tuple}. Expected last element + or -. fname: {bw_fname}")
                        raise Exception
                bw_pos.close(); bw_neg.close()
                
                read_bw_fnames |= set([bw_pos_fname, bw_neg_fname])
                
        print(f"Finished reading coverages.")
        return coverages

    def norm_len(self, a, target_size=50) -> list:
        
        if len(a) == 0:
            return [np.nan for x in target_size]
        
        if len(a) >= target_size:
            y = np.array_split(a, target_size)

        else:
            f = int(np.ceil(target_size/len(a)))  # Round up.
            new = np.zeros(f * len(a))
            for n, v in enumerate(a):
                new[f*n:(f*n)+f] = v
            y = np.array_split(new, target_size)     
            
        b = [np.mean(x) for x in y]
        #print('b=', b)
#        b = [np.mean(a[int(k*binsize):int(k*binsize)+int(binsize)]) for k in range(0,int(len(a)/binsize))]
        return b

    def is_iv_ok(interval, bw):
        if interval[2] == interval[1]:
            return False
        try:
            sig_arr = np.nan_to_num(bw.values(*interval[:-1]))
        except:
            # Invalid interval bounds.
            print("Invalid interval:", interval[:-1])
            return False
        
        return True
    
    def signal_one_interval(self, bw, interval, target_size=20) -> list:
        if self.bw_uses_chr_in_chrom_names:
            interval = ('chr' + interval[0].split('chr')[-1], interval[1], interval[2], interval[3])
            
        if not self.is_iv_ok(interval, bw):
            return [np.nan] * target_size
        
        if not(np.any(np.isnan(sig_arr))):
            
            # Flip if negative strand.
            if interval[-1] == '-':
                return np.flip(sig_arr)
            elif interval[-1] == '+':
                return sig_arr
            else:
                print(f"Error, expected strand info in interval {interval}")
                
        return [np.nan] * target_size
    
    def signal(self, bw, interval, target_size=20, regions_are_split=False) -> list:
        if not regions_are_split:
            return self.signal_one_interval(bw, interval, target_size=target_size)
        else:
            arrs = []
            if interval[0][-1] == '-':
                for iv in interval[::-1]:
                    arrs.append(self.signal_one_interval(bw, iv, target_size=target_size))
            else:
                for iv in interval:
                    arrs.append(self.signal_one_interval(bw, iv, target_size=target_size))
            return np.concatenate(arrs)
