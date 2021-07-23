import scipy, os, pysam, importlib
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scripts.deletionCounter
importlib.reload(scripts.deletionCounter)

def get_coverage_in_an_interval_for_bam_file(bam_fname, interval):
    #count_coverage(self, contig, start=None, stop=None, region=None, quality_threshold=15, read_callback='all', reference=None, end=None)
    samfile = pysam.AlignmentFile(bam_fname, "rb" )
    # Get np.array of coverage depth across interval.

    if type(interval) == type(''):
        b = interval.strip("()[]").split(', ')
        interval = [b[0].strip("''"), int(b[1]), int(b[2])]

    iv = interval[:3]

    arr = np.sum(samfile.count_coverage(*iv), axis=0)
    samfile.close()
    return arr

# Get the average correlations.
def get_corr_matrix_for_iv(
    data1, data2, replicates=['easy1', 'easy2', 'easy3', 'eclip1', 'eclip2'],
    plot=False, verbose=False, method='spearman'):
    # Create a cross-correlation matrix:
    # Rows: easyCLIP replicates then eCLIP replicates.
    # Columns: easyCLIP replicates then eCLIP replicates (same order).
    
    corr_matrix = np.zeros(shape=(len(replicates), len(replicates)))

    done_a = set()
    for a, name_a in zip(data1 + data2, replicates):

        for b, name_b in zip(data1 + data2, replicates):
            if name_b in done_a:
                continue

            #print(f"Correlating len {len(a)} with {len(b)}.")
            if method == 'spearman':
                r = scipy.stats.spearmanr(a, b).correlation
            elif method == 'pearson':
                r = scipy.stats.pearsonr(a, b)[0]
            
            verbose and print(f"{name_a} v {name_b}", r)
            #print("..-> r =", r)
            corr_matrix[replicates.index(name_a)][replicates.index(name_b)] = r
        done_a.add(name_a)

    verbose and print(corr_matrix)

    if plot:
        sns.heatmap(corr_matrix, cmap='binary')
        plt.show();plt.clf();plt.close()
    
    return corr_matrix


class correlations():
    
    def __init__(self):
        self.matrices = {}
        
    def add(self, matrix, name):
        if np.any([np.isnan(x) for x in matrix]):
            print("Skipping due to NaN values.")
        else:
            self.matrices[name] = matrix
    
    def average(self):
        mm = [v for v in self.matrices.values()]
        self.mean = np.average(mm, axis=0)
        return np.average(mm, axis=0)

    
def smooth(y, box_pts):
    box_pts = int(box_pts)
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def get_corr_matrix_of_deletions_vs_rt_stops_for_list_of_ivs(
    fnames1: list, fnames2: list, reps1: list, reps2: list,
    intervals: list, apply_a_cutoff=False):


    dc = scripts.deletionCounter.deletionCounter()

    corrs = correlations()
    replicates = reps1 + reps2
    for n, iv in enumerate(intervals):
        if n % 20 == 0:
            print(f"{n}: {len(corrs.matrices)}", end=', ')
            
        data1 = [list(dc._deletion_locations_for_intervals(f, [iv]).values())[0] for f in fnames1]
        data2 = [np.abs(get_coverage_in_an_interval_for_bam_file(f, iv)) for f in fnames2]
        
        if apply_a_cutoff:
            if np.max(data1) < apply_a_cutoff:
                continue
        
        window = 200
        #data1 = [smooth(_, window) for _ in data1]
        data2 = [smooth(_, window) for _ in data2]
        
        corrs.add(get_corr_matrix_for_iv(
            data1, data2, replicates=replicates, verbose=False), str(iv))        


    return corrs

def bin_array(a, bin_width=50):
    if type(a) != type(np.array([1])):
        a = np.array(a)
    bin_width = 50
    out_arr = []
    #print(a)
    for b in range(0, len(a), bin_width):
        #print(a[b:b+50].sum())
        out_arr.append(a[b:b+bin_width].sum())
    return np.array(out_arr)


def get_corr_matrix_for_list_of_ivs(
    fnames1: list, fnames2: list, reps1: list, reps2: list,
    intervals: list, apply_a_cutoff=False, smoothing=False, method='spearman',
    bin_window=50):
    
    corrs = correlations()
    replicates = reps1 + reps2
    for n, iv in enumerate(intervals):
        if n % 20 == 0:
            print(f"{n}: {len(corrs.matrices)}", end=', ')

        data1 = [np.abs(get_coverage_in_an_interval_for_bam_file(f, iv)) for f in fnames1]
        data2 = [np.abs(get_coverage_in_an_interval_for_bam_file(f, iv)) for f in fnames2]
        
        if apply_a_cutoff:
            if np.max(data1) < apply_a_cutoff:
                continue
        
        if bin_window:
            data1 = [bin_array(arr, bin_width=bin_window) for arr in data1]
            data2 = [bin_array(arr, bin_width=bin_window) for arr in data2]
            #print(data1, data2)
        elif smoothing:
            window = int(smoothing)
            data1 = [smooth(_, window) for _ in data1]
            data2 = [smooth(_, window) for _ in data2]
        
        corrs.add(get_corr_matrix_for_iv(
            data1, data2, replicates=replicates, method=method, verbose=False), str(iv))        


    return corrs
