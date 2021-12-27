from typing import Mapping, List
import os, glob


def for_split_bedgraph_lines(fname: str) -> List:
    with open(fname) as f:
        try:
            next(f)  # Skip header.
        except:
            yield ('', 0, 1, 0)
        for li in f:
            sp = li.rstrip('\n').split('\t')
            # (_, start, end, value)
            yield (sp[0], int(sp[1]) , int(sp[2]), sp[3])


def get_auc(fname: str) -> float:
    auc = 0
    for s in for_split_bedgraph_lines(fname):
        # (end - start) * value
        try:
            auc += (s[2] - s[1]) * float(s[3])
        except:
            auc += 0
            print(f"scripts.total_read_numbers.get_auc({fname}):\nError processing values {s}.")
            print(f"Entered zero signal for the line. Current total auc for this file is {auc}.")
    return auc


def total_read_numbers(folder: str, outfile='./data/total_read_numbers.txt') -> Mapping[str, float]:
    """Total area under the curve for bedgraphs in a folder.
    Assuming each read is a single point with value 1, this is the total read number.
    Write the results to a file.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    
    aucs = {}
    for fname in glob.glob(f"{folder}/*.+.wig"):
        print(f"{fname}...", end=" ")
        aucs[fname] = get_auc(fname)
        minus_fname = fname.split('.+.wig')[0] + '.-.wig'
        aucs[fname] += get_auc(minus_fname)
        print(f"AUC={aucs[fname]:,}")
    
    original_aucs = {k:v for k,v in aucs.items()}
    
    # For good measure, write the same values with some different versions of the filename.
    aucs.update({os.path.basename(k).rstrip('.+.wig'): v for k,v in original_aucs.items()})
    aucs.update({os.path.basename(k): v for k,v in original_aucs.items()})
    
    with open(outfile, 'w') as f:
        f.write("Dataset\tTotal read number\n")
        [f.write(f"{name}\t{val}\n") for name, val in aucs.items()]
    
    print(f"Wrote total read numbers to {outfile}.")
    return aucs