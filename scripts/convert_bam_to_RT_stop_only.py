import pysam, sys, os
from itertools import groupby


bam = sys.argv[1]
out = sys.argv[2]

bam = pysam.AlignmentFile(bam, 'rb')

sam = open(out, 'w')
sam.write(str(bam.header))

for r in bam.fetch():
    
    # BAM: 0-based, SAM: 1-based.
    if r.is_reverse: 
        query_sequence = r.query_sequence[-1] 
        query_length = 1
        ref_start = r.reference_start + 1 
    else: 
        query_sequence = r.query_sequence[0] 
        query_length = 1
        ref_start = r.reference_end 
        
    read_str = f"{r.query_name}\t{r.flag}\t{r.reference_name}\t{ref_start}\t{r.mapping_quality}\t1M\t*\t0\t1\t{query_sequence}\tF\n"
    sam.write(read_str)


sam.close()
