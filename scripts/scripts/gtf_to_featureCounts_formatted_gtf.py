import re, sys
"""
Gencode format (e.g., gencode.v37.annotation.gtf):
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";

featureCounts format:
GeneID	Chr	Start	End	Strand
497097	chr1	3204563	3207049	-
"""


def reformat_gtf(input_file, output_file):
    outf = open(output_file, 'w')

    outf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    
    with open(input_file) as f:
        for li in f:
            
            if li[0] == '#':
                continue
            
            if (('transcript_support_level "1"' in li) or ('transcript_support_level "NA"' in li)):
                
                s = li.split("\t")
                
                if s[2] == 'transcript':# or s[2] == 'intron'

                    gene_name = re.search('gene_name "([^"]+)"', s[-1])

                    if gene_name is not None:
                        outf.write(f"{gene_name.group(1)}::{s[2]}\t{s[0]}\t{s[3]}\t{s[4]}\t{s[6]}\n")
                    
    outf.close()

    
if __name__ == '__main__':
    input_file, output_file = sys.argv[1], sys.argv[2]
    reformat_gtf(input_file, output_file)