import pysam, os, subprocess, argparse


def edit_read_names_in_bam(
    input_bam: str, output_bam: str,
    l5_inline_pattern="BBBBBBNNNNNNN", l3_inline_pattern="BBBBBFFFFF"):
    """
    Remove fixed sequences that had been moved to the read name as "UMIs" but
    which are fixed. This is to be done before removing duplicates.
    
    input_bam: bam filename.
    
    The original l5_inline_pattern/l3_inline_pattern: The config lists L5 strings of 
    barcode (B) and UMI (N) in the order it would be sequenced in read 1 
    (that is, 6 barcode bases and then 7 UMI bases). The L3_inline entry is the same,
    but in the order it would be sequenced from read 2. These values are from the
    configfile.
    
    Which bases are to be removed are given in string format. For an input L3 inline
    sequence of B=barcode and N="UMI":
    BBBBBNNNNN
    
    If the N's were actually fixed, we pass the following string to have all fixed 
    bases removed:
    FFFFFFFFFF
    
    That is, an "F" denotes a fixed base, which is cut, and all other letters are kept.
    """
    
    if '/' in output_bam:
        os.makedirs(os.path.dirname(output_bam), exist_ok=True)

    alignment = pysam.AlignmentFile(input_bam)
    
    l5_bases_to_keep = [x for x in range(len(l5_inline_pattern)) if l5_inline_pattern[x]!='F']
    l3_bases_to_keep = [x for x in range(len(l3_inline_pattern)) if l3_inline_pattern[x]!='F']
    
    writer = pysam.AlignmentFile(output_bam, "wb", template=alignment)
    
    recs = alignment.fetch()    
        
    for rec in recs:
        try:
            try:
                l5, l3 = rec.query_name.split('__')[-1].split('|')[0].split('-')
            except:
                print(rec.query_name)
                sys.exit()
            pcr_index = rec.query_name.split('__')[-1].split('|')[-1]
            
            if len(l5_bases_to_keep):
                l5 = ''.join([l5[x] for x in l5_bases_to_keep])
            else:
                l5 = 'N'
            if len(l3_bases_to_keep):
                l3 = ''.join([l3[x] for x in l3_bases_to_keep])
            else:
                l3 = 'N'
            
            _name = rec.query_name.split('__')[0]
            rec.query_name = f"{_name}__{l5}-{l3}"
            
            writer.write(rec)
            
        except KeyError:
            pass

    writer.close()

    cmd = f"samtools index {output_bam}"
    print(cmd)
    res = subprocess.check_output(cmd.split(' '))
        
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Edit read names in a bam file to remove fixed sequences.')
    parser.add_argument('-i', help='Input bam file.')
    parser.add_argument('-o', help='Output bam file.')
    parser.add_argument('--l5', help='String denoting what bases to remove from the L5 inline barcode.' +\
                       " 'F' denotes a fixed base to remove. Example: FFFFFFFNNN keeps only the final 3 bases.")
    parser.add_argument('--l3', help='String denoting what bases to remove from the L3 inline barcode.' +\
                       " 'F' denotes a fixed base to remove. Example: FFFFFFFNNN keeps only the final 3 bases.")
    args = parser.parse_args()
    print(args.i, args.o, args.l5, args.l3)
    
    edit_read_names_in_bam(
        input_bam=args.i, output_bam=args.o,
        l5_inline_pattern=args.l5, l3_inline_pattern=args.l3)
