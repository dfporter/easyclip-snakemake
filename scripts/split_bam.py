"""
Heavily modified from https://bioinformatics.stackexchange.com/questions/12868/is-there-a-command-line-tool-to-split-a-sam-bam-file-by-cb-cell-barcode-tag

MIT License
Copyright (c) 2020 Warren W. Kretzschmar
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import pysam, os, sys, subprocess
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import scripts
import scripts.scheme
import random

class BamWriter:
    def __init__(self, alignment, barcodes, barcode_to_fname, prefix):
        self.alignment = alignment
        self.prefix = prefix
        self.barcodes = set(barcodes)
        self.barcode_to_fname = barcode_to_fname
        self.outfiles = {}
        files_opened = {}  # Avoid errors when multiple raw input fastqs have the same output.
        
        for barcode in barcodes:
            out_bam_name = f"{self.prefix}{barcode_to_fname[barcode]}.bam"
            if barcode not in self.outfiles:
                if out_bam_name not in files_opened:
                    self.outfiles[barcode] = self.open_bam(out_bam_name)
                    files_opened[out_bam_name] = barcode
                else:
                    self.outfiles[barcode] = self.outfiles[files_opened[out_bam_name]]
        
        #print(f"bamWriter.__init__(): outfiles = {self.outfiles}")
        
        self.outfile_names = list(set([
            f"{self.prefix}{barcode_to_fname[barcode]}.bam" for barcode in barcodes]))

        # Make copies of the barcode with all possible 1 nt SNP that all
        # output to the same file. If two barcodes differ by only 1 nt then this fails,
        # but such a design should never be used anyway.
        bases = ['A', 'T', 'C', 'G', 'N']

        for barcode in list(self.barcodes):
            l5, l3 = barcode.split('|')[0].split('__')
            pcr_index = barcode.split('|')[-1]
            
            bc_list = list(l5)

            for pos in range(len(l5)):
                for base in bases:
                    new_l5_bc = bc_list[:]
                    new_l5_bc[pos] = base
                    new_l5_bc = ''.join(new_l5_bc)
                    
                    l3_bc_list = list(l3)
                    for pos in range(len(l3)):
                        for l3_base in bases:
                            new_l3_bc = l3_bc_list[:]
                            new_l3_bc[pos] = l3_base
                            new_l3_bc = ''.join(new_l3_bc)
                            new_bc = f"{new_l5_bc}__{new_l3_bc}|{pcr_index}"
                            if new_bc not in self.outfiles:
                                self.outfiles[new_bc] = self.outfiles[barcode]

    def write_record_to_barcode(self, rec, barcode):
        if barcode not in self.outfiles.keys():
            barcode = 'unrecognized'

        if barcode not in self.outfiles:
            self.outfiles[barcode] = pysam.AlignmentFile(
                f"{self.prefix}{barcode}.bam", "wb", template=self.alignment)
            
        self.outfiles[barcode].write(rec)

    def open_bam(self, name):
        return pysam.AlignmentFile(
            name, "wb", template=self.alignment)
        
def split_bam(
    input_bam: str, barcodes: list, barcode_to_fname: dict, output_prefix: str,
    l5_inline_pattern="BBBBBBNNNNNNN", l3_inline_pattern="BBBBBNNNNN"):
    """Split a bam into separate bams based on barcodes in the read name.
    input_bam: bam filename.
    barcodes: list of strings, [L5BC_L3BC, L5BC_L3BC, ect.], found in the samples file.
    barcode_to_fname: assigns each string of L5/L3 barcode pairs to a filename.
    output_prefix: prefix for output filenames: {self.prefix}{barcode_to_fname[barcode]}.bam
    
    l5_inline_pattern/l3_inline_pattern: The L5_inline entry represents the 
    barcode (B) and UMI (N) in the order it would be sequenced in read 1 
    (that is, 6 barcode bases and then 7 UMI bases). The L3_inline entry is the same,
    but in the order it would be sequenced from read 2. These values are from the
    configfile.
    """
    
    if '/' in output_prefix:
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    alignment = pysam.AlignmentFile(input_bam)
    
    print(f"barcodes: {barcodes}")
    writer = BamWriter(
        alignment=alignment, barcodes=barcodes,
        barcode_to_fname=barcode_to_fname, prefix=output_prefix)
    
    recs = alignment.fetch()
    # Read name suffix format: 3773_AGCTAGAAAATCG_AGT
    # AGCTAGAAAATCG_AGT -> AGCTAG=L5 BC. AGT=L3 BC. AAAATCG=UMI.
    
    l5_bc_len = len(barcodes[0].split('|')[0].split('__')[0])
    l3_bc_len = len(barcodes[0].split('|')[0].split('__')[1])
    
    l5_bc_pos = [pos for pos, base in enumerate(l5_inline_pattern) if base == 'B']
    l3_bc_pos = [pos for pos, base in enumerate(l3_inline_pattern) if base == 'B']
    
    for rec in recs:
        try:
            try:
                l5, l3 = rec.query_name.split('__')[-1].split('|')[0].split('-')
            except:
                print(rec.query_name)
                sys.exit()
            pcr_index = rec.query_name.split('__')[-1].split('|')[-1]
            
            l5 = ''.join([l5[i] for i in l5_bc_pos])
            l3 = ''.join([l3[i] for i in l3_bc_pos])
            barcode = f"{l5}__{l3}|{pcr_index}"

            writer.write_record_to_barcode(rec=rec, barcode=barcode)
        except KeyError:
            pass

    for fh in writer.outfiles.values():
        fh.close()
    
    for bamfilename in writer.outfile_names:
        
        if not os.path.exists(bamfilename):
            continue
           
        base = os.path.splitext(bamfilename)[0]
        try:
            
            cmd = f"samtools sort -o {base}.sorted.bam {bamfilename}"
            print(cmd)
            res = subprocess.check_output(cmd.split(' '))

            cmd = f"mv {base}.sorted.bam {bamfilename}"
            print(cmd)
            res = subprocess.check_output(cmd.split(' '))

            cmd = f"samtools index {bamfilename}"
            print(cmd)
            res = subprocess.check_output(cmd.split(' '))
        except:
            print(f"Error processing {bamfilename}. Tried commands: ")
            print(f"samtools sort -o {base}.sorted.bam {bamfilename}")
            print(f"mv {base}.sorted.bam {bamfilename}")
            print(f"samtools index {bamfilename}")
        
if __name__ == '__main__':
    import pandas, re
    df = pandas.read_csv('testdata/samples.txt', sep='\t')
    barcodes = ['__'.join([x,y]) for x,y in zip(df['L5_BC'], df['L3_BC'])]

    scheme = scripts.scheme.scheme('testdata/samples.txt')
    barcode_to_fname = {
        '_'.join(l5_l3): os.path.splitext(scheme.p6p3_to_long_filename_r1[l5_l3])[0] \
        for l5_l3 in zip(scheme.scheme_df['L5_BC'], scheme.scheme_df['L3_BC'])}

    df = pandas.read_csv('testdata/samples.txt', sep='\t')
    df['Gene'] = [re.sub(' ', '-', x) for x in df['Gene']]  # Get rid of spaces.
    df = df.loc[[type(x)==type('') for x in df['Gene']], :]

    # Define the samples list used throughout the workflow.
    samples = [f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}" for exp,protein,l5_bc,l3_bc,rep in zip(
        df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'])]

    # Read name suffix format example: 3773_AGCTAGAAAATCG_AGT
    # AGCTAGAAAATCG_AGT -> AGCTAG=L5 BC. AGT=L3 BC. AAAATCG=UMI.

    print('=' * 140)
    print(df.head())

    barcodes = []
    barcode_to_fname = {}
    for exp,protein,l5_bc,l3_bc,rep,r1_fastq in zip(
        df['Experiment'], df['Gene'], df['L5_BC'], df['L3_BC'], df['Replicate'], df['R1_fastq']):

        pcr_prefix = r1_fastq.split('R1.fastq')[0].split('R1.fq')[0].split('1.fq')[0].split('1.fastq')[0]
        barcode = f"{l5_bc}__{l3_bc}|{pcr_prefix}"
        fname = f"{exp}_{protein}_{rep}_{l5_bc}_{l3_bc}"
        barcode_to_fname[barcode] = fname
        barcodes.append(barcode)

    print(barcodes)
    print(barcode_to_fname)
    split_bam('testdata/sams/all_reads.bam', barcodes, barcode_to_fname, 'testdata/sams/split/',
             l5_inline_pattern="BBBBBBNNNNNNN", l3_inline_pattern="BBB")