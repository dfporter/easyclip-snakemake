import pandas, os, re

"""
if not os.path.exists("assets/reference/Dfam_curatedonly.embl"):
    os.makedirs("assets/reference/", exist_ok=True)
    os.chdir("assets/reference/")
    os.system("wget https://www.dfam.org/releases/Dfam_3.3/families/Dfam_curatedonly.embl.gz")
    os.system("gunzip Dfam_curatedonly.embl.gz")
    os.chdir("../../")
"""

class emblParser():
    
    def __init__(self, fname=''):
        self.fname = fname
        
    def parse_entry(self, li):
        entry = {'id': li.split(';')[0].split(' ')[-1], 'seq': []}
        
        nm_pat = 'NM   '
        sq_pat = 'SQ   '
        os_pat = 'OS   '
        in_sequence = False
        for li in self.embl:
            if li[:2] == '//':
                entry['seq'] = ''.join(entry['seq'])
                return entry
            if in_sequence:
                entry['seq'].append(li.rstrip('\n').replace(' ', '').rstrip('0123456789'))
            if li[:len(nm_pat)] == nm_pat:
                entry['nm'] = li.rstrip('\n').split(' ')[-1]
            if li[:len(sq_pat)] == sq_pat:
                in_sequence = True
            if li[:len(os_pat)] == os_pat:
                entry['os'] = re.search(os_pat + '(.+)$', li).group(1).rstrip('\n')
    def parse(self, fname=''):
        if fname == '':
            fname = self.fname
        
        includes_humans = set([
            'Eutheria', 'Primates', 'Simiiformes', 'Mammalia', 'Vertebrata <vertebrates>', 
            'Theria <mammals>', 'Hominoidea', 'Hominidae', 'Homo sapiens', 'Homininae', 'Homo',
            'Metazoa',
        ])        
        
        self.entries = []
        self.embl = open(fname, 'r')
        
        id_pat = 'ID   '
        for li in self.embl:
            if li[:len(id_pat)] == id_pat:
                entry = self.parse_entry(li)
                if entry['os'] in includes_humans \
                and ('rRNA_' not in entry['nm']):
                    self.entries.append(entry)
        self.embl.close()
        
        entries_df = pandas.DataFrame(self.entries)
        print(entries_df)

        entries_df = entries_df.loc[[x in includes_humans for x in entries_df['os']], :]
        """
print(entries_df['os'].value_counts()) ->
Halyomorpha halys            1941  # Stinkbug
Danio rerio                  1740  # Fish
Heliconius                    538  # Butterfly
Eutheria                      460
Primates                      245
Mus <genus>                   219
Caenorhabditis elegans        180
Drosophila <flies,genus>      180
Muridae                       171
Danio                         149
Simiiformes                   121
Mus musculus                  121
Amniota                        97
Ficedula albicollis            91
Heliconiini                    76
Catarrhini                     71
Mammalia                       66
Vertebrata <vertebrates>       66
Uraeginthus cyanocephalus      66
Theria <mammals>               53
Euarchontoglires               33
Hominoidea                     29
Rodentia                       29
Drosophila melanogaster        25
Testudinoidea                  23
Murinae                        23
Hominidae                      21
Tetrapoda                      19
Homo sapiens                    8
Testudines                      6
Metazoa                         5
Durocryptodira                  5
Chrysemys                       5
Chrysochloris asiatica          5
Sauropsida                      4
Glires                          4
Homininae                       3
Haplorrhini                     3
Cyprinidae                      2
Boreoeutheria                   2
Protostomia                     2
Chrysemys picta bellii          1
Teleostei                       1
Sciuromorpha                    1
Drosophila simulans             1
Afrotheria                      1
Artiodactyla                    1
Homo                            1
Euteleostomi                    1
"""
        return entries_df
    
    def write_as_chromosome(
        self, outfasta='assets/reference/repeats_chrom.fa', outgtf='assets/reference/repeats.gtf'):
        """Example gtf line:
chr1	HAVANA	transcript	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; 
transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; 
gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; 
level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2";
havana_transcript "OTTHUMT00000362751.1";
"""
        fasta, gtf = open(outfasta, 'w'), open(outgtf, 'w')
        spacer_len = 100
        spacer = 'N' * spacer_len
        chrom_len = 0
        fasta.write('>repeats\n')
        for entry in self.entries:
            fasta.write(spacer)
            start = chrom_len + spacer_len
            end = chrom_len + spacer_len + len(entry['seq'])
            fasta.write(entry['seq'].upper())
            _id, nm = entry['id'], entry['nm']
            gtf.write("\t".join([
                'repeats', 'repeats', 'gene', str(start), str(end), '.', '+', '.',
                f'gene_id "{_id}"; transcript_id "{_id}"; gene_type "repeat"; gene_name "{nm}";']) + '\n')
            # Identical line for transcript except labelled transcript.
            gtf.write("\t".join([
                'repeats', 'repeats', 'transcript', str(start), str(end), '.', '+', '.',
                f'gene_id "{_id}"; transcript_id "{_id}"; gene_type "repeat"; gene_name "{nm}";']) + '\n')
            gtf.write("\t".join([
                'repeats', 'repeats', 'exon', str(start), str(end), '.', '+', '.',
                f'gene_id "{_id}"; transcript_id "{_id}"; gene_type "repeat"; gene_name "{nm}"; exon_number 1;']) + '\n')
            chrom_len = end
        gtf.close()
        fasta.close()
    
        
    def write_fasta_separate_entries(self, outfname=''):
        outf = open(outfname, 'w')
        for entry in self.entries:
            outf.write(f">{entry['nm']}-{entry['id']}\n{entry['seq']}\n")
        outf.close()

if __name__ == '__main__':
    fname = "assets/reference/Dfam_curatedonly.embl"
    ep = emblParser(fname)
    entries_df = ep.parse()
    ep.write_as_chromosome()

