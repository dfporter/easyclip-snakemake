import HTSeq, collections, pandas, os, re, importlib, csv
import numpy as np
from pathlib import Path
from typing import List, Mapping, Union

import scripts
import scripts.scheme
#import scripts.biotypeUtils
import scripts.biotypeLookupFileMaker
importlib.reload(scripts.biotypeLookupFileMaker)
importlib.reload(scripts.scheme)
#importlib.reload(scripts.biotypeUtils)

class readsPerGene():

    def __init__(
        self,
        filename: Union[str, None],  # The counts.txt file of reads per gene.
        scheme_filename: Union[str, None]=None,
        indexing_col_name: str='gene_name', index_col: Union[int, None]=0,
        exclude_unmapped: bool=True, ignore_biotypes: List[str]=[]):
        """Reads per gene counts, and methods to load to and process.

        The counts filename passed to readsPerGene is expected to be raw read counts,
        without normalization.

        counts.txt and ann_counts.txt files are always raw read counts.
        """
        self.log = f"Instantiated with {locals()}\n"

        self.data_folder = Path(os.path.dirname(filename), 'data')

        self.df = {}
        self.indexing_col_name = indexing_col_name
        self.filename = filename

        if filename.lower().endswith(('xls', 'xlsx')):
            self.df = pandas.read_excel(filename, index_col=index_col)
        else:
            self.df = pandas.read_csv(filename, sep='\t', index_col=index_col)

        self.scheme = None
        if scheme_filename is not None:
            self.scheme = scripts.scheme.scheme(scheme_filename)
            self.edit()

    def __add__(self, rpg_to_add):
        """Combine dataframe with another rpg object and return the combination."""
        self.df = self.df.merge(rpg_to_add.df, left_index=True, right_index=True, how='inner')
        if self.scheme is not None:
            self.scheme = self.scheme + rpg_to_add.scheme
        return self
    
    def edit(self, black_list=None, verbose=False, save=False,
        save_to=None, xl_rate_fname='percentCrosslinked.xlsx',
        mapping_fname='enst_transcript_id_name_biotype_map.txt',
        include_percent_xl=True):

        # Remove directories from column names.
        self.df.columns = [x.split('/')[-1].split('.+.wig')[0] for x in self.df.columns]

        # Remove unwanted datasets (appear contaminated, too small, ...).
        self.set_blacklist(black_list=black_list)
        self.df = self.df.loc[:, [x for x in self.df.columns if not any(b in x for b in self.black_list)]]

        # Remove columns we can't match to proteins in the scheme.
        if verbose: print(f"Columns before matching with scheme: {self.df.columns}")
        self.df = self.df.loc[:, [x for x in self.df.columns if (
            x in ['Gene type', 'gene_name'] or self.scheme.gene_from_fname(x))]]  # '' is False
        if verbose: print(f"After filtering: {self.df.columns}")
        self.df.fillna(0, inplace=True)

    def add_biotypes_column(
        self, mapping_fname='./temp/enst_transcript_id_name_biotype_map.txt',
        df=None, reload=False, gtf_filename='combined.gtf'):
        """
        Expect a mapping file with this format:
        gene_id   gene_name   gene_type
        ENST00000456328 DDX11L1 processed_transcript
        ENST00000450305 DDX11L1 transcribed_unprocessed_pseudogene
        ENST00000488147 WASH7P  unprocessed_pseudogene
        """
        
        #if (not os.path.exists(mapping_fname)) or reload:
            
        #    scripts.biotypeLookupFileMaker.make_enst_gene_name_biotype_map_file(
        #        gtf_filename=gtf_filename, out_filename=mapping_fname)

        # Create a lookup dictionary from the mapping text file.
        if (not hasattr(self, 'to_biotypes')) or reload:
            print(f"Loading {os.path.realpath(mapping_fname)} for biotype lookup.")
            to_biotypes_df = pandas.read_csv(mapping_fname, sep='\t', index_col=None)
            if 'gene_type' in to_biotypes_df.columns:
                self.to_biotypes = dict(zip(to_biotypes_df['gene_name'], to_biotypes_df['gene_type']))      
            elif ('transcript_biotype' in to_biotypes_df.columns):
                print(f"Warning: added biotype to genes from {mapping_fname} but found only transcript_biotype column.")
                self.to_biotypes = dict(zip(to_biotypes_df['gene_name'], to_biotypes_df['transcript_biotype']))
            else:
                raise ValueError(
                    f"Tried to add biotypes from {mapping_fname} but did not find expected column names *_biotype.")

        print("Adding biotypes.")

        def fix(name, biotype):
            if re.search('\ASNHG\d*::intron', name):
                biotype = 'snoRNA'
            if name=='_no_feature':
                biotype = 'Intergenic'
            if name=='_ambiguous':
                biotype = 'Target ambiguous'
            return biotype
            
        if df is None:
            df = self.df

        df['Gene type'] = [
            fix(x, self.to_biotypes.get(x.split('::')[0], 'Unknown')) for x in df.index]

        return df
        

    def proteins(self, df=None):

        if df is None:
            df = self.df

        _proteins = set([self.scheme.gene_from_fname(x) for x in df.columns])

        if not len(_proteins):
            print(f"No proteins found in columns: {df.columns}")

        return _proteins - set(['', None])

    def numeric_columns(self, df=None):

        if df is None:
            df = self.df

        return [x for x in df.columns if (df[x].dtype.kind in 'bifc')]

    def set_blacklist(self, black_list: Union[None, List[str]] = None) -> None:
        """Define column names in counts.txt/ann_counts.txt reads per gene files to remove.
        """
        if black_list is None:
            self.black_list = [
                'Exp61_HCT116_GTAGCC_TCA',  # <1,000 reads
                'Exp61_HCT116_GTAGCC_AGT',  # <10,000 reads
                'Exp15_AURKA_TCTGAG_TCA',  # <10,000 reads
                #'CCIN',
                #'EPB41L5',
                'Exp16_FBL_AGCTAG_CAG',  # Very small dataset.
                'Exp16_hnRNPC_TGAGTG_AGT',
                'Exp16_hnRNPC_TGAGTG_CAG',
                'Exp28_hnRNPC_CGATTA_AAC',
                'Exp31_UBA2_TGAGTG_AAC',# Too correlated with CDK4 AAC Exp31.
                'Exp31_UBA2_TGAGTG_CAG', # Empty.
                'Exp31_UBA2_TGAGTG_CAG',
                'Exp31_CDK4_GCCATG_AAC', # Too correlated with UBA2 AAC Exp31.
                'Exp31_CDK4_GCCATG_CAG', # Empty.
                'Exp33_CDK4_GCCATG_AGT',
                'Exp31_CDK4_GCCATG_CAG',
                'Exp33_CAPNS6_CACTGT_TCA', # Empty.
                'Exp61_PCBP1-100P_GCCATG_TCA',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1-100P_GCCATG_AGT',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1-100Q_AGCTAG_TCA',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1-100Q_AGCTAG_AGT',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1-dKH_ATCGTG_TCA',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1-dKH_ATCGTG_AGT',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1_CGATTA_AGT',  # Miseq run replaced by Hiseq (hp).
                'Exp61_PCBP1_CGATTA_TCA',  # Miseq run replaced by Hiseq (hp).
                'Exp00_Unknown',  # Like Null Island.
                'CSRP2',  # Might be a good dataset; just not using.
                '100Qox',  # Good quality data. Just no space in paper.
                '100Pox',  # Good quality data. Just no space in paper.
                'PCBP1ox',  # Good quality data. Just no space in paper.
                'Exp32_100P_CGAAAC_TCA',  # I believe this was OX.
                'YBX',
                'AURKA'
            ]
            self.black_list.extend([
                'Exp16_FBL_24AGCTAG_CAG',  # Very small dataset.
                'Exp16_hnRNPC_24TGAGTG_AGT',
                'Exp16_hnRNPC_24TGAGTG_CAG',
                'Exp28_hnRNPC_17CGATTA_AAC',
                'Exp31_UBA2_15TGAGTG_AAC',# Too correlated with CDK4 AAC Exp31.
                'Exp31_UBA2_15TGAGTG_CAG', # Empty.
                'Exp31_UBA2_05TGAGTG_CAG',
                'Exp31_CDK4_15GCCATG_AAC', # Too correlated with UBA2 AAC Exp31.
                'Exp31_CDK4_15GCCATG_CAG', # Empty.
                'Exp33_CDK4_17GCCATG_AGT',
                'Exp31_CDK4_05GCCATG_CAG',
                'Exp33_CAPNS6_17CACTGT_TCA', # Empty.
                'Exp61_PCBP1-100P_hpGCCATG_TCA',  # Too small.
                'Exp61_PCBP1-100P_hpGCCATG_AGT',  # Too small.
                'YBX',
                'AURKA'
            ])
        else:
            self.black_list = black_list

    def columns_for_a_protein(self, protein, df=None):
        if df is None:
            df = self.df
        return [x for x in self.numeric_columns(df) \
            if (self.scheme.gene_from_fname(x)==protein)]

    def average_by_protein(self) -> pandas.DataFrame:

        aves = pandas.DataFrame(index=self.df.index)
        for protein in self.proteins():
            cols = self.columns_for_a_protein(protein)
            aves[protein] = self.df[cols].mean(axis=1)

        return aves

class rawReadsPerGene(readsPerGene):
    """For if we want a readsPerGene object with a more explicit type name.
    """
    
    def foo(self):
        pass

class readsPerMillion(readsPerGene):
    """
    Expect to be used as:
    raw = readsPerGene('counts.txt', scheme='scheme.xlsx')
    rpm = readsPerMillion(raw)
    """

    def __init__(self,
        raw_reads_per_gene: Union[readsPerGene, rawReadsPerGene],
        #raw_counts_filename: Union[None, str] = None,
        use_bed_dir=None, load_total_read_numbers=False):

        self.raw_reads_per_gene = raw_reads_per_gene
        self.use_bed_dir = use_bed_dir
        self.load_total_read_numbers = load_total_read_numbers
        self.scheme = raw_reads_per_gene.scheme

        self.normalize_to_per_million_reads(raw_reads_per_gene.df)
        
    def __add__(self, rpm_to_add):
        """Combine dataframe with another rpm object and return the combination."""
        self.df = self.df.merge(rpm_to_add.df, left_index=True, right_index=True, how='inner')
        if self.scheme is not None:
            self.scheme = self.scheme + rpm_to_add.scheme
        return self
    
    def get_total_reads(self, folder_name, load=False):
        
        if load:            
            print(f"Loading total read counts from {load}.")
            with open(load) as fh:
                next(fh)  # Skip header.
                self.total_counts = dict((r[0],float(r[1])) for r in csv.reader(fh, delimiter='\t'))

            for protein in self.proteins(df=self.raw_reads_per_gene.df):
                cols = [x for x in self.total_counts if self.scheme.gene_from_fname(x)==protein]
                self.total_counts[protein] = sum([self.total_counts[x] for x in cols])

    def normalize_to_per_million_reads(self, df: pandas.DataFrame) -> None:
        """Set self.counts_per_million_df from self.raw_counts_df and
        a folder of bed files (for total read count).
        """
        # Make a copy, which we edit to be per million.
        df = df.copy()

        # Determine total read count for each replicate bed file.
        if self.use_bed_dir:  # Count the lines in a bed file.
            self.get_total_reads(self.use_bed_dir)#, load=self.load_total_read_numbers)
        elif self.load_total_read_numbers:  # Get AUC for a bedgraph.
            self.get_total_reads('', load=self.load_total_read_numbers)
        else:  # Just use the column totals.
            self.total_counts = df.loc[:, self.numeric_columns(df)].sum(axis=0).to_dict()

        # Remove columns with zero total reads.
        for column in [x for x in self.total_counts if x != 'gene_name']:
            if self.total_counts[column] == 0:
                print(f"No read counts in {column}. Removing.")
                df = df.loc[:, [x for x in df.columns if x!=column]]

        # Convert to reads per million.
        for col in self.numeric_columns(df):
            df[col] = (1E6/self.total_counts[col]) * df[col]

        self.df = df
        self.reads_per_million_df = self.df

class xlsPerProtein(readsPerGene):
    """
    Expect to be used as:
    raw = readsPerGene('counts.txt', scheme='scheme.xlsx')
    rpm = readsPerMillion(raw)
    xpp = xlsPerProtein(rpm, xl_rate_fname='percentCrosslinked.xlsx')
    """

    def __init__(
        self, reads_per_million: readsPerMillion,
        xl_rate_fname: str = 'percentCrosslinked.xlsx'):
        """Sets counts_per_protein_df from counts_per_million_df and xl_rate_fname.
        """

        self.scheme = reads_per_million.scheme
        self.df = reads_per_million.df.copy()
        
        print("Multiplying by XL rates and outputing rates as PER 1E6 reads, PER 10,000 proteins...")
        print("(so if XL rate=0.01%, 3 reads per million stays 3. XL=1% -> 3 reads per million becomes 300.)")
        print('Proteins:', self.proteins())

        for col in self.numeric_columns(self.df):
            protein = self.scheme.gene_from_fname(col)
            xl_rate = self.scheme.percent_xl_of_protein_from_file(
                    protein, xl_rate_fname=xl_rate_fname)
            #print(f'Multiplying by {xl_rate} for {protein}.')
            self.df[col] = 100. * float(xl_rate) * self.df[col].astype(float) 

        self.counts_per_protein_df = self.df
