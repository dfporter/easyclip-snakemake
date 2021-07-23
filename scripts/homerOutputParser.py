import re, pandas

class homerOutputParser():

    def parse(self, fname: str='homerResults.html') -> pandas.DataFrame:
        # Write the % of peaks with motifs.
        prefix = '/Users/dp/homer/used/top10Ktargets_5_min_height_no_Pvalue_cutoffs_20k_controlSeqs_mRNA_only_old_mapping_eg_repeats_first_Thursday/'
        a = open(fname).readlines()
        a = ''.join(a)

        n_pos = int(re.search('Total target sequences = (\d+)<BR/>', a).group(1))
        n_neg = int(re.search('Total background sequences = (\d+)<BR/>', a).group(1))

        order = "<TR><TD>Rank</TD><TD>Motif</TD><TD>P-value</TD><TD>log P-pvalue</TD><TD>% of Targets</TD><TD>% of Background</TD>".split('<TD>')
        order = [x.split('</TD>')[0] for x in order]
        _full_dict = {}
        for chunk in a.split('<TR>')[1:]:
            chunk = [x.split('</TD>')[0] for x in chunk.split('<TD>')]
            _dict = dict(zip(order, chunk))
            del _dict['Motif']
            if _dict['Rank'] != 'Rank':  # Skip header.
                if '<FONT color="red">*</FONT>' not in _dict['Rank']:  # Skip what Homer thinks are false positives (it marks with red *).
                    _dict['Rank'] = int(_dict['Rank'].rstrip('\n'))
                    del _dict['<TR>']  # Clean HTML marker.
                    _full_dict[_dict['Rank']] = _dict
                    
        df = pandas.DataFrame(_full_dict).T
        df['N pos'] = n_pos
        df['N background'] = n_neg

        return df