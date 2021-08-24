import sys, os, re

# Example line
#>GCUGAGAGUGUAGGAUGUUUACA        hsa-miR-30c MIMAT0000244 Homo sapiens miR-30c Targets (miRBase) 10.1864640240426

def is_human(name):
    if re.search('Homo.sapiens', name) is not None:
        return True
    return False

def is_miRNA(name):
    if 'hsa-miR' in name or 'hsa-let' in name:
        return True
    return False

if __name__ == '__main__':
    to_keep = {}
    keeping_current = False
    with open("data/external/known.rna.motifs") as f:
        for li in f:
            if li[0] == '>':
                if is_human(li) and not is_miRNA(li):
                    name = li.rstrip('\n')
                    to_keep[name] = []
                    to_keep[name].append(li)
                    keeping_current = True
                else:
                    keeping_current = False
            elif keeping_current:
                to_keep[name].append(li)
    with open("data/processed/filtered_homer_known.rna.motifs", 'w') as f:     
        for lines in to_keep.values():
            f.write(''.join(lines))
        
