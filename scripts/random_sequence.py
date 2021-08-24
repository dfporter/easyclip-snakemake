import random

def write_random_seqs(seqs, outfname, minimum=10000):

    alphabet = ['A', 'T', 'C', 'G']

    n_seqs = len(seqs) if len(seqs) > minimum else minimum

    if len(seqs) == 0:
        return ''

    with open(outfname, 'w') as outfh:
        for n, seq in enumerate(seqs):
            outlines = ">{}\n".format(n)
            outlines += "".join(random.choices(alphabet, k=len(seq))) + '\n'
            outfh.write(outlines)

    return outlines
