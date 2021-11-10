# The input file is the ClonalFrameML output '_importation_status.txt' (better if NODE* rows are removed) and the original .fasta alignment.
# The output is splitted: 1) STDOUT: Cleaned fasta alignment 2) STDERR: fasta alignment with homoplasies
# USAGE: python3 clean_homoplasy_from_fasta.py <file.importation_status.txt> <original.fasta> > cleaned.fasta 2> homoplasies.fasta

import pandas as pd
from sys import argv, stderr

homopl_regions = argv[1]
fasta = argv[2]

m = pd.read_csv(homopl_regions, delimiter='\t', header=0)
l = max(m.End)

remove = set(range(1,l+1))
for row in range(m.shape[0]):
    lower = m.iloc[row,1]
    upper = m.iloc[row,2]+1
    k = remove.difference(set(range(lower,upper)))
    remove = k
total = set(range(1,l+1))
keep = total.difference(remove)

with open(fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            print(line.strip())
            print(line.strip(), file = stderr)
        else:
            seq = line.strip()
            print(''.join(seq[i-1] for i in keep)) # i-1 to match python's 0-based coordiate system
            print(''.join([seq[i-1] for i in remove]), file = stderr) # Same

