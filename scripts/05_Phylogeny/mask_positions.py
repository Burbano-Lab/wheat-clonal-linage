# USAGE: python3 mask_positions.py <input.fasta> <output_clonal_frame.importation_status.txt>
from sys import argv

fasta = argv[1]
to_mask = argv[2]

# Create dictionary with isolate names as keys and empty lists as values
mask = {}
with open(to_mask, 'r') as f:
    f.readlines(1)
    for line in f:
        if line.split('\t')[0] in mask.keys():
            continue
        else:
            mask[line.split('\t')[0]] = []
# Store masking targets per isolate
with open(to_mask, 'r') as f:
    f.readlines(1)
    for line in f:
        isolate, start, end = line.rsplit()
        mask[isolate] = mask[isolate] + [i-1 for i in range(int(start), int(end)+1)]
# Iterate through the fasta by replacing the ,asking targets with '?'
with open(fasta, 'r') as f:
    for line in f:
        print(line.strip())
        isol = line.split('>')[1].strip()
        line = list(f.readlines(1)[0].strip())
        for i in mask[isol]:
            line[i] = '?'
        print(''.join(line))
