from Bio import SeqIO
from sys import argv

blastout = argv[1]
refinp = argv[2]

reference = {}
for contig in SeqIO.parse(refinp, 'fasta'):
    reference[contig.id] = contig.seq

with open(blastout, 'r') as f:
    for line in f:
        name = line.split('\t')[0]
        chrom = line.split('\t')[1]
        val1 = int(line.split('\t')[8])
        val2 = int(line.split('\t')[9])
        corr = int(line.split('\t')[6])
        if val1 < val2:
            loc = val1 + 100 - corr
            #REF = reference[chrom][loc]
            REF = (reference[chrom][loc:loc+1]).reverse_complement() # To match PacBio REF version with genotypes
        else:
            loc = val1 - 100 - corr
            #REF = reference[chrom][loc]
            REF = (reference[chrom][loc:loc+1]).reverse_complement() # The same
        print(name, chrom, loc+1, REF, sep = '\t')
