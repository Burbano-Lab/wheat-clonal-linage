# A pandemic clonal lineage of the wheat blast fungus
# 4. Recombination Analyses

Program                | Location
---------------------- | -----------
*PLINK v.1.9*          | (https://www.cog-genomics.org/plink)
*R*                    | (https://cran.r-project.org/)
*VCFtools v.0.1.16*    | (https://github.com/vcftools/vcftools)
*FFPopSim*             | (https://github.com/neherlab/ffpopsim)
*simulations.py*       | This repository
*RminCutter*           | (https://github.com/RILAB/rmin_cut)
*samtools v.1.11*      | (https://github.com/samtools/samtools)
*bcftools v.1.11*      | (https://github.com/samtools/bcftools)


### Inference of recombination via Linkage Desequilibrium (LD) decay patterns
First, we grouped the isolates in genetic groups based on the previously described PCA-based population structure analyses (e.g. B71 cluster, PY0925 cluster, other Brazilian isolates). For each genetic group we used *VCFtools* to compute several measures of LD. Since *VCFtools* is designed to handle diploid organisms, we transformed the haploid *VCF* files into "phased double haploid" *VCFs*
```bash
# Create a VCF as diplod
plink --allow-extra-chr --vcf wheat-blast.snps.filtered.vcf.gz --recode vcf \
--out wheat-blast.snps.filtered.as_diploid

# Artificially phase the VCF file
sed 's/\//\|/g' wheat-blast.snps.filtered.as_dip.vcf |
bgzip > wheat-blast.snps.filtered.as_dip_phased.vcf.gz

# Remove intermediate files
rm wheat-blast.snps.filtered.as_diploid.*
```

Next, using *VCFtools* we computed pairwise SNP correlations as *r<sup>2</sup>*, as well as Lewontin's *D* and *D'*
```bash
# Clusters of isolates were grouped in files named as "cluster_X.list"
vcftools --keep cluster_X.list --gzvcf wheat-blast.snps.filtered.as_dip_phased.vcf.gz \
--max-alleles 2 --min-alleles 2 --min-r2 0.1 --hap-r2 --phased --stdout |
gzip > cluster_X.LD.gz
```
Resulting files: [B71 cluster](/data/04_Recombination/B71_cluster.LD.gz) ; [PY0925 cluster](/data/04_Recombination/PY0925_cluster.LD.gz) ; [Brazilian samples](/data/04_Recombination/Brazilian_cluster.thinned.LD.gz).  
Finally, for each of the genetic clusters,  using *R*, we calculated the average of each LD measure (*r<sup>2</sup>*, Lewontin's *D* and *D'*) in bins of physical genomic distance.
```{r}
#R
bin_size <- 1000 # We used a bin size of 1000. Smaller sizes will result in a higher number of measures

# Load the dataset and transform the LD measures for their absolute values
cl <- read.table('cluster_X.LD.gz', header = TRUE)
cl$R.2 <- abs(cl$R.2)
cl$D <- abs(cl$D)
cl$Dprime <- abs(cl$Dprime)

# Restrict to distaces <= 2Mb
distances <- cl$POS2 - cl$POS1
cl <- cbind(cl, distances)
cl <- cl[(cl$distances <= 2000000), ]

# Create bins and empty vectors for final results
bins <- seq(1, max(cl$distances), bin_size)
r2_out <- c()
D_out <- c()
Dprime_out <- c()

# Iterate through bins and compute mean values for each of the LD measures
for(bin in bins){
    r2_mval <- mean(cl[(cl$distances >= bin) & (cl$distances < bin + bin_size), 5])
    r2_out <- c(r2_out, r2_mval)
    D_mval <- mean(cl[(cl$distances >= bin) & (cl$distances < bin + bin_size), 6])
    D_out <- c(D_out, D_mval)
    Dprime_mval <- mean(cl[(cl$distances >= bin) & (cl$distances < bin + bin_size), 7])
    Dprime_out <- c(Dprime_out, Dprime_mval)
}

# Generate plots
par(mfrow=c(1,3))
plot(bins, r2_out, main = 'r^2')
plot(bins, D_out, main = 'D')
plot(bins, Dprime_out, main = 'D prime')
```
![LD](/data/04_Recombination/LD.png)


## Forward simulations of LD decay patterns
To understand the impact of the probabability of sexual reproduction per generation and the population size on the patterns of LD decay, we used the Python interface of the forward simulator *FFPopSim*.

```python
# python2.7

# Load required libraries
import FFPopSim as h
import matplotlib.pyplot as plt
import numpy as np
import itertools

# Main function that accepts as input a python directory with the different parameters
def simul(params):
    L = params['L']
    N = params['N']
    ### set up
    pop = h.haploid_highd(L) # Use high dimensional function
    pop.set_wildtype(N)         # start with N wildtype individuals
    pop.mutation_rate = params['mu']
    pop.outcrossing_rate = params['r']
    pop.crossover_rate = params['cr']
    # Half of the population with genotype: 00000... and half with genotype: 111111...
    pop.set_genotypes([[0]*L, [1]*L],[N*0.5, N*0.5])
    pop.evolve(100) # The population will evolve for 100 generations

    # Iterate through all possible pairwise loci and measure LD
    LDs = []
    for pair in itertools.combinations(range(L), 2):
        l1 = pair[0]
        l2 = pair[1]
        diff = l2 - l1
        LDs.append([diff, abs(pop.get_LD(l1, l2))])
    LDs = np.array(LDs)

    # Compute median values per distace bin
    out = []
    for i in range(1,int(max(LDs[:,0]))+1):
        med = np.median(LDs[LDs[:,0] == i, 1])
        out.append([med])
    out = (LDs, np.array(out))
    return(out)

#   PARAMETERS
#   N: Population size
#   L: Number of loci
#   mu: mutation rate per site per generation
#   r: Probability of sexual reproduction per generation
#   cr: Crossover rate per site per generation
parameters = {'N':100000, 'L':300, 'mu':0.0003 , 'r':0.01, 'cr':0.01} # Example of set of parameters

# Run the simulation
result = simul(parameters)

# Generate plots
indx = np.array([[i] for i in range(1,parameters['L'])]) # X-axis
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(result[0][:,0], result[0][:,1], alpha=0.1, marker='.')
plt.xlabel('Distance between loci')
plt.ylabel('LD')
plt.show()
```
![LD simulation example](/data/04_Recombination/LD_simulation_example.png)  
A dedicated *Python* script can be found in the file [LD_simulations.py](/scripts/04_Recombination/LD_simulations.py)

## Detection of recombination events using the four-gamete test
We used the four-gamete test to detect recombination events in each of the defined clusters. For this purpose, we used *RminCutter*. Since this program takes alignments as input, we prepared the data using a combination of *samtools* and *bcftools*.

```bash
'''
The following snippet iterates throught 8 chromosomes of the 70-15 reference genome. 
Inside each loop, a new iteration over isolates of a cluster in a file 'cluster_N.list' is performed. 
Within each sub-loop, all the variants from each isolate are applied to the reference genome
and the outputs are pseudo-fasta files per chromosome
'''
for chr in {1..8}; do
    (while read sample; do
        samtools faidx 70-15.fasta $chr: |
	bcftools consensus -s $sample -p $sample\_ wheat-blast.snps.filtered.vcf.gz |
	tr "\n" " " | sed 's/ //g' | sed 's/$/\n/g' |
	sed 's/\:/\:\n/g' >> cluster_N.CHR$chr.fasta
        done < cluster_N.list )
done
```

We used the formatted pseudo-fasta files as input for *RminCutter*
```bash
# The following loop will run RminCutter on a given cluster of samples per chromosome
for chr in {1..8}; do
    perl RminCutter.pl -g -i $cluster.Chr$chr.fasta
done
```
The output consists of multiple fasta files that correspond to genomic regions without violations of the four-gamete test.

---
[Main README](/README.md) | [Previous - 03. Population Structure Analyses](/03_Population_Structure.md) | [Next - 05. Phylogenetic Analyses](/05_Phylogenetic_Analyses.md)
