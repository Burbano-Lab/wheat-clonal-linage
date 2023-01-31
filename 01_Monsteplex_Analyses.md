# A pandemic clonal lineage of the wheat blast fungus
# 1. Phylogenetic placement of wheat-infecting isolates from the Zambian wheat blast outbreak (2018-2020) using 84 SNPs

To establish the genetic makeup and the phylogenetic placement of the wheat-infecting lineage that caused a wheat blast outbreak in Zambia (2018-2020) we analyzed a set of 84 SNPs, which were designed  to distinguish  between  the  pandemic  clonal  lineage  of  the  wheat  blast fungus (*Magnaporthe oryzae*) that reached South East Asia in 2016 from other *M. oryzae* genotypes. This set of 84 SNPs was: (i) genotyped by multiplex-PCR (Monsterplex) derived from blast fungus  isolates (N=186) from 13 different grasses; and (ii) extracted from whole-genome sequecing data derired from *M. oryzae* isolates (N=XXX) from XX different grasses. The total dataset consisted of 576 *M. oryzae* isolates. 

Program                                  | Location
---------------------------------------- | --------------------------------------------------
*AdapterRemoval v2*                      | (https://github.com/mikkelschubert/adapterremoval)
*BWA v.0.7.12*                           | (https://github.com/lh3/bwa)
*samtools v.1.11*                        | (https://github.com/samtools/samtools)
*blast+ v.2.12*                          | (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
*extract_positions_from_blast_result.py* | (/scripts/01_Monsterplex_84SNPs_analyses/extract_positions_from_blast_result.py)
*extract_genotypes_from_BAM.sh*          | (/scripts/01_Monsterplex_84SNPs_analyses/extract_genotypes_from_BAM.sh)

## Alignment of short reads to the *M. oryzae* B71 reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*.
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the wheat-infecting *Magnaporthe oryzae* B71 assembly as the reference genome and indexed this genome using *BWA*.
```bash
bwa index B71.fa
```

*BWA mem* was used to map the trimmed reads to the reference genome, and *samtools* to discard non-mapped reads, sort and convert the mapped file into a binary format (BAM).
```bash
bwa mem B71.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R1.fastq.gz > sample1_mapped.sam
samtools view -F 4 -Sbh sample1_mapped.sam | samtools sort - > sample1_mapped.bam
samtools index sample1_mapped.bam
```

## Identification of the genomic location of the 84 diagnostic SNPs in the *M. oryzae* B71 reference genome

We used *blastn* to align the the 84 SNPs and their surrounding regions to the B71 reference genome and subsequently identify the exact genomic coordinates of all 84 SNPs.
First, we indexed the B71 reference genome with *makeblastdb*
```bash
makeblastdb -dbtype nucl -in B71.fa -out B71

```
We then used *blastn* to align the [regions around the diagnostic SNPs]([/data/01_Monsterplex_84SNPs_analyses/Diagnostic_SNPs_regions.fasta) to the B71 reference genome
```bash
blastn -db B71 -query Diagnostic_SNPs_regions.fasta -outfmt 6 > diagnostic_SNPs_regions_against_B71.out
```
We used the python custom script *extract_positions_from_blast_result.py* to extract the genomic coordinates of the B71 reference genome corresponding to each diagnostic SNP
```bash
python3 extract_positions_from_blast_result.py diagnostic_SNPs_regions_against_B71.out B71.fa > positions_diagnostic_SNPs.out
```

## Extraction of genotypes for each *M. oryzae* isolate
Finally, we extracted the genotypes for each of the per isolate mapped BAM files with a custom bash script *extract_genotypes_from_BAM.sh*, which uses samtools mpileup to retrieve the genotypes at specific genomic locations
```bash
bash extract_genotypes_from_BAM.sh sample1_mapped.bam B71.fa positions_diagnostic_SNPs.out
```
The output format is a [fasta-like file](/data/01_Monsterplex_84SNPs_analyses/all.maxmiss0.05.fasta) that can be used for subsequent phylogenetic analyses.

## Phylogenetic placement of wheat-infecting isolates
We generated a [Neighbor-Joining](/data/01_Monsterplex_84SNPs_analyses/NJ_all.maxmiss0.05_Bootstrap1K.nexus) tree of 532 worldwide distributed *M. oryzae* isolates based on the 84 concatenated SNPs.

## Benchmark of the 84 set of SNPs
To show that the set of 84 SNPs are informative, we subset the common isolates for which there is whole genomic data available. We computed a [WGS-based Neighbor-Joining tree](/data/01_Monsterplex_84SNPs_analyses/WGS_Benchmark_tree.newick) and [a pruned tree from the 84 SNPs based tree with the common isolates](/data/01_Monsterplex_84SNPs_analyses/84SNPs_Benchmark_tree.newick).

We further calculated genetic (Hamming) distances between each pair of blast isolates using the set of [84 SNPs](/data/01_Monsterplex_84SNPs_analyses/Hamming_Distances_Benchmark_84SNPs.tsv), and the [genome-wide SNPs](/data/01_Monsterplex_84SNPs_analyses/Hamming_Distances_Benchmark_WGS.tsv)

```R
# R

# Load the pairwise matrix
SNPsdist <- read.table('Hamming_Distances_Benchmark_84SNPs.tsv', row.names = 1, header = T)
WGSdist <- read.table('Hamming_Distances_Benchmark_WGS.tsv', row.names = 1, header = T)

# Normalize Hamming Distances
SNPsdist <- SNPsdist / max(SNPsdist)
WGSdist <- WGSdist / max(WGSdist)

# Subset the upper triangles
SNPsvals <- SNPsdist[upper.tri(SNPsdist)]
WGSsvals <- WGSdist[upper.tri(WGSdist)]

# Produce plot and compute correlation
plot(SNPsvals, WGSsvals, xlab = 'Monsterplex SNP dataset', ylab = 'WGS SNP dataset', main = 'Pairwise Hamming Distances')
reg <- lm(WGSsvals ~ SNPsvals)
abline(reg, lty = 2)
corrcoef <- cor.test(SNPsvals, WGSsvals, method = 'spe')
mtext(paste("Spearman's p = ", round(corrcoef$estimate, 2), '\tpval = ', corrcoef$p.value, sep = ''), 3, at = 0.2)
```
![WGS_Vs_Monsterplex](/data/01_Monsterplex_84SNPs_analyses/WGS_Vs_Monsterplex.png)

Finally, we created a random sampled distribution and compared with a permuted distribution to assess the robustness of the previous estimation
```R
# R

set.seed(12345)
corrcoef <- c()
corrcoefpermut <- c()
for(i in 1:100){
  # Random subsample
  subs <- sample(length(SNPsvals), length(SNPsvals)*0.1)
  SNPsubset <- SNPsvals[subs]
  WGSsubset <- WGSsvals[subs]
  r2 <- cor(SNPsubset, WGSsubset, method = 'spe')
  corrcoef <- c(corrcoef, r2)
  
  # Permutation
  SNPsubsetPermut <- sample(SNPsubset, length(SNPsubset))
  r2 <- cor(SNPsubsetPermut, WGSsubset, method = 'spe')
  corrcoefpermut <- c(corrcoefpermut, r2)
}

# Plot
boxplot(list(corrcoef, corrcoefpermut), names = c('Subsampling','Permutation'), ylab = "Spearman's p")
```
![Subsampling_and_Permutation](/data/01_Monsterplex_84SNPs_analyses/WGS_Vs_Monsterplex_Subsampling_and_Permutation.png)

---
[Main README](/README.md) | [Next - 02. Preprocessing and Variant Calling Analyses](/02_Preprocessing_and_Variant_Calling.md)
