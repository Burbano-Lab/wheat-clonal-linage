# A pandemic clonal lineage of the wheat blast fungus
# 6. Mating Type analyses

Program                         | Location
------------------------------- | --------------------------------------
*Bwa-mem2 v.2.1*                | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*               | (https://github.com/samtools/samtools)

The mating type of each isolate was classified using the breadth of coverage of the MAT1-1-1 and MAT1-2-1 loci as proxy.
We used the nucleotides codifying for MAT1-1-1 (GenBank: BAC65091.1) and MAT1-2-1 (GenBank: BAC65094.1) to create a [fasta file used as reference](/data/06_Mating_Type/MAT_genes.fasta). We generated an index using  *Bwa-mem2*.

```bash
bwa-mem2 index MAT_genes.fasta
```
Per-isolate sorted BAM files were created, following the same steps described in [section 2](/02_Preprocessing_and_Variant_Calling.md). Briefly:
```bash
bwa-mem2 mem MAT_genes.fasta $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R2.fastq.gz > sample1.sam
samtools view -SbhF 4 > sample1_mapped.bam
sambamba sort -o sample1_mapped_sorted.bam sample1_mapped.bam
sambamba markdup sample1_mapped_sorted.bam sample1_mapped_sorted.dd.bam
```

Finally, we used a combination of *samtools view* and  *samtools depth* to compute the breadth of coverage as well as the average depth for both MAT1-1-1 and MAT1-2-1 per BAM file.
```bash
# Breadth of coverage
samtools view -bh -q 30 sample1_mapped_sorted.dd.bam MAT1-1-1: | # Filter all reads mapped to MAT1-1-1 with Mapping Quality > 30
samtools depth -a - | # Print out the depth per position
awk '{if($3 > 0){sum += 1}}END{print (sum / 984)*100}' # Filter positions covered at least with 1X and divide the count over 984 (length of MAT1-1-1)

samtools view -bh -q 30 sample1_mapped_sorted.dd.bam MAT1-2-1: | # Filter all reads mapped to MAT1-2-1 with Mapping Quality > 30
samtools depth -a - | # Print out the depth per position
awk '{if($3 > 0){sum += 1}}END{print (sum / 1330)*100}' # Filter positions covered at least with 1X and divide the count over 1330 (length of MAT1-2-1)

```
A [summary of the computed breadth of coverage for MAT1-1-1 and MAT1-2-1](/data/06_Mating_Type/breadth_coverage_MAT.txt) allowed to classify the mating type per isolate

---
[Main README](/README.md) | [Previous - 05. Phylogenetic Analyses](/05_Phylogenetic_Analyses.md) | [Next - 07. AVR Rmg8 Analyses](/07_AVR_Rmg8.md)
