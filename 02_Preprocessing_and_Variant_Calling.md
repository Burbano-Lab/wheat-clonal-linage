# A pandemic clonal lineage of the wheat blast fungus
# 2. Preprocessing and Variant Calling

Program                  | Location
------------------------ | ----------------------------
*AdapterRemoval v2*      | (https://github.com/mikkelschubert/adapterremoval)
*Bwa-mem2 v.2.1*         | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*        | (https://github.com/samtools/samtools)
*sambamba v0.8.0*        | (https://github.com/biod/sambamba)
*GATK v4.2*              | (https://github.com/broadinstitute/gatk/releases)
*bcftools v.1.11*        | (https://github.com/samtools/bcftools)

## Alignment of short reads to the rice-infecting *M. oryzae* reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the rice-infecting *Magnaporthe oryzae* 70-15 assembly as the reference genome and indexed this genome using *Bwa-mem2*.
```bash
bwa index 70-15.fa
```

*BWA mem2* was used to map the trimmed reads to the reference genome, *samtools* to discard non-mapped reads, and *sambamba* to sort and mark PCR optical duplicates.
```bash
bwa-mem2 mem -R "@RG\tID:$sample\tSM:$sample" 70-15.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R2.fastq.gz > sample1.sam
samtools view -SbhF 4 > sample1_mapped.bam
sambamba sort -o sample1_mapped_sorted.bam sample1_mapped.bam
sambamba markdup sample1_mapped_sorted.bam sample1_mapped_sorted.dd.bam
```

Depth and coverage statistics were computed using samtools depth
```bash
# Average depth
samtools depth sample1_mapped_sorted.dd.bam | awk '{sum += $3}END{print sum / NR}'

# Percentage of genome covered. 70-15 genome reference length = 41027733
samtools depth sample1_mapped_sorted.dd.bam | awk '{NR / 41027733}'

```
A summary of these statistics was organized in a [table](/data/02_Preprocessing_and_Variant_Calling/coverage_summary.tsv)


## Variant calling
We used the *HaplotypeCaller* from *GATK* to generate genomic haplotype calls per individual using the duplicate-marked BAM file as input.
```bash
gatk HaplotypeCaller -R 70-15.fa -I sample1_mapped_sorted.dd.bam -O sample1.g.vcf.gz
```

We used *CombineGVCFs*, *GenotypeGVCFs* and *SelectVariants* from *GATK* to combine the individual genomic VCFs, call genotypes and filter SNPs, respectively.
```bash
gatk CombineGVCFs -R 70-15.fa -V sample1.g.vcf.gz -V sample2.g.vcf.gz -V sampleN.g.vcf.gz -O wheat-blast.g.vcf.gz
gatk GenotypeGVCFs -R 70-15.fa -ploidy 1 -V wheat-blast.g.vcf.gz -O wheat-blast.raw.vcf.gz
gatk SelectVariants -select-type SNP -V wheat-blast.raw.vcf.gz -O wheat-blast.raw.snps.vcf.gz
```

We extracted all Quality-by-Depth (QD) values
```bash
bcftools view -H wheat-blast.raw.snps.vcf.gz | cut -f8 |
awk -F "QD=" '{print $2}' | cut -f1 -d ";" | gzip >  wheat-blast.raw.snps.QD.gz
```

Based on the distribution of Quality-by-Depth values, we set filters of one standard deviation around the median value.
```python
# Python
import pandas as pd
QD = pd.read_csv('wheat-blast.raw.snps.QD.gz', header = None, compression = 'gzip')
med = QD.median()
lower = med - QD.std()
upper = med + QD.std()
print(lower, upper)
```

Finally, using the above-mentioned scheme, we filtered SNPs using *GATK VariantFiltration* and created a new VCF file, keeping non-missing positions, using *bcftools*.
```bash
gatk VariantFiltration --filter-name "QD" \
--filter-expression "QD <= $lower || QD >= $upper" \
-V wheat-blast.raw.snps.QD.gz \
-O wheat-blast.snps.filter.vcf.gz

bcftools view -g ^miss wheat-blast.snps.filter.vcf.gz | bgzip > wheat-blast.snps.filtered.vcf.gz
```

---
[Main README](/README.md) | [Previous - 01. Monsterplex Analyses](/01_Monsteplex_Analyses.md) | [Next - 03. Population Structure Analyses](/03_Population_Structure.md)
