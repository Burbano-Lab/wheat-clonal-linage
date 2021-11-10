# XX General Title XX
# 2. Population Structure Analyses

Program                  | Location
------------------------ | ----------------------------
*AdapterRemoval v2*      | (https://github.com/mikkelschubert/adapterremoval)
*Bwa-mem2 v.2.1*         | (https://github.com/bwa-mem2/bwa-mem2)
*samtools v.1.11*        | (https://github.com/samtools/samtools)
*sambamba v0.8.0*        | (https://github.com/biod/sambamba)
*GATK v4.2*              | (https://github.com/broadinstitute/gatk/releases)
*PLINK v.1.9*            | (https://www.cog-genomics.org/plink)
*FFPopSim*               | (https://github.com/neherlab/ffpopsim)
*R*                      | (https://cran.r-project.org/)

## Alignment of short reads to reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the rice-infecting *Magnaporthe oryzae* 70-15 assembly as the reference genome and used *Bwa-mem2* to create the index.
```bash
bwa index B71.fa
```

*BWA mem2* was used to map the trimmed reads to the reference genome, *samtools* was used to to discard non-mapped reads, and *sambamba* was used to sort and mark PCR optical duplicates.
```bash
bwa-mem2 mem -R "@RG\tID:$sample\tSM:$sample" 70-15.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R2.fastq.gz > sample1.sam
samtools view -SbhF 4 > sample1_mapped.bam
sambamba sort -o sample1_mapped_sorted.bam sample1_mapped.bam
sambamba markdup sample1_mapped_sorted.bam sample1_mapped_sorted.dd.bam
```

## Genotype calling
We used *HaplotypeCaller* from *GATK* to generate genomic haplotype calls per individual using the duplicate-marked resulting BAM file as input
```bash
gatk HaplotypeCaller -R 70-15.fa -I sample1_mapped_sorted.dd.bam -O sample1.g.vcf.gz
```

We used *CombineGVCFs*, *GenotypeGVCFs* and *SelectVariants* from *GATK* to combine the individual genomic VCFs, call genotypes and filter for SNPs, respectively
```bash
gatk CombineGVCFs -R 70-15.fa -V sample1.g.vcf.gz -V sample2.g.vcf.gz -V sampleN.g.vcf.gz -O wheat-blast.g.vcf.gz
gatk GenotypeGVCFs -R 70-15.fa -ploidy 1 -V wheat-blast.g.vcf.gz -O wheat-blast.raw.vcf.gz
gatk SelectVariants -select-type SNP -V wheat-blast.raw.vcf.gz -O wheat-blast.raw.snps.vcf.gz
```

We extracted all Quality by Depth values
```bash
bcfools view -H wheat-blast.raw.snps.vcf.gz | cut -f8 |
awk -F "QD=" '{print $2}' | cut -f1 -d ";" | gzip >  wheat-blast.raw.snps.QD.gz
```

We assesed the distribution of the Quality by Depth, to set a filters of one standard deviation around the median value
```python
# Python
import pandas as pd
QD = pd.read_csv('wheat-blast.raw.snps.QD.gz', header = None, compression = 'gzip')
med = QD.median()
lower = med - QD.std()
upper = med + QD.std()
print(lower, upper)
```

Finally, we filtered accordingly using *GATK VariantFiltration* and created a new VCF keeping non-missing positions using *bcftools*
```bash
gatk VariantFiltration --filter-name "QD" \
--filter-expression "QD <= $lower || QD >= $upper" \
-V wheat-blast.raw.snps.QD.gz \
-O wheat-blast.snps.filter.vcf.gz

bcftools view -g ^miss wheat-blast.snps.filter.vcf.gz | bgzip > wheat-blast.snps.filtered.vcf.gz
```

## PCA analyses
We computed pairwise Hamming distances on non-missing SNP positions.
```bash
plink --allow-extra-chr --vcf wheat-blast.snps.filtered.vcf.gz --out wheat-blast.snps.filtered
```

And performed a PCA using *R*
```R
# R
m <- read.table('wheat-blast.snps.filtered.dist', header = False)
npca <- prcomp(m, scale.=True)
plot(npca$x[,1], npca$x[,2])
```

