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

We used *blastn* to align the the 84 SNPs and their sorrounding regions to the B71 reference genome and subsequently identify the exact genomic coordinates of all 84 SNPs.
First, we indexed the B71 reference genome with *makeblastdb*
```bash
makeblastdb -dbtype nucl -in B71.fa -out B71

```
We then used *blastn* to align the [regions around the diagnostic SNPs]([/data/01_Monsterplex_84SNPs_analyses/Diagnostic_SNPs_regions.fasta) to the B71 reference genome
```bash
blastn -db B71 -query Diagnostic_SNPs_regions.fasta -outfmt 6 > diagnostic_SNPs_regions_against_B71.out
```
We used the python custom script *extract_positions_from_blast_result.py* to extract the genomic coordintes of the B71 reference genome corresponding to each diagnostic SNP
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

---
[Main README](/README.md) | [Next - 02. Preprocessing and Variant Calling Analyses](/02_Preprocessing_and_Variant_Calling.md)
