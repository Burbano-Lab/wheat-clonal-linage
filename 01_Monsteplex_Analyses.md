# General Title
# 1. Monsterplex analyses

An impotant source of genetic information come from rapid screening analyses on the fileds, from which few SNPs can be easily assessed with the Monsterplex methodology.

To show how few SNP positions (N = 84) can discriminate the most important infecting lineages in *Pyricularia oryzae*, we gathered whole genomic information from 576 samples and reconstructed a quick phylogenetic analyses based only on those 84 positions.



Program                                  | Location
---------------------------------------- | --------------------------------------------------
*AdapterRemoval v2*                      | (https://github.com/mikkelschubert/adapterremoval)
*BWA v.0.7.12*                           | (https://github.com/lh3/bwa)
*samtools v.1.11*                        | (https://github.com/samtools/samtools)
*blast+ v.2.12*                          | (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
*extract_positions_from_blast_result.py* | (./scripts/01_Monsterplex_84SNPs_analyses/extract_positions_from_blast_result.py)
*extract_genotypes_from_BAM.sh*          | (./scripts/01_Monsterplex_84SNPs_analyses/extract_genotypes_from_BAM.sh)

## Alignment of short reads to reference genome

Raw .fastq sequences were trimmed with *AdapterRemoval2*.
```bash
AdapterRemoval --file1 $sample1.R1.fastq.gz --file2 $sample1.R2.fastq.gz --gzip --basename $sample.trimmed
```

We used the wheat-infecting *Magnaporthe oryzae* B71 assembly as the reference genome and used *BWA* to create the index.
```bash
bwa index B71.fa
```

*BWA mem* was used to map the trimmed reads to the reference genome and *samtools* was used to discard non-mapped reads, sort and convert the mapped file into a binary format (BAM).
```bash
bwa mem B71.fa $sample1.trimmed.R1.fastq.gz $sample1.trimmed.R1.fastq.gz > sample1_mapped.sam
samtools view -F 4 -Sbh sample1_mapped.sam | samtools sort - > sample1_mapped.bam
samtools index sample1_mapped.bam
```

## Finding the position of the diagnostic SNPs in the B71 alignment

The general approach is to use *blastn* to align the regions where the SNPs are located to the B71 alignment and extract the exact coordinates.
First, we index the B71 with *makeblastdb*
```bash
makeblastdb -dbtype nucl -in B71.fa -out B71

```
We then use *blastn* to align the [regions around the diagnostic SNPs]([./data/01_Monsterplex_84SNPs_analyses/Diagnostic_SNPs_regions.fasta) to the B71 reference genome
```bash
blastn -db B71 -query Diagnostic_SNPs_regions.fasta -outfmt 6 > diagnostic_SNPs_regions_against_B71.out
```
With the python custom script *extract_positions_from_blast_result.py* we now extract the positions of the B71 reference genome associated with each diagnostic SNP
```bash
python3 extract_positions_from_blast_result.py diagnostic_SNPs_regions_against_B71.out B71.fa > positions_diagnostic_SNPs.out
```

## Extracting genotypes per sample
Finally, we extract the genotypes on each of the mapped BAM files with a custom bash script *extract_genotypes_from_BAM.sh* which uses samtools mpileup to sample the genotypes at specific positions
```bash
bash extract_genotypes_from_BAM.sh sample1_mapped.bam B71.fa positions_diagnostic_SNPs.out
```
The generaed format is a fasta-like file which can be used for further phylogenetic analyses
