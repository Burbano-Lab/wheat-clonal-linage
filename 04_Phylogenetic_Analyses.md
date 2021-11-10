# XX General Title XX
# 4. Phylogenetic Analyses

Program                | Location
---------------------- | --------------------------------------
*bcftools v.1.11*      | (https://github.com/samtools/bcftools)
*PLINK v.1.9*          | (https://www.cog-genomics.org/plink)
*tped2fasta*           | (https://github.com/smlatorreo/misc_tools)
*RAxML-NG v.1.0.3*     | (https://github.com/amkozlov/raxml-ng)
*ClonalFrameML v.1.12* | (https://github.com/xavierdidelot/clonalframeml)
*BactDating*           | (https://github.com/xavierdidelot/BactDating)
*R*                    | (https://cran.r-project.org/)

After removing those samples where recombination is present (see 3. Recombination analyses), we thinned the dataset to the B71 clonal lineag and the PY0925 cluster (also clonal) which can be used as an outgroup. We kept positions with no-missing data.

```bash
bcftools view -a -S B71clust_PY0925clust.list wheat-blast.snps.filtered.vcf.gz | bcftools view -m2 -M2 -g ^miss | bgzip > B71clust_PY0925clust.snps.filtered.fullinfo.vcf.gz
```

## Maximum-Likelihood (ML) phylogeny
We converted the VCF file into a pseudo-fasta format to have whole-genome concatenated SNPs per isoalte as a suitable input for the phylogenetic analyses.

```bash
plink --allow-extra-chr --vcf B71clust_PY0925clust.snps.filtered.fullinfo.vcf.gz --recode transpose --out B71clust_PY0925clust.snps.filtered.fullinfo
tped2fasta B71clust_PY0925clust.snps.filtered.fullinfo > B71clust_PY0925clust.snps.filtered.fullinfo.fasta

```

Then, we generated a ML phylogeny using RAxML-NG with a GTR+G substituion model and 1,000 bootstrap replicates.
```bash
raxml-ng --all --msa B71clust_PY0925clust.snps.filtered.fullinfo.fasta --msa-format FASTA --data-type DNA --model GTR+G --bs-trees 1000
```

## Removing homoplasy
As low bootstrap values can be caused due to presence homoplasies, we used *ClonalFrameML* in order to detect and remove those events from the calculated phylogeny.

```bash
ClonalFrameML B71clust_PY0925clust.snps.filtered.fullinfo.fasta.raxml.bestTree B71clust_PY0925clust.snps.filtered.fullinfo.fasta
```

