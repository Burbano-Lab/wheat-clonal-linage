# A pandemic clonal lineage of the wheat blast fungus
# 3. Population Structure Analyses

Program                  | Location
------------------------ | ----------------------------
*PLINK v.1.9*            | (https://www.cog-genomics.org/plink)
*R*                      | (https://cran.r-project.org/)


## PCA analyses
We computed pairwise Hamming distances on non-missing SNP positions (full information).
```bash
plink --allow-extra-chr --vcf wheat-blast.snps.filtered.vcf.gz --out wheat-blast.snps.filtered
```

The resulting files ([wheat-blast.snps.filtered.dist](/data/03_Population_Structure/wheat-blast.snps.filtered.dist) ; [wheat-blast.snps.filtered.dist.id](/data/03_Population_Structure/wheat-blast.snps.filtered.dist.id)) were used as input to performed a PCA using *R*
```R
# R
m <- read.table('wheat-blast.snps.filtered.dist', header = FALSE)
npca <- prcomp(m, scale. = TRUE)
plot(npca$x[,1], npca$x[,2], xlab = 'PC1', ylab = 'PC2')
```
![Wheat blast PCA](/data/03_Population_Structure/Wheat_blast_PCA.png)

---
[Main README](/README.md) | [Previous - 02. Preprocessing and Variant Calling Analyses](/02_Preprocessing_and_Variant_Calling.md) | [Next - 04. Recombination Analyses](/04_Recombination_Analyses.md)
