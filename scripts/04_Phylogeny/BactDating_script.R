#### Importing from ClonalFrame-ML ####
library(ape)
library(scales)
setwd('~/Work/wheat_moryzae/clonal_frame_ML/mapped_to_70-15/')


# Temporal Signal #
t=read.tree('B71_and_PY0925_clust.snps.filtered.fullinfo.labelled_tree.newick')

distances <- cophenetic(t)
d_to_PY0925 <- distances[colnames(distances) == 'PY0925', ]

# Just keep distances of the B71 lineage (remove those from isolates: '053i','PY0925','117','37','12.1.037')
dist_to_PY0925 <- dist_to_PY0925[! names(d_to_PY0925) %in% c('053i','PY0925','117','37','12.1.037')]




dt <- read.table('~/Work/wheat_moryzae/ML-trees/mapped_to_70-15/B71_and_PY0925_clust.dates', header = F)
dts <- c()
for(n in names(d_to_PY0925)){
	  dts <- c(dts, dt[dt[,1] == n, 2])
}

m  <- data.frame(Coll_Year=dts, Patr_dist_to_PY0925=d_to_PY0925)

plot(m)
abline(lm(m$Distance ~ m$Date))
cor.test(m$Coll_Year, m$Patr_dist_to_PY0925, method = 'pear')

