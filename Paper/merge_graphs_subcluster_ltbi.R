# File: merge_graphs_subcluster_ltbi.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 23/11/2016
# Desc: generate subclusters in the ltbi dataset after merging two datasets


######### load libraries and set global variables
library(org.Hs.eg.db)
library(downloader)
source('CGraphClust.R')
# plotting parameters
p.old = par()

## load reactome data
### load the uniprot2reactome mapping obtained from
# http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
# convert enterez ids to uniprot as Reactome database file uses UNIPROT ids
dfMap = AnnotationDbi::select(org.Hs.eg.db, dfReactome$V1, 'ENTREZID', 'UNIPROT')
dfMap = na.omit(dfMap)

## map reactome ids to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
rm(dfReactome, dfMap)
gc()
###

########### utility functions
get.reactome.name = function(csNames){
  i = which(dfReactome.sub$V2 %in% csNames)
  dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
  dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
  rownames(dfCluster.name) = NULL
  return(dfCluster.name)
}

## load the two tb datasets to merge
load('Paper/Results/ltb_atb_data.rds')
load('Paper/Results/tb_data.rds')

## genes that show association in both graphs are selected by
## intersecting the two graphs
# intersect the 2 graphs and recalibrate the graph object with the common vertices
igi = graph.intersection(getFinalGraph(ltb_atb_data$graph), getFinalGraph(tb_data$graph))
d = degree(igi)
igi = delete.vertices(igi, which(d == 0))
ecount(igi)
vcount(igi)
## ~ 224 genes are coexpressed in both the datasets
cVertID = V(igi)$name
## remake the graph using the ltb_atb graph
oGr = CGraphClust.recalibrate(ltb_atb_data$graph, cVertID)

mCounts = ltb_atb_data$matrix
fGroups = ltb_atb_data$groups

## general graph structure
## we would like to see how does the graph look like, are the clusters connected or in subgraphs
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))

# get the expression matrix for the clusters
mMarginal = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters

csClust = rownames(mMarginal)
length(csClust)
pdf('Paper/Results/cluster_variance_merged_dataset.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, mMarginal, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(mMarginal[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())


# plot pca biplot for this data using posterior mean
set.seed(1)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
pdf('Paper/Results/pca_clusters_merged_dataset.pdf')
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
dev.off(dev.cur())

## as we have reduced dimensions we try and cluster ltbi into 2 groups
## this is based on biological knowledge that LTBI can be Recently exposed and remotely exposed
## and LTBI is a more heterogenous group and some clusters e.g.
## 1280215 and 372790 which mostly contain genes from Cytokine signalling and GPCR signalling pathways
## have a higher posterior variance when compared with Healthy and active tb groups.
## subcluster the ltb group into two subgroups
ivLtb = which(fGroups == 'LTB')
mLtb = mMarginal[,ivLtb]
mLtb = t(mLtb)
set.seed(123)
# cluster along the row vectors
km.out = kmeans(mLtb, centers = 2, nstart = 20)
km.out
table(km.out$cluster)
fGroups = as.character(fGroups)
fGroups[ivLtb] = paste0(fGroups[ivLtb], km.out$cluster)
fGroups = factor(fGroups, levels = c('HC', 'LTB1', 'LTB2', 'ATB'))
names(fGroups) = names(ltb_atb_data$groups)
# pca with new groupings
set.seed(1)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
pdf('Paper/Results/pca_clusters_merged_dataset_after_LTBI_subcluster.pdf')
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
dev.off(dev.cur())

# reorder the groups and the matrix
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]
rownames(mCounts) = fGroups
# plot the heatmaps and expressions
## you can see that in some clusters LTB2 are closer to HC while in others they are closer to LTB1 and HC
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2, bStabalize = T, cex.axis=0.7)

# plot summary heatmaps which show a similar pattern that in some clusters
# LTB2 are like HC and in others they are like ATB
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)

mMarginal = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters

csClust = rownames(mMarginal)
length(csClust)
pdf('Paper/Results/cluster_variance_merged_dataset_after_LTBI_subcluster.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, mMarginal, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(mMarginal[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())

# save the new grouping
n = V(getFinalGraph(oGr))$name
head(n)
mCounts = mCounts[,n]
ltb2_atb_data = list(graph=oGr, groups=fGroups, matrix=mCounts)
save(ltb2_atb_data, file='Paper/Results/ltb2_atb_data.rds')

