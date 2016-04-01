# File: cluster_ltbi_subgroups.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 31/03/2016
# Desc: generate subclusters in the ltbi dataset


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
load('Objects/ltb_atb_data.rds')
load('Objects/tb_data.rds')
ig.merge = CGraphClust.intersect.union(ltb_atb_data$graph, tb_data$graph)

oGr = ig.merge
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
pdf('Temp/var.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, mMarginal, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(mMarginal[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())

## plot heatmap of this matrix
library(NMF)
m1 = t(mMarginal)
m1 = t(scale(m1))
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)

# plot pca biplot for this
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

## subcluster the ltb group into two subgroups
ivLtb = which(fGroups == 'LTB')
mLtb = mMarginal[,ivLtb]
mLtb = t(mLtb)
# cluster along the row vectors
km.out = kmeans(mLtb, centers = 2, nstart = 20)
km.out
table(km.out$cluster)
fGroups = as.character(fGroups)
fGroups[ivLtb] = paste0(fGroups[ivLtb], km.out$cluster)
fGroups = factor(fGroups, levels = c('HC', 'LTB2', 'LTB1', 'ATB'))

# pca with new groupings
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

# reorder the groups and the matrix
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]
rownames(mCounts) = fGroups
# plot the heatmaps and expressions
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)

# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)

mMarginal = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters

csClust = rownames(mMarginal)
length(csClust)
pdf('Temp/var_sub.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, mMarginal, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(mMarginal[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())


#boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust))

# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

# # plot a cluster of choice as heatmap
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280218')
#### plot a graph of clusters 
#m = getSignificantClusters(oGr, t(mCounts), fGroups, bStabalize = T)
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
data.frame(sort(table(dfCluster$cluster)))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))


# plot the graphs at each contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
#   # plot the largest compoenent only
#   com = components(ig)
#   com.lar = which.max(com$csize)
#   ig = induced_subgraph(ig, vids = V(ig)[which(com$membership == com.lar)])
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}


df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Results/Clusters_ltb_unique.csv')

# # save the graph and data objects
# ltb_atb_data = list(graph=oGr, matrix=mCounts, groups=fGroups)
# save(ltb_atb_data, file='Objects/ltb_atb_data.rds')

# saving graph object to visualize in cytoscape or other graph viewers
lev = levels(fGroups)[-1]
m = mCounts
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  nm = paste('Results/ltb_unique_', lev[i], 'vs', levels(fGroups)[1], '.graphml', sep='')
  write.graph(ig, file = nm, format = 'graphml')
}

#######################################################################################
### selection of plots for various clusters
#######################################################################################


# Various plots for one cluster of choice
csClust = '212436'

lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=60)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

# heatmap of the genes
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)


# if we want to plot variance of one gene at a time
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
rn = rownames(mC)
length(rn)
for (i in seq_along(rn)){
  temp = t(as.matrix(mC[rn[i],]))
  rownames(temp) = rn[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE)
}

# saving graph object to visualize in cytoscape or other graph viewers
csClust = as.character(unique(dfCluster$cluster))
lev = levels(fGroups)[-1]
m = mCounts
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  nm = paste('Results/', lev[i], 'vs', levels(fGroups)[1], '.graphml', sep='')
  write.graph(ig, file = nm, format = 'graphml')
}
# ig = getFinalGraph(oGr)
# n = f_dfGetGeneAnnotation(V(ig)$name)
# V(ig)[n$ENTREZID]$label = n$SYMBOL
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
# write.graph(ig, file = 'Results/graph.graphml', format = 'graphml')
