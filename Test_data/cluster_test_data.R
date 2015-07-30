# File: cluster_test_data.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 24/06/2015
# Desc: generate clusters for data


library(reactome.db)
library(org.Hs.eg.db)
source('CGraphClust.R')
p.old = par()
# load the test data
dfData = read.csv(file.choose(), header=T, row.names=1)
n = gsub('X(\\d+)', replacement = '\\1', x = colnames(dfData))
colnames(dfData) = n

# separate the factor and the count matrix
fGroups = factor(dfData$fSamples)
mCounts = as.matrix(dfData[,1:(ncol(dfData)-1)])

dfGraph = AnnotationDbi::select(reactome.db, colnames(mCounts), 'REACTOMEID', 'ENTREZID')
dfGraph = na.omit(dfGraph)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]

# create a correlation matrix
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='')

# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7)

# order the count matrix before making heatmaps or plots
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

# plot the main communities and centrality graphs
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
plot(getCommunity(oGr), ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color='grey')
par(p.old)
# look at the graph centrality properties
plot.centrality.graph(oGr)
# plot largest connected graph - clique
plot.graph.clique(oGr)

# get the genes of importance based on graph properties
# largest clique
n = getLargestCliques(oGr)
n = names(unlist(n))
f_dfGetGeneAnnotation(cvEnterezID = n)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)
# top  5% of the vertices from each category and largest clique
l = lGetTopVertices(oGr)
# top genes based on centrality parameters
f_dfGetGeneAnnotation(l$degree)
f_dfGetGeneAnnotation(l$hub)
f_dfGetGeneAnnotation(l$betweenness)

## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster')
# only significant clusters
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2)

# plot summary heatmaps
# expression in all clusters
plot.heatmap.all(oGr, t(mCounts))
# marginal expression level in each cluster
plot.heatmap.marginal(oGr, t(mCounts))
# plot selected clusters
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280215')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '168249')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '3247509')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '72203')

# plot PCA on components of clusters
# with some variance stabalization
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

# get clusters of choice to make subgraphs
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
dfCluster = dfCluster[order(dfCluster$cluster),]

ig = getFinalGraph(oGr)
# look at subgraph properties in clusters
# 1280215
n = as.character(dfCluster[dfCluster$cluster == '1280215', 'gene'])
ig.sub = induced_subgraph(ig, vids = n)
par(mar=c(1,1,1,1)+0.1)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
plot(ig.sub, vertex.label=NA, vertex.size=2, layout=layout_with_fr, main='Cluster 1280215')
v.l = largest_cliques(ig.sub)
ig.sub = induced_subgraph(ig.sub, unlist(v.l))
plot(ig.sub, vertex.label.cex=0.8, vertex.size=20, layout=layout_with_fr, main='Cluster 1280215')

# 168249
n = as.character(dfCluster[dfCluster$cluster == '168249', 'gene'])
ig.sub = induced_subgraph(ig, vids = n)
par(mar=c(1,1,1,1)+0.1)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
plot(ig.sub, vertex.label=NA, vertex.size=2, layout=layout_with_fr, main='Cluster 168249')
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
v.l = largest_cliques(ig.sub)
ig.sub = induced_subgraph(ig.sub, unlist(v.l))
plot(ig.sub, vertex.label.cex=0.8, vertex.size=20, layout=layout_with_fr, main='Cluster 168249')

# 3247509
n = as.character(dfCluster[dfCluster$cluster == '3247509', 'gene'])
ig.sub = induced_subgraph(ig, vids = n)
par(mar=c(1,1,1,1)+0.1)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
plot(ig.sub, vertex.label=NA, vertex.size=2, layout=layout_with_fr, main='Cluster 3247509')
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
v.l = largest_cliques(ig.sub)
ig.sub = induced_subgraph(ig.sub, unlist(v.l))
plot(ig.sub, vertex.label.cex=0.8, vertex.size=20, layout=layout_with_fr, main='Cluster 3247509')

# 72203
n = as.character(dfCluster[dfCluster$cluster == '72203', 'gene'])
ig.sub = induced_subgraph(ig, vids = n)
par(mar=c(1,1,1,1)+0.1)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
plot(ig.sub, vertex.label=NA, vertex.size=2, layout=layout_with_fr, main='Cluster 72203')
V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
v.l = largest_cliques(ig.sub)
ig.sub = induced_subgraph(ig.sub, unlist(v.l))
plot(ig.sub, vertex.label.cex=0.8, vertex.size=20, layout=layout_with_fr, main='Cluster 72203')

# saving graph object to visualize in cytoscape or other graph viewers
ig = getFinalGraph(oGr)
n = f_dfGetGeneAnnotation(V(ig)$name)
V(ig)[n$ENTREZID]$label = n$SYMBOL

dir.create('Test_data', showWarnings = F)
write.csv(dfCluster, 'Test_data/clusters.csv')
write.graph(ig, file = 'Test_data/graph.graphml', format = 'graphml')
