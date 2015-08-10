# File: cluster_test_data.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 24/06/2015
# Desc: generate clusters for data


#library(reactome.db)
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

# convert enterez ids to uniprot
dfMap = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'UNIPROT', 'ENTREZID')
dfMap = na.omit(dfMap)
# load the uniprot2reactome mapping obtained from
# http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
dfReactome = read.csv(file.choose(), header = F, stringsAsFactors=F, sep='\t')
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

# dfGraph = AnnotationDbi::select(reactome.db, colnames(mCounts), 'REACTOMEID', 'ENTREZID')
# dfGraph = na.omit(dfGraph)

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
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')

# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')
par(p.old)
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)
# top  5% of the vertices from each category and largest clique
l = lGetTopVertices(oGr, iQuantile = 0.90)
# top genes based on centrality parameters
dfClique = f_dfGetGeneAnnotation(l$clique)
cvSum = f_csGetGeneSummaryFromGenbank(dfClique$ENTREZID)
cvSum.2 = dfClique$SYMBOL
dfClique$Summary = cvSum[cvSum.2]
write.csv(dfClique, file='Temp/dfClique.csv')

dfDegree = f_dfGetGeneAnnotation(l$degree)
cvSum = f_csGetGeneSummaryFromGenbank(dfDegree$ENTREZID)
cvSum.2 = dfDegree$SYMBOL
dfDegree$Summary = cvSum[cvSum.2]
write.csv(dfDegree, file='Temp/dfDegree.csv')

dfHub = f_dfGetGeneAnnotation(l$hub)
cvSum = f_csGetGeneSummaryFromGenbank(dfHub$ENTREZID)
cvSum.2 = dfHub$SYMBOL
dfHub$Summary = cvSum[cvSum.2]
write.csv(dfHub, file='Temp/dfHub.csv')

dfBetweenness = f_dfGetGeneAnnotation(l$betweenness)
cvSum = f_csGetGeneSummaryFromGenbank(dfBetweenness$ENTREZID)
cvSum.2 = dfBetweenness$SYMBOL
dfBetweenness$Summary = cvSum[cvSum.2]
write.csv(dfBetweenness, file='Temp/dfBetweenness.csv')

dfCloseness = f_dfGetGeneAnnotation(l$closeness)
cvSum = f_csGetGeneSummaryFromGenbank(dfCloseness$ENTREZID)
cvSum.2 = dfCloseness$SYMBOL
dfCloseness$Summary = cvSum[cvSum.2]
write.csv(dfCloseness, file='Temp/dfCloseness.csv')


# plot largest connected graph - clique
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
plot(ig, layout=layout_with_fr)
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
par(p.old)

# plot the graphs of centrality parameter genes 
par(mar=c(1,1,1,1)+0.1)
ig = induced_subgraph(getFinalGraph(oGr), vids = dfDegree$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Degree')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

ig = induced_subgraph(getFinalGraph(oGr), vids = dfCloseness$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Closeness')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

ig = induced_subgraph(getFinalGraph(oGr), vids = dfBetweenness$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Betweenness')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

ig = induced_subgraph(getFinalGraph(oGr), vids = dfHub$ENTREZID)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.5, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey', main='Hub')
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
par(p.old)

## we can look at the problem from the other direction and look at clusters instead of genes
# sample plots
# mean expression of groups in every cluster
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster')
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2, bStabalize = T, cex.axis=0.7)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# expression in all clusters
plot.heatmap.all(oGr, t(mCounts))
# marginal expression level in each cluster
plot.heatmap.marginal(oGr, t(mCounts))
# # plot selected clusters
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280215')
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '168249')
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '3247509')
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '72203')


# get clusters of choice to make subgraphs
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
dfCluster = dfCluster[order(dfCluster$cluster),]

# plot the subgraphs of the significant clusters
l = getSignificantClusters(oGr, mCounts = t(mCounts), fGroups, bStabalize = T)
csClust = rownames(l$clusters)
pdf('Temp/graphs.pdf')
par(mar=c(1,1,1,1)+0.1)
sapply(seq_along(csClust), function(x){
  set.seed(1)
  ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust[x])
  ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T)
  n = f_dfGetGeneAnnotation(V(ig.sub)$name)
  V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
  plot(ig.sub, vertex.label.cex=0.5, layout=layout_with_fr, main=csClust[x])
  ig.sub = getLargestCliqueInCluster(oGr, csClustLabel = csClust[x])
  ig.sub = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts), fGroups, bColor = T)
  n = f_dfGetGeneAnnotation(V(ig.sub)$name)
  V(ig.sub)[n$ENTREZID]$label = n$SYMBOL
  plot(ig.sub, vertex.label.cex=0.8, layout=layout_with_fr, main=csClust[x])
  plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = csClust[x])
})
sym = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
cluster = as.character(dfCluster$cluster)
dfCluster = cbind(sym, cluster)
cvSum = f_csGetGeneSummaryFromGenbank(iID = as.character(dfCluster$ENTREZID))
cvSum.2 = as.character(dfCluster$SYMBOL)
dfCluster$Summary = cvSum[cvSum.2]

write.csv(dfCluster, file='Temp/dfCluster.csv')

# saving graph object to visualize in cytoscape or other graph viewers
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
set.seed(1)
par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color='grey')
set.seed(1)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color='grey')

n = f_dfGetGeneAnnotation(V(ig)$name)
V(ig)[n$ENTREZID]$label = n$SYMBOL
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
dir.create('Test_data', showWarnings = F)
write.csv(dfCluster, 'Test_data/clusters.csv')
write.graph(ig, file = 'Test_data/graph.graphml', format = 'graphml')
