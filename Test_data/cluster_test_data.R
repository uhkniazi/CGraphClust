# File: cluster_test_data.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 24/06/2015
# Desc: generate clusters for data


library(reactome.db)
library(org.Hs.eg.db)
source('CGraphClust.R')

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
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7)#, iCorCut = 0.7)
#oGr = CGraphClust(dfGraph, (mCor), iCorCut = 0.8)#, iCorCut = 0.7)

# order the count matrix before making heatmaps or plots
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

# sample plots
# mean expression in every cluster
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster')
# only significant clusters
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=2)
# only one cluster
plot.cluster.expressions(oGr, t(mCounts), fGroups, csClustLabel = '3247509', main='cluster')
# plot summary heatmaps
plot.heatmap.all(oGr, t(mCounts))
plot.heatmap.means(oGr, t(mCounts))
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '392499', main='392499 - Metabolism of proteins')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '168249', main='cluster')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280215', main='cluster')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '3247509', main='cluster')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '72203', main='cluster')
# without data stabalization
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = F)
biplot(pr.out)

# with some stabalization
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

# graph properties
n = getLargestCliques(oGr)
n = names(unlist(n))
f_dfGetGeneAnnotation(cvEnterezID = n)
plot.graph.clique(obj = oGr)
# final graph
plot.final.graph(oGr)

# top vertices based on centrality
mCent = mPrintCentralitySummary(oGr)
# top 2% of the vertices from each category and largest clique
l = lGetTopVertices(oGr)
l = unique(unlist(l))
f_dfGetGeneAnnotation(cvEnterezID = l)



# plotting of the igraph object and saving for cytoscape
ig = getFinalGraph(oGr)
p.old = par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$ob_to_ex), vertex.frame.color=NA)
plot(getCommunity(oGr), ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
par(p.old)
df = getClusterMapping(oGr)
colnames(df) = c('gene', 'cluster')
df = df[order(df$cluster),]

write.csv(df, 'Test_data/clusters.csv')
write.graph(ig, file = 'Temp/graph.gml', format = 'graphml')
