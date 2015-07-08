# File: cluster_test_data.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 24/06/2015
# Desc: generate clusters for data


library(igraph)
library(reactome.db)
source('CGraphClust.R')

# load the test data
dfData = read.csv(file.choose(), header=T, row.names=1)
n = gsub('X(\\d+)', replacement = '\\1', x = colnames(dfData))
colnames(dfData) = n

# separate the factor and the count matrix
fGroups = dfData$fSamples
mCounts = as.matrix(dfData[,1:(ncol(dfData)-1)])

dfGraph = AnnotationDbi::select(reactome.db, colnames(mCounts), 'REACTOMEID', 'ENTREZID')
dfGraph = na.omit(dfGraph)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]

# create a correlation matrix
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T)

# create the graph cluster object
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)#, iCorCut = 0.7)

# order the count matrix before making heatmaps
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft')
plot.heatmap.means(oGr, t(mCounts))
plot.heatmap.all(oGr, t(mCounts))

ig = getFinalGraph(oGr)
p.old = par(mar=c(1,1,1,1)+0.1)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
plot(getCommunity(oGr), ig, vertex.label=NA, vertex.size=2, layout=layout.fruchterman.reingold, vertex.frame.color=NA)
par(p.old)
df = getClusterMapping(oGr)
colnames(df) = c('gene', 'cluster')
df = df[order(df$cluster),]

write.csv(df, 'Test_data/clusters.csv')
