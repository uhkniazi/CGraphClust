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
hist(sample(mCor, 1000, replace = F), prob=T)

# create the graph cluster object
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.6)#, iCorCut = 0.7)

# order the count matrix before making heatmaps
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft')
plot.significant.expressions(oGr, t(mCounts), fGroups)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)

n = getLargestCliques(oGr)
n = names(unlist(n))
f_dfGetGeneAnnotation(cvEnterezID = n)

plot.graph.clique(obj = oGr)

l = lGetTopVertices(oGr)

f_dfGetGeneAnnotation(cvEnterezID = l)

# try on stabalized data
# mCounts.st = t(mCounts)
# mCounts.st = t(apply(mCounts.st, 1, function(x) f_ivStabilizeData(x, fGroups)))
# colnames(mCounts.st) = fGroups
# plot.significant.expressions(oGr, mCounts.st, fGroups)
# 
# plot.heatmap.means(oGr, t(mCounts))
# plot.heatmap.all(oGr, t(mCounts))

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
