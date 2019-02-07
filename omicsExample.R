# File: omicsExample.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 17/11/2018
# Desc: analysis of omics example data sets

source('CGraphClust.R')
p.old = par()

# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

lData.train = f_LoadObject('testData/lData.train.rds')
names(lData.train)

lData.test = f_LoadObject('testData/lData.test.rds')
names(lData.test)

mCounts.train = t(lData.train$data)
mCounts.test = t(lData.test$data)

str(mCounts.test)
str(mCounts.train)

load('testData/dfPathways.rds')

table(dfPathways$Database)

## separate the pathways and do one database at a time
dfGraph = dfPathways[dfPathways$Database == 'Reactome',c('Gene', 'Pathway')]
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)

library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys=dfGraph$Gene, keytype = 'ENTREZID', columns =  'SYMBOL')
dfGraph$Gene = df$SYMBOL
# select subset of genes from our 2 matrices
dfGraph = dfGraph[dfGraph$Gene %in% unique(c(colnames(mCounts.test), colnames(mCounts.train))), ]
dfGraph = na.omit(dfGraph)
dim(dfGraph)
length(unique(dfGraph$Gene))
length(unique(c(colnames(mCounts.test), colnames(mCounts.train))))

oCGbp.reactome = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.reactome))$weight)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red', 'yellow'), bDropOrphans = T)


#### use the GO database
dfGraph = AnnotationDbi::select(org.Hs.eg.db, 
                                unique(c(colnames(mCounts.test), colnames(mCounts.train)))
                                , 'GO', 'SYMBOL')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'BP',]
dfGraph = dfGraph[,c('SYMBOL', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$SYMBOL))

oCGbp.gobp = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.gobp))$weight)
plot.projected.graph(oCGbp.gobp, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

### use disgenet database
load('testData/dfDisgenet.rds')
dfGraph = dfDisgenet[dfDisgenet$score > 0.1,c('SYMBOL', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$SYMBOL)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[,-1]
dfGraph = dfGraph[,c(2,1)]

dfGraph = dfGraph[dfGraph$Gene %in% unique(c(colnames(mCounts.test), colnames(mCounts.train))), ]
length(unique(dfGraph$Gene))

oCGbp.disgenet = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.disgenet))$weight)
plot.projected.graph(oCGbp.disgenet, bDropOrphans = T, cDropEdges = c('red', 'yellow'))

# create a template graph for making correlation graphs
ig = CGraph.union(getProjectedGraph(oCGbp.reactome),
                  getProjectedGraph(oCGbp.gobp),
                  getProjectedGraph(oCGbp.disgenet))
table(E(ig)$weight)

ig = delete.edges(ig, which(E(ig)$weight < 1))
vcount(ig)
ecount(ig)
ig.p = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig.p)
plot(ig.p, vertex.label=NA, vertex.size=2, layout=layout_with_fr, vertex.frame.color=NA)

# create 2 correlation graphs
oCGcor.train = CGraph.cor(ig.template = ig.p, mCor = cor(mCounts.train), ivWeights = c(1, 0, 0)) 
table(E(getProjectedGraph(oCGcor.train))$weight)


oCGcor.test = CGraph.cor(ig.template = ig.p, mCor = cor(mCounts.test), ivWeights = c(1, 0, 0))
table(E(getProjectedGraph(oCGcor.test))$weight)

ig = CGraph.union(getProjectedGraph(oCGcor.train),
                  getProjectedGraph(oCGcor.test),
                  getProjectedGraph(oCGbp.reactome),
                  getProjectedGraph(oCGbp.gobp),
                  getProjectedGraph(oCGbp.disgenet))
table(E(ig)$weight)

## save the graph object in graphml format to use in cytoscape
write.graph(ig, file= 'temp/tb_graph_weights.graphml', format='graphml')


# genes at different weights
ig.p = delete.edges(ig, which(E(ig)$weight < 4))
vcount(ig.p)
ecount(ig.p)
ig.p = delete.vertices(ig.p, which(degree(ig.p) == 0))
vcount(ig.p)
#write.csv(V(ig.p)$name, file='temp/names.csv')
#write.csv(AnnotationDbi::select(org.Hs.eg.db, keys=V(ig.p)$name, keytype = 'SYMBOL', columns =  'ENTREZID'), file='temp/names.csv')
cvSeed = V(ig.p)$name

dfKegg = dfPathways[dfPathways$Database == 'KEGG',c('Gene', 'Pathway')]
dfKegg$Gene = as.character(dfKegg$Gene)
dfKegg$Pathway = as.character(dfKegg$Pathway)
str(dfKegg)

df = AnnotationDbi::select(org.Hs.eg.db, keys=dfKegg$Gene, keytype = 'ENTREZID', columns =  'SYMBOL')
dfKegg$Gene = df$SYMBOL

dfKegg.sub = dfKegg[dfKegg$Gene %in% cvSeed,]
dim(dfKegg)
dim(dfKegg.sub)
sort(table(dfKegg.sub$Pathway))
table(dfKegg.sub$Pathway %in% 'hsa:05152')





## calculate modularity at different cutoffs
ivMod = rep(NA, times=5)
for (i in 0:4){
  ig.p = delete.edges(ig, which(E(ig)$weight < i))
  com = cluster_louvain(ig.p)
  ivMod [i+1] = modularity(ig.p, membership(com))
}

## import the list of genes after analysis in cytoscape
dfGenes = read.csv(file.choose(), header=T, stringsAsFactors = F)
ig.p = induced.subgraph(ig, V(ig)[dfGenes$name])
ig.p = delete.edges(ig.p, which(E(ig.p)$weight < 1))

## grouping factor to calculate fold changes
fGroups = lData.test$grouping
levels(fGroups)
ig.draw = f_igCalculateVertexSizesAndColors(ig.p, mCounts = lData.test$data, fGroups = fGroups, bColor = T, iSize = 40)
plot(ig.draw, vertex.label=names(V(ig.draw)), vertex.label.cex=1,  
     edge.color='darkgrey', edge.width=1, layout=layout_with_fr)

# 
# 
# ig = delete.edges(ig, which(E(ig)$weight < 2))
# vcount(ig)
# ecount(ig)
# ig.p = delete.vertices(ig, which(degree(ig) == 0))
# vcount(ig.p)
# plot(ig.p, vertex.label=NA, vertex.size=2, layout=layout_with_fr, vertex.frame.color=NA)
# 
# pdf('temp/graph2.pdf')
# par(mar=c(1,1,1,1)+0.1, family='Helvetica')
# set.seed(123)
# plot(ig.p, vertex.label=names(V(ig.p)), vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr)
# dev.off(dev.cur())
# 
# ## calculate clusters and modularity
# # m = max((E(ig.p)$weight))
# # com = cluster_edge_betweenness(ig.p, weights = abs((E(ig.p)$weight)-m)+1)
# com = cluster_louvain(ig.p)
# table(membership(com))
# modularity(ig.p, membership(com))
# 
# ## recalculate after removing smaller communities
# cl = clusters(ig.p)
# cl$csize
# i = which(cl$csize < 5)
# v = which(cl$membership %in% i)
# # delete the components that are small
# ig.p = delete.vertices(ig.p, v = v)
# 
# ## calculate clusters and modularity
# # m = max((E(ig.p)$weight))
# # com = cluster_edge_betweenness(ig.p, weights = abs((E(ig.p)$weight)-m)+1)
# com = cluster_spinglass(ig.p)
# table(membership(com))
# modularity(ig.p, membership(com))
# 
# hc = as.hclust(com)
# plot(hc)
# 
# pdf('temp/graph.pdf')
# par(mar=c(1,1,1,1)+0.1, family='Helvetica')
# set.seed(123)
# plot(com, ig.p, vertex.label=names(V(ig.p)), vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr)
# dev.off(dev.cur())
# 
# 
# df = data.frame(n=com$names, c=com$membership)
# df
# membership(com)
# df = df[order(df$c),]
# write.csv(df, file='temp/df.csv')
# 
