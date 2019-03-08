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

# lData.test = f_LoadObject('testData/lData.test.rds')
# names(lData.test)

mCounts.train = t(lData.train$data)
# mCounts.test = t(lData.test$data)

# str(mCounts.test)
str(mCounts.train)

load('testData/dfPathways.rds')

table(dfPathways$Database)

## separate the pathways and do one database at a time
dfGraph = dfPathways[dfPathways$Database == 'Reactome',c('Gene', 'Pathway')]
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)

# convert gene ids to symbols
library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys=dfGraph$Gene, keytype = 'ENTREZID', columns =  'SYMBOL')
dfGraph$Gene = df$SYMBOL
# select subset of genes from our data set
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts.train), ]
dfGraph = na.omit(dfGraph)
dim(dfGraph)
length(unique(dfGraph$Gene))
# 1075 genes have annotations

oCGbp.reactome = CGraph.bipartite(dfGraph, ivWeights = c(2, 1, 0))
table(E(getProjectedGraph(oCGbp.reactome))$weight)
plot.projected.graph(oCGbp.reactome, cDropEdges = c(''), bDropOrphans = T)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red'), bDropOrphans = T)
set.seed(123)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

ig = getProjectedGraph(oCGbp.reactome)
table(E(ig)$weight)

ig = delete.edges(ig, which(E(ig)$weight < 2))
vcount(ig)
ecount(ig)
ig = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig)

set.seed(123)
plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$green))

set.seed(123)
plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$weight))

set.seed(123)
plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$red))

set.seed(123)
plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$yellow))


## save the graph object in graphml format to use in cytoscape
write.graph(ig, file= 'temp/tb_graph_weights.graphml', format='graphml')



## calculate modularity at different cutoffs
ivMod = rep(NA, times=length(iIndex))
for (i in seq_along(iIndex)){
  ig.p = delete.edges(ig.2, which(E(ig)$weight < iIndex[i]))
  print(table(E(ig.p)$weight))
  com = cluster_louvain(ig.p)
  ivMod[i] = modularity(ig.p, membership(com))
}

plot(ivMod, xaxt='n', main='Modularity')
axis(1, at = 1:length(iIndex), labels = iIndex)

## get the clusters at a particular cutoff
## select 2 large clusters, preferebly similar sizes
## assign go terms to the cluster, and keep go terms that are more common
## fit a model to the data cluster ~ GO terms
## this should be able to identify/predict clusters if the go terms are concentrated in clusters
## repeat this section multiple times manually to select optimal cluster sizes to compare
library(lme4)
iErrorRate = rep(NA, times=length(iIndex))
names(iErrorRate) = iIndex
iAIC = rep(NA, times=length(iIndex))
names(iAIC) = iIndex
length(iIndex)
cutoff = 8
## different cutoffs
ecount(ig.2)
table(E(ig)$weight)
ig.p = delete.edges(ig.2, which(E(ig)$weight < iIndex[cutoff]))
ecount(ig.p)
ig.p = delete.vertices(ig.p, which(degree(ig.p) == 0))
ecount(ig.p)
com = cluster_louvain(ig.p)
dfCom = data.frame(gene=com$names, com=com$membership)
i = sort(table(dfCom$com), decreasing = T)
i
# choose clusters of comparable sizes
i = names(i)[c(1,2)]
dfCom = dfCom[dfCom$com %in% i,]
dfCom$cluster = factor(dfCom$com)
table(dfCom$cluster)
# assign go terms to the genes in the clusters
df = AnnotationDbi::select(org.Hs.eg.db, keys=as.character(dfCom$gene), keytype='SYMBOL', columns='GO')
df = df[df$ONTOLOGY == 'BP', ]
#df = df[df$EVIDENCE != 'TAS', ]
df = na.omit(df)
i = match(df$SYMBOL, as.character(dfCom$gene))
dfCom = dfCom[i,]
identical(as.character(dfCom$gene), df$SYMBOL)
dfCom$GO = factor(df$GO)
# choose more frequent go terms regardless of which cluster they belong to
i = sort(table(dfCom$GO), decreasing = T)
quantile(i, 0:10/10)
i = i[i >= quantile(i, 0.90)]
dfCom = dfCom[dfCom$GO %in% names(i), ]
dfCom = droplevels.data.frame(dfCom)
str(dfCom)
# fit the model
fit.cluster = glmer(cluster ~ 1 + (1|GO), data=dfCom, family='binomial')
summary(fit.cluster)
p = predict(fit.cluster, type='response')
l = levels(dfCom$cluster)
pred = ifelse(p > 0.5, l[2], l[1])
table(pred, actual=dfCom$cluster)
mean(pred != dfCom$cluster)
iErrorRate[cutoff] = mean(pred != dfCom$cluster)
iAIC[cutoff] = AIC(fit.cluster)

plot(iErrorRate, xaxt='n', main='Cluster ~ GO Terms', ylab='Prediction Error', xlab='Edge Weight Cutoff', type='b')
axis(1, at = 1:length(iIndex), labels = iIndex)

plot(iAIC, xaxt='n', main='Cluster Purity: Cluster ID ~ GO Terms', ylab='Akaike Information Criterion', xlab='Edge Weight Cutoff',
     type='b')
axis(1, at = 1:length(iIndex), labels = iIndex)


