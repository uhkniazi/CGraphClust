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


# #### use the GO database
# dfGraph = AnnotationDbi::select(org.Hs.eg.db, 
#                                 unique(c(colnames(mCounts.test), colnames(mCounts.train)))
#                                 , 'GO', 'SYMBOL')
# dfGraph = dfGraph[dfGraph$ONTOLOGY == 'BP',]
# dfGraph = dfGraph[,c('SYMBOL', 'GO')]
# dfGraph = na.omit(dfGraph)
# str(dfGraph)
# length(unique(dfGraph$SYMBOL))
# 
# oCGbp.gobp = CGraph.bipartite(dfGraph)
# table(E(getProjectedGraph(oCGbp.gobp))$weight)
# plot.projected.graph(oCGbp.gobp, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

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
                  #getProjectedGraph(oCGbp.gobp),
                  getProjectedGraph(oCGbp.disgenet))
table(E(ig)$weight)

ig = delete.edges(ig, which(E(ig)$weight < 1))
vcount(ig)
ecount(ig)
ig.p = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig.p)
plot(ig.p, vertex.label=NA, vertex.size=2, layout=layout_with_fr, vertex.frame.color=NA)

# create 2 correlation graphs
oCGcor.train = CGraph.cor(ig.template = ig.p, mCor = cor(mCounts.train), ivWeights = c(1, 0, -1)) 
table(E(getProjectedGraph(oCGcor.train))$weight)


oCGcor.test = CGraph.cor(ig.template = ig.p, mCor = cor(mCounts.test), ivWeights = c(1, 0, -1))
table(E(getProjectedGraph(oCGcor.test))$weight)

ig = CGraph.union(getProjectedGraph(oCGcor.train),
                  getProjectedGraph(oCGcor.test),
                  getProjectedGraph(oCGbp.reactome),
                  #getProjectedGraph(oCGbp.gobp),
                  getProjectedGraph(oCGbp.disgenet))
table(E(ig)$weight)

## save the graph object in graphml format to use in cytoscape
write.graph(ig, file= 'temp/tb_graph_weights.graphml', format='graphml')


# genes at different weights in kegg pathway
# belonging to TB
iIndex = sort(unique(E(ig)$weight))
iIndex = iIndex[-c(length(iIndex))]
mRes = matrix(NA, length(iIndex), ncol = 2)
rownames(mRes) = iIndex
for (i in seq_along(iIndex)){
  ig.p = delete.edges(ig, which(E(ig)$weight < iIndex[i]))
  ig.p = delete.vertices(ig.p, which(degree(ig.p) == 0))
  cvSeed = V(ig.p)$name
  dfKegg = dfPathways[dfPathways$Database == 'KEGG',c('Gene', 'Pathway')]
  dfKegg$Gene = as.character(dfKegg$Gene)
  dfKegg$Pathway = as.character(dfKegg$Pathway)
  df = AnnotationDbi::select(org.Hs.eg.db, keys=dfKegg$Gene, keytype = 'ENTREZID', columns =  'SYMBOL')
  dfKegg$Gene = df$SYMBOL
  dfKegg.sub = dfKegg[dfKegg$Gene %in% cvSeed,]
  mRes[i,] = as.numeric(table(dfKegg.sub$Pathway %in% 'hsa:05152'))
}
#mRes[7,2] = 0
x = (mRes[,2]+0.01) / mRes[,1]
plot(x, xaxt='n', main='TB in Kegg Annotation', ylab='Proportion Genes Detected', xlab='Edge Weight Cutoff', type='b')
axis(1, at = 1:length(iIndex), labels = iIndex)
## number of significant go terms 
library(GOstats)
# get the universe of genes with go terms
univ = keys(org.Hs.eg.db, 'ENTREZID')
dfUniv = AnnotationDbi::select(org.Hs.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
dim(dfUniv)
dfUniv = na.omit(dfUniv)
dim(dfUniv)
univ = unique(dfUniv$ENTREZID)
length(univ)
iIndex = sort(unique(E(ig)$weight))
iIndex = iIndex[-c(length(iIndex))]
mRes = matrix(NA, length(iIndex), ncol = 2)
for (i in seq_along(iIndex)){
  ig.p = delete.edges(ig, which(E(ig)$weight < iIndex[i]))
  ig.p = delete.vertices(ig.p, which(igraph::degree(ig.p) == 0))
  cvSeed = V(ig.p)$name
  df = AnnotationDbi::select(org.Hs.eg.db, keys=cvSeed, keytype = 'SYMBOL', columns =  'ENTREZID')
  
  ## make hypergeometric test object for each type, CC, BP and MF
  params = new('GOHyperGParams', geneIds=unique(df$ENTREZID),
               annotation='org.Hs.eg.db',
               universeGeneIds=univ,
               ontology='BP',
               pvalueCutoff= 0.01,
               conditional=FALSE,
               testDirection='over')
  
  oGOStat = hyperGTest(params) 
  # get pvalues
  ivPGO = pvalues(oGOStat)
  # fdr
  ivPGO.adj = p.adjust(ivPGO, 'BH')
  
  mRes[i,] = table(ivPGO.adj < 0.01)
}
rownames(mRes) = iIndex
plot(mRes[,2], xaxt='n', main='Enrichment of GO Terms', xlab='Edge Weight Cutoff', ylab='Count', type='b')
axis(1, at = 1:length(iIndex), labels = iIndex)

rs = rowSums(mRes)
plot(mRes[,2]/rs, xaxt='n', main='Proportion of GO Enrichment', xlab='Edge Weight Cutoff', ylab='Proportion', type='b')
axis(1, at = 1:length(iIndex), labels = iIndex)

## random list of genes at cutoff 3
ig.p = delete.edges(ig, which(E(ig)$weight < 3))
table(E(ig.p)$weight)
ig.p = delete.vertices(ig.p, which(igraph::degree(ig.p) == 0))
ecount(ig.p)
vcount(ig.p)
cvGenes = V(ig)$name

f1 = function(){
  cvSeed = sample(cvGenes, 136)
  df = AnnotationDbi::select(org.Hs.eg.db, keys=cvSeed, keytype = 'SYMBOL', columns =  'ENTREZID')
  ## make hypergeometric test object for each type, CC, BP and MF
  params = new('GOHyperGParams', geneIds=unique(df$ENTREZID),
               annotation='org.Hs.eg.db',
               universeGeneIds=univ,
               ontology='BP',
               pvalueCutoff= 0.01,
               conditional=FALSE,
               testDirection='over')
  
  oGOStat = hyperGTest(params) 
  # get pvalues
  ivPGO = pvalues(oGOStat)
  # fdr
  ivPGO.adj = p.adjust(ivPGO, 'BH')
  ret = sum((ivPGO.adj < 0.01))
  return((ret+0.5)/(length(ivPGO.adj)+1))
}

go.sim = replicate(30, f1())

## random list of genes at cutoff -1
ig.p = delete.edges(ig, which(E(ig)$weight < -1))
table(E(ig.p)$weight)
ig.p = delete.vertices(ig.p, which(igraph::degree(ig.p) == 0))
ecount(ig.p)
vcount(ig.p)
cvGenes = V(ig)$name

f1 = function(){
  cvSeed = sample(cvGenes, 1171)
  df = AnnotationDbi::select(org.Hs.eg.db, keys=cvSeed, keytype = 'SYMBOL', columns =  'ENTREZID')
  ## make hypergeometric test object for each type, CC, BP and MF
  params = new('GOHyperGParams', geneIds=unique(df$ENTREZID),
               annotation='org.Hs.eg.db',
               universeGeneIds=univ,
               ontology='BP',
               pvalueCutoff= 0.01,
               conditional=FALSE,
               testDirection='over')
  
  oGOStat = hyperGTest(params) 
  # get pvalues
  ivPGO = pvalues(oGOStat)
  # fdr
  ivPGO.adj = p.adjust(ivPGO, 'BH')
  ret = sum((ivPGO.adj < 0.01))
  return((ret+0.5)/(length(ivPGO.adj)+1))
}

go.sim.minus1 = replicate(30, f1())

library(lattice)
df = data.frame(go.sim, go.sim.minus1)
colnames(df) = c('Weight 3', 'Weight -1')
df = stack(df)
histogram( ~ values | ind, data=df)
hist(go.sim)

# house keeping
detach("package:GOstats", unload=T)
detach("package:igraph", unload=T)
library(igraph)
ig.2 = ig
E(ig.2)$weight = E(ig.2)$weight + abs(min(E(ig.2)$weight))
table(E(ig.2)$weight)
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

cutoff = 6
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


