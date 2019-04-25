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

## load tb data
lData.tb = f_LoadObject('testData/lData.tb.rds')
names(lData.tb)

## load sepsis data
lData.sepsis = f_LoadObject('testData/lData.sepsis.rds')
names(lData.sepsis)

mCounts.tb = t(lData.tb$data)
mCounts.sepsis = t(lData.sepsis$data)

str(mCounts.sepsis)
str(mCounts.tb)

load('testData/dfPathways.rds')

table(dfPathways$Database)

## separate the pathways and do one database
dfGraph = dfPathways[dfPathways$Database == 'Reactome',c('Gene', 'Pathway')]
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)

# convert gene ids to symbols
library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys=dfGraph$Gene, keytype = 'ENTREZID', columns =  'SYMBOL')
dfGraph$Gene = df$SYMBOL
# select subset of genes from our data set
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts.tb), ]
dfGraph = na.omit(dfGraph)
dim(dfGraph)
length(unique(dfGraph$Gene))
# 1080 genes have annotations

oCGbp.tb = CGraph.bipartite2(dfGraph, ivWeights = c(2, 1, 0))
table(E(getProjectedGraph(oCGbp.tb))$weight)

# 0     1     2 
# 46445 19059  8037 

# some figures
plot.projected.graph(oCGbp.tb, cDropEdges = c(''), bDropOrphans = T)
plot.projected.graph(oCGbp.tb, cDropEdges = c('red'), bDropOrphans = T)
set.seed(123)
plot.projected.graph(oCGbp.tb, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

# extract the igraph object
ig.tb= getProjectedGraph(oCGbp.tb)
table(E(ig.tb)$weight)

# drop all the red edges
ig.tb= delete.edges(ig.tb, which(E(ig.tb)$weight < 2))
vcount(ig.tb)
ecount(ig.tb)
ig.tb= delete.vertices(ig.tb, which(degree(ig.tb) == 0))
vcount(ig.tb)
# 761 genes left

# some plots
set.seed(123)
plot(ig.tb, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig.tb, weights = E(ig.tb)$green))

pdf('temp/omicsGraphs_tb.pdf')
par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
fGroups.tb = lData.tb$grouping
ig.plot.tb = f_igCalculateVertexSizesAndColors(ig.tb, t(mCounts.tb), fGroups.tb, bColor = T, iSize = 10)
set.seed(123)
plot(ig.plot.tb, vertex.label.cex=0.2, layout=layout_with_fr(ig.plot.tb, weights=E(ig.plot.tb)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste(levels(fGroups.tb)[nlevels(fGroups.tb)], 'vs', levels(fGroups.tb)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
dev.off(dev.cur())

## location of the largest clique
set.seed(123)
plot.graph.clique(oCGbp.tb)

## plot subgraph of largest clique
pdf('temp/omicsGraphs_clique_tb.pdf')
ig.plot.tb = induced_subgraph(ig.tb, unlist(largest_cliques(ig.tb)))
ecount(ig.plot.tb)
vcount(ig.plot.tb)
par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
ig.plot.tb = f_igCalculateVertexSizesAndColors(ig.plot.tb, t(mCounts.tb), fGroups.tb, bColor = T, iSize = 20)
set.seed(123)
plot(ig.plot.tb, vertex.label.cex=0.5, layout=layout_with_fr(ig.plot.tb, weights=E(ig.plot.tb)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste(levels(fGroups.tb)[nlevels(fGroups.tb)], 'vs', levels(fGroups.tb)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
dev.off(dev.cur())

# library(lattice)
# df = data.frame(mCounts.tb[,V(ig.plot)$name])
# cSel = c('GBP5', 'GBP6', 'IFI30', 'OAS2', 'FCGR1B', 'RSAD2', 'MT2A')
# df = stack(df[,cSel])
# df$fGroups = fGroups
# df$hiv = lData.tb$adjust
# 
# xyplot( values ~ fGroups | ind, groups=hiv, data=df, type='p', pch=20, scales=list(relation='free'))
# bwplot( values ~ fGroups | ind, data=df[df$hiv == 'HIV-',], type='p', pch=20, scales=list(relation='free'))
# bwplot( values ~ fGroups | ind, data=df[df$hiv == 'HIV+',], type='p', pch=20, scales=list(relation='free'))


###################################################################
### sepsis data set
###################################################################
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
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts.sepsis), ]
dfGraph = na.omit(dfGraph)
dim(dfGraph)
length(unique(dfGraph$Gene))
# 1148 genes have annotations

oCGbp.sepsis = CGraph.bipartite2(dfGraph, ivWeights = c(2, 1, 0))
table(E(getProjectedGraph(oCGbp.sepsis))$weight)
# 0     1     2 
# 55433 23291  9428

# some figures
plot.projected.graph(oCGbp.sepsis, cDropEdges = c(''), bDropOrphans = T)
plot.projected.graph(oCGbp.sepsis, cDropEdges = c('red'), bDropOrphans = T)
set.seed(123)
plot.projected.graph(oCGbp.sepsis, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

# extract the igraph object
ig.sepsis = getProjectedGraph(oCGbp.sepsis)
table(E(ig.sepsis)$weight)

# drop all the red edges
ig.sepsis = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 2))
vcount(ig.sepsis)
ecount(ig.sepsis)
ig.sepsis = delete.vertices(ig.sepsis, which(degree(ig.sepsis) == 0))
vcount(ig.sepsis)
# 836 genes left

# some plots
set.seed(123)
plot(ig.sepsis, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig.sepsis, weights = E(ig.sepsis)$green))

pdf('temp/omicsGraphs_sepsis.pdf')
par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
fGroups.sep = lData.sepsis$grouping
ig.plot.sepsis = f_igCalculateVertexSizesAndColors(ig.sepsis, t(mCounts.sepsis), fGroups.sep, bColor = T, iSize = 5)
set.seed(123)
plot(ig.plot.sepsis, vertex.label.cex=0.2, layout=layout_with_fr(ig.plot.sepsis, weights=E(ig.plot.sepsis)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste(levels(fGroups.sep)[nlevels(fGroups.sep)], 'vs', levels(fGroups.sep)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
dev.off(dev.cur())

## location of the largest clique
set.seed(123)
plot.graph.clique(oCGbp.sepsis)

## plot subgraph of largest clique
pdf('temp/omicsGraphs_clique_sepsis.pdf')
ig.plot.sepsis = induced_subgraph(ig.sepsis, unlist(largest_cliques(ig.sepsis)))
ecount(ig.plot.sepsis)
vcount(ig.plot.sepsis)
par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
ig.plot.sepsis = f_igCalculateVertexSizesAndColors(ig.plot.sepsis, t(mCounts.sepsis), fGroups.sep, bColor = T, iSize = 5)
set.seed(123)
plot(ig.plot.sepsis, vertex.label.cex=0.5, layout=layout_with_fr(ig.plot.sepsis, weights=E(ig.plot.sepsis)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste(levels(fGroups.sep)[nlevels(fGroups.sep)], 'vs', levels(fGroups.sep)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
dev.off(dev.cur())

# df = data.frame(mCounts.sepsis[,V(ig.plot.sepsis)$name])
# df = stack(df)
# df$fGroups = fGroups
# 
# xyplot( values ~ fGroups | ind, data=df, type='smooth', pch=20, scales=list(relation='free'))

########## end sepsis data set

## coreness test tb dataset
# extract the igraph object
ig.tb= getProjectedGraph(oCGbp.tb)
table(E(ig.tb)$weight)

# create sub graphs after edge pruning
ig.tb.y = delete.edges(ig.tb, which(E(ig.tb)$weight < 1))
ig.tb.y = delete.vertices(ig.tb.y, which(degree(ig.tb.y) == 0))

ig.tb.g = delete.edges(ig.tb, which(E(ig.tb)$weight < 2))
ig.tb.g = delete.vertices(ig.tb.g, which(degree(ig.tb.g) == 0))

vcount(ig.tb); vcount(ig.tb.y); vcount(ig.tb.g)

mtb = sapply(c(degree, coreness), function(x){
  return(x(ig.tb))
})

mtb.y = sapply(c(degree, coreness), function(x){
  return(x(ig.tb.y))
})

mtb.g = sapply(c(degree, coreness), function(x){
  return(x(ig.tb.g))
})

plot(mtb, pch=20); cor(mtb)
plot(mtb.y, pch=20); cor(mtb.y)
plot(mtb.g, pch=20); cor(mtb.g)

itb.deg = c(max(mtb[,2]), max(mtb.y[,2]), max(mtb.g[,2]))
itb.tri = c(sum(count_triangles(ig.tb)), sum(count_triangles(ig.tb.y)), sum(count_triangles(ig.tb.g)))

plot(log(itb.tri), log(itb.deg), pch=20)
fit.tb = lm(log(itb.deg) ~ log(itb.tri))
abline(fit.tb)

## coreness test sepsis dataset
# extract the igraph object
ig.sepsis= getProjectedGraph(oCGbp.sepsis)
table(E(ig.sepsis)$weight)

# create sub graphs after edge pruning
ig.sepsis.y = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 1))
ig.sepsis.y = delete.vertices(ig.sepsis.y, which(degree(ig.sepsis.y) == 0))

ig.sepsis.g = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 2))
ig.sepsis.g = delete.vertices(ig.sepsis.g, which(degree(ig.sepsis.g) == 0))

vcount(ig.sepsis); vcount(ig.sepsis.y); vcount(ig.sepsis.g)

msepsis = sapply(c(degree, coreness), function(x){
  return(x(ig.sepsis))
})

msepsis.y = sapply(c(degree, coreness), function(x){
  return(x(ig.sepsis.y))
})

msepsis.g = sapply(c(degree, coreness), function(x){
  return(x(ig.sepsis.g))
})

plot(msepsis, pch=20); cor(msepsis)
plot(msepsis.y, pch=20); cor(msepsis.y)
plot(msepsis.g, pch=20); cor(msepsis.g)

isepsis.deg = c(max(msepsis[,2]), max(msepsis.y[,2]), max(msepsis.g[,2]))
isepsis.tri = c(sum(count_triangles(ig.sepsis)), sum(count_triangles(ig.sepsis.y)), sum(count_triangles(ig.sepsis.g)))

plot(log(isepsis.tri), log(isepsis.deg), pch=20)
fit.sepsis = lm(log(isepsis.deg) ~ log(isepsis.tri))
abline(fit.sepsis)

### combine the data for triangles and degeneracy
iDegeneracy = c(itb.deg, isepsis.deg)
iTriangles = c(itb.tri, isepsis.tri)

plot(log(iTriangles), log(iDegeneracy), pch=20, xlim=c(0, 20), ylim=c(0,6))
fit.joined = lm(log(iDegeneracy) ~ log(iTriangles) - 1)
summary(fit.joined)
abline(fit.joined)
abline(0, 0.3, col=2)

### generate some random graphs
mRandom = matrix(NA, nrow = 20, ncol=2)
colnames(mRandom) = c('tri', 'dege')
for (i in 1:20){
  ig.r1 = erdos.renyi.game(vcount(ig.tb), p.or.m = ecount(ig.tb), type='gnm')
  # count coreness
  mr1 = sapply(c(degree, coreness), function(x){
    return(x(ig.r1))
  })
  mRandom[i,] = c(log(sum(count_triangles(ig.r1))), log(max(mr1[,2])))
}
points(mRandom, pch=20, col=2)

mRandom = matrix(NA, nrow = 20, ncol=2)
colnames(mRandom) = c('tri', 'dege')
for (i in 1:20){
  ig.r1 = erdos.renyi.game(vcount(ig.tb.y), p.or.m = ecount(ig.tb.y), type='gnm')
  # count coreness
  mr1 = sapply(c(degree, coreness), function(x){
    return(x(ig.r1))
  })
  mRandom[i,] = c(log(sum(count_triangles(ig.r1))), log(max(mr1[,2])))
}
points(mRandom, pch=20, col=2)

mRandom = matrix(NA, nrow = 20, ncol=2)
colnames(mRandom) = c('tri', 'dege')
for (i in 1:20){
  ig.r1 = erdos.renyi.game(vcount(ig.tb.g), p.or.m = ecount(ig.tb.g), type='gnm')
  # count coreness
  mr1 = sapply(c(degree, coreness), function(x){
    return(x(ig.r1))
  })
  mRandom[i,] = c(log(sum(count_triangles(ig.r1))), log(max(mr1[,2])))
}
points(mRandom, pch=20, col=2)




########### test 1 - GO Stats

goTest = function(cvSeed, univ = keys(org.Hs.eg.db, 'ENTREZID')){
  library(GOstats)
  ## set up universe background
  dfUniv = AnnotationDbi::select(org.Hs.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
  dfUniv = na.omit(dfUniv)
  univ = unique(dfUniv$ENTREZID)
  
  ## perform test
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
  
  # clean up
  detach("package:GOstats", unload=T)
  detach("package:igraph", unload=T)
  library(igraph)
  return(table(ivPGO.adj < 0.01))
}

## tb graph at cutoffs
ig.tb= getProjectedGraph(oCGbp.tb)
table(E(ig.tb)$weight)

# drop all the red edges
ig.tb= delete.edges(ig.tb, which(E(ig.tb)$weight < 2))
ig.tb= delete.vertices(ig.tb, which(degree(ig.tb) == 0))

goTest(names(unlist(largest_cliques(ig.tb))))






# 
# set.seed(123)
# plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
#      vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5,
#      layout=layout_with_fr(ig, weights = E(ig)$weight))
# 
# set.seed(123)
# plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
#      vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5,
#      layout=layout_with_fr(ig, weights = E(ig)$red))
# 
# set.seed(123)
# plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
#      vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5,
#      layout=layout_with_fr(ig, weights = E(ig)$yellow))


## save the graph object in graphml format to use in cytoscape
write.graph(ig, file= 'temp/tb_graph_weights.graphml', format='graphml')

## import the gene list from cytoscape analysis
cvClusterOne.genes = scan(what=character())
mCounts.sub = mCounts.tb[,cvClusterOne.genes]
dim(mCounts.sub)

ig.sub = induced_subgraph(ig, cvClusterOne.genes)
ecount(ig.sub)
vcount(ig.sub)
fGroups = lData.tb$grouping
fAdjust = lData.tb$adjust

ig.plot = f_igCalculateVertexSizesAndColors(ig.sub, t(mCounts.sub), fGroups, bColor = T, iSize = 20)
set.seed(123)
plot(ig.plot, vertex.label.cex=0.5, layout=layout_with_fr(ig.plot, weights=E(ig.plot)$green),
     vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste(levels(fGroups)[nlevels(fGroups)], 'vs', levels(fGroups)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# community detection
ecount(ig)
#com.r = edge.betweenness.community(ig, weights=E(ig)$yellow)
com.g = cluster_(ig, weights=E(ig)$green)

table(com.g$membership)
mCom = mCompressMatrixByRow(t(mCounts.tb), ig, com.g)

# set.seed(123)
# plot(com.y, ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
#      vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5,
#      layout=layout_with_fr(ig, weights = E(ig)$yellow))

set.seed(123)
plot(com.g, ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$green))


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


