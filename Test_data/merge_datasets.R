# File: merge_datasets.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 11/03/2016
# Desc: merge the sepsis and tb datasets

######### load libraries and set global variables
library(org.Hs.eg.db)
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

get.reactome.name = function(csNames){
  i = which(dfReactome.sub$V2 %in% csNames)
  dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
  dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
  rownames(dfCluster.name) = NULL
  return(dfCluster.name)
}

# load the tb and sepsis datasets
load('Objects/tb_data.rds')
load('Objects/sepsis_data.rds')
load('Objects/ltb_atb_data.rds')

ig.tb = getFinalGraph(tb_data$graph)
ig.sep = getFinalGraph(sepsis_data$graph)
ig.lt = getFinalGraph(ltb_atb_data$graph)

# merge the 3 graphs
ig.merge = CGraphClust.union(sepsis_data$graph, ltb_atb_data$graph)
ig.merge = CGraphClust.union(ig.merge, sepsis_data$graph)
plot.final.graph(ig.merge)

par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(ig.merge, t(tb_data$matrix), tb_data$groups, main='TB Significant Clusters', 
                             lwd=1, bStabalize = T, cex.axis=0.7)

par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(ig.merge, t(ltb_atb_data$matrix), ltb_atb_data$groups, main='TB LTBI Significant Clusters', 
                             lwd=1, bStabalize = T, cex.axis=0.7)

par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(ig.merge, t(sepsis_data$matrix), sepsis_data$groups, main='Sepsis Significant Clusters', 
                             lwd=1, bStabalize = T, cex.axis=0.7)


# principal component plots
pr.out = plot.components(ig.merge, t(tb_data$matrix), tb_data$groups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='TB')# principal component plots

pr.out = plot.components(ig.merge, t(ltb_atb_data$matrix), ltb_atb_data$groups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='LTB')# principal component plots

pr.out = plot.components(ig.merge, t(sepsis_data$matrix), sepsis_data$groups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='sepsis')



## create graphs with common genes
iVertID.tb = which(V(ig.tb)$name %in% V(ig.sep)$name)
iVertID.sep = which(V(ig.sep)$name %in% V(ig.tb)$name)

gr.tb = CGraphClust.recalibrate(tb_data$graph, iVertID.tb)
gr.sep = CGraphClust.recalibrate(sepsis_data$graph, iVertID.sep)

set.seed(1)
plot.final.graph(gr.tb)

set.seed(1)
plot.final.graph(gr.sep)

par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(gr.tb, t(tb_data$matrix), tb_data$groups, main='TB Significant Clusters', 
                             lwd=1, bStabalize = T, cex.axis=0.7)

par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(gr.sep, t(sepsis_data$matrix), sepsis_data$groups, main='Sepsis Significant Clusters', 
                             lwd=1, bStabalize = T, cex.axis=0.7)


# principal component plots
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)


# create graphs with different genes
#iVertID.tb = which(!(V(ig.tb)$name %in% V(ig.sep)$name))
ig.merge = CGraphClust.union(ltb_atb_data$graph, tb_data$graph)
ig.lt = getFinalGraph(ig.merge)
ig.sep = getFinalGraph(sepsis_data$graph)
iVertID.lt = which(!(V(ig.lt)$name %in% V(ig.sep)$name))
iVertID.sep = which(!(V(ig.sep)$name %in% V(ig.lt)$name))

#gr.tb = CGraphClust.recalibrate(tb_data$graph, iVertID.tb)
#gr.lt = CGraphClust.recalibrate(ltb_atb_data$graph, iVertID.lt)
gr.lt = CGraphClust.recalibrate(ig.merge, iVertID.lt)
gr.sep = CGraphClust.recalibrate(sepsis_data$graph, iVertID.sep)

set.seed(1)
plot.final.graph(gr.lt)

set.seed(1)
plot.final.graph(gr.sep)
par(mar=c(7, 3, 2, 2)+0.1)
m = plot.significant.expressions(gr.lt, t(ltb_atb_data$matrix), ltb_atb_data$groups, main='LTB difference sepsis', 
                             lwd=1, bStabalize = T, cex.axis=0.7)
m2 = getSignificantClusters(gr.lt, t(ltb_atb_data$matrix), ltb_atb_data$groups)$clusters
get.reactome.name(rownames(m$means))
pr.out = plot.components(gr.lt, t(ltb_atb_data$matrix), ltb_atb_data$groups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='LTB difference sepsis')# principal component plots
par(mfrow=c(2,2))
boxplot.cluster.variance(gr.lt, m2, ltb_atb_data$groups, log=T, iDrawCount = nrow(m2), las=2)
par(p.old)
i = 1
temp = t(as.matrix(m2[rownames(m2)[i],]))
rownames(temp) = rownames(m2)[i]
plot.cluster.variance(gr.lt, temp, ltb_atb_data$groups, log=FALSE); i = i+1

par(p.old)
par(mar=c(7, 3, 2, 2)+0.1)
m = plot.significant.expressions(gr.sep, t(sepsis_data$matrix), sepsis_data$groups, main='Sepsis difference ltb', 
                             lwd=1, bStabalize = T, cex.axis=0.7)
m2 = getSignificantClusters(gr.sep, t(sepsis_data$matrix), sepsis_data$groups)$clusters
get.reactome.name(rownames(m$means))
pr.out = plot.components(gr.sep, t(sepsis_data$matrix), sepsis_data$groups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='sepsis difference ltb')# principal component plots
par(mfrow=c(2,2))
boxplot.cluster.variance(gr.sep, m2, sepsis_data$groups, log=T, iDrawCount = nrow(m2), las=2)
par(p.old)
i = 1
temp = t(as.matrix(m2[rownames(m2)[i],]))
rownames(temp) = rownames(m2)[i]
plot.cluster.variance(gr.sep, temp, sepsis_data$groups, log=FALSE); i = i+1






# intersec the 2 gaphs
ig.merge = CGraphClust.union(ltb_atb_data$graph, tb_data$graph)
ig.sep = getFinalGraph(sepsis_data$graph)
#ig.lt = getFinalGraph(ltb_atb_data$graph)

igi = graph.intersection(getFinalGraph(ig.merge), ig.sep)
d = degree(igi)
igi = delete.vertices(igi, which(d == 0))
par(p.old)
## plotting the intersected graphs
# plot the graph for tb and sepsis
m = tb_data$matrix
fGroups = tb_data$groups
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(tb_data$matrix)
ig = f_igCalculateVertexSizesAndColors(igi, t(m), fGroups, bColor = T, iSize=70)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste('TB', levels(fGroups)[length(levels(fGroups))], 'vs', levels(fGroups)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# plot the graph for sepsis
m = sepsis_data$matrix
fGroups = sepsis_data$groups
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(sepsis_data$matrix)
# relevel 
fGroups = factor(fGroups, levels = c('D1', 'D2', 'D4', 'D5', 'D3'))
ig = f_igCalculateVertexSizesAndColors(igi, t(m), fGroups, bColor = T, iSize=120)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
set.seed(1)
plot(ig, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste('Sepsis', levels(fGroups)[length(levels(fGroups))], 'vs', levels(fGroups)[1]))
legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
##

## select a certain vertex of interest
# get vertex names
n = V(igi)$name
lab = f_dfGetGeneAnnotation(n)
V(igi)$label = as.character(lab$SYMBOL)
plot(igi)

i = which(V(igi)$label == 'SIGLEC5')
i = which(V(igi)$label == 'MT2A')
i = which(V(igi)$label == 'CRCP')
i = which(V(igi)$label == 'MT2A')

cID = V(igi)$name[i]
iNeighSize = 2
## neighbours of selected gene
ig.tb.s = graph.neighborhood(ig.tb, iNeighSize, cID)[[1]]
ig.sep.s = graph.neighborhood(ig.sep, iNeighSize, cID)[[1]]

# ## set the names and colours for plotting
# n = V(ig.tb.s)$name
# lab = f_dfGetGeneAnnotation(n)
# V(ig.tb.s)$label = as.character(lab$SYMBOL)
# V(ig.tb.s)$color = 'grey'
# V(ig.tb.s)[cID]$color = 'red'
# set.seed(1)
# plot(ig.tb.s, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
#      main=paste(cID, 'TB'))
# 
# 
# n = V(ig.sep.s)$name
# lab = f_dfGetGeneAnnotation(n)
# V(ig.sep.s)$label = as.character(lab$SYMBOL)
# V(ig.sep.s)$color = 'grey'
# V(ig.sep.s)[cID]$color = 'red'
# set.seed(1)
# plot(ig.sep.s, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
#      main=paste(cID, 'Sepsis'))

## set names colours and sizes for plotting
# TB
m = tb_data$matrix
fGroups = tb_data$groups
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(tb_data$matrix)
ig = f_igCalculateVertexSizesAndColors(ig.tb.s, t(m), fGroups, bColor = T, iSize=70)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
V(ig)[cID]$color = 'red'

set.seed(1)
plot(ig, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste('TB', levels(fGroups)[length(levels(fGroups))], 'vs', levels(fGroups)[1]))
#legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# sepsis
m = sepsis_data$matrix
fGroups = sepsis_data$groups
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(sepsis_data$matrix)
# relevel 
fGroups = factor(fGroups, levels = c('D1', 'D2', 'D4', 'D5', 'D3'))
ig = f_igCalculateVertexSizesAndColors(ig.sep.s, t(m), fGroups, bColor = T, iSize=120)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)
V(ig)[cID]$color = 'red'

set.seed(1)
plot(ig, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste('Sepsis', levels(fGroups)[length(levels(fGroups))], 'vs', levels(fGroups)[1]))
#legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

#### compare top centrality genes
## tb merge
dfTopGenes.cent = dfGetTopVertices(ig.merge, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(ig.merge)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.tb = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)
rm(dfTopGenes.cent)

## sepsis
dfTopGenes.cent = dfGetTopVertices(sepsis_data$graph, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(sepsis_data$graph)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.sepsis = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)
rm(dfTopGenes.cent)

table(rownames(dfTopGenes.sepsis) %in% rownames(dfTopGenes.tb))
i = which(rownames(dfTopGenes.sepsis) %in% rownames(dfTopGenes.tb))
dfTopGenes.sepsis[i,]

com = compare(getCommunity(tb_data$graph), getCommunity(sepsis_data$graph))


### compare clusters
csClust = '194315'
m = tb_data$matrix
fGroups = tb_data$groups
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(tb_data$matrix)
ig = f_igCalculateVertexSizesAndColors(getClusterSubgraph(tb_data$graph, csClust), t(m), fGroups, bColor = T, iSize=70)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)

set.seed(1)
plot(ig, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste('TB', levels(fGroups)[length(levels(fGroups))], 'vs', levels(fGroups)[1]))
#legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

# sepsis
m = sepsis_data$matrix
fGroups = sepsis_data$groups
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(sepsis_data$matrix)
# relevel 
fGroups = factor(fGroups, levels = c('D1', 'D2', 'D4', 'D5', 'D3'))
ig = f_igCalculateVertexSizesAndColors(getClusterSubgraph(sepsis_data$graph, csClust), t(m), fGroups, bColor = T, iSize=120)
n = V(ig)$name
lab = f_dfGetGeneAnnotation(n)
V(ig)$label = as.character(lab$SYMBOL)

set.seed(1)
plot(ig, vertex.label.cex=0.8, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
     main=paste('Sepsis', levels(fGroups)[length(levels(fGroups))], 'vs', levels(fGroups)[1]))


## heatmaps of the selected cluster
# tb
ig.sub = getClusterSubgraph(tb_data$graph, csClustLabel = csClust)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(tb_data$matrix)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
fGroups = tb_data$groups
# stabalize data
mC = t(apply(mC, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(mC) = fGroups
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)
# sepsis
ig.sub = getClusterSubgraph(sepsis_data$graph, csClustLabel = csClust)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(sepsis_data$matrix)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
fGroups = sepsis_data$groups
# stabalize data
mC = t(apply(mC, 1, function(x) f_ivStabilizeData(x, fGroups)))
colnames(mC) = fGroups
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)

## boxplot
mTb = getClusterMarginal(tb_data$graph, t(tb_data$matrix))
boxplot(mTb[csClust, ] ~ tb_data$groups)

mSep = getClusterMarginal(sepsis_data$graph, t(sepsis_data$matrix))
boxplot(mSep[csClust, ] ~ sepsis_data$groups)
