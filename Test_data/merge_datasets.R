# File: merge_datasets.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 11/03/2016
# Desc: merge the sepsis and tb datasets

######### load libraries and set global variables
library(org.Hs.eg.db)
source('CGraphClust.R')
# plotting parameters
p.old = par()


# load the tb and sepsis datasets
load('Objects/tb_data.rds')
load('Objects/sepsis_data.rds')

ig.tb = getFinalGraph(tb_data$graph)
ig.sep = getFinalGraph(sepsis_data$graph)

# intersec the 2 gaphs
igi = graph.intersection(ig.tb, ig.sep)
d = degree(igi)
igi = delete.vertices(igi, which(d == 0))

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
## tb
dfTopGenes.cent = dfGetTopVertices(tb_data$graph, iQuantile = 0.80)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(tb_data$graph)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.tb = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)
rm(dfTopGenes.cent)

## sepsis
dfTopGenes.cent = dfGetTopVertices(sepsis_data$graph, iQuantile = 0.80)
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
