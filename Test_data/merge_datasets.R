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

########### utility functions
get.reactome.name = function(csNames){
  i = which(dfReactome.sub$V2 %in% csNames)
  dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
  dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
  rownames(dfCluster.name) = NULL
  return(dfCluster.name)
}

## collection of steps to plot some graph diagnostics
f_plot.diagnostics = function(oGr, mCounts, fGroups, ...){
  set.seed(1)
  plot.final.graph(oGr)
  ecount(getFinalGraph(oGr))
  vcount(getFinalGraph(oGr))
  set.seed(1)
  ig = plot.centrality.graph(oGr)
  par(p.old)
  plot.centrality.diagnostics(oGr)
  dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
  rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
  # assign metadata annotation to these genes and clusters
  dfCluster = getClusterMapping(oGr)
  colnames(dfCluster) = c('gene', 'cluster')
  rownames(dfCluster) = dfCluster$gene
  df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
  dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
  dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
  dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)
  library(NMF)
  m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
  m1 = scale(m1)
  m1 = t(m1)
  # threshhold the values
  m1[m1 < -3] = -3
  m1[m1 > 3] = 3
  rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
  # draw the heatmap  color='-RdBu:50'
  aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
           annColors=NA, Colv=NA)
  
  m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
  m1 = apply(m1, 2, f_ivStabilizeData, fGroups)
  rownames(m1) = fGroups
  m1 = scale(m1)
  m1 = t(m1)
  # threshhold the values
  m1[m1 < -3] = -3
  m1[m1 > 3] = 3
  rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
  # draw the heatmap  color='-RdBu:50'
  aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
           annColors=NA, Colv=NA)
  lev = levels(fGroups)[-1]
  m = mCounts
  m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
  rownames(m) = rownames(mCounts)
  par(mar=c(1,1,1,1)+0.1)
  for(i in 1:length(lev)){
    ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
    fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
    ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
    n = V(ig)$name
    lab = f_dfGetGeneAnnotation(n)
    V(ig)$label = as.character(lab$SYMBOL)
    set.seed(1)
    plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
         main=paste(lev[i], 'vs', levels(fGroups)[1]))
    legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
  }
  set.seed(1)
  ig = plot.graph.clique(oGr)
  lev = levels(fGroups)[-1]
  m = mCounts
  #m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
  #rownames(m) = rownames(mCounts)
  par(mar=c(1,1,1,1)+0.1)
  for(i in 1:length(lev)){
    ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
    fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
    ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
    n = V(ig)$name
    lab = f_dfGetGeneAnnotation(n)
    V(ig)$label = as.character(lab$SYMBOL)
    set.seed(1)
    plot(ig, layout=layout_with_fr, main=paste(lev[i], 'vs', levels(fGroups)[1]))
    legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
  }
  par(p.old)
  plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
  # only significant clusters
  par(mar=c(7, 3, 2, 2)+0.1)
  plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
  # principal component plots
  pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
  par(mar=c(4,2,4,2))
  biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
  # plot summary heatmaps
  # marginal expression level in each cluster
  plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
  plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)
  m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
  csClust = rownames(m)
  print(paste('number of significant clusters', length(csClust)))
  print(get.reactome.name(csClust))
  # make graph of clusters
  dfCluster = getClusterMapping(oGr)
  colnames(dfCluster) = c('gene', 'cluster')
  rownames(dfCluster) = dfCluster$gene
  # how many genes in each cluster
  print(data.frame(sort(table(dfCluster$cluster))))
  #csClust = rownames(m$clusters)
  csClust = as.character(unique(dfCluster$cluster))
  # plot the graphs at each contrast
  lev = levels(fGroups)[-1]
  m = mCounts
  #m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
  #rownames(m) = rownames(mCounts)
  par(mar=c(1,1,1,1)+0.1)
  for(i in 1:length(lev)){
    ig = getClusterSubgraph(oGr, csClust)
#     # plot the largest compoenent only
#     com = components(ig)
#     com.lar = which.max(com$csize)
#     ig = induced_subgraph(ig, vids = V(ig)[which(com$membership == com.lar)])
    fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
    ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=50)
    n = V(ig)$name
    lab = f_dfGetGeneAnnotation(n)
    V(ig)$label = as.character(lab$SYMBOL)
    set.seed(1)
    plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
         main=paste(lev[i], 'vs', levels(fGroups)[1]))
    legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
  }
}
f_plot.cluster.diagnostics = function(oGr, mCounts, fGroups, csClust){
  lev = levels(fGroups)[-1]
  m = mCounts
  #m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
  #rownames(m) = rownames(mCounts)
  par(mar=c(1,1,1,1)+0.1)
  for(i in 1:length(lev)){
    ig = getClusterSubgraph(oGr, csClust)
    fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
    ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=60)
    n = V(ig)$name
    lab = f_dfGetGeneAnnotation(n)
    V(ig)$label = as.character(lab$SYMBOL)
    set.seed(1)
    plot(ig, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
         main=paste(lev[i], 'vs', levels(fGroups)[1]))
    legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
  }
  # heatmap of the genes
  ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
  n = f_dfGetGeneAnnotation(V(ig.sub)$name)
  mC = t(mCounts)
  mC = mC[n$ENTREZID,]
  rownames(mC) = n$SYMBOL
  mC = t(scale(t(mC)))
  # threshhold the values
  mC[mC < -3] = -3
  mC[mC > +3] = +3
  # draw the heatmap
  hc = hclust(dist(mC))
  aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
           annColors=NA, Colv=NA)
}

# load the tb and sepsis datasets
load('Objects/tb_data.rds')
load('Objects/sepsis_data.rds')
load('Objects/ltb_atb_data.rds')
load('Objects/sepsis_ns_data.rds')
# ig.tb = getFinalGraph(tb_data$graph)
# ig.sep = getFinalGraph(sepsis_data$graph)
# ig.lt = getFinalGraph(ltb_atb_data$graph)

# merge the sepsis datasets
ig.merge = CGraphClust.intersect.union(sepsis_data$graph, sepsis_ns_data$graph)

f_plot.diagnostics(ig.merge, sepsis_ns_data$matrix, sepsis_ns_data$groups)
f_plot.diagnostics(ig.merge, sepsis_data$matrix, sepsis_data$groups)

# remove genes common between sepsis ns and sepsis
ig.ns = getFinalGraph(sepsis_ns_data$graph)
ig.sur = getFinalGraph(sepsis_data$graph)
# get vertices not common in ltbi graph
iVertID.ns = which(!(V(ig.ns)$name %in% V(ig.sur)$name))
# get the graph
ig.ns.unique = CGraphClust.recalibrate(sepsis_ns_data$graph, iVertID.ns)
f_plot.diagnostics(ig.ns.unique, sepsis_ns_data$matrix, sepsis_ns_data$groups)


# merge sepsis ns and ltbi data
ig.merge = CGraphClust.intersect.union(ltb_atb_data$graph, sepsis_ns_data$graph)
f_plot.diagnostics(ig.merge, sepsis_ns_data$matrix, sepsis_ns_data$groups)
f_plot.diagnostics(ig.merge, ltb_atb_data$matrix, ltb_atb_data$groups)


# merge the 3 graphs
ig.merge = CGraphClust.intersect.union(tb_data$graph, ltb_atb_data$graph)
f_plot.diagnostics(ig.merge, ltb_atb_data$matrix, ltb_atb_data$groups, main='LTBI merged')
f_plot.diagnostics(ig.merge, tb_data$matrix, tb_data$groups, main='long tb and LTBI merged')
# remove genes common between tb longitudinal and ltbi dataset 
ig.ltbi = getFinalGraph(ltb_atb_data$graph)
ig.tblong = getFinalGraph(tb_data$graph)
# get vertices not common in ltbi graph
iVertID.ltbi = which(!(V(ig.ltbi)$name %in% V(ig.tblong)$name))
# get the graph
ig.ltbi.unique = CGraphClust.recalibrate(ltb_atb_data$graph, iVertID.ltbi)
f_plot.diagnostics(ig.ltbi.unique, ltb_atb_data$matrix, ltb_atb_data$groups, main='LTBI unique')
# remove vertices common between ltb and sepsis groups
ig.ltbi = getFinalGraph(ltb_atb_data$graph)
ig.sep = getFinalGraph(sepsis_data$graph)
# get vertices not common in sepsis graph
iVertID.sep = which(!(V(ig.sep)$name %in% V(ig.ltbi)$name))
# get the graph
ig.sep.unique = CGraphClust.recalibrate(sepsis_data$graph, iVertID.sep)
f_plot.diagnostics(ig.sep.unique, sepsis_data$matrix, sepsis_data$groups, main='sepsis unique')



ig.merge.sep = CGraphClust.intersect.union(ig.merge, sepsis_data$graph)
f_plot.diagnostics(ig.merge.sep, ltb_atb_data$matrix, ltb_atb_data$groups, main='LTBI sepsis merged')
f_plot.diagnostics(ig.merge.sep, tb_data$matrix, tb_data$groups, main='TB long sepsis merged')
f_plot.diagnostics(ig.merge.sep, sepsis_data$matrix, sepsis_data$groups, main='Sepsis TB merged')

# plot.final.graph(ig.merge)
# 
# par(mar=c(7, 3, 2, 2)+0.1)
# plot.significant.expressions(ig.merge, t(tb_data$matrix), tb_data$groups, main='TB Significant Clusters', 
#                              lwd=1, bStabalize = T, cex.axis=0.7)
# 
# par(mar=c(7, 3, 2, 2)+0.1)
# plot.significant.expressions(ig.merge, t(ltb_atb_data$matrix), ltb_atb_data$groups, main='TB LTBI Significant Clusters', 
#                              lwd=1, bStabalize = T, cex.axis=0.7)
# 
# par(mar=c(7, 3, 2, 2)+0.1)
# plot.significant.expressions(ig.merge, t(sepsis_data$matrix), sepsis_data$groups, main='Sepsis Significant Clusters', 
#                              lwd=1, bStabalize = T, cex.axis=0.7)
# 
# 
# # principal component plots
# pr.out = plot.components(ig.merge, t(tb_data$matrix), tb_data$groups, bStabalize = T)
# par(mar=c(4,2,4,2))
# biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='TB')# principal component plots
# 
# pr.out = plot.components(ig.merge, t(ltb_atb_data$matrix), ltb_atb_data$groups, bStabalize = T)
# par(mar=c(4,2,4,2))
# biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='LTB')# principal component plots
# 
# pr.out = plot.components(ig.merge, t(sepsis_data$matrix), sepsis_data$groups, bStabalize = T)
# par(mar=c(4,2,4,2))
# biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, main='sepsis')



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
