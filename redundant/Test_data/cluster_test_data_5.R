# File: cluster_test_data_5.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 23/03/2016
# Desc: generate clusters for integrated data from ltbi and longitudinal data


######### load libraries and set global variables
library(org.Hs.eg.db)
library(downloader)
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

load('Objects/ltb2_atb_data.rds')
mCounts = ltb2_atb_data$matrix
fGroups = ltb2_atb_data$groups
oGr = ltb2_atb_data$graph


## general graph structure
## we would like to see how does the graph look like, are the clusters connected or in subgraphs
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))

## community structure
## overview of how the commuinties look like
# plot the main communities in 2 different ways
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')

# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
# remove communities smaller than 5 members or choose a size of your liking
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

# # make a pdf output for publication
# dir.create('Results', showWarnings = F)
# pdf('Results/Graph_structure.pdf')
# par(mar=c(1,1,1,1)+0.1, family='Helvetica')
# ig = getFinalGraph(oGr)
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
# set.seed(1)
# plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
#      vertex.frame.color=NA, edge.color='darkgrey')
# set.seed(1)
# ig = plot.centrality.graph(oGr)
# dev.off(dev.cur())

## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)


## top vertices based on centrality scores
## get a table of top vertices 
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

####### NOTE: This section of code is very slow, use only if you need data from genbank
# loop and get the data from genbank
n = rep(NA, length=nrow(dfTopGenes.cent))
names(n) = as.character(dfTopGenes.cent$VertexID)
for (i in seq_along(n)){
  n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
  # this wait time is required to stop sending queries to ncbi server very quickly
  Sys.sleep(time = 3)
}
cvSum.2 = as.character(dfTopGenes.cent$VertexID)
dfTopGenes.cent$Summary = n[cvSum.2]
####### Section ends

dir.create('Results', showWarnings = F)
write.csv(dfTopGenes.cent, file='Results/Top_Centrality_Genes_ltb2_atb.csv')

## if we want to look at the expression profiles of the top genes
# plot a heatmap of these top genes
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

## in addition to heatmaps the graphs can be plotted
# plot a graph of these top genes
# plot for each contrast i.e. base line vs other level
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


### Looking at the largest clique can be informative in the graph
# plot the graph with location of the clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique at each grouping contrast
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


## instead of looking at individual genes we can look at clusters
## we can look at the problem from the other direction and look at clusters instead of genes
# some sample plots
# mean expression of groups in every cluster
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
# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters

csClust = rownames(m)
length(csClust)
pdf('Results/cluster_variance_ltb2_tb.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(m[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())
#boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust))

# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name

# # plot a cluster of choice as heatmap
# plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280218')
#### plot a graph of clusters 
#m = getSignificantClusters(oGr, t(mCounts), fGroups, bStabalize = T)
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
data.frame(sort(table(dfCluster$cluster)))
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
  #   # plot the largest compoenent only
  #   com = components(ig)
  #   com.lar = which.max(com$csize)
  #   ig = induced_subgraph(ig, vids = V(ig)[which(com$membership == com.lar)])
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


#dfCluster = dfCluster[dfCluster$cluster %in% csClust,]
df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Results/Clusters_ltb2_atb.csv')

# saving graph object to visualize in cytoscape or other graph viewers
lev = levels(fGroups)[-1]
m = mCounts
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  nm = paste('Results/ltb_atb_data', lev[i], 'vs', levels(fGroups)[1], '.graphml', sep='')
  write.graph(ig, file = nm, format = 'graphml')
}

#######################################################################################
### selection of plots for various clusters
#######################################################################################


# Various plots for one cluster of choice
csClust = '109582'

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


# if we want to plot variance of one gene at a time
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
rn = rownames(mC)
length(rn)

for (i in seq_along(rn)){
  temp = t(as.matrix(mC[rn[i],]))
  rownames(temp) = rn[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}


# saving graph object to visualize in cytoscape or other graph viewers
csClust = as.character(unique(dfCluster$cluster))
lev = levels(fGroups)[-1]
m = mCounts
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  nm = paste('Results/', lev[i], 'vs', levels(fGroups)[1], '.graphml', sep='')
  write.graph(ig, file = nm, format = 'graphml')
}
# ig = getFinalGraph(oGr)
# n = f_dfGetGeneAnnotation(V(ig)$name)
# V(ig)[n$ENTREZID]$label = n$SYMBOL
# ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T)
# write.graph(ig, file = 'Results/graph.graphml', format = 'graphml')
