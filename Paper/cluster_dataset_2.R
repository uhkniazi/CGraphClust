# File: cluster_test_data_2.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 23/11/2016
# Desc: generate graph and clusters for the hc, ltb and atb dataset


######### load libraries and set global variables
library(org.Hs.eg.db)
library(downloader)
source('CGraphClust.R')
# plotting parameters
p.old = par()

# load the test data
dfData = read.csv('Paper/Test_data/test_data_GSE19491_atb_ltb.csv', header=T, row.names=1)
n = gsub('X(\\d+)', replacement = '\\1', x = colnames(dfData))
colnames(dfData) = n

# separate the factor and the count matrix
fGroups = factor(dfData$fSamples, levels = c('HC', 'LTB', 'ATB'))
names(fGroups) = rownames(dfData)
mCounts = as.matrix(dfData[,1:(ncol(dfData)-1)])

# convert enterez ids to uniprot as Reactome database file uses UNIPROT ids
dfMap = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'UNIPROT', 'ENTREZID')
dfMap = na.omit(dfMap)

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

## map reactome ids to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

# select genes that have a reactome term attached
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
# order the count matrix based on grouping factor
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]
levels(fGroups)

# create a correlation matrix to decide cor cutoff
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 10000, replace = T), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation to cluster positively and negatively correlated genes
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)

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
head(iSizes)
# remove communities smaller than 5 members or choose a size of your liking
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

# make a pdf output for publication
dir.create('Paper/Results', showWarnings = F)
pdf('Paper/Results/community_atb_ltb_hc.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
set.seed(1)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr,
     vertex.frame.color=NA, edge.color='darkgrey')
set.seed(1)
ig = plot.centrality.graph(oGr)
dev.off(dev.cur())

## centrality diagnostics
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
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

write.csv(dfTopGenes.cent, file='Paper/Results/Top_Centrality_Genes_ltb_atb.csv')

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

## if we model the data in a hierarchical manner, and sample values from the posterior means for each group
## this plot may give a clearer/smoother view of the data
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


### Looking at the largest clique can be informative in the graph
# plot the graph with location of the clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

pdf('Paper/Results/largest_clique_atb_ltb_hc.pdf')
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
dev.off(dev.cur())

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
set.seed(1)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
pdf('Paper/Results/pca_clusters_ltb_hc_atb.pdf')
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0, ylim=c(-0.3, 0.35))
dev.off(dev.cur())
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)
# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
csClust = rownames(m)

pdf('Paper/Results/cluster_variance_ltb_tb.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(m[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())

dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene

df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Paper/Results/Clusters_ltb_atb.csv')

# save the graph and data objects
ltb_atb_data = list(graph=oGr, matrix=mCounts, groups=fGroups)
save(ltb_atb_data, file='Paper/Results/ltb_atb_data.rds')

