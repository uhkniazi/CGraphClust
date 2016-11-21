# File: cluster_dataset_1.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 21/11/2016
# Desc: generate graph and clusters for the tb longitudinal dataset


######### load libraries and set global variables
library(org.Hs.eg.db)
library(downloader)
source('CGraphClust.R')
# plotting parameters
p.old = par()

# load the longitudinal dataset
dfData = read.csv('Paper/Test_data/test_data_GSE19491.csv', header=T, row.names=1)
n = gsub('X(\\d+)', replacement = '\\1', x = colnames(dfData))
colnames(dfData) = n
dim(dfData)

# separate the factor and the count matrix
fGroups = factor(dfData$fSamples, levels = c('12', '2', '0'))
names(fGroups) = rownames(dfData)
mCounts = as.matrix(dfData[,1:(ncol(dfData)-1)])

# convert enterez ids to uniprot as Reactome database file uses UNIPROT ids
dfMap = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'UNIPROT', 'ENTREZID')
dfMap = na.omit(dfMap)

### load the uniprot2reactome mapping obtained from
## any other databsae of choice can be used as long as there are 2 columns
## where Column 1 = Gene and Column 2 = Database ID (label)

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
## where dfMap is the dataframe created earlier with gene ids mapped to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

# select genes that have a reactome term attached
# this will reduce your initial gene list as a lot of genes actually have 
# no annotations, you may use a low level database like GO if you think 
# that you lose too many genes
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
# order the count matrix based on grouping factor
# this will serve useful later for plotting but is not necessary
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

head(fGroups)

# create a correlation matrix to decide cor cutoff
# in a random collection of genes the correlations should follow rougly a 
# normal distribution centered at 0, as this is a subset of differentially expressed
# genes, they should show a bimodal pattern
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation to cluster positively and negatively correlated genes
# the cutoff of 0.7 and absolute is your choice, we would base this cutoff after looking
# at the histogram, if you move this too low then you may end up with a graph with too many
# connections and it may be dominated by noise/random associations, and it will likely crash you machine
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7, bSuppressPlots = F)

# you may not want to suppress the 2 plots, to see the distribution fits
# 1- degree distribution of type 2 vertices, i.e. connections of IDs and Genes
# bipartite graph is created and some cleaning performed, i.e. 
# those type 2 vertices with too many connections (degree) with type 1 vertices are removed.
# this is done by looking at the distribution of the degrees on a log scale and approximating
# a negative binomial and poisson distributions on top. the cutoff is by default 0.95 quantile under
# a negative binomial model, anything over that is removed. The idea is that type vertices that have a
# lot of connections, are not very interesting and they tend to hide the more interesting connections.

# 2- distribution of observed to expected ratio weights for connections between genes
# graph is projected on to one dimension (type 1 vertices) and connections between type 1
# vertices are assigned interestingness or observed to expected probability ratios as weights.
# the weights on a square root scale follow a negative binomial distribution and only the highest
# weights are chosen, as they are the most interesting. So anything over 0.95 quantile under a neg bin model
# are chosen and the other connections are discarded.

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

## if there are small sized communities/clusters then you may want to remove them?
## we normally remove clusters with 5 or less genes
# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
iSizes
# remove communities smaller than 5 members
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

# make a pdf output for publication
dir.create('Paper/Results', showWarnings = F)
pdf('Paper/Results/community_tb_longitudinal.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
set.seed(1)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr,
     vertex.frame.color=NA, edge.color='darkgrey')
set.seed(1)
ig = plot.centrality.graph(oGr)
dev.off(dev.cur())

## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
## centrality in social networks is simpler to understand intuitively 
## e.g. degree = people with lots of followers 
## hub = people with lots of followers connected to people with lots of followers
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
# these plots are bar plots of ordered statistics - sorted centrality parameter
# usually the last 10% (on the right) will show a sudden rise and may be the interesting genes/nodes
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

## top vertices/nodes/genes based on centrality scores
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
## not run
# # loop and get the data from genbank
# n = rep(NA, length=nrow(dfTopGenes.cent))
# names(n) = as.character(dfTopGenes.cent$VertexID)
# for (i in seq_along(n)){
#   n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
#   # this wait time is required to stop sending queries to ncbi server very quickly
#   Sys.sleep(time = 3)
# }
# cvSum.2 = as.character(dfTopGenes.cent$VertexID)
# dfTopGenes.cent$Summary = n[cvSum.2]
####### Section ends

dir.create('Paper/Results', showWarnings = F)
write.csv(dfTopGenes.cent, file='Paper/Results/Top_Centrality_Genes_tb_longitudinal.csv')

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


## in addition to heatmaps the graphs can be plotted
# plot a graph of these top genes
# plot for each contrast i.e. base line vs other level
lev = levels(fGroups)[-1]
m = mCounts
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
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

## instead of looking at individual genes we can look at clusters
## we can look at the problem from the other direction and look at clusters instead of genes
## each cluster is labelled by the most common group ID term shared by the genes in that cluster
# some sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
# the p-values for comparisons are calculated by comparing with the baseline level in the factor fGroups
# and are calculated via simulation from the posterior using a prior parameterized using ebayes approach
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
# principal component plots
# 7 values for each group 0, 2 and 12 are sampled from their posterior mean
# after a hierarchical fit, where we use the full data to calculate a weakly informative prior
# if seed not set then each time the simulated draws may differ slightly / or significantly (if variance very high and noisy 
# set of genes in cluster)
set.seed(123)
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)

# plot summary heatmaps 
# marginal expression level in each cluster
# see methods for details on how this is calculated
# set bStabalize to F, if we do not want to calculate posterior paramters and just plot the data
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = T)

# plot variance of cluster, posterior variances in each cluster may
# highlight clusters that show different variances in each group
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters

csClust = rownames(m)

pdf('Paper/Results/cluster_variance_tb_longitudinal.pdf')
par(mfrow=c(1,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust), las=2)
for (i in seq_along(csClust)){
  temp = t(as.matrix(m[csClust[i],]))
  rownames(temp) = csClust[i]
  plot.cluster.variance(oGr, temp, fGroups, log=FALSE);
}
dev.off(dev.cur())
# examine this figure to see clusters of interest e.g. clusters 109582 and 1280215
# show a marked reduction in variance suggesting that the genes in this cluster may
# be under pressure of drug treatment and are responding similarly as treatment progresses?
# can we find possible biomarkers in this cluster?

# print the labels of the clusters from reactome table
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name


# plot the genes in a cluster of choice as heatmap
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = '1280218')

# write a csv file of gene names and cluster ids
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene

df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Paper/Results/Clusters_tb_longitudinal.csv')

# save the graph and data objects
tb_data = list(graph=oGr, matrix=mCounts, groups=fGroups)
save(tb_data, file='Paper/Results/tb_data.rds')

# saving graph object to visualize in cytoscape or other graph viewers
# not run
# lev = levels(fGroups)[-1]
# m = mCounts
# for(i in 1:length(lev)){
#   ig = getClusterSubgraph(oGr, csClust)
#   fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
#   ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
#   n = V(ig)$name
#   lab = f_dfGetGeneAnnotation(n)
#   V(ig)$label = as.character(lab$SYMBOL)
#   nm = paste('Paper/Results/tb_data', lev[i], 'vs', levels(fGroups)[1], '.graphml', sep='')
#   write.graph(ig, file = nm, format = 'graphml')
# }

