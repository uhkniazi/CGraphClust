# scratch.R
# Desc: used for testing and rough work


### create a data set to use for bipartite graph
library(org.Hs.eg.db)
library(downloader)
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

lData.train = f_LoadObject('workflow/results/lData.train.rds')
names(lData.train)

oExp = lData.train$data

fGroups = oExp$fSamples
table(fGroups)
table(oExp$fSamples)

## assign annotation names
mCounts = t(exprs(oExp))
dim(mCounts)

## sub select top genes sorted on p-values
cvTopGenes = rownames(lData.train$results)[1:2000]
mCounts = mCounts[,cvTopGenes]

load('workflow/results/dfPathways.rds')

table(dfPathways$Database)

## separate the pathways and do one database at a time

dfGraph = dfPathways[dfPathways$Database == 'Reactome',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)

# select subset of genes
# this will reduce your initial gene list as a lot of genes actually have 
# no annotations, you may use a low level database like GO if you think 
# that you lose too many genes
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))
#mCounts = mCounts[,n]
#print(paste('Total number of genes with terms', length(n)))
#dim(mCounts)

dfGraph = na.omit(dfGraph)

oCGbp.reactome = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.reactome))$weight)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

dfGraph = dfPathways[dfPathways$Database == 'KEGG',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.kegg = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.kegg, bDropOrphans = T)

dfGraph = dfPathways[dfPathways$Database == 'PANTHER',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.panther = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.panther, bDropOrphans = F)

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'BioCarta',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.biocarta = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.biocarta)

ig = CGraph.union(getProjectedGraph(oCGbp.biocarta), oCGbp.kegg@ig.p, oCGbp.panther@ig.p, oCGbp.reactome@ig.p)
table(E(ig)$weight)

# create a correlation graph
oCGcor = CGraph.cor(ig, cor(mCounts))
table(E(getProjectedGraph(oCGcor))$weight)
plot.projected.graph(oCGcor)
ig = CGraph.union(getProjectedGraph(oCGcor), getProjectedGraph(oCGbp.biocarta), oCGbp.kegg@ig.p, oCGbp.panther@ig.p, oCGbp.reactome@ig.p)
table(E(ig)$weight)


ig = delete.edges(ig, which(E(ig)$weight < 4))
vcount(ig)
ig = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$weight), vertex.frame.color=NA)

plot(ig, vertex.label=f_dfGetGeneAnnotation(names(V(ig)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, edge.color='darkgrey')

## add fold changes for infomap method
n = names(V(ig))
V(ig)$weight = abs(lData.train$results[n, 'logFC'])
#E(ig)$weight = E(ig)$weight + abs(min(E(ig)$weight))
# com = cluster_leading_eigen(ig)
# com = cluster_label_prop(ig)
com = cluster_infomap(ig, nb.trials = 100)
#com = cluster_louvain(ig)
table(membership(com))
pdf('temp/graph.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
set.seed(123)
plot(ig, vertex.label=f_dfGetGeneAnnotation(names(V(ig)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig, weights = E(ig)$weight)*10)
set.seed(123)
plot(com, ig, vertex.label=f_dfGetGeneAnnotation(names(V(ig)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.1, layout=layout_with_fr(ig, weights = E(ig)$weight))



m = membership(com)
t = table(membership(com))
t = names(which(t > 11))
m = m[m %in% as.numeric(t)]

ig.s = induced.subgraph(ig, V(ig)[names(m)])

set.seed(123)
plot(ig.s, vertex.label=f_dfGetGeneAnnotation(names(V(ig.s)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig.s, weights = E(ig.s)$weight))

m = membership(com)
t = table(membership(com))
t = names(which.max(t))
m = m[m %in% as.numeric(t)]

ig.s = induced.subgraph(ig, V(ig)[names(m)])

set.seed(123)
plot(ig.s, vertex.label=f_dfGetGeneAnnotation(names(V(ig.s)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig.s, weights = E(ig.s)$weight))


dev.off(dev.cur())
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
###################################################### scratch for distance based graph
dfResults = lData.train$results[unique(dfGraph$Gene),]
plot(dfResults$logFC, scale(-1*log10(dfResults$P.Value)))
mData = cbind(dfResults$logFC, dfResults$P.Value)#-1*log10(dfResults$P.Value))


library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(t(mData), 'FC, PValue')

plot.PCA(oDiag.1, u, legend.pos = 'topright')
p = oDiag.1@lData$PCA

whiten = function(mData){
  ## center the data
  ivMeans = colMeans(mData)
  # centered data
  mData.s = sweep(mData, 2, ivMeans, '-')
  ## calculate covariance matrix
  mCov = cov(mData)
  ## see bishop 2006 chapter 12 page 568 for formula
  # y = 1/sqrt(L) * t(U) * centered data
  ## get the eigen vectors and values
  lEigens = eigen(mCov)
  L = diag(lEigens$values)
  U = lEigens$vectors
  # invert after taking square root
  Z = solve(sqrt(L))
  Z = Z %*% t(U)
  yn = Z %*% t(mData.s)
  rownames(yn) = colnames(mData)
  return(t(yn))
}

plot(whiten(mData))























## make the bipartite graph using the section of script in the main CGraph script
dfGraph = na.omit(dfGraph)
# some error checks
if (ncol(dfGraph) != 2) {
  stop(paste('data frame dfGraph should have 2 columns only', 
             'column 1 for vertex of type 1, and column 2 for vertex of',
             'type 2'))
}
# create bipartite graph
oIGbp = graph.data.frame(dfGraph, directed = F)
# set the vertex type variable to make graph bipartite
f = rep(c(T, F), times = c(length(unique(dfGraph[,1])),length(unique(dfGraph[,2]))))
V(oIGbp)$type = f
# sanity check - is graph bipartite
if (!is.bipartite(oIGbp)) {
  stop(paste('Graph is not bipartite'))
}

## plot this bipartie graph
fType = V(oIGbp)$type
V(oIGbp)[fType]$shape = 'circle'
V(oIGbp)[!fType]$shape = 'square'

#pdf('Temp/Figures/graphs.pdf')
#par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
#plot(oIGbp, layout=layout_as_bipartite, vertex.size=10)


## get alpha and beta parameters for the beta distribution from the data
## see gelman [2] p 583
getalphabeta = function(m, v){
  al.be = (m * (1-m) / v) - 1
  al = al.be * m
  be = al.be * (1-m)
  return(c(alpha=al, beta=be))
}

g.mix = function(theta) mix.prior[1]*exp(m1(theta)) + mix.prior[2]*exp(m2(theta)) + mix.prior[3]*exp(m3(theta))

