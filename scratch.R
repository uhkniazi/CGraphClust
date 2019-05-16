# scratch.R
# Desc: used for testing and rough work

# https://github.com/uhkniazi/Scratch/blob/0bc2dffc459aa621adebc1076f2f2608ac81dd99/goTermsForJens.R
my_scale = function(x){
  return((x-min(x))/(max(x)-min(x)))
}


############ poisson model
## utility functions
# calculates the gamma prior parameters for a poisson sampling distribution
# see page 5 in notes here: https://www.evernote.com/shard/s288/res/659dc700-ccb9-457e-aea6-aa54bc3abbb9
# and for an example see page 154, chapter on Hierarchical Modeling Bayesian Computation with R.
## DESC: using the poisson sampling model, the data vector is used to count values of alpha (shape), beta (rate)
## parameters for the gamma prior
getalphabeta.poisson = function(lambda){
  m = mean(lambda)
  v = var(lambda)
  alpha = (m^2)/v
  beta = alpha/m
  return(c(alpha=alpha, beta=beta))
}

simGammaPost = function(data, prior){
  alpha = prior['alpha']+sum(data)
  beta = prior['beta']+length(data)
  return(rgamma(1000, shape = alpha, rate = beta))
}

simPostPredictPoisson = function(post, len, nc=20){
  mDraws = matrix(NA, nrow = len, ncol=nc)
  for (i in 1:nc){
    p = sample(post, size = 1)
    mDraws[,i] = rpois(len, p)
  }
  return(mDraws)
}

pr.tb = getalphabeta.poisson(mtb[,1])
mPost.tb = sapply(mtb[,1], function(x) simGammaPost(x, pr.tb))
post.tb = colMeans(mPost.tb)

simPostPredictPoisson2 = function(post, len, nc=20){
  mDraws = matrix(NA, nrow = len, ncol=nc)
  for (i in 1:nc){
    p = sample(1:nrow(post), size = 1)
    mDraws[,i] = rpois(len, post[p,])
  }
  return(mDraws)
}








plot(ig, vertex.label=NA, vertex.size=2, 
     layout=layout_with_fr, vertex.frame.color=NA)

plot(ig, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
     vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5,
     layout=layout_with_fr(ig, weights = E(ig)$green))

plot(ig.s, vertex.label=f_dfGetGeneAnnotation(names(V(ig.s)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig.s, weights = E(ig.s)$weight))



temp = get.edgelist(ig.p)
temp
E(ig.p)
E(ig.p)$weight
which(E(ig.p)$weight > 4)
get.edge.ids(ig.p, temp)
apply(temp, 1, function(x) get.edge.ids(ig.p, x))
ig.p = delete.edges(ig.p, 30)
plot(ig.p)
temp = get.edgelist(ig.p)
dim(temp)
ig.p = delete.edges(ig.p, c(30, 31))



write.graph(ig, file= 'temp/mansoor.graphml', format='graphml')


el = as.matrix(get.edgelist(ig.tb.g))
head(el)
el = el[i,]
el
E(ig.tb.g)[i]
apply(el, 1, function(x) get.edge.ids(ig.tb.g, x))
i
i2 = apply(el, 1, function(x) get.edge.ids(ig.tb, x))
i2
lBet[[1]][i2]
lBet[[3]][i]
xg = lBet[[3]][i]
x = lBet[[1]][i2]
x-xg
log(x)-log(xg)
log(xg)-log(x)
plot(log(xg)-log(x))
i
E(ig.tb.g)[7059]


#################### go
library(GOstats)
# get the universe of genes with go terms
univ = keys(org.Hs.eg.db, 'ENTREZID')
dfUniv = AnnotationDbi::select(org.Hs.eg.db, keys = univ, columns = c('GO'), keytype = 'ENTREZID')
dim(dfUniv)
dfUniv = na.omit(dfUniv)
dim(dfUniv)
univ = unique(dfUniv$ENTREZID)
length(univ)

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

table(ivPGO.adj < 0.01)
p = sort(ivPGO.adj[ivPGO.adj < 0.01], decreasing = F)
columns(GO.db)
temp.3 = select(GO.db, keys=names(p), keytype='GOID', columns=columns(GO.db))


####################

mPrintCentralitySummary = function(ig){
  # calculate 3 measures of centrality i.e. degree, closeness and betweenness
  deg = degree(ig)
  clo = closeness(ig)
  bet = betweenness(ig, directed = F)
  # calculate the page rank and authority_score
  aut = authority_score(ig, scale = F)
  aut = aut$vector
  # print the summary of the data
  print('Degree summary')
  print(summary(deg))
  print('Closeness summary')
  print(summary(clo))
  print('Betweenness summary')
  print(summary(bet))
  print('Hub score summary')
  print(summary(aut))
  # print correlation summary
  m = cbind(degree=deg, closeness=clo, betweenness=bet, hub=aut)
  print('Correlation between scores')
  print(round(cor(m), 3))
  return(m)  
}




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

load('workflow/results/dfPathways_2.rds')

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

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'KEGG',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.kegg = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.kegg, bDropOrphans = T)

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'panther',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.panther = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.panther, bDropOrphans = F)

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'biocarta',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.biocarta = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.biocarta)

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'nci',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.nci = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.nci)

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'PharmGKB',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.pharm = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.pharm)


table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'humancyc',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.humancyc = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.humancyc)

table(dfPathways$Database)
dfGraph = dfPathways[dfPathways$Database == 'SMPDB',c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
length(unique(dfGraph$Gene))

oCGbp.smpdb = CGraph.bipartite(dfGraph)
plot.projected.graph(oCGbp.smpdb)

ig = CGraph.union(getProjectedGraph(oCGbp.biocarta), oCGbp.kegg@ig.p, oCGbp.panther@ig.p, oCGbp.reactome@ig.p, 
                  getProjectedGraph(oCGbp.humancyc), getProjectedGraph(oCGbp.nci), getProjectedGraph(oCGbp.pharm),
                  getProjectedGraph(oCGbp.smpdb))
table(E(ig)$weight)


########### add the go term graphs
dfGraph = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'GO', 'ENTREZID')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'BP',]
dfGraph = dfGraph[,c('ENTREZID', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$ENTREZID))

oCGbp.gobp = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.gobp))$weight)
plot.projected.graph(oCGbp.gobp, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

dfGraph = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'GO', 'ENTREZID')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'CC',]
dfGraph = dfGraph[,c('ENTREZID', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$ENTREZID))

oCGbp.gocc = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.gocc))$weight)
plot.projected.graph(oCGbp.gocc, cDropEdges = c('red', 'yellow'), bDropOrphans = T)


dfGraph = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'GO', 'ENTREZID')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'MF',]
dfGraph = dfGraph[,c('ENTREZID', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$ENTREZID))

oCGbp.gomf = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.gomf))$weight)
plot.projected.graph(oCGbp.gomf, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

ig = CGraph.union(getProjectedGraph(oCGbp.biocarta), oCGbp.kegg@ig.p, oCGbp.panther@ig.p, oCGbp.reactome@ig.p, 
                  getProjectedGraph(oCGbp.humancyc), getProjectedGraph(oCGbp.nci), getProjectedGraph(oCGbp.pharm),
                  getProjectedGraph(oCGbp.smpdb), getProjectedGraph(oCGbp.gobp), getProjectedGraph(oCGbp.gocc),
                  getProjectedGraph(oCGbp.gomf))
table(E(ig)$weight)

# create a correlation graph
oCGcor = CGraph.cor(mCor = cor(mCounts))
table(E(getProjectedGraph(oCGcor))$weight)
table(E(getProjectedGraph(oCGbp.reactome))$weight)
#plot.projected.graph(oCGcor)
# ig = CGraph.union(getProjectedGraph(oCGcor), getProjectedGraph(oCGbp.biocarta), oCGbp.kegg@ig.p, oCGbp.panther@ig.p, oCGbp.reactome@ig.p, 
#                   getProjectedGraph(oCGbp.humancyc), getProjectedGraph(oCGbp.nci), getProjectedGraph(oCGbp.pharm),
#                   getProjectedGraph(oCGbp.smpdb), getProjectedGraph(oCGbp.gobp), getProjectedGraph(oCGbp.gocc),
#                   getProjectedGraph(oCGbp.gomf))

ig = CGraph.union(getProjectedGraph(oCGcor), oCGbp.reactome@ig.p, oCGbp.gobp@ig.p, oCGbp.gomf@ig.p, oCGbp.gocc@ig.p)
table(E(ig)$weight)




ig = delete.edges(ig, which(E(ig)$weight < 3))
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
com = cluster_infomap(ig)#, nb.trials = 100)
com = cluster_walktrap(ig)
m = max((E(ig)$weight))
com = cluster_edge_betweenness(ig, weights = abs((E(ig)$weight)-m)+1)
#com = cluster_louvain(ig)
table(membership(com))
mCent = na.omit(mCompressMatrixByRow((t(mCounts)), ig, com))
library(lattice)
df = data.frame(t(mCent))
df = stack(df)
df$fGroups = fGroups
df$hiv = lData.train$adjust
  
xyplot( values ~ fGroups | ind, groups=hiv, data=df, type='p', pch=20, scales=list(relation='free'))

pdf('temp/graph.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
set.seed(123)
plot(ig, vertex.label=f_dfGetGeneAnnotation(names(V(ig)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig, weights = E(ig)$weight)*10)
set.seed(123)
plot(com, ig, vertex.label=f_dfGetGeneAnnotation(names(V(ig)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.1, layout=layout_with_fr(ig, weights = E(ig)$weight))



m = membership(com)
table(membership(com))
#t = names(which(t >= 10))
m = m[m %in% rownames(mCent)]

ig.s = induced.subgraph(ig, V(ig)[names(m)])

set.seed(123)
plot(ig.s, vertex.label=f_dfGetGeneAnnotation(names(V(ig.s)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig.s, weights = E(ig.s)$weight))

m = membership(com)
t = table(membership(com))
t = names(which.max(t))
m = m[m %in% '24']#  as.numeric(t)]

ig.s = induced.subgraph(ig, V(ig)[names(m)])

set.seed(123)
plot(ig.s, vertex.label=f_dfGetGeneAnnotation(names(V(ig.s)))$SYMBOL, vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr(ig.s, weights = E(ig.s)$weight))


dev.off(dev.cur())
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
###################################################### scratch for distance based graph
dfResults = lData.test$data[V(ig.p)$name,]
plot(dfResults$logFC, scale(-1*log10(dfResults$P.Value)))
mData = dfResults

fGroups = lData.test$grouping
colnames(mData) = as.character(fGroups)
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')
#colnames(mData) = fGroups
oDiag.1 = CDiagnosticPlots(mData, 'selected')

plot.PCA(oDiag.1, fGroups, legend.pos = 'topright', csLabels = '')
plot.dendogram(oDiag.1, fGroups)
plot.mean.summary(oDiag.1, fGroups)
plot.sigma.summary(oDiag.1, fGroups)
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

plot(whiten(t(mCent)), col=c(2,3,4)[as.numeric(fGroups)])
legend('topright', legend = levels(fGroups), fill=c(2,3,4))























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

