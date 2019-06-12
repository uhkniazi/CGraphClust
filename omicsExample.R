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
# 
# # some figures
# plot.projected.graph(oCGbp.sepsis, cDropEdges = c(''), bDropOrphans = T)
# plot.projected.graph(oCGbp.sepsis, cDropEdges = c('red'), bDropOrphans = T)
# set.seed(123)
# plot.projected.graph(oCGbp.sepsis, cDropEdges = c('red', 'yellow'), bDropOrphans = T)
# 
# # extract the igraph object
# ig.sepsis = getProjectedGraph(oCGbp.sepsis)
# table(E(ig.sepsis)$weight)
# 
# # drop all the red edges
# ig.sepsis = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 2))
# vcount(ig.sepsis)
# ecount(ig.sepsis)
# ig.sepsis = delete.vertices(ig.sepsis, which(degree(ig.sepsis) == 0))
# vcount(ig.sepsis)
# # 836 genes left
# 
# # some plots
# set.seed(123)
# plot(ig.sepsis, vertex.label=NA, vertex.label.cex=0.1, vertex.size=2, 
#      vertex.frame.color=NA, 
#      edge.color='darkgrey', edge.width=0.5,
#      layout=layout_with_fr(ig.sepsis, weights = E(ig.sepsis)$green))
# 
# pdf('temp/omicsGraphs_sepsis.pdf')
# par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
# fGroups.sep = lData.sepsis$grouping
# ig.plot.sepsis = f_igCalculateVertexSizesAndColors(ig.sepsis, t(mCounts.sepsis), fGroups.sep, bColor = T, iSize = 5)
# set.seed(123)
# plot(ig.plot.sepsis, vertex.label.cex=0.2, layout=layout_with_fr(ig.plot.sepsis, weights=E(ig.plot.sepsis)$green),
#      vertex.frame.color='darkgrey', edge.color='lightgrey', 
#      main=paste(levels(fGroups.sep)[nlevels(fGroups.sep)], 'vs', levels(fGroups.sep)[1]))
# legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
# dev.off(dev.cur())
# 
# ## location of the largest clique
# set.seed(123)
# plot.graph.clique(oCGbp.sepsis)
# 
# ## plot subgraph of largest clique
# pdf('temp/omicsGraphs_clique_sepsis.pdf')
# ig.plot.sepsis = induced_subgraph(ig.sepsis, unlist(largest_cliques(ig.sepsis)))
# ecount(ig.plot.sepsis)
# vcount(ig.plot.sepsis)
# par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
# ig.plot.sepsis = f_igCalculateVertexSizesAndColors(ig.plot.sepsis, t(mCounts.sepsis), fGroups.sep, bColor = T, iSize = 5)
# set.seed(123)
# plot(ig.plot.sepsis, vertex.label.cex=0.5, layout=layout_with_fr(ig.plot.sepsis, weights=E(ig.plot.sepsis)$green),
#      vertex.frame.color='darkgrey', edge.color='lightgrey', 
#      main=paste(levels(fGroups.sep)[nlevels(fGroups.sep)], 'vs', levels(fGroups.sep)[1]))
# legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
# dev.off(dev.cur())

# df = data.frame(mCounts.sepsis[,V(ig.plot.sepsis)$name])
# df = stack(df)
# df$fGroups = fGroups
# 
# xyplot( values ~ fGroups | ind, data=df, type='smooth', pch=20, scales=list(relation='free'))

########## end sepsis data set

################ add some tests and comparisons for results
################ coreness test tb dataset
# extract the igraph object
ig.tb= getProjectedGraph(oCGbp.tb)
table(E(ig.tb)$weight)

# create sub graphs after edge pruning
ig.tb.y = delete.edges(ig.tb, which(E(ig.tb)$weight < 1))
ig.tb.y = delete.vertices(ig.tb.y, which(degree(ig.tb.y) == 0))

ig.tb.g = delete.edges(ig.tb, which(E(ig.tb)$weight < 2))
ig.tb.g = delete.vertices(ig.tb.g, which(degree(ig.tb.g) == 0))

vcount(ig.tb); vcount(ig.tb.y); vcount(ig.tb.g)

## calculate degree and coreness
mtb = sapply(c(degree, coreness), function(x){
  return(x(ig.tb))
})

mtb.y = sapply(c(degree, coreness), function(x){
  return(x(ig.tb.y))
})

mtb.g = sapply(c(degree, coreness), function(x){
  return(x(ig.tb.g))
})

# visualise results
plot(mtb, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Red'); cor(mtb)
plot(mtb.y, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Yellow'); cor(mtb.y)
plot(mtb.g, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Green'); cor(mtb.g)

################### atypical patterns in coreness vs degree plots
## largest cliques location
# red graph
i = names(unlist(largest_cliques(ig.tb)))
i2 = which(rownames(mtb) %in% i)
plot(mtb, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Red')
points(mtb[i2,], pch=20, col=2)
# yellow graph
i = names(unlist(largest_cliques(ig.tb.y)))
i2 = which(rownames(mtb.y) %in% i)
plot(mtb.y, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Yellow')
points(mtb.y[i2,], pch=20, col=2)
# green graph
i = names(unlist(largest_cliques(ig.tb.g)))
i2 = which(rownames(mtb.g) %in% i)
plot(mtb.g, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Green', cex=0.7, col='grey')
points(mtb.g[i2,], pch=20, col='green')
# 3 graphs together
i = names(unlist(largest_cliques(ig.tb)))
i2 = which(rownames(mtb) %in% i)
plot(mtb, pch=20, xlab='Degree', ylab='Coreness', main='TB data - Red', cex=0.7, col='grey')
points(mtb[i2,], pch=20, col=2)
i = names(unlist(largest_cliques(ig.tb.y)))
i2 = which(rownames(mtb) %in% i)
points(mtb[i2,], pch=20, col='yellow')
i = names(unlist(largest_cliques(ig.tb.g)))
i2 = which(rownames(mtb) %in% i)
points(mtb[i2,], pch=20, col='green')

## extract the green graph at top and selected core cutoffs to make figures
par(mar=c(1,1,1,1)+0.1, mfrow=c(2,2))
sort(unique(mtb.g[,2]))
# repeat with appropriate choices of cutoffs
iCut = 44
i = names(which(mtb.g[,2] == iCut))
ig.plot.tb = induced_subgraph(ig.tb.g, i)
ecount(ig.plot.tb)
vcount(ig.plot.tb)
# plot the graph
ig.plot.tb = f_igCalculateVertexSizesAndColors(ig.plot.tb, t(mCounts.tb), fGroups.tb, bColor = T, iSize = 20)
set.seed(123)
plot(ig.plot.tb, vertex.label=NA, layout=layout_with_fr(ig.plot.tb, weights=E(ig.plot.tb)$green),
     vertex.frame.color='grey', edge.color='grey',edge.width=0.3, 
     main=paste(iCut, levels(fGroups.tb)[nlevels(fGroups.tb)], 'vs', levels(fGroups.tb)[1]))
#legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))

##################### Poisson Models for degree distribution from Gelman book and paper example
## degree distribution and poisson models
iDeg = degree(ig.tb.g)

## create a random graph for comparison and checking of the model 
ig.ran = erdos.renyi.game(vcount(ig.tb.g), p.or.m = ecount(ig.tb.g), type='gnm')
iDeg.ran = degree(ig.ran)

## stan does not like this package so unload annotation packages
detach("package:org.Hs.eg.db", unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='testData/poissonRegressionNoPooling.stan')
stanDso.2 = rstan::stan_model(file='testData/poissonRegressionPartialPooling_1.stan')

## y is observed data, e is exposure or offset, 1 in this case as log 1 = 0 in the model
# dfData = data.frame(y=iDeg.ran, e=1)
# # complete pooling model - ER model
# m = model.matrix(y ~ 1, data=dfData)
# 
# lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
#                  y=dfData$y)
# 
# fit.stan = sampling(stanDso, data=lStanData, iter=5000, chains=2, pars=c('betas', 'mu'),
#                     cores=2)
# print(fit.stan, c('betas'), digits=3)
# # extract fitted values and generate some posterior predictive samples
# mFitted = extract(fit.stan)$mu
# dim(mFitted)
# mFitted = mFitted[sample(1:nrow(mFitted), 20, replace = F),]
# mSim.fit.ran = apply(mFitted, 1, function(x){
#   return(rpois(length(x), exp(x)))
# })
# # plots to compare with observed data
# hist(dfData$y, prob=T)
# temp = apply(mSim.fit.ran, 2, function(x) lines(density(x)))

### apply erdos renyi model to observed data in green graph
dfData = data.frame(y=iDeg, e=1, gene=factor(names(iDeg)))
m = model.matrix(y ~ 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('betas', 'mu'),
                    cores=4)
print(fit.stan, c('betas'), digits=3)

mFitted.er.green = extract(fit.stan)$mu
dim(mFitted.er.green)
mFitted = mFitted.er.green[sample(1:nrow(mFitted.er.green), 100, replace = F),]
mSim.er.green = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})
rm(mFitted)
# hist(dfData$y, prob=T)
# temp = apply(mSim.fit.ran, 2, function(x) lines(density(x)))

## apply the second model, partial pooling model, by adding gene information
m = model.matrix(y ~ gene - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso.2, data=lStanData, iter=1000, chains=4, pars=c('populationMean', 'sigmaRan', 'betas', 'mu'),
                    cores=4)
print(fit.stan, c('populationMean', 'sigmaRan'), digits=3)

pairs(fit.stan, pars = c("sigmaRan", "populationMean", "lp__"))
# some diagnostics for stan
traceplot(fit.stan, c('sigmaRan'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('populationMean'), ncol=1, inc_warmup=F)

## quick comparison with lme4 results
library(lme4)
fit.lme = glmer(y ~ 1 + (1 | gene), offset=log(dfData$e), data=dfData, family=poisson(link = "log"))
## compare lme4 and stan
summary(fit.lme)

## extract the coefficients and fitted values of interest
mCoef = extract(fit.stan)$betas
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$gene))
colnames(mCoef) =  as.character(d$mods)
mCoef.green = mCoef; rm(mCoef)

s = extract(fit.stan)$sigmaRan
mCoef.green = sweep(mCoef.green, 1, s, FUN = '*')

p = extract(fit.stan)$populationMean
mCoef.green = sweep(mCoef.green, 1, p, '+')

mFitted.pp.green = extract(fit.stan)$mu
dim(mFitted.pp.green)
mFitted = mFitted.pp.green[sample(1:nrow(mFitted.pp.green), 100, replace = F),]
mSim.pp.green = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})
rm(mFitted)

########## apply the 2 models again to the red graph
ig.tb= getProjectedGraph(oCGbp.tb)
table(E(ig.tb)$weight)
iDeg = degree(ig.tb)
summary(iDeg)
dfData = data.frame(y=iDeg, e=1, gene=factor(names(iDeg)))
m = model.matrix(y ~ 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('betas', 'mu'),
                    cores=4)
print(fit.stan, c('betas'), digits=3)

mFitted.er.red = extract(fit.stan)$mu
dim(mFitted.er.red)
mFitted = mFitted.er.red[sample(1:nrow(mFitted.er.red), 100, replace = F),]
mSim.er.red = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})
rm(mFitted)
# hist(dfData$y, prob=T)
# temp = apply(mSim.pp.red, 2, function(x) lines(density(x)))

## apply the second model, partial pooling model, by adding gene information
m = model.matrix(y ~ gene - 1, data=dfData)

lStanData = list(Ntotal=nrow(dfData), Ncol=ncol(m), X=m, offset=log(dfData$e),
                 y=dfData$y)

fit.stan = sampling(stanDso.2, data=lStanData, iter=1000, chains=4, pars=c('populationMean', 'sigmaRan', 'betas', 'mu'),
                    cores=4)
print(fit.stan, c('populationMean', 'sigmaRan'), digits=3)

pairs(fit.stan, pars = c("sigmaRan", "populationMean", "lp__"))
# some diagnostics for stan
traceplot(fit.stan, c('sigmaRan'), ncol=1, inc_warmup=F)
traceplot(fit.stan, c('populationMean'), ncol=1, inc_warmup=F)

## quick comparison with lme4 results
library(lme4)
fit.lme = glmer(y ~ 1 + (1 | gene), offset=log(dfData$e), data=dfData, family=poisson(link = "log"))
## compare lme4 and stan
summary(fit.lme)

## extract the coefficients and fitted values of interest
mCoef = extract(fit.stan)$betas
dim(mCoef)
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$gene))
colnames(mCoef) =  as.character(d$mods)
mCoef.red = mCoef; rm(mCoef)

s = extract(fit.stan)$sigmaRan
mCoef.red = sweep(mCoef.red, 1, s, FUN = '*')

p = extract(fit.stan)$populationMean
mCoef.red = sweep(mCoef.red, 1, p, '+')

mFitted.pp.red = extract(fit.stan)$mu
dim(mFitted.pp.red)
mFitted = mFitted.pp.red[sample(1:nrow(mFitted.pp.red), 100, replace = F),]
mSim.pp.red = apply(mFitted, 1, function(x){
  return(rpois(length(x), exp(x)))
})
rm(mFitted)

#### figures with data and posterior predictive simulations
par(mfrow=c(3,2))
hist(degree(ig.tb), xlim=c(0, 700), prob=T, xlab='', main='Observed Degree Distribution - Red')
hist(degree(ig.tb.g), xlim=c(0, 90), xaxt='n', prob=T, xlab='', main='Observed Degree Distribution - Green')
axis(1, c(0, 20, 40, 60, 80, 90), labels = c(0, 20, 40, 60, 80, 90))
hist(mSim.er.red[,1], xlim=c(0, 700), prob=T, xlab='', main='Erdos-Renyi Model - Red')
apply(mSim.er.red, 2, function(x) lines(density(x), lwd=0.8, col='grey'))
hist(mSim.er.green[,1], xlim=c(0, 90), xaxt='n', prob=T, xlab='', main='Erdos-Renyi Model - Green')
axis(1, c(0, 20, 40, 60, 80, 90), labels = c(0, 20, 40, 60, 80, 90))
apply(mSim.er.green, 2, function(x) lines(density(x), lwd=0.8, col='grey'))
hist(mSim.pp.red[,1], xlim=c(0, 700), prob=T, xlab='', main='Partial Pooling Model - Red')
apply(mSim.pp.red, 2, function(x) lines(density(x), lwd=0.8, col='grey'))
hist(mSim.pp.green[,1], xlim=c(0, 90), xaxt='n', prob=T, xlab='', main='Partial Pooling Model - Green')
axis(1, c(0, 20, 40, 60, 80, 90), labels = c(0, 20, 40, 60, 80, 90))
apply(mSim.pp.green, 2, function(x) lines(density(x), lwd=0.8, col='grey'))

##### comparisons of the extracted coefficients
g = colMeans(mCoef.green)
r = colMeans(mCoef.red)
# 
# i = match(names(g), names(r))
# r = r[i]
# identical(names(g), names(r))
# plot(g, r, pch=20)

# #### clustering of networks to assign groupings to genes
# com.r = cluster_louvain(ig.tb, weights=NULL)
# table(com.r$membership)
# head(com.r$membership)
# rc = com.r$membership
# names(rc) = com.r$names
# i = match(names(r), names(rc))
# rc = rc[i]
# identical(names(r), names(rc))
# head(r); head(rc)
# rc = factor(rc)
# i = which(table(rc) > 10)
# table(rc)[i]
# df = data.frame(r, rc)
# df = df[df$rc %in% i, ]
# library(lattice)
# bwplot(r ~ rc, data=df)
# 
# ## green
# com.g = cluster_louvain(ig.tb.g, weights=NULL)
# gc = com.g$membership
# names(gc) = com.g$names
# i = match(names(g), names(gc))
# gc = gc[i]
# identical(names(g), names(gc))
# gc = factor(gc)
# i = which(table(gc) > 10)
# table(gc)[i]
# df = data.frame(g, gc)
# df = df[df$gc %in% i, ]
# bwplot(g ~ gc, data=df)
# summary(lm(g ~ gc, data=df))

############# checking for cluster purity
## get the clusters at red or green graphs
## select 2 large clusters, preferebly similar sizes
## assign go terms to the cluster, and keep go terms that are more common in frequency
## fit a model to the data cluster_id ~ GO terms
## this should be able to identify/predict clusters if the go terms are concentrated in clusters
## repeat this section multiple times manually to select optimal cluster sizes to compare
library(lme4)

getScore = function(ids){
  ig.p = delete.edges(ig.tb, ids)
  ig.p = delete.vertices(ig.p, which(degree(ig.p) == 0))
  com = cluster_louvain(ig.p, weight=NULL)
  ## map the cluster id to the gene name
  dfCom = data.frame(gene=com$names, com=com$membership)
  # get the frequency/size of the clusters
  i = sort(table(dfCom$com), decreasing = T)
  # choose clusters of comparable sizes
  i = names(i)[c(1,2)]
  # subset the data to the 2 largest clusters
  dfCom = dfCom[dfCom$com %in% i,]
  dfCom$cluster = factor(dfCom$com)
  table(dfCom$cluster)
  str(dfCom)
  # assign go terms to the genes in the clusters, which means go terms are assigned to clusters
  df = AnnotationDbi::select(org.Hs.eg.db, keys=as.character(dfCom$gene), keytype='SYMBOL', columns='GO')
  df = df[df$ONTOLOGY == 'BP', ]
  #df = df[df$EVIDENCE != 'TAS', ]
  df = na.omit(df)
  i = match(df$SYMBOL, as.character(dfCom$gene))
  dfCom = dfCom[i,]
  identical(as.character(dfCom$gene), df$SYMBOL)
  dfCom$GO = factor(df$GO)
  ## reduce the number of go terms i.e. drop rare terms or rare factor levels
  # choose more frequent go terms regardless of which cluster they belong to
  i = sort(table(dfCom$GO), decreasing = T)
  quantile(i, 0:10/10)
  # choose the most frequent
  i = i[i >= quantile(i, 0.90)]
  dfCom = dfCom[dfCom$GO %in% names(i), ]
  dfCom = droplevels.data.frame(dfCom)
  str(dfCom)
  # fit the model and calculate AIC
  fit.cluster = glmer(cluster ~ 1 + (1|GO), data=dfCom, family='binomial')
  summary(fit.cluster)
  return(c(aic=AIC(fit.cluster), edges=ecount(ig.p), vertices=vcount(ig.p)))
}

table(E(ig.tb)$weight)
iAIC = rep(NA, times=2)
names(iAIC) = c('Red', 'Green')
iAIC['Red'] = getScore(which(E(ig.tb)$weight > 2))['aic']
iAIC['Green'] = getScore(which(E(ig.tb)$weight < 2))['aic']

# 
# # weights in the graph
# iIndex = c(0, 2)
# iErrorRate = rep(NA, times=length(iIndex))
# names(iErrorRate) = c('Red', 'Green')
# iAIC = rep(NA, times=length(iIndex))
# names(iAIC) = names(iErrorRate)
# length(iIndex)
# ## weight cutoff index
# cutoff = 2
# ## different cutoffs
# ecount(ig.tb)
# table(E(ig.tb)$weight)
# ig.p = delete.edges(ig.tb, which(E(ig.tb)$weight < iIndex[cutoff]))
# ecount(ig.p)
# ig.p = delete.vertices(ig.p, which(degree(ig.p) == 0))
# vcount(ig.p)
# com = cluster_louvain(ig.p, weight=NULL)
# ## map the cluster id to the gene name
# dfCom = data.frame(gene=com$names, com=com$membership)
# # get the frequency/size of the clusters
# i = sort(table(dfCom$com), decreasing = T)
# i
# # choose clusters of comparable sizes
# i = names(i)[c(1,2)]
# # subset the data to the 2 largest clusters
# dfCom = dfCom[dfCom$com %in% i,]
# dfCom$cluster = factor(dfCom$com)
# table(dfCom$cluster)
# str(dfCom)
# # assign go terms to the genes in the clusters, which means go terms are assigned to clusters
# df = AnnotationDbi::select(org.Hs.eg.db, keys=as.character(dfCom$gene), keytype='SYMBOL', columns='GO')
# df = df[df$ONTOLOGY == 'BP', ]
# #df = df[df$EVIDENCE != 'TAS', ]
# df = na.omit(df)
# i = match(df$SYMBOL, as.character(dfCom$gene))
# dfCom = dfCom[i,]
# identical(as.character(dfCom$gene), df$SYMBOL)
# dfCom$GO = factor(df$GO)
# ## reduce the number of go terms i.e. drop rare terms or rare factor levels
# # choose more frequent go terms regardless of which cluster they belong to
# i = sort(table(dfCom$GO), decreasing = T)
# quantile(i, 0:10/10)
# # choose the most frequent
# i = i[i >= quantile(i, 0.90)]
# dfCom = dfCom[dfCom$GO %in% names(i), ]
# dfCom = droplevels.data.frame(dfCom)
# str(dfCom)
# # fit the model and calculate AIC
# fit.cluster = glmer(cluster ~ 1 + (1|GO), data=dfCom, family='binomial')
# summary(fit.cluster)
# p = predict(fit.cluster, type='response')
# l = levels(dfCom$cluster)
# pred = ifelse(p > 0.5, l[2], l[1])
# table(pred, actual=dfCom$cluster)
# mean(pred != dfCom$cluster)
# iErrorRate[cutoff] = mean(pred != dfCom$cluster)
# iAIC[cutoff] = AIC(fit.cluster)
# 
# df = data.frame(y=c(iErrorRate, iAIC), z=factor(c(1,1,2,2), labels=c('Prediction Error Rate', 'AIC')),
#                 x=factor(c(1,2,1,2), labels=c('Red', 'Green')))
# barchart(y ~ x | z, data=df, scales=list(relation='free'), ylab='Model Predictive Score', main='GO Term Purity')

#### repeat with random sizes of edges to delete
mEdges = get.edgelist(ig.tb)
table(E(ig.tb)$weight)
dim(mEdges)

iRange = round(seq(5000, 70000, length.out = 20),0)
# create 20 subsets of these random graphs
lSub = lapply(iRange, function(x) sample(1:nrow(mEdges), x))
mScores = sapply(lSub, getScore)

plot(c(mScores['edges',], 8000, 75000), c(mScores['aic',], iAIC), type='n',
     xlab='No. of Edges', ylab='AIC', cex.axis=0.8)

lines(mScores['edges',], mScores['aic',])
points(ecount(ig.tb), iAIC['Red'], pch=20, cex=2, col='red')
points(ecount(ig.tb.g), iAIC['Green'], pch=20, cex=2, col='green')

## maximal cliques 
# iCliques.green = sapply(max_cliques(ig.tb.g, 3, length(largest_cliques(ig.tb.g)[[1]])), length)
# iCliques.red = sapply(max_cliques(ig.tb, 3, length(largest_cliques(ig.tb)[[1]])), length)

getCliques = function(ids){
  ig.p = delete.edges(ig.tb, ids)
  ig.p = delete.vertices(ig.p, which(degree(ig.p) == 0))
  return(sapply(max_cliques(ig.p, 3, length(largest_cliques(ig.p)[[1]])), length))
}
iCliques.green = getCliques(which(E(ig.tb)$weight < 2))
iCliques.yellow = getCliques(which(E(ig.tb)$weight < 1))
iCliques.red = getCliques(which(E(ig.tb)$weight > 2))
iCliques.ran = getCliques(lSub[[20]])
iCliques.ran.2 = getCliques(lSub[[19]])
iCliques.ran.3 = getCliques(lSub[[18]])
iCliques.ran.4 = getCliques(lSub[[11]])
iCliques.ran.5 = getCliques(lSub[[15]])
iCliques.ran.6 = getCliques(lSub[[13]])
par(mfrow=c(3,3))
plot(density(iCliques.green), xlab='Cliques', main=paste0('Green E-count ', ecount(ig.tb.g)), ylab='')
plot(density(iCliques.yellow), xlab='Cliques', main=paste0('Yellow E-count ', ecount(ig.tb.y)), ylab='')
plot(density(iCliques.red), xlab='Cliques', main=paste0('Red E-count ', ecount(ig.tb)), ylab='')
plot(density(iCliques.ran), xlab='Cliques', main=paste0('Random E-count ', ecount(ig.tb) - length(lSub[[20]])), ylab='')
plot(density(iCliques.ran.2), xlab='Cliques', main=paste0('Random E-count ', ecount(ig.tb) - length(lSub[[19]])), ylab='')
plot(density(iCliques.ran.3), xlab='Cliques', main=paste0('Random E-count ', ecount(ig.tb) - length(lSub[[18]])), ylab='')
plot(density(iCliques.ran.5), xlab='Cliques', main=paste0('Random E-count ', ecount(ig.tb) - length(lSub[[15]])), ylab='')
plot(density(iCliques.ran.6), xlab='Cliques', main=paste0('Random E-count ', ecount(ig.tb) - length(lSub[[13]])), ylab='')
plot(density(iCliques.ran.4), xlab='Cliques', main=paste0('Random E-count ', ecount(ig.tb) - length(lSub[[11]])), ylab='')

## create a random graph for comparison and checking of the model 
ig.ran = erdos.renyi.game(vcount(ig.tb), p.or.m = ecount(ig.tb), type='gnm')
ig.ran.y = erdos.renyi.game(vcount(ig.tb.y), p.or.m = ecount(ig.tb.y), type='gnm')
ig.ran.g = erdos.renyi.game(vcount(ig.tb.g), p.or.m = ecount(ig.tb.g), type='gnm')
ER.1 = sapply(max_cliques(ig.ran, 3, length(largest_cliques(ig.ran)[[1]])), length)
ER.2 = sapply(max_cliques(ig.ran.y, 3, length(largest_cliques(ig.ran.y)[[1]])), length)
ER.3 = sapply(max_cliques(ig.ran.g, 3, length(largest_cliques(ig.ran.g)[[1]])), length)



plot(density(iCliques.green[iCliques.green <= 20]), xlab='Cliques', main=paste0('Green E-count ', ecount(ig.tb.g)), ylab='')
plot(density(iCliques.yellow[iCliques.yellow <= 20]), xlab='Cliques', main=paste0('Yellow E-count ', ecount(ig.tb.y)), ylab='')
plot(density(iCliques.red[iCliques.red <= 20]), xlab='Cliques', main=paste0('Red E-count ', ecount(ig.tb)), ylab='')
plot(density(ER.1), xlab='Cliques', main=paste0('ER E-count ', ecount(ig.ran)), ylab='')
plot(density(ER.2), xlab='Cliques', main=paste0('ER E-count ', ecount(ig.ran.y)), ylab='')
plot(density(ER.3), xlab='Cliques', main=paste0('ER E-count ', ecount(ig.ran.g)), ylab='')

#########################################################
####### Graph union of tb and sepsis
#########################################################
ig.tb = getProjectedGraph(oCGbp.tb)
ig.sepsis = getProjectedGraph(oCGbp.sepsis)
ig.tb.g = delete.edges(ig.tb, which(E(ig.tb)$weight < 2))
ig.tb.g = delete.vertices(ig.tb.g, which(degree(ig.tb.g) == 0))
ig.sepsis.g = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 2))
ig.sepsis.g = delete.vertices(ig.sepsis.g, which(degree(ig.sepsis.g) == 0))
table(E(ig.sepsis.g)$weight); table(E(ig.tb.g)$weight)
vcount(ig.sepsis.g); vcount(ig.tb.g)

ig.un = CGraph.union(ig.tb.g, ig.sepsis.g)
table(E(ig.un)$weight)
## find genes only present in TB and sepsis
cTb = V(ig.tb.g)$name
cSepsis = V(ig.sepsis.g)$name
table(cTb %in% cSepsis)
V(ig.un)$color = 'lightgrey'
V(ig.un)[cTb]$color = 'lightgreen' # tb color
V(ig.un)[cSepsis]$color = 'moccasin'
## find common genes
cCommon = V(ig.un)$name
f = (cCommon %in% cTb) & (cCommon %in% cSepsis)
table(f)
cCommon = cCommon[f]
head(cCommon)
V(ig.un)[cCommon]$color = 'lightskyblue1'
# keep only the largest connected component
c = components(ig.un)
cComponent = names(c$membership[c$membership != 1])
length(cComponent)
ig.un = delete.vertices(ig.un, cComponent)
## find largest cliques in tb and sepsis
c = names(unlist(largest_cliques(ig.tb.g)[[1]]))
V(ig.un)[c]$color = 'pink'
c = names(unlist(largest_cliques(ig.sepsis.g)[[1]]))
V(ig.un)[c]$color = 'yellow'


set.seed(123)
par(mar=c(1,1,1,1)+0.1)
plot(ig.un, vertex.label.cex=0.2, layout=layout_with_fr(ig.un, weights=rep(1, times=ecount(ig.un))),
     vertex.frame.color='darkgrey', edge.color='lightgrey', vertex.size=5)
legend('topright', legend = c('TB', 'Sepsis', 'Common'), fill = c('lightgreen', 'moccasin', 'lightskyblue1'))
## save the graph object in graphml format to use in cytoscape
write.graph(ig.un, file= 'temp/union.graphml', format='graphml')





# 
# plot(density(iCliques.green), xlab='Cliques', main=paste0('Green E-count ', ecount(ig.tb.g)), ylab='')
# plot(density(iCliques.yellow), xlab='Cliques', main=paste0('Yellow E-count ', ecount(ig.tb.y)), ylab='')
# plot(density(iCliques.red), xlab='Cliques', main=paste0('Red E-count ', ecount(ig.tb)), ylab='')


# iRange = rep(ecount(ig.tb) - 8037, times=30)
# lSub = lapply(iRange, function(x) sample(1:nrow(mEdges), x))
# 
# lCliques.ran = lapply(lSub, getCliques)

# hist(iCliques.green, prob=T); hist(iCliques.red, prob=T); hist(iCliques.ran, prob=T)
#iCliques.all.green = cliques(ig.tb.g, 3, 7)
## centrality measures
# ig.tb= getProjectedGraph(oCGbp.tb)
# table(E(ig.tb)$weight)
# 
# # create sub graphs after edge pruning
# ig.tb.y = delete.edges(ig.tb, which(E(ig.tb)$weight < 1))
# ig.tb.g = delete.edges(ig.tb, which(E(ig.tb)$weight < 2))
# vcount(ig.tb); vcount(ig.tb.y); vcount(ig.tb.g)
# 
# # extract edge betweenness
# # mBet = lapply(list(ig.tb, ig.tb.y, ig.tb.g), function(x) { 
# #   x = delete_edge_attr(x, 'weight')
# #   betweenness(x, directed = F, weights = NULL)})
# # 
# # mBet = do.call(cbind, mBet)
# 
# lBet = lapply(list(ig.tb, ig.tb.y, ig.tb.g), function(x) { 
#   x = delete_edge_attr(x, 'weight')
#   edge_betweenness(x, directed = F, weights = NULL)})
# 
# i = sapply(lBet, which.max)
# E(ig.tb)[i[1]]
# E(ig.tb.y)[i[2]]
# E(ig.tb.g)[i[3]]
# 
# lBet[[1]][get.edge.ids(ig.tb, c('BCL6', 'BATF'))]
# lBet[[3]][get.edge.ids(ig.tb.g, c('EPHA1', 'FOXP1'))]
# head(sort(lBet[[3]], decreasing = T), 10)
# i = which(lBet[[3]] >= 8468)
# E(ig.tb.g)[i]

## degeneracy and triangles
itb.deg = c(max(mtb[,2]), max(mtb.y[,2]), max(mtb.g[,2]))
itb.tri = c(sum(count_triangles(ig.tb)), sum(count_triangles(ig.tb.y)), sum(count_triangles(ig.tb.g)))

plot(log(itb.tri), log(itb.deg), pch=20)
fit.tb = lm(log(itb.deg) ~ log(itb.tri) - 1)
abline(fit.tb)

## coreness test sepsis dataset
# extract the igraph object
# ig.sepsis= getProjectedGraph(oCGbp.sepsis)
# table(E(ig.sepsis)$weight)
# 
# # create sub graphs after edge pruning
# ig.sepsis.y = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 1))
# ig.sepsis.y = delete.vertices(ig.sepsis.y, which(degree(ig.sepsis.y) == 0))
# 
# ig.sepsis.g = delete.edges(ig.sepsis, which(E(ig.sepsis)$weight < 2))
# ig.sepsis.g = delete.vertices(ig.sepsis.g, which(degree(ig.sepsis.g) == 0))
# 
# vcount(ig.sepsis); vcount(ig.sepsis.y); vcount(ig.sepsis.g)
# 
# msepsis = sapply(c(degree, coreness), function(x){
#   return(x(ig.sepsis))
# })
# 
# msepsis.y = sapply(c(degree, coreness), function(x){
#   return(x(ig.sepsis.y))
# })
# 
# msepsis.g = sapply(c(degree, coreness), function(x){
#   return(x(ig.sepsis.g))
# })
# 
# plot(msepsis, pch=20); cor(msepsis)
# plot(msepsis.y, pch=20); cor(msepsis.y)
# plot(msepsis.g, pch=20); cor(msepsis.g)
# 
# ### atypical patterns
# ## largest cliques location
# i = names(unlist(largest_cliques(ig.sepsis)))
# i2 = which(rownames(msepsis) %in% i)
# plot(msepsis, pch=20, xlab='Degree', ylab='Coreness', main='Sepsis data - Red')
# points(msepsis[i2,], pch=20, col=2)
# 
# i = names(unlist(largest_cliques(ig.sepsis.y)))
# i2 = which(rownames(msepsis.y) %in% i)
# plot(msepsis.y, pch=20, xlab='Degree', ylab='Coreness', main='Sepsis data - Yellow')
# points(msepsis.y[i2,], pch=20, col=2)
# 
# i = names(unlist(largest_cliques(ig.sepsis.g)))
# i2 = which(rownames(msepsis.g) %in% i)
# plot(msepsis.g, pch=20, xlab='Degree', ylab='Coreness', main='Sepsis data - Green')
# points(msepsis.g[i2,], pch=20, col=2)
# 
# i = names(unlist(largest_cliques(ig.sepsis)))
# i2 = which(rownames(msepsis) %in% i)
# plot(msepsis, pch=20, xlab='Degree', ylab='Coreness', main='Sepsis data - Red')
# points(msepsis[i2,], pch=20, col=2)
# i = names(unlist(largest_cliques(ig.sepsis.y)))
# i2 = which(rownames(msepsis) %in% i)
# points(msepsis[i2,], pch=20, col='yellow')
# i = names(unlist(largest_cliques(ig.sepsis.g)))
# i2 = which(rownames(msepsis) %in% i)
# points(msepsis[i2,], pch=20, col='green')
# 
# ## degeneracy and triangles
# isepsis.deg = c(max(msepsis[,2]), max(msepsis.y[,2]), max(msepsis.g[,2]))
# isepsis.tri = c(sum(count_triangles(ig.sepsis)), sum(count_triangles(ig.sepsis.y)), sum(count_triangles(ig.sepsis.g)))
# 
# plot(log(isepsis.tri), log(isepsis.deg), pch=20)
# fit.sepsis = lm(log(isepsis.deg) ~ log(isepsis.tri) - 1)
# abline(fit.sepsis)

### combine the data for triangles and degeneracy
iDegeneracy = c(itb.deg, isepsis.deg)
iTriangles = c(itb.tri, isepsis.tri)

plot(log(iTriangles), log(iDegeneracy), pch=20, xlim=c(0, 20), ylim=c(0,6))
fit.joined = lm(log(iDegeneracy) ~ log(iTriangles) - 1)
summary(fit.joined)
abline(fit.joined)
abline(0, 0.3, col=2)

### generate some random graphs
mGenerateRandomGraph = function(v){
  mRandom = matrix(NA, nrow = 20, ncol=2)
  colnames(mRandom) = c('tri', 'dege')
  p = runif(20, 0.01, 0.2)
  for (i in 1:20){
    ig.r1 = erdos.renyi.game(v, p.or.m = p[i], type='gnp')
    # count coreness
    mr1 = sapply(c(degree, coreness), function(x){
      return(x(ig.r1))
    })
    mRandom[i,] = c(log(sum(count_triangles(ig.r1))), log(max(mr1[,2])))
  }
  return(mRandom)
}

mRan.ig.tb = mGenerateRandomGraph(vcount(ig.tb))
points(mRan.ig.tb, pch=20, col=3)

mRan.ig.tb.y = mGenerateRandomGraph(vcount(ig.tb.y))
points(mRan.ig.tb.y, pch=20, col=3)

mRan.ig.tb.g = mGenerateRandomGraph(vcount(ig.tb.g))
points(mRan.ig.tb.g, pch=20, col=3)

mRan.ig.sepsis = mGenerateRandomGraph(vcount(ig.sepsis))
points(mRan.ig.sepsis, pch=20, col=2)

mRan.ig.sepsis.y = mGenerateRandomGraph(vcount(ig.sepsis.y))
points(mRan.ig.sepsis.y, pch=20, col=3)

mRan.ig.sepsis.g = mGenerateRandomGraph(vcount(ig.sepsis.g))
points(mRan.ig.sepsis.g, pch=20, col=3)

iTriangles.ran = c((mRan.ig.tb[,'tri']),
                   (mRan.ig.tb.y[,'tri']),
                   (mRan.ig.tb.g[,'tri']),
                   (mRan.ig.sepsis[,'tri']),
                   (mRan.ig.sepsis.y[,'tri']),
                   (mRan.ig.sepsis.g[,'tri']))

iDegeneracy.ran = c((mRan.ig.tb[,'dege']),
                   (mRan.ig.tb.y[,'dege']),
                   (mRan.ig.tb.g[,'dege']),
                   (mRan.ig.sepsis[,'dege']),
                   (mRan.ig.sepsis.y[,'dege']),
                   (mRan.ig.sepsis.g[,'dege']))
iTriangles = log(iTriangles)
iDegeneracy = log(iDegeneracy)

plot(iTriangles, iDegeneracy, pch=20, 
     xlim=c(min(c(iTriangles.ran, iTriangles)), max(c(iTriangles.ran, iTriangles))),
     ylim=c(min(c(iDegeneracy.ran, iDegeneracy)),max(c(iDegeneracy.ran, iDegeneracy))))
fit.joined = lm(iDegeneracy ~ iTriangles - 1)
summary(fit.joined)
abline(fit.joined)
abline(0, 0.3, col=1, lty=2)

fit.ran = lm(iDegeneracy.ran ~ iTriangles.ran - 1)
summary(fit.ran)
points(iTriangles.ran, iDegeneracy.ran, pch=20, col=2)
abline(fit.ran, col=2)








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


