# File: toy_example.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 27/07/2015
# Desc: create a toy example with figures

source('CGraphClust.R')
p.old = par()
set.seed(123)
# universe of possible type 2 vertices
type.2.universe = LETTERS[1:26]
type.1.universe = 1:7
# graph data frame
dfGraph = NULL
# randomly assign labels to type 1 vertices
for (i in seq_along(type.1.universe)){
  s = sample(type.2.universe, runif(1, 1, length(type.2.universe)), replace = F)
  df = data.frame(i, s)
  dfGraph = rbind(dfGraph, df)
}
## assign some labels non randomly
i = 8:10
for (x in sample(1:26, 5, replace = F)){
  s = LETTERS[x]
  df = data.frame(i, s)
  dfGraph = rbind(dfGraph, df)
}

####### create bipartite graph and projected graph
oGr.pathways = CGraph.bipartite2(dfGraph, F, F)
table(E(getProjectedGraph(oGr.pathways))$weight)

oIGbp = getBipartiteGraph(oGr.pathways)
# sanity check - is graph bipartite
if (!is.bipartite(oIGbp)) {
  stop(paste('Graph is not bipartite'))
}

# make the type 2 vertices square
fType = V(oIGbp)$type
V(oIGbp)[fType]$shape = 'circle'
V(oIGbp)[!fType]$shape = 'square'

pdf('temp/toyGraphs.pdf')
par(mar=c(1,1,1,1)+0.1)#, mfrow=c(2,2))
set.seed(123)
plot(oIGbp, layout=layout_as_bipartite, vertex.size=10)

oIG.proj = getProjectedGraph(oGr.pathways)
#E(oIG.proj)$weight = E(oIG.proj)$ob_to_ex
set.seed(123)
plot(oIG.proj, vertex.size=10, edge.label=round(E(oIG.proj)$ob_to_ex, 2), edge.label.cex=0.8, 
     layout=layout_with_fr(oIG.proj, weights = E(oIG.proj)$green), edge.width=2,
     edge.color=E(oIG.proj)$weight_cat)


iWeight = sort(E(oIG.proj)$ob_to_ex)

############# weighting functions
generate.weights = function(iCor){
  # scale value for log posterior function
  iScale = sd(iCor)
  # starting value for search - initial value
  start = c('mu'=mean(iCor))
  
  ## break the weight vector into 3 quantiles
  q1 = quantile(iCor, 0.975)
  
  q2 = quantile(iCor, 0.75)
  
  q3 = quantile(iCor, 0.5)
  
  
  ## define 3 functions with a prior for each of the quantiles
  ## model m1 - green
  m1 = function(th) dcauchy(th, q1, log = T)
  
  ## model m2 - yellow
  m2 = function(th) dcauchy(th, q2, log = T)
  
  ## model m3 - red
  m3 = function(th) dcauchy(th, q3, log = T)
  
  ## define an array that represents number of models in our parameter space
  ## each index has a prior weight/probability of being selected
  ## this can be thought of coming from a categorical distribution 
  ## moved to the arguments section of the constructor
  ##mix.prior = c(m1=3/9 ,m2= 3/9 ,m3= 3/9)
  library(LearnBayes)
  library(car)
  
  lp1 = function(theta, data){
    m = theta['mu']
    d = data # data observed
    log.lik = sum(dnorm(d, m, iScale, log=T))
    log.prior = m1(m)
    log.post = log.lik + log.prior
    return(log.post)
  }
  
  lp2 = function(theta, data){
    m = theta['mu']
    d = data # data observed
    log.lik = sum(dnorm(d, m, iScale, log=T))
    log.prior = m2(m) 
    log.post = log.lik + log.prior
    return(log.post)
  }
  
  lp3 = function(theta, data){
    m = theta['mu']
    d = data # data observed
    log.lik = sum(dnorm(d, m, iScale, log=T))
    log.prior = m3(m) 
    log.post = log.lik + log.prior
    return(log.post)
  }
  
  mMixs = sapply(seq_along(iCor), function(x){
    data = iCor[x]
    fit_m1 = laplace(lp1, start, data)
    fit_m2 = laplace(lp2, start, data)
    fit_m3 = laplace(lp3, start, data)
    mix.post = mix.prior
    mix.post[1] = exp(fit_m1$int) * mix.prior[1] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2]
                                                    + exp(fit_m3$int) * mix.prior[3])
    
    mix.post[2] = exp(fit_m2$int) * mix.prior[2] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2]
                                                    + exp(fit_m3$int) * mix.prior[3])
    
    mix.post[3] = exp(fit_m3$int) * mix.prior[3] / (exp(fit_m1$int) * mix.prior[1] + exp(fit_m2$int) * mix.prior[2]
                                                    + exp(fit_m3$int) * mix.prior[3])
    return(mix.post)
  })
  return(mMixs)
}
mix.prior = c(m1=3/9 ,m2= 3/9 ,m3= 3/9)
mPost = generate.weights(iWeight)
i = apply(mPost, 2, which.max)
cat = c('green', 'yellow', 'red')[i]
dfPost = data.frame(iWeight, cat)

## break the weight vector into 3 quantiles
q1 = quantile(iWeight, 0.975)

q2 = quantile(iWeight, 0.75)
  
q3 = quantile(iWeight, 0.5)
  
## define 3 functions with a prior for each of the quantiles
## model m1 - green
m1 = function(th) dcauchy(th, q1, log = F)
  
## model m2 - yellow
m2 = function(th) dcauchy(th, q2, log = F)
  
## model m3 - red
m3 = function(th) dcauchy(th, q3, log = F)
  
range(iWeight)
par(p.old)
iGrid = seq(-3, 3, length.out = 500)
y = seq(0, 1, length.out = 500)
d = density(iWeight)
d$y = d$y/max(d$y)
plot(iGrid, y, type='n', ylim=c(0,1), xlab='Weight', ylab='Normalised Density',
     main='Models supported by Observed Weight')
# green
m1.d = sapply(iGrid, m1)
m1.d = m1.d / max(m1.d)
lines(iGrid, m1.d, col='green')
# yellow
m2.d = sapply(iGrid, m2)
m2.d = m2.d / max(m2.d)
lines(iGrid, m2.d, col='yellow2')
# red
m3.d = sapply(iGrid, m3)
m3.d = m3.d / max(m3.d)
lines(iGrid, m3.d, col='red')
lines(d$x, d$y)

## polygons for red
d = density(iWeight)
d$y = d$y/max(d$y)

#dx = which(d$x <= max(dfPost$iWeight[dfPost$cat == 'red']))
dx = which(d$x <= 0.15)
py = c(d$y[dx], c(rep(0, length(d$x[dx]))))
px = c(d$x[dx], rev(d$x[dx]))
polygon(px, py, col='red')

## polygons for yellow
# dx = which(d$x > max(dfPost$iWeight[dfPost$cat == 'red']) & 
#                        d$x <= max(dfPost$iWeight[dfPost$cat == 'yellow']))
# dx = which(d$x > max(dfPost$iWeight[dfPost$cat == 'red']) & 
#              d$x <= 1)
dx = which(d$x > 0.15 & 
             d$x <= 1)


py = c(d$y[dx], c(rep(0, length(d$x[dx]))))
px = c(d$x[dx], rev(d$x[dx]))
polygon(px, py, col='yellow2')

## polygon for green
# dx = which(d$x >= max(dfPost$iWeight[dfPost$cat == 'green']))
dx = which(d$x > 1)
py = c(d$y[dx], c(rep(0, length(d$x[dx]))))
px = c(d$x[dx], rev(d$x[dx]))
polygon(px, py, col='green')


# d = density(dfPost$iWeight[dfPost$cat == 'yellow'])
# d$y = d$y/max(d$y)
# py = c(d$y, c(rep(0, length(d$x))))
# px = c(d$x, rev(d$x))
# polygon(px, py, col='yellow2')


# py = c(d$y, c(rep(0, length(d$x))))
# px = c(d$x, rev(d$x))
# polygon(px, py, col='green')


dev.off(dev.cur())

# # remove the low weight edges
# w2 = E(oIG.proj)$weight
# c = qnbinom(0.05, size = median(w2), mu=median(w2), lower.tail = F)
# f = which(w2 < c)
# oIGProj = delete.edges(oIG.proj, edges = f)
# # plot the new graph
# #plot(oIGProj, vertex.size=10, edge.label=round(E(oIGProj)$ob_to_ex, 2), edge.label.cex=0.7, layout=layout_with_fr)

##################################################
## section for correlation and union 
# # groups control and treatment
# fGroups = gl(2, k = 5, labels = c('con', 'treat'))
# # generate some test data
# mCounts = matrix(NA, nrow = 10, ncol = 10, dimnames = list(1:10, fGroups))
# for (r in 1:(nrow(mCounts))){
#   theta.var = 1
#   # calculate possible values of prior
#   theta.mean = rnorm(10000, 0, sd = sqrt(theta.var))
#   # data values
#   dat.1 = 0
#   # treatment will randomly be higher or lower
#   dat.2 = sample(c(-4, 4), 1)
#   # likelihoods for each prior mean  
#   lik.1 = dnorm(dat.1, mean = theta.mean, sqrt(theta.var))
#   lik.2 = dnorm(dat.2, mean = theta.mean, sqrt(theta.var))
#   # use rejection sampling to get posterior possible means
#   post.mean.1 = sample(theta.mean, 10000, replace = T, prob = lik.1)
#   post.mean.2 = sample(theta.mean, 10000, replace = T, prob = lik.2)
#   # create new data points by sampling from the posterior mean
#   # create posterior predictive distribution
#   iSize = 5
#   dat.1.new = rep(NA, iSize)
#   dat.2.new = rep(NA, iSize)
#   for (i in 1:iSize) {
#     dat.1.new[i] = rnorm(1, sample(post.mean.1, 1), sqrt(theta.var))
#     dat.2.new[i] = rnorm(1, sample(post.mean.2, 1), sqrt(theta.var))
#   }
#   mCounts[r,] = c(dat.1.new, dat.2.new)
# }
# 
# # create correlation matrix
# mCor = cor(t(mCounts))
# hist(mCor)



### create correlation graph and intersect
# mCor = abs(mCor)
# diag(mCor) = 0
# create the graph of correlations
# m = c(m1=10, m2=9, m3=8)
# m = m/sum(m)
# oGr.cor = CGraph.cor(mCor = mCor, mix.prior = m)
# oIGcor = getProjectedGraph(oGr.cor)  #graph.adjacency(mCor, mode='min', weighted=T)
# # c = E(oIGcor)$weight
# # E(oIGcor)$cor = E(oIGcor)$weight
# # iCorCut = 0.5
# # f = which(c < iCorCut)
# # oIGcor = delete.edges(oIGcor, edges = f)
# set.seed(123)
# plot(oIGcor, vertex.size=10, edge.label=round(E(oIGcor)$cor, 2), edge.label.cex=1, layout=layout_with_fr,
#      edge.color=E(oIGcor)$weight_cat)
# 
# plot(oIGcor, vertex.size=10, edge.label=E(oIGcor)$weight, edge.label.cex=1, layout=layout_with_fr,
#      edge.color=E(oIGcor)$weight_cat)


# # intersect the 2 graphs
# ig.1 = igraph::graph.intersection(oIGProj, oIGcor)
# # set observed to expected ratio as weight
# E(ig.1)$weight = E(ig.1)$ob_to_ex
# d = degree(ig.1)
# plot(ig.1, vertex.size=10, edge.label=round(E(ig.1)$ob_to_ex, 2), edge.label.cex=0.7, layout=layout_with_fr)

# dev.off(dev.cur())
# 
# ig = CGraph.union(oIG.proj, oIGcor)
# 
# plot(ig, vertex.size=10, edge.label=E(ig)$weight, edge.label.cex=1, layout=layout_with_fr)
# table(E(ig)$weight)
# ig = delete.edges(ig, which(E(ig)$weight < 1))
# vcount(ig)
# ig = delete.vertices(ig, which(degree(ig) == 0))
# vcount(ig)
