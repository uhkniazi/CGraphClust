# Copyright (C) 2015  Umar Niazi
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

library(methods)
if (!require(igraph)) stop('CGraphClust.R: library igraph required')

##### Class CGraph
# Name: Class CGgraph
# Desc: assigns weights to one mode projection of graphs based on observed to expected probabilities of 
#       vertices of the first kind i.e. with value TRUE using the igraph library
#       Zweig, K. A., & Kaufmann, M. (2011). A systematic approach to the one-mode projection of 
#       bipartite graphs. Social Network Analysis and Mining (Vol. 1, pp. 187–218). 
#       doi:10.1007/s13278-011-0021-0

# declaration
# f = flag to identify type 1 or type 2 vertices
# r = total number of type 2 vertices
# ig = bipartite igraph object
# ig.p = projected and weighted igraph object
setClass('CGraph', slots=list(ig='ANY', r='numeric', f='logical', ig.p='ANY'))

# object constructor
CGraph.bipartite = function(dfGraph, bFilterLowDegreeType2Edges=T, bFilterWeakLinks=T,  ivWeights=c(1, 0, -1)){
  # check if igraph library present
  if (!require(igraph)) stop('R library igraph required')
  if (!require(LearnBayes)) stop('R library LearnBayes required')

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
  
  ## graph cleaning 
  if (bFilterLowDegreeType2Edges){
    # remove  type 2 terms that have low degrees, 
    # these are rare terms that add little to the association scores
    f = V(oIGbp)$type
    # degree vector of type 2 vertices
    ivDegGo = degree(oIGbp, V(oIGbp)[!f])
    # setting this at less than 2, if 2 type 1 terms share a lot of weak terms, those should be preserved
    c = names(which(ivDegGo < 2))
    v = V(oIGbp)[c]
    oIGbp = delete.vertices(oIGbp, v)
    # check if graph is bipartite
    if (!is.bipartite(oIGbp)) stop('Graph is not bipartite')
  }
  
  ############################### internal private functions
  ### called by constructor

  ############# weighting functions
  generate.weights = function(suc, trials){
    ## break the weight vector into 3 quantiles
    m1.suc = quantile(suc, 0.975)
    m1.fail = trials - m1.suc
    
    m2.suc = quantile(suc, 0.75)
    m2.fail = trials - m2.suc
    
    m3.suc = quantile(suc, 0.5)
    m3.fail = trials - m3.suc
    
    
    ## define 3 functions with a prior for each of the quantiles
    ## model m1 - green
    m1 = function(th) dbeta(th, m1.suc, m1.fail, log = T)
    
    ## model m2 - yellow
    m2 = function(th) dbeta(th, m2.suc, m2.fail, log = T)
    
    ## model m3 - red
    m3 = function(th) dbeta(th, m3.suc, m3.fail, log = T)
    
    ## define an array that represents number of models in our parameter space
    ## each index has a prior weight/probability of being selected
    ## this can be thought of coming from a categorical distribution 
    mix.prior = c(m1=3/9 ,m2= 3/9 ,m3= 3/9)
    
    library(LearnBayes)
    library(car)
    logit.inv = function(p) {exp(p)/(exp(p)+1) }
    
    mylogpost_m1 = function(theta, data){
      ## theta contains parameters we wish to track
      th = logit.inv(theta['theta'])
      success = data['suc']
      fail = data['fail']
      
      # define likelihood function
      lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
      
      # calculate log posterior
      val = lf(success, fail, th) + m1(th) 
      return(val)
    }
    
    mylogpost_m2 = function(theta, data){
      ## theta contains parameters we wish to track
      th = logit.inv(theta['theta'])
      success = data['suc']
      fail = data['fail']
      
      # define likelihood function
      lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
      
      # calculate log posterior
      val = lf(success, fail, th) + m2(th)
      return(val)
    }
    
    mylogpost_m3 = function(theta, data){
      ## theta contains parameters we wish to track
      th = logit.inv(theta['theta'])
      success = data['suc']
      fail = data['fail']
      
      # define likelihood function
      lf = function(s, f, t) return(dbinom(s, s+f, t, log=T))
      
      # calculate log posterior
      val = lf(success, fail, th) + m3(th)
      return(val)
    }
    
    # starting value for search - initial value
    start = c(theta=logit(median(suc)))
    mMixs = sapply(seq_along(suc), function(x){
      data = c(suc=suc[x], fail=trials-suc[x])
      fit_m1 = laplace(mylogpost_m1, start, data)
      fit_m2 = laplace(mylogpost_m2, start, data)
      fit_m3 = laplace(mylogpost_m3, start, data)
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
  
  # Name: CGraph.project
  # Desc: assigns a score to each edge based on the model scores
  #       the score for each model is calculated using a mixture of binomial models with beta priors
  # Args: called internally no need to do it externally, 
  #       will project on vertex with TYPE=TRUE
  CGraph.project = function(obj){
    # project the graph in one dimension and assign weights
    g.p = bipartite.projection(obj@ig, which = 'TRUE')
    w = E(g.p)$weight
    if (bFilterWeakLinks) {
      # remove low weight edges i.e. if two type 1 vertices share only one type 2 vertex
      f = which(w < 2)
      g.p = delete.edges(g.p, edges=f)
      w = E(g.p)$weight
    }
    ## assign weights and categories to weights
    mWeights = generate.weights(w, obj@r)
    i = apply(mWeights, 2, which.max)
    cat = c('green', 'yellow', 'red')[i]
    num = ivWeights[i]
    E(g.p)$weight_cat = cat
    E(g.p)$weight_projection = E(g.p)$weight
    E(g.p)$weight = num
    obj@ig.p = g.p
    return(obj)
  }
  
  ######## end internal functions called by constructor
  # create the object
  g = new('CGraph', ig=oIGbp, r = 0, f= F, ig.p=NULL)
  f = V(g@ig)$type
  # r is the total numbers of vertices of the second kind
  g@r = sum(!f)
  g@f = f
  # assign weights on one mode projection
  g = CGraph.project(g)
  return(g)
}


# object constructor 2 for correlation matrix
CGraph.cor = function(ig.template, mCor, ivWeights=c(2, 1, -2)){
  # check if igraph library present
  if (!require(igraph)) stop('R library igraph required')
  if (!require(LearnBayes)) stop('R library LearnBayes required')
  library(car)
  
  ############################### internal private functions
  ### called by constructor
  
  ############# weighting functions
  generate.weights = function(iCor){
    # scale value for log posterior function
    iScale = sd(iCor)
    # starting value for search - initial value
    start = c('mu'=mean(iCor))
    # define a grid for calculating model scores and then
    # bin the actual correlations in those grid bins
    r = range(iCor)
    iGrid = seq(floor(r[1]), ceiling(r[2]), length.out = 100)
    ## break the weight vector into 3 quantiles
    q1 = quantile(iCor, 0.975)
    
    q2 = quantile(iCor, 0.75)
    
    q3 = quantile(iCor, 0.5)
    
    
    ## define 3 functions with a prior for each of the quantiles
    ## model m1 - green
    m1 = function(th) dnorm(th, q1, log = T)
    
    ## model m2 - yellow
    m2 = function(th) dnorm(th, q2, log = T)
    
    ## model m3 - red
    m3 = function(th) dnorm(th, q3, log = T)
    
    
    ## define an array that represents number of models in our parameter space
    ## each index has a prior weight/probability of being selected
    ## this can be thought of coming from a categorical distribution 
    mix.prior = c(m1=3/9 ,m2= 3/9 ,m3= 3/9)
    
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
    
    mMixs = sapply(seq_along(iGrid), function(x){
      data = iGrid[x]
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
    # bin the actual correlation vector on the grid bin
    i = apply(mMixs, 2, which.max)
    cat = c('green', 'yellow', 'red')[i]
    cat = factor(cat, levels = c('red', 'yellow', 'green'))
    # break points for bins
    p = tapply(iGrid, cat, range)
    p = c(p$red, p$yellow[2], p$green[2])
    cat.all = as.character(cut(iCor, breaks = p, include.lowest = F, labels = levels(cat)))
    return(cat.all)
  }
  
  ######## end internal functions called by constructor
  ## create correlation matrix graph, by treating it as an adjacency matrix
  mCor = round(mCor, 3)
  diag(mCor) = 0  
  # create the graph of correlations
  oIGcor = graph.adjacency(mCor, mode='min', weighted=T)
  ## house keeping and cleaning template graph to drop edge attributes
  m = as_adjacency_matrix(ig.template)
  ig.template = graph.adjacency(m, mode = 'min', weighted = NULL)
  oIGcor = graph.intersection(ig.template, oIGcor)
  c = E(oIGcor)$weight
  E(oIGcor)$cor = c
  c = logit(abs(c))
  cat = generate.weights(c)
  i = rep(NA, length(cat))
  i[cat == 'green'] = ivWeights[1]
  i[cat == 'yellow'] = ivWeights[2]
  i[cat == 'red'] = ivWeights[3]
  E(oIGcor)$weight_cat = cat
  E(oIGcor)$weight = i
  g = new('CGraph', ig=NULL, r = 0, f= F, ig.p=oIGcor)
  return(g)
}

# data acccessor functions
setGeneric('getBipartiteGraph', function(obj)standardGeneric('getBipartiteGraph'))
setMethod('getBipartiteGraph', signature = 'CGraph', definition = function(obj){
  return(obj@ig)
})

setGeneric('getProjectedGraph', function(obj)standardGeneric('getProjectedGraph'))
setMethod('getProjectedGraph', signature = 'CGraph', definition = function(obj){
  return(obj@ig.p)
})

# simple plotting function for the graph
setGeneric('plot.projected.graph', function(obj, cDropEdges='red', bDropOrphans=T)standardGeneric('plot.projected.graph'))
setMethod('plot.projected.graph', signature = 'CGraph', definition = function(obj, cDropEdges='red', bDropOrphans=T){
  ig = getProjectedGraph(obj)
  i = which(E(ig)$weight_cat %in% cDropEdges)
  if (length(i) > 0) ig = delete.edges(ig, i)
  if (bDropOrphans) {
    c = which(degree(ig) == 0)
    ig = delete.vertices(ig, c)
  }
  p.old = par(mar=c(1,1,1,1)+0.1)
  plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$weight), vertex.frame.color=NA)
  par(p.old)
})

# function to perform unions of the graph objects
CGraph.union = function(g1, g2, ...){
  u = graph.union(g1, g2, ...)
  i = grep('weight_\\d+', edge_attr_names(u))
  ## add the weights together
  w = do.call(cbind, args = edge.attributes(u)[i])
  E(u)$weight = rowSums(w, na.rm = T)
  #g = new('CGraph', ig=NULL, r = 0, f= F, ig.p=u)
  return(u)
}

## utility functions
f_dfGetGeneAnnotation = function(cvEnterezID = NULL) {
  if (!require(org.Hs.eg.db)) stop('org.Hs.eg.db annotation library required')
  return(AnnotationDbi::select(org.Hs.eg.db, cvEnterezID, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID'))  
}


############### end class CGraph

# Name: Class CGgraphClust
# Decs: class to create a igraph and hclust object based on 2 criteria: 1) shared 
#       properties or connections with type 2 vertices in a bipartite graph.
#       2) positive correlation value
# Inherits properties from CGraph class

# declaration
# contains: hc = hclust object
#           com = community object
#           labels = the most frequent type 2 vertex shared with the members of the 
#                   community (cluster) of type 1 vertices
#           ig.p2 = projected graph after cutoffs of obs to exp frequencies
#           ig.c = correlation matrix graph
#           ig.i = intersection graph
#           CGraph object with original graph
setClass('CGraphClust', slots=list(hc='ANY', com='ANY', 
                                   labels='character', ig.p2 = 'ANY',
                                   ig.c = 'ANY', ig.i = 'ANY'), 
         contains='CGraph')


# constructor
CGraphClust = function(dfGraph, mCor, iCorCut=0.5, bSuppressPlots = T, iMinComponentSize=6, clusterMethod=cluster_walktrap){
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
  
  ## graph cleaning
  # remove  type 2 terms that have low degrees, 
  # these are rare terms that add little to the association scores
  f = V(oIGbp)$type
  # degree vector of type 2 vertices
  ivDegGo = degree(oIGbp, V(oIGbp)[!f])
  # # on a log scale it follows a poisson or negative binomial dist
  # t = log(ivDegGo)
  # r = range(t)
  # s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
  # r[1] = floor(r[1])
  # r[2] = ceiling(r[2])
  # if (!bSuppressPlots){
  #   # which distribution can approximate the frequency of reactome terms
  #   hist(t, prob=T, main='degree distribution of type 2 vertices', breaks=s,
  #        xlab='log degree', ylab='')
  #   # try negative binomial and poisson distributions
  #   # parameterized on the means
  #   dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
  #   dp = dpois(r[1]:r[2], mean(t))
  #   lines(r[1]:r[2], dn, col='black', type='b')
  #   lines(r[1]:r[2], dp, col='red', type='b')
  #   legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
  # }
  # # a poisson distribution with mean(t) fits well - use this as cutoff
  # # however a negative binomial will adjust for overdispertion, try both perhaps
  # i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
  # #i = round(exp(qnbinom(0.05, size = mean(t), mu = mean(t), lower.tail = F)))
  # c = names(which(ivDegGo>i))
  c = names(which(ivDegGo<=2))
  v = V(oIGbp)[c]
  oIGbp = delete.vertices(oIGbp, v)
  # delete any orphan type 1 vertices left behind
  # d = degree(oIGbp)
  # oIGbp = delete.vertices(oIGbp, which(d == 0))

  ## graph projection to one dimension
  # create the CGraph object and calculate obs to exp weights after projection
  obj = CGraph(oIGbp)
  # create a projection of the graph 
  # oIGProj = getProjectedGraph(obj)
  oIGProj = obj@ig.p
  ## some type 1 vertices are orphans as they don't share
  # type 2 vertices with other type 1 and will now be orphans after projection,
  # remove those
  d = degree(oIGProj)
  oIGProj = delete.vertices(oIGProj, which(d == 0))
  # switch the weights with obs to exp ratio
  E(oIGProj)$weight_old = E(oIGProj)$weight
  w = rep(0, length=ecount(oIGProj))
  w[E(oIGProj)$weight_cat == 'red'] = -1
  w[E(oIGProj)$weight_cat == 'green'] = 1
  E(oIGProj)$weight = w #E(oIGProj)$ob_to_ex
  
  ## remove low observed to expected probabilities
  w = E(oIGProj)$weight_cat
  # # choose a cutoff by modelling the distribution shape
  # # it appears that the distribution follows a power law?
  # # taking square root means we can fit a poisson or neg bin distribution
  # w2 = sqrt(w)
  # r = range(w2)
  # s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by = 1)
  # r[1] = floor(r[1])
  # r[2] = ceiling(r[2])
  # if (!bSuppressPlots){
  #   hist(w2, prob=T, breaks=s, main='distribution of obs to exp ratios', 
  #        xlab='square root obs to exp ratio', ylab='')
  #   r = round(r)
  #   dp = dpois(r[1]:r[2], lambda = median(w2))
  #   dn = dnbinom(r[1]:r[2], size = median(w2), mu = median(w2))
  #   lines(r[1]:r[2], dp, col='red', type='b')
  #   lines(r[1]:r[2], dn, col='blue', type='b')
  #   legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))
  # }
  # # NOTE: this cutoff can be changed, the lower it is the more edges in the graph
  # # use negative binomial to choose cutoff
  # c = qnbinom(0.05, size = median(w2), mu=median(w2), lower.tail = F)
  # f = which(w2 < c)
  f = which(w == 'red')
  oIGProj = delete.edges(oIGProj, edges = f)
  
  ## reweight the edges
  reweight_edges = function(g.p, r=obj@r){
    w = E(g.p)$weight_old
    ob = (w+1e-6) / r
    ## assign categories to weights
    cat = cut(ob, quantile(jitter(ob), prob=c(0, 0.5, 0.75, 1)), labels = c('red', 'yellow', 'green'))
    # # calculate expected 
    # mExp = cbind(V(g.p)[m[,1]]$prob_marginal, V(g.p)[m[,2]]$prob_marginal)
    # ex = mExp[,1] * mExp[,2]
    # E(g.p)$observed = ob
    # E(g.p)$expected = ex
    # E(g.p)$ob_to_ex = ob / ex
    E(g.p)$weight_cat = as.character(cat)
    
  }
  
  ## create correlation matrix graph, by treating it as an adjacency matrix
  diag(mCor) = 0  
  # create the graph of correlations
  oIGcor = graph.adjacency(mCor, mode='min', weighted=T)
  ## house keeping and cleaning graph
  c = E(oIGcor)$weight
  E(oIGcor)$cor = E(oIGcor)$weight
  # keep only positively correlated genes connected
  # above chosen cutoff
  f = which(c < iCorCut)
  oIGcor = delete.edges(oIGcor, edges = f)
  
  ### graph intersection
  # this function causing problems?
  #l = list(oIGProj, oIGcor)  
  #ig.1 = graph.intersection(l)
  # use non list version of function
  ig.1 = igraph::graph.intersection(oIGProj, oIGcor)
  # set observed to expected ratio as weight
  E(ig.1)$weight = E(ig.1)$ob_to_ex
  d = degree(ig.1)
  # delete any orphan edges
  ig.1 = delete.vertices(ig.1, which(d == 0))  
  
  ## remove small components
  cl = clusters(ig.1)
#   t = log(cl$csize)
#   r = range(t)
#   s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
#   r[1] = floor(r[1])
#   r[2] = ceiling(r[2])
#   # which distribution can approximate the distribution of cluster sizes
#   hist(t, prob=T, main='distribution of cluster sizes', breaks=s,
#        xlab='log size', ylab='')
#   # try negative binomial and poisson distributions
#   # parameterized on the means
#   dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
#   dp = dpois(r[1]:r[2], mean(t))
#   lines(r[1]:r[2], dn, col='black', type='b')
#   lines(r[1]:r[2], dp, col='red', type='b')
#   legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
  # a poisson distribution with mean(t) fits well - use this as cutoff
  # however a negative binomial will adjust for overdispertion, try both perhaps
  ## EDIT HERE to get larger clusters i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
  #i = round(exp(qnbinom(0.05, size = mean(t), mu = mean(t), lower.tail = F)))
  i = iMinComponentSize
  i = which(cl$csize < i)
  v = which(cl$membership %in% i)
  # delete the components that are small
  ig.1 = delete.vertices(ig.1, v = v)
  
  ## clean up the bipartite graph by removing type 2 nodes
  ## that are now redundant, as the intersected final graph has less type 1
  ## vertices than the original building of the bipartite graph. 
  # recreate the bipartite graph but only with nodes that are in our final graph
  oIGbp = getBipartiteGraph(obj)
  # get the indices for the vertices of type 1
  f = V(oIGbp)$type
  n = V(oIGbp)[f]$name
  # get names of genes present in last graph i.e. intersected graph
  n2 = V(ig.1)$name
  # intersect the names to select those not present in the bipartite graph
  i = !(n %in% n2)
  n = n[i]
  # delete these type 1 vertices from the bipartite graph
  oIGbp = delete.vertices(oIGbp, v = n)
  d = degree(oIGbp)
  # delete orphan nodes left behind (which will include some type 2 vertices)
  oIGbp = delete.vertices(oIGbp, which(d == 0))
  # reset the type flag
  f = V(oIGbp)$type
  # create communities in the graph
  # NOTE: if number of edges in the graph larger than 5000 or so then
  # it may take too long or crash the system, so put in a safety check here
  # and choose a different community finding algorithm
  com = NULL
  if (ecount(ig.1) > 5000) {
    print('Too many edges in graph using cluster_walktrap')
    com = walktrap.community(ig.1)
  } else com = clusterMethod(ig.1)
  # get the hclust object 
  hc = as.hclust(com)
  memb = membership(com)
  # variable to hold the type 2 vertex common between 
  # members of a community
  rv.g = rep(NA, length=vcount(ig.1))
  rn = V(ig.1)$name
  for (i in 1:length(unique(memb))){
    # get the type 2 names names
    nei = graph.neighborhood(oIGbp, order = 1, nodes = rn[memb == i])
    # this neighbourhood graph is a list of graphs with each 
    # graph consisting of type 2 vertices that are connected to the 
    # corresponding type 1 vertex in condition rn[memb == i]
    # go through list to get the names
    pw = sapply(seq_along(nei), function(x) V(nei[[x]])$name)
    pw = unlist(pw)
    pw = as.data.frame(table(pw))
    # assign the most frequent type 2 vertex
    rv.g[memb == i] = as.character(pw[which.max(pw$Freq), 1])
  }
  # we are ready to create the object
  return(new('CGraphClust', hc=hc, com=com, labels=rv.g, ig.p2=oIGProj,
             ig.c = oIGcor, ig.i = ig.1, obj))  
  
} # constructor


#### constructor 2 to create a new CGraphClust object based on subset of vertices
# constructor
CGraphClust.recalibrate = function(obj, ivVertexID.keep, iMinComponentSize=6, clusterMethod=cluster_walktrap){
  ig.1 = getFinalGraph(obj)
  # only keep the desired vertices
  ig.1 = induced.subgraph(ig.1, V(ig.1)[ivVertexID.keep])
  d = degree(ig.1)
  # delete any orphan edges
  ig.1 = delete.vertices(ig.1, which(d == 0))  
  
  ## remove small components
  cl = clusters(ig.1)
  i = iMinComponentSize
  i = which(cl$csize < i)
  v = which(cl$membership %in% i)
  # delete the components that are small
  ig.1 = delete.vertices(ig.1, v = v)
  
  ## clean up the bipartite graph by removing type 2 nodes
  ## that are now redundant, as the intersected final graph has less type 1
  ## vertices than the original building of the bipartite graph. 
  # recreate the bipartite graph but only with nodes that are in our final graph
  oIGbp = getBipartiteGraph(obj)
  # get the indices for the vertices of type 1
  f = V(oIGbp)$type
  n = V(oIGbp)[f]$name
  # get names of genes present in last graph i.e. intersected graph
  n2 = V(ig.1)$name
  # intersect the names to select those not present in the bipartite graph
  i = !(n %in% n2)
  n = n[i]
  # delete these type 1 vertices from the bipartite graph
  oIGbp = delete.vertices(oIGbp, v = n)
  d = degree(oIGbp)
  # delete orphan nodes left behind (which will include some type 2 vertices)
  oIGbp = delete.vertices(oIGbp, which(d == 0))
  # reset the type flag
  f = V(oIGbp)$type
  # create communities in the graph
  # NOTE: if number of edges in the graph larger than 5000 or so then
  # it may take too long or crash the system, so put in a safety check here
  # and choose a different community finding algorithm
  com = NULL
  if (ecount(ig.1) > 5000) {
    print('Too many edges in graph using cluster_walktrap')
    com = walktrap.community(ig.1)
  } else com = clusterMethod(ig.1)
  # get the hclust object 
  hc = as.hclust(com)
  memb = membership(com)
  # variable to hold the type 2 vertex common between 
  # members of a community
  rv.g = rep(NA, length=vcount(ig.1))
  rn = V(ig.1)$name
  for (i in 1:length(unique(memb))){
    # get the type 2 names names
    nei = graph.neighborhood(oIGbp, order = 1, nodes = rn[memb == i])
    # this neighbourhood graph is a list of graphs with each 
    # graph consisting of type 2 vertices that are connected to the 
    # corresponding type 1 vertex in condition rn[memb == i]
    # go through list to get the names
    pw = sapply(seq_along(nei), function(x) V(nei[[x]])$name)
    pw = unlist(pw)
    pw = as.data.frame(table(pw))
    # assign the most frequent type 2 vertex
    rv.g[memb == i] = as.character(pw[which.max(pw$Freq), 1])
  }
  # we are ready to create the object
  return(new('CGraphClust', hc=hc, com=com, labels=rv.g, ig.p2=getProjectedGraph(obj),
             ig.c = getCorrelationGraph(obj), ig.i = ig.1, ig=obj@ig, r=obj@r, f=obj@f, ig.p=obj@ig.p))  
  
} # constructor 2

####

####### constructor 3
#### constructor 3 to create a new CGraphClust object by combining two CGraphClust objects
# constructor
CGraphClust.intersect.union = function(obj1, obj2, iMinOverlap=10, iMinComponentSize=6, clusterMethod=cluster_walktrap){
  # create a union graph based on common vertices
  iVertID.obj1 = which(V(getFinalGraph(obj1))$name %in% V(getFinalGraph(obj2))$name)
  iVertID.obj2 = which(V(getFinalGraph(obj2))$name %in% V(getFinalGraph(obj1))$name)
  
  # get 2 subgraphs
  ig.s1 = induced.subgraph(getFinalGraph(obj1), V(getFinalGraph(obj1))[iVertID.obj1])
  ig.s2 = induced.subgraph(getFinalGraph(obj2), V(getFinalGraph(obj2))[iVertID.obj2])
  
  ## sanity check
  f = V(ig.s1)$name %in% V(ig.s2)$name
  if (sum(f) < iMinOverlap) stop('CGraphClust.union: overlapping vertices less than cutoff')
  
  ## create a union of the graph
  ig.1 = graph.union(ig.s1, ig.s2)
  d = degree(ig.1)
  # delete any orphan edges
  ig.1 = delete.vertices(ig.1, which(d == 0))  
  
  # get new weight 
  wt = cbind(E(ig.1)$weight_1, E(ig.1)$weight_2)
  wt = rowMeans(wt, na.rm=T)
  E(ig.1)$weight = wt
  E(ig.1)$ob_to_ex = wt
  ## remove small components
  cl = clusters(ig.1)
  i = iMinComponentSize
  i = which(cl$csize < i)
  v = which(cl$membership %in% i)
  # delete the components that are small
  ig.1 = delete.vertices(ig.1, v = v)
  
  ## clean up the bipartite graph by removing type 2 nodes
  ## that are now redundant, as the intersected final graph has less type 1
  ## vertices than the original building of the bipartite graph. 
  # recreate the bipartite graph but only with nodes that are in our final graph
  bp1 = as_data_frame(getBipartiteGraph(obj1))
  bp2 = as_data_frame(getBipartiteGraph(obj2))
  dfGraph = rbind(bp1, bp2)
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
  # get the indices for the vertices of type 1
  f = V(oIGbp)$type
  n = V(oIGbp)[f]$name
  # get names of genes present in last graph i.e. intersected graph
  n2 = V(ig.1)$name
  # intersect the names to select those not present in the bipartite graph
  i = !(n %in% n2)
  n = n[i]
  # delete these type 1 vertices from the bipartite graph
  oIGbp = delete.vertices(oIGbp, v = n)
  d = degree(oIGbp)
  # delete orphan nodes left behind (which will include some type 2 vertices)
  oIGbp = delete.vertices(oIGbp, which(d == 0))
  # reset the type flag
  f = V(oIGbp)$type
  # create communities in the graph
  # NOTE: if number of edges in the graph larger than 5000 or so then
  # it may take too long or crash the system, so put in a safety check here
  # and choose a different community finding algorithm
  com = NULL
  if (ecount(ig.1) > 5000) {
    print('Too many edges in graph using cluster_walktrap')
    com = walktrap.community(ig.1)
  } else com = clusterMethod(ig.1)
  # get the hclust object 
  hc = as.hclust(com)
  memb = membership(com)
  # variable to hold the type 2 vertex common between 
  # members of a community
  rv.g = rep(NA, length=vcount(ig.1))
  rn = V(ig.1)$name
  for (i in 1:length(unique(memb))){
    # get the type 2 names names
    nei = graph.neighborhood(oIGbp, order = 1, nodes = rn[memb == i])
    # this neighbourhood graph is a list of graphs with each 
    # graph consisting of type 2 vertices that are connected to the 
    # corresponding type 1 vertex in condition rn[memb == i]
    # go through list to get the names
    pw = sapply(seq_along(nei), function(x) V(nei[[x]])$name)
    pw = unlist(pw)
    pw = as.data.frame(table(pw))
    # assign the most frequent type 2 vertex
    rv.g[memb == i] = as.character(pw[which.max(pw$Freq), 1])
  }
  # we are ready to create the object
  return(new('CGraphClust', hc=hc, com=com, labels=rv.g, ig.p2=getProjectedGraph(obj1),
             ig.c = getCorrelationGraph(obj1), ig.i = ig.1, ig=obj1@ig, r=obj1@r, f=obj1@f, ig.p=obj1@ig.p))  
  
} # constructor 3
### end constructor 3



# data acccessor functions

# projected one dimension graph after removing low obs to exp ratios
setMethod('getProjectedGraph', signature = 'CGraphClust', definition = function(obj){
  return(obj@ig.p2)
})

# graph of correlation matrix after cleanup
setGeneric('getCorrelationGraph', function(obj)standardGeneric('getCorrelationGraph'))
setMethod('getCorrelationGraph', signature = 'CGraphClust', definition = function(obj){
  return(obj@ig.c)
})

# final graph after intersection of correlation and projected graph
setGeneric('getFinalGraph', function(obj)standardGeneric('getFinalGraph'))
setMethod('getFinalGraph', signature = 'CGraphClust', definition = function(obj){
  return(obj@ig.i)
})


# labels for type 2 vertices common in clusters
setGeneric('getClusterLabels', function(obj)standardGeneric('getClusterLabels'))
setMethod('getClusterLabels', signature = 'CGraphClust', definition = function(obj){
  return(obj@labels)
})


# the hclust object for the clusters 
setGeneric('getHclust', function(obj)standardGeneric('getHclust'))
setMethod('getHclust', signature = 'CGraphClust', definition = function(obj){
  return(obj@hc)
})

# the community object after performing community searching
setGeneric('getCommunity', function(obj)standardGeneric('getCommunity'))
setMethod('getCommunity', signature = 'CGraphClust', definition = function(obj){
  return(obj@com)
})

# get mapping of type 1 vertices to cluster
setGeneric('getClusterMapping', function(obj)standardGeneric('getClusterMapping'))
setMethod('getClusterMapping', signature = 'CGraphClust', definition = function(obj){
  hc = getHclust(obj)
  df = data.frame(type.1=hc$labels, type.2=getClusterLabels(obj))
  return(df)
})


# get the largest cliques in the graph
setGeneric('getLargestCliques', function(obj)standardGeneric('getLargestCliques'))
setMethod('getLargestCliques', signature = 'CGraphClust', definition = function(obj){
  ig = getFinalGraph(obj)
  return(largest_cliques(ig))
})

# get the matrix for the cluster 
setGeneric('getClusterMatrix', def = function(obj, mCounts, csClustLabel) standardGeneric('getClusterMatrix'))
setMethod('getClusterMatrix', signature='CGraphClust', definition = function(obj, mCounts, csClustLabel){
  n = V(getClusterSubgraph(obj, csClustLabel))$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('getClusterMatrix: Row names of count matrix do not match with genes')
  mCounts = mCounts[rownames(mCounts) %in% n,]
  return(mCounts)
})


# simple plotting function for the graph
setGeneric('plot.final.graph', function(obj)standardGeneric('plot.final.graph'))
setMethod('plot.final.graph', signature = 'CGraphClust', definition = function(obj){
  ig = getFinalGraph(obj)
  p.old = par(mar=c(1,1,1,1)+0.1)
  plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$ob_to_ex), vertex.frame.color=NA)
  par(p.old)
})

setGeneric('plot.centrality.graph', function(obj, iQuantile=0.95)standardGeneric('plot.centrality.graph'))
setMethod('plot.centrality.graph', signature = 'CGraphClust', definition = function(obj, iQuantile=0.95){
  ig = getFinalGraph(obj)
  mCent = mPrintCentralitySummary(obj)
  V(ig)$color = 'lightgrey'
  # get vertices with highest 5% of degrees
  d = mCent[,'degree']
  c = quantile(d, iQuantile)
  d = d[d >= c]
  n = names(d)
  V(ig)[n]$color = 'blue'
  
  # betweenness
  d = mCent[,'betweenness']
  c = quantile(d, iQuantile)
  d = d[d >= c]
  n = names(d)
  V(ig)[n]$color = 'red'
  
  # hub
  d = mCent[,'hub']
  c = quantile(d, iQuantile)
  d = d[d >= c]
  n = names(d)
  V(ig)[n]$color = 'blue'
  
  # closeness
  d = mCent[,'closeness']
  c = quantile(d, iQuantile)
  d = d[d >= c]
  n = names(d)
  V(ig)[n]$color = 'orange'
  # plot the graph
  p.old = par(mar=c(1,1,1,1)+0.1)
  plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$ob_to_ex), vertex.frame.color=NA)
  legend('topright', legend = c('Degree', 'Hub', 'Betweenness', 'Closeness'), fill = c('blue', 'blue', 'red', 'orange'))
  par(p.old)
  return(ig)
})

# simple plotting function for the graph to highlight cliques
setGeneric('plot.graph.clique', function(obj)standardGeneric('plot.graph.clique'))
setMethod('plot.graph.clique', signature = 'CGraphClust', definition = function(obj){
  cl = getLargestCliques(obj)
  ig = getFinalGraph(obj)
  c = rainbow(length(cl)+1)
  V(ig)$color = c[length(c)]
  # cliques are a list, so choose colours based on number of cliques
  for (i in 1:length(cl)) V(ig)[cl[[i]]]$color = c[i]
  p.old = par(mar=c(1,1,1,1)+0.1)
  plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$ob_to_ex), vertex.frame.color=NA)
  par(p.old)
  return(ig)
})


# get the marginal of each cluster based on the count matrix
setGeneric('getClusterMarginal', def = function(obj, mCounts, bScaled=FALSE) standardGeneric('getClusterMarginal'))
setMethod('getClusterMarginal', signature='CGraphClust', definition = function(obj, mCounts, bScaled=FALSE){
  # internal function
  f_edge.score = function(cluster, mCl){
    g = getClusterSubgraph(obj, cluster)
    # get the list of edges
    el = get.edgelist(g)
    # for each edge list row, get the 
    # expression matrix rows
    mRet = apply(el, 1, function(x){
      m = colSums(mCl[x,])
      fac = 10*log10(E(g)[get.edge.ids(g, x)]$ob_to_ex)
      # multiply the weight by the weight factor
      return(m * fac)
#        return(m)
    })
    mRet = t(mRet)
    return(mRet)    
  }
  n = V(getFinalGraph(obj))$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('getClusterMarginal: Row names of count matrix do not match with genes')
  mCounts = mCounts[rownames(mCounts) %in% n,]
  hc = getHclust(obj)
  l = hc$labels
  memb = getClusterLabels(obj)
  # reorder genes according to their sequence in hc object
  mCounts = mCounts[l,]
  # center and scale the data across genes before adding
  # this highlights the effects on the bar plots better showing directions of effect
  if (bScaled) mCounts = t(scale(t(mCounts)))
  # get marginal data for each cluster
  mCent = matrix(NA, nrow=length(unique(memb)), ncol = ncol(mCounts))
  rownames(mCent) = unique(memb)
  colnames(mCent) = colnames(mCounts)
  # loop and calculate marginal for each cluster
  for(a in 1:nrow(mCent)){
    i = rownames(mCent)[a]
    # if cluster has only one or two members then remove from analysis
    # issue1 - if 2 members in cluster then only one edge, which crashes the 
    # colSumns function in f_edge.score function, and a subgraph can have no edges
    if (sum(memb == i) <= 2 || ecount(getClusterSubgraph(obj, i)) <= 1) {
      # mCent[i,] = mCounts[memb == i,]
      next;
    } else {
      # else if more than one member, we can use mean 
      ## make change here for using edges 
      #mCent[i,] = colMeans(mCounts[memb == i,])}
      mCent[i,] = colMeans(f_edge.score(i, mCounts[memb==i,]))}
  }
  # remove rows with NA in it, i.e. clusters with only one member
  mCent = na.omit(mCent)
  return(mCent)
})

# plot heatmap of cluster
setGeneric('plot.heatmap.all', def = function(obj, mCounts, ivScale = c(-3, 3), ...) standardGeneric('plot.heatmap.all'))
setMethod('plot.heatmap.all', signature='CGraphClust', definition = function(obj, mCounts, ivScale = c(-3, 3), ...){
  if (!require(NMF)) stop('R package NMF needs to be installed.')
  n = V(getFinalGraph(obj))$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('plot.heatmap.all: Row names of count matrix do not match with genes')
  mCounts = mCounts[rownames(mCounts) %in% n,]
  # scale across the rows
  mCounts = t(mCounts)
  mCounts = scale(mCounts)
  mCounts = t(mCounts)
  # threshhold the values
  mCounts[mCounts < ivScale[1]] = ivScale[1]
  mCounts[mCounts > ivScale[2]] = ivScale[2]
  # draw the heatmap  color='-RdBu:50'
  aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = as.dendrogram(getHclust(obj)),
           annRow=as.factor(getClusterLabels(obj)), annColors=NA, Colv=NA)
})


# plot heatmap of cluster means
setGeneric('plot.heatmap.marginal', def = function(obj, mCounts, ivScale = c(-3, 3), ...) standardGeneric('plot.heatmap.marginal'))
setMethod('plot.heatmap.marginal', signature='CGraphClust', definition = function(obj, mCounts, ivScale = c(-3, 3), ...){
  if (!require(NMF)) stop('R package NMF needs to be installed.')
  mCent = getClusterMarginal(obj, mCounts, bScaled = F)
  mCounts = mCent  
  # scale across the rows
  mCounts = t(mCounts)
  mCounts = scale(mCounts)
  mCounts = t(mCounts)
  # threshhold the values
  mCounts[mCounts < ivScale[1]] = ivScale[1]
  mCounts[mCounts > ivScale[2]] = ivScale[2]
  # draw the heatmap
  hc = hclust(dist(mCounts))
  #   aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=as.factor(hc$labels), 
  #            annColors=NA, Colv=NA)
  # removed annColors = 'Set1'
  aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, Colv=NA)
})

# plot heatmap of significant clusters only
setGeneric('plot.heatmap.significant.clusters', def = function(obj, mCounts, fGroups, bStabalize = F, ivScale = c(-3, 3), ...) standardGeneric('plot.heatmap.significant.clusters'))
setMethod('plot.heatmap.significant.clusters', signature='CGraphClust', definition = function(obj, mCounts, fGroups, bStabalize = F, ivScale = c(-3, 3), ...){
  if (!require(NMF)) stop('R package NMF needs to be installed.')
  #   # stabalize the data before performing DE
  #   if (bStabalize){
  #     mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
  #     colnames(mCounts) = fGroups
  #   }  
  # get significant clusters which also gives the marginals
  mCent = getSignificantClusters(obj, mCounts, fGroups, ...)$clusters
  mCounts = mCent  
  # stabalize the significant cluster marginals
  if (bStabalize){
    mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
    colnames(mCounts) = fGroups
  }
  # scale across the rows
  mCounts = t(mCounts)
  mCounts = scale(mCounts)
  mCounts = t(mCounts)
  # threshhold the values
  mCounts[mCounts < ivScale[1]] = ivScale[1]
  mCounts[mCounts > ivScale[2]] = ivScale[2]
  # draw the heatmap
  hc = hclust(dist(mCounts))
  #   aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=as.factor(hc$labels), 
  #            annColors=NA, Colv=NA)
  # removed annColors = 'Set1'
  aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, Colv=NA)
})

# plots heatmap for one cluster
setGeneric('plot.heatmap.cluster', def = function(obj, mCounts, csClustLabel, ivScale = c(-3, 3), ...) standardGeneric('plot.heatmap.cluster'))
setMethod('plot.heatmap.cluster', signature='CGraphClust', definition = function(obj, mCounts, csClustLabel, ivScale = c(-3, 3), ...){
  if (!require(NMF)) stop('R package NMF needs to be installed.')
  n = V(getClusterSubgraph(obj, csClustLabel))$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('plot.heatmap.cluster: Row names of count matrix do not match with genes')
  mCounts = mCounts[rownames(mCounts) %in% n,]
  # scale before plotting, across genes i.e. rows
  mCounts = t(scale(t(mCounts)))
  # threshhold the values
  mCounts[mCounts < ivScale[1]] = ivScale[1]
  mCounts[mCounts > ivScale[2]] = ivScale[2]
  # draw the heatmap
  hc = hclust(dist(mCounts))
  aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
           annColors=NA, Colv=NA)
  # removed annColors = 'Set1'
})

# plot line graph of mean expressions in each cluster and each group
setGeneric('plot.mean.expressions', def = function(obj, mCounts, fGroups, legend.pos='topright', ...) standardGeneric('plot.mean.expressions'))
setMethod('plot.mean.expressions', signature='CGraphClust', definition = function(obj, mCounts, fGroups, legend.pos='topright', ...){
  mCent = getClusterMarginal(obj, mCounts, bScaled = FALSE)
  #mCent = t(scale(t(mCent)))
  # plot the means for each level of the factor fGroups
  mPlot = matrix(NA, nrow = nrow(mCent), ncol = length(unique(fGroups)), 
                 dimnames = list(rownames(mCent), as.character(unique(fGroups))) )
  mSD = matrix(NA, nrow = nrow(mCent), ncol = length(unique(fGroups)), 
               dimnames = list(rownames(mCent), as.character(unique(fGroups))) )
  for(i in 1:nrow(mPlot)){
    # get the mean for each factor level and assign to the plot matrix
    mPlot[i,] = tapply(mCent[i,], INDEX = fGroups, FUN = mean)
    mSD[i,] = tapply(mCent[i,], INDEX = fGroups, FUN = sd)
  }
  # select number of colours
  c = c('black', 'red', 'darkgreen')
  if (ncol(mPlot) > 3)  c = rainbow(ncol(mPlot))
  # plot the matrix
  matplot(mPlot, type='b', lty=1, pch=20, xaxt='n', col=c, ylab='Mean Expression', ...)
  axis(1, at=1:nrow(mPlot), labels = rownames(mPlot), las=2)
  # if there are more than 3 factors then plot legend separately
  if (ncol(mPlot) > 3){
    plot.new()
    legend('center', legend = colnames(mPlot), col=c, lty=1, lwd=2)
  } else legend(legend.pos, legend = colnames(mPlot), col=c, lty=1, lwd=2)
  lRet = list(means=mPlot, sd=mSD)  
  return(lRet)
})


# plot line graph of significant expressions in each cluster and each group
setGeneric('plot.significant.expressions', def = function(obj, mCounts, fGroups, legend.pos='topright', bStabalize=FALSE, ...) standardGeneric('plot.significant.expressions'))
setMethod('plot.significant.expressions', signature='CGraphClust', definition = function(obj, mCounts, fGroups, legend.pos='topright', bStabalize=FALSE, ...){
  #   # stabalize the data before performing DE
  #   if (bStabalize){
  #     mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
  #     colnames(mCounts) = fGroups
  #   }  
  # get significant clusters
  mCent = getSignificantClusters(obj, mCounts, fGroups, ...)$clusters
  #   mCent = getClusterMarginal(obj, mCounts)
  #   # check which cluster shows significant p-values
  #   #p.vals = na.omit(apply(mCent, 1, function(x) pairwise.t.test(x, fGroups, p.adjust.method = 'BH')$p.value))
  #   #fSig = apply(p.vals, 2, function(x) any(x < 0.01))
  #   p.val = apply(mCent, 1, function(x) anova(lm(x ~ fGroups))$Pr[1])
  #   p.val = p.adjust(p.val, method = 'BH')
  #   fSig = p.val < 0.01
  #   mCent = mCent[fSig,]
  #   p.val = p.val[fSig]
  #   # reorder the matrix based on range of mean
  #   rSort = apply(mCent, 1, function(x){ m = tapply(x, fGroups, mean); r = range(m); diff(r)}) 
  #   mCent = mCent[order(rSort, decreasing = T),]
  # plot the means for each level of the factor fGroups
  # stabalize the significant clusters marginal
  if (bStabalize){
    mCent = t(apply(mCent, 1, function(x) f_ivStabilizeData(x, fGroups)))
    colnames(mCent) = fGroups
  }  
  mPlot = matrix(NA, nrow = nrow(mCent), ncol = length(unique(fGroups)), 
                 dimnames = list(rownames(mCent), as.character(unique(fGroups))) )
  mSD = matrix(NA, nrow = nrow(mCent), ncol = length(unique(fGroups)), 
               dimnames = list(rownames(mCent), as.character(unique(fGroups))) )
  for(i in 1:nrow(mPlot)){
    # get the mean for each factor level and assign to the plot matrix
    mPlot[i,] = tapply(mCent[i,], INDEX = fGroups, FUN = mean)
    mSD[i,] = tapply(mCent[i,], INDEX = fGroups, FUN = sd)
  }
  # select number of colours
  c = c('black', 'red', 'darkgreen')
  if (ncol(mPlot) > 3)  c = rainbow(ncol(mPlot))
  # plot the matrix
  matplot(mPlot, type='b', lty=1, pch=20, xaxt='n', col=c, ylab='Mean Expression', ...)
  axis(1, at=1:nrow(mPlot), labels = rownames(mPlot), las=2)
  # if there are more than 3 factors then plot legend separately
  if (ncol(mPlot) > 3){
    plot.new()
    legend('center', legend = colnames(mPlot), col=c, lty=1, lwd=2)
  } else legend(legend.pos, legend = colnames(mPlot), col=c, lty=1, lwd=2)
  lRet = list(means=mPlot, sd=mSD)  
  return(lRet)
})

# perform PCA on clusters (cor matrix of clusters) and return the PCA object
setGeneric('plot.components', def = function(obj, mCounts, fGroups, legend.pos='topright', bStabalize=TRUE, ...) standardGeneric('plot.components'))
setMethod('plot.components', signature='CGraphClust', definition = function(obj, mCounts, fGroups, legend.pos='topright', bStabalize=TRUE, ...){
#   # stabalize the data before pca
#   if (bStabalize){
#     mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
#     colnames(mCounts) = fGroups
#   }
  # compress the data for each cluster i.e. get marginals
#   mCent = getClusterMarginal(obj, mCounts, bScaled = F)
#   # center data across clusters i.e. rows
#   mCent = t(scale(t(mCent)))
#   # plot only significant clusters
#   l = getSignificantClusters(obj, mCounts, fGroups)
#   csClust = rownames(l$clusters)
#   mCent = mCent[csClust,]  
  mCent = getSignificantClusters(obj, mCounts, fGroups, ...)$clusters
  # stabalize the significant clusters marginal
  if (bStabalize){
    mCent = t(apply(mCent, 1, function(x) f_ivStabilizeData(x, fGroups)))
    colnames(mCent) = fGroups
  }  
  # center data across clusters i.e. rows
  mCent = t(scale(t(mCent)))
  pr.out = prcomp(mCent, scale=T)
  plot(pr.out$x[,1:2], pch=19, xlab='Z1', ylab='Z2')
  text(pr.out$x[,1:2], labels=rownames(pr.out$x), cex=0.6, pos=2)
  return(pr.out)
})

# plot heatmap of cluster means
setGeneric('plot.cluster.variance', def = function(obj, mCent, fGroups, log=TRUE, iDrawCount=4, ...) standardGeneric('plot.cluster.variance'))
setMethod('plot.cluster.variance', signature='CGraphClust', definition = function(obj, mCent, fGroups, log=TRUE, iDrawCount=4, ...){
  if (!require(lattice)) stop('R package lattice needs to be installed.')
  # for each cluster calculate the posterior variance
  fac = sapply(seq_along(levels(fGroups)), function(x) {
    rep(levels(fGroups)[x], times=100)
  })
  fac.1 = as.vector(fac)
  rn = rownames(mCent)
  if (length(rn) > iDrawCount) rn = rn[1:iDrawCount]
  # plot 4 clusters per panel
  dfVar = sapply(seq_along(rn), function(x){
    l = f_lpostVariance(mCent[rn[x],], fGroups)
    l2 = lapply(l, function(x) sample(x, 100, replace = T))
#     l2 = lapply(l, function(x){
#       # get the mean and se using central limit theorem and simulation
#       m = mean(x); se = sd(x)/sqrt(length(x))
#       return(rnorm(100, m, se))      
#     })
    # return log of the variance
    return(log(unlist(l2)))
    
  })
  colnames(dfVar) = rn
  dfVar = stack(as.data.frame(dfVar))
  dfVar$fac = factor(fac.1,levels = levels(fGroups)) 
  # if plot on log scale or exp scale
  if (log == TRUE){
    p = bwplot(~values | ind+fac, data=dfVar, do.out=TRUE, xlab='Log Variance') } else {
      p = bwplot(~ exp(values) | ind+fac, data=dfVar, do.out=TRUE, xlab='Variance')
    }
  print(p)
})


setGeneric('boxplot.cluster.variance', def = function(obj, mCent, fGroups, log=TRUE, iDrawCount=4, ...) standardGeneric('boxplot.cluster.variance'))
setMethod('boxplot.cluster.variance', signature='CGraphClust', definition = function(obj, mCent, fGroups, log=TRUE, iDrawCount=4, ...){
  # for each cluster calculate the posterior variance
  fac = sapply(seq_along(levels(fGroups)), function(x) {
    rep(levels(fGroups)[x], times=1000)
  })
  fac.1 = as.vector(fac)
  rn = rownames(mCent)
  if (length(rn) > iDrawCount) rn = rn[1:iDrawCount]
  # plot 4 clusters per panel
  dfVar = sapply(seq_along(rn), function(x){
    l2 = f_lpostVariance(mCent[rn[x],], fGroups)
    #l2 = lapply(l, function(x) sample(x, 100, replace = F))
    #     l2 = lapply(l, function(x){
    #       # get the mean and se using central limit theorem and simulation
    #       m = mean(x); se = sd(x)/sqrt(length(x))
    #       return(rnorm(100, m, se))      
    #     })
    # return log of the variance
    return(log(unlist(l2)))
    
  })
  colnames(dfVar) = rn
  fac = factor(fac.1,levels = levels(fGroups)) 
  for (i in seq_along(rn)){
    if (log == TRUE){
      boxplot( dfVar[,i] ~  fac, do.out=TRUE, ylab='Log Variance', main=rn[i], ...)} else {
        boxplot( exp(dfVar[,i]) ~  fac, do.out=TRUE, ylab='Variance', main=rn[i], ...)
      }
  }
})


# extract subgraph for a given cluster
setGeneric('getClusterSubgraph', def = function(obj, csClustLabel = NULL) standardGeneric('getClusterSubgraph'))
setMethod('getClusterSubgraph', signature='CGraphClust', definition = function(obj, csClustLabel = NULL){
  # check if label has been provided
  if (is.null(csClustLabel)) stop('Provide a cluster label')
  # check if label actually exists in the cluster labels list
  dfClust = getClusterMapping(obj)
  f = dfClust$type.2 %in% csClustLabel
  if (sum(f) == 0) stop('Cluster label does not match with actual labels')
  # if labels match correctly get the corresponding vertex labels
  v = as.character(dfClust$type.1[f])
  ig.f = getFinalGraph(obj)
  ig.ret = induced.subgraph(ig.f, V(ig.f)[v])
  return(ig.ret)
})

setGeneric('mPrintCentralitySummary', def = function(obj) standardGeneric('mPrintCentralitySummary'))
setMethod('mPrintCentralitySummary', signature='CGraphClust', definition = function(obj){
  # calculate 3 measures of centrality i.e. degree, closeness and betweenness
  ig.f = getFinalGraph(obj)
  deg = degree(ig.f)
  clo = closeness(ig.f)
  bet = betweenness(ig.f, directed = F)
  # calculate the page rank and authority_score
  aut = authority_score(ig.f, scale = F)
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
})

# get the top gene list
setGeneric('lGetTopVertices', function(obj, iQuantile=0.95)standardGeneric('lGetTopVertices'))
setMethod('lGetTopVertices', signature = 'CGraphClust', definition = function(obj, iQuantile=0.95){
  cl = getLargestCliques(obj)
  ig.f = getFinalGraph(obj)
  cvTop.cl = names(unlist(cl))
  cvTop.cl = unique(cvTop.cl)
  # 2 % of the top genes
  top.2 = vcount(ig.f) * (1-iQuantile)
  deg = degree(ig.f)
  clo = closeness(ig.f)
  bet = betweenness(ig.f, directed = F)
  # calculate the page rank and authority_score
  aut = authority_score(ig.f, scale = F)
  aut = aut$vector
  lTop = list()
  lTop$clique = cvTop.cl
  lTop$degree = names(sort(deg, decreasing = T)[1:top.2])
  lTop$closeness = names(sort(clo, decreasing = T)[1:top.2])
  lTop$betweenness = names(sort(bet, decreasing = T)[1:top.2])
  lTop$hub = names(sort(aut, decreasing = T)[1:top.2])
  return(lTop)
})


# # get the significant clusters matrix and the p.values
# setGeneric('getSignificantClusters', def = function(obj, mCounts, fGroups, ...) standardGeneric('getSignificantClusters'))
# setMethod('getSignificantClusters', signature='CGraphClust', definition = function(obj, mCounts, fGroups, ...){
# #   # stabalize the data before performing DE
# #   if (bStabalize){
# #     mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
# #     colnames(mCounts) = fGroups
# #   }  
#   # get the marginal of each cluster
#   mCent = getClusterMarginal(obj, mCounts, bScaled = F)
#   # check which cluster shows significant p-values
#   #p.vals = na.omit(apply(mCent, 1, function(x) pairwise.t.test(x, fGroups, p.adjust.method = 'BH')$p.value))
#   #fSig = apply(p.vals, 2, function(x) any(x < 0.01))
#   p.val = apply(mCent, 1, function(x) anova(lm(x ~ fGroups))$Pr[1])
#   #p.val = apply(mCent, 1, function(x) oneway.test(x ~ fGroups)$p.value)
#   p.val = p.adjust(p.val, method = 'bonferroni')
#   fSig = p.val < 0.01
#   mCent = mCent[fSig,]
#   p.val = p.val[fSig]
#   # reorder the matrix based on range of mean
#   rSort = apply(mCent, 1, function(x){ m = tapply(x, fGroups, mean); r = range(m); diff(r)}) 
#   mCent = mCent[order(rSort, decreasing = T),]
#   p.val = p.val[order(rSort, decreasing = T)]
#   lRet = list(clusters=mCent, p.val=p.val)  
#   return(lRet)
# })


# get the significant clusters matrix and the p.values
setGeneric('getSignificantClusters', def = function(obj, mCounts, fGroups, p.cut=0.001, ...) standardGeneric('getSignificantClusters'))
setMethod('getSignificantClusters', signature='CGraphClust', definition = function(obj, mCounts, fGroups, ...){
  # get the marginal of each cluster
  mCent = getClusterMarginal(obj, mCounts, bScaled = F)
  # check if grouping variable a factor
  if (!is.factor(fGroups)) stop('getSignificantClusters: Grouping variable not a factor')
  # create a new factor to represent sampled groups
  fac = sapply(seq_along(levels(fGroups)), function(x) {
    rep(levels(fGroups)[x], times=100000)
  })
  fac.1 = as.vector(fac)
  fac = factor(fac.1,levels = levels(fGroups)) 
  # calculate the posterior mean for each of the groups for each row of matrix
  rn = rownames(mCent)
  dfMean = sapply(seq_along(rn), function(x){
    l2 = f_lpostMean(mCent[rn[x],], fGroups)
    return(unlist(l2))
  })
  colnames(dfMean) = rn
  # get prob for comparison with base
  get.prob = function(x, f){
    # get the baseline
    cBaseline = levels(f)[1]
    ivBaseline = x[which(f == cBaseline)]
    # remaining levels/groups to compare against
    cvLevels = levels(f)[-1]
    pret = sapply(cvLevels, function(l){
      x.l = x[which(f == l)]
      x.l.m = mean(x.l)
      # calculate two sided p-value
      return(min(c(sum(ivBaseline <= x.l.m)/length(ivBaseline), sum(ivBaseline >= x.l.m)/length(ivBaseline))) * 2)
    })
   return(pret) 
  }
  # check which cluster shows significant p-values
  p.vals = apply(dfMean, 2, function(x) get.prob(x, fac))
  # error checking, if there were only two levels in the factor 
  # then p.vals will be a vector instead of a matrix and the apply method will not work
  # check if its matrix or vector
  fSig = NULL
  if (is.matrix(p.vals)) {
    fSig = apply(p.vals, 2, function(x) any(x < p.cut))
    p.vals = p.vals[,fSig]
  } else {
    fSig = p.vals < p.cut
    p.vals = p.vals[fSig]
  }
  mCent = mCent[fSig,]
  #p.vals = p.vals[,fSig]
  # reorder the matrix based on range of mean
  rSort = apply(mCent, 1, function(x){ m = tapply(x, fGroups, mean); r = range(m); diff(r)}) 
  mCent = mCent[order(rSort, decreasing = T),]
  # do another check for p.val being matrix or vector
  if (is.matrix(p.vals)) {
    p.vals = p.vals[,order(rSort, decreasing = T)]
  } else {
    p.vals = p.vals[order(rSort, decreasing = T)]
  }
  #p.vals = p.vals[,order(rSort, decreasing = T)]
  ## add an error check
  if (nrow(mCent) == 0) stop('getSignificantClusters: no significant p-values')
  lRet = list(clusters=mCent, p.val=p.vals)  
  return(lRet)
})




# returns igraph object with largest clique in given cluster
setGeneric('getLargestCliqueInCluster', def = function(obj, csClustLabel, ...) standardGeneric('getLargestCliqueInCluster'))
setMethod('getLargestCliqueInCluster', signature='CGraphClust', definition = function(obj, csClustLabel, ...){
  ig.sub = getClusterSubgraph(obj, csClustLabel)
  v.l = largest_cliques(ig.sub)
  return(induced_subgraph(ig.sub, unlist(v.l)))  
})

# makes bar plots of ordered statistics of centrality measures coloured by clusters
setGeneric('plot.centrality.diagnostics', def = function(obj, ...) standardGeneric('plot.centrality.diagnostics'))
setMethod('plot.centrality.diagnostics', signature='CGraphClust', definition = function(obj, ...){
  p.old = par()
  # get the cluster to gene mapping
  dfCluster = getClusterMapping(obj)
  colnames(dfCluster) = c('gene', 'cluster')
  # get the centrality parameters
  mCent = getCentralityMatrix(obj)
  rownames(dfCluster) = dfCluster$gene
  mCent = mCent[rownames(dfCluster),]
  dfCluster = cbind(dfCluster, mCent)
  # plot bar plots
  col = rainbow(length(unique(dfCluster$cluster)))
  csPlots = c('degree', 'closeness', 'betweenness', 'hub')
  lRet = vector('list', length(csPlots))
  names(lRet) = csPlots
  for (x in seq_along(csPlots)){#sapply(seq_along(csPlots), function(x){
    dfPlot = dfCluster[,c('cluster', csPlots[x])]
    dfPlot = dfPlot[order(dfPlot[,csPlots[x]], decreasing = F),]
    col.p = col[as.numeric(dfPlot$cluster)]
    barplot(dfPlot[,csPlots[x]], col=col.p, main=csPlots[x])
    cut.pts = quantile(dfPlot[,csPlots[x]], probs = c(0, 0.9, 1))
    groups = cut(dfPlot[,csPlots[x]], breaks = cut.pts, labels = c('[0,90]', '(90,100]'), include.lowest = T)
    dfPlot$groups = groups
    dfPlot.sub = dfPlot[dfPlot$groups == '(90,100]',]
    # plot the subplot of the last 90% quantile
    col.p = col[as.numeric(dfPlot.sub$cluster)]
    barplot(dfPlot.sub[,csPlots[x]], col=col.p, main=paste(csPlots[x], '(90,100] Quantile'))
    legend('topleft', legend = unique(dfPlot.sub$cluster), fill = col[as.numeric(unique(dfPlot.sub$cluster))])
    lRet[[csPlots[x]]] = dfPlot.sub
  }
  # return the list of data frames with the top genes and the associated clusters
  return(lRet)
})


# similar to mPrintCentralitySummary, just returns matrix instead of printing
setGeneric('getCentralityMatrix', def = function(obj) standardGeneric('getCentralityMatrix'))
setMethod('getCentralityMatrix', signature='CGraphClust', definition = function(obj){
  # calculate 3 measures of centrality i.e. degree, closeness and betweenness
  ig.f = getFinalGraph(obj)
  deg = degree(ig.f)
  clo = closeness(ig.f)
  bet = betweenness(ig.f, directed = F)
  # calculate the page rank and authority_score
  aut = authority_score(ig.f, scale = F)
  aut = aut$vector
  m = cbind(degree=deg, closeness=clo, betweenness=bet, hub=aut)
  return(m)  
})


# returns a data frame with the top vertices
setGeneric('dfGetTopVertices', def = function(obj, iQuantile=0.95) standardGeneric('dfGetTopVertices'))
setMethod('dfGetTopVertices', signature='CGraphClust', definition = function(obj, iQuantile=0.95){
  # get the top vertices
  l = lGetTopVertices(obj, iQuantile)
  top.genes = unique(unlist(l))
  dfRet = data.frame(VertexID = as.character(top.genes))
  # create a true / false vector of matches
  f = sapply(seq_along(1:length(l)), function(x) dfRet$VertexID %in% l[[x]])
  colnames(f) = names(l)
  dfRet = cbind(dfRet, f)
  return(dfRet)
})



## utility functions for data stabilization
## NOTE: Revisit these stabalization functions
f_ivStabilizeData = function(ivDat, fGroups){
  #set.seed(123)
  # if fGroups is not a factor
  if (!is.factor(fGroups)) stop('f_ivStabalizeData: Grouping variable not a factor')
  # calculate prior parameters
  #prior.ssd = sum((ivDat - mean(ivDat))^2)
  # uncomment these lines if prior parameters calculated by pooled data
#     sigma.0 = var(ivDat)
#     k.0 = length(ivDat)
#     v.0 = k.0 - 1
#     mu.0 = mean(ivDat)
  # comment these lines if using not using non-informataive prior
  sigma.0 = var(ivDat)
  k.0 = 2
  v.0 = k.0 - 1
  mu.0 = mean(ivDat)
  
  ## look at page 68 of Bayesian Data Analysis (Gelman) for formula
  sim.post = function(dat.grp){
    # calculate conjugate posterior
    n = length(dat.grp)
    k.n = k.0 + n
    v.n = v.0 + n
    y.bar = mean(dat.grp)
    s = sd(dat.grp)
    mu.n = (k.0/k.n * mu.0) + (n/k.n * y.bar)
    sigma.n = (( v.0*sigma.0 ) + ( (n-1)*(s^2) ) + ( (k.0*n/k.n)*((y.bar-mu.0)^2) )) / v.n
    #post.scale = ((prior.dof * prior.scale) + (var(dat.grp) * (length(dat.grp) - 1))) / post.dof
    ## simulate variance
    sigma = (sigma.n * v.n)/rchisq(1000, v.n)
    mu = rnorm(1000, mu.n, sqrt(sigma)/sqrt(k.n))
    return(list(mu=mu, var=sigma))
  }
  
  #   post.pred = function(theta, n){
  #     return(rnorm(n, mean(theta$mu), sd = mean(sqrt(theta$var))))    
  #   }
  
  # get a sample of the posterior means for each group
  ivDat.new = tapply(ivDat, fGroups, function(x) sample(sim.post(x)$mu, size = length(x), replace = T))
  ivDat.ret = vector('numeric', length=length(ivDat))
  # put the data in the same order as the groups
  l = levels(fGroups)
  for (i in seq_along(l)){
    temp = ivDat.new[l[[i]]]
    ivDat.ret[fGroups == l[i]] = unlist(temp)
  }
  return(ivDat.ret)  
}

f_lpostVariance = function(ivDat, fGroups){
  #set.seed(123)
  # if fGroups is not a factor
  if (!is.factor(fGroups)) stop('f_lpostVariance: Grouping variable not a factor')
  # calculate using non-informative prior parameters  
  sigma.0 = 0
  k.0 = 0
  v.0 = k.0 - 1
  mu.0 = 0
  
  ## look at page 68 of Bayesian Data Analysis (Gelman) for formula
  sim.post = function(dat.grp){
    # calculate conjugate posterior
    n = length(dat.grp)
    k.n = k.0 + n
    v.n = v.0 + n
    y.bar = mean(dat.grp)
    s = sd(dat.grp)
    mu.n = (k.0/k.n * mu.0) + (n/k.n * y.bar)
    sigma.n = (( v.0*sigma.0 ) + ( (n-1)*(s^2) ) + ( (k.0*n/k.n)*((y.bar-mu.0)^2) )) / v.n
    #post.scale = ((prior.dof * prior.scale) + (var(dat.grp) * (length(dat.grp) - 1))) / post.dof
    ## simulate variance
    sigma = (sigma.n * v.n)/rchisq(1000, v.n)
    mu = rnorm(1000, mu.n, sqrt(sigma)/sqrt(k.n))
    return(list(mu=mu, var=sigma))
  }
    
  # get a sample of the posterior variance for each group
  return(tapply(ivDat, fGroups, function(x) sim.post(x)$var))
}

f_lpostMean = function(ivDat, fGroups, sim.size=100000){
  #set.seed(123)
  # if fGroups is not a factor
  if (!is.factor(fGroups)) stop('f_lpostMean: Grouping variable not a factor')
  # calculate using non-informative prior parameters  
  sigma.0 = 0
  k.0 = 0
  v.0 = k.0 - 1
  mu.0 = 0
  
  ## look at page 68 of Bayesian Data Analysis (Gelman) for formula
  sim.post = function(dat.grp){
    # calculate conjugate posterior
    n = length(dat.grp)
    k.n = k.0 + n
    v.n = v.0 + n
    y.bar = mean(dat.grp)
    s = sd(dat.grp)
    mu.n = (k.0/k.n * mu.0) + (n/k.n * y.bar)
    sigma.n = (( v.0*sigma.0 ) + ( (n-1)*(s^2) ) + ( (k.0*n/k.n)*((y.bar-mu.0)^2) )) / v.n
    #post.scale = ((prior.dof * prior.scale) + (var(dat.grp) * (length(dat.grp) - 1))) / post.dof
    ## simulate variance
    sigma = (sigma.n * v.n)/rchisq(sim.size, v.n)
    mu = rnorm(sim.size, mu.n, sqrt(sigma)/sqrt(k.n))
    return(list(mu=mu, var=sigma))
  }
  
  # get a sample of the posterior mean for each group
  return(tapply(ivDat, fGroups, function(x) sim.post(x)$mu))
}




# f_igCalculateVertexSizes = function(ig, mCounts, fGroups, bStabalize=FALSE, iSize=NULL){
#   n = V(ig)$name
#   # sanity check
#   if (sum(rownames(mCounts) %in% n) == 0) stop('f_igCalculatevertexSizes: Row names of count matrix do not match with genes')
#   
#   # calculate fold changes function
#   lf_getFC = function(x, f, bS=FALSE){
#     # if data stabalization required
#     if (bS) x = f_ivStabilizeData(x, f)
#     r = range(tapply(x, f, mean))
#     fc = log10(r[2]) - log10(r[1])
#     return(fc)
#   }
#   if (is.null(iSize)) iSize = 4000/vcount(ig)
#   mCounts = mCounts[n,]
#   s = apply(mCounts, 1, function(x) lf_getFC(x, fGroups, bStabalize))
#   V(ig)[n]$size = s * iSize
#   return(ig)
# }

# utility function to assign colours and sizes to vertices
f_igCalculateVertexSizesAndColors = function(ig, mCounts, fGroups, bColor = FALSE, bStabalize=FALSE, iSize=NULL){
  n = V(ig)$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('f_igCalculatevertexSizesAndColors: Row names of count matrix do not match with genes')
  
  # calculate fold changes function
  lf_getFC = function(x, f, bS=FALSE){
    # if data stabalization required
    l = levels(f)
    if (bS) x = f_ivStabilizeData(x, f)
    r = tapply(x, f, mean)
    fc = log2(r[l[length(l)]]) - log2(r[l[1]])
    return(fc)
  }
  # calculate colour function
  lf_getDirection = function(x, f, bS=FALSE){
    l = levels(f)
    if (bS) x = f_ivStabilizeData(x, f)
    r = tapply(x, f, mean)
    c = ifelse(r[l[1]] < r[l[length(l)]], 'pink', 'lightblue')
    return(c)
  }
  
  if (is.null(iSize)) iSize = 4000/vcount(ig)
  mCounts = mCounts[n,]
  s = apply(mCounts, 1, function(x) lf_getFC(x, fGroups, bStabalize))
  V(ig)[n]$size = abs(s * iSize)
  # assign colours if required
  if (bColor){c = sapply(seq_along(n), function(x) lf_getDirection(mCounts[n[x], ], fGroups, bStabalize))
              V(ig)[n]$color = c
  }  
  return(ig)
}

# utility function to obtain gene information from genbank
f_csGetGeneSummaryFromGenbank = function(iID){
  if (!require(annotate) || !require(XML)) stop('CRAN packge XML and Bioconductor library annotate required')
  warning(paste('The powers that be at NCBI have been known to ban the', 
                'IP addresses of users who abuse their servers (currently', 
                'defined as less then 2 seconds between queries). Do NOT',
                'put this function in a tight loop or you may find your access revoked',
                'see help annotate::genbank for details.'))
  gene = genbank(as.numeric(iID))
  r = xmlRoot(gene)
  # count number of nodes
  iNode = length(xmlChildren(r))
  csRet = rep(NA, iNode)
  #names(csRet) = as.character(iID)
  names(csRet) = sapply(seq_along(1:iNode), function(x) r[[x]][['Entrezgene_gene']][[1]][[1]][[1]]$value)
  for (i in 1:iNode){
    x = r[[i]][['Entrezgene_summary']][[1]]
    if (is.null(x$value)) next;
    csRet[i] = x$value
  }
  return(csRet)
}
