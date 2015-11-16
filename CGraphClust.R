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

# CGraphClust.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 25/4/2015
# Desc: class to create a igraph and hclust object based on 2 criteria: 1) shared 
#       properties or connections with type 2 vertices in a bipartite graph.
#       2) positive correlation value

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
setClass('CGraph', slots=list(ig='ANY', r='numeric', f='logical', ig.p='ANY'))

# object constructor
CGraph = function(oGraph){
  # check if igraph library present
  if (!require(igraph)) stop('R library igraph required')
  # check if graph is bipartite
  if (!is.bipartite(oGraph)) stop('Graph is not bipartite')
  #### internal private functions
  # processing steps - called by constructor  
  # assign probabilities to vertex of first kind
  # Name: CGraph.assign.marginal.probabilities
  # Desc: assigns probabilities to each vertex of the first kind (TRUE) 
  #       based on how many times it is connected to the vertex of the 
  #       second kind i.e. degree(V1) / (total number of V-type2)
  # Args: internal function - object of CGraph class
  CGraph.assign.marginal.probabilities = function(obj){
    # vertex of the first kind will be assigned probabilities
    # based on their relations with the vertices of the second kind
    # flag to identify vertex types
    f = V(obj@ig)$type
    d = degree(obj@ig)
    d = d[f]
    # r is the total numbers of vertices of the second kind
    r = sum(!f)
    p = d/r
    V(obj@ig)[f]$prob_marginal = p
    obj@r = r
    obj@f = f
    return(obj)
  }
  
  # Name: CGraph.project
  # Desc: assigns a level of interestingness/leverage or observed to expected ratio to 
  #       each edge after graph projection on the vertex of first kind i.e. type = TRUE 
  #       Observed frequency = weight of edge / (total number of vertices of second type)
  #       i.e. how many shared vertices of type 2 are between the 2 type 1 vertices
  #       Expected frequency = how many times we expect to see them based on their 
  #       joint probability under assumption of independence. 
  #       (marginal.prob of V1 * marginal.prob of V2)
  # Args: called internally no need to do it externally, 
  #       will project on vertex with TYPE=TRUE
  CGraph.project = function(obj){
    # project the graph in one dimension and
    # assign weights based on observed to expected ratios
    g.p = bipartite.projection(obj@ig, which = 'TRUE')
    # get the matrix with rows representing each edge
    m = get.edgelist(g.p)
    w = E(g.p)$weight
    # calculate observed ratio
    # weight / r
    ob = w / obj@r
    # calculate expected 
    mExp = cbind(V(g.p)[m[,1]]$prob_marginal, V(g.p)[m[,2]]$prob_marginal)
    ex = mExp[,1] * mExp[,2]
    E(g.p)$observed = ob
    E(g.p)$expected = ex
    E(g.p)$ob_to_ex = ob / ex
    obj@ig.p = g.p
    return(obj)
  }
  
  ####
  # create the object
  g = new('CGraph', ig=oGraph, r = 0, f= F, ig.p=NULL)
  # assign marginal probabilities
  g = CGraph.assign.marginal.probabilities(g)
  # assign weights on one mode projection
  g = CGraph.project(g)
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
CGraphClust = function(dfGraph, mCor, iCorCut=0.5, bSuppressPlots = T, iMinComponentSize=6){
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
  # remove very frequent type 2 terms, as they create too many edges
  # and may hide real relationships
  f = V(oIGbp)$type
  # degree vector of type 2 vertices
  ivDegGo = degree(oIGbp, V(oIGbp)[!f])
  # on a log scale it follows a poisson or negative binomial dist
  t = log(ivDegGo)
  r = range(t)
  s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
  r[1] = floor(r[1])
  r[2] = ceiling(r[2])
  if (!bSuppressPlots){
    # which distribution can approximate the frequency of reactome terms
    hist(t, prob=T, main='degree distribution of type 2 vertices', breaks=s,
         xlab='log degree', ylab='')
    # try negative binomial and poisson distributions
    # parameterized on the means
    dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
    dp = dpois(r[1]:r[2], mean(t))
    lines(r[1]:r[2], dn, col='black', type='b')
    lines(r[1]:r[2], dp, col='red', type='b')
    legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
  }
  # a poisson distribution with mean(t) fits well - use this as cutoff
  # however a negative binomial will adjust for overdispertion, try both perhaps
  #i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
  i = round(exp(qnbinom(0.05, size = mean(t), mu = mean(t), lower.tail = F)))
  c = names(which(ivDegGo>i))
  v = V(oIGbp)[c]
  oIGbp = delete.vertices(oIGbp, v)
  # delete any orphan type 1 vertices left behind
  d = degree(oIGbp)
  oIGbp = delete.vertices(oIGbp, which(d == 0))
  
  ## graph projection to one dimension
  # create the CGraph object and calculate obs to exp weights after projection
  obj = CGraph(oIGbp)
  # create a projection of the graph 
  oIGProj = getProjectedGraph(obj)
  ## some type 1 vertices are orphans as they don't share
  # type 2 vertices with other type 1 and will now be orphans after projection,
  # remove those
  d = degree(oIGProj)
  oIGProj = delete.vertices(oIGProj, which(d == 0))
  # switch the weights with obs to exp ratio
  E(oIGProj)$weight_old = E(oIGProj)$weight
  E(oIGProj)$weight = E(oIGProj)$ob_to_ex
  
  ## remove low observed to expected probabilities
  w = E(oIGProj)$weight
  # choose a cutoff by modelling the distribution shape
  # it appears that the distribution follows a power law?
  # taking square root means we can fit a poisson or neg bin distribution
  w2 = sqrt(w)
  r = range(w2)
  s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by = 1)
  r[1] = floor(r[1])
  r[2] = ceiling(r[2])
  if (!bSuppressPlots){
    hist(w2, prob=T, breaks=s, main='distribution of obs to exp ratios', 
         xlab='square root obs to exp ratio', ylab='')
    r = round(r)
    dp = dpois(r[1]:r[2], lambda = median(w2))
    dn = dnbinom(r[1]:r[2], size = median(w2), mu = median(w2))
    lines(r[1]:r[2], dp, col='red', type='b')
    lines(r[1]:r[2], dn, col='blue', type='b')
    legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))
  }
  # NOTE: this cutoff can be changed, the lower it is the more edges in the graph
  # use negative binomial to choose cutoff
  c = qnbinom(0.05, size = median(w2), mu=median(w2), lower.tail = F)
  f = which(w2 < c)
  oIGProj = delete.edges(oIGProj, edges = f)
  
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
    print('Too many edges in graph to use edge.betweenness communities')
    com = walktrap.community(ig.1)
  } else com = edge.betweenness.community(ig.1)
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
      fac = E(g)[get.edge.ids(g, x)]$ob_to_ex
#       # multiply the weight by the weight factor
#       return(m * fac)
        return(m)
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
    # if cluster has only one member then remove from analysis
    if (sum(memb == i) == 1) {
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
  mCent = getSignificantClusters(obj, mCounts, fGroups)$clusters
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
  mCent = getSignificantClusters(obj, mCounts, fGroups)$clusters
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
  mCent = getSignificantClusters(obj, mCounts, fGroups)$clusters
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
setGeneric('plot.cluster.variance', def = function(obj, mCent, fGroups, iDrawCount=4, ...) standardGeneric('plot.cluster.variance'))
setMethod('plot.cluster.variance', signature='CGraphClust', definition = function(obj, mCent, fGroups, iDrawCount=4, ...){
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
  bwplot(~values | ind+fac, data=dfVar, do.out=TRUE, xlab='Log Variance') 
})

# plot the expression of all members of the given cluster
setGeneric('plot.cluster.expressions', def = function(obj, mCounts, fGroups, csClustLabel, ...) standardGeneric('plot.cluster.expressions'))
setMethod('plot.cluster.expressions', signature='CGraphClust', definition = function(obj, mCounts, fGroups, csClustLabel, ...){
  n = V(getClusterSubgraph(obj, csClustLabel))$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('plot.cluster.expressions: Row names of count matrix do not match with genes')
  mCounts = mCounts[rownames(mCounts) %in% n,]  
  # scale before plotting
  mPlot = scale(t(mCounts))
  c = rainbow(ncol(mPlot))
  # plot the matrix
  matplot(mPlot, type='l', lty=1, pch=20, xaxt='n', col=c, ylab='Expression', ...)
  axis(1, at=1:nrow(mPlot), labels = rownames(mPlot), las=2, cex.axis=0.5)
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


# get the significant clusters matrix and the p.values
setGeneric('getSignificantClusters', def = function(obj, mCounts, fGroups, ...) standardGeneric('getSignificantClusters'))
setMethod('getSignificantClusters', signature='CGraphClust', definition = function(obj, mCounts, fGroups, ...){
#   # stabalize the data before performing DE
#   if (bStabalize){
#     mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
#     colnames(mCounts) = fGroups
#   }  
  # get the marginal of each cluster
  mCent = getClusterMarginal(obj, mCounts, bScaled = F)
  # check which cluster shows significant p-values
  #p.vals = na.omit(apply(mCent, 1, function(x) pairwise.t.test(x, fGroups, p.adjust.method = 'BH')$p.value))
  #fSig = apply(p.vals, 2, function(x) any(x < 0.01))
  p.val = apply(mCent, 1, function(x) anova(lm(x ~ fGroups))$Pr[1])
  p.val = p.adjust(p.val, method = 'BH')
  fSig = p.val < 0.01
  mCent = mCent[fSig,]
  p.val = p.val[fSig]
  # reorder the matrix based on range of mean
  rSort = apply(mCent, 1, function(x){ m = tapply(x, fGroups, mean); r = range(m); diff(r)}) 
  mCent = mCent[order(rSort, decreasing = T),]
  p.val = p.val[order(rSort, decreasing = T)]
  lRet = list(clusters=mCent, p.val=p.val)  
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
  l = lGetTopVertices(oGr, iQuantile)
  top.genes = unique(unlist(l))
  dfRet = data.frame(VertexID = as.character(top.genes))
  # create a true / false vector of matches
  f = sapply(seq_along(1:length(l)), function(x) dfRet$VertexID %in% l[[x]])
  colnames(f) = names(l)
  dfRet = cbind(dfRet, f)
  return(dfRet)
})



## utility functions for data stabilization
# f_mCalculateLikelihoodMatrix = function(ivDat, fGroups){
#   # if fGroups is not a factor
#   if (!is.factor(fGroups)) stop('f_mCalculateLikelihoodMatrix: Grouping variable not a factor')
#   # get the parameters for the data of the data
#   v.dat = var(ivDat)
#   var.dat = tapply(ivDat, fGroups, var)
#   l.dat = tapply(ivDat, fGroups, length)
#   m.dat = tapply(ivDat, fGroups, mean)
#   se.dat = sqrt(var.dat/l.dat)
#   # sample possible values of the means from normal distribution
#   theta.mean = rnorm(10000, m.dat, sqrt(v.dat))
#   mLik = matrix(NA, nrow=length(theta.mean), ncol=length(levels(fGroups)))
#   colnames(mLik) = levels(fGroups)
#   
#   for (i in 1:length(levels(fGroups))){
#     # calculate likelihood function for this mean
#     y = dnorm(theta.mean, m.dat[i], se.dat[i])
#     # use rejection sampling to sample from posterior based on likelihood
#     mLik[,i] = sample(theta.mean, 10000, replace = T, prob=y)
#   }
#   return(mLik)
# }
# 
# f_ivStabilizeData = function(ivDat, fGroups){
#   # set seed 
#   set.seed(123)
#   mNew = f_mCalculateLikelihoodMatrix(ivDat, fGroups)
#   var.dat = tapply(ivDat, fGroups, var)
#   l.dat = tapply(ivDat, fGroups, length)
#   se.dat = sqrt(var.dat/l.dat)
#   sd.dat = tapply(ivDat, fGroups, sd)
#   ivDat.new = NULL
#   for (i in 1:length(levels(fGroups))){
#     # sample of means from the posterior means
#     m = sample(mNew[,i], size = l.dat[i], replace = T)
#     ivSam = sapply(seq_along(m), function(x) rnorm(1, m[x], se.dat[i]))
#     ivDat.new = c(ivDat.new, ivSam)
#   }
#   return(ivDat.new)
# }

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
  if (!is.factor(fGroups)) stop('f_ivStabalizeData: Grouping variable not a factor')
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


f_dfGetGeneAnnotation = function(cvEnterezID = NULL) {
  if (!require(org.Hs.eg.db)) stop('org.Hs.eg.db annotation library required')
  return(AnnotationDbi::select(org.Hs.eg.db, cvEnterezID, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID'))  
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
