# CGraphClust.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 25/4/2015
# Desc: class to create a igraph and hclust object based on 2 criteria: 1) shared 
#       properties or connections with type 2 vertices in a bipartite graph.
#       2) positive correlation value

library(methods)

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
CGraphClust = function(dfGraph, mCor, iCorCut=0.5){
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
  hist(w2, prob=T, breaks=s, main='distribution of obs to exp ratios', 
       xlab='square root obs to exp ratio', ylab='')
  r = round(r)
  dp = dpois(r[1]:r[2], lambda = median(w2))
  dn = dnbinom(r[1]:r[2], size = median(w2), mu = median(w2))
  lines(r[1]:r[2], dp, col='red', type='b')
  lines(r[1]:r[2], dn, col='blue', type='b')
  legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))
  
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
  l = list(oIGProj, oIGcor)  
  ig.1 = graph.intersection(l)
  # set observed to expected ratio as weight
  E(ig.1)$weight = E(ig.1)$ob_to_ex
  d = degree(ig.1)
  # delete any orphan edges
  ig.1 = delete.vertices(ig.1, which(d == 0))  
  
  ## remove small components
  cl = clusters(ig.1)
  t = log(cl$csize)
  r = range(t)
  s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
  r[1] = floor(r[1])
  r[2] = ceiling(r[2])
  # which distribution can approximate the distribution of cluster sizes
  hist(t, prob=T, main='distribution of cluster sizes', breaks=s,
       xlab='log size', ylab='')
  # try negative binomial and poisson distributions
  # parameterized on the means
  dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
  dp = dpois(r[1]:r[2], mean(t))
  lines(r[1]:r[2], dn, col='black', type='b')
  lines(r[1]:r[2], dp, col='red', type='b')
  legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
  # a poisson distribution with mean(t) fits well - use this as cutoff
  # however a negative binomial will adjust for overdispertion, try both perhaps
  ## EDIT HERE to get larger clusters i = round(exp(qpois(0.05, mean(t), lower.tail = F)))
  #i = round(exp(qnbinom(0.05, size = mean(t), mu = mean(t), lower.tail = F)))
  i = 6
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
  # NOTE: if number of edges in the graph larger than 3000 or so then
  # it may take too long or crash the system, so put in a safety check here
  # and choose a different community finding algorithm
  com = NULL
  if (ecount(ig.1) > 3000) {
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


# plot heatmap of cluster
setGeneric('plot.heatmap', def = function(obj, mCounts, ivScale = c(-3, 3), ...) standardGeneric('plot.heatmap'))
setMethod('plot.heatmap', signature='CGraphClust', definition = function(obj, mCounts, ivScale = c(-3, 3), ...){
  if (!require(NMF)) stop('R package NMF needs to be installed.')
  n = V(getFinalGraph(obj))$name
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
setGeneric('plot.heatmap.means', def = function(obj, mCounts, ivScale = c(-3, 3), ...) standardGeneric('plot.heatmap.means'))
setMethod('plot.heatmap.means', signature='CGraphClust', definition = function(obj, mCounts, ivScale = c(-3, 3), ...){
  if (!require(NMF)) stop('R package NMF needs to be installed.')
  n = V(getFinalGraph(obj))$name
  mCounts = mCounts[rownames(mCounts) %in% n,]
  hc = getHclust(obj)
  l = hc$labels
  memb = getClusterLabels(obj)
  # reorder genes according to their sequence in hc object
  mCounts = mCounts[l,]
  mCent = matrix(NA, nrow=length(unique(memb)), ncol = ncol(mCounts))
  rownames(mCent) = unique(memb)
  # loop and calculate means for each cluster
  for(a in 1:nrow(mCent)){
    i = rownames(mCent)[a]
    # if cluster has only one member
    if (sum(memb == i) == 1) {
      mCent[i,] = mCounts[memb == i,]
    } else {
      # else if more than one member, we can use mean 
      mCent[i,] = colMeans(mCounts[memb == i,])}
  }
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
  aheatmap(mCounts, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=as.factor(hc$labels), 
           annColors=NA, Colv=NA)
  # removed annColors = 'Set1'
})


# plot line graph of mean expressions in each cluster and each group
setGeneric('plot.mean.expressions', def = function(obj, mCounts, fGroups, legend.pos='topright', ...) standardGeneric('plot.mean.expressions'))
setMethod('plot.mean.expressions', signature='CGraphClust', definition = function(obj, mCounts, fGroups, legend.pos='topright', ...){
  # get the names of the genes present in the final graph
  n = V(getFinalGraph(obj))$name
  # sanity check
  if (sum(rownames(mCounts) %in% n) == 0) stop('Row names of count matrix do not match with genes')
  # subset the rows of the count matrix based on the genes
  mCounts = mCounts[rownames(mCounts) %in% n,]  
  # get the cluster labels from the cluster object
  hc = getHclust(obj)
  l = hc$labels
  memb = getClusterLabels(obj)
  # reorder genes according to their sequence in hc object
  mCounts = mCounts[l,]
  mCent = matrix(NA, nrow=length(unique(memb)), ncol = ncol(mCounts))
  rownames(mCent) = unique(memb)
  # loop and calculate means for each cluster
  for(a in 1:nrow(mCent)){
    i = rownames(mCent)[a]
    # if cluster has only one member
    if (sum(memb == i) == 1) {
      mCent[i,] = mCounts[memb == i,]
    } else {
      # else if more than one member, we can use mean 
      mCent[i,] = colMeans(mCounts[memb == i,])}
  }
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
  matplot(mPlot, type='b', lty=1, pch=20, xaxt='n', ylab='Mean Expression', ...)
  axis(1, at=1:nrow(mPlot), labels = rownames(mPlot), las=2)
  legend(legend.pos, legend = colnames(mPlot), col=1:ncol(mPlot), lty=1)
  lRet = list(means=mPlot, sd=mSD)  
  return(lRet)
})