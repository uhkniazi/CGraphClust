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


## this function removes small communities, but affects how the heatmaps are made at the end
# not very appropriate, add later if required.
# remove small communities below a certain size
setGeneric('oCGdeleteSmallCommunities', function(obj, size=3)standardGeneric('oCGdeleteSmallCommunities'))
setMethod('oCGdeleteSmallCommunities', signature ='CGraphClust', definition = function(obj, size=3){
  # get the community object and cut out small communities
  com = getCommunity(obj)
  ig.1 = getFinalGraph(obj)
  i = which(sizes(com) <= size)
  m = membership(com)
  # get community member vertex names
  n = names(m[m %in% i])
  # delete these vertices from the graph and recalculate the communities
  ig.1 = delete.vertices(ig.1, n)  
  # find communities again and reassign labels
  if (ecount(ig.1) > 3000) {
    stop('Too many edges in graph to safely find communities')
  }
  com = edge.betweenness.community(ig.1)
  # get the hclust object 
  hc = as.hclust(com)
  memb = membership(com)
  # variable to hold the type 2 vertex common between 
  # members of a community
  rv.g = rep(NA, length=vcount(ig.1))
  rn = V(ig.1)$name
  for (i in 1:length(unique(memb))){
    # get the type 2 names names
    nei = graph.neighborhood(getBipartiteGraph(obj), order = 1, nodes = rn[memb == i])
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
  obj@hc = hc
  obj@com = com
  obj@labels = rv.g
  obj@ig.i = ig.1
  return(obj)
})


### code to add later if we want to do a odds ratio test for terms
f_getSignificantTerm = function(dfDat.t, dfDat.r){
  f.t = sum(dfDat.t$Freq)
  f.r = sum(dfDat.r$Freq)
  p.val = rep(NA, length=nrow(dfDat.t))
  names(p.val) = dfDat.t[,1]
  for (i in 1:length(p.val)){
    # how many times term is present and absent in test data
    t.p = dfDat.t[i,2]
    t.a = f.t - t.p
    # in reference data
    r.p = dfDat.r[as.character(dfDat.r[,1]) == as.character(dfDat.t[i,1]), 2]
    if (length(r.p) == 0) r.p = 1
    r.a = f.r - r.p
    m = matrix(c(t.p, t.a, r.p, r.a), nrow = 2, byrow = T)
    p.val[i] = fisher.test(m)$p.value
  }
  return(names(which.min(p.val)))
}

## this code will fit in the loop where we assign common terms to communities 
for (i in 1:length(unique(memb))){
  # get the type 2 names from members of cluster
  nei = graph.neighborhood(getBipartiteGraph(obj), order = 1, nodes = rn[memb == i])
  # get type 2 names from members not in cluster
  nei.all = graph.neighborhood(getBipartiteGraph(obj), order = 1, nodes = rn[memb != i])
  # this neighbourhood graph is a list of graphs with each 
  # graph consisting of type 2 vertices that are connected to the 
  # corresponding type 1 vertex in condition rn[memb == i]
  # go through list to get the names
  pw = sapply(seq_along(nei), function(x) V(nei[[x]])$name)
  pw = unlist(pw)
  pw = as.data.frame(table(pw))
  ## get names for all the rest not in the cluster
  pw.all = sapply(seq_along(nei.all), function(x) V(nei.all[[x]])$name)
  pw.all = unlist(pw.all)
  pw.all = as.data.frame(table(pw.all))
  # assign the most frequent type 2 vertex
  #rv.g[memb == i] = as.character(pw[which.max(pw$Freq), 1])
  rv.g[memb == i] = as.character(f_getSignificantTerm(pw, pw.all))
}
