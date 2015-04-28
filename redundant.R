
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
