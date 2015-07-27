# File: toy_example.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 27/07/2015
# Desc: create a toy example with figures



source('CGraphClust.R')

set.seed(123)
# universe of possible type 2 vertices
type.2.universe = LETTERS[1:26]
type.1.universe = 1:7

dfGraph = NULL
p.old = par()

for (i in seq_along(type.1.universe)){
  s = sample(type.2.universe, 5, replace = F)
  df = data.frame(i, s)
  dfGraph = rbind(dfGraph, df)
}

## assign some labels non randomly
i = 8:10#sample(1:10, size = 3, replace = F)
for (x in sample(1:26, 5, replace = F)){
  s = LETTERS[x]
  df = data.frame(i, s)
  dfGraph = rbind(dfGraph, df)
}

fGroups = gl(2, k = 5, labels = c('con', 'treat'))

# generate some test data
mCounts = matrix(NA, nrow = 10, ncol = 10, dimnames = list(1:10, fGroups))
for (r in 1:(nrow(mCounts))){
  theta.var = 1
  # calculate possible values of prior
  theta.mean = rnorm(10000, 0, sd = sqrt(theta.var))
  # data values
  dat.1 = 0
  # treatment will randomly be higher or lower
  dat.2 = sample(c(-4, 4), 1)
  # likelihoods for each prior mean  
  lik.1 = dnorm(dat.1, mean = theta.mean, sqrt(theta.var))
  lik.2 = dnorm(dat.2, mean = theta.mean, sqrt(theta.var))
  # use rejection sampling to get posterior possible means
  post.mean.1 = sample(theta.mean, 10000, replace = T, prob = lik.1)
  post.mean.2 = sample(theta.mean, 10000, replace = T, prob = lik.2)
  # create new data points by sampling from the posterior mean
  # create posterior predictive distribution
  iSize = 5
  dat.1.new = rep(NA, iSize)
  dat.2.new = rep(NA, iSize)
  for (i in 1:iSize) {
    dat.1.new[i] = rnorm(1, sample(post.mean.1, 1), sqrt(theta.var))
    dat.2.new[i] = rnorm(1, sample(post.mean.2, 1), sqrt(theta.var))
  }
  mCounts[r,] = c(dat.1.new, dat.2.new)
}

# create correlation matrix
mCor = cor(t(mCounts))
hist(mCor)
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)

####### create bipartite graph
oIGbp = graph.data.frame(dfGraph, directed = F)
# set the vertex type variable to make graph bipartite
f = rep(c(T, F), times = c(length(unique(dfGraph[,1])),length(unique(dfGraph[,2]))))
V(oIGbp)$type = f
# sanity check - is graph bipartite
if (!is.bipartite(oIGbp)) {
  stop(paste('Graph is not bipartite'))
}

# make the type 2 vertices square
fType = V(oIGbp)$type
V(oIGbp)[fType]$shape = 'circle'
V(oIGbp)[!fType]$shape = 'square'

par(mar=c(1,1,1,1)+0.1)
plot(oIGbp, layout=layout_as_bipartite, vertex.size=10)

oGr = CGraph(oIGbp)
oIG.proj = getProjectedGraph(oGr)
E(oIG.proj)$weight = E(oIG.proj)$ob_to_ex
plot(oIG.proj, vertex.size=10, edge.label=round(E(oIG.proj)$ob_to_ex, 2), edge.label.cex=0.7, layout=layout_with_fr)


