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

dfGraph = dfPathways[,c('Gene', 'Pathway')]
dfGraph = na.omit(dfGraph)
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)

# select subset of genes
# this will reduce your initial gene list as a lot of genes actually have 
# no annotations, you may use a low level database like GO if you think 
# that you lose too many genes
dfGraph = dfGraph[dfGraph$Gene %in% colnames(mCounts), ]
n = unique(dfGraph$Gene)
mCounts = mCounts[,n]
print(paste('Total number of genes with terms', length(n)))
dim(mCounts)




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



