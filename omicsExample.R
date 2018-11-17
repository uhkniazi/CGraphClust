# File: omicsExample.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 17/11/2018
# Desc: analysis of omics example data sets

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

lData.train = f_LoadObject('testData/lData.train.rds')
names(lData.train)

lData.test = f_LoadObject('testData/lData.test.rds')
names(lData.test)

mCounts.train = t(lData.train$data)
mCounts.test = t(lData.test$data)

str(mCounts.test)
str(mCounts.train)

load('testData/dfPathways.rds')

table(dfPathways$Database)

## separate the pathways and do one database at a time
dfGraph = dfPathways[dfPathways$Database == 'Reactome',c('Gene', 'Pathway')]
dfGraph$Gene = as.character(dfGraph$Gene)
dfGraph$Pathway = as.character(dfGraph$Pathway)
str(dfGraph)

library(org.Hs.eg.db)
df = AnnotationDbi::select(org.Hs.eg.db, keys=dfGraph$Gene, keytype = 'ENTREZID', columns =  'SYMBOL')
dfGraph$Gene = df$SYMBOL
# select subset of genes from our 2 matrices
dfGraph = dfGraph[dfGraph$Gene %in% unique(c(colnames(mCounts.test), colnames(mCounts.train))), ]
dfGraph = na.omit(dfGraph)
dim(dfGraph)
length(unique(dfGraph$Gene))
length(unique(c(colnames(mCounts.test), colnames(mCounts.train))))

oCGbp.reactome = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.reactome))$weight)
plot.projected.graph(oCGbp.reactome, cDropEdges = c('red', 'yellow'), bDropOrphans = T)

dfGraph = AnnotationDbi::select(org.Hs.eg.db, 
                                unique(c(colnames(mCounts.test), colnames(mCounts.train)))
                                , 'GO', 'SYMBOL')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'BP',]
dfGraph = dfGraph[,c('SYMBOL', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$SYMBOL))

oCGbp.gobp = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.gobp))$weight)
plot.projected.graph(oCGbp.gobp, cDropEdges = c('red', 'yellow'), bDropOrphans = T)


dfGraph = AnnotationDbi::select(org.Hs.eg.db, 
                                unique(c(colnames(mCounts.test), colnames(mCounts.train)))
                                , 'GO', 'SYMBOL')
dfGraph = dfGraph[dfGraph$ONTOLOGY == 'MF',]
dfGraph = dfGraph[,c('SYMBOL', 'GO')]
dfGraph = na.omit(dfGraph)
str(dfGraph)
length(unique(dfGraph$SYMBOL))

oCGbp.gomf = CGraph.bipartite(dfGraph)
table(E(getProjectedGraph(oCGbp.gomf))$weight)
plot.projected.graph(oCGbp.gomf, cDropEdges = c('red', 'yellow'), bDropOrphans = T)


# create 2 correlation graphs
oCGcor.train = CGraph.cor(mCor = cor(mCounts.train))
table(E(getProjectedGraph(oCGcor.train))$weight)

oCGcor.test = CGraph.cor(mCor = cor(mCounts.test))
table(E(getProjectedGraph(oCGcor.test))$weight)

ig = CGraph.union(getProjectedGraph(oCGcor.train),
                  getProjectedGraph(oCGcor.test),
                  getProjectedGraph(oCGbp.reactome),
                  getProjectedGraph(oCGbp.gobp),
                  getProjectedGraph(oCGbp.gomf))
table(E(ig)$weight)

ig = delete.edges(ig, which(E(ig)$weight < 3))
vcount(ig)
ecount(ig)
ig = delete.vertices(ig, which(degree(ig) == 0))
vcount(ig)
plot(ig, vertex.label=NA, vertex.size=2, layout=layout_with_fr(ig, weights = E(ig)$weight), vertex.frame.color=NA)

pdf('temp/graph.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
set.seed(123)
plot(ig, vertex.label=names(V(ig)), vertex.label.cex=0.1, vertex.size=2, vertex.frame.color=NA, 
     edge.color='darkgrey', edge.width=0.5, layout=layout_with_fr)
dev.off(dev.cur())
