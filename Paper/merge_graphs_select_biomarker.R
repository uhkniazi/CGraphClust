# File: merge_graphs_select_biomarker.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 23/11/2016
# Desc: merge the two networks and select top genes based on centrality


######### load libraries and set global variables
library(org.Hs.eg.db)
library(downloader)
source('CGraphClust.R')
# plotting parameters
p.old = par()

## load reactome data
### load the uniprot2reactome mapping obtained from
# http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x
# convert enterez ids to uniprot as Reactome database file uses UNIPROT ids
dfMap = AnnotationDbi::select(org.Hs.eg.db, dfReactome$V1, 'ENTREZID', 'UNIPROT')
dfMap = na.omit(dfMap)

## map reactome ids to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
rm(dfReactome, dfMap)
gc()
###

########### utility functions
get.reactome.name = function(csNames){
  i = which(dfReactome.sub$V2 %in% csNames)
  dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
  dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
  rownames(dfCluster.name) = NULL
  return(dfCluster.name)
}

## load the two tb datasets to merge
load('Paper/Results/ltb_atb_data.rds')
load('Paper/Results/tb_data.rds')

## genes that show association in both graphs are selected by
## intersecting the two graphs
# intersect the 2 graphs and recalibrate the graph object with the common vertices
igi = graph.intersection(getFinalGraph(ltb_atb_data$graph), getFinalGraph(tb_data$graph))
d = degree(igi)
igi = delete.vertices(igi, which(d == 0))
ecount(igi)
vcount(igi)
## ~ 224 genes are coexpressed in both the datasets
cVertID = V(igi)$name
## remake the graph using the ltb_atb graph
oGr = CGraphClust.recalibrate(ltb_atb_data$graph, cVertID)

mCounts = ltb_atb_data$matrix
fGroups = ltb_atb_data$groups

## general graph structure
## we would like to see how does the graph look like, are the clusters connected or in subgraphs
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))

## look at the top centrality genes
set.seed(1)
ig = plot.centrality.graph(oGr)

dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

dir.create('Paper/Results', showWarnings = F)
write.csv(dfTopGenes.cent, file='Paper/Results/Top_Centrality_Genes_tb_combined.csv')

## interestingly 1280215 cluster has genes from 
## the largest clique and have highest degrees
## see table for dfTopGenes.cent
## get the list of these genes
l = getLargestCliqueInCluster(oGr, '1280215')
cvTopGenes = names(V(l))

## second cluster with highest hub score
l = getLargestCliqueInCluster(oGr, '6798695')
cvTopGenes = c(cvTopGenes, names(V(l)))

# select only 30 genes for this test, as exhaustive selection doesnt work well with more than 30 genes
cvTopGenes = cvTopGenes[1:30]

## are any of these diagnostic/biomarkers
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

dfDat = data.frame(mCounts[,cvTopGenes])
fSamples = fGroups
colnames(dfDat) = f_dfGetGeneAnnotation(cvTopGenes)$SYMBOL
head(dfDat)

# merge the LTB and HC together as Healthy and ATB as disease
i = grep('LTB|HC', fSamples)
#dfDat = dfDat[i,]
rownames(dfDat) = names(fSamples)
fSamples = as.character(fSamples)
fSamples[i] = 'Healthy'
fSamples[-i] = 'Disease'
fSamples = factor(fSamples, levels = c('Healthy', 'Disease'))

par(p.old)
## variable combinations
# create a test set variable
set.seed(123)
test = sample(1:nrow(dfDat), size = nrow(dfDat)*0.20, replace = F)
# warning: on slow machines this can take 10 to 20 minutes to finish, reduce boot number in that case
oVar = CVariableSelection.ReduceModel(dfDat[-test,], fSamples[-test], 100)

pdf('Paper/Results/crossvalidation_variable_selection.pdf')
plot.var.selection(oVar)
table(fSamples[test])

## 10 fold nested cross validation with various variable combinations

par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:5){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar, i)
  dfData.train = as.data.frame(dfDat[-test, cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  dfData.test = data.frame(dfDat[test, cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fSamples[test],
                             train.groups = fSamples[-test], level.predict = 'Disease', boot.num = 100)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

dev.off(dev.cur())
sapply(1:5, function(x) {cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar, x)
})


