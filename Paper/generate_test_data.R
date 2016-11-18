# File: generate_test_data.R
# Auth: Umar Niazi umar.niazi@kcl.ac.uk
# Date: 18/11/2016
# Desc: generate 2 TB datasets to analyse

library(GEOquery)
library(Biobase)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(annotate)
library(limma)
library(sva)

## internal function
# Function: f_lGetPCAClusterCount
# Desc: Takes first two components of the PCA and counts possible clusters. The function does this by 
#       binning the vector of data into X bins, assigning a class label to each bin, counting how many
#       observations in each bin, any total number of bins with at least one observations. This is calculated
#       for both the components of the pca matrix, and the max number of bins with at least one observation, in
#       first or second dimension is reported along with a data.frame with cluster labels.
# Args: pr.out = principal component object returned by prcomp function
# Rets: returns list with 2 elements: 
#       1 - cluster.count = possible number of clusters in the data
#       2 - cluster.label = data.frame with cluster labels
f_lGetPCAClusterCount = function(pr.out){
  # how many clusters in data, using first 2 components
  x1 = pr.out$x[,1]
  x2 = pr.out$x[,2]
  # bin the data from the 2 components
  h1 = hist(x1, plot=F)
  # give a class label to each bin
  c1 = cut(x1, h1$breaks, labels = 1:(length(h1$mids)))
  h2 = hist(x2, plot=F)
  c2 = cut(x2, h2$breaks, labels = 1:(length(h2$mids)))
  # labels for vectors and the class labels
  dfClust = data.frame(lab=names(x1), c1, c2)
  # get contingency table
  mClust = as.matrix(table(c1 = dfClust$c1, c2 = dfClust$c2))
  # count the max of row and col sums that are not zero
  ir = length(which(rowSums(mClust) != 0))
  ic = length(which(colSums(mClust) != 0))
  iClust.count = ifelse(ir > ic, ir, ic)
  lRet = list(cluster.count=iClust.count, cluster.label=dfClust)
  return(lRet)
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

### end internal functions

# global variables
p.old = par()
setwd('Paper/')
###############################################################################
#### dataset 1 TB longitudinal
###############################################################################
## data loading
# load the data, clean and create factors
dir.create('Data_external', showWarnings = F)
gse =  getGEO('GSE19491', GSEMatrix = T, destdir = 'Data_external/')
oExp = gse$GSE19491_series_matrix.txt.gz

# print samples
as.data.frame(table(oExp$source_name_ch1))

# get the whole blood data
i = grep('treatment', x = oExp$source_name_ch1, ignore.case = F, perl = T)
oExp = oExp[,i]

## data normalization
# normalize and log2 transform the data using lumi
oExp.lumi = lumiT(oExp, 'log2')
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))
fSamples = rep(NA, 21)
f = as.character(oExp$source_name_ch1)
# create factors
i = grep('before treatment', f)
fSamples[i] = '0'

i = grep('2 months', f)
fSamples[i] = '2'

i = grep('12 months', f)
fSamples[i] = '12'

fSamples = factor(fSamples, levels = c('12', '2', '0'))
levels(fSamples)
oExp.lumi$fSamples = fSamples

## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

# add lumi nuIDs 
oExp = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )
# remove any NA data
exprs(oExp) = na.omit(exprs(oExp))
fSamples = oExp$fSamples
table(fSamples)
levels(fSamples)

### perform DE analysis
mDat = exprs(oExp)
# map the nuID to human symbols
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]

design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get annotation
df = select(lumiHumanAll.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# sanity check
nrow(mDat) == nrow(df)
# add annotation to limma object
fit$genes = df
topTable(fit, adjust='BH')

# look at top tables for each comparison
for (i in 2:length(levels(fSamples))){
  print(paste(levels(fSamples)[i], i))
  print(topTable(fit, coef=i, adjust='BH'))
}

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.adj = p.adjust(fit$p.value[,i], method = 'BH')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.01]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

######### Volcano plots
# plot volcano plots
par(p.old)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  f_plotVolcano(dfGenes, paste(names(n[i])), fc.lim = c(-3, 3))
}

# get the common genes
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)

## export the result for cluster analysis
dfData = t(mDat[cvCommonGenes,])
cn = select(lumiHumanAll.db, keys = colnames(dfData), columns = c('ENTREZID'), keytype = 'PROBEID')
colnames(dfData) = cn$ENTREZID
# remove duplicate probes
f = !duplicated(colnames(dfData))
dfData = dfData[,f]
dfData = data.frame(dfData)
dfData$fSamples = fSamples


dir.create('Test_data', showWarnings = F)

write.csv(dfData, file='Test_data/test_data_GSE19491.csv')



#################################################################################
## dataset 2, atb, ptb and hc
################################################################################
## data loading
# load the data, clean and create factors
dir.create('Data_external', showWarnings = F)
gse =  getGEO('GSE19491', GSEMatrix = T, destdir = 'Data_external/')
oExp = gse$GSE19491_series_matrix.txt.gz

# add lumi nuIDs 
oExp = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )
# remove any NA data
i = which(is.na(rowSums(exprs(oExp))))
oExp = oExp[-i,]

# print samples
as.data.frame(table(oExp$source_name_ch1))
as.data.frame(table(oExp$characteristics_ch1.3))

# get the whole blood data
i = grep('^Whole blood from healthy control$|^Whole Blood from patient with Active TB$|^Whole Blood from patient with Latent TB$', 
         x = oExp$source_name_ch1, ignore.case = T, perl = T)
oExp = oExp[,i]

# sanity check
as.data.frame(table(oExp$source_name_ch1))
as.data.frame(table(oExp$characteristics_ch1.3))

## data normalization
# normalize and log2 transform the data using lumi
# remove negative values first and set minimum value to 1
exprs(oExp) = exprs(oExp) + abs(min(exprs(oExp))) + 1
oExp.lumi = lumiT(oExp, 'log2')
fSamples = rep(NA, times=nrow(pData(oExp.lumi)))
f = as.character(oExp$source_name_ch1)
# create factors
i = grep('healthy control', f)
fSamples[i] = 'HC'

i = grep('Active TB', f)
fSamples[i] = 'ATB'

i = grep('Latent TB', f)
fSamples[i] = 'LTB'

fSamples = factor(fSamples, levels = c('HC', 'LTB', 'ATB'))
table(fSamples)
levels(fSamples)

oExp.lumi$fSamples = fSamples

## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

## check for outliers
# check quality
m = exprs(oExp)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = oExp$fSamples

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)
i = which(l$cluster.label$c1 %in% as.character(6:11))

# sanity check for the outlier
c = col
c[i] = 'black'
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=c, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=c, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)

# create a factor for this outlier
l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)
i = (l$cluster.label$c1 %in% as.character(6:11))
fBatch = factor(as.numeric(i)) 

# check what are these batches
temp = pData(oExp)
xtabs(~ (temp$characteristics_ch1.4) + fBatch + fSamples)

i = grep('London|South Africa', oExp$characteristics_ch1.4)
oExp = oExp[,i]
temp = pData(oExp)
xtabs(~ (temp$characteristics_ch1.4) + temp$fSamples)

fBatch.1 = as.character(oExp$characteristics_ch1.4)
fBatch.1 = gsub('[\\w ]*:[ ]?(\\w+)', replacement = '\\1', fBatch.1, perl = T)

fBatch.2 = as.character(oExp$title)
fBatch.2 = (gsub('^[A-Z]{2,}_(\\w+)', '\\1', fBatch.2, perl = T))
fBatch.2 = gsub('^(LON|SA)[-_]([A-Za-z]+)\\d+', '\\1.\\2', fBatch.2, perl=T)

fBatch = paste0(fBatch.1, fBatch.2)
oExp$fBatch = factor(fBatch)

## add the combat step for batch covariate
modcombat = model.matrix(~1, data=pData(oExp))
oCexp = ComBat(exprs(oExp), batch = oExp$fBatch, mod=modcombat)
exprs(oExp) = oCexp

# # remove the outlier
# oExp = oExp[,-i]

## test quality again
## check for outliers
m = exprs(oExp)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = oExp$fSamples
fSamples = oExp$fBatch
col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')

l = f_lGetPCAClusterCount(pr.out)
l$cluster.count
table(c1 = l$cluster.label$c1, c2 = l$cluster.label$c2)
i = which(l$cluster.label$c1 %in% as.character(6:11))

# sanity check for the outlier
c = col
c[i] = 'black'
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=c, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=c, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=c, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')
par(p.old)

# remove the outlier
oExp = oExp[,-i]

m = exprs(oExp)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = oExp$fSamples

col.p = rainbow(length(unique(fSamples)))
col = col.p[as.numeric(fSamples)]
# plot the pca components
par(mfrow=c(2,2))
plot.new()
legend('center', legend = unique(fSamples), fill=col.p[as.numeric(unique(fSamples))])
plot(pr.out$x[,1:2], col=col, pch=19, xlab='Z1', ylab='Z2',
     main='PCA comp 1 and 2')
plot(pr.out$x[,c(1,3)], col=col, pch=19, xlab='Z1', ylab='Z3',
     main='PCA comp 1 and 3')
plot(pr.out$x[,c(2,3)], col=col, pch=19, xlab='Z2', ylab='Z3',
     main='PCA comp 2 and 3')


### perform DE analysis
mDat = exprs(oExp)
# map the nuID to human symbols
cvSym = getSYMBOL(rownames(mDat), 'lumiHumanAll.db')
# number of genes not annotated
table(is.na(cvSym))
# remove unannotated genes
mDat = mDat[!is.na(cvSym),]
fSamples = oExp$fSamples
design = model.matrix(~fSamples)
colnames(design) = levels(fSamples)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

# get annotation
df = select(lumiHumanAll.db, keys = rownames(mDat), columns = c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
# sanity check
nrow(mDat) == nrow(df)
# add annotation to limma object
fit$genes = df
topTable(fit, adjust='BH')

# look at top tables for each comparison
for (i in 2:length(levels(fSamples))){
  print(paste(levels(fSamples)[i], i))
  print(topTable(fit, coef=i, adjust='BH'))
}

# get the list of genes for each comparison i.e. each coefficient compared to base line
lSigGenes.adj = vector('list', length = length(levels(fSamples))-1)
names(lSigGenes.adj) = levels(fSamples)[2:length(levels(fSamples))]

for (i in 2:length(levels(fSamples))){
  p.adj = p.adjust(fit$p.value[,i], method = 'BH')
  lSigGenes.adj[[i-1]] = names(p.adj)[p.adj < 0.01]
}

#cvSigGenes.adj = unique(cvSigGenes.adj)
sapply(lSigGenes.adj, length)

######### Volcano plots
# plot volcano plots
par(p.old)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1

for (i in seq_along(n)) {
  dfGenes = topTable(fit, coef = n[i], number = Inf)
  f_plotVolcano(dfGenes, paste(names(n[i])), fc.lim = c(-2, 2))
}

# get the common genes
dfRes = topTable(fit, adjust='BH', number=Inf, p.value=0.1)
n = (which(sapply(lSigGenes.adj, length) >= 10)) + 1
cvCommonGenes = NULL
for (i in seq_along(n)) {
  cvCommonGenes = append(cvCommonGenes, lSigGenes.adj[[names(n[i])]])
}
cvCommonGenes = unique(cvCommonGenes)

## export the result for cluster analysis
dfData = t(mDat[cvCommonGenes,])
cn = select(lumiHumanAll.db, keys = colnames(dfData), columns = c('ENTREZID'), keytype = 'PROBEID')
colnames(dfData) = cn$ENTREZID
# remove duplicate probes
f = !duplicated(colnames(dfData))
dfData = dfData[,f]
dfData = data.frame(dfData)
dfData$fSamples = fSamples

dir.create('Test_data', showWarnings = F)

write.csv(dfData, file='Test_data/test_data_GSE19491_atb_ltb.csv')

