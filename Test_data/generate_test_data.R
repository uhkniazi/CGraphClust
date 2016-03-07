# generate_test_data.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 23/04/2015
# Desc: analysis of dataset to generate test data

library(GEOquery)
library(Biobase)
library(lumi)

# global variables
p.old = par()

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

fSamples = factor(fSamples)

oExp.lumi$fSamples = fSamples

## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

## extract data for significant genes
fSamples = oExp$fSamples
pano = apply(exprs(oExp), 1, function(x) anova(lm(x~fSamples))$Pr[1])
pano.adj = p.adjust(pano, method = 'BH')
n = which(pano.adj < 0.1)
n = names(n)
# count matrix
dfData = data.frame(t(exprs(oExp)[n,]))
n = colnames(dfData)
dfAnnotation = fData(oExp)
# replace names by enterez ids
dfAnnotation = dfAnnotation[n,c('ID', 'Entrez_Gene_ID')]
# remove empty enterez ids and NAs
i = as.numeric(dfAnnotation$Entrez_Gene_ID)
i = which(is.na(i))
dfAnnotation = dfAnnotation[-i,]
# remove duplicated enterez ids
i = duplicated(dfAnnotation$Entrez_Gene_ID)
dfAnnotation = dfAnnotation[!i,]
# get the names of these annotation ids
n = as.character(dfAnnotation$ID)
dfData = dfData[,n]
# replace names by enterez ids
n = as.character(dfAnnotation$Entrez_Gene_ID)
colnames(dfData) = n
# assign sample ids
fSamples = oExp$fSamples
dfData$fSamples = fSamples
dir.create('Test_data', showWarnings = F)

write.csv(dfData, file='Test_data/test_data_GSE19491.csv')


############################################################################
## generate second dataset
############################################################################

library(GEOquery)
library(Biobase)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(annotate)
library(limma)

# global variables
p.old = par()

## data loading
# load the data, clean and create factors
dir.create('Data_external', showWarnings = F)
gse =  getGEO('GSE54514', GSEMatrix = T, destdir = 'Data_external/', AnnotGPL = T, getGPL = T)
oExp = gse$GSE54514_series_matrix.txt.gz

# add lumi nuIDs 
oExp.lumi = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))

# get the grouping factor
fSamples = as.character(pData(oExp.lumi)$source_name_ch1)
# comparing survivors over days
# remove control and nonsurvivors
i = grep('Control|sepsis_nonsurvivor', fSamples)
oExp.lumi = oExp.lumi[,-i]
dfSamples = pData(oExp.lumi)

fSamples = as.character(dfSamples$source_name_ch1)
# get the covariates
fDays = gsub('[\\w+ ,]+Day (\\d+).+', '\\1', x = fSamples, perl = T)
fGender = as.character(dfSamples$characteristics_ch1.4)
fGender = gsub('\\gender: (\\w)', '\\1', fGender)
fAge = as.character(dfSamples$characteristics_ch1.5)
iAge = as.numeric(gsub('age \\(years\\): (\\d+)', '\\1', fAge))
fAge = cut(iAge, quantile(iAge, 0:5/5), include.lowest = T, labels = paste0('Age', 1:5))

fBatch = paste0(fGender, fAge)
oExp.lumi$fDays = fDays
oExp.lumi$fBatch = fBatch
# normalize the data
lumi.n = lumiN(oExp.lumi, method = 'rsn')

# check quality
m = exprs(lumi.n)
# pca on samples i.e. covariance matrix of m
pr.out = prcomp(t(m), scale=T)
## choose appropriate factor
fSamples = as.factor(lumi.n$fBatch)
fSamples = as.factor(lumi.n$fDays)

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


## add the combat step for batch covariate
library(sva)
modcombat = model.matrix(~1, data=pData(lumi.n))
oCexp = ComBat(exprs(lumi.n), batch = lumi.n$fBatch, mod=modcombat)
exprs(lumi.n) = oCexp

oExp.lumi = lumi.n
## select grouping and data for DE analysis
dfSamples = pData(oExp.lumi)

## create grouping factors by data
fSamples = factor(paste0('D', dfSamples$fDays))
#fSamples = factor(dfSamples$fBatch)
table(fSamples)

### perform DE analysis
mDat = exprs(oExp.lumi)
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


###### internal function
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
##########

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

## save the data
dir.create('Test_data', showWarnings = F)

write.csv(dfData, file='Test_data/test_data_GSE54514.csv')

