# Name: generateTrainingData.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 23/04/2018
# Desc: generate training data

### libraries to load
library(GEOquery)
library(Biobase)
library(lumi)
library(lumiHumanAll.db)
library(lumiHumanIDMapping)
library(annotate)
library(limma)
library(downloader)

setwd('workflow/')
dir.create('dataExternal')
## download the data set
## go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37250

# url = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE37nnn/GSE37250/matrix/GSE37250_series_matrix.txt.gz'
# download(url, destfile='dataExternal/GSE37250_series_matrix.txt.gz')
oExp = getGEO(filename = file.choose())
## ignore the warning as it is a bug
# add lumi nuIDs - converting probe ids to identify genes
oExp = addNuID2lumi(oExp, lib.mapping = 'lumiHumanIDMapping' )

## check the metadata
df = pData(oExp)
str(df)
head(df)

## the main grouping of interest and control variable
xtabs(~ df$characteristics_ch1 + df$characteristics_ch1.1)

## create short names for these factors
fSamples = rep(NA, times=nrow(df))
f = as.character(df$characteristics_ch1)
# create factors
i = grep('active tuberculosis', f)
fSamples[i] = 'ATB'

i = grep('latent TB', f)
fSamples[i] = 'LTB'

i = grep('other', f)
fSamples[i] = 'HC'

fSamples = factor(fSamples, levels = c('HC', 'LTB', 'ATB'))
table(fSamples, f)
table(fSamples)
levels(fSamples)

df$fSamples = fSamples

## second factor for hiv status
fSamples = rep(NA, times=nrow(df))
f = as.character(df$characteristics_ch1.1)

# create factors
i = grep('HIV positive', f)
fSamples[i] = 'HIV+'

i = grep('HIV negative', f)
fSamples[i] = 'HIV-'

fSamples = factor(fSamples, levels = c('HIV-', 'HIV+'))
table(fSamples, f)
table(fSamples)
levels(fSamples)

df$fHIV = fSamples

## update the table
pData(oExp) = df

#### examine the data matrix
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

mData = exprs(oExp)
range(mData)
dim(mData)
i = sample(1:nrow(mData), 2000, replace = F)
mData = mData[i,]
oDiag.1 = CDiagnosticPlots(mData, 'Raw Matrix')
fBatch = df$fSamples
## check data quality
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)

#################################### data normalization
# remove negative values first and set minimum value to 1
exprs(oExp) = exprs(oExp) + abs(min(exprs(oExp))) + 1
range(exprs(oExp))

mData = log(exprs(oExp))
range(mData)
dim(mData)
i = sample(1:nrow(mData), 2000, replace = F)
mData = mData[i,]
## data has negative values, make positive before further analysis
oDiag.1 = CDiagnosticPlots(mData, 'Log Raw Matrix')
fBatch = df$fSamples

## check normalisation
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F
## this should give an error if scaling can't be done
## if all the vector 0 for PCA
oDiag.1.2 = CDiagnosticPlotsSetParameters(oDiag.1, l)
plot.PCA(oDiag.1.2, fBatch, legend.pos = 'topright')
plot.dendogram(oDiag.1.2, fBatch, labels_cex = 0.8, cex.main=0.7)

# normalize and log2 transform the data using lumi
oExp.lumi = lumiT(oExp, 'log2')
## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

## check data quality after normalisation
mData = exprs(oExp)
range(mData)
dim(mData)
i = sample(1:nrow(mData), 2000, replace = F)
mData = mData[i,]
## data has negative values, make positive before further analysis
oDiag.2 = CDiagnosticPlots(mData, 'Normalised Matrix')
fBatch = df$fSamples

## check normalisation
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)
plot.PCA(oDiag.2, fBatch, cex.main=1)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)
## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.2)
l$PCA.jitter = F
l$HC.jitter = F
## this should give an error if scaling can't be done
## if all the vector 0 for PCA
oDiag.2.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.2.2, fBatch, legend.pos = 'topright')
plot.dendogram(oDiag.2.2, fBatch, labels_cex = 0.8, cex.main=0.7)
rm(mData)

##################################### perform DE analysis
mDat = exprs(oExp)
df = select(lumiHumanAll.db, keys = rownames(mDat), columns=c('ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'PROBEID')
## drop columns with NA
df = na.omit(df)
## drop columns with duplicate genes
i = duplicated(df$SYMBOL)
table(i)
df = df[!i, ]
head(df)
## put both tables in same order
i = match(df$PROBEID, rownames(mDat))
mDat = mDat[i,]
# sanity check
identical(rownames(mDat), df$PROBEID)
rownames(mDat) = df$ENTREZID

##### perform DE analysis
## first merge the groups LTB and HC 
fSamples = as.character(oExp$fSamples)
fSamples[fSamples != 'ATB'] = 'Other'
fSamples = factor(fSamples, levels = c('Other', 'ATB'))
table(fSamples, oExp$fSamples)

design = model.matrix(~ fSamples + oExp$fHIV)
head(design)

fit = lmFit(mDat, design)
fit = eBayes(fit)

dfLimmma = topTable(fit, coef=2, adjust='BH', number = Inf)
head(dfLimmma)

hist(dfLimmma$logFC)
hist(dfLimmma$adj.P.Val)
table(dfLimmma$adj.P.Val < 0.01)

dfLimmma = dfLimmma[order(dfLimmma$adj.P.Val, decreasing = F), ]

### save for further analysis
oExp.save = ExpressionSet(mDat)
pData(oExp.save) = pData(oExp)

## save this data 
lData.train = list(data=oExp.save, results=dfLimmma)
dir.create('results')
save(lData.train, file='results/lData.train.rds')
