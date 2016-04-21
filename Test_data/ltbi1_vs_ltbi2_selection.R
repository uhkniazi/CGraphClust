# Name: ltbi1_vs_ltbi2_selection.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 05/04/2016
# Desc: script for variable selection

if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')


# utility function to load object
f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}

f_dfGetGeneAnnotation = function(cvEnterezID = NULL) {
  if (!require(org.Hs.eg.db)) stop('org.Hs.eg.db annotation library required')
  return(AnnotationDbi::select(org.Hs.eg.db, cvEnterezID, columns = c('SYMBOL', 'GENENAME'), keytype = 'ENTREZID'))  
}


## load the test data
lData = f_LoadObject('Objects/ltb2_atb_data.rds')
# it is a list with various components
names(lData)
# setup data for variable selection
dfDat = data.frame(lData$matrix)
fSamples = lData$groups
i = grep('LTB', fSamples)
dfDat = dfDat[i,]
rownames(dfDat) = names(fSamples[i])
fSamples = as.character(fSamples[i])
fSamples = factor(fSamples, levels = c('LTB2', 'LTB1'))


oRan = CVariableSelection.RandomForest(dfDat, fSamples, 100, big.warn = F)
plot.var.selection(oRan)

par(p.old)
dfRF = CVariableSelection.RandomForest.getVariables(oRan)
hist(dfRF$ivMean)

# get the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]
f_dfGetGeneAnnotation(gsub('^X(\\d+)', '\\1', cvTopGenes))
## variable combinations
dfDat.top = dfDat[,cvTopGenes]

oVar = CVariableSelection.ReduceModel(dfDat.top, fSamples, 100)
plot.var.selection(oVar)


## get the combination that produces the best result
set.seed(1234)
test = sample(1:nrow(dfDat), size = nrow(dfDat) * 0.2)
table(fSamples[test])

## 10 fold nested cross validation with various variable combinations
#par(mfrow=c(2,2))
# try models of various sizes with CV
for (i in 1:5){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar, i)
  dfData.train = as.data.frame(dfDat[-test, cvTopGenes.sub])
  colnames(dfData.train) = cvTopGenes.sub
  dfData.test = data.frame(dfDat[test, cvTopGenes.sub])
  colnames(dfData.test) = cvTopGenes.sub
  
  oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fSamples[test],
                             train.groups = fSamples[-test], level.predict = 'LTB1', boot.num = 100)
  
  plot.cv.performance(oCV)
  # print variable names and 95% confidence interval for AUC
  temp = oCV@oAuc.cv
  x = as.numeric(temp@y.values)
  print(paste('Variable Count', i))
  print(cvTopGenes.sub)
  print(signif(quantile(x, probs = c(0.025, 0.975)), 2))
}

sapply(1:5, function(x) {cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar, x)
    print(f_dfGetGeneAnnotation(gsub('^X(\\d+)', '\\1', cvTopGenes.sub)))})


## make plots for the top genes
cvTopGenes.sub = gsub('^X(\\d+)', '\\1', cvTopGenes)
mCounts = lData$matrix
mCounts = mCounts[,cvTopGenes.sub]
df = f_dfGetGeneAnnotation(cvTopGenes.sub)
colnames(mCounts) = df$SYMBOL
fSamples = lData$groups

par(mfrow=c(1,2))
sapply(seq_along(1:ncol(mCounts)), function(x) boxplot(mCounts[,x] ~ fSamples, main=colnames(mCounts)[x]))

library(NMF)
m1 = mCounts
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = colnames(mCounts)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)



