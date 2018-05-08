# Name: biomarkers.R
# Auth: uhkniazi
# Date: 02/05/2018
# Desc: variable selection for biomarker discovery in TB dataset

library(org.Hs.eg.db)
library(downloader)
source('CGraphClust.R')

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

## combine the groups LTB and HC
fGroups = rep('ATB', times=length(oExp$fSamples))
fGroups[oExp$fSamples != 'ATB'] = 'Other'
fGroups = factor(fGroups, levels = c('Other', 'ATB'))
table(fGroups)
table(oExp$fSamples)

## assign annotation names
mCounts = t(exprs(oExp))
dim(mCounts)

## sub select top genes sorted on p-values
cvTopGenes = rownames(lData.train$results)[1:2000]
mCounts = mCounts[,cvTopGenes]

# convert enterez ids to uniprot as Reactome database file uses UNIPROT ids
dfMap = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'UNIPROT', 'ENTREZID')
dfMap = na.omit(dfMap)

### load the uniprot2reactome mapping obtained from
## any other databsae of choice can be used as long as there are 2 columns
## where Column 1 = Gene and Column 2 = Database ID (label)

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
head(dfReactome)

## map reactome ids to uniprot ids
## where dfMap is the dataframe created earlier with gene ids mapped to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

rm(dfReactome)
gc(reset = T)

# select genes that have a reactome term attached
# this will reduce your initial gene list as a lot of genes actually have 
# no annotations, you may use a low level database like GO if you think 
# that you lose too many genes
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
# order the count matrix based on grouping factor
# this will serve useful later for plotting but is not necessary
# rownames(mCounts) = fGroups
# mCounts = mCounts[order(fGroups),]
# fGroups = fGroups[order(fGroups)]

head(fGroups)

# create a correlation matrix to decide cor cutoff
# in a random collection of genes the correlations should follow roughly a 
# normal distribution centered at 0, as this is a subset of differentially expressed
# genes, they should show a bimodal pattern
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# create the graph cluster object
# using absolute correlation to cluster positively and negatively correlated genes
# the cutoff of 0.6 and absolute is your choice, we would base this cutoff after looking
# at the histogram, if you move this too low then you may end up with a graph with too many
# connections and it may be dominated by noise/random associations, and it will likely crash you machine
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.4, bSuppressPlots = F)

# you may not want to suppress the 2 plots, to see the distribution fits
# 1- degree distribution of type 2 vertices, i.e. connections of IDs and Genes
# bipartite graph is created and some cleaning performed, i.e. 
# those type 2 vertices with too many connections (degree) with type 1 vertices are removed.
# this is done by looking at the distribution of the degrees on a log scale and approximating
# a negative binomial and poisson distributions on top. the cutoff is by default 0.95 quantile under
# a negative binomial model, anything over that is removed. The idea is that type vertices that have a
# lot of connections, are not very interesting and they tend to hide the more interesting connections.

# 2- distribution of observed to expected ratio weights for connections between genes
# graph is projected on to one dimension (type 1 vertices) and connections between type 1
# vertices are assigned interestingness or observed to expected probability ratios as weights.
# the weights on a square root scale follow a negative binomial distribution and only the highest
# weights are chosen, as they are the most interesting. So anything over 0.95 quantile under a neg bin model
# are chosen and the other connections are discarded.

## general graph structure
## we would like to see how does the graph look like, are the clusters connected or in subgraphs
set.seed(1)
plot.final.graph(oGr)
ecount(getFinalGraph(oGr))
vcount(getFinalGraph(oGr))

## community structure
## overview of how the commuinties look like
# plot the main communities in 2 different ways
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 20)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 30)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')


## if there are small sized communities/clusters then you may want to remove them?
## we normally remove clusters with 5 or less genes
# get community sizes
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
iSizes
# remove communities smaller than 5 members
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}

# make a pdf output for publication
dir.create('Paper/Results', showWarnings = F)
pdf('Paper/Results/community_tb.pdf')
par(mar=c(1,1,1,1)+0.1, family='Helvetica')
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
set.seed(1)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr,
     vertex.frame.color=NA, edge.color='darkgrey')
set.seed(1)
ig = plot.centrality.graph(oGr)
dev.off(dev.cur())

## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
## centrality in social networks is simpler to understand intuitively 
## e.g. degree = people with lots of followers 
## hub = people with lots of followers connected to people with lots of followers
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 20)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
# these plots are bar plots of ordered statistics - sorted centrality parameter
# usually the last 10% (on the right) will show a sudden rise and may be the interesting genes/nodes
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

## try the genes with various centrality scores
colnames(mCent)
mCent = mCent[order(mCent[,'closeness'], decreasing = T), ]
head(mCent)
cvTopGenes = rownames(mCent)[1:30]
dfTopGenes.cent = f_dfGetGeneAnnotation(cvTopGenes)

## top vertices/nodes/genes based on centrality scores
## get a table of top vertices 
# dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.80)
# rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# # assign metadata annotation to these genes and clusters
# dfCluster = getClusterMapping(oGr)
# colnames(dfCluster) = c('gene', 'cluster')
# rownames(dfCluster) = dfCluster$gene
# df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
# dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
# dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
# dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

####### NOTE: This section of code is very slow, use only if you need data from genbank
## not run
# # loop and get the data from genbank
# n = rep(NA, length=nrow(dfTopGenes.cent))
# names(n) = as.character(dfTopGenes.cent$VertexID)
# for (i in seq_along(n)){
#   n[i] = f_csGetGeneSummaryFromGenbank(names(n)[i])
#   # this wait time is required to stop sending queries to ncbi server very quickly
#   Sys.sleep(time = 3)
# }
# cvSum.2 = as.character(dfTopGenes.cent$VertexID)
# dfTopGenes.cent$Summary = n[cvSum.2]
####### Section ends

# dir.create('Paper/Results', showWarnings = F)
# write.csv(dfTopGenes.cent, file='Paper/Results/Top_Centrality_Genes_tb.xls')

## perform variable selection
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## subset the data matrix
## assign annotation names
mCounts = t(exprs(oExp))
dim(mCounts)
mCounts = mCounts[, dfTopGenes.cent$ENTREZID]
colnames(mCounts) = as.character(dfTopGenes.cent$SYMBOL)
dim(mCounts)

## prepare data for modelling in stan
dfData = data.frame(scale(mCounts))
dim(dfData)
dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == 'ATB', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

detach("package:org.Hs.eg.db", unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='workflow/binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

save(fit.stan, file='workflow/results/fit.stan.binom_reverse.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]

## function to calculate statistics for a coefficient
getDifference = function(ivData){
  # get the difference vector
  d = ivData
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(p)
}

ivPval = apply(mCoef, 2, getDifference)
hist(ivPval)
plot(colMeans(mCoef), ivPval, pch=19)
m = colMeans(mCoef)
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
m = abs(m)
m = sort(m, decreasing = T)
cvTopGenes.binomial = names(m)[1:30] #names(m[m > 0.25])

par(mar=c(6,3,4,2)+0.1)
l2 = barplot(m[cvTopGenes.binomial], 
             las=2, xaxt='n', col='grey', main='Top Variables')
axis(1, at = l2, labels = cvTopGenes.binomial, tick = F, las=2, cex.axis=0.7 )

# use the top 30 genes to find top combinations of genes
dfData = data.frame(mCounts)#[,cvTopGenes.binomial])
head(dfData)
oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 50)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

cvGenes.hub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 4)
cvGenes.degree = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 1)
cvGenes.between = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 4)
cvGenes.closeness = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 4)
cvTopGenes = unique(c(cvGenes.hub, cvGenes.degree, cvGenes.closeness, cvGenes.between))
## choose the 1 variable model
# dfData.train = data.frame(mCounts)
# dfData = data.frame(dfData.train[,CVariableSelection.ReduceModel.getMinModel(oVar.sub, 2)])
# colnames(dfData) = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 2)
mCounts = t(exprs(oExp))
dim(mCounts)
df = AnnotationDbi::select(org.Hs.eg.db, keys = cvTopGenes, keytype = 'SYMBOL', columns = 'ENTREZID')
mCounts = mCounts[, df$ENTREZID]
colnames(mCounts) = as.character(df$SYMBOL)
dim(mCounts)
dfData = data.frame(mCounts)

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 50)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

dfData = data.frame(dfData[,CVariableSelection.ReduceModel.getMinModel(oVar.sub, 7)])
head(dfData)
dfData$fGroups = fGroups
### fit a model to this data
fit.1 = lda(fGroups ~ ., data=dfData)

## check the prediction on this data
ivPredict = predict(fit.1, newdata = dfData)$posterior[, 'ATB']

## how good is this predictor
library(lattice)
densityplot(~ logit(ivPredict), groups=fGroups, type='n', xlab='Predicted Score', main='Model predicted score', auto.key = list(columns=2))
xyplot(logit(ivPredict) ~ lData.train$data$fSamples, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)',
       main='Predicted scores vs Actual groups')

### plot the score performance
cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, 7)
dfData.train = data.frame(mCounts)
dfData.train = data.frame(dfData.train[,cvTopGenes.sub])
colnames(dfData.train) = cvTopGenes.sub
head(dfData.train)
oCV = CCrossValidation.LDA(test.dat = dfData.train, train.dat = dfData.train, test.groups = fGroups,
                           train.groups = fGroups, level.predict = 'ATB', boot.num = 100)

plot.cv.performance(oCV)

df = getCutoffTprFprCrossValidation(oCV)

fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.6] = 'ATB'
table(fPredict, fGroups)







## perform nested random forest
## adjust boot.num as desired
oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100, big.warn = F)
save(oVar.r, file='oVar.r_reverse.rds')

plot.var.selection(oVar.r)

############## use a binomial regression approach this time to rank variables
dfData = data.frame(scale(t(lData.train$data[cvTopVariables, ])))
dim(dfData)
dfData$fGroups = fGroups

lData = list(resp=ifelse(dfData$fGroups == 'ATB', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 11))

save(fit.stan, file='fit.stan.binom_reverse.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
# ## get the intercept at population level
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]

## function to calculate statistics for a coefficient
getDifference = function(ivData){
  # get the difference vector
  d = ivData
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(p)
}

ivPval = apply(mCoef, 2, getDifference)
hist(ivPval)
plot(colMeans(mCoef), ivPval, pch=19)
m = colMeans(mCoef)
names(m) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
m = abs(m)
m = sort(m, decreasing = T)
cvTopGenes.binomial = names(m)[1:30] #names(m[m > 0.25])

p.old = par(mar=c(6,3,4,2)+0.1)
y = c(0, 2.5)
l2 = barplot(m[cvTopGenes.binomial], 
             las=2, xaxt='n', col='grey', main='Top Variables', ylim=y)
axis(1, at = l2, labels = cvTopGenes.binomial, tick = F, las=2, cex.axis=0.7 )


#cvTopGenes.binomial[11] = "HLA-DPA1"

### test performance of both results, from Random Forest and Binomial Regression
dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# select the top 30 variables
cvTopGenes = rownames(dfRF)[1:30]

# use the top 30 genes to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes, ]))

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

# use the top 30 genes to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes.binomial, ]))

oVar.sub2 = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub2)

# print variable combinations
for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

for (i in 1:7){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

### try a combination of top genes from both models
cvTopGenes.comb = NULL;
for (i in 1:15){
  cvTopGenes.comb = append(cvTopGenes.comb, CVariableSelection.ReduceModel.getMinModel(oVar.sub, i)) 
  cat(i)
}
cvTopGenes.comb = unique(cvTopGenes.comb)

for (i in 1:15){
  cvTopGenes.comb = append(cvTopGenes.comb, CVariableSelection.ReduceModel.getMinModel(oVar.sub2, i))
  cat(i)
}

cvTopGenes.comb = unique(cvTopGenes.comb)
length(cvTopGenes.comb)

# use these combined variables to find top combinations of genes
dfData = data.frame(t(lData.train$data[cvTopGenes.comb, ]))

oVar.subComb = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)

# plot the number of variables vs average error rate
plot.var.selection(oVar.subComb)

for (i in 1:10){
  cvTopGenes.sub = CVariableSelection.ReduceModel.getMinModel(oVar.subComb, i)
  cat('Variable Count', i, paste(cvTopGenes.sub), '\n')
  #print(cvTopGenes.sub)
}

################################################# binomial regression with mixture model section
################ fit a binomial model on the chosen model size based on previous results
## this can be another classifier as well e.g. LDA. Using this model check how is the performance 
## and using this make some calibration curves to select decision boundary

library(LearnBayes)
logit.inv = function(p) {exp(p)/(exp(p)+1) }
dfData = data.frame(t(lData.train$data[cvTopGenes.comb, ]))

## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  return(iFitted)
}

## lets write a custom glm using a bayesian approach
## write the log posterior function
mylogpost = function(theta, data){
  ## parameters to track/estimate
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  resp = data$resp # resp
  mModMatrix = data$mModMatrix
  
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  iFitted = logit.inv(iFitted)
  # write the priors and likelihood
  lp = dnorm(betas[1], 0, 10, log=T) + sum(dnorm(betas[-1], 0, 10, log=T))
  lik = sum(dbinom(resp, 1, iFitted, log=T))
  val = lik + lp
  return(val)
}

dfData = data.frame(dfData[ , CVariableSelection.ReduceModel.getMinModel(oVar.subComb, 7)])
colnames(dfData) = CVariableSelection.ReduceModel.getMinModel(oVar.subComb, 7)
dim(dfData)
head(dfData)
dfData = data.frame(dfData, fGroups=fGroups)

lData = list(resp=ifelse(dfData$fGroups == 'ATB', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
start = c(rep(0, times=ncol(lData$mModMatrix)))
mylogpost(start, lData)

fit.2 = laplace(mylogpost, start, lData)
fit.2

fit.1 = glm(fGroups ~ ., data=dfData, family='binomial')
data.frame(coef(fit.1), fit.2$mode)
# se = sqrt(diag(fit.2$var))
# summary(fit.1)
# ## lets take a sample from this
# # parameters for the multivariate t density
# tpar = list(m=fit.2$mode, var=fit.2$var*2, df=4)
# ## get a sample directly and using sir (sampling importance resampling with a t proposal density)
# s = sir(mylogpost, tpar, 5000, lData)
# colnames(s) = colnames(lData$mModMatrix)
# apply(s, 2, mean)
# apply(s, 2, sd)
# pairs(s, pch=20)
# fit.2$sir = s

stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=3, pars=c('tau', 'betas2'), init=initf, 
                    control=list(adapt_delta=0.99, max_treedepth = 13))

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', CVariableSelection.ReduceModel.getMinModel(oVar.subComb, 7))
pairs(mCoef, pch=20)


### once we have results from the classifier we can make some plots to see
### the performance
library(lattice)
library(car)
## get the predicted values
dfData.new = dfData
str(dfData.new)
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData.new)), dfData.new[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
ivPredict = mypred(colMeans(mCoef), list(mModMatrix=X))[,1]
xyplot(ivPredict ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)')
xyplot(ivPredict ~ lData.train$grouping, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)',
       main='Predicted scores vs Actual groups')
densityplot(~ ivPredict, data=dfData, type='n')
densityplot(~ ivPredict | fGroups, data=dfData, type='n', xlab='Predicted Score', main='Actual Scale')
densityplot(~ ivPredict, groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Actual Scale', auto.key = list(columns=2))

## lets check on a different scale of the score
densityplot(~ logit(ivPredict), data=dfData)
xyplot(logit(ivPredict) ~ lData.train$grouping, xlab='Actual Group', ylab='Predicted Probability of Being ATB (1)')
densityplot(~ logit(ivPredict), groups=fGroups, data=dfData, type='n', 
            xlab='Predicted Score', main='Logit Scale', auto.key = list(columns=2))

# convert to logit scale for model fitting
ivPredict = logit(ivPredict)
################################ section for mixture model
######## this mixture model will help us decide an appropriate cutoff for the decision rule
######## see Gelman 2013 around P18 for an example of record linking score calibration
stanDso = rstan::stan_model(file='normResponseFiniteMixture_2.stan')

## take a subset of the data
lStanData = list(Ntotal=length(ivPredict), y=ivPredict, iMixtures=2)

## give initial values if you want, look at the density plot 
initf = function(chain_id = 1) {
  list(mu = c(-5, 5), sigma = c(1, 1), iMixWeights=c(0.5, 0.5))
} 

## give initial values function to stan
# l = lapply(1, initf)
fit.stan = sampling(stanDso, data=lStanData, iter=4000, chains=4, cores=4, init=initf)
print(fit.stan, digi=3)
traceplot(fit.stan)

## check if labelling degeneracy has occured
## see here: http://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
params1 = as.data.frame(extract(fit.stan, permuted=FALSE)[,1,])
params2 = as.data.frame(extract(fit.stan, permuted=FALSE)[,2,])
params3 = as.data.frame(extract(fit.stan, permuted=FALSE)[,3,])
params4 = as.data.frame(extract(fit.stan, permuted=FALSE)[,4,])

## check if the means from different chains overlap
## Labeling Degeneracy by Enforcing an Ordering
par(mfrow=c(2,2))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
plot(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)
plot(params3$`mu[1]`, params3$`mu[2]`, pch=20, col=4)
plot(params4$`mu[1]`, params4$`mu[2]`, pch=20, col=5)

par(mfrow=c(1,1))
plot(params1$`mu[1]`, params1$`mu[2]`, pch=20, col=2)
points(params2$`mu[1]`, params2$`mu[2]`, pch=20, col=3)
points(params3$`mu[1]`, params3$`mu[2]`, pch=20, col=4)
points(params4$`mu[1]`, params4$`mu[2]`, pch=20, col=5)


# model checks
############# extract the mcmc sample values from stan
mStan = do.call(cbind, extract(fit.stan))
mStan = mStan[,-(ncol(mStan))]
colnames(mStan) = c('mu1', 'mu2', 'sigma1', 'sigma2', 'mix1', 'mix2')
dim(mStan)
## get a sample for this distribution
########## simulate 200 test quantities
mDraws = matrix(NA, nrow = length(ivPredict), ncol=200)

for (i in 1:200){
  p = sample(1:nrow(mStan), size = 1)
  mix = mean(mStan[,'mix1'])
  ## this will take a sample from a normal mixture distribution
  sam = function() {
    ind = rbinom(1, 1, prob = mix)
    return(ind * rnorm(1, mStan[p, 'mu1'], mStan[p, 'sigma1']) + 
             (1-ind) * rnorm(1, mStan[p, 'mu2'], mStan[p, 'sigma2']))
  }
  mDraws[,i] = replicate(length(ivPredict), sam())
}

mDraws.normMix = mDraws

yresp = density(ivPredict)
yresp$y = yresp$y/max(yresp$y)
plot(yresp, xlab='', main='Fitted distribution', ylab='scaled density', lwd=2)
temp = apply(mDraws, 2, function(x) {x = density(x)
x$y = x$y/max(x$y)
lines(x, col='darkgrey', lwd=0.6)
})
lines(yresp, lwd=2)


print(fit.stan)

range(ivPredict)
## reconvert back to inverse logit scale i.e. 0 to 1 range
ivPredict = plogis(ivPredict)
# temp = rnorm(1000, 0.02, 0.02)
# sapply(seq(0, 1, length.out = 10), function(x) sum(temp >= x)/1000)
# round(pnorm(seq(0, 1, length.out = 10), 0.02, 0.02, lower.tail = F), 3)
# round(seq(0, 1, length.out = 10), 3)

## draw a ROC curve first for calibration performance test
ivTruth = fGroups == 'ATB'
p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
plot(perf.alive)

## draw the simulation lines
## these are p-values from the mixture components
## create posterior smatter lines
grid = seq(-7.5, 7, length.out = 100)
f_getSmatterLines = function(m, s, g){
  return(pnorm(g, m, s, lower.tail = F))
}
y = f_getSmatterLines(2.27, 1.96, grid)
x = f_getSmatterLines(-2.57, 1.58, grid)
lines(x, y, col=2)

## holders for the simulated p-values
mTP = matrix(NA, nrow = length(grid), ncol = 2000)
mFP = matrix(NA, nrow = length(grid), ncol = 2000)

for (i in 1:2000){
  p = sample(1:nrow(mStan), size = 1)
  x = pnorm(grid, mStan[p, 'mu1'], mStan[p, 'sigma1'], lower.tail = F) 
  y = pnorm(grid, mStan[p, 'mu2'], mStan[p, 'sigma2'], lower.tail=F)
  lines(x, y, col='darkgrey', lwd=0.5)
  mFP[,i] = x
  mTP[,i] = y
}

plot(perf.alive, add=T, col='blue')

### check the p-values from the model to select a cutoff point to reduce false positive rate
### i.e. to test the class at certain cutoff of the score and then declare it as false match
### the others will go for rejection category
### see gelman page 19 for such a calibration method.

## third way of using cross validation to draw the ROC curves
## this method uses LDA but results are pretty comparable to binomial model
dfData = dfData.new[,-8]
#colnames(dfData) = CVariableSelection.ReduceModel.getMinModel(oVar.subComb, 7)
dim(dfData)

oCV = CCrossValidation.LDA(test.dat = dfData, train.dat = dfData, test.groups = fGroups,
                           train.groups = fGroups, level.predict = 'ATB', boot.num = 100)

plot.cv.performance(oCV)

## compare these 3 tables to decide on a cutoff
## the simulation suggests somewhere around 0.6 should be the decision cutoff
df = getCutoffTprFprCrossValidation(oCV)
dfSim = round(data.frame(logit.inv(grid), fp=rowMeans(mFP), tp=rowMeans(mTP)), 3)
head(dfPerf.alive)

fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.6] = 'ATB'
table(fPredict, fGroups)

## draw these accept reject points
xyplot(ivPredict ~ fGroups, xlab='Actual Group', ylab='Predicted Probability of ATB (1)', groups=fPredict,
       auto.key = list(columns=2))

xyplot(ivPredict ~ lData.train$grouping, xlab='Actual Group', ylab='Predicted Probability of ATB (1)', groups=fPredict,
       auto.key = list(columns=2))

p = sample(1:nrow(mStan), size = 2000)
x = pnorm(logit(0.6), mStan[p, 'mu1'], mStan[p, 'sigma1'], lower.tail = F)
y = pnorm(logit(0.6), mStan[p, 'mu2'], mStan[p, 'sigma2'], lower.tail=F)

par(p.old)
par(mfrow=c(1,2))
hist(x, xlab='Expected Rate', main='False Positive Rate at cutoff')
hist(y, xlab='Expected Rate', main='True Positive Rate at cutoff')
# at this cutoff what is the expect false positive and true positive rate according to the model
cvTopGenes.chosen = colnames(dfData)

################################ end section of choosing model on training data


############################################ Test data set
## load the test data and try these combinations
lData.test = f_LoadObject('GSE19491_normalised_train.RList')
fGroups.test = rep('ATB', times=length(lData.test$grouping))
fGroups.test[lData.test$grouping != 'ATB'] = 'Other'
fGroups.test = factor(fGroups.test, levels = c('Other', 'ATB'))
table(fGroups.test)
table(lData.test$grouping)

dfData = data.frame(t(lData.test$data))
dim(dfData)
dfData = dfData[,cvTopGenes.chosen]
rm(fGroups)
head(dfData)

oCV.test = CCrossValidation.LDA(dfData, dfData, fGroups.test, fGroups.test, level.predict = 'ATB', 100)
plot.cv.performance(oCV.test)

d = data.frame(dfData, fGroups.test)
f = glm( fGroups.test ~ ., data=d, family='binomial')
ivPredict = predict(f, type = 'response')

fPredict = rep('reject', times=length(ivPredict))
fPredict[ivPredict >= 0.6] = 'ATB'
table(fPredict, fGroups.test)

## draw these accept reject points
xyplot(ivPredict ~ fGroups.test, xlab='Actual Group', ylab='Predicted Probability of ATB (1)', groups=fPredict,
       auto.key = list(columns=2))

xyplot(ivPredict ~ lData.test$grouping, xlab='Actual Group', ylab='Predicted Probability of ATB (1)', groups=fPredict,
       auto.key = list(columns=2))





