# generate_test_data.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 23/04/2015
# Desc: analysis of dataset to generate test data

library(GEOquery)
library(Biobase)
#library(annotate)
#library(org.Hs.eg.db)
library(lumi)
#source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')

# global variables
p.old = par()

## data loading
# load the data, clean and create factors
gse =  getGEO('GSE63881', GSEMatrix = T, destdir = 'Data_external/')
oExp = gse$GSE63881_series_matrix.txt.gz

# print samples
dfSam = as.data.frame(table(oExp$characteristics_ch1.1))

## data normalization
# normalize and log2 transform the data using lumi
oExp.lumi = oExp; #lumiT(oExp, 'log2')
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))
fSamples = factor(as.character(oExp$characteristics_ch1.1), levels = c("phase: Acute",
                                                                 "phase: Convalescent"),
                  labels = c('AC', 'CONV'))  
oExp.lumi$fSamples = fSamples

## data normalization
oExp = lumiN(oExp.lumi, method='rsn')
rm(oExp.lumi)
gc()

## extract data for significant genes
fSamples = oExp$fSamples
pano = apply(exprs(oExp), 1, function(x) anova(lm(x~fSamples))$Pr[1])
pano.adj = p.adjust(pano, method = 'BH')
n = which(pano.adj < 0.01)
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

write.csv(dfData, file='Test_data/test_data_GSE63881.csv')