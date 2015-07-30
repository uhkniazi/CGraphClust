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
gse =  getGEO('GSE19491', GSEMatrix = T, destdir = 'Data_external/')
oExp = gse$GSE19491_series_matrix.txt.gz

# print samples
as.data.frame(table(oExp$source_name_ch1))

# get the whole blood data
i = grep('Whole', x = oExp$source_name_ch1, ignore.case = F, perl = T)
oExp = oExp[,i]

## data normalization
# normalize and log2 transform the data using lumi
oExp.lumi = lumiT(oExp, 'log2')
# remove any NA data
exprs(oExp.lumi) = na.omit(exprs(oExp.lumi))
fSamples = rep(NA, 180)
f = as.character(oExp$source_name_ch1)
# create factors
i = grep('healthy', f)
fSamples[i] = 'HC'

i = grep('Latent', f)
fSamples[i] = 'LTBI'

i = grep('Active TB', f, ignore.case = F)
fSamples[i] = 'ATB'

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

#write.csv(dfData, file='Test_data/test_data_GSE19491_full.csv')

# save only the 0, 2, 12 month for cluster analysis test data
i = which(fSamples %in% c('0', '2', '12'))

dfData = dfData[i,]
write.csv(dfData, file='Test_data/test_data_GSE19491.csv')