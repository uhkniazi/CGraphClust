# CGraphClust
R Class to create graphs from Omics data, e.g. Gene expression matrix and a Table of Genes and Associated pathway terms.  

# Quick Start Guide  
The following are required to be present on the system:  
+ R version 3.4.0 or higher  
+ *igraph* package in R  
+ *downloader* package in R (if downloading annotation files required)  
+ Bioconductor, (version 3.0 or higher) with packages  
  * org.Hs.eg.db (if annotations for human genes is required)  
 
 
## Running Example  
The usage is explained via an example.  
### Data Preparation
The data matrix consists of 21 samples and 1362 genes, while the grouping factor shows 7 samples in each group.

```R
> dim(mCounts)
[1]   21 1362
> table(fGroups)
fGroups
12  2  0 
 7  7  7 
```
The second piece of information we need is a mapping table that has Variable names (Gene IDs) in the first column and Grouping terms (e.g. REACTOME ids for that gene).  
```R
> mCounts[1:5, 1:10]
               2597     9792    65980     1592     8720     2790     3779     2000     6398   222223
GSM484380 10.158587 8.371591 7.812608 6.667878 7.604684 8.327295 7.709530 9.491233 5.543011 4.165178
GSM484381  9.939359 8.877825 7.948275 6.030240 7.894243 8.174115 7.152345 8.806329 5.266879 4.492998
GSM484382  9.491808 9.333348 8.542408 5.298431 8.084001 6.839688 6.253693 8.620620 3.947724 4.729314
GSM484383 10.513044 8.368436 7.524928 6.451396 7.919592 8.666990 7.588009 9.189281 4.962316 4.237957
GSM484384  9.568857 9.054925 8.204923 5.436528 8.069337 8.066015 6.575776 8.679335 4.415582 4.252351
> head(dfGraph, 10)
      ENTREZID      V2
30887    57488 1430728
30888    57488 1660662
30889    57488  428157
30890    57488  556833
33498   338440 2672351
33499   338440  382551
33500   338440  983712
51717    23019  212436
51718    23019 3700989
51719    23019  429914
```
In the code snippet below, we subset the count matrix, by the number of genes with actual REACTOME terms attached.
```R
# select genes that have a reactome term attached
# this will reduce your initial gene list as a lot of genes actually have 
# no annotations, you may use a low level database like GO if you think 
# that you lose too many genes
n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
[1] "Total number of genes with Reactome terms 823"
# order the count matrix based on grouping factor
# this will serve useful later for plotting but is not necessary
rownames(mCounts) = fGroups
mCounts = mCounts[order(fGroups),]
fGroups = fGroups[order(fGroups)]

head(fGroups)
GSM484382 GSM484385 GSM484392 GSM484393 GSM484394 GSM484397 
       12        12        12        12        12        12 
Levels: 12 2 0
```
### Building the Graph  
The constructor for the CGraphClust object requires a matrix of absolute correlations of genes and the mapping data frame with Gene Ids and Pathway labels.  
```R
mCor = cor(mCounts)
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.7, bSuppressPlots = F)
```
If everything goes well you should have the graph object in the variable called oGr. This can now be used for further analysis.  

### Graph Pruning  
Due to the nature of the data, sometimes you will have graphs with very small clusters, which will contain 2 to 5 genes, however most clusters will be pretty large e.g. 10 to 90 plus genes. Small clusters can be removed by identifying these genes and dropping them from the graph.
```R
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
iSizes = sort(table(dfCluster$cluster))
iSizes

 170834   71387  112316  397014   72306 1266738  194315 2454202   15869 1852241  379724 1428517   72312  109582 
      2       2       3       4       4       6       6       6       8       8       8       9      10      12 
 382551 2262752  196849  372790   72203 5653656  392499  556833 1280218 1280215  212436 6798695 
     13      14      15      18      28      35      36      36      44      51      66      87 
# remove communities smaller than 5 members
i = which(iSizes <= 5)
if (length(i) > 0) {
  cVertRem = as.character(dfCluster[dfCluster$cluster %in% names(i),'gene'])
  iVertKeep = which(!(V(getFinalGraph(oGr))$name %in% cVertRem))
  oGr = CGraphClust.recalibrate(oGr, iVertKeep)
}
```
### Analysis  
The types of analyses that can be performed on the CGraphClust object include:  
1- Individual types of igraph objects can be extracted using the functions: getBipartiteGraph, getCorrelationGraph and getFinalGraph.  
2- Centrality summaries of Degree, Closeness, Betweenness and Hub centralities can be obtained using plot.centrality.graph, plot.centrality.diagnostics, mPrintCentralitySummary, dfGetTopVertices.  
3- Information regarding clusters/communities in the graph can be obtained by getClusterMapping, plot.graph.clique, getCommunity, getLargestCliques, getSignificantClusters, getHclust.  
4- Multiple plotting functions, e.g. plot.mean.expressions, plot.significant.expressions, plot.components etc.  

For typical usage where we want the list of genes that are the top 10% with respect to at least one centrality score, we can use the following approach:  
```R
## top vertices/nodes/genes based on centrality scores
## get a table of top vertices 
dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
head(dfTopGenes.cent)
       VertexID clique degree closeness betweenness   hub
4940       4940   TRUE  FALSE     FALSE       FALSE FALSE
2633       2633   TRUE   TRUE     FALSE       FALSE FALSE
2634       2634   TRUE   TRUE     FALSE       FALSE FALSE
163351   163351   TRUE   TRUE     FALSE       FALSE FALSE
10346     10346   TRUE   TRUE     FALSE       FALSE FALSE
115362   115362   TRUE   TRUE     FALSE       FALSE FALSE
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)
head(dfTopGenes.cent)
       VertexID clique degree closeness betweenness   hub SYMBOL                                  GENENAME Cluster
4940       4940   TRUE  FALSE     FALSE       FALSE FALSE   OAS3         2'-5'-oligoadenylate synthetase 3 1280215
2633       2633   TRUE   TRUE     FALSE       FALSE FALSE   GBP1               guanylate binding protein 1 1280215
2634       2634   TRUE   TRUE     FALSE       FALSE FALSE   GBP2               guanylate binding protein 2 1280215
163351   163351   TRUE   TRUE     FALSE       FALSE FALSE   GBP6 guanylate binding protein family member 6 1280215
10346     10346   TRUE   TRUE     FALSE       FALSE FALSE TRIM22            tripartite motif containing 22 1280215
115362   115362   TRUE   TRUE     FALSE       FALSE FALSE   GBP5               guanylate binding protein 5 1280215
```
The scripts in the [Paper folder](https://github.com/uhkniazi/CGraphClust/tree/master/Paper) contain a working example to generate test data and build the respective graphs; while the [Test_data](https://github.com/uhkniazi/CGraphClust/tree/master/Test_data) folder contains various other test examples.  
  
  
The following section contains details about each function in the class.


***
# Inheritence Structure
CGraph - Parent Class --> CGraphClust - Child Class

# CGraph Class
assigns weights to one mode projection of graphs based on observed to expected probabilities of vertices of the first kind
i.e. with value TRUE using the igraph library
Zweig, K. A., & Kaufmann, M. (2011). A systematic approach to the one-mode projection of bipartite graphs. 
Social Network Analysis and Mining (Vol. 1, pp. 187–218). doi:10.1007/s13278-011-0021-0

# Constructor
##CGraph 
Arguments - bipartite graph object of igraph library, with vertices of the first kind assigned id TRUE and of second type
assigned a FALSE. 

# Internal Functions
# CGraph.assign.marginal.probabilities
vertices of the first kind are assigned probabilities
# CGraph.project
assigns a level of interestingness/leverage or observed to expected ratio to each edge after graph projection on the vertex of first kind i.e. type = TRUE.

# data acccessor functions
getBipartiteGraph & getProjectedGraph

# CGraphClust Class
class to create a igraph and hclust object based on 2 criteria: 1) shared properties or connections with type 2 vertices in a bipartite graph and 2) positive correlation value.

# Constructor - CGraphClust
## Arguments
dfGraph = a data frame of 2 columns and n rows, where column 1 represents vertices of type 1 (e.g. Gene symbols) and column 2
represents vertices of type 2 (e.g. Pathways related to the corresponding gene symbol).  
mCor = A correlation matrix of the expression values of e.g. the Genes with column and rownames being the same as the entries in column 1 of the data frame.  
iCorCut = 0.5 DEFAULT, a cutoff value of correlations (default is 0.5) - i.e. anything below 0.5 is not considered to be related.  
bSuppressPlots = TRUE DEFAULT, should the plots be drawn while building graph.  
iMinComponentSize = 6 DEFAULT, Before finding the communities, the smallest strongly connected components less than the this size 
are removed, this can be adjusted, if a graph has many small communities.
### step 1
bipartite graph is created and some cleaning performed, i.e. those type 2 vertices with too many connections (degree) with type 1
vertices are removed. this is done by looking at the distribution of the degrees on a log scale and approximating a negative
binomial and poisson distributions on top. the cutoff is by default 0.95 quantile under a negative binomial model, anything over that
is removed. The idea is that type vertices that have a lot of connections, are not very interesting and they tend to hide the
more interesting connections.
### step 2
graph is projected on to one dimension (type 1 vertices) and connections between type 1 vertices are assigned interestingness or 
observed to expected probability ratios as weights. the weights on a square root scale follow a negative binomial distribution
and only the highest weights are chosen, as they are the most interesting. So anything over 0.95 quantile under a neg bin model
are chosen and the other connections are discarded.
### step 3
a correlation matrix is converted to an adjacency matrix with a 0 diagonal and converted to a graph. Only those edges greater than
or equal to the cutoff value (default 0.5) are chosen.
### step 4
the 2 graphs are intersected and only those connections are remaining that are positively correlated and have a very high observed 
to expected ratios. We remove small clusters by plotting the distribution of log cluster sizes and keep only the large sized clusters.
edge betweeness community finding algorithm is used to detect clusters in this graph. NOTE: if the number of edges is over 3000 a 
simpler community finding algorithm is used and a message is written out.

### step 5
the community is converted to a hclust object which can be used for plotting and heatmaps. each cluster is also assigned a 
label based on which is the most common type 2 vertex (e.g. pathway) in that cluster and can be accessed via the
getClusterLabels function.


# Constructor - CGraphClust.recalibrate
##ARGS
obj = object of CGraphClust class  
ivVertexID.keep = integer vector of vertex IDs that will not be removed from the obj  
iMinComponentSize=6 If any components are smaller than this, those will be removed, used to remove small sized subgraphs not connected to the main graph  

##DESC
This constructor is not called directly first time when creating a new CGraph object. It will create a new CGraph object that only contains the subset of the vertices. It works in the following steps:  
1- get the subgraph (igraph object) and delete any orphan vertices i.e. with degree 0.  
2- remove any small components.  
3- get communities in the new graph, assign type 2 vertex names to each cluster and create a new CGraph object.  
Typically this function can be used to remove those communities which have a small number of members, or to create a new graph after intersecting two graphs from different datasets (used in data integration).  

##RETS
CGraph object.

# Constructor - CGraphClust.intersect.union
##ARGS
obj1 = object of CGraphClust class  
obj2 = object of CGraphClust class  
iMinOverlap=10 minimum overlap between nodes of obj1 and obj2
iMinComponentSize=6 If any components are smaller than this, those will be removed, used to remove small sized subgraphs not connected to the main graph  

##DESC
A direct union using the igraph library function will produce a union of two graphs and even if the vertices are not common, they still be included in the final graph. This constructor is useful typically in situations where a direct intersection and CGraphClust.recalibrate does not have many nodes remaining in the graph. In that case, this function may be useful, however it should be used with care and the results not over interpreted. It however allows to compare two very different datasets, e.g. a tuberculosis and sepsis datasets.  
1- Get vertices common between the two graphs.  
2- create two induced subgraphs.  
3- create a union of the two graphs.  
4- recalculate the edge weights (observed to expected) as the mean of the weights provided by each graph. however some edges may not be present in one graph and that will contribute NA edge weight.  
5- Recreate the bipartite graph using the steps from the original constructor.  
NOTE: most of the components of this new CGraph object are meaningless, and only the following components are used:  
hc = hclust community object  
com = community object  
labels = the community labels (most common type 2 vertex) assigned to each vertex  
ig.i = the graph after the union, use getFinalGraph function to get this.  


# data accessor functions
getProjectedGraph, getCorrelationGraph, getBipartiteGraph, getFinalGraph (to get the intersected graph), 
getClusterLabels, getHclust, getCommunity

# oCGdeleteSmallCommunities (REMOVED)
if we want to remove vertices/communities below a certain size. this can be useful to remove very small clusters which may
not really give useful information. it takes as an argument the object of CGraphCluster and a size (default = 3), and removes
any clusters with number of vertices equal to or less than size. It returns a new object of CGraphCluster size, which has
the bipartite, correlated and projected graphs similar to the original objects, however the final graph (intersected graph),
community object, labels and hclust object are different due to removing smaller communities.

# getClusterMapping
returns a data.frame with mapping type 1 vertex names to corresponding type 2 vertex name that has been used to assign cluster
name.

# getClusterMatrix
returns the subset of the matrix that contains only the data from the named cluster

# getLargestCliques
returns the largest number/s of cliques in a list with vertex names

# plot.final.graph and plot.graph.clique
simple graph plotting functions to give an overview of the graph, uses the Fruchterman-Reingold layout algorithm.

# plot.centrality.graph
plots the graph with the nodes colored on centrality. Returns the igraph object.

# plot.heatmap.all
plots an annotated heatmap using the NMF library. Requires the object CGraphClust, count matrix with type 1 vertices in the 
rows and type 2 vertices in the columns and a default cutoff expression values (-3, 3) - all values below or above these cutoffs
are thredholded to this value as the extreme values affect the heatmap colouring. The heatmap matrix is scaled by sd and centered 
to a zero mean along the rows (type 1 vertex) before plotting.

# plot.heatmap.significant.clusters
similar to plot heatmap all, but uses the getSignificantClusters function to get the data matrix for the marginals of clusters
showing significant expression.

# plot.heatmap.marginal
similar to plot heatmap but will plot the marginal profile of type 1 vertices in each cluster, i.e. by calling getClusterMarginal
on the data matrix. 

# plot.mean.expressions
plots the mean expression of each group (based on factor fGroups) in each cluster. The function takes the object of CGraphClust type
a count matrix with rows as type 1 vertices and columns as components of the row vectors (they would usually be your samples), 
a grouping factor defining the groups, a legend position (default it topright). The function checks the final projected graph with
type 1 vertices, and uses the names of those vertices to subset the count matrix rows. It assigns cluster labels to each gene and 
uses that cluster label as a grouping factor to calculate the mean values in that group (similar to function plot.heatmap.means) - i.e. it is taking the rowmeans of the count matrix based on each cluster grouping factor. It then takes that matrix of cluster means and
takes the mean and standard deviations for each cluster vector - grouped on the second factor (fGroups). Effectively we are compressing the data using getClusterMarginal on the count matrix (grouped on clusters) and second rowMeans of the matrix (grouped on samples e.g. control vs treatment).
The data is plotted and a list with the matrix of mean and standard deviations is returned.

NOTE: Ideally one would want to reassign cluster labels, instead of using the most common label in a group. This can be done
by using the function getClusterMapping, and looking for the type 1 vertices associated with a cluster label. Those type 1 vertex
names (e.g. gene names in that cluster) can be put into reactome database to see which is the most sensible label to assign to that
cluster.

# plot.significant.expressions
Very similar to plot.mean.expressions, but it will optionally stabalize the data before plotting (default FALSE). See the function 
f_ivStabilizeData for details. The main difference to the previous function is that it will call getSignificantClusters and use the
names of significant clusters to plot only those.

# plot.heatmap.cluster
very similar to the previous function, but plots heatmap instead of line graph


# getClusterMarginal
## ARGS:
obj = object of CGraphClust class.  
mCounts = count matrix with rows as type 1 vertices (genes) and columns as samples.  
bScaled = FALSE (default) - will scale the mCounts matrix (i.e. z scaled with zero center) before calculations. This can be useful 
as it highlights the direction of effect.
## DESC:
The function does some error checking first to see if the matrix row names match with the type 1 vertex names. It then reorders the
count matrix based on its order in the hclust object (using the type 1 vertex names). Each cluster was assigned a label earlier 
which was the most common shared type 2 vertex in that cluster (e.g. the most common REACTOME term for a cluster of genes). If a cluster has only one member then remove it from calculation.  For each set of genes in the cluster based on the cluster label (we haven't seen it yet, but it may be possible that a cluster label may be assigned to more than one clusters?) - we use the internal function f_edge.score to calculate the score vector for each cluster and the column means of this matrix to make a vector, where the components of the vector are the number of columns in the matrix mCounts (i.e. samples). 
### f_edge.score (internal function)
NOTE: before calling this function make sure that the graph has more than one edges
Gets the subgraph of each cluster, gets the edge list of this subgraph, gets 2 expression vectors from mCounts for the 2 nodes connecting the edges, adds these together. (Considering weighting them by observed to expected probability ratio as well). Returns a matrix of scores where rows = number of edges and columns = number of samples   
## RETS:
The expression in each cluster


# plot.components
Performs a principal component analysis on the correlation matrix of the clusters. The count matrix is marginalized by getClusterMarginal to get marginal data for each cluster. The significant clusters are selected after using the function getSignificantClusters (to select only clusters that show a significant change between the conditions using ANOVA) and this matrix is used for PCA where the row vectors are clusters and the columns are the components (i.e. sample measurements). The first 2 components are plotted and the prcomp object is returned. The data is stabalized (optional) using f_ivStabalizeData (set bStabalize=FALSE to not do this) before performing PCA, if the data has a high variance.

# getClusterSubgraph
Takes a cluster label and returns the igraph object for that subgraph

# mPrintCentralitySummary 
Centrality can be measured by various methods, the function returns a matrix with various centrality measures for each vertex.  
1- Degree: the number of connections for each vertex.  
2- Closeness: how close a vertex is to other vertices in the graph.  
3- Betweenness: how many shortest paths between other vertices pass through the given vertex.  
4- Hub: hub score or authority score: http://www.cs.cornell.edu/home/kleinber/auth.pdf - these are nodes that have a large in-degree
 and a lot of overlap in the sets of nodes that point to them. Techically they are calculated by http://igraph.org/r/doc/authority_score.html  
The correlation should be small between these metrics, and the function reports a correlation matrix and returns the measures in
a matrix.
  
some examples can be seen here: http://cs.brynmawr.edu/Courses/cs380/spring2013/section02/slides/05_Centrality.pdf  
http://www.evernote.com/l/ASAPToyPhMRNM6235zEpyVJdU1KPSISX5Do/  

# getCentralityMatrix
similar to mPrintCentralitySummary but just returns the matrix


# plot.centrality.diagnostics
The function creates an ordered statistic of each centrality measure, and makes a bar plot, while colouring each bar with the
colour of the cluster of that value. It then groups the vector into 0-90% and 91-100% quantiles, subsets the data over 90% quantile
and draws another bar plot with cluster wise colours, but only for the top 10% quantile. This kind of plot can highlight which
centrality measure may be important.  
RETS: returns the list of data frames with the top 10% vertices and its associated cluster label.

# plot.cluster.variance
## ARGS:
obj = object of CGraphClust.  
mCent = matrix of data for which variance is to be calculated. Typically this would be Row Vectors, where columns are components.  
fGroups = grouping factor for components of mCent columns.  
log=TRUE default value, for drawing on log scale
iDrawCount=4 Default value for how many rows will be drawn (as too many rows will not be visible on the lattice box plots).  
## DESC:
Uses the package lattice::bwplot to draw the distribution of the posterior variance of each group (columns) in the data matrix, grouped into sub panels by Row vectors.


# lGetTopVertices
Simple summary function, reports the names of vertices in:  
1- The vertices in the largest clique.  
2- The top Quantile of the vertices from the centrality scores - default is 0.95, change iQuantile argument.  
The data is returned in a list format.  

# dfGetTopVertices
Very similar to previous function lGetTopVertices, but data is returned in a data frame format.  
This can be useful for quickly getting a summary of genes and printing them out.

# getSignificantClusters
## ARGS:
the function will take the graph object, count matirx (rows are type 1 vertices e.g. genes), grouping factor (representing 
columns), p.cut = 0.001 (default), the p-value to use as cutoff
## DESC: 
The function calculates the marginal for each cluster using getClusterMarginal function. For each cluster and the groups in that cluster, using a non-informative prior distribution, 100000 samples are generated for each group (using f_lpostMean function). Using these simulated samples, a p-value is calculated by comparing the posterior distribution of the mean of the baseline group and a two sided test for means of the other groups. If a simulation p-value for any of the comparison is less than p.cut, than that cluster is marked as showing significant differences from the baseline. Each group mean for a cluster vector is used to calculate the group means and then the difference in the maximum and miminum mean to rank/sort the clusters on that. The return value is a list with the marginal matrix sorted on the ranking of clusters and the p values.  
## RETS:
list with significant clusters in matrix and p values

# getLargestCliqueInCluster
Takes a cluster label, finds the largest clique in that cluster and returns an igraph object

# f_igCalculateVertexSizesAndColors
Utility function to calculate vertex sizes of nodes of the graph, using the expression data, grouping factor and number of nodes
in the graph. Node sizes are calculated by log2 fold change between the mean of 2 groups (the 2 extreme levels are used i.e. the base line level of the factor and the last level of the factor to calculate the difference). This size factor is multiplied by a iSize (default NULL) argument and absolute value is saved. If it is provided (set value to around 20 or 30 for large graphs - depending on drawing area) - for small graphs leave this area blank and a multiplicative factor is calculated using 4000/vcount(ig). 
If the bColors variable is set (default FALSE), then factor with levels(factor)[1] is used as baseline and levels(factor)[last]
is used as the second level, and if mean of levels(factor)[1] is < mean levels(factor)[length(levels(factor))], it suggests 
a positive change with a colour of red (pink) or else negative change with blue colour

# f_csGetGeneSummaryFromGenbank
## NOTE
This function does not seem to work well in vectorized format, needs fixing at some point.
Utility function using the bioconductor package annoate and cran package XML to query genbank for a gene, download the xml file
and extracts the summary for the gene returning it in a character string format. It can take a vectorized form of the gene list and
returns a vectorized named character string


# f_ivStabalizeData
## ARGS:
ivDat = integer vector for data  
fGroups = corresponding grouping factor for ivDat  
## DESC:
Using a weakly informative conjugate prior, calculated from ivDat, where the prior n=2 and prior degrees of freedom = 1. Look at page 68 of Bayesian Data Analysis (Gelman 2008) for the formula. The posterior mean is simulated using the data for each group and the prior, where the posterior variance is inverse chi-square and posterior mean given the posterior variance/n1. 
## RETS:
Returns a data vector of size ivDat, and order fGroups, where each value is a sample from the respective posterior mean.  
This function is typically used to stabalize the variance in the data by pushing the measurements towards the mean.

# f_lpostVariance
## ARGS:
ivDat = integer vector for data  
fGroups = corresponding grouping factor for ivDat  
## DESC:
very similar to f_ivStabalizeData.
## RETS:
Vector of size 1000 with simulated values of posterior variance.

# f_lpostMean
## ARGS:
ivDat = integer vector for data  
fGroups = corresponding grouping factor for ivDat  
sim.size = 100000
## DESC:
very similar to f_ivStabalizeData, but will use a large sample size, as it is used to calculate p-values
## RETS:
List of components = levels(fGroups) and length of each component is = sim.size with simulated values of posterior mean for each group



# REDUNDANT - REMOVED
## plot.cluster.expressions
Plot a line graph of expressions of all the members of the cluster, in order of the given
columns of the count matrix. These patterns make more sense if the graph was built on the
correlation matrix instead of absolute correlation values.

## f_ivStabilizeData and f_mCalculateLikelihoodMatrix
Utility functions for count matrix stabalization, if you expect that the data groups have outliers. The function sets the seed before doing anything else in order to produce reproducible results. This will use the overall 
data variance as a fixed prior variance, and create a prior distribution of the means by sampling from the normal distribution with
the parameters fixed data variance (prior) and X means (number of means = number of levels). Using this prior distribution of 
possible means, a likelihood vector is calculated for each possible mean, where the data parameter is the mean for the group, and
the standard error of mean is the SD parameter of normal density function. Using this likelihood vector as a probability of 
sampling (similar to rejection sampling) - we sample from the prior to create a posterior set of means. To create a new set of 
data (posterior predictive) for each group - we sample from the posterior mean, and using that and the SE (not SD) for a rnorm
function we sample new data.









