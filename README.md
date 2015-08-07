# CGraphClust
Class to create a Hclust object based on shared properties (e.g. pathways) and positive correlation

# Inheritence Structure
CGraph - Parent Class --> CGraphClust - Child Class

# CGraph Class
assigns weights to one mode projection of graphs based on observed to expected probabilities of vertices of the first kind
i.e. with value TRUE using the igraph library
Zweig, K. A., & Kaufmann, M. (2011). A systematic approach to the one-mode projection of bipartite graphs. 
Social Network Analysis and Mining (Vol. 1, pp. 187–218). doi:10.1007/s13278-011-0021-0

# Constructor - CGraph 
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
Arguments - a data frame of 2 columns and n rows, where column 1 represents vertices of type 1 (e.g. Gene symbols) and column 2
represents vertices of type 2 (e.g. Pathways related to the corresponding gene symbol). A correlation matrix of the expression
values of e.g. the Genes with column and rownames being the same as the entries in column 1 of the data frame. The third argument
is a cutoff value of correlations (default is 0.5) - i.e. anything below 0.5 is not considered to be related. 
# step 1
bipartite graph is created and some cleaning performed, i.e. those type 2 vertices with too many connections (degree) with type 1
vertices are removed. this is done by looking at the distribution of the degrees on a log scale and approximating a negative
binomial and poisson distributions on top. the cutoff is by default 0.95 quantile under a negative binomial model, anything over that
is removed. The idea is that type vertices that have a lot of connections, are not very interesting and they tend to hide the
more interesting connections.
# step 2
graph is projected on to one dimension (type 1 vertices) and connections between type 1 vertices are assigned interestingness or 
observed to expected probability ratios as weights. the weights on a square root scale follow a negative binomial distribution
and only the highest weights are chosen, as they are the most interesting. So anything over 0.95 quantile under a neg bin model
are chosen and the other connections are discarded.
# step 3
a correlation matrix is converted to an adjacency matrix with a 0 diagonal and converted to a graph. Only those edges greater than
or equal to the cutoff value (default 0.5) are chosen.
# step 4
the 2 graphs are intersected and only those connections are remaining that are positively correlated and have a very high observed 
to expected ratios. We remove small clusters by plotting the distribution of log cluster sizes and keep only the large sized clusters.
edge betweeness community finding algorithm is used to detect clusters in this graph. NOTE: if the number of edges is over 3000 a 
simpler community finding algorithm is used and a message is written out.

# step 5
the community is converted to a hclust object which can be used for plotting and heatmaps. each cluster is also assigned a 
label based on which is the most common type 2 vertex (e.g. pathway) in that cluster and can be accessed via the
getClusterLabels function.

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

# plot.heatmap.marginal
similar to plot heatmap but will plot the marginal profile of type 1 vertices in each cluster, i.e. by performing a column sum
on the data matrix, where the rows represent the type 1 vertices in the cluster and the columns are the individual samples (components) of the vector.  

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
Very similar to plot.mean.expressions, but it will first optionally stabalize the data (default FALSE). See the function 
f_ivStabilizeData for details. The main difference to the previous function is that it will call getSignificantClusters and use the
names of significant clusters to plot only those.

# plot.cluster.expressions
Plot a line graph of expressions of all the members of the cluster, in order of the given
columns of the count matrix. These patterns make more sense if the graph was built on the
correlation matrix instead of absolute correlation values.

# plot.heatmap.cluster
very similar to the previous function, but plots heatmap instead of line graph


# getClusterMarginal
ARGS: object of CGraphClust class, mCounts = count matrix with rows as type 1 vertices (genes) and columns as samples;  
bScaled = TRUE (default) - will scale the mCounts matrix (i.e. z scaled with zero center) before calculations. This can be useful 
as it highlights the direction of effect.  
The function does some error checking first to see if the matrix row names match with the type 1 vertex names. It then reorders the
count matrix based on its order in the hclust object (using the type 1 vertex names). Each cluster was assigned a label earlier 
which was the most common shared type 2 vertex in that cluster (e.g. the most common REACTOME term for a cluster of genes). For each
set of genes in the cluster based on the cluster label (we haven't seen it yet, but it may be possible that a cluster label may 
be assigned to more than one clusters?) - we take the mean of the scaled (or not scaled) expression values for each gene together, i.e. taking the marginal of the contingency table, averaging across the columns (compressing the columns) to create one total expression row for each cluster.


# f_ivStabilizeData and f_mCalculateLikelihoodMatrix
Utility functions for count matrix stabalization, if you expect that the data groups have outliers. The function sets the seed before doing anything else in order to produce reproducible results. This will use the overall 
data variance as a fixed prior variance, and create a prior distribution of the means by sampling from the normal distribution with
the parameters fixed data variance (prior) and X means (number of means = number of levels). Using this prior distribution of 
possible means, a likelihood vector is calculated for each possible mean, where the data parameter is the mean for the group, and
the standard error of mean is the SD parameter of normal density function. Using this likelihood vector as a probability of 
sampling (similar to rejection sampling) - we sample from the prior to create a posterior set of means. To create a new set of 
data (posterior predictive) for each group - we sample from the posterior mean, and using that and the SE (not SD) for a rnorm
function we sample new data.

# plot.components
Performs a principal component analysis on the correlation matrix of the clusters. The data is first stabalized using 
f_ivStabalizeData (set bStabalize=FALSE to not do this), then getClusterMarginal is used to get marginal data for each cluster. The 
significant clusters are selected after using the function getSignificantClusters (to select only clusters that show a significant 
change between the conditions using ANOVA) and this matrix is used for PCA where the row vectors are clusters and the columns are
the components (i.e. sample measurements). The first 2 components are plotted and the prcomp object is returned.

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


# plot.centrality.diagnostics'
The function creates an ordered statistic of each centrality measure, and makes a bar plot, while colouring each bar with the
colour of the cluster of that value. It then groups the vector into 0-90% and 91-100% quantiles, subsets the data over 90% quantile
and draws another bar plot with cluster wise colours, but only for the top 10% quantile. This kind of plot can highlight which
centrality measure may be important.  
RETS: returns the list of data frames with the top 10% vertices and its associated cluster label.


# lGetTopVertices
Simple summary function, reports the names of vertices in:  
1- The vertices in the largest clique.  
2- The top Quantile of the vertices from the centrality scores - default is 0.95, change iQuantile argument.  
The data is returned in a list format.  

# getSignificantClusters
ARGS: the function will take the graph object, count matirx (rows are type 1 vertices e.g. genes), grouping factor (representing 
columns) and a boolean value (if to perform data stabalization).  
The function calculates the marginal for each cluster using getClusterMarginal function. Then performs an anova on each vector of cluster values grouping them on the fGroups factor. Multiple testing adjustment method used is 'BH' and at a FDR of 0.01 the significant clusters are selected. Each group mean for a cluster vector is used to calculate the group means 
and then the difference in the maximum and miminum mean to rank the clusters on that. The return value is the marginal matrix 
sorted on the ranking of clusters and the p values.  
RETS: list

# getLargestCliqueInCluster
Takes a cluster label, finds the largest clique in that cluster and returns an igraph object

# f_igCalculateVertexSizesAndColors
Utility function to calculate vertex sizes of nodes of the graph, using the expression data, grouping factor and number of nodes
in the graph. Node sizes are calculated by log10 fold change between the mean of 2 groups (in case of more than 2 groups, the 
means of 2 groups at the extremes are used using the range function). This size factor is multiplied by a iSize (default NULL) 
argument. If it is provided (set value to around 20 or 30 for large graphs - depending on drawing area) - for small graphs leave 
this area blank and a multiplicative factor is calculated using 4000/vcount(ig). 
If the bColors variable is set (default FALSE), then factor with levels(factor)[1] is used as baseline and levels(factor)[last]
is used as the second level, and if mean of levels(factor)[1] is < mean levels(factor)[length(levels(factor))], it suggests 
a positive change with a colour of red (pink) or else negative change with blue colour


# f_csGetGeneSummaryFromGenbank
Utility function using the bioconductor package annoate and cran package XML to query genbank for a gene, download the xml file
and extracts the summary for the gene returning it in a character string format. It can take a vectorized form of the gene list and
returns a vectorized named character string













