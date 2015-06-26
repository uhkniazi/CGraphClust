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

# plot.heatmap.all
plots an annotated heatmap using the NMF library. Requires the object CGraphClust, count matrix with type 1 vertices in the 
rows and type 2 vertices in the columns and a default cutoff expression values (-3, 3) - all values below or above these cutoffs
are thredholded to this value as the extreme values affect the heatmap colouring. The heatmap matrix is scaled by sd and centered 
to a zero mean along the rows (type 1 vertex) before plotting.

# plot.heatmap.mean
similar to plot heatmap but will plot the mean profile of type 1 vertices in each cluster. 

# plot.mean.expressions
plots the mean expression of each group (based on factor fGroups) in each cluster. The function takes the object of CGraphClust type
a count matrix with rows as type 1 vertices and columns as components of the row vectors (they would usually be your samples), 
a grouping factor defining the groups, a legend position (default it topright). The function checks the final projected graph with
type 1 vertices, and uses the names of those vertices to subset the count matrix rows. It assigns cluster labels to each gene and 
uses that cluster label as a grouping factor to calculate the mean values in that group (similar to function plot.heatmap.means) - i.e. it is taking the rowmeans of the count matrix based on each cluster grouping factor. It then takes that matrix of cluster means and
takes the mean and standard deviations for each cluster vector - grouped on the second factor (fGroups). Effectively we are taking
two means of the data, one colMeans of the count matrix (grouped on clusters) and second rowMeans of the matrix (grouped on samples).
The data is plotted and a list with the matrix of mean and standard deviations is returned.

NOTE: Ideally one would want to reassign cluster labels, instead of using the most common label in a group. This can be done
by using the function getClusterMapping, and looking for the type 1 vertices associated with a cluster label. Those type 1 vertex
names (e.g. gene names in that cluster) can be put into reactome database to see which is the most sensible label to assign to that
cluster.
