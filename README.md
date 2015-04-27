# CGraphClust
Class to create a Hclust object based on shared properties (e.g. pathways) and positive correlation

# Inheritence Structure
CGraph - Parent Class
|
| 
CGraphClust - Child Class

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
to expected ratios. edge betweeness community finding algorithm is used to detect clusters in this graph. NOTE: if the number of 
edges is over 3000 then the program stops as this may not give a solution without running out of memory or a reasonable amount of 
time. A different algorithm is recommended.
# step 5
the community is converted to a hclust object which can be used for plotting and heatmaps. each cluster is also assigned a 
label based on which is the most common type 2 vertex (e.g. pathway) in that cluster and can be accessed via the
getClusterLabels function.

# data accessor functions
getProjectedGraph, getCorrelationGraph, getBipartiteGraph, getFinalGraph (to get the intersected graph), 
getClusterLabels, getHclust, getCommunity


