# Title
Clustering and integrating differentially expressed and co-expressed genes to stratify Latent TB into Recents and Remotes

## Scripts to Generate Data and Results  
1. generate_test_data.R  
  * Downloads two different microarray datasets for TB from GEO database and saves the resulting data for the genes showing differential expression.
2. cluster_dataset_1.R  
  * Clustering of Longitudinal dataset using CGraphClust class.  
3. cluster_dataset_2.R  
  * Clustering of Cross-sectional dataset using CGraphClust class.  
4. merge_graphs_subcluster_ltbi.R  
  * Merge the two graphs for the two datasets, selecting genes that are co-expressed in both experiments. Used unsupervised clustering like kmeans to look for 2 sub-groups within the LTB groups using the reduced dimension cluster profile matrix.
  




  
  

