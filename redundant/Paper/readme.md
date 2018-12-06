# Title
Clustering and integrating differentially expressed and co-expressed genes to select genes with high association with the disease.

## Scripts to Generate Data and Results  
1. generate_test_data.R  
  * Downloads two different microarray datasets for TB from GEO database and saves the resulting data for the genes showing differential expression.
2. cluster_dataset_1.R  
  * Clustering of Longitudinal dataset using CGraphClust class.  
3. cluster_dataset_2.R  
  * Clustering of Cross-sectional dataset using CGraphClust class.  
4. merge_graphs_select_biomarker.R
  * Merge the two graphs for the two datasets, selecting genes that are co-expressed in both experiments. Select genes with high hub and degree scores from two clusters, perform variable steps to show these are also good biomarkers.  
  

  




  
  

