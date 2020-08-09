- **1_ChangeAndUsePreferredKinaseNameInAllAnalysis.ipynb**: set and update the preferred kinase for the downstream analysis on the KinPred
- **2_ThresholdsAndFilteringPredictionData.ipynb**: Filter the final predictions by different stringencies
- **3_PredictionNetworkSizeComparison.ipynb**: The size of the networks predicted by different predictors are compared in the following ways:
  - *Total Number of Predicted Edge* by each predictor
  - *Unique Substrates* covered in each predictor
  - *Substrates Overlapped Among Predictors*
  - *Unique Kinase* covered in each predictor
  - *Kinases Overlapped Among Predictors*
  - *Substrate Degree*
  - *Kinase Degree*
- **4_KinaseNetworkSimilarityComparisonWithinPredictor.ipynb**: This notebook compares the similarity of the different-kinase networks predicted by the same predictor. The similarity is measured by Jaccard's Index. The Jaccard’s index is the ratio of the intersection to the union of the two sets. The values range from 0 for no similarity to 1 for the complete overlap of the two compared sets.
- **5.1_KinaseNetworkSimilarityComparisonBetweenPredictors.ipynb**: 
  1. compares the similarity of the kinase networks predicted by the different predictors. The similarity is measured by Jaccard's Index. The Jaccard's index is the ratio of the intersection to the union of the two sets. The values range from 0 for no similarity to 1 for the complete overlap of the two compared sets.
  2. plots the distribution of the same-kinase similarity ranks.
- **5.2_RandomlizedNetworkSimilarity.ipynb**:
  1. randomized the predicted networks with preserved kinase degree
  2. compares the similarity of the kinase networks predicted by the different predictors. The similarity is measured by the Jaccard's Index. The Jaccard's index is the ratio of the intersection to the union of the two sets. The values range from 0 for no similarity to 1 for the complete overlap of the two compared sets.
  3. plots the distribution of the same-kinase similarity ranks.
- **6_StudyBias.ipynb**: This notebook accesses if the extensiveness of the studies on the known phosphosite effects the distribution of:
  1. prediction score
  2. predicted substrate network degrees
- **meltMatrix.py**: create a list-like df given matrix dataframe
- **df_by_kin.py**: divided the given prediction dataframe by kinase type, return two dataframes only contain either Y kinases or S/T kinases
- **kinase_mutual_information.py**:
  - *between_predictors*: calculate the similarity among the preicidted kinase networks from two given predictors by Jaccard's Index.  
  - *within_predictors*: calculate the similarity among the preicidted kinase networks within the given predictor by Jaccard's Index.
- **plot_kinase_heatmap.py**: plot heapmap of given Jaccard's Index matrix
  - *sort_heatmap* sort the kinases in decreasing order of the same-kinase Jaccard's Index, plot the heatmap in the sorted order
  - *reset_label* for S/T kinases, only show label if the kinases:
    - within predictor: have at lease one Jaccard's Index, other than the comparison with itself, higher than 0.49 
    - between predictors: have at lease one Jaccard's Index is at the top 1% 
- **sameKinaseRank**: Sort the Jaccard's Index for each kinases from the two compared predictors (index, columns) in given heatmap.  
- **plot_sameKinase_rank.py**: plot the same-kinase Jaccard's Index rank between the two compared predictors
  - *plot_rank_relation*: scatter plot and histogram 
  - *plot_rank_density*: plot multiple CDFs in as subplots
  - *get_cumelative_count*: get the cumelative fraction kinases of the same-kinase rank
  - *set_up_subplot*: setup and plot subplot
    
