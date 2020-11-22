Welcome to the Github repository for the KinPred v1.0 project, which seeks to create sustainable and usable formats of whole human proteome predictions of kinases with site-specific substrates. Details of code are available within the Code and Data directory readme files. For orientation, the likely workflows you may be interested in, include:

# Users of KinPred-formatted Prediction Algorithm Resources


# Developers of kinase-substrate prediction algorithms
KinPred specifically recommends that to create sustainable and maximally usable results, you perform the following steps:
1. Run the algorithm for your set of kinases on the KinPred reference human proteome. Publish this proteome file and *ALL* predictions in a permanent location, using the same Substrate ID and Kinase ID unified ontology. For KinPred v1.0 publication, these whole-proteome formats are published for NetworKIN, GPS, and PhosphoPICK [here](https://figshare.com/projects/KinPred_v1_0/86885 "FigShare KinPred"). See example code for these projects in Code/PreProcessingPredictionData
2. Use ProteomeScout's current reference phosphoproteome to filter raw predictions for known phosphorylation sites. See Code/CrossReferenceWithProteomeScout
3. Provide an easy-to-use matrix format of these outputs, include all edges (as we noted in our [publication](https://www.biorxiv.org/content/10.1101/2020.08.10.244426v1), removing edges based on a stringency cutoff will cause serious limitations in the flexible use of your predictions). Definitions of matrix files can be seen for [Final Data, here](https://figshare.com/projects/KinPred_v1_0/86885 "FigShare KinPred")

# KinPred and Prediction Algorithm Developers: Updating prediction data
There are two types of updates that a kinase-substrate prediction resource will undergo.
### Phosphoproteome Update
When the phosphoproteome has been updated, then use the UpdatingResourceData/UpdateNewPhosphosites.py script to take the full list of whole-proteome predictions and filter on the new phosphoproteome.
### Underlying Reference Proteome Update
The underlying Uniprot database of the human reference proteome often undergoes changes, such that the substrate and kinase relationships of an old reference are now invalid. Since running predictions anew on the entire reference proteome is costly, we have provided scripts to identify changes in the reference proteome that have an effect on predictions (i.e. the sequence or position of a phosphorylation site has changed), and that site position in the ontology is updated, or a subset of sequences must be run and these replace the existing predictions for that protein in the Formatted, list output of prediction resources. See Code/UpdatingResourceData

All the code used in the study "Unifying and Comparing Multiple Kinase-Substrate Prediction Networks for the Human Phosphoproteome" are stored under "Code".  Please see the README.md files inside "Code" for more detail. Generating updates for current phosphoproteome requires the ProteomeScoutAPI available at https://github.com/NaegleLab/ProteomeScoutAPI. 

All the raw data and output files of the study should store under "Data".  Please see the README.md files inside "Data" for more detail on the directory structure.
