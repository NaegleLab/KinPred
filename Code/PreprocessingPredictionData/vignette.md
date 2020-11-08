## Introduction

This vignette shows how to use KinPred code to:
1) replicate the preprocessed prediction data for PhosphoPICK, NetworKIN, and GPS
2) incorporate the provided python modules to process the prediction data of other kinase-substrate predictor to the standarlized format. 

## Replicate KinPred preprocessed data for PhosphoPICK, NetworKIN, and GPS
Please see [FormattingPhosphoPICK.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingPhosphoPICK.ipynb), 
[FormattingNetworKIN.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingNetworKIN.ipynb), 
and [FormattingGPS.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingGPS.ipynb), 
and follow the steps to replicate the KinPred preprocessed data.

## Preprocess data from other kinase-substrate predictor to the standarlized format
The goal of the KinPred preprocessing step is to standerlized the data format across all predictors. The standard formatted file includes all 
predictions on human kinase-substrate interaction in which the predicted sites are mapped to the fixed version of human proteome. 
The standard formatted file contains information of unique IDs for the predicted phosphorylation site (substrate protein accession + position in protein seq),
gene name for the substrates, Uniprot accessions for the substrates, site (aa + position in protein seq), peptide sequences around the sites, scores, 
common kinase names use across all predictors.

Two python modules are used to achieve the goal: [getUniprotID.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/getUniprotID.py) 
and [checkSite.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/checkSite.py)
