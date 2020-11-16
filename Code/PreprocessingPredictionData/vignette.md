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

Two python modules are used to achieve the goal: 
- [getUniprotID.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/getUniprotID.py): to get the Uniprot accession by query the given `id` (and `id_type`) in Uniprot.org. The Uniprot accession is used as the common asscession for substrate (protein) and kinase identifier.  
  ```def getUniprotID(id, id_type):```\
   Parameters:   
   - id (str): the search term 
   - id_type (str): the search type:
     - gene name (gene/protein name)
     - other (other type of accession)
   
   Return:
   - Uniport accession 
    
- [checkSite.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/checkSite.py): positionally map the given substrate site to the given reference humnan proteome. \
  ```def checkSite(id, site, pep, pos_in_pep, HP_df):```\
    Parameters:
    - id (str):  substrate accession (uniprot accession)
    - site (str or int):
        - str : aa + postion in protein sequence
        - int: postion in protein sequence
    - pep (str): given peptide sequence around the site
    - pos_in_pep (int): position in the peptide
    - HP_df (df): the dataframe of the reference human proteome sequence
    
    Returns:
    - new_site (str): site mapped to the given reference human proteome in the format of amino acid + postion in protein sequence
    - site_confirm (bool)
    

Please see [gps_convert.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/gps_convert.py), [phosphoPick_convert.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/phosphoPick_convert.py), and [networKin_convert.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/networKin_convert.py) as example to apply the above modules while iterate thought predicted substrates and sites in multiple files. 
