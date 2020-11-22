## Introduction

This vignette shows how to use KinPred code to:
1) replicate the preprocessed prediction data for PhosphoPICK, NetworKIN, and GPS
2) incorporate the provided python modules to process the prediction data of other kinase-substrate predictor to the standardized format. 

## Replicate KinPred preprocessed data for PhosphoPICK, NetworKIN, and GPS
Please see [FormattingPhosphoPICK.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingPhosphoPICK.ipynb), 
[FormattingNetworKIN.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingNetworKIN.ipynb), 
and [FormattingGPS.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/FormattingGPS.ipynb), 
and follow the steps to replicate the KinPred preprocessed data.

## Preprocess data from other kinase-substrate predictor to the standardized format
The goal of the KinPred preprocessing step is to standardize the data format across all prediction algorithms. The standard formatted file includes all 
predictions on human kinase-substrate relationships where the predicted sites are mapped to the fixed version of human proteome. 
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
   - Uniprot accession 
    
- [checkSite.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/checkSite.py): positionally map the given substrate site to the given reference human proteome. \
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
    
Use the following code/template to process data from other predictors. 
1. create a python helper module with the following template (Please see [gps_convert.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/gps_convert.py), [phosphoPick_convert.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/phosphoPick_convert.py), and [networKin_convert.py](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/networKin_convert.py) for example):
```
import pandas as pd
import re
import os
import glob
import checkSite, getUniprotID

def convert_acc (filename):
    """
    Given a filename of raw prediction results ,the file is converted to a dataframe and returned 
    
    Parameters
    ----------
    filename : str
        
    Returns
    ----------
    pandas df: input columns +
                'substrate_acc'     substrate uniprotID
                'kinase_acc'        kinase uniprotID
        
    """
    ### load data
    df_raw = pd.read_csv(filename, sep='\t') 
    ### retrieve the current uniprotID for the predicted substrates
    # sub_col = the substrarte col name of your data, it can be the gene name or accession 
    df_unique_sub = df_raw[[sub_col]]
    df_unique_sub = df_unique_sub.drop_duplicates()
    # type = 'gene name' or 'other
    df_unique_sub['substrate_acc'] = df_unique_sub.apply(lambda row :  getUniprotID.getUniprotID(str(row['Uniprot-Acc']),type), axis = 1) 
    ### merge the 'substrate_acc' col to df_raw
    df_raw = df_raw.merge(df_unique_sub, left_on=[sub_col], right_on=[sub_col], how = 'right')

    ### retrieve the current uniprotID for the predicted kinases
    # kin_col = the kinase col name of your data, it can be the gene name or accession  
    df_unique_kin = df_raw[[kin_col]]
    df_unique_kin = df_unique_kin.drop_duplicates()
    # type = 'gene name' or 'other
    df_unique_kin['kinase_acc'] = df_unique_kin.apply(lambda row :  getUniprotID.getUniprotID(str(row['kinase']), type), axis = 1) 
    ### merge the 'kinase_acc' col to df_raw
    df_raw = df_raw.merge(df_unique_kin, left_on=[kin_col], right_on=[kin_col], how = 'left')

    return df_raw

def map_site (filename,  human_proteome):
    """
    Given prodiction dataframe file, positionally map the site to the given human proteome sequence dataframe
    
    Parameters
    ----------
    filename : str
    human_proteome : str
        given human proteome sequence dataframe
        
    Returns
    ----------
    pandas df: input columns +
                'substrate_id'      unique IDs for the substrate phosphorylation site (substrate_acc + position)
                'site'              aa + position in protein sequence
                'pep'               +/- 7 AA
        
    """
    
    ### load file
    df_temp = pd.read_csv(filename)

    # get a list of unique substrate/site 
    # pos_col = predicted site (can be position or aa + position )
    df_unique_sub = df_temp[['substrate_acc', pos_col]]
    df_unique_sub = df_unique_sub.drop_duplicates()
    
    for index, row in df_unique_sub.iterrows():
        id = df_unique_sub.at[index, 'substrate_acc']
        try: 
            # get the sequence of the given uniprotID
            seq = checkSite.getSequence(id, human_proteome)
            pos = int(df_unique_sub.at[index, pos_col])
            # get the substrate site (aa + pos), add the "site" column
            # if the "site" column name already exist, rename the orginal "site" column  
            df_unique_sub.at[index, 'site'] = seq[pos - 1] + str(pos)
            # get the +/- 7 AA around the site from the given reference human proteom 
            pep = checkSite.getPep(pos, 8, seq)
            # add the "pep" column, if the input file doesn't have 
            df_unique_sub.at[index, 'pep'] = pep
            # create unique id for the predicted sites
            df_unique_sub.at[index, 'substrate_id'] = id + '_' + str(pos)
        # except any outdated records such as: 
        # . deleted uniprotID
        # . updated UniprotID
        # . sequence of the given uniprotID changed causing out of range error
        except:
            df_unique_sub.at[index, 'substrate_id'] = 'outdated'

    # add the 'substrate_id' df_temp 
    df_temp = df_temp.merge(df_unique_sub, left_on=['substrate_acc', pos_col], right_on=['substrate_acc', col], how = 'left')
          
    return df_temp

def pick_convert_directory(load_directory, reference_filename, save_directory, convert_type):
    """
    Pulls all files from loading directory and save the output dataframe to the given saving directory
    
    Parameters
    ----------
    load_directory : str
        location to pull files from 
        (directory of prediction data that need to process)
    reference_filename : str
        filename of the referece human proteome sequence dataframe
        OR, 'na'
    save_directory : str
        location to save the processed files result dataframe
    convert_type : str
        acc (convert substrate/kinase accession)
        OR, site (mapping site to reference sequence)
    """
    # for converting substrate and kinase accession
    if convert_type == 'acc':
        all_files = glob.glob(load_directory + "*.txt") # load files (replace with you file extention)
        for filename in all_files:
            name = re.search(r'.+\/(.+)$', filename).group(1) # get the file name
            df = convert_acc(filename) # convert accession
            df.to_csv(save_directory + name + '_mappedAcc.csv', index=False) # save converted file

    # for mapping the site to the reference human proteome seq
    elif convert_type == 'site':
        HP_df = pd.read_csv(reference_filename, sep = '\t') # load referece human proteome sequence dataframe

        all_files = glob.glob(load_directory + "*.txt") # load files (replace with you file extention)
        for filename in all_files:
            name = re.search(r'.+\/(.+)\.txt$', filename).group(1) # get the file name
            df = map_site(filename, HP_df) # map the site to the given reference human proteome seq
            df.to_csv(save_directory + name + '_mappedSite.csv', index=False) # save converted file
 ```
       
2. use `convert_directory(load_directory, save_directory, "acc")` to retrieve the uniprot accession of all substrate protein accession and kinase accession. Make sure your prediction data only contain predictions on human kinase-substrate interaction with at least the following columns
  - substrate (id or protein name)
  - kinase (id or protein name)
  - site (position in the protein sequence, or aa + position)
  - score

3. use `convert_directory(load_directory, reference_filename, save_directory, "site")` to map the sites to the given reference human proteom sequence dataframe. You can convert the reference human proteom sequence fasta file to a dataframe using `humanProteomesReference.fastaToCSV(HP_fasta, HP_csv)` (see [CreateHumanProteomeDF.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/PreprocessingPredictionData/CreateHumanProteomeDF.ipynb) for detail)

4. use pandas to format the dataframe, the preprecessed file should have the following columns:
- substrate_id - unique IDs for the substrate phosphorylation site (substrate_acc + position)
- substrate - gene name for the substrates (if the raw prediction doesn't contain gene name, get the gene name from the reference human proteom sequence dataframe - using uniprot accession)
- substrate_acc - mapped UniprotIDs for the substrates
- site - phosphorylation site
- pep - +/- 7 AA peptide sequence around the site
- kinase - Kinase name (get the common name used across different predictor from the `globalKinaseMap.csv`, available at: https://doi.org/10.6084/m9.figshare.12749333.v1, Kinase Ontology,  using uniprot accession)


