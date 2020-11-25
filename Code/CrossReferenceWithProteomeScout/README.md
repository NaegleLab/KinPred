This section requires the ProteomeScoutAPI, please make sure the ProteomeScoutAPI repository is in the same directory as this repository (KinaseSubstratePredictions
)  

- **[CrossReferenceWithProteomeScout.ipynb](https://github.com/NaegleLab/KinPred/blob/master/Code/CrossReferenceWithProteomeScout/CrossReferenceWithProteomeScout.ipynb)**:Filter the pre-processed predictions to keep only the known human phosphosites based on ProteomeScout(2020-02-07) annotation. The final prediction data is in matrix form. 
- **[createSubKinMatrix.py](https://github.com/NaegleLab/KinPred/blob/master/Code/CrossReferenceWithProteomeScout/createSubKinMatrix.py)**:convert the list-like dataframe into a matrix
- **[XRefProteomeScout.py](https://github.com/NaegleLab/KinPred/blob/master/Code/CrossReferenceWithProteomeScout/XRefProteomeScout.py)**: 
  - ```getPScoutData()```: Download, unzip and save the current ProteomeScout Data.\
   Output:   
   download and unzipped ProteomeScout data are saved in `/Data/ProteomeScout_YYYY-MM-DD` directory:
      - data.tsv
      - citations.tsv 
  - ```getHumanPTMs(pscout_data, ref_proteome)```: get the known human phosphosites from downloaded and map to the given reference human proteome.\
    Parameters: 
      - pscout_data : str
          - path to the ProteomeScout data
      - ref_proteome : str
          - path to the reference human proteome (the current KinPred reference human proteome is available at: https://doi.org/10.6084/m9.figshare.12749342.v2, Raw(unfiltered) Data Files for KinPred)

    Return: 
      - df: dataframe of the ProteomeScout data that mapped to the reference human proteome
      
   - ```XRefProteomeScout(pscout_data, ref_proteome, old_version)```:cross referece the (unfiltered) formated prediction data with the mapped human phosphosites. The final formatted file only contains data with confirmed phosphorylation site. Create final substrate/kinase matrices for PhosphoPICK, GPS, and NetworKIN.\
      Parameters:
        - pscout_data : str
            - path to the ProteomeScout data
        - ref_proteome : str
            - path to the reference human proteome
        - old_version: str (YYYY-MM-DD)
            - input (unfiltered) file version

      Output:\
      final standard formated files saved in `../../Data/Final/` directory
