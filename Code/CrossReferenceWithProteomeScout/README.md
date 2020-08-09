This section requires the ProteomeScoutAPI, please make sure the ProteomeScoutAPI repository is in the same directory as this repository (KinaseSubstratePredictions
)  

- **CrossReferenceWithProteomeScout.ipynb**:Filter the pre-processed predictions to keep only the known human phosphosites based on ProteomeScout(2020-02-07) annotation. The final prediction data is in matrix form. 
- **createSubKinMatrix.py**:convert the list-like dataframe into a matrix
- **proteome_scout_phosphorylation.py**: get the known human phosphosite from ProteomeScout
